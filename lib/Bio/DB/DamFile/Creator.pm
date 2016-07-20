package Bio::DB::DamFile::Creator;
use strict;

use Bio::DB::DamFile::Common;
use IO::Compress::Bzip2;
use Cwd qw(abs_path);

sub new {
    my $class   = shift;
    my $damfile = shift or die "Usage: Bio::DB::DamFile::Creator->new(\$output_dam_file)";
    my $tmpdirs = shift;
    return bless {damfile=>$damfile,
		  tmpdirs=>$tmpdirs
    },ref $class || $class;
}

sub damfile { shift->{damfile} }

sub dessicate {
    my $self   = shift;
    my $infile = shift;

    $infile or die "usage: \$bd->dessicate(\$sam_or_bamfile_in)";

    my $damfile      = $self->damfile;

    my $outfh        = $self->create_damfile($infile,$damfile);
    my $header_start = HEADER;

    # Transcribe the SAM/BAM @header information. 
    # The two formats need to be treated slightly differently.
    my $is_bam = $infile =~ /\.bam$/;
    if ($is_bam) {
	$self->transcribe_bam_header($infile,$outfh);
    } else {
	$self->transcribe_sam_header($infile,$outfh)
    }

    # open a read name-sorted filehandle
    my $infh  = $self->open_sam_or_bam($infile);

    # now start dessicating blocks
    my $block_start = tell($outfh);
    my $index       = $self->transcribe_blocks($infh,$outfh);

    # now we write out the name index
    my $index_start = tell($outfh);
    $self->write_dam_index($outfh,$index);

    # and update the header to indicate correct starting points
    $self->update_dam_header($outfh,$header_start,$block_start,$index_start);
}

sub hydrate {
    my $self = shift;
    my ($readfile,$damfile) = @_;
    $readfile && $damfile or die "usage: \$bd->hydrate(\$sam_bam_or_fastqfile_in,\$damfile_in,\$bamfile_out)";
}

# See Bio::DB::DamFile::Common for more information on
# the header format.
sub create_damfile {
    my $self                  = shift;
    my ($sourcefile,$outfile) = @_;

    my $abspath  = abs_path($sourcefile);
    my $limit    = HEADER - ( 4   # magic number
			     +4   # version number
			     +8   # offset to beginning of BAM/SAM header data
			     +8   # offset to beginning of gzip block data
			     +8   # offset to beginning of read name index
			     +1 );# room for zero-terminated string containing pathname

    my $pathlen  = length($abspath);
    die "Absolute path '$abspath' is too long (path length=$pathlen bytes; limit=$limit bytes)"
	unless $pathlen <= $limit;  # (this would be a VERY long pathname, but might happen)

    open my $outfh,'>',$outfile or die "$outfile: $!";

    # magic number, offsets to beginning of BAM header data, gzip data, index
    print $outfh pack(HEADER_STRUCT,MAGIC,FORMAT_VERSION*100,HEADER,0,0,$abspath); 

    seek($outfh,HEADER,0);                        # start writing at beginning of BAM header area
    return $outfh;
}

sub transcribe_bam_header {
    my $self             = shift;
    my ($bamfile,$damfh) = @_;
    open my $infh,"samtools view -H $bamfile | " or die "samtools view $bamfile: $!";
    print $damfh $_ while <$infh>;
    close $infh;
}

sub transcribe_sam_header {
    my $self             = shift;
    my ($samfile,$damfh) = @_;

    open my $infh,'<',$samfile or die "$samfile: $!";
    while (<$infh>) {
	last unless /^\@/;
	print $damfh $_;
    }
    close $infh;
}

sub tmpdir_string {
    my $self = shift;
    my $tmpdirs = $self->{tmpdirs} or return '';
    if (ref $tmpdirs) {
	return join ' ',map {"-T $_"} @$tmpdirs;
    } else {
	return "-T $tmpdirs" 
    }
}

sub open_sam_or_bam {
    my $self   = shift;
    my $infile = shift;

    my $infh;
    my $tmpdir = $self->tmpdir_string;
    if ($infile =~ /\.bam$/) {
	open $infh,"samtools view $infile | sort $tmpdir -k1,1 | "  or die "samtools view $infile: $!";
    } else {
	open $infh,"sort $tmpdir -k1,1 $infile                 |"  or die "sort $infile: $!";
    }
    return $infh;
}

sub transcribe_blocks {
    my $self          = shift;
    my ($infh,$outfh) = @_;

    my ($index_id,$block_buffer,$offset,$block_index);
    $offset = tell($outfh);

    while (<$infh>) {
	next if /^@/; # ignore headers

	my @fields = split "\t";
	my $line   = join("\t",@fields[0,1,2,3,4,5,6,7,8,11..$#fields]); # everything but the read and quals

	# the key keeps track of the first ID in the block
	# we constrain IDs to always fall in the same block
	$index_id    = $fields[0] if !defined $index_id;

	if ( ($fields[0] ne $index_id) 
	     && (length($block_buffer) + length($line) > BLOCKSIZE)) {

	    $self->update_read_index($index_id,$offset,\$block_index);
	    $self->write_compressed_block($outfh,\$block_buffer);
	    $offset       = tell($outfh);
	    $block_buffer = $line;
	    $index_id     = $fields[0];
	} else {
	    $block_buffer .= $line;
	}
    }

    # last block
    $self->update_read_index($index_id,$offset,\$block_index);
    $self->write_compressed_block($outfh,\$block_buffer);
    
    # provide dummy key at the very end to enable length retrieval
    $self->update_read_index('~',tell($outfh),\$block_index);

    return \$block_index;
}

sub write_dam_index {
    my $self           = shift;
    my ($outfh,$index) = @_;
    $self->write_compressed_block($outfh,$index);
}

sub update_dam_header {
    my $self = shift;
    my ($outfh,$header_start,$block_start,$index_start) = @_;
    seek($outfh,8,0);
    print $outfh pack('QQQ',$header_start,$block_start,$index_start) 
	or die "write of output file failed: $!";
    close $outfh or die "close of output file failed: $!";
}

sub update_read_index {
    my $self    = shift;
    my ($read_name,$offset,$index) = @_;
    $$index .= pack('Z*Q',$read_name,$offset);
}

sub write_compressed_block {
    my $self          = shift;
    my ($outfh,$data) = @_;   # $data is a scalar reference
    my $compressed_stream = IO::Compress::Bzip2->new($outfh);
    $compressed_stream->print($$data) or die "write compressed stream failed: $!";
    $compressed_stream->close()       or die "compression failed: $!";

}
