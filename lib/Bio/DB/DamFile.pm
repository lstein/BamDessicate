package Bio::DB::DamFile;

our $VERSION = '1.01';

=head1 NAME

Bio::DB::DamFile -- Create and manage dessicated BAM sequence read files

=head1 SYNOPSIS

 # Create a new DAM file from a BAM file.  "reads.bam" is the existing BAM
 # file.  "dessicated.dam" is the new dessicated file.  There is
 # typically an 8x reduction in size.

 $dam = Bio::DB::DamFile->new;
 $dam->dessicate('reads.bam','dessicated.dam');

 # Open an existing DAM file:
 $dam = Bio::DB::DamFile->new('dessicated.dam');

 # Restore a BAM file from a DAM file. The source of reads can be a
 # BAM, TAM or FASTQ file. The mapping information from the DAM file is
 # retained.
 # "reads.bam" is a source of read and quality information; can be
 #             BAM, TAM or FASTQ
 # "hydrated.bam" is the reconstituted BAM file

 $dam->rehydrate('reads.bam','hydrated.bam')

 # Fetch SAM lines one at a time from the DAM file.
 # Reads will come out in alphabetic order by read_id.

 while ($sam_line = $dam->next_read) {
    print $sam_line,"\n";
 }

 # More sophisticated: use an iterator to seek into
 # the file at the location of a particular read.
 # Reads will come out in alphabetic order by read_id,
 # starting with the provided read id.
 $iterator = $dam->read_iterator('NA06984-SRR006041.1000244');
 while ($sam_line = $iterator->next_read) {
    print $sam_line,"\n";
 }

 # Return all SAM lines that match a particular id
 $sam_lines = $dam->fetch_read('NA06984-SRR006041.1000244');
 print join ("\n",@$sam_lines),"\n";

 # Get the SAM header information
 $header    = $dam->sam_header;  # BUG - needs to be written

 # Accessing information in the DAM header
 $magic     = $dam->header_magic;  # DAM magic number (DAM1)
 $path      = $dam->source_path;   # Path to the original BAM file that was used to create the DAM file
 $offset    = $dam->header_offset; # Offset (in bytes) to where the SAM header starts
 $offset    = $dam->block_offset;  # Offset (in bytes) to where the bzip2-compressed alignmenet data starts
 $offset    = $dam->index_offset;  # Offset (in bytes) to where the name-sorted index of reads begins

 # COMMAND LINE TOOLS:
 $ dessicate.pl in.bam out.dam           # dessicates in.bam, stores it in out.dam
 $ hydrate.pl   in.dam in.bam   out.bam  # rehydrates in.dam from read data stored in in.bam and stores in ou t.ba
 $ hydrate.pl   in.dam in.fastq out.bam  # same idea, but takes read data from in.fastq (works on fastq.gz too)
 $ dam_view.pl  in.dam                   # extracts SAM lines from in.dam and displays to stdout
 $ dam_view.pl  in.dam read0010 read0020 # displays SAM lines from start read to end read

=head1 DESCRIPTION

This module was created to solve the issue of maintaining multiple
aligned BAM files from the same set of reads. This happens when a BAM
file is remapped onto different genome builds or using different
alignment software/settings. Rather than have multiple copies of the
same read and quality score information, one wishes to maintain a
single BAM (or FASTQ) file with the read information, and store the
alternative alignments in separate data files.

The DAM ("dessicated BAM") format is very simple, and consists of the
standard SAM header followed by a series of bzip2-compressed chunks of
alpha-sorted SAM (text format) alignment lines from which the read and
quality score fields have been removed. This is followed by a
alphabetically-sorted index of the first read name in each block
followed by the offset to that block. An individual read line can be
found by performing a binary search on the read name index, fetching
and uncompressing the corresponding block, and then performing a
binary search on the SAM lines contained within the uncompressed
block.

The DAM file is typically reduced by a factor of 6-8 relative to the
size of the original BAM file, making it an economical alternative to
storing multiple remapped BAM files.

Here is an ASCII text representation of the DAM file:

 --------------------------------------------------------------------------------------
 |                       HEADER  (512 bytes)                                          |
 |DAM1                -- magic number,                                      4 bytes   |
 |SAM_header_offset   -- byte position of the start of the SAM header,      16 bytes  |
 |block_offset        -- byte position of the first compressed block,       16 bytes  |
 |index_offset        -- byte position of the start of the read name index, 16 bytes  |
 |path_name           -- full path to the source BAM/SAM file                variable |
 |                        (this is a zero-terminated string)                          |
 --------------------------------------------------------------------------------------
 |                      SAM HEADER (variable length)                                  |
 | Uncompressed text version of the SAM/BAM header                                    |
 --------------------------------------------------------------------------------------
 |  BLOCK 1 (variable length) -bzip2 compressed SAM lines, newline terminated,        |
 |                              up to 1 MB in length prior to compression.            |
 --------------------------------------------------------------------------------------
 |  BLOCK 2 (variable length)                                                         |
 --------------------------------------------------------------------------------------
 ~
 --------------------------------------------------------------------------------------
 |  BLOCK N (variable length)                                                         |
 --------------------------------------------------------------------------------------
 |  READ Index (variable length) -- 
              bzip2-compressed data. When uncompressed looks like this:               |  
 | ReadName1\0Offset
 | ReadName2\0Offset
 |...
 | ReadNameN\Offset
 |
 | ReadName is a null terminated string. Offset is the position of the compressed block
 | that begins with the ReadName SAM line. Offset is 8 bytes in length.
 --------------------------------------------------------------------------------------
 

=head1 METHODS

=cut

use strict;

use Bio::DB::DamFile::Common;
use Bio::DB::DamFile::Iterator;

use IO::Uncompress::Bunzip2 qw(bunzip2);
use List::BinarySearch      qw(binsearch_pos binsearch);
use Carp                    qw(croak);
use Tie::Cache;             # SHOULD USE Tie::Cache::LRU FOR SPEED

=head2 $dam = Bio::DB::DamFile->new(['path_to_damfile'] [,\%options])

This is the constructor for the class. When called without arguments
it can be used to create a new DAM file using the dessicate() method:

  my $dam = Bio::DB::DamFile->new();
  $dam->dessicate('input.bam','output.dam');

This will create the dessicated BAM file, "output.dam", and then
open up the DamFile for reading and further manipulation.

Alternatively, when called with the path to an existing .dam file, the
method will open it for reading:

 my $dam = Bio::DB::DamFile->new('/path/to/file.dam');

The optional second argument is a hashref to a series of named
options. Currently only one option, 'cache_size' is recognized:

 $dam = Bio::DB::DamFile->new('/path/to/file.dam',{cache_size => 1_048_576_000})

This controls the size of the least-recently-used cache of
uncompressed SAM data blocks. Roughly 200 data blocks can be held in
memory. For some access patterns (fetching reads that are close
together alphabetically), performance may be affected by the cache
size, but in most cases it won't make a difference.

=cut

sub new {
    my $class    = shift;
    my $damfile  = shift;
    my $options  = shift;
    return bless {
	damfile       => $damfile,
	header_data   => undef,
	damfh         => undef,
	options       => $options || {},
    },ref $class || $class;
}

=head2 $file = $dam->damfile

Read-only accessor that returns the path to the .dam file that was
passed to new().

=cut

sub damfile     { shift->{damfile} }

=head2 $cach_size = $dam->block_cache_size

Read-only accessor that returns the size of the LRU cache in which
uncompressed blocks are stored.

=cut

sub block_cache_size {
    my $self = shift;
    return $self->{options}{cache_size} ||= DEFAULT_BLOCK_CACHE_SIZE;
}

sub tmpdir {
    my $self = shift;
    return $self->{options}{tmpdir};
}

=head2 $magic = $dam->header_magic

Returns the "magic number" that begins the header of the .dam file.

=cut

sub header_magic  {
    my $self = shift;
    $self->{header_data} ||= $self->_get_dam_header;
    $self->{header_data}{magic};
}

=head2 $version = $dam->format_version

Returns the version of the file format. Currently 1.01

=cut

sub format_version  {
    my $self = shift;
    $self->{header_data} ||= $self->_get_dam_header;
    $self->{header_data}{format_version};
}

=head2 $bytes = $dam->header_offset

Returns the number of bytes from the beginning of the file to the
beginning of the (uncompressed) SAM header data.

=cut

sub header_offset  {
    my $self = shift;
    $self->{header_data} ||= $self->_get_dam_header;
    $self->{header_data}{header_offset};
}

=head2 $bytes = $dam->block_offset

Returns the number of bytes from the beginning of the file to the
beginning of the compressed read data.

=cut


sub block_offset {
    my $self = shift;
    $self->{header_data} ||= $self->_get_dam_header;
    $self->{header_data}{block_offset};
}

=head2 $bytes = $dam->index_offset

Returns the number of bytes from the beginning of the file to the
beginning of the compressed read index data.

=cut

sub index_offset {
    my $self = shift;
    $self->{header_data} ||= $self->_get_dam_header;
    $self->{header_data}{index_offset};
}

=head2 $path = $dam->source_path

Returns the absolute path to the source BAM/SAM file from which the
.dam file was derived.

=cut

sub source_path {
    my $self = shift;
    $self->{header_data} ||= $self->_get_dam_header;
    $self->{header_data}{original_path};
}

=head2 $header = $dam->read_header

This method returns the entire SAM header as a string. You must split
and parse it yourself:

    my $header  = $dam->read_header;
    my @entries = split "\n",$h;
    for my $l (@entries) {
       my ($tag,@fields) = split "\t",$l;
       # now do something
    }

=cut

sub sam_header {
    my $self = shift;
    my $sam_start = $self->header_offset;
    my $sam_length= $self->block_offset-$sam_start;
    my $buffer;
    my $fh = $self->_open_damfile;
    seek($fh,$sam_start,0)        or die "Can't seek on dam fh: $!";
    read($fh,$buffer,$sam_length) or die "Can't read on dam fh: $!";
    return $buffer;
}

=head2 $dam->dessicate('in.bam','out.dam')

The dessicate() method is used to create a .dam file from an input BAM
or SAM file. It takes two arguments: the path to the input BAM (or
SAM) file, and the path to the output .dam file. On success, it will
return the Bio::DB::DamFile object attached to the newly-created .dam
file.

It can be called as a class method like so:

  $dam = Bio::DB::DamFile->dessicate('in.bam','out.dam')

The input file can be a compressed or uncompressed .sam file, or a
.bam file.

Internally, the method uses Bio::DB::DamFile::Creator to do the dirty
work by making this call:

     Bio::DB::DamFile::Creator->new('out.bam')->dessicate('in.bam');

=cut

sub dessicate {
    my $self   = shift;
    my ($infile,$damfile) = @_;
    my $obj = ref $self ? $self : $self->new();
    $damfile            ||= eval {$obj->damfile};
    $infile && $damfile   or die "usage: \$bd->dessicate(\$sam_or_bamfile_in,\$damfile_out)";

    eval 'require Bio::DB::DamFile::Creator' 
	unless Bio::DB::DamFile::Creator->can('new');

    Bio::DB::DamFile::Creator->new($damfile,eval{$self->tmpdir})->dessicate($infile);
    $obj->{damfile} = $damfile;
    return $obj;
}

=head2 $dam->rehydrate('in.bam','out.bam')

The rehydrate() method takes the mapping information in its associated
.dam file, adds in read and quality score information from a provided
BAM, SAM or FASTQ file, and produces a BAM file that is equivalent to
the original.

The first argument is a BAM, SAM or FASTQ file that contains read and
quality score information. The SAM and FASTQ files may be compressed
or uncompressed. The second argument is the destination BAM file.

Here is an example of roundtripping from BAM to DAM to BAM again:

    my $dam = Bio::DB::DamFile->dessicate('original.bam','test.dam');
    $dam->rehydrate('original.bam','output.bam');

This will create a file "output.bam" that contains the same
information as "original.bam". The dessicated data will be in
"test.dam".

=cut

sub rehydrate {
    my $self   = shift;
    my ($infile,$outfile) = @_;

    $infile && $outfile or croak "Usage: Bio::DB::DamFile->rehydrate(\$bam_sam_or_fastq_file_in,\$bamfile_out)";

    open my $outfh,"| samtools view -b -S - > $outfile" or die "Can't open samtools pipe to write $outfile";

    # write out the SAM header
    my $fh = $self->_open_damfile;
    my $offset = $self->header_offset;
    my $len    = $self->block_offset - $offset;
    seek($fh,$self->header_offset,0);
    my $buffer;
    read($fh,$buffer,$len) or die "Couldn't read $len header bytes from DAM file at offset $offset";
    print $outfh $buffer;

    # write out the contents
    # infile can be any of:
    # (1) a bam file (.bam)
    # (2) a sam file (.sam or .tam)
    # (3) a fastq file (.fastq)

    if ($infile =~ /\.bam$/) {
	$self->_rehydrate_bam($infile,$outfh);
    } elsif ($infile =~ /\.[st]am$/) {
	$self->_rehydrate_sam($infile,$outfh);
    } elsif ($infile =~ /\.fastq(?:\.gz|\.bz2)?$/) {
	$self->_rehydrate_fastq($infile,$outfh);
    } else {
	croak "$infile has an unknown extension (must be one of .bam, .sam, .tam or .fastq)";
    }

    close $outfh or die "error closing $outfile: $!";
}

=head2 $read = $dam->next_read(['starting_read'] [,'ending_read'])

Called repeatedly, the next_read() method will return an
alphabetically-sorted list of the SAM read lines in the attached .dam
file. You may optionally provide the ID(s) of the read(s) to start and
end at. When no more reads are available, the method returns undef.

This code fragment iterates through all reads in the file:

  while ( my $read = $dam->next_read() ) {
        my @fields = split "\t",$r;
        print "read id = $fields[0], reference chromosome = $fields[2]\n";
  }


This code fragment iterates through all the reads in the file,
starting with the read named "na12345" (inclusive):

  while ( my $read = $dam->next_read('na12345') ) {
        my @fields = split "\t",$r;
        print "read id = $fields[0], reference chromosome = $fields[2]\n";
  }

And this fragment starts at na12345 and ends at na90000 (inclusive):

  while ( my $read = $dam->next_read('na12345','na90000') ) {
        my @fields = split "\t",$r;
        print "read id = $fields[0], reference chromosome = $fields[2]\n";
  }

If the start and/or end read are not present in the .dam file, then
iteration will span the alphabetic range defined by the read names.

=cut

sub next_read {
    my $self  = shift;
    $self->{iterator} ||= $self->read_iterator(@_);
    my $read = $self->{iterator}->next_read();
    undef $self->{iterator} if !defined $read;
    return $read;
}

=head2 $iterator = $dam->read_iterator(['starting_read'] [,'ending_read'])

The read_iterator() method operates similarly to next_read() except
that it returns an intermediate object, the
Bio::DB::DamFile::Iterator. From this you can call next_read()
repeatedly to retrieve all, or a range of, reads in the file. The main
advantage of this method over calling next_read() directly is that you
can set up multiple simultaneous iterators to return reads across
different alphabetic ranges.  Example:

 my $iterator = $dam->read_iterator('na12345','na90000');
 while (my $read = $iterator->next_read()) {
      # do something
 }

=cut

sub read_iterator {
    my $self = shift;
    return Bio::DB::DamFile::Iterator->new($self,@_);
}

=head2 $read_lines = $dam->fetch_read('read_id')

This method takes a read id and returns a reference to an array of all
the SAM read lines that match. If the read ID is not found in the .dam
file, then this method will throw an exception. You may wish to put it
inside an eval{}:

 my $reads = eval {$dam->fetch_read('ABC')} // croak "read ABC not found";
 for my $r (@$reads) {
    my @fields = split "\t",$r;
    print "read id = $fields[0], reference chromosome = $fields[2]\n";
 }

Note that fetch_read() causes a binary search lookup on the block
index, followed by a binary search on the uncompressed block. An
internal least-recently-used cache speeds up access, but it will only
be effective if you are accessing the reads in alphabetic or
near-alphabetic order. Consider using next_read() or a read iterator
instead.

=cut

sub fetch_read {
    my $self    = shift;
    my $read_id = shift or die "Usage \$readline = Bio::DB::DamFile->fetch_read(\$read_id)";

    my ($block_index) = $self->_lookup_block($read_id)    or croak "Read $read_id not found in DB.";
    my $lines         = $self->_fetch_block($block_index) or croak "Block fetch error while retrieving block at $block_index: $!";

    my $key   = "$read_id\t"; # terminate match at tab

    my $len   = length($key);
    my $i     = binsearch {$a cmp substr($b,0,$len)} $key,@$lines;

    croak "Read $read_id not found in DB." unless defined $i;
    
    # there may be more than one matching read line, but they will be adjacent!
    my @matches;
    while (substr($lines->[$i],0,$len) eq $key) {
	my $line = _interpolate_stars($lines->[$i++]);
	push @matches,$line;
    }
    
    return \@matches;
}

sub _tmpdir_string {
    my $self = shift;

    my $tmpdir_string = '';
    if (my $tmpdir = $self->tmpdir) {
	$tmpdir_string = ref $tmpdir ? (join ' ',map{"-T $_"} @$tmpdir) : "-T $tmpdir";
    }

    return $tmpdir_string;
}

sub _rehydrate_bam {
    my $self = shift;
    my ($infile,$outfh) = @_;

    my $tmpdir = $self->_tmpdir_string;
    open my $infh,"samtools view $infile | sort $tmpdir -k1,1 |" 
	or die "Can't open samtools to read from $infile: $!";
    warn "Sorting input BAM by read name. This may take a while...\n";
    $self->_rehydrate_stream($infh,$outfh);
}

sub _rehydrate_sam {
    my $self = shift;
    my ($infile,$outfh) = @_;
    my $tmpdir = $self->_tmpdir_string;
    open my $infh,"grep -v '\@' $infile | sort $tmpdir -k1,1 |" or die "Can't open $infile: $!";
    warn "Sorting input SAM by read name. This may take a while...\n";
    $self->_rehydrate_stream($infh,$outfh);
}

sub _rehydrate_fastq {
    my $self = shift;
    my ($infile,$outfh) = @_;

    # need to create a stream of sorted SAM-like data for passing to _rehydrate_stream()
    my $pid = open(my $infh,"-|") // die "Can't fork!";

    if ($pid) { # in parent process
	$self->_rehydrate_stream($infh,$outfh);
	return;
    } 

    else {  # in child process
	my $fastq_fh;

	my $tmpdir = $self->_tmpdir_string;

	# open appropriate unzipper
	if ($infile =~ /\.gz$/) {
	    open $fastq_fh,"gunzip -c $infile |"  or die "gunzip   -c $infile: $!";
	} elsif ($infile =~ /\.bz2$/) {
	    open $fastq_fh,"bunzip2 -c $infile |" or die "bunzip2 -c $infile: $!";
	} else {
	    open $fastq_fh,'<',$infile            or die "failed opening $infile for reading: $!";
	}
	open my $sort_fh,"| sort $tmpdir -k1,1"   or die "failed opening output pipe to sort: $!";

	local $/ = '@';
	while (<$fastq_fh>) {
	    chomp;
	    next unless $_;
	    my ($read_id,$dna,undef,$quality) = split "\n";
	    #             field      0        1     2     3     4     5     6     7     8     9    10
	    print $sort_fh join("\t",$read_id,undef,undef,undef,undef,undef,undef,undef,undef,$dna,$quality),"\n";
	}
	close $sort_fh                            or die "An error occurred while closing the sort pipe filehandle: $!";
	exit 0;
    }
}

sub _rehydrate_stream {
    my $self = shift;
    my ($infh,$outfh) = @_;

    my $iterator = $self->read_iterator();

    my @sam_fields = ('');
    my $sam_done;

    while (my $dam_line = $iterator->next_read) {
	my @dam_fields = split "\t",$dam_line;

	while (!$sam_done && ($sam_fields[0] lt $dam_fields[0])) { # read from sam file until we match
	    chomp (my $sam_line = <$infh>);
	    $sam_done++ unless $sam_line;
	    @sam_fields = split "\t",$sam_line;
	    last if $sam_fields[0] ge $dam_fields[0];
	}

	if ($sam_done) {
	    # sequence missing
	    print $outfh join("\t",@dam_fields),"\n";
	} 

	elsif ($dam_fields[0] eq $sam_fields[0]) { #match
	    print $outfh join("\t",@dam_fields[0..8],
			           @sam_fields[9,10],
			           @dam_fields[11..$#dam_fields]),"\n";
	}
    }
}

# This is not a method, but a subroutine call.
# It adds stars to the SAM line for the sequence and quality scores
sub _interpolate_stars {
    my $line   = shift;
    my @fields = split "\t",$line;
    return join ("\t",@fields[0..8],'*','*',@fields[9..$#fields]);
}

sub _open_damfile {
    my $self = shift;
    return $self->{damfh} 
       if defined $self->{damfh} && defined fileno($self->{damfh});

    open my $fh,'<',$self->damfile or die $self->damfile,": $!";
    $self->{damfh} = $fh;
    
    return $self->{damfh};
}

sub _get_dam_header {
    my $self = shift;
    my $fh = $self->_open_damfile;
    my $buffer;
    seek($fh,0,0);
    read($fh,$buffer,HEADER) or die "Couldn't read from ",$self->damfile,": $!";
    my @data = unpack(HEADER_STRUCT,$buffer);

    my %fields;
    @fields{'magic','format_version','header_offset','block_offset','index_offset','original_path'} = @data;

    $fields{magic}          eq MAGIC          or croak $self->damfile," doesn't have the right magic number";

    $fields{format_version} /= 100;
    $fields{format_version} == FORMAT_VERSION 
	or croak $self->damfile," doesn't have the correct format version (got $fields{format_version} but expect ",FORMAT_VERSION,")";
    return \%fields;
}

sub _lookup_block {
    my $self = shift;
    my $key  = shift;  # a read id

    my $index = $self->_get_read_index;

    # find the first block that might contain the key
    my $i     = binsearch_pos {$a cmp $b->[0]} $key,@$index;
    return if $i < 0 or $i > $#$index;

    $i-- unless $index->[$i][0] eq $key;  # distinguish between an exact match and an insert position
    return $i;
}

sub _fetch_block {
    my $self  = shift;
    my $i     = shift;

    my $cache = $self->_block_cache;
    return $cache->{$i} if defined $cache->{$i};

    my $index = $self->_get_read_index;
    my $offset = $index->[$i][1];
    my $length = $index->[$i+1][1]-$offset;
    return unless $length > 0;

    my $fh     = $self->_open_damfile;

    my $block;
    seek($fh,$offset,0)      or die "seek failed: $!";
    read($fh,$block,$length) or die "read failed: $!";

    my $uncompressed;
    bunzip2(\$block,\$uncompressed);
    my @lines           = split "\n",$uncompressed;
    return $cache->{$i} = \@lines;
}

sub _get_read_index {
    my $self = shift;
    return $self->{read_index} if $self->{read_index};

    my $fh   = $self->_open_damfile;
    seek($fh,$self->index_offset,0);

    my $data = '';
    do {1} while read($fh, $data, 8192, length $data);

    my $index;
    bunzip2(\$data,\$index);
    my @flat = unpack('(Z*Q)*',$index);

    # turn into list of lists
    my @index;
    while (my ($key,$offset) = splice(@flat,0,2)) {
	push @index,[$key,$offset];
    }
    return $self->{read_index} = \@index;
}

sub _block_cache {
    my $self = shift;

    return $self->{block_cache} if defined $self->{block_cache};

    my %c;
    tie %c,'Tie::Cache',{MaxBytes => $self->block_cache_size};
    return $self->{block_cache} = \%c;
}

=head1 SEE ALSO

L<Bio::Perl>, L<Bio::DB::Bam>

=head1 AUTHOR

Lincoln Stein E<lt>lincoln.stein@oicr.on.caE<gt>.
E<lt>lincoln.stein@bmail.comE<gt>

Copyright (c) 2016 Ontario Institute for Cancer Research.

This package and its accompanying libraries are free software; you can
redistribute it and/or modify it under the terms of the Artistic
License 2.0, the Apache 2.0 License, or the GNU General Public License
(version 1 or higher).  Refer to LICENSE for the full license text.

=cut


1;
