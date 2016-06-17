package Bio::DB::DamFile::Iterator;
use strict;

use Bio::DB::DamFile;
use List::BinarySearch qw(binsearch);
use Carp 'croak';

sub new {
    my $self          = shift;
    my $damfile       = shift or die "Usage: Bio::DB::DamFile::Iterator->new(\$damfile [,$\starting_read])";
    my $starting_read = shift;
    my $ending_read   = shift;
    my $object        = bless { damfile     => $damfile,
				lines       => undef,
	                        block_index => 0,
				read_index  => 0,
				ending_read => $ending_read,
                              },ref $self || $self;
    $object->seek($starting_read) or croak "Read `$starting_read` not found"
	if defined $starting_read;
    return $object;
}

sub damfile { shift->{damfile} }

sub seek {
    my $self = shift;
    my $read = shift;

    my $dam  = $self->damfile;

    # default to an invalid seek
    $self->{block_index} = -1;
    $self->{read_index}  = -1;

    my $block_index = $dam->_lookup_block($read);
    defined $block_index or return;

    my $lines       = $dam->_fetch_block($block_index) or return;
    $self->{lines}  = $lines; # cache

    my $key         = "$read\t"; # terminate match at tab
    my $len         = length($key);
    my $read_index  = binsearch {$a cmp substr($b,0,$len)} $key,@$lines;

    defined $read_index or return;

    $self->{block_index} = $block_index;
    $self->{read_index}  = $read_index;

    return 1;
}

sub next_read {
    my $self  = shift;

    $self->{block_index} ||= 0;
    $self->{read_index}  ||= 0;

    if ($self->{block_index} < 0) {
	$self->{block_index} = 0;
	return;
    }

    my $dam  = $self->damfile;

    # actually a two-tier cache here, one in memory, and one on disk
    my $lines = $self->{lines} 
            ||= $dam->_fetch_block($self->{block_index}) or return;
    
    if (@$lines <= $self->{read_index}) { # refresh cache
	    # otherwise we get next block
	$self->{block_index}++;
	$self->{read_index} = 0;
	$lines   = $self->{lines} = $dam->_fetch_block($self->{block_index}) or return;
    }

    my $next = $lines->[$self->{read_index}++];
    if (defined $self->{ending_read}) {
	my ($id) = $next =~ /^(\S+)\t/;
	if ($id gt $self->{ending_read}) {
	    $self->{block_index} = -1;
	    return;
	}
    }
    return $next;
}

sub reset {
    my $self = shift;
    $self->{block_index} = 0;
    $self->{read_index} = 0;
}

1;
