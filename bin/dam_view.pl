#!/usr/bin/perl

use strict;
use lib './lib','../lib';
use Bio::DB::DamFile;

if ($ARGV[0] =~ /^--?h/) { # help!
    die <<END;
Usage: dam_view.pl <infile.dam> [start_read] [end_read]

Prints the SAM representation of the specified DAM file. If desired,
you may provide an inclusive range of read IDs to display.

END
}

my $damfile = shift or die "Usage: dam_view.pl <infile.dam> [start_read] [end_read]";
my $dam     = Bio::DB::DamFile->new($damfile);

print $dam->sam_header;

my $iterator = $dam->read_iterator(@ARGV);
while (my $read = $iterator->next_read) {
    print $read,"\n";
}

