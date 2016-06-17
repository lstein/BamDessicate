#!/usr/bin/perl

use strict;
use lib './lib','../lib';
use Bio::DB::DamFile;

if ($ARGV[0] =~ /^--?h/) { # help!
    die <<END;
Usage: dessicate.pl <infile.{bam,tam}> <outfile.dam>

This script uses Bio::DB::DamFile to strip the sequence and quality
information from each read of a standard BAM or TAM (text BAM) file
and write it out in compressed "dessicated BAM" format to
<outfile.dam>.

The sequence and quality information can later be added back to the
dessicated BAM file by using the hydrate.pl script, restoring the
original BAM file.
END
}

my $in  = shift;
my $out = shift;

$in && $out or die "Usage: dessicate.pl <infile.{bam,tam}> <outfile.dam>";

my $bd = Bio::DB::DamFile->dessicate($in,$out);

exit 0;

