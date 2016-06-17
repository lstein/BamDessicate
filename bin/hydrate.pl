#!/usr/bin/perl

use strict;
use lib './lib','../lib';
use Bio::DB::DamFile;

if ($ARGV[0] =~ /^--?h/) { # help!
    die <<END;
Usage: hydrate.pl <infile.dam> <readfile.{bam,tam,fastq} <outfile.bam>

This script uses Bio::DB::DamFile to add back the sequence and quality
information to a "dessicated BAM" (.dam) file, creating a proper BAM
file as the result. The read and quality information is provided by
the BAM binary (.bam), BAM text (.tam), or FASTQ file given in the
second command-line argument. The FASTQ file may be gzip or bzip2
compressed.

The output BAM file will be sorted by read name, not by map position.
END
}

@ARGV == 3 or die "Usage: hydrate.pl <in.dam> <reads.{sam,bam,fastq}> <out.bam>";

my $dam_in   = shift;
my $reads_in = shift;
my $bam_out  = shift;

my $bd = Bio::DB::DamFile->new($dam_in);
$bd->rehydrate($reads_in,$bam_out);

0;



