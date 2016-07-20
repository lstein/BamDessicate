#!/usr/bin/perl

use strict;
use lib './lib','../lib';
use Bio::DB::DamFile;
use Getopt::Long;

my (@tmpdirs);
GetOptions('tmpdir|temporary-directory|T=s' => \@tmpdirs
    ) or die <<END;
Usage: dessicate.pl [options] <infile.{bam,tam}> <outfile.dam>

This script uses Bio::DB::DamFile to strip the sequence and quality
information from each read of a standard BAM or TAM (text BAM) file
and write it out in compressed "dessicated BAM" format to
<outfile.dam>.

The sequence and quality information can later be added back to the
dessicated BAM file by using the hydrate.pl script, restoring the
original BAM file.

Options:

   -T,--tmpdir,--temporary-directory=DIR    Use DIR for temporary files. Multiple
                                             options specificy multiple directories.
END

my $in  = shift;
my $out = shift;

$in && $out or die "Usage: dessicate.pl <infile.{bam,tam}> <outfile.dam>";

@tmpdirs = split(/,/,join(',',@tmpdirs));
my $dam = Bio::DB::DamFile->new(undef,{tmpdir=>\@tmpdirs});

$dam->dessicate($in,$out);

exit 0;

