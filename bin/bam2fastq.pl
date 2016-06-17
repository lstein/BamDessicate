#!/usr/bin/perl

use strict;
open OUT,"|gzip -c" or die $!;

while (@ARGV) {
    my $file = shift;

    my $fh;
    if ($file =~ /\.bam$/) {
	open $fh,"samtools view $file | " or die $!;
    } else {
	open $fh,$file or die $!;
    }
    while (<$fh>) {
	chomp;
	print OUT "$_\n" if /^@/;
	my @fields = split "\t";
	my ($id,$seq,$qual) = @fields[0,9,10];
	print OUT "\@$id\n$seq\n+\n$qual\n";
    }
}
