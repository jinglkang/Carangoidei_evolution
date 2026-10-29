#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

my $bed=$ARGV[0]; my (%hash1, %hash2);
open BED, $bed or die "can not open $bed\n";
while (<BED>) {
        chomp;
        my @a=split;
        $hash1{$a[3]}=$a[0];
}

my $simp=$ARGV[1];
my @chrs=split /\,/, $ARGV[2];
foreach my $chr (@chrs) {
        $hash2{$chr}++;
}

open SIMP, $simp or die "can not open $simp\n";
while (<SIMP>) {
        chomp;
        my @a=split;
        my $chr=$hash1{$a[0]};
        if ($hash2{$chr}) {
                print "m*$_\n";
        } else {
                print "$_\n";
        }
}
