#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

my $genome=$ARGV[0]; my $nb=$ARGV[1]; my $gff=$ARGV[2];
my (%hash, %hash1); my $name;
open GENOME, $genome or die "can not open $genome\n";
while (<GENOME>) {
        chomp;
        if (/>/) {
                s/>//;
                $name=$_;
        } else {
                $hash{$name}.=$_;
        }
}

my @sorted = sort {length($hash{$b}) <=> length($hash{$a})} keys %hash;
my ($i, $j);
foreach my $id (@sorted) {
        my @a=split /\s+/, $id;
        $i++;
        if ($i<=$nb) {
                my $chr="Chr"."$i";
                $hash1{$a[0]}=$chr;
                #print ">$chr\n$hash{$id}\n";
        } else {
                $j++;
                my $sca="Scaffold"."$j";
                $hash1{$a[0]}=$sca;
                #print ">$sca\n$hash{$id}\n";
        }
}

open GFF, $gff or die "can not open $gff\n";
while (<GFF>) {
        chomp;
        my @a=split /\t/;
        if (/^#/) {
                print "$_\n";
        } else {
                my $info;
                for (my $i = 1; $i < @a; $i++) {
                        $info.=$a[$i]."\t";
                }
                $info=~s/\s+$//;
                print "$hash1{$a[0]}\t$info\n";
        }
}
