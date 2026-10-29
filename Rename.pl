#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

my $genome=$ARGV[0]; my $nb=$ARGV[1];
my %hash; my $name;
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
        $i++;
        if ($i<=$nb) {
                my $chr="Chr"."$i";
                print ">$chr\n$hash{$id}\n";
        } else {
                $j++;
                my $sca="Scaffold"."$j";
                print ">$sca\n$hash{$id}\n";
        }
}
