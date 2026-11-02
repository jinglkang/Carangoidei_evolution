#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Parallel::ForkManager;

# gffread -g $genome -x $genome.cds -y $genome.pep $genome.gff3
my @cmds;
my @gffs=<*.gff3>;
foreach my $gff (@gffs) {
	(my $spe)=$gff=~/(.*)\.gff3/;
	my $new=$spe.".longest.gff3";
	my $genome=$spe.".fasta";
	my $cmd="gffread -g $genome -x longest_nuc/$spe.cds -y longest_pep/$spe.pep longest_gff/$new";
	# print "$cmd\n";
	push @cmds, $cmd;
}

my $manager = new Parallel::ForkManager(32);
foreach my $cmd (@cmds) {
        $manager->start and next;
        system($cmd);
        $manager->finish;
}
$manager -> wait_all_children;
