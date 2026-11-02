#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Parallel::ForkManager;

# fission
my @cmds;
my @gffs=<*.gff3>;
# agat_sp_keep_longest_isoform.pl --gff Alectis_indicus.gff3 -o longest_gff/Alectis_indicus.longest.gff3
foreach my $gff (@gffs) {
	(my $spe)=$gff=~/(.*)\.gff3/;
	my $new=$spe.".longest.gff3";
	my $cmd="agat_sp_keep_longest_isoform.pl --gff $gff -o longest_gff/$new";
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
