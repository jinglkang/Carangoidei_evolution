# Single copy genes evolution
## 1. Devergence time estimation
```bash
# 1.1 把所有的单拷贝基因cds序列提取出来
# (base) shichuang001@login01 Fri Nov 06 2026 14:17:46 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03
mkdir Single_Copy_Orthologue_Sequences_cds
perl extract_scgs_cds.pl
```

- extract_scgs_cds.pl
```perl
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

my %seq;
my @nucfas=<~/jlkang/Carangoidei_evolution/Revise/longest_nuc/*.fas>;
foreach my $fa (@nucfas) {
	my $seqid;
	open FA, $fa or die "can not open $fa\n";
	while (<FA>) {
		chomp;
		if (/>/) {
			s/>//;
			$seqid=$_;
		} else {
			$seq{$seqid}.=$_;
		}
	}
}

my @scgfas=<Single_Copy_Orthologue_Sequences/*.fa>;
foreach	my $fa (@scgfas) {
	my ($newfa, $seqid);
	$newfa="Single_Copy_Orthologue_Sequences_cds/".basename($fa);
	open NEW, ">$newfa" or die "can not create $newfa\n";
	open FA, $fa or die "can not open $fa\n";
	while (<FA>) {
		chomp;
		if (/>/) {
			s/>//;
			$seqid=$_;
			my $seqs=$seq{$seqid};
			print NEW ">$seqid\n$seqs\n";
		}
	}
}
```

```bash
# 1.2 比对和剪切cds序列
# (base) shichuang001@login01 Fri Nov 06 2026 14:42:22 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Single_Copy_Orthologue_Sequences_cds
mkdir muscle; mkdir trimal
vi Align_trim.pl; vi Align_trim.sh
sbatch Align_trim.sh
# Submitted batch job 911730
```
**等待Align_trim.pl运行结束**

## 2. prepare input for paml and Hyphy
```bash
# (base) shichuang001@login01 Fri Nov 06 2026 15:50:17 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Orthogroups
mkdir paml_input; cd paml_input

# (base) shichuang001@login01 Fri Nov 06 2026 15:57:30 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Orthogroups/paml_input
vi prepare_input_paml.pl; vi prepare_input_paml_parallel.pl

# 准备ortho_list.txt, correlation.txt
# (base) shichuang001@login01 Fri Nov 06 2026 16:20:30 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Orthogroups/paml_input
perl generate_list_corr.pl
sbatch run_prepare.sh
# Submitted batch job 911781
```
**等待prepare_input_paml.pl运行结束**

- **prepare_input_paml.pl**
```perl
use strict;
use warnings;

use File::Basename;
use Getopt::Long 'HelpMessage';

my $usage=<<_EOH_;;
--------------------------------------------------------------------------------------------------------------------------------
This script is used to prepare the input for PAML
Usage: 
perl prepare_input_paml.pl --input ortho_list.txt --seq_dir . --cor_list correlation.txt --output .
Example:
1. --input:
the first column is the protein sequence id of reference, the other columns are the nucleotide sequences of each species
# pep and nuc are with same id
# orth_id	spe1_nuc	spe2_nuc	spe3_nuc	spe4_nuc	spe5_nuc	spe6_nuc
OG0000014	Apoly_6299	Padel_40993	Daru_99354	Acura_116661	Ocomp_39102	Pmol_166022
OG0000021	Apoly_3749	Padel_10792	Daru_24893	Acura_7918	Ocomp_30861	Pmol_46346
OG0000035	Apoly_14462	Padel_8374	Daru_143421	Acura_125242	Ocomp_155787	Pmol_97813
OG0000047	Apoly_17254	Padel_33286	Daru_13087	Acura_32119	Ocomp_171369	Pmol_47987

2. --nuc_dir:
directory of nuc sequences

3. --pep_dir:
directory of pep sequences

4. --cor_list:
spe1 	spe1.pep.fasta	spe1.nuc.fasta	# the sec and third col of cor_list are pep and nuc sequences of spe1
spe2 	spe2.pep.fasta	spe2.nuc.fasta
... 	...							...
Options:
	--input			the list of orthologous genes, which should include the reference protein sequences id
					and the nucleotide sequences per species
	--nuc_dir		the directory of nucleotide sequences
	--pep_dir		the directory of protein sequences
	--cor_list		the corresponding list between species and sequence fasta file
	--output,-o 	the directory of the output results
	--help,-h 		Print this help
								
					Kang 2021-7-12
--------------------------------------------------------------------------------------------------------------------------------
_EOH_
;

GetOptions('input:s', \ my $input,
	'nuc_dir:s', \ my $nuc_dir,
	'pep_dir:s', \ my $pep_dir,
	'cor_list:s', \ my $cor_list,
	'output:s', \ my $output,
	'help', \ my $help
	);

if ($help || (! $input) || (! $nuc_dir) || (! $pep_dir) || (! $cor_list) || (! $output)) {
	die $usage; # all of these options are mandatory requirements
}

my $pwd=`pwd`;
chomp($pwd);

open COR, "$cor_list" or die "can not open $cor_list\n";
my (%seq_pep, %seq_nuc, %spe);
my @spe; my $i;
while (<COR>) {
	chomp; $i++;
	my @a=split; my ($seq_id, $seq);
	my $species=$a[0];
	push @spe, $species;
	$spe{$i}=$species;
	
	open FIL, "$pep_dir/$a[1]" or die "can not open $pep_dir/$a[1]\n";
	while (<FIL>) {
		chomp;
		if (/>/) {
			s/>//; my @a=split;
			$seq_id=$a[0];
		} else {
			$seq_pep{$i}->{$seq_id} .= $_;
		}
	}

	open FIL, "$nuc_dir/$a[2]" or die "can not open $nuc_dir/$a[1]\n";
	while (<FIL>) {
		chomp;
		if (/>/) {
			s/>//; my @a=split;
			$seq_id=$a[0];
		} else {
			$seq_nuc{$i}->{$seq_id} .= $_;
		}
	}
}

mkdir $output unless (-e $output);

open INPUT, "$input" or die "can not open $input\n";
while (<INPUT>) {
	chomp;
	my @a=split;
	my $orth_id=$a[0];
	for (my $i = 0; $i < @a; $i++) {
		if ($i==0) {
			mkdir "$output/$a[0]" unless (-e "$output/$a[0]");
		} else {
			my $species=$spe{$i};
			my $seq_id=$a[$i];
			my $seq_pep=$seq_pep{$i}->{$seq_id};
			my $seq_nuc=$seq_nuc{$i}->{$seq_id};
			my $file_pep="pep.fasta";
			my $file_nuc="nuc.fasta";
			open PEP, ">>$output/$a[0]/$file_pep" or die "can not create $output/$a[0]/$file_pep\n";
			print PEP ">$species\n$seq_pep\n";
			open PEP, ">>$output/$a[0]/$file_nuc" or die "can not create $output/$a[0]/$file_nuc\n";
			print PEP ">$species\n$seq_nuc\n";
		}
	}
	
	chdir "$output/$a[0]/";
	system("clustalo -i pep.fasta -t Protein -o pep.aln --outfmt=fa"); # clustalo
	system("pal2nal.pl pep.aln nuc.fasta -output fasta >pal2nal.fasta"); # pal2nal: pal2nal.fasta
	system("Gblocks pal2nal.fasta -b4=10 -b5=n -b3=5 -t=c"); # Gblocks: pal2nal.fasta-gb

    open FIL2, "pal2nal.fasta-gb" or die "can not open pal2nal.fasta-gb\n";
    open FIL3, ">final_alignment.phy" or die "can not create final_alignment.phy\n";
    open FIL4, ">final_alignment.fa" or die "can not create final_alignment.fa\n";
    my $SPE; my (%hash1, %hash2);
    while (<FIL2>) {
    	chomp;
    	if (/>/) {
    		s/>//; $SPE=$_;
    	} else {
    		s/ //g; s/\d/-/g;
    		$hash1{$SPE}.=$_;
    	}
    }
    my $i=0;
    foreach my $spe (@spe) {
    	my $seq=$hash1{$spe};
    	for (my $i = 0; $i < length($seq); $i+=3) {
    		my $codon=substr($seq,$i,3);
    		if ($codon=~/taa/i||$codon=~/tga/i||$codon=~/tag/i) {
    			$codon="---"; # replace stop codon to gap
    		}
    		$hash2{$spe}.=$codon;
    	}
    	$i++; if ($i==1) {
    		my $num=@spe; my $len=length($hash2{$spe});
    		print FIL3 "$num  $len\n";
    		print FIL3 "$spe  $hash2{$spe}\n";
    		print FIL4 ">$spe\n$hash2{$spe}\n";
    	} else {
    		print FIL3 "$spe  $hash2{$spe}\n";
    		print FIL4 ">$spe\n$hash2{$spe}\n";
    	}
    }
    # check gaps in the sequences
    open FIL3, "final_alignment.phy";
    open FIL5, ">>$pwd/final_orth_input_paml.txt";
    my $j=-1; my $len;
    while (<FIL3>) {
    	chomp;
    	my @a=split; $j++;
    	$len=$a[1] if $j==0 && $a[1];
    	if ($a[1]=~/n|-/i) {
    		last;
    	}
    	if ($j==@spe && $len && $len>=150) {
    		print FIL5 "$orth_id\n";
    	}
    }
    system("rm pal2nal.fasta-gb");
    chdir "$pwd";
}
```

- **prepare_input_paml_parallel.pl**
```perl
use strict;
use warnings;
use Parallel::ForkManager;

# prepare_input_paml.pl
my $fm=$ARGV[0];
system("mkdir temp");

my @cmds;
open FM, $fm or die "can not open $fm\n";
while (<FM>) {
	chomp; my @a=split;
	my $temp=$a[0].".txt";
	open TEMP, ">temp/$temp" or die "can not create >temp/$temp\n";
	print TEMP "$_\n";
	my $cmd="perl prepare_input_paml.pl --input temp/$temp --nuc_dir ~/jlkang/Carangoidei_evolution/Revise/longest_nuc/ --pep_dir ~/jlkang/Carangoidei_evolution/Revise/longest_pep/ --cor_list correlation.txt --output .";
	push @cmds, $cmd;
}

my $manager = new Parallel::ForkManager(100);
foreach my $cmd (@cmds) {
	$manager->start and next;
	system($cmd);
	$manager->finish;
}
$manager -> wait_all_children;
system("rm -rf temp/");
```

- **generate_list_corr.pl**
```perl
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Parallel::ForkManager;

my @fas=</jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Single_Copy_Orthologue_Sequences/*fa>;
my @spes; my $i;

my $list="ortho_list.txt";
open LIST, ">$list" or die "can not creat $list\n";

my $corr="correlation.txt";
open CORR, ">$corr" or die "can not creat $corr\n";

foreach my $fa (@fas) {
	$i++;
	(my $id)=basename($fa)=~/(.*?)\.fa/;
	my $info=$id."\t";
	open FA, $fa or die "can not open $fa\n";
	while (<FA>) {
		chomp;
		if (/^>/) {
			s/>//;
			(my $spe)=$_=~/(.*?)\|.*/;
			push @spes, $spe if $i==1;
			$info.=$_."\t";
		}
	}
	$info=~s/\s+$//;
	print LIST "$info\n";
}

foreach my $spe (@spes) {
	my $file=$spe.".fas";
	print CORR "$spe\t$file\t$file\n";
}
```

- **run_prepare.sh**
```bash
#!/bin/bash
#SBATCH -J prepare # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 100 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.prepare # 标准输出文件
#SBATCH -e logs/%j.prepare # 错误输出文件
#SBATCH --mail-user=jinglkang@outlook.com

perl prepare_input_paml_parallel.pl ortho_list.txt
```





