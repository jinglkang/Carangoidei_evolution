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
sbatch Align_trim.sh # trimal -gt 0.9 :1 - (fraction of sequences with a gap allowed). Range: [0 - 1]
# Submitted batch job 912039
```
**等待Align_trim.pl运行结束**
**运行结束**

```bash
# 1.3 把single copy genes的cds序列连在一块
# (base) shichuang001@login01 Sat Nov 07 2026 11:44:05 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Single_Copy_Orthologue_Sequences_cds/trimal
perl Conca_seq.pl > ../../Phylogeny/Single_copy_conca_cds.fa

# 1.4 把fasta格式转换成phylip格式
# (base) shichuang001@login01 Sat Nov 07 2026 11:53:15 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Phylogeny
perl ~/bin/fasta2phy.pl Single_copy_conca_cds.fa Single_copy_conca_cds.phy # fasta2phy.pl有问题，第一行应该把制表符(\t)改成一个空格(\s)

# 1.5 run mcmctree
# (base) shichuang001@login01 Sat Nov 07 2026 20:41:43 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Phylogeny
mkdir mcmctree;cd mcmctree
# (base) shichuang001@login01 Sat Nov 07 2026 20:42:16 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Phylogeny/mcmctree
cp ../Carangoidei_28sp.contree ./

# 1.5.1 baseml估算substitution rate
# 修改树
# (base) shichuang001@login01 Sat Nov 07 2026 20:58:53 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Phylogeny/mcmctree
less mcmc.tree |perl -alne 's/\:\d+\.\d+//g;s/\d+//g;s/\[//g;s/\]//g;print' > mcmc.tree.1; mv mcmc.tree.1 mcmc.tree
# 需要编辑mcmc.tree: 第一行表明有多少个物种，以及多少个树
# 两个化石时间点: O. latipes–T. nigroviridis (~96.9–150.9 Ma), D. rerio–G. aculeatus (~149.85–165.2 Ma)
# 从itol导出以斑马鱼为root的树，命名为mcmc.tree
vi mcmc.tree
cp ~/software/paml-4.10.10-linux-x86_64/src/baseml.ctl ./
mv ../Single_copy_conca_cds.phy ./
baseml baseml.ctl
```

- **baseml.ctl**
```baseml.ctl
      seqfile = Single_copy_conca_cds.phy * sequence data filename
     treefile = mcmc.tree * tree Structure file name
      outfile = baseml.out           * main result file name

        noisy = 3  * 0,1,2,3,9: how much rubbish on the screen
      verbose = 0  * 0: concise; 1: detailed, 2: too much
      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

        noisy = 3   * 0,1,2,3: how much rubbish on the screen
      verbose = 0   * 1: detailed output, 0: concise output
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI

*        model = 10 [5 (AC CA) (AG GA) (AT TA) (CG GC) (CT TC)] /* SYM */
        model = 7   * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85, 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu * REV = GTR

        ndata = 1  # run ML using subtrees generated from maintree
        Mgene = 0   * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff

        clock = 0   * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
    fix_kappa = 0   * 0: estimate kappa; 1: fix kappa at value below; 2: kappa for branches
        kappa = 5   * initial or fixed kappa

    fix_alpha = 0   * 0: estimate alpha; 1: fix alpha at value below
        alpha = 0   * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * 1: different alpha's for genes, 0: one alpha
        ncatG = 4   * # of categories in the dG, AdG, or nparK models of rates
        nparK = 0   * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK

        nhomo = 0   * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2
        getSE = 1   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 0   * (0,1,2): rates (alpha>0) or ancestral states

   Small_Diff = 1e-8
    cleandata = 1   * remove sites with ambiguity data (1:yes, 0:no)?
*        icode = 0   * (with RateAncestor=1. try "GC" in data,model=4,Mgene=4)
  fix_blength = 0  * 0: ignore, -1: random, 1: initial, 2: fixed, 3: proportional
       method = 0   * Optimization method 0: simultaneous; 1: one branch a time
```

- **mcmc.tree**
```mcmc.tree
32 1
(Danio_rerio,(Oryzias_latipes,((Gasterosteus_aculeatus,Takifugu_rubripes),(Polydactylus_sextarius,(((Xiphias_gladius,Toxotes_jaculatrix),(Lates_calcarifer,(Sphyraena_barracuda,Sphyraena_novaehollandiae))),(((Echeneis_naucrates,Remorina_albescens),(Scomberoides_lysan,(Trachinotus_ovatus,(Trachinotus_anak,(Trachinotus_carolinus,Trachinotus_blochii))))),(Elagatis_bipinnulata,(((Decapterus_macrosoma,Decapterus_maruadsi),(Trachurus_trachurus,Trachurus_japonicus)),((Parastromateus_niger,(Caranx_equula,Alectis_indicus)),((Selaroides_leptolepis,(Alepes_djedaba,Megalaspis_cordyla)),(Caranx_crysos,(Caranx_sexfasciatus,(Caranx_ignobilis,Caranx_heberi))))))))))))'@(.969,1.509)')'@(1.4885,1.652)';
```

- **Conca_seq.pl**
```perl
#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

my (%hash, %hash1); my @spes;
my @fas=<*fa.trimal>;
foreach my $fa (@fas) {
    my (%hash2); my $spe;
    open FA, $fa or die "can not open $fa\n";
    while (<FA>) {
        chomp;
        if (/^>/) {
            s/>//;
            ($spe)=$_=~/(.*?)\|.*/;
            $hash1{$spe}++;
            push @spes, $spe if $hash1{$spe}==1;
        } else {
            $hash2{$spe}.=$_;
            }
    }

    my $len=length($hash2{$spe});
    if ($len >= 150) {
        open FA, $fa or die "can not open $fa\n";
        while (<FA>) {
            chomp;
            if (/^>/) {
                s/>//;
                ($spe)=$_=~/(.*?)\|.*/;
            } else {
                $hash{$spe}.=$_;
            }
        }
    }
}

foreach my $spe (@spes) {
    print ">$spe\n$hash{$spe}\n";
}
```


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
# Submitted batch job 911845
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
	system("~/software/Gblocks_0.91b/Gblocks pal2nal.fasta -b4=10 -b5=n -b3=5 -t=c"); # Gblocks: pal2nal.fasta-gb

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

## 3. Run PAML free-ratio model for evolutionary rate
**codeml.pl: https://github.com/jinglkang/small_tools/blob/master/codeml.pl**
```bash
# (base) shichuang001@login01 Sun Nov 08 2026 16:04:13 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Orthogroups/paml_input
cp ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Phylogeny/Carangoidei_28sp.contree ./
less Carangoidei_28sp.contree|perl -alne 's/\:\d+\.\d+//g;s/\d+//g;s/\[//g;s/\]//g;print' > spe.tree

# 3.1 Run PAML free-ratio
# perl codeml_parallel.pl final_orth_input_paml.txt temp1 free-ratio spe.tree
sbatch run_free-ratio.sh
# Submitted batch job 913416

# 3.2 Run Hyphy for positive selection analysis
# (eggnog) shichuang001@login01 Mon Nov 09 2026 13:43:24 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Orthogroups/paml_input
vi spe_hyphy.tree

# 3.3 逐个谱系进行BUSTED分析
# (base) shichuang001@login01 Mon Nov 09 2026 14:04:28 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Orthogroups/paml_input
vi run_hyphy.pl
# 把三类具有扁平性状的鱼类设置为前景枝: FG1, FG2, FG3; 详情见spe_hyphy.tree
sbatch run_hyphy.sh
# Submitted batch job 913882
```

- **run_hyphy.sh**
```bash
#!/bin/bash
#SBATCH -J hyphy_busted # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 100 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.hyphy_busted # 标准输出文件
#SBATCH -e logs/%j.hyphy_busted # 错误输出文件
#SBATCH --mail-user=jinglkang@outlook.com

perl run_hyphy.pl final_orth_input_paml.txt spe_hyphy.tree
```

- **run_hyphy.pl**
```perl
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;

# BUSTED
# hyphy busted --alignment final_alignment.fa --tree spe_hyphy.tree --branches FG1 --output OG000001.FG1.BUSTED.json
# hyphy busted --alignment final_alignment.fa --tree spe_hyphy.tree --branches FG2 --output OG000001.FG2.BUSTED.json
# hyphy busted --alignment final_alignment.fa --tree spe_hyphy.tree --branches FG3 --output OG000001.FG3.BUSTED.json

my $list=$ARGV[0]; # The list
my $tree=$ARGV[1]; # The tree

my @cmds;
open LIST, $list or die "can not open $list\n";
while (<LIST>) {
	chomp;
	my @a=split;
	my $orth=$a[0];
	my $alig="$orth/final_alignment.fa"; # the alignment

	my $FG1out=$orth.".FG1.BUSTED.json";
	my $FG2out=$orth.".FG2.BUSTED.json";
	my $FG3out=$orth.".FG3.BUSTED.json";

	my $cmd ="hyphy busted --alignment $alig --tree $tree --branches FG1 --output $orth/$FG1out";
	$cmd.="; hyphy busted --alignment $alig --tree $tree --branches FG2 --output $orth/$FG2out";
	$cmd.="; hyphy busted --alignment $alig --tree $tree --branches FG3 --output $orth/$FG3out";

#	print "$cmd\n";
	push @cmds, $cmd;
}

my $manager = new Parallel::ForkManager(100);
foreach my $cmd (@cmds) {
	$manager->start and next;
    system($cmd);
    $manager->finish;
}
$manager -> wait_all_children;
```

- **spe_hyphy.tree**
```txt
(((((((Scomberoides_lysan,((Trachinotus_anak{FG2},(Trachinotus_blochii{FG2},Trachinotus_carolinus{FG2})),Trachinotus_ovatus{FG2})),(Remorina_albescens,Echeneis_naucrates)),((((Trachurus_japonicus,Trachurus_trachurus),(Decapterus_maruadsi,Decapterus_macrosoma)),(((Selaroides_leptolepis,(Megalaspis_cordyla,Alepes_djedaba)),((Caranx_sexfasciatus{FG3},(Caranx_heberi{FG3},Caranx_ignobilis{FG3})),Caranx_crysos{FG3})),((Alectis_indicus{FG3},Caranx_equula),Parastromateus_niger{FG3}))),Elagatis_bipinnulata)),((Lates_calcarifer,(Sphyraena_novaehollandiae,Sphyraena_barracuda)),(Toxotes_jaculatrix{FG1},Xiphias_gladius))),Polydactylus_sextarius),(Takifugu_rubripes,Gasterosteus_aculeatus)),Danio_rerio,Oryzias_latipes);
```

- **codeml_parallel.pl**
```perl
use strict;
use warnings;
use Parallel::ForkManager;

# codeml.pl
my $fm   =$ARGV[0]; # input
my $dir  =$ARGV[1]; # make a temporary directory to save the temporary files
my $model=$ARGV[2]; # select the model
my $tree =$ARGV[3]; # the input phylogenetic tree

system("mkdir $dir");

my @cmds;
open FM, $fm or die "can not open $fm\n";
while (<FM>) {
	chomp; my @a=split;
	my $temp=$a[0].".txt";
	open TEMP, ">$dir/$temp" or die "can not create $dir/$temp\n";
	print TEMP "$_\n";
	my $cmd="perl codeml.pl --input $dir/$temp --model $model --dir ./ --tree $tree --icode 0 --omega 1.2";
	push @cmds, $cmd;
}

my $manager = new Parallel::ForkManager(100);
foreach my $cmd (@cmds) {
	$manager->start and next;
	system($cmd);
	$manager->finish;
}
$manager -> wait_all_children;
system("rm -rf $dir/");
```

- **run_free-ratio.sh**
```bash
#!/bin/bash
#SBATCH -J free-ratio # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 100 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.free-ratio # 标准输出文件
#SBATCH -e logs/%j.free-ratio # 错误输出文件
#SBATCH --mail-user=jinglkang@outlook.com

perl codeml_parallel.pl final_orth_input_paml.txt temp1 free-ratio spe.tree
```


