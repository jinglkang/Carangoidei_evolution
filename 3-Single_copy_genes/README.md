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

# 可以参考~/software/paml-4.10.10-linux-x86_64/examples/TipDate.HIV2/baseml.ctl
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
# baseml baseml.ctl
sbatch run_baseml.sh
# Submitted batch job 920110
# Substitution rate is per time unit: 0.133443

# 1.5.2 run mcmctree first time
mcmctree mcmctree_first.ctl # 只是设置了"usedata = 3"

# 1.5.3 run mcmctree second time
# 基于baseml的结果: Substitution rate is per time unit: 0.133443 (m=0.133443, s=0.133443)
# a = (m/s)^2 and b = m/s^2
# a = 1 && b = 7.493836
# 相对应的把mcmctree.ctl改为: rgene_gamma = 1 7.493836; RootAge = '<7'
mcmctree mcmctree_second.ctl
```

- **mcmctree_second.ctl**
```
       seqfile = Single_copy_conca_cds.phy
      treefile = mcmc.tree
      mcmcfile = mcmc_2.txt
       outfile = mcmc_2.txt

         ndata = 1
       seqtype = 0 * 0 : nucleotides; 1:codons; 2:AAs
       usedata = 2 * 0: no data; 1:seq Like; 2:normal approximation; 3:out.BV(in.BV)
         clock = 2 * 1: Global clock; 2: independent rates; 3: correlated rates

         model = 7 * 0 : JC69, 1: K80, 2: F81, 3: F84, 4: HKY85
         alpha = 0    * alpha for gamma rates at sites
         ncatG = 5    * No. categories In discrete gamma
   duplication = 1
     cleandata = 0    * remove sites With ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0.1 M   * birth, death, sampling
   kappa_gamma = 6 2       * gamma prior for kappa
   alpha_gamma = 1 1       * gamma prior for alpha

   rgene_gamma = 1 7.493836    * gammaDir prior for rate for genes
  sigma2_gamma = 1 10 1    * gammaDir prior for sigma^2     (for clock=2 or 3)

      finetune = 1: .1 .1 .1 .1 .1 .1 * auto (0 or 1): times, musigma2, rates, mixing, paras, FossilErr
	   RootAge = '<7'

         print = 1   * 0: no mcmc sample; 1: everything except branch rates 2: everything
        burnin = 20000
      sampfreq = 100
       nsample = 2000
```


- **mcmctree_first.ctl**
```mcmctree_first.ctl
       seqfile = Single_copy_conca_cds.phy
      treefile = mcmc.tree
      mcmcfile = mcmc_1.txt
       outfile = mcmc_1.txt

         ndata = 1
       seqtype = 0 * 0 : nucleotides; 1:codons; 2:AAs
       usedata = 3 * 0: no data; 1:seq Like; 2:normal approximation; 3:out.BV(in.BV)
         clock = 2 * 1: Global clock; 2: independent rates; 3: correlated rates

         model = 7 * 0 : JC69, 1: K80, 2: F81, 3: F84, 4: HKY85
         alpha = 0    * alpha for gamma rates at sites
         ncatG = 5    * No. categories In discrete gamma
   duplication = 1
     cleandata = 0    * remove sites With ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0.1 M   * birth, death, sampling
   kappa_gamma = 6 2       * gamma prior for kappa
   alpha_gamma = 1 1       * gamma prior for alpha

   rgene_gamma = 2 20 1    * gammaDir prior for rate for genes
  sigma2_gamma = 1 10 1    * gammaDir prior for sigma^2     (for clock=2 or 3)

      finetune = 1: .1 .1 .1 .1 .1 .1 * auto (0 or 1): times, musigma2, rates, mixing, paras, FossilErr

         print = 1   * 0: no mcmc sample; 1: everything except branch rates 2: everything
        burnin = 20000
      sampfreq = 100
       nsample = 2000
```

- **run_baseml.sh**
```bash
#!/bin/bash
#SBATCH -J baseml # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 1 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.baseml # 标准输出文件
#SBATCH -e logs/%j.baseml_error # 错误输出文件
#SBATCH --mail-user=jinglkang@outlook.com

baseml baseml.ctl
```

- **baseml.ctl**
```baseml.ctl
      seqfile = Single_copy_conca_cds.phy * sequence data filename
     treefile = baseml.tree * tree Structure file name [Divergence date symbol @, Fossil calibration information is specified using the symbol @; 第一行的两个数字分别代表树有多少个物种和多少棵树, @ 主要用于 BASEML 和 CODEML 的 clock/local clock 模型, B用于MCMCtree]
      outfile = baseml.out           * main result file name

        noisy = 3  * 0,1,2,3,9: how much rubbish on the screen
      verbose = 0  * 0: concise; 1: detailed, 2: too much
      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

        noisy = 3   * 0,1,2,3: how much rubbish on the screen
      verbose = 0   * 1: detailed output, 0: concise output
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI

        model = 7   * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85, 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu * REV = GTR [7是GTR模型]

        ndata = 1   * run ML using subtrees generated from maintree [产生多少个dataset，用于模拟，例如用evolver产生200个重复的数据级，则设置 ndata = 200]
        Mgene = 0   * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff [是否分区，还是认定为好几个区；假如序列是由很多基因连接起来的，则设置 Mgene = 0]

        clock = 1   * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis [clock = 0, 没有时间且不同枝具有不同的rate; clock = 1, 有一个整体的时间且所有枝都是一样的rate; clock = 2依次类推]
    
    fix_kappa = 0   * 0: estimate kappa; 1: fix kappa at value below; 2: kappa for branches [kappa 表示转换（transition）与颠换（transversion）的相对替换率参数]
        kappa = 5   * initial or fixed kappa

    fix_alpha = 0   * 0: estimate alpha; 1: fix alpha at value below [alpha 通常表示伽马分布（Gamma distribution）的形状参数（shape parameter），用于描述不同核苷酸位点或密码子位点之间的进化速率异质性]
        alpha = 0.5   * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * 1: different alpha's for genes, 0: one alpha

        ncatG = 5   * # of categories in the dG, AdG, or nparK models of rates [将 Gamma 分布划分为多少个离散类别]
        nparK = 0   * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK

        nhomo = 0   * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 
 RateAncestor = 0   * (0,1,2): rates (alpha>0) or ancestral states [RateAncestor = 1 强制计算序列每个位点的rate]

   Small_Diff = 1e-8
    cleandata = 1   * remove sites with ambiguity data (1:yes, 0:no)?
*        icode = 0   * (with RateAncestor=1. try "GC" in data,model=4,Mgene=4)
  fix_blength = 0  * 0: ignore, -1: random, 1: initial, 2: fixed, 3: proportional
       method = 0   * Optimization method 0: simultaneous; 1: one branch a time
```

- **baseml.tree**
```mcmc.tree
32 1
(Danio_rerio,(Oryzias_latipes,((Gasterosteus_aculeatus,Takifugu_rubripes),(Polydactylus_sextarius,(((Xiphias_gladius,Toxotes_jaculatrix),(Lates_calcarifer,(Sphyraena_barracuda,Sphyraena_novaehollandiae))),(((Echeneis_naucrates,Remorina_albescens),(Scomberoides_lysan,(Trachinotus_ovatus,(Trachinotus_anak,(Trachinotus_carolinus,Trachinotus_blochii))))'@(.193,.509)')'@(.49,.591)',(Elagatis_bipinnulata,(((Decapterus_macrosoma,Decapterus_maruadsi),(Trachurus_trachurus,Trachurus_japonicus)),((Parastromateus_niger,(Caranx_equula,Alectis_indicus)),((Selaroides_leptolepis,(Alepes_djedaba,Megalaspis_cordyla)),(Caranx_crysos,(Caranx_sexfasciatus,(Caranx_ignobilis,Caranx_heberi))))))))))))'@(.969,1.509)');
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
# (base) shichuang001@login01 Sun Nov 08 2026 16:04:13  
cp ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Phylogeny/Carangoidei_28sp.contree ./
less Carangoidei_28sp.contree|perl -alne 's/\:\d+\.\d+//g;s/\d+//g;s/\[//g;s/\]//g;print' > spe.tree

# 3.1 Run PAML free-ratio
# perl codeml_parallel.pl final_orth_input_paml.txt temp1 free-ratio spe.tree
sbatch run_free-ratio.sh
# Submitted batch job 913416
# (base) shichuang001@login01 Tue Nov 10 2026 23:25:59 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Orthogroups/paml_input
cat free-ratio/*.txt|perl -alne 's/\s+/\t/g;print' > Total_free-ratio_result.txt
# 画图
less Total_free-ratio_result.txt|perl -alne '$i++;print "Orth_id\tspecies\tlength\tbranch\tt\tN\tS\tdN/dS\tdN\tdS\tN*dN\tS*dS" if $i==1;print if $F[1]=~/[a-z]+/ && $i > 1' > Total_free-ratio_result_onlyspe.txt

# 3.2 Run Hyphy busted for positive selection analysis
# hyphy busted --alignment final_alignment.fa --tree spe_hyphy.tree --branches FG1 --output OG000001.FG1.BUSTED.json
# (eggnog) shichuang001@login01 Mon Nov 09 2026 13:43:24 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Orthogroups/paml_input
vi spe_hyphy.tree

# 3.3 逐个谱系进行BUSTED分析
# (base) shichuang001@login01 Mon Nov 09 2026 14:04:28 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Orthogroups/paml_input
# perl run_hyphy.pl final_orth_input_paml.txt spe_hyphy.tree
vi run_hyphy.pl
# 把三类具有扁平性状的鱼类设置为前景枝: FG1, FG2, FG3; 详情见spe_hyphy.tree
sbatch run_hyphy.sh
# Submitted batch job 914530

# 任务非正常结束
# 3.3.1 查看已经运行了并且产生了非空文件的orthID
# (base) shichuang001@login01 Tue Nov 10 2026 15:40:30 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Orthogroups/paml_input
find . -mindepth 2 -type f -name '*.BUSTED.json' -print0 |
while IFS= read -r -d '' f; do
    if [ -s "$f" ]; then
        echo "NON-EMPTY  $f"
    else
        echo "EMPTY      $f"
    fi
done | \
perl -alne 'if ($F[0] eq "NON-EMPTY"){(my $orth)=$F[1]=~/\.\/(.*?)\//;$hash{$orth}++};END{foreach $key (sort keys %hash){print "$key\t$hash{$key}" if $hash{$key}==3}}' > Finish_hyphy_orthID.txt
# 1463 orth finished

# 3.3.2 除去已经运行的orthID，从而整理出未运行hyphy的orthID
# (base) shichuang001@login01 Tue Nov 10 2026 15:47:30 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Orthogroups/paml_input
perl Extract_unfinish_orth.pl > unfinish_orth.txt

# 3.3.3 把run_hyphy.sh中的"final_orth_input_paml.txt"替换成"unfinish_orth.txt"
# (base) shichuang001@login01 Tue Nov 10 2026 15:51:04 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Orthogroups/paml_input
sbatch run_hyphy.sh
# Submitted batch job 916635

# 3.3.4 整理hyphy busted结果
# FG1
# (base) shichuang001@login01 Wed Nov 11 2026 21:49:19 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Orthogroups/paml_input
perl extract_hyphy_busted.pl FG1.BUSTED.json FG1 > hyphy_FG1_BUSTED_summary.tsv

# 3.3.5 BH-FDR对p-value进行校正
# (base) shichuang001@login01 Wed Nov 11 2026 22:36:12 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Orthogroups/paml_input
perl Calculate_FDR.pl hyphy_FG1_BUSTED_summary.tsv yes 6 > hyphy_FG1_BUSTED_summary_FDR.tsv

# 3.3.6 挑选出正选择基因
less hyphy_FG1_BUSTED_summary_FDR.tsv|perl -alne 'if (/^Gene/){print}else{print if $F[-1]<=0.05}' > hyphy_FG1_PSGs.tsv

# 3.3.7 把注释信息添加进去
perl anno.pl hyphy_FG1_PSGs.tsv > hyphy_FG1_PSGs_anno.tsv # 5 PSGs

# FG2
perl extract_hyphy_busted.pl FG2.BUSTED.json FG2 > hyphy_FG2_BUSTED_summary.tsv
perl Calculate_FDR.pl hyphy_FG2_BUSTED_summary.tsv yes 6 > hyphy_FG2_BUSTED_summary_FDR.tsv
less hyphy_FG2_BUSTED_summary_FDR.tsv|perl -alne 'if (/^Gene/){print}else{print if $F[-1]<=0.05}' > hyphy_FG2_PSGs.tsv
perl anno.pl hyphy_FG2_PSGs.tsv > hyphy_FG2_PSGs_anno.tsv # 66 PSGs

# FG3
perl extract_hyphy_busted.pl FG3.BUSTED.json FG3 > hyphy_FG3_BUSTED_summary.tsv
perl Calculate_FDR.pl hyphy_FG3_BUSTED_summary.tsv yes 6 > hyphy_FG3_BUSTED_summary_FDR.tsv
less hyphy_FG3_BUSTED_summary_FDR.tsv|perl -alne 'if (/^Gene/){print}else{print if $F[-1]<=0.05}' > hyphy_FG3_PSGs.tsv
perl anno.pl hyphy_FG3_PSGs.tsv > hyphy_FG3_PSGs_anno.tsv # 394 PSGs


# 3.4 逐个谱系进行aBSREL分析
# (base) shichuang001@login01 Tue Nov 10 2026 16:35:09 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Orthogroups/paml_input
# perl run_hyphy_update.pl final_orth_input_paml.txt spe_hyphy.tree aBSREL
vi run_hyphy_update.pl; vi run_hyphy_aBSREL.sh
sbatch run_hyphy_aBSREL.sh
# Submitted batch job 916653
# 终止全部基因的aBSREL分析，着重分析正选择基因的aBSREL分析
# 已经完成的基因
find . -mindepth 2 -type f -name '*.aBSREL.json' -print0 |
while IFS= read -r -d '' f; do
    if [ -s "$f" ]; then
        echo "NON-EMPTY  $f"
    else
        echo "EMPTY      $f"
    fi
done | \
perl -alne 'if ($F[0] eq "NON-EMPTY"){(my $orth)=$F[1]=~/\.\/(.*?)\//;$hash{$orth}++};END{foreach $key (sort keys %hash){print "$key\t$hash{$key}" if $hash{$key}==3}}' > Finish_hyphy_aBSREL_orthID.txt

cat *PSGs_anno.tsv|perl -alne 'next if /^Gene/;print $F[0]' > PSGs_id.txt

# 3.4.1 查看哪些PSGs没有运行aBSREL
perl left_aBSREL.pl > left_aBSREL_PSGs.tsv

# 3.4.2 对left_aBSREL_PSGs.tsv里面的基因运行aBSREL
sbatch run_hyphy_aBSREL.sh
# Submitted batch job 919595

# 3.4.3 aBSREL的json文件出错误，可能是因为上次结束任务，导致某些基因结果没跑完，把所有之前认为跑完的PSGs再重新跑一下
# last_aBSREL_PSGs.tsv:剩下123个
sbatch run_hyphy_aBSREL.sh
# Submitted batch job 925392


# 3.4.3 FG1, FG2, FG3
# (base) shichuang001@login01 Thu Nov 12 2026 16:20:04 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Orthogroups/paml_input
perl extract_aBSREL.pl OG*/*FG1.aBSREL.json > FG1_aBSREL_summary.tsv
perl Extract_PSGs_aBSREL.pl hyphy_FG1_PSGs_id.tsv FG1_aBSREL_summary.tsv > FG1_PSGs_aBSREL.tsv

perl extract_aBSREL.pl OG*/*FG2.aBSREL.json > FG2_aBSREL_summary.tsv
perl Extract_PSGs_aBSREL.pl hyphy_FG2_PSGs_id.tsv FG2_aBSREL_summary.tsv > FG2_PSGs_aBSREL.tsv # 60个，少了6个

perl extract_aBSREL.pl OG*/*FG3.aBSREL.json > FG3_aBSREL_summary.tsv
perl Extract_PSGs_aBSREL.pl hyphy_FG3_PSGs_id.tsv FG3_aBSREL_summary.tsv > FG3_PSGs_aBSREL.tsv # 392个，少了2个


# 3.5 MEME分析正选择位点
# (base) shichuang001@login01 Tue Nov 10 2026 16:40:19 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Orthogroups/paml_input
vi run_hyphy_meme.pl; vi run_hyphy_MEME.sh
sbatch run_hyphy_MEME.sh
# Submitted batch job 916658

# 3.5.1 提取MEME结果(只包括p-value <= 0.05的基因)
# (base) shichuang001@login01 Thu Nov 12 2026 00:07:45 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Orthogroups/paml_input
perl extract_MEME.pl OG*/*.MEME.json # 结果文件为MEME_candidates.tsv

# 3.5.2 提取出PSGs的MEME结果(只包括p-value <= 0.05的基因)
perl Extract_PSGs_MEME.pl PSGs_id.txt MEME_candidates.tsv > PSGs_MEME.tsv

# 3.5.3 再运行一下MEME，设定前景分支，只运行正选择基因
# perl run_hyphy_update.pl PSGs_id.txt spe_hyphy.tree MEME
# (base) shichuang001@login01 Thu Nov 12 2026 01:02:21 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Orthogroups/paml_input
sbatch run_hyphy_MEME_PSGs.sh
# Submitted batch job 919615
# 修改extract_MEME.pl为extract_MEME_update.pl: 自行定义输出文件名字

# 3.5.4 FG1
# (base) shichuang001@login01 Thu Nov 12 2026 15:12:37 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Orthogroups/paml_input
perl extract_MEME_update.pl OG*/*.FG1.MEME.json > FG1_MEME_summary.tsv
less hyphy_FG1_PSGs.tsv|perl -alne 'next if /^Gene/;print $F[0]' > hyphy_FG1_PSGs_id.tsv
perl Extract_PSGs_MEME.pl hyphy_FG1_PSGs_id.tsv FG1_MEME_summary.tsv > FG1_PSGs_MEME.tsv # 5个

# 3.5.5 FG2
perl extract_MEME_update.pl OG*/*.FG2.MEME.json > FG2_MEME_summary.tsv
less hyphy_FG2_PSGs.tsv|perl -alne 'next if /^Gene/;print $F[0]' > hyphy_FG2_PSGs_id.tsv
perl Extract_PSGs_MEME.pl hyphy_FG2_PSGs_id.tsv FG2_MEME_summary.tsv > FG2_PSGs_MEME.tsv # 64个，少了2个

# 3.5.6 FG3
perl extract_MEME_update.pl OG*/*.FG3.MEME.json > FG3_MEME_summary.tsv
less hyphy_FG3_PSGs.tsv|perl -alne 'next if /^Gene/;print $F[0]' > hyphy_FG3_PSGs_id.tsv
perl Extract_PSGs_MEME.pl hyphy_FG3_PSGs_id.tsv FG3_MEME_summary.tsv > FG3_PSGs_MEME.tsv # 380个，少了14个
```

- **extract_aBSREL.pl**
```perl
#!/usr/bin/perl
use strict;
use warnings;
use JSON::PP;
use File::Basename;

print join("\t",
    "Gene", "Branch", "LRT", "P_uncorrected", "P_corrected",
    "Baseline_omega", "Rate_classes", "Omega_classes",
    "N_sequences", "N_sites"
), "\n";

foreach my $file (@ARGV) {

    open(my $IN, "<", $file) or die "Cannot open $file: $!\n";
    local $/;
    my $json_text = <$IN>;
    close($IN);

    # 检查文件是否为空
	unless (defined $json_text && $json_text =~ /\S/) {
    	warn "WARNING: Empty JSON file, skipped: $file\n";
    	next;
	}

	# 检查JSON格式是否正确
	my $data = eval { decode_json($json_text) };

	if ($@) {
    	warn "WARNING: Invalid JSON, skipped: $file\n";
    	warn "Reason: $@\n";
    	next;
	}

    (my $gene) = basename($file) =~ /(OG\d+?)\./;

    my $tested = $data->{"tested"}{"0"};
    my $attrs  = $data->{"branch attributes"}{"0"};
    my $nseq   = $data->{"input"}{"number of sequences"};
    my $nsites = $data->{"input"}{"number of sites"};

    foreach my $branch (sort keys %$tested) {

        # 只提取被检验的分支
        next unless $tested->{$branch} eq "test";

        my $r = $attrs->{$branch};
        next unless defined $r;

        my $lrt = defined $r->{"LRT"} ? $r->{"LRT"} : "NA";
        my $p1  = defined $r->{"Uncorrected P-value"}
                    ? $r->{"Uncorrected P-value"} : "NA";
        my $p2  = defined $r->{"Corrected P-value"}
                    ? $r->{"Corrected P-value"} : "NA";

        my $omega = defined $r->{"Baseline MG94xREV omega ratio"}
                    ? $r->{"Baseline MG94xREV omega ratio"} : "NA";

        my $classes = defined $r->{"Rate classes"}
                    ? $r->{"Rate classes"} : "NA";

        # 提取各个omega类别及其对应比例
        my @omega_info;
        if (defined $r->{"Rate Distributions"}) {
            foreach my $class (@{$r->{"Rate Distributions"}}) {
                my ($w, $prop) = @$class;
                push @omega_info, "omega=$w,prop=$prop";
            }
        }
        my $omega_str = @omega_info ? join(";", @omega_info) : "NA";

        print join("\t",
            $gene, $branch, $lrt, $p1, $p2,
            $omega, $classes, $omega_str,
            $nseq, $nsites
        ), "\n";
    }
}
```

- **extract_MEME_update.pl**
```perl
#!/usr/bin/perl
use strict;
use warnings;
use JSON::PP;
use File::Basename;

# Usage:
# perl extract_MEME.pl OG0010263.MEME.json
# perl extract_MEME.pl *.MEME.json


print join("\t",
    qw(Gene Site Alpha Beta_minus P_minus Beta_plus P_plus LRT P_value FDR Branches)
), "\n";

my @results;

my $k;
for my $file (@ARGV) {
    (my $gene) = basename($file) =~ /(OG\d+?)\./;

    # Read JSON
    open my $fh, "<", $file or die "Cannot open $file: $!\n";
    local $/;
    my $json = <$fh>;
    close $fh;

    my $data = decode_json($json);

    # Extract MEME site results
    my $content = $data->{MLE}->{content}->{"0"};

    die "Cannot find MLE content in $file\n"
        unless defined $content && ref($content) eq "ARRAY";

    my @gene_results;

    for my $i (0 .. $#$content) {

        my $row = $content->[$i];

        # MLE.headers column order:
        # 0 alpha
        # 1 beta-
        # 2 p-
        # 3 beta+
        # 4 p+
        # 5 LRT
        # 6 p-value
        # 7 # branches under selection

        my $alpha   = $row->[0];
        my $beta_m  = $row->[1];
        my $p_m     = $row->[2];
        my $beta_p  = $row->[3];
        my $p_p     = $row->[4];
        my $lrt     = $row->[5];
        my $pvalue  = $row->[6];
        my $branches = $row->[7];

        next unless defined $pvalue && $pvalue =~ /^\d*\.?\d+(?:[eE][+-]?\d+)?$/;

        push @gene_results, {
            gene     => $gene,
            site     => $i + 1,
            alpha    => $alpha,
            beta_m   => $beta_m,
            p_m      => $p_m,
            beta_p   => $beta_p,
            p_p      => $p_p,
            lrt      => $lrt,
            pvalue   => $pvalue,
            branches => $branches,
        };
    }

    # BH-FDR correction within this gene
    my @sorted = sort {
        $a->{pvalue} <=> $b->{pvalue}
    } @gene_results;

    my $m = scalar @sorted;
    my $min_q = 1;

    for (my $i = $m - 1; $i >= 0; $i--) {

        my $rank = $i + 1;
        my $q = $sorted[$i]->{pvalue} * $m / $rank;

        $q = 1 if $q > 1;
        $min_q = $q if $q < $min_q;

        $sorted[$i]->{fdr} = $min_q;
    }

    # Output sites with raw P < 0.05
    for my $r (@gene_results) {

        next unless $r->{pvalue} < 0.05;

        printf "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%.8g\t%.8g\t%s\n",
            $r->{gene},
            $r->{site},
            $r->{alpha},
            $r->{beta_m},
            $r->{p_m},
            $r->{beta_p},
            $r->{p_p},
            $r->{lrt},
            $r->{pvalue},
            $r->{fdr},
            $r->{branches};
    }
}
```

- **Extract_PSGs_MEME.pl**
```perl
#!/usr/bin/perl
use strict;
use warnings;

my %hash;
my $psg =$ARGV[0];
open PSG, $psg or die "can not open $psg:$!\n";
while (<PSG>) {
	chomp;
	$hash{$_}++;
}

my $meme=$ARGV[1];
open MEME, $meme or die "can not open $meme\n";
while (<MEME>) {
	chomp;
	my @a=split;
	if (/^Gene/) {
		print "$_\n";
	} elsif ($hash{$a[0]} && $a[-3] <= 0.05) {
		print "$_\n";
	}
}
```

- **extract_MEME.pl**
```perl
#!/usr/bin/perl
use strict;
use warnings;
use JSON::PP;
use File::Basename;

# Usage:
# perl extract_MEME.pl OG0010263.MEME.json
# perl extract_MEME.pl *.MEME.json

my $outfile = "MEME_candidates.tsv";

open my $out, ">", $outfile or die "Cannot write $outfile: $!\n";

print $out join("\t",
    qw(Gene Site Alpha Beta_minus P_minus Beta_plus P_plus LRT P_value FDR Branches)
), "\n";

my @results;

for my $file (@ARGV) {

    my $gene = basename($file);
    $gene =~ s/\.MEME\.json$//;

    # Read JSON
    open my $fh, "<", $file or die "Cannot open $file: $!\n";
    local $/;
    my $json = <$fh>;
    close $fh;

    my $data = decode_json($json);

    # Extract MEME site results
    my $content = $data->{MLE}->{content}->{"0"};

    die "Cannot find MLE content in $file\n"
        unless defined $content && ref($content) eq "ARRAY";

    my @gene_results;

    for my $i (0 .. $#$content) {

        my $row = $content->[$i];

        # MLE.headers column order:
        # 0 alpha
        # 1 beta-
        # 2 p-
        # 3 beta+
        # 4 p+
        # 5 LRT
        # 6 p-value
        # 7 # branches under selection

        my $alpha   = $row->[0];
        my $beta_m  = $row->[1];
        my $p_m     = $row->[2];
        my $beta_p  = $row->[3];
        my $p_p     = $row->[4];
        my $lrt     = $row->[5];
        my $pvalue  = $row->[6];
        my $branches = $row->[7];

        next unless defined $pvalue && $pvalue =~ /^\d*\.?\d+(?:[eE][+-]?\d+)?$/;

        push @gene_results, {
            gene     => $gene,
            site     => $i + 1,
            alpha    => $alpha,
            beta_m   => $beta_m,
            p_m      => $p_m,
            beta_p   => $beta_p,
            p_p      => $p_p,
            lrt      => $lrt,
            pvalue   => $pvalue,
            branches => $branches,
        };
    }

    # BH-FDR correction within this gene
    my @sorted = sort {
        $a->{pvalue} <=> $b->{pvalue}
    } @gene_results;

    my $m = scalar @sorted;
    my $min_q = 1;

    for (my $i = $m - 1; $i >= 0; $i--) {

        my $rank = $i + 1;
        my $q = $sorted[$i]->{pvalue} * $m / $rank;

        $q = 1 if $q > 1;
        $min_q = $q if $q < $min_q;

        $sorted[$i]->{fdr} = $min_q;
    }

    # Output sites with raw P < 0.05
    for my $r (@gene_results) {

        next unless $r->{pvalue} < 0.05;

        printf $out "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%.8g\t%.8g\t%s\n",
            $r->{gene},
            $r->{site},
            $r->{alpha},
            $r->{beta_m},
            $r->{p_m},
            $r->{beta_p},
            $r->{p_p},
            $r->{lrt},
            $r->{pvalue},
            $r->{fdr},
            $r->{branches};
    }
}

close $out;

print "Finished! Results saved to $outfile\n";
```

- **extract_hyphy_busted.pl**
```perl
#!/usr/bin/perl
use strict;
use warnings;
use JSON::PP;

my $list  ="final_orth_input_paml.txt";
my $suffix=$ARGV[0];
my $fg    =$ARGV[1];
print "Gene\tForeground\tLRT\tω3\tω3_proportion\tP-value\n";

open LIST, $list or die "can not open $list:$!\n";
while (<LIST>) {
	chomp;
	my $file="$_/$_.$suffix";
	open FILE, $file or die "can not open $file:$!\n";
	local $/; # 不再寻找分隔符，直接读取文件句柄中剩余的所有内容，直到文件末尾（EOF）
	my $json=<FILE>;
	close FILE;
	my $data = decode_json($json);
	
	# LRT
	my $lrt = $data->{"test results"}->{LRT};

	# P-value
	my $pvalue = $data->{"test results"}->{"p-value"};

	# Unconstrained model
	my $model = $data->{fits}->{"Unconstrained model"};

	# Rate distributions
	my $rate = $model->{"Rate Distributions"};

	# Foreground: Test
	my $test = $rate->{Test};

	# omega0
	my $omega0 = $test->{0}->{omega};

	# omega1
	my $omega1 = $test->{1}->{omega};

	# omega3
	my $omega3 = $test->{2}->{omega};

	# omega3 proportion
	my $proportion = $test->{2}->{proportion};

	print "$_\t$fg\t$lrt\t$omega3\t$proportion\t$pvalue\n";
}
```

- **run_hyphy_MEME.sh**
```bash
#!/bin/bash
#SBATCH -J hyphy_MEME # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 50 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.hyphy_MEME # 标准输出文件
#SBATCH -e logs/%j.hyphy_MEME # 错误输出文件
#SBATCH --mail-user=jinglkang@outlook.com

perl run_hyphy_meme.pl final_orth_input_paml.txt spe_hyphy.tree
```

- **run_hyphy_meme.pl**
```perl
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;

# MEME
# hyphy meme --alignment final_alignment.fa --tree spe_hyphy.tree --output OG000001.MEME.json

my $list=$ARGV[0]; # The list
my $tree=$ARGV[1]; # The tree

my @cmds;
open LIST, $list or die "can not open $list\n";
while (<LIST>) {
	chomp;
	my @a=split;
	my $orth=$a[0];
	my $alig="$orth/final_alignment.fa"; # the alignment

	my $out=$orth.".MEME.json";
	my $cmd ="hyphy meme --alignment $alig --tree $tree --output $orth/$out";

#	print "$cmd\n";
	push @cmds, $cmd;
}

my $manager = new Parallel::ForkManager(50);
foreach my $cmd (@cmds) {
	$manager->start and next;
    system($cmd);
    $manager->finish;
}
$manager -> wait_all_children;
```

- **run_hyphy_aBSREL.sh**
```bash
#!/bin/bash
#SBATCH -J hyphy_aBSREL # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 50 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.hyphy_aBSREL # 标准输出文件
#SBATCH -e logs/%j.hyphy_aBSREL # 错误输出文件
#SBATCH --mail-user=jinglkang@outlook.com

perl run_hyphy_update.pl final_orth_input_paml.txt spe_hyphy.tree aBSREL
```

- **run_hyphy_update.pl**
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

# aBSREL

# hyphy absrel --alignment final_alignment.fa --tree spe_hyphy.tree --branches FG1 --output OG000001.FG1.aBSREL.json
# hyphy absrel --alignment final_alignment.fa --tree spe_hyphy.tree --branches FG2 --output OG000001.FG2.aBSREL.json
# hyphy absrel --alignment final_alignment.fa --tree spe_hyphy.tree --branches FG3 --output OG000001.FG3.aBSREL.json

my $list=$ARGV[0]; # The list
my $tree=$ARGV[1]; # The tree
my $type=$ARGV[2]; # run BUSTED or aBSREL

my @cmds;
open LIST, $list or die "can not open $list\n";
while (<LIST>) {
	chomp;
	my @a=split;
	my $orth=$a[0];
	my $alig="$orth/final_alignment.fa"; # the alignment

	my $FG1out=$orth.".FG1.$type.json";
	my $FG2out=$orth.".FG2.$type.json";
	my $FG3out=$orth.".FG3.$type.json";

	my $typeNew=lc($type);

	my $cmd ="hyphy $typeNew --alignment $alig --tree $tree --branches FG1 --output $orth/$FG1out";
	$cmd.="; hyphy $typeNew --alignment $alig --tree $tree --branches FG2 --output $orth/$FG2out";
	$cmd.="; hyphy $typeNew --alignment $alig --tree $tree --branches FG3 --output $orth/$FG3out";

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

- **Extract_unfinish_orth.pl**
```perl
#!/usr/bin/perl
use strict;
use warnings;

my %hash;
my $name="Finish_hyphy_orthID.txt";
open NAME, $name or die "can not open $name\n";
while (<NAME>) {
	chomp;
	my @a=split /\t/;
	$hash{$a[0]}++;
}

my $Torth="final_orth_input_paml.txt";
open TORTH, $Torth or die "can not open $Torth\n";
while (<TORTH>) {
	chomp;
	print "$_\n" unless $hash{$_};	
}
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


