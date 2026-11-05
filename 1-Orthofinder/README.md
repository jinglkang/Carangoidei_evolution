# obtain the orthogroups for evolutionary analysis
## 1. orthofinder
### 1.1 obtain the longest pep gff3 for each species
```bash
# 1. keep_longest_isoform: keep the one with the longest CDS
# (agat) shichuang001@login01 Sun Nov 01 2026 23:53:53 ~/jlkang/Carangoidei_evolution/rename
conda activate agat
# agat_sp_keep_longest_isoform.pl --gff Alectis_indicus.gff3 -o longest_gff/Alectis_indicus.longest.gff3

sbatch keep_longest_isoform.sh
# Submitted batch job 909345
```

- run_agat_sp_keep_longest_isoform.pl
```run_agat_sp_keep_longest_isoform.pl
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
```

- keep_longest_isoform.sh
```keep_longest_isoform.sh
#!/bin/bash
#SBATCH -J agat # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 32 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.agat # 标准输出文件
#SBATCH -e logs/%j.agat # 错误输出文件

perl run_agat_sp_keep_longest_isoform.pl
```

### 1.2 obtain the sequences of nuc and pep
```bash
# (base) shichuang001@login01 Mon Nov 02 2026 00:48:52 ~/jlkang/Carangoidei_evolution/rename
mkdir longest_nuc; mkdir longest_pep
# gffread -g $genome -x longest_nuc/$spe.cds -y longest_pep/$spe.pep longest_gff/$new
sbatch run_obtain_longest_nuc_pep.sh
# Submitted batch job 909347
```

- obtain_longest_nuc_pep.pl
```obtain_longest_nuc_pep.pl
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
```

```run_obtain_longest_nuc_pep.sh
#!/bin/bash
#SBATCH -J longest_nuc_pep # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 32 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.longest_nuc_pep # 标准输出文件
#SBATCH -e logs/%j.longest_nuc_pep # 错误输出文件

perl obtain_longest_nuc_pep.pl
```

**把判定80%位点为ATCG的序列去除掉，把序列中的"."去掉**
```bash
# (base) shichuang001@login01 Mon Nov 02 2026 01:06:57 ~/jlkang/Carangoidei_evolution/rename/longest_pep
perl remove_suspicious_proteins.pl

# 改后缀
# (base) shichuang001@login01 Mon Nov 02 2026 01:11:34 ~/jlkang/Carangoidei_evolution/rename/longest_pep/longest_pep_2
for file in *.fa; do mv "$file" "${file%.fa}.fas"; done

mkdir Revise; perl temp1.pl

# (orthofinder_env) shichuang001@login01 Mon Nov 02 2026 01:28:12 ~/jlkang/Carangoidei_evolution/rename/longest_pep/longest_pep_2/Revise
sbatch orthofinder.job
# Submitted batch job 909352
```

### 1.2 根据结果所得的单拷贝基因进行系统发育树构建
**共有3689个单拷贝基因**
```bash
# 1. 先在每个文件的序列header加上相对应物种的名称，以便于后续序列根据物种进行连接
# (base) shichuang001@login01 Mon Nov 02 2026 14:46:17 ~/jlkang/Carangoidei_evolution/rename/longest_pep/longest_pep_2/Revise/OrthoFinder/Results_Nov02
mkdir Rename_Single_Copy_Orthologue_Sequences; perl rename_single_copy_fa.pl

# 2. muscle先比对，然后根据相似性进行剪切
# (base) shichuang001@login01 Mon Nov 02 2026 15:08:51 ~/jlkang/Carangoidei_evolution/rename/longest_pep/longest_pep_2/Revise/OrthoFinder/Results_Nov02/Rename_Single_Copy_Orthologue_Sequences
mkdir muscle; mkdir trimal
# perl Align_trim.pl
# (base) shichuang001@login01 Mon Nov 02 2026 16:27:47 ~/jlkang/Carangoidei_evolution/rename/longest_pep/longest_pep_2/Revise/OrthoFinder/Results_Nov02/Rename_Single_Copy_Orthologue_Sequences
sbatch Align_trim.sh
# Submitted batch job 909528

# 3. 每个fasta最终的trimal文件中，只要存在某个蛋白序列长度小于50，则去除；然后按照物种依次连接剩下trimal文件的相对应的序列
# (base) shichuang001@login01 Mon Nov 02 2026 16:53:01 ~/jlkang/Carangoidei_evolution/rename/longest_pep/longest_pep_2/Revise/OrthoFinder/Results_Nov02/Rename_Single_Copy_Orthologue_Sequences/trimal
perl Conca_seq.pl > ../../Phylogeny/Single_copy_conca.fa
# 共2,071,585个蛋白

# 4. Single_copy_conca.fa，用iqtree进行构树
# 以Oryzias_latipes,Takifugu_rubripes,Danio_rerio,Gasterosteus_aculeatus为外类群
# iqtree -s Single_copy_conca.fa -o Oryzias_latipes,Takifugu_rubripes,Danio_rerio,Gasterosteus_aculeatus --seqtype AA -T 32 -B 1000 -m MFP
# (base) shichuang001@login01 Mon Nov 02 2026 17:00:50 ~/jlkang/Carangoidei_evolution/rename/longest_pep/longest_pep_2/Revise/OrthoFinder/Results_Nov02/Phylogeny
sbatch iqtree.sh
# Submitted batch job 909539

# 以Danio_rerio为外类群
# iqtree3 -s Single_copy_conca.fa -o Danio_rerio -st AA -m MFP -B 1000 -T 50 -pre Carangoidei_28sp
sbatch iqtree2.sh
# Submitted batch job 909610
``` 

- rename_single_copy_fa.pl
```perl
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

my $nb="Orthogroups/Orthogroups.GeneCount.tsv";
my @headers;
open NB, $nb or die "can not open $nb\n";
while (<NB>) {
	chomp;
	my @a=split;
	if (/^Orthogroup/) {
		@headers=@a;
	}
}

my @fas=<Single_Copy_Orthologue_Sequences/*.fa>;
foreach my $fa (@fas) {
	my %hash;
	my $i=0;
	my $name=basename($fa);
	open FA, $fa or die "can not open $fa\n";
	my $newfa="Rename_Single_Copy_Orthologue_Sequences/$name";
	open NEW, ">$newfa" or die "can not create $newfa\n";
	while (<FA>) {
		chomp;
		if (/^>/) {
			$i++;
			#print NEW ">$headers[$i]\n";
		} else {
			#print NEW "$_\n";
			$hash{$headers[$i]}.=$_;
		}
	}
	for (my $j = 1; $j < @headers-1; $j++) {
		print NEW ">$headers[$j]\n$hash{$headers[$j]}\n";
	}
}
```

- Align_trim.pl
```perl
#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

my @fasta=<*.fa>;
my @cmds;
foreach my $fasta (@fasta) {
    my $align=$fasta.".align";
    my $cmd1="muscle -align $fasta -output muscle/$align";
    #system("mv $align $fasta");
    my $new=$fasta.".trimal";
    $cmd1.=";trimal -in muscle/$align -out trimal/$new -gt 0.8 -st 0.001 -cons 60";
    push @cmds, $cmd1;
}

my $manager = new Parallel::ForkManager(100);
foreach my $cmd (@cmds) {
        $manager->start and next;
        system($cmd);
        $manager->finish;
}
$manager -> wait_all_children;
```

- Align_trim.sh
```Align_trim.sh
#!/bin/bash
#SBATCH -J align_trim # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 100 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.align_trim # 标准输出文件
#SBATCH -e logs/%j.align_trim # 错误输出文件

perl Align_trim.pl
```

- Conca_seq.pl
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
			$spe=$_;
			$hash1{$spe}++;
			push @spes, $spe if $hash1{$spe}==1;
		} else {
			$hash2{$spe}.=$_;
		}
	}

	my $len=length($hash2{$spe});
	if ($len >= 50) {
		open FA, $fa or die "can not open $fa\n";
		while (<FA>) {
			chomp;
			if (/^>/) {
				s/>//;
				$spe=$_;
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

- iqtree.sh
```bash
#!/bin/bash
#SBATCH -J iqtree # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 32 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.iqtree # 标准输出文件
#SBATCH -e logs/%j.iqtree # 错误输出文件
#SBATCH --mail-user=jinglkang@outlook.com

iqtree -s Single_copy_conca.fa -o Oryzias_latipes,Takifugu_rubripes,Danio_rerio,Gasterosteus_aculeatus --seqtype AA -T 32 -B 1000 -m MFP
```

- iqtree2.sh
```bash
#!/bin/bash
#SBATCH -J iqtree # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 50 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.iqtree # 标准输出文件
#SBATCH -e logs/%j.iqtree # 错误输出文件
#SBATCH --mail-user=jinglkang@outlook.com

iqtree3 -s Single_copy_conca.fa -o Danio_rerio -st AA -m MFP -B 1000 -T 50 -pre Carangoidei_28sp
```

# 重新跑Orthofinder
```bash
# (base) shichuang001@login01 Wed Nov 04 2026 05:18:01 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Single_Copy_Orthologue_Sequences
ll *.fa|wc -l # 共3472个单拷贝同源基因
mkdir muscle; mkdir trimal
sbatch Align_trim.sh
# (base) shichuang001@login01 Wed Nov 04 2026 05:32:24 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Single_Copy_Orthologue_Sequences/trimal
perl Conca_seq.pl > ../../Phylogeny/Single_copy_conca.fa
# 共1,989,887个AA
# (base) shichuang001@login01 Wed Nov 04 2026 05:51:05 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Phylogeny
sbatch iqtree2.sh
# Submitted batch job 909917
```

- Conca_seq.pl
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
    if ($len >= 50) {
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
