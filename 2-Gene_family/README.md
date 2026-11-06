# Gene family evolution
## 1. Remove gene families with ≥ 100 gene copies
```bash
# (base) shichuang001@login01 Mon Nov 02 2026 23:09:43 ~/jlkang/Carangoidei_evolution/rename/longest_pep/longest_pep_2/Revise/OrthoFinder/Results_Nov02/Orthogroups
less Orthogroups.GeneCount.tsv|perl -alne 'my $info;for (my $i = 0; $i < @F-1; $i++){$info.=$F[$i]."\t"};$info=~s/\s+$//;if(/^Orthogroup/){print "Desc\t$info"}else{print "null\t$info"}' > Carangoidei_genefamilies.tab
python clade_and_size_filter.py -i Carangoidei_genefamilies.tab -o filtered_cafe_input.txt -s # 34599 -> 33249
```

## 2. Making the species tree ultrametric
```bash
# 1. 以斑马鱼和三刺鱼的分化时间to scale branch lengths into time units
# (base) shichuang001@login01 Tue Nov 03 2026 20:57:03 ~/jlkang/Carangoidei_evolution/rename/longest_pep/longest_pep_2/Revise/OrthoFinder/Results_Nov02/Phylogeny
python prep_r8s.py -i Carangoidei_28sp.contree -o r8s_ctl_file.txt -s 638853 -p 'Danio_rerio,Gasterosteus_aculeatus' -c '156.53'
r8s -b -f r8s_ctl_file.txt > r8s_tmp.txt
tail -n 1 r8s_tmp.txt | cut -c 16- > Carangoidei28spe_tree.txt
cafe5 -i filtered_cafe_input.txt -t Carangoidei28spe_tree.txt
```

# 重新开始，因为之前的三刺鱼基因组和注释文件拷贝错误
## CAFE5结果文件解读
**model_change.txt** - A tab-separated file listing, for each family and clade, the difference between it and its parent clade in the reconstruction that was performed.
**model_count.txt** - A tab-separated file listing, for each family and clade, the reconstructed value in that clade.
**model_clade_results.txt** - It consists of a header line, “Taxon Increase Decrease”, and for each node in the tree, a tab-separated count of the number of families which have increased and decreased for that node.
**model_asr.tre** -  the reconstructed states of the families in the Nexus file format, "\*" indicate significant or not.
**每个基因家族在特定节点是否显著扩张和收缩将基于三个文件**：1. model_asr.tre; 2. model_change.txt; 3. model_count.txt.

```bash
# 1. 筛选出显著扩张或者收缩的树
echo $'#nexus\nbegin trees;'>Significant_trees.tre
grep "*" Base_asr.tre >>Significant_trees.tre
echo "end;">>Significant_trees.tre

# 2. 基于Base_count.tab, Significant_trees.tre, Base_change.tab三个文件整合信息
perl extract_sig_fm.pl > Sig_family.txt
```

- **extract_sig_fm.pl**
```bash
#!/usr/bin/perl
use strict;
use warnings;

my $fanb  ="Base_count.tab"; # extract the gene copy number
my $sigtre="Significant_trees.tre"; # extract the significant infor for specific node
my $fach  ="Base_change.tab"; # extract the info of expasion or contraction

# Final output: gene_copy_nb(change_nb;significant_info)
my (%sig, %nbc);
open SIG, $sigtre or die "can not open $sigtre\n";
while (<SIG>) {
	chomp;
	s/^\s+//;
	my @a=split;
	while (/\*(?:_[^:,()]+)?/g) {
		my $star_pos = $-[0];
		my $before = substr($_, 0, $star_pos);
		if ($before =~ /<(\d+)>(?!.*<\d+>)/s) {
     		#print "$1\n";
     		my $id=$1;
     		$sig{$a[1]}->{$id}++;
    	}
	}
}

my @headers;
open FACH, $fach or die "can not open $fach\n";
while (<FACH>) {
	chomp;
	my @a=split;
	if (/^FamilyID/) {
		@headers=@a;
	} else {
		for (my $i = 1; $i < @a; $i++) {
			(my $id)=$headers[$i]=~/<(\d+?)>/;
			$nbc{$a[0]}->{$id}=$a[$i];
		}
	}
}

open FANB, $fanb or die "can not open $fanb\n";
while (<FANB>) {
	chomp;
	my @a=split;
	if (/^FamilyID/) {
		print "$_\n";
	} elsif	($sig{$a[0]}) {
		my $info;
		for (my $i = 1; $i < @a; $i++) {
			(my $id)=$headers[$i]=~/<(\d+?)>/;
			my ($change, $signif);
			$change=$nbc{$a[0]}->{$id};
			my $new;
			if ($sig{$a[0]}->{$id}) {
				$signif="sig";
				$new=$a[$i]."($change;$signif)";
			} else {
				$new=$a[$i]."($change)";
			}
			$info.=$new."\t";
		}
		$info=~s/\s+$//;
		print "$a[0]\t$info\n";
	}
}
```

## Gene family analysis
### 1. Remove gene families with ≥ 100 gene copies
```bash
# (base) shichuang001@login01 Thu Nov 05 2026 21:21:47 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Orthogroups
less Orthogroups.GeneCount.tsv|perl -alne 'my $info;for (my $i = 0; $i < @F-1; $i++){$info.=$F[$i]."\t"};$info=~s/\s+$//;if(/^Orthogroup/){print "Desc\t$info"}else{print "null\t$info"}' > Carangoidei_genefamilies.tab
python clade_and_size_filter.py -i Carangoidei_genefamilies.tab -o filtered_cafe_input.txt -s # 34169 -> 32722
```

### 2. Making the species tree ultrametric
```bash
# 1. 以斑马鱼和三刺鱼的分化时间to scale branch lengths into time units
# (base) shichuang001@login01 Thu Nov 05 2026 21:49:12 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Phylogeny
python prep_r8s.py -i Carangoidei_28sp.contree -o r8s_ctl_file.txt -s 638853 -p 'Danio_rerio,Gasterosteus_aculeatus' -c '156.53'
r8s -b -f r8s_ctl_file.txt > r8s_tmp.txt
tail -n 1 r8s_tmp.txt | cut -c 16- > Carangoidei28spe_tree.txt
sbatch cafe5.sh
# Submitted batch job 910699
```

- cafe5.sh
```bash
#!/bin/bash
#SBATCH -J cafe5 # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 32 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.cafe5 # 标准输出文件
#SBATCH -e logs/%j.cafe5 # 错误输出文件
#SBATCH --mail-user=jinglkang@outlook.com

cafe5 -i filtered_cafe_input.txt -t Carangoidei28spe_tree.txt -p -c 32 -o results_single_lambda
```

**有20个Orthogroups gene copy差异过大，导致cafe5运行停止，需要去除掉**
```bash
# (base) shichuang001@login01 Thu Nov 05 2026 23:20:30 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Phylogeny
perl temp1.pl logs/910699.cafe5 > filtered_cafe_input_2.txt

# 把输入文件替换成filtered_cafe_input_2.txt
sbatch cafe5.sh
# Submitted batch job 910707
```

- temp1.pl
```perl
#!/usr/bin/perl
use strict;
use warnings;

my $log  =$ARGV[0];
my %hash;
open LOG, $log or die "can not open $log\n";
while (<LOG>) {
    chomp;
    if (/^OG/) {
        my @a=split;
        (my $id)=$a[0]=~/(.*?)\:/;
        $hash{$id}++;
    }
}

my $fami="filtered_cafe_input.txt";
open FAMI, $fami or die "can not open $fami\n";
while (<FAMI>) {
    chomp;
    my @a=split;
    unless ($hash{$a[1]}) {
    	print "$_\n";
    }
}
```

**筛选显著扩张或者收缩的基因家族**
```bash
# 1. 筛选出显著扩张或者收缩的树
# (base) shichuang001@login01 Thu Nov 05 2026 23:33:16 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Phylogeny/results_single_lambda
echo $'#nexus\nbegin trees;'>Significant_trees.tre
grep "*" Base_asr.tre >>Significant_trees.tre
echo "end;">>Significant_trees.tre

# 2. 基于Base_count.tab, Significant_trees.tre, Base_change.tab三个文件整合信息
perl extract_sig_fm.pl > Sig_family.txt

# 3. 总结全部node显著扩张和收缩的基因家族数量
# (base) shichuang001@login01 Fri Nov 06 2026 00:03:57 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Phylogeny/results_single_lambda
perl temp1.pl > Sig_family_summary.txt

# 4. R显示所有nodes的数字ID(包括内部和末端物种ID)
# 从Base_report.cafe中复制"# IDs of nodes"后面的那棵树
# 然后用print_nodes_id.R来画出一棵包括所有nodes的数字ID的系统发育树
print_nodes_id.R
```

- temp1.pl
```bash
#!/usr/bin/perl
use strict;
use warnings;

my $fmnb  ="Sig_family.txt";
my @headers;
my %hash;
open FMNB, $fmnb or die "can not open $fmnb\n";
while (<FMNB>) {
	chomp;
	my @a=split;
	if (/^FamilyID/) {
		@headers=@a;
	} else {
		for (my $i = 1; $i < @a; $i++) {
			if ($a[$i]=~/sig/) {
				(my $chan)=$a[$i]=~/\d+\((.*?)\;sig\)/;
				($chan > 0)?($hash{$headers[$i]}->{"expansion"}++):($hash{$headers[$i]}->{"constraction"}++);
			}
		}
	}
}

foreach my $id (1..@headers-1) {
	my ($expan, $const);
	($hash{$headers[$id]}->{"expansion"})?($expan=$hash{$headers[$id]}->{"expansion"}):($expan=0);
	($hash{$headers[$id]}->{"constraction"})?($const=$hash{$headers[$id]}->{"constraction"}):($const=0);
	print "$headers[$id]\t+$expan/-$const\n";
}
```

- print_nodes_id.R
```R
library(ape)
library(ggtree)
library(ggplot2)
library(stringr)

setwd("~/Documents/2026/鲳鲹进化/Gene_family")

# 1. 读取 Newick 树
newick <- "((((((((((Trachinotus_blochii<1>,Trachinotus_carolinus<2>)<42>,Trachinotus_anak<3>)<41>,Trachinotus_ovatus<4>)<40>,Scomberoides_lysan<5>)<39>,(Remorina_albescens<6>,Echeneis_naucrates<7>)<43>)<38>,((((Trachurus_japonicus<8>,Trachurus_trachurus<9>)<47>,(Decapterus_maruadsi<10>,Decapterus_macrosoma<11>)<48>)<46>,((((Megalaspis_cordyla<12>,Alepes_djedaba<13>)<52>,Selaroides_leptolepis<14>)<51>,(((Caranx_heberi<15>,Caranx_ignobilis<16>)<55>,Caranx_sexfasciatus<17>)<54>,Caranx_crysos<18>)<53>)<50>,((Alectis_indicus<19>,Caranx_equula<20>)<57>,Parastromateus_niger<21>)<56>)<49>)<45>,Elagatis_bipinnulata<22>)<44>)<37>,(((Sphyraena_novaehollandiae<23>,Sphyraena_barracuda<24>)<60>,Lates_calcarifer<25>)<59>,(Toxotes_jaculatrix<26>,Xiphias_gladius<27>)<61>)<58>)<36>,Polydactylus_sextarius<28>)<35>,(Takifugu_rubripes<29>,Gasterosteus_aculeatus<30>)<62>)<34>,Danio_rerio<31>,Oryzias_latipes<32>)<33>;"

# 2. 将尖括号中的 ID 转换为标准节点标签
newick <- gsub("<([0-9]+)>", "NID\\1", newick)

# 3. 读取树
tree <- read.tree(text = newick)

# 4. 绘制树
p <- ggtree(tree, layout = "rectangular")

# 5. 提取所有节点（包括末端和内部节点）
node_df <- p$data

# 6. 提取所有节点的 CAFE5 ID
node_df$CAFE_ID <- str_extract(
  node_df$label,
  "(?<=NID)\\d+"
)

# 7. 只保留有 CAFE5 ID 的节点
node_df <- node_df[!is.na(node_df$CAFE_ID), ]

# 8. 绘图
p +
  
  # 显示物种名称
  geom_tiplab(size = 2.8, offset = 0.15, align = FALSE) +
  
  # 显示所有节点 ID
  geom_text(data = node_df, aes(x = x, y = y, label = CAFE_ID),
    inherit.aes = FALSE, color = "red", size = 2.8,
    hjust = -0.5, vjust = -0.5) +
  
  # 扩大树右侧空间
  xlim_tree(4) +
  
  # 允许文字超出绘图区
  coord_cartesian(clip = "off") + theme_tree2() +
  theme(plot.margin = margin(
      t = 10, r = 100, b = 10, l = 10, unit = "pt"))

```
