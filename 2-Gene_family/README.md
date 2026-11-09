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

## All genes annotation
```bash
# 1. 将所有物种的pep序列合并到一个文件中
# (base) shichuang001@login01 Fri Nov 06 2026 20:29:55 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03
cat ~/jlkang/Carangoidei_evolution/Revise/longest_pep/*.fas > all.fas
# 2. 运行
# diamond blastp -q ~/jlkang/Carangoidei_evolution/Revise/longest_pep/Alectis_indicus.fas -e 1e-5 --sensitive -k 1 -d ~/jlkang/database/SwissProt.fasta --threads 100 --out Alectis_indicus.diamond.result
# (base) shichuang001@login01 Fri Nov 06 2026 20:37:42 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03
mkdir Annotation; cd Annotation
perl temp1.pl # 然后把输入内容拷贝至run_diamond.sh
# 建库
# (base) shichuang001@login01 Sat Nov 07 2026 11:15:46 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Annotation
diamond makedb --in uniprot_sprot.fasta -d uniprot_sprot.fasta
sbatch run_diamond.sh
# Submitted batch job 913209

# 3. 选出注释出现频率最高的基因名字作为Orthogroup的代表性名字
# (base) shichuang001@login01 Sat Nov 07 2026 16:48:15 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Annotation
perl generate_nm_for_orth.pl
# 产生3个文件: Orthogroup_representative_gene.txt, Orthogroup_gene_annotation.txt, Orthogroup_annotation.log

# 4. GO/KEGG注释
# (base) shichuang001@login01 Mon Nov 09 2026 13:18:14 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Annotation
conda activate eggnog
cat ~/jlkang/Carangoidei_evolution/Revise/longest_pep/*.fas > all_species.protein.fas
# emapper.py -i all_species.protein.fas --itype proteins -o fish_eggnog --cpu 100 --data_dir ~/yichengyu/eggnog_data/eggnog.db
sbatch run_emapper.sh
# Submitted batch job 913800

# 4.1 提取出每个基因的GO和KEGG注释，结果文件为Orthogroups_GOs.txt和Orthogroups_KEGG_pathways.txt
# (base) shichuang001@login01 Mon Nov 09 2026 22:25:18 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Orthogroups
perl Extract_GOs_KEGGs_orthogroups.pl

# 4.2 下载GO_id和GO_name信息
# (base) kangjingliang@KangdeMBP-2 一 11 09 2026 23:35:15 ~/Documents/2026/鲳鲹进化/Gene_family
wget https://current.geneontology.org/ontology/go-basic.obo?utm_source=chatgpt.com # 改名字为go-basic.obo
python3 extract_GOname.py # 得到对应关系表GO_ID_Name.txt
# clusterProfiler进行GO enrichment
GO_enrichment.R

# 4.3 clusterProfiler进行KEGG pathway enrichment
# 下载KEGG pathway id与name的对应表
# (base) kangjingliang@KangdeMBP-2 二 11 10 2026 00:01:02 ~/Documents/2026/鲳鲹进化/Gene_family
less Orthogroups_KEGG_pathways.txt|grep 'map' > Orthogroups_KEGG_pathways.txt.1; mv Orthogroups_KEGG_pathways.txt.1 Orthogroups_KEGG_pathways.txt
wget https://rest.kegg.jp/list/pathway?utm_source=chatgpt.com # 改名字为KEGG_Pathway_Name.txt
KEGG_enrichment.R
```

**KEGG_enrichment.R**
```R
# ============================================================
# KEGG enrichment analysis based on Orthogroups
# ============================================================

library(clusterProfiler)
library(dplyr)
library(readr)

setwd("~/Documents/2026/鲳鲹进化/Gene_family")

# 1. Read candidate OGs
candidate <- readLines(
  "Trachinotus_yes_OGs.txt",
  warn = FALSE
)

candidate <- unique(trimws(candidate))
candidate <- candidate[candidate != ""]

# 2. Read OG-KEGG Pathway annotation
anno <- read_tsv(
  "Orthogroups_KEGG_pathways.txt",
  col_names = c("OG", "Pathway"),
  show_col_types = FALSE,
  col_types = cols(
    OG = col_character(),
    Pathway = col_character()
  )
)

# 3. Clean annotation
anno <- anno %>%
  mutate(
    OG = trimws(OG),
    Pathway = trimws(Pathway)
  ) %>%
  filter(
    !is.na(OG),
    !is.na(Pathway),
    OG != "",
    Pathway != ""
  ) %>%
  distinct(OG, Pathway)

# 4. Build TERM2GENE
TERM2GENE <- anno %>%
  dplyr::select(Pathway, OG) %>%
  distinct()

# 5. Define background universe
universe <- unique(anno$OG)

# 6. Match candidate OGs with background
candidate_use <- intersect(candidate, universe)

# 7. Print statistics
cat("\n========== Basic statistics ==========\n")
cat("Total candidate OGs:", length(candidate), "\n")
cat("Total annotated OGs:", length(universe), "\n")
cat("Matched candidate OGs:", length(candidate_use), "\n")
cat("Total OG-Pathway associations:", nrow(TERM2GENE), "\n")
cat("Total unique pathways:", length(unique(TERM2GENE$Pathway)), "\n")

# Export unmatched OGs
unmatched <- setdiff(candidate, universe)

write.table(
  unmatched,
  "Unmatched_KEGG_candidate_OGs.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

if (length(candidate_use) == 0) {
  stop("No candidate OGs matched the KEGG annotation background!")
}

# 8. KEGG enrichment
ekegg <- enricher(
  gene = candidate_use,
  universe = universe,
  TERM2GENE = TERM2GENE,
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  minGSSize = 1,
  maxGSSize = 10000
)

# 9. Convert results to data frame
ekegg_df <- as.data.frame(ekegg)

cat("\nNumber of KEGG pathways returned:", nrow(ekegg_df), "\n")

if (nrow(ekegg_df) == 0) {
  
  warning("No KEGG pathways could be tested.")
  
} else {
  
  # 10. Read KEGG pathway names
  pathway_names <- read_tsv(
    "KEGG_Pathway_Name.txt",
    col_names = c("Pathway_ID", "Pathway_Name"),
    show_col_types = FALSE
  )
  
  pathway_names <- pathway_names %>%
    mutate(
      Pathway_ID = sub("^path:", "", trimws(Pathway_ID)),
      Pathway_Name = trimws(Pathway_Name)
    )
  
  # 11. Merge pathway names
  ekegg_df_named <- ekegg_df %>%
    dplyr::left_join(
      pathway_names,
      by = c("ID" = "Pathway_ID")
    ) %>%
    mutate(
      Description = ifelse(
        is.na(Pathway_Name),
        ID,
        Pathway_Name
      )
    ) %>%
    dplyr::select(
      ID,
      Description,
      everything(),
      -Pathway_Name
    ) %>%
    dplyr::arrange(p.adjust)
  
  # 12. Export all KEGG results
  write.csv(
    ekegg_df_named,
    "KEGG_all_enrichment_with_names.csv",
    row.names = FALSE
  )
  
  # 13. Export significant pathways (FDR < 0.05)
  sig_KEGG <- ekegg_df_named %>%
    dplyr::filter(p.adjust < 0.05)
  
  write.csv(
    sig_KEGG,
    "KEGG_significant_enrichment.csv",
    row.names = FALSE
  )
  
  # 14. Export pathways with FDR < 0.10
  sig_KEGG_010 <- ekegg_df_named %>%
    dplyr::filter(p.adjust < 0.10)
  
  write.csv(
    sig_KEGG_010,
    "KEGG_enrichment_FDR_010.csv",
    row.names = FALSE
  )
  
  # 15. Summary
  cat("\n========== KEGG enrichment summary ==========\n")
  cat("Total pathways:", nrow(ekegg_df_named), "\n")
  cat("Significant pathways (FDR < 0.05):",
      nrow(sig_KEGG), "\n")
  cat("Pathways with FDR < 0.10:",
      nrow(sig_KEGG_010), "\n")
}

cat("\nKEGG enrichment analysis completed!\n")

```

**GO_enrichment.R**
```R
library(clusterProfiler)
library(dplyr)
library(readr)

setwd("~/Documents/2026/鲳鲹进化/Gene_family")

# 1. Read candidate orthogroups
candidate <- readLines("Trachinotus_yes_OGs.txt",warn = FALSE)

candidate <- unique(trimws(candidate))
candidate <- candidate[candidate != ""]

# 2. Read all orthogroup GO annotations
anno <- read_tsv(
  "Orthogroups_GOs.txt",
  col_names = c("OG", "GO"),
  show_col_types = FALSE
)

anno <- anno %>%
  mutate(
    OG = trimws(OG),
    GO = trimws(GO)
  ) %>%
  filter(
    !is.na(OG),
    !is.na(GO),
    OG != "",
    GO != ""
  ) %>%
  distinct(OG, GO)

# 3. Build TERM2GENE
TERM2GENE <- anno %>%
  dplyr::select(GO, OG) %>%
  dplyr::distinct()

# 4. Background
universe <- unique(anno$OG)

# 5. Keep candidates present in annotation background
candidate_use <- intersect(candidate, universe)

# Print statistics
cat("Total candidate OGs:", length(candidate), "\n")
cat("Matched candidate OGs:", length(candidate_use), "\n")
cat("Background OGs:", length(universe), "\n")

if (length(candidate_use) == 0) {
  stop("No candidate OGs matched the GO annotation background!")
}

# 6. GO enrichment
ego <- enricher(
  gene = candidate_use,
  universe = universe,
  TERM2GENE = TERM2GENE,
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 1,
  maxGSSize = 10000
)

# 7. Convert enrichment results to data frame
ego_df <- as.data.frame(ego)

if (nrow(ego_df) == 0) {
  
  warning("No GO terms could be tested.")
  
} else {
  
  # 8. Read GO ID - GO Name mapping file
  go_names <- read_tsv(
    "GO_ID_Name.txt",
    show_col_types = FALSE
  )
  
  # Check and clean GO mapping
  go_names <- go_names %>%
    mutate(
      GO_ID = trimws(GO_ID),
      GO_Name = trimws(GO_Name)
    ) %>%
    filter(
      !is.na(GO_ID),
      !is.na(GO_Name),
      GO_ID != "",
      GO_Name != ""
    ) %>%
    distinct(GO_ID, .keep_all = TRUE)
  
  # 9. Merge GO names into enrichment results
  ego_df_named <- ego_df %>%
    dplyr::left_join(
      go_names,
      by = c("ID" = "GO_ID")
    ) %>%
    mutate(
      Description = ifelse(
        is.na(GO_Name),
        ID,
        GO_Name
      )
    ) %>%
    dplyr::select(
      ID,
      Description,
      everything(),
      -GO_Name
    ) %>%
    dplyr::arrange(p.adjust)
  
  # 10. Export all GO enrichment results
  write.csv(
    ego_df_named,
    "GO_all_enrichment_with_names.csv",
    row.names = FALSE
  )
  
  # 11. Export significant GO terms (FDR < 0.05)
  sig_GO <- ego_df_named %>%
    dplyr::filter(p.adjust < 0.05)
  
  write.csv(
    sig_GO,
    "GO_significant_enrichment.csv",
    row.names = FALSE
  )
  # 12. Print summary
  cat("\n========== GO enrichment summary ==========\n")
  cat("Total GO terms:", nrow(ego_df_named), "\n")
  cat("Significant GO terms (FDR < 0.05):", nrow(sig_GO), "\n")
}

cat("\nGO enrichment analysis completed!\n")
```

- **extract_GOname.py**
```python3
#!/usr/bin/env python3

obo_file = "go-basic.obo"
out_file = "GO_ID_Name.txt"

with open(obo_file, "r") as fin, \
     open(out_file, "w") as fout:

    fout.write("GO_ID\tGO_Name\n")

    go_id = None
    go_name = None

    for line in fin:

        line = line.strip()

        if line == "[Term]":
            go_id = None
            go_name = None

        elif line.startswith("id: GO:"):
            go_id = line.split("id: ")[1]

        elif line.startswith("name:"):
            go_name = line.split("name: ")[1]

        elif line == "" and go_id and go_name:
            fout.write(f"{go_id}\t{go_name}\n")
            go_id = None
            go_name = None

    if go_id and go_name:
        fout.write(f"{go_id}\t{go_name}\n")
```

- **Extract_GOs_KEGGs_orthogroups.pl**
```perl
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;

my (%GO, %KEGG);
my $egg="../Annotation/fish_eggnog.emapper.annotations";
open EGG, $egg or die "can not open $egg\n";
while (<EGG>) {
	chomp;
	my @a=split /\t/;
	next if /^#/;
	# GOs:9; KEGG_Pathway:12;
	$GO{$a[0]}  =$a[9];
	$KEGG{$a[0]}=$a[12];
}

my $orthgo  ="Orthogroups_GOs.txt";
my $orthkegg="Orthogroups_KEGG_pathways.txt";
open ORTHGO, ">$orthgo" or die "can not create $orthgo\n";
open ORTHKEGG, ">$orthkegg" or die "can not create $orthkegg\n";


my $orth="Orthogroups.tsv";
open ORTH, $orth or die "can not open $orth\n";
while (<ORTH>) {
	chomp;
	s/\,//g;
	my @a=split;
	if (/^Orthogroup/) {
		next;
	} else {
		my (@GOs, @KEGGs);
		for (my $i = 1; $i < @a; $i++) {
			if ($GO{$a[$i]} && $GO{$a[$i]} ne "-") {
				my @bb=split /\,/,$GO{$a[$i]};
				push @GOs, @bb;
			}

			if ($KEGG{$a[$i]} && $KEGG{$a[$i]} ne "-") {
				my @bb=split /\,/,$KEGG{$a[$i]};
				push @KEGGs, @bb;
			}
		}
		
		my (%hash1, %hash2);
		my @GOsU=grep { !$hash1{$_}++ } @GOs;
		my @KEGGsU=grep { !$hash2{$_}++ } @KEGGs;
		
		foreach my $id (@GOsU) {
			print ORTHGO "$a[0]\t$id\n";
		}

		foreach my $id (@KEGGsU) {
			print ORTHKEGG "$a[0]\t$id\n";
		}
	}
}
```

- **run_emapper.sh**
```bash
#!/bin/bash
#SBATCH -J emapper # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 100 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.emapper # 标准输出文件
#SBATCH -e logs/%j.emapper # 错误输出文件
#SBATCH --mail-user=jinglkang@outlook.com
emapper.py -i all_species.protein.fas --itype proteins -o fish_eggnog --cpu 100 --data_dir ~/yichengyu/eggnog_data
```

- **temp1.pl**
```perl
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

# diamond blastp -q ~/jlkang/Carangoidei_evolution/Revise/longest_pep/Alectis_indicus.fas -e 1e-5 --sensitive -k 1 -d ~/jlkang/database/SwissProt.fasta --threads 100 --out Alectis_indicus.diamond.result
my @fas=<~/jlkang/Carangoidei_evolution/Revise/longest_pep/*.fas>;
foreach my $fa (@fas) {
        (my $name)=basename($fa)=~/(.*)\.fa/;
        my $result=$name.".diamond.result";
        my $cmd="diamond blastp -q $fa -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out $result";
        print "$cmd\n";
}
```

- **run_diamond.sh**
```bash
#!/bin/bash
#SBATCH -J diamond # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 100 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.diamond # 标准输出文件
#SBATCH -e logs/%j.diamond # 错误输出文件
#SBATCH --mail-user=jinglkang@outlook.com

diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Alectis_indicus.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Alectis_indicus.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Alepes_djedaba.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Alepes_djedaba.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Caranx_crysos.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Caranx_crysos.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Caranx_equula.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Caranx_equula.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Caranx_heberi.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Caranx_heberi.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Caranx_ignobilis.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Caranx_ignobilis.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Caranx_sexfasciatus.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Caranx_sexfasciatus.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Danio_rerio.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Danio_rerio.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Decapterus_macrosoma.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Decapterus_macrosoma.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Decapterus_maruadsi.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Decapterus_maruadsi.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Echeneis_naucrates.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Echeneis_naucrates.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Elagatis_bipinnulata.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Elagatis_bipinnulata.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Gasterosteus_aculeatus.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Gasterosteus_aculeatus.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Lates_calcarifer.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Lates_calcarifer.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Megalaspis_cordyla.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Megalaspis_cordyla.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Oryzias_latipes.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Oryzias_latipes.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Parastromateus_niger.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Parastromateus_niger.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Polydactylus_sextarius.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Polydactylus_sextarius.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Remorina_albescens.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Remorina_albescens.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Scomberoides_lysan.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Scomberoides_lysan.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Selaroides_leptolepis.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Selaroides_leptolepis.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Sphyraena_barracuda.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Sphyraena_barracuda.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Sphyraena_novaehollandiae.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Sphyraena_novaehollandiae.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Takifugu_rubripes.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Takifugu_rubripes.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Toxotes_jaculatrix.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Toxotes_jaculatrix.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Trachinotus_anak.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Trachinotus_anak.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Trachinotus_blochii.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Trachinotus_blochii.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Trachinotus_carolinus.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Trachinotus_carolinus.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Trachinotus_ovatus.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Trachinotus_ovatus.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Trachurus_japonicus.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Trachurus_japonicus.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Trachurus_trachurus.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Trachurus_trachurus.diamond.result
diamond blastp -q /jinxianstor/home/shichuang001/jlkang/Carangoidei_evolution/Revise/longest_pep/Xiphias_gladius.fas -e 1e-5 --sensitive -k 1 -d uniprot_sprot.fasta --threads 100 --out Xiphias_gladius.diamond.result
```

- **generate_nm_for_orth.pl**
```perl
#!/usr/bin/perl
use strict;
use warnings;

# ============================================================
# Configuration
# ============================================================

my $unip_file = "uniprot_sprot.fasta";
my $orth_file = "../Orthogroups/Orthogroups.tsv";

my $out_file = "Orthogroup_gene_annotation.txt";
my $rep_file = "Orthogroup_representative_gene.txt";
my $log_file = "Orthogroup_annotation.log";


# ============================================================
# 1. Read DIAMOND results
# ============================================================

my %blastp;

my @resus = sort glob("*diamond.result");

die "ERROR: No *diamond.result files found!\n"
    unless @resus;

print "Reading DIAMOND results...\n";

foreach my $resu (@resus) {

    open my $RESU, "<", $resu
        or die "Cannot open $resu: $!\n";

    while (<$RESU>) {

        chomp;
        next if /^\s*$/;

        my @a = split /\t/;

        next unless defined $a[0] && defined $a[1];

        my $query_id   = $a[0];
        my $subject_id = $a[1];

        # Keep the first hit for each query
        $blastp{$query_id} //= $subject_id;
    }

    close $RESU;
}

print "  DIAMOND files: ", scalar(@resus), "\n";
print "  Unique query IDs: ", scalar(keys %blastp), "\n";


# ============================================================
# 2. Read UniProt FASTA annotations
# ============================================================

my %uniprot;

print "\nReading UniProt annotations...\n";

open my $UNIP, "<", $unip_file
    or die "Cannot open $unip_file: $!\n";

while (<$UNIP>) {

    chomp;
    next unless /^>/;

    my $header = substr($_, 1);

    # Example:
    # sp|Q9BWW7|SCRT1_HUMAN Transcriptional repressor scratch 1
    # OS=Homo sapiens OX=9606 GN=SCRT1 PE=1 SV=2

    my ($full_id) = $header =~ /^(\S+)/;
    next unless defined $full_id;

    # Extract accession
    my ($accession) =
        $full_id =~ /^(?:sp|tr)\|([^|]+)\|/;

    $accession = $full_id
        unless defined $accession;

    # Extract entry name
    my ($entry_name) =
        $full_id =~ /^(?:sp|tr)\|[^|]+\|([^|]+)/;

    # Extract complete protein description
    my ($description) =
        $header =~ /^\S+\s+(.*?)\s+OS=/;

    $description = "NA"
        unless defined $description && $description ne "";

    # Extract gene name
    my ($gene_name) =
        $header =~ /\bGN=([^\s]+)/;

    # Fallback if GN is absent
    unless (defined $gene_name && $gene_name ne "") {

        $gene_name = $entry_name // "Unknown";
        $gene_name =~ s/_.*$//;
        $gene_name = uc($gene_name);
    }

    # Normalize gene name to uppercase
    $gene_name = uc($gene_name);

    # Store using accession, full ID and entry name
    foreach my $id ($accession, $full_id, $entry_name) {

        next unless defined $id && $id ne "";

        $uniprot{$id}->{name} = $gene_name;
        $uniprot{$id}->{description} = $description;
    }
}

close $UNIP;

print "  UniProt annotation IDs: ", scalar(keys %uniprot), "\n";


# ============================================================
# 3. Read OrthoFinder Orthogroups.tsv
# ============================================================

print "\nReading OrthoFinder orthogroups...\n";

open my $ORTH, "<", $orth_file
    or die "Cannot open $orth_file: $!\n";

open my $OUT, ">", $out_file
    or die "Cannot create $out_file: $!\n";

open my $REP, ">", $rep_file
    or die "Cannot create $rep_file: $!\n";

open my $LOG, ">", $log_file
    or die "Cannot create $log_file: $!\n";

print $OUT
    "Orthogroup\tGeneName\tCount\tTotalGenes\tRatio\tDescription\n";

print $REP
    "Orthogroup\tRepresentativeGene\tCount\tTotalGenes\tRatio\tDescription\n";


# Summary statistics
my $total_genes = 0;
my $matched_genes = 0;
my $no_diamond_hit = 0;
my $no_uniprot_annotation = 0;
my $orthogroup_count = 0;


# ============================================================
# 4. Process each orthogroup
# ============================================================

while (<$ORTH>) {

    chomp;

    next if /^Orthogroup/;
    next if /^\s*$/;

    # Preserve empty columns
    my @a = split /\t/, $_, -1;

    my $orth_id = shift @a;

    next unless defined $orth_id && $orth_id ne "";

    $orthogroup_count++;

    my (%gene_count, %gene_description);

    my $total_in_orthogroup = 0;

    # --------------------------------------------------------
    # Process genes in this orthogroup
    # --------------------------------------------------------

    foreach my $cell (@a) {

        next unless defined $cell && $cell ne "";

        # Multiple genes in one species cell are comma-separated
        my @genes = split /,\s*/, $cell;

        foreach my $gene (@genes) {

            next unless defined $gene && $gene ne "";

            $total_in_orthogroup++;
            $total_genes++;

            # Find DIAMOND hit
            my $subject_id = $blastp{$gene};

            unless (defined $subject_id) {

                $no_diamond_hit++;

                print $LOG
                    "NO_DIAMOND_HIT\t$orth_id\t$gene\n";

                next;
            }

            # Normalize DIAMOND subject ID
            my $lookup_id = $subject_id;

            # sp|Q6GZX4|001R_FRG3G -> Q6GZX4
            if ($subject_id =~ /^(?:sp|tr)\|([^|]+)\|/) {
                $lookup_id = $1;
            }

            # BL_ORD_ID cannot be directly mapped
            if ($lookup_id =~ /^BL_ORD_ID:/) {

                $no_uniprot_annotation++;

                print $LOG
                    "UNMAPPED_DIAMOND_ID\t$orth_id\t$gene\t$subject_id\n";

                next;
            }

            # Find UniProt annotation
            unless (exists $uniprot{$lookup_id}) {

                $no_uniprot_annotation++;

                print $LOG
                    "NO_UNIPROT_ANNOTATION\t$orth_id\t$gene\t$subject_id\n";

                next;
            }

            my $gene_name =
                $uniprot{$lookup_id}->{name};

            my $description =
                $uniprot{$lookup_id}->{description};

            unless (defined $gene_name && $gene_name ne "") {

                $no_uniprot_annotation++;

                print $LOG
                    "MISSING_GENE_NAME\t$orth_id\t$gene\t$subject_id\n";

                next;
            }

            # Normalize gene name before counting
            $gene_name = uc($gene_name);

            $description = "NA"
                unless defined $description && $description ne "";

            # Count gene names
            $gene_count{$gene_name}++;

            # Keep first description for each gene name
            $gene_description{$gene_name} //= $description;

            $matched_genes++;
        }
    }


    # --------------------------------------------------------
    # 5. Output all gene names for this orthogroup
    # --------------------------------------------------------

    my @sorted_names = sort {
        $gene_count{$b} <=> $gene_count{$a}
        ||
        $a cmp $b
    } keys %gene_count;

    foreach my $gene_name (@sorted_names) {

        my $count = $gene_count{$gene_name};

        my $description =
            $gene_description{$gene_name} // "NA";

        my $ratio = sprintf(
            "%.4f",
            $count / $total_in_orthogroup
        );

        print $OUT join(
            "\t",
            $orth_id,
            $gene_name,
            $count,
            $total_in_orthogroup,
            "$count/$total_in_orthogroup",
            $description
        ), "\n";
    }


    # --------------------------------------------------------
    # 6. Output representative gene name
    # --------------------------------------------------------

    if (@sorted_names) {

        # The first name has the highest frequency.
        # Ties are resolved alphabetically.
        my $representative = $sorted_names[0];

        my $count = $gene_count{$representative};

        my $description =
            $gene_description{$representative} // "NA";

        my $ratio = sprintf(
            "%.4f",
            $count / $total_in_orthogroup
        );

        print $REP join(
            "\t",
            $orth_id,
            $representative,
            $count,
            $total_in_orthogroup,
            "$count/$total_in_orthogroup",
            $description
        ), "\n";
    }
    else {

        # No annotated genes in this orthogroup
        print $REP join(
            "\t",
            $orth_id,
            "NA",
            0,
            $total_in_orthogroup,
            "0/$total_in_orthogroup",
            "NA"
        ), "\n";
    }
}

close $ORTH;
close $OUT;
close $REP;
close $LOG;


# ============================================================
# 7. Summary
# ============================================================

print "\n========================================\n";
print "Analysis finished!\n";
print "========================================\n";

print "Orthogroups processed: $orthogroup_count\n";
print "Total genes: $total_genes\n";
print "Genes with UniProt annotation: $matched_genes\n";
print "Genes without DIAMOND hits: $no_diamond_hit\n";
print "Genes without UniProt annotation: $no_uniprot_annotation\n";

if ($total_genes > 0) {

    my $annotation_rate =
        $matched_genes / $total_genes * 100;

    printf "Annotation rate: %.2f%%\n", $annotation_rate;
}

print "\nOutput files:\n";
print "  $out_file\n";
print "  $rep_file\n";
print "  $log_file\n";
```

## 筛选出高扁平和明显非高扁平性状鱼类之间的基因家族数量差异
**四个比较组:**
```txt
Type1:
高扁平组(Toxotes jaculatrix:26);
明显非高扁平组(Sphyraena novaehollandiae:23; Sphyraena barracuda:24; Xiphias gladius:27);

Type2:
高扁平组(Trachinotus blochii:1; Trachinotus carolinus:2; Trachinotus anak:3; Trachinotus ovatus:4);
明显非高扁平组(Remorina albescens:6; Echeneis naucrates:7);

Type3:
高扁平组(Caranx heberi:15; Caranx ignobilis:16; Caranx sexfasciatus:17; Caranx crysos:18; Alectis indicus:19; Parastromateus niger:21);
明显非高扁平组(Trachurus japonicus:8; Trachurus trachurus:9; Decapterus maruadsi:10; Decapterus macrosoma:11; Megalaspis cordyla:12; Elagatis bipinnulata:22);
```

**if Group B species have more than two species, only absent for once**
**if Group B species have less two species, no absent allow**

### 1. 高扁平存在，而非高扁平完全没有
```bash
# (base) shichuang001@login01 Mon Nov 09 2026 01:04:41 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Phylogeny/results_single_lambda
# Type1: GroupA 0 gene copy number
perl filter_gene_families.pl Base_count.tab 23,24,27 26 filtered_families_Toxotes_yes.tsv
perl anno.pl filtered_families_Toxotes_yes.tsv > filtered_families_Toxotes_yes_anno.tsv

# Type2
perl filter_gene_families.pl Base_count.tab 6,7 1,2,3,4 filtered_families_Trachinotus_yes.tsv
perl anno.pl filtered_families_Trachinotus_yes.tsv > filtered_families_Trachinotus_yes_anno.tsv

# Type3
perl filter_gene_families.pl Base_count.tab 8,9,10,11,12,22 15,16,17,18,19,21 filtered_families_Caranginae_yes.tsv # 0 number
```

### 2. 高扁平完全没有，而非扁平存在
```bash
# (base) shichuang001@login01 Mon Nov 09 2026 01:04:41 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Phylogeny/results_single_lambda
# Type1
perl filter_gene_families.pl Base_count.tab 26 23,24,27 filtered_families_Toxotes_no.tsv
perl anno.pl filtered_families_Toxotes_no.tsv > filtered_families_Toxotes_no_anno.tsv

# Type2
perl filter_gene_families.pl Base_count.tab 1,2,3,4 6,7 filtered_families_Trachinotus_no.tsv
perl anno.pl filtered_families_Trachinotus_no.tsv > filtered_families_Trachinotus_no_anno.tsv

# Type3
perl filter_gene_families.pl Base_count.tab 15,16,17,18,19,21 8,9,10,11,12,22 filtered_families_Caranginae_no.tsv
perl anno.pl filtered_families_Caranginae_no.tsv > filtered_families_Caranginae_no_anno.tsv
```

### 3. 一致性状的node同时跟祖先枝node相比
```bash
# (base) shichuang001@login01 Mon Nov 09 2026 10:39:54 ~/jlkang/Carangoidei_evolution/Revise/longest_pep/OrthoFinder/Results_Nov03/Phylogeny/results_single_lambda
# Type1: 高扁平组(26); 非高扁平组(23,24,27)
perl filter_cafe_contrast_groups.pl -i Base_count.tab -a 58 -f 26 -n 23,24,27 -o Toxotes
perl anno.pl Toxotes.flat_contracted_nonflat_expanded.tsv > Toxotes.flat_contracted_nonflat_expanded_anno.tsv # 1个
perl anno.pl Toxotes.flat_expanded_nonflat_contracted.tsv > Toxotes.flat_expanded_nonflat_contracted_anno.tsv # 23个

# Type2: 高扁平组(1,2,3,4); 非高扁平组(6,7)
perl filter_cafe_contrast_groups.pl -i Base_count.tab -a 38 -f 1,2,3,4 -n 6,7 -o Trachinotus
perl anno.pl Trachinotus.flat_expanded_nonflat_contracted.tsv > Trachinotus.flat_expanded_nonflat_contracted_anno.tsv # 41个
perl anno.pl Trachinotus.flat_contracted_nonflat_expanded.tsv > Trachinotus.flat_contracted_nonflat_expanded_anno.tsv # 7个

# Type3: 高扁平组(15,16,17,18,19,21); 非高扁平组(15,16,17,18,19,21)
perl filter_cafe_contrast_groups.pl -i Base_count.tab -a 44 -f 15,16,17,18,19,21 -n 8,9,10,11,12,22 -o Caranginae # 没有找到
```


- **filter_gene_families.pl**
```perl
#!/usr/bin/perl
use strict;
use warnings;

# ============================================================
# Filter gene families based on two groups of CAFE copy numbers
#
# Group A: all selected nodes must have copy number == 0
# Group B: at least 75% of selected nodes must have copy number > 0
#
# Input:  CAFE family count matrix (tab-delimited)
# Output: FamilyID + selected Group A and Group B columns
# ============================================================

my $input = $ARGV[0];

my $group_a_input = $ARGV[1];

my $group_b_input = $ARGV[2];

my $output = $ARGV[3];

# ------------------------------------------------------------
# Parse node IDs
# ------------------------------------------------------------

sub parse_nodes {
    my ($input, $group_name) = @_;

    my @nodes = split /,/, $input;
    my %seen;

    for my $node (@nodes) {
        $node =~ s/^\s+|\s+$//g;

        die "ERROR: Invalid node ID '$node' in $group_name\n"
            unless $node =~ /^\d+$/;

        die "ERROR: Duplicate node ID '$node' in $group_name\n"
            if $seen{$node}++;

    }

    die "ERROR: $group_name is empty!\n" unless @nodes;

    return @nodes;
}

my @group_a = parse_nodes($group_a_input, "Group A");
my @group_b = parse_nodes($group_b_input, "Group B");

# Check that the two groups do not overlap
my %group_a_hash = map { $_ => 1 } @group_a;

for my $node (@group_b) {
    die "ERROR: Node $node occurs in both Group A and Group B!\n"
        if $group_a_hash{$node};
}

# if Group B species have more than two species, only absent for once
# if Group B species have less two species, no absent allow 
my $min_present;
if (@group_b > 3) {
	$min_present = @group_b - 1;
} else {
	 $min_present = @group_b;
}

print "\nGroup A nodes: @group_a\n";
print "Group B nodes: @group_b\n";
print "Minimum Group B species with copies: $min_present / "
      . scalar(@group_b) . "\n\n";

# ------------------------------------------------------------
# Read input matrix
# ------------------------------------------------------------

open my $IN, "<", $input
    or die "Cannot open input file '$input': $!\n";

my $header = <$IN>;
defined $header or die "ERROR: Input file is empty!\n";

$header =~ s/\r?\n$//;

my @headers = split /\t/, $header, -1;

# Find node IDs from column headers, e.g.
# Trachinotus_blochii<1>
# Trachinotus_carolinus<2>

my %node_column;

for my $i (0 .. $#headers) {
    if ($headers[$i] =~ /<(\d+)>\s*$/) {
        my $node = $1;

        die "ERROR: Duplicate node ID $node in header!\n"
            if exists $node_column{$node};

        $node_column{$node} = $i;
    }
}

# Check FamilyID column
die "ERROR: First column is not FamilyID!\n"
    unless $headers[0] =~ /^FamilyID$/i;

# Check whether all requested nodes exist
for my $node (@group_a, @group_b) {
    die "ERROR: Node ID $node was not found in the header!\n"
        unless exists $node_column{$node};
}

# Obtain column indices
my @a_indices = map { $node_column{$_} } @group_a;
my @b_indices = map { $node_column{$_} } @group_b;

# ------------------------------------------------------------
# Prepare output
# ------------------------------------------------------------

open my $OUT, ">", $output
    or die "Cannot create output file '$output': $!\n";

my @selected_indices = (0, @a_indices, @b_indices);

print $OUT join("\t", @headers[@selected_indices]), "\n";

# ------------------------------------------------------------
# Filter gene families
# ------------------------------------------------------------

my ($total, $selected) = (0, 0);

while (my $line = <$IN>) {

    $line =~ s/\r?\n$//;

    next if $line =~ /^\s*$/;

    my @data = split /\t/, $line, -1;

    $total++;

    die "ERROR: Incorrect number of columns at input line $total + 1\n"
        unless @data == @headers;

    # Check all Group A nodes have zero copies
    my $a_all_zero = 1;

    for my $idx (@a_indices) {

        my $value = $data[$idx];

        die "ERROR: Non-numeric copy number in Group A, family $data[0]\n"
            unless defined $value && $value =~ /^\d+$/;

        if ($value != 0) {
            $a_all_zero = 0;
            last;
        }
    }

    next unless $a_all_zero;

    # Count Group B species with copy number > 0
    my $b_present = 0;

    for my $idx (@b_indices) {

        my $value = $data[$idx];

        die "ERROR: Non-numeric copy number in Group B, family $data[0]\n"
            unless defined $value && $value =~ /^\d+$/;

        $b_present++ if $value > 0;
    }

    next if $b_present < $min_present;

    # Output FamilyID + selected Group A and Group B columns
    print $OUT join("\t", @data[@selected_indices]), "\n";

    $selected++;
}

close $IN;
close $OUT;

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

print "========================================\n";
print "Filtering completed!\n";
print "Total gene families: $total\n";
print "Selected gene families: $selected\n";
print "Group A: all copy numbers = 0\n";
print "Group B: at least $min_present / " . scalar(@group_b)
      . " species have copies\n";
print "Output file: $output\n";
print "========================================\n";

```

- **filter_cafe_contrast_groups.pl**
```perl
#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

# ============================================================
# CAFE5 contrasting gene family copy number analysis
#
# Compare flat and nonflat groups against one ancestral node.
#
# Pattern A:
#   Flat group:    expansion
#   Nonflat group: contraction
#
# Pattern B:
#   Flat group:    contraction
#   Nonflat group: expansion
#
# Automatic threshold:
#   Group size <= 3: all nodes must agree
#   Group size > 3 : allow at most one inconsistent node
#
# Header formats supported:
#   FamilyID
#   Species_name<1>
#   <32>
#   <44>
# ============================================================

my ($input, $ancestor, $flat_nodes, $nonflat_nodes, $out_prefix);

GetOptions(
    "input|i=s"    => \$input,
    "ancestor|a=s" => \$ancestor,
    "flat|f=s"     => \$flat_nodes,
    "nonflat|n=s"  => \$nonflat_nodes,
    "out|o=s"      => \$out_prefix,
) or die usage();

die usage() unless
    defined $input &&
    defined $ancestor &&
    defined $flat_nodes &&
    defined $nonflat_nodes;

$out_prefix //= "CAFE_contrast";

# ------------------------------------------------------------
# 1. Validate parameters
# ------------------------------------------------------------

die "Error: ancestor node ID must be numeric\n"
    unless $ancestor =~ /^\d+$/;

sub parse_nodes {

    my ($string, $group_name) = @_;

    my @nodes = split /,/, $string;

    die "Error: empty $group_name node list\n"
        unless @nodes;

    my %seen;

    for my $node (@nodes) {

        $node =~ s/^\s+|\s+$//g;

        die "Error: invalid node ID '$node' in $group_name\n"
            unless $node =~ /^\d+$/;

        die "Error: duplicated node ID $node in $group_name\n"
            if $seen{$node}++;
    }

    return @nodes;
}

my @flat = parse_nodes($flat_nodes, "flat group");

my @nonflat = parse_nodes($nonflat_nodes, "nonflat group");

# Check overlap between groups
my %flat_seen = map { $_ => 1 } @flat;

for my $node (@nonflat) {

    die "Error: node $node occurs in both groups\n"
        if $flat_seen{$node};
}

# Ancestor cannot be included in either group
for my $node (@flat, @nonflat) {

    die "Error: ancestor node $ancestor is included in a group\n"
        if $node eq $ancestor;
}

# ------------------------------------------------------------
# 2. Automatically determine thresholds
# ------------------------------------------------------------

my $flat_threshold =
    (@flat > 3) ? scalar(@flat) - 1 : scalar(@flat);

my $nonflat_threshold =
    (@nonflat > 3) ? scalar(@nonflat) - 1 : scalar(@nonflat);

print "========================================\n";
print "CAFE5 Contrasting Group Analysis\n";
print "========================================\n";

print "Ancestor node: $ancestor\n";

print "Flat group nodes: ",
    join(",", @flat), "\n";

print "Nonflat group nodes: ",
    join(",", @nonflat), "\n";

print "Flat group size: ", scalar(@flat), "\n";
print "Flat group threshold: $flat_threshold\n";

print "Nonflat group size: ", scalar(@nonflat), "\n";
print "Nonflat group threshold: $nonflat_threshold\n";

print "\nThreshold rules:\n";
print "  Group size <= 3: all nodes must agree\n";
print "  Group size > 3 : allow one inconsistent node\n\n";

# ------------------------------------------------------------
# 3. Read Base_count.tab
# ------------------------------------------------------------

open my $IN, "<", $input
    or die "Error: cannot open $input: $!\n";

my $header_line = <$IN>;

die "Error: input file is empty\n"
    unless defined $header_line;

chomp $header_line;
$header_line =~ s/\r$//;

# Remove UTF-8 BOM
$header_line =~ s/^\xEF\xBB\xBF//;

my @header = split /\t/, $header_line, -1;

# ------------------------------------------------------------
# 4. Parse header and extract node IDs
# ------------------------------------------------------------

my %col;

for my $i (0 .. $#header) {

    my $name = $header[$i];

    $name =~ s/^\s+|\s+$//g;

    if ($i == 0) {

        $col{"FamilyID"} = $i;
        next;
    }

    # Extract numeric node ID from:
    # Species_name<1>
    # <44>

    if ($name =~ /<(\d+)>$/) {

        my $node_id = $1;

        die "Error: duplicated node ID $node_id in header\n"
            if exists $col{$node_id};

        $col{$node_id} = $i;
    }
}

# ------------------------------------------------------------
# 5. Check requested nodes
# ------------------------------------------------------------

my @required_nodes = (
    $ancestor,
    @flat,
    @nonflat
);

my @missing;

for my $node (@required_nodes) {

    push @missing, $node
        unless exists $col{$node};
}

if (@missing) {

    print STDERR "ERROR: Missing node IDs: ";
    print STDERR join(",", @missing), "\n";

    print STDERR "\nAvailable node IDs:\n";

    my @available = sort { $a <=> $b } keys %col;

    print STDERR join(",", @available), "\n";

    die "\nPlease check Base_count.tab header.\n";
}

# ------------------------------------------------------------
# 6. Prepare output files
# ------------------------------------------------------------

my $file_A =
    "${out_prefix}.flat_expanded_nonflat_contracted.tsv";

my $file_B =
    "${out_prefix}.flat_contracted_nonflat_expanded.tsv";

my $file_all =
    "${out_prefix}.all_families.tsv";

my $file_summary =
    "${out_prefix}.summary.txt";

open my $OUT_A, ">", $file_A
    or die "Cannot write $file_A: $!\n";

open my $OUT_B, ">", $file_B
    or die "Cannot write $file_B: $!\n";

open my $OUT_ALL, ">", $file_all
    or die "Cannot write $file_all: $!\n";

open my $OUT_SUM, ">", $file_summary
    or die "Cannot write $file_summary: $!\n";

# ------------------------------------------------------------
# 7. Construct output headers
# ------------------------------------------------------------

my @out_header = (
    "FamilyID",
    "Ancestor_$ancestor",

    map { "Flat_node_$_" } @flat,

    "Flat_expanded_votes",
    "Flat_contracted_votes",
    "Flat_unchanged_votes",

    map { "Nonflat_node_$_" } @nonflat,

    "Nonflat_expanded_votes",
    "Nonflat_contracted_votes",
    "Nonflat_unchanged_votes",

    "Classification"
);

my $header = join("\t", @out_header) . "\n";

print $OUT_A $header;
print $OUT_B $header;
print $OUT_ALL $header;

# ------------------------------------------------------------
# 8. Initialize statistics
# ------------------------------------------------------------

my $total = 0;
my $valid = 0;
my $missing_count = 0;

my $pattern_A = 0;
my $pattern_B = 0;
my $other = 0;

# ------------------------------------------------------------
# 9. Analyze each gene family
# ------------------------------------------------------------

while (my $line = <$IN>) {

    chomp $line;
    $line =~ s/\r$//;

    next if $line =~ /^\s*$/;
    next if $line =~ /^\s*#/;

    my @data = split /\t/, $line, -1;

    $total++;

    my $family = $data[0];

    unless (defined $family && $family ne "") {

        $missing_count++;
        next;
    }

    # --------------------------------------------------------
    # Ancestral copy number
    # --------------------------------------------------------

    my $anc_count = $data[$col{$ancestor}];

    unless (
        defined $anc_count &&
        $anc_count =~ /^\d+(?:\.\d+)?$/
    ) {

        $missing_count++;
        next;
    }

    # --------------------------------------------------------
    # Flat group copy numbers
    # --------------------------------------------------------

    my @flat_counts;
    my $ok = 1;

    for my $node (@flat) {

        my $value = $data[$col{$node}];

        unless (
            defined $value &&
            $value =~ /^\d+(?:\.\d+)?$/
        ) {

            $ok = 0;
            last;
        }

        push @flat_counts, $value;
    }

    unless ($ok) {

        $missing_count++;
        next;
    }

    # --------------------------------------------------------
    # Nonflat group copy numbers
    # --------------------------------------------------------

    my @nonflat_counts;

    for my $node (@nonflat) {

        my $value = $data[$col{$node}];

        unless (
            defined $value &&
            $value =~ /^\d+(?:\.\d+)?$/
        ) {

            $ok = 0;
            last;
        }

        push @nonflat_counts, $value;
    }

    unless ($ok) {

        $missing_count++;
        next;
    }

    $valid++;

    # --------------------------------------------------------
    # Count expansion / contraction / unchanged
    # --------------------------------------------------------

    my ($flat_exp, $flat_con, $flat_unch) = (0, 0, 0);

    for my $count (@flat_counts) {

        if ($count > $anc_count) {

            $flat_exp++;
        }
        elsif ($count < $anc_count) {

            $flat_con++;
        }
        else {

            $flat_unch++;
        }
    }

    my ($nonflat_exp, $nonflat_con, $nonflat_unch) = (0, 0, 0);

    for my $count (@nonflat_counts) {

        if ($count > $anc_count) {

            $nonflat_exp++;
        }
        elsif ($count < $anc_count) {

            $nonflat_con++;
        }
        else {

            $nonflat_unch++;
        }
    }

    # --------------------------------------------------------
    # Classify contrasting patterns
    # --------------------------------------------------------

    my $classification = "Other";

    # Pattern A:
    # Flat expansion + Nonflat contraction

    if (
        $flat_exp >= $flat_threshold &&
        $nonflat_con >= $nonflat_threshold
    ) {

        $classification = "Flat_expanded_Nonflat_contracted";

        $pattern_A++;
    }

    # Pattern B:
    # Flat contraction + Nonflat expansion

    elsif (
        $flat_con >= $flat_threshold &&
        $nonflat_exp >= $nonflat_threshold
    ) {

        $classification = "Flat_contracted_Nonflat_expanded";

        $pattern_B++;
    }

    else {

        $other++;
    }

    # --------------------------------------------------------
    # Write output
    # --------------------------------------------------------

    my @out = (
        $family,
        $anc_count,

        @flat_counts,

        $flat_exp,
        $flat_con,
        $flat_unch,

        @nonflat_counts,

        $nonflat_exp,
        $nonflat_con,
        $nonflat_unch,

        $classification
    );

    my $out_line = join("\t", @out) . "\n";

    print $OUT_ALL $out_line;

    if (
        $classification eq
        "Flat_expanded_Nonflat_contracted"
    ) {

        print $OUT_A $out_line;
    }

    elsif (
        $classification eq
        "Flat_contracted_Nonflat_expanded"
    ) {

        print $OUT_B $out_line;
    }
}

# ------------------------------------------------------------
# 10. Close output files
# ------------------------------------------------------------

close $IN;
close $OUT_A;
close $OUT_B;
close $OUT_ALL;

# ------------------------------------------------------------
# 11. Write summary
# ------------------------------------------------------------

print $OUT_SUM "========================================\n";
print $OUT_SUM "CAFE5 Contrasting Group Analysis\n";
print $OUT_SUM "========================================\n\n";

print $OUT_SUM "Input file: $input\n";
print $OUT_SUM "Ancestor node: $ancestor\n";

print $OUT_SUM "\nFlat group nodes:\n";
print $OUT_SUM join(",", @flat), "\n";

print $OUT_SUM "Flat group size: ", scalar(@flat), "\n";
print $OUT_SUM "Flat group threshold: $flat_threshold\n";

print $OUT_SUM "\nNonflat group nodes:\n";
print $OUT_SUM join(",", @nonflat), "\n";

print $OUT_SUM "Nonflat group size: ", scalar(@nonflat), "\n";
print $OUT_SUM "Nonflat group threshold: $nonflat_threshold\n";

print $OUT_SUM "\nAutomatic threshold rule:\n";
print $OUT_SUM "Group size <= 3: all nodes must agree\n";
print $OUT_SUM "Group size > 3: allow one inconsistent node\n";

print $OUT_SUM "\nTotal families: $total\n";
print $OUT_SUM "Valid families: $valid\n";
print $OUT_SUM "Missing/invalid families: $missing_count\n";

print $OUT_SUM "\nPattern A:\n";
print $OUT_SUM "Flat expanded + Nonflat contracted: $pattern_A\n";

print $OUT_SUM "\nPattern B:\n";
print $OUT_SUM "Flat contracted + Nonflat expanded: $pattern_B\n";

print $OUT_SUM "\nOther patterns: $other\n";

close $OUT_SUM;

# ------------------------------------------------------------
# 12. Final report
# ------------------------------------------------------------

print "\n========================================\n";
print "Analysis completed!\n";
print "========================================\n";

print "Total families: $total\n";
print "Valid families: $valid\n";

print "\nPattern A:\n";
print "Flat expanded + Nonflat contracted: $pattern_A\n";

print "\nPattern B:\n";
print "Flat contracted + Nonflat expanded: $pattern_B\n";

print "\nOther patterns: $other\n";
print "Missing/invalid: $missing_count\n";

print "\nOutput files:\n";
print "  $file_A\n";
print "  $file_B\n";
print "  $file_all\n";
print "  $file_summary\n";

# ------------------------------------------------------------
# 13. Usage
# ------------------------------------------------------------

sub usage {

    return <<'USAGE';

Usage:

perl filter_cafe_contrast_groups.pl \
    -i Base_count.tab \
    -a 44 \
    -f 32,33,34,35 \
    -n 36,37,38,39 \
    -o flat_vs_nonflat

Parameters:

  -i, --input       CAFE5 Base_count.tab
  -a, --ancestor    Ancestral node ID
  -f, --flat        Flat fish node IDs
  -n, --nonflat     Nonflat fish node IDs
  -o, --out         Output prefix (optional)

Threshold rules:

  Group size <= 3:
      All nodes must show the same direction.

  Group size > 3:
      Allow at most one inconsistent node.

USAGE
}

```
