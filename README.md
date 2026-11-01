# 鲳鲹进化
## 下载数据
- 自测基因组数据
Carangiformes(鲹形目);Carangoidei(鲹亚目);
Carangidae(鲹科): Naucratinae(鰤亚科);Caranginae(鲹亚科)
Trachinotidae(鲳鲹科);Echeneidae(鮣科)

| 中文名 | 拉丁学名 | 目（Order） | 科（Family） | 属（Genus） |
|---|---|---|---|---|
| 印度丝鲹 | *Alectis indica* | 鲹形目 Carangiformes | 鲹科 Carangidae | 丝鲹属 *Alectis* |
| 六带鲹 | *Caranx sexfasciatus* | 鲹形目 Carangiformes | 鲹科 Carangidae | 鲹属 *Caranx* |
| 纺锤鰤 | *Elagatis bipinnulata* | 鲹形目 Carangiformes | 鲹科 Carangidae | 纺锤鰤属 *Elagatis* |
| 乌鲹 | *Parastromateus niger* | 鲹形目 Carangiformes | 鲹科 Carangidae | 乌鲹属 *Parastromateus* |
| 长颌似鲹 | *Scomberoides lysan* | 鲹形目 Carangiformes | 鲹科 Carangidae | 似鲹属 *Scomberoides* |
| 金带细鲹 | *Selaroides leptolepis* | 鲹形目 Carangiformes | 鲹科 Carangidae | 细鲹属 *Selaroides* |
| 大甲鲹 | *Megalaspis cordyla* | 鲹形目 Carangiformes | 鲹科 Carangidae | 大甲鲹属 *Megalaspis* |
| 长体圆鲹 | *Decapterus macrosoma* | 鲹形目 Carangiformes | 鲹科 Carangidae | 圆鲹属 *Decapterus* |
| 高体若鲹 | *Caranx equula*（建议核对） | 鲹形目 Carangiformes | 鲹科 Carangidae | 若鲹属（旧组合） |
| 吉打副叶鲹 | *Alepes djedaba* | 鲹形目 Carangiformes | 鲹科 Carangidae | 副叶鲹属 *Alepes* |
| 珍鲹? | *Caranx ignobilis* | 鲹形目 Carangiformes | 鲹科 Carangidae | 鲹属 *Caranx* |
| 蓝鳍鲹? | *Caranx melampygus* | 鲹形目 Carangiformes | 鲹科 Carangidae | 鲹属 *Caranx* |

- Naucratinae(鰤亚科);Caranginae(鲹亚科);Trachinotidae(鲳鲹科);Echeneidae(鮣科)形态上的差异总结

| 比较特征 | 鰤亚科 | 鲹亚科 | 鲳鲹科 | 鮣科 |
|:---|:---:|:---:|:---:|:---:|
| 体形 | 多为纺锤形 | 多样，纺锤形或侧扁 | 典型种类体高、侧扁 | 细长、近圆筒形 |
| 头部 | 通常无特殊附着结构 | 通常无特殊附着结构 | 头部较圆，背部轮廓明显 | 头顶具有特化吸盘 |
| 尾鳍 | 通常深叉，部分呈新月形 | 通常深叉 | 通常深叉 | 较浅叉、截形或略圆 |
| 背鳍与臀鳍 | 形态多样，部分种类呈镰刀形 | 形态多样 | 常较突出，部分呈镰刀形 | 基底长，通常无硬棘 |
| 侧线棱鳞 | 视种类而异 | 许多种类明显 | 视种类而异 | 不具有典型鲹类式棱鳞特征 |
| 主要鉴别点 | 流线型体形 | 形态多样、尾柄发达 | 体高侧扁、鳍形突出 | 头部吸盘 |

- 用于比较的系统发育关系较近的鱼类
Leptobramidae（单鳍鲳科）：单鳍鲳属 Leptobrama
Toxotidae（射水鱼科）：射水鱼属 Toxotes
Nematistiidae（长吻鲹科）：长吻鲹属 Nematistius
Menidae（月鱼科）：月鱼属 Mene
Xiphiidae（剑鱼科）：剑鱼属 Xiphias
Istiophoridae（旗鱼科）：旗鱼属 Istiophorus；枪鱼属 Makaira；四鳍旗鱼属 Kajikia；四鳍枪鱼属 Tetrapturus

## 选取基因组进行染色体进化分析
**阿纳鲳鲹_Trachinotus_anak, 卵形鲳鲹_trachinotus_ovatus, 纺锤鰤_Elagatis_bipinnulata, 长体圆鲹_Decapterus_macrosoma, 珍鲹_Caranx_ignobilis, 长印鱼_Echeneis_naucrates, 白短印鱼Remorina_albescens**
以阿纳鲳鲹_Trachinotus_anak为参考基因组
```bash
# (base) shichuang001@login01 Sun Oct 25 2026 22:46:31 ~/jlkang/Carangoidei_evolution
mkdir synteny
cp trachinotus_anna_GCA_977005145.1.fasta trachinotus_ovatus_GCA_022709315.2.fasta Elagatis_bipinnulata_assembly_final.fasta Decapterus_macrosoma_assembly_final.fasta Caranx_ignobilis.Genome.fasta synteny/
# 长印鱼_Echeneis_naucrates, 白短印鱼Remorina_albescens的基因组还在下载
# (base) shichuang001@login01 Sun Oct 25 2026 23:09:13 ~/jlkang/Carangoidei_evolution/synteny/logs
conda deactivate # 用自己安装的MUMmer
# shichuang001@login01 Sun Oct 25 2026 23:10:15 ~/jlkang/Carangoidei_evolution/synteny
sbatch nucmer.sh
# Submitted batch job 900785
```

```nucmer.sh
#!/bin/bash
#SBATCH -J nucmer # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 100 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.nucmer # 标准输出文件
#SBATCH -e logs/%j.nucmer # 错误输出文件

nucmer -b 500 -t 100 Trachinotus_anak_GCA_977005145.1.fasta Trachinotus_ovatus_GCA_022709315.2.fasta -p Tanak-Tovatus
delta-filter -1 -l 500 Tanak-Tovatus.delta > Tanak-Tovatus.delta.filt
show-coords -H -c -l -o -r -T Tanak-Tovatus.delta.filt > Tanak-Tovatus.delta.filt.coords

cat Tanak-Tovatus.delta.filt.coords |  awk '$5>2000' | awk '$3<$4{print $1"\t"$3"\tf\t"$(NF-1)"-"$NF"\n"$2"\t"$4"\tf\t"$(NF-1)"-"$NF}$3>$4{print $1"\t"$3"\tr\t"$(NF-1)"-"$NF"\n"$2"\t"$4"\tr\t"$(NF-1)"-"$NF}'  > Tanak-Tovatus.delta.filt.coords.rf
```

## 重新命名基因组序列，同时替换gff文件中的原序列名字
```bash
# (base) shichuang001@login01 Tue Oct 27 2026 16:38:30 ~/jlkang/Carangoidei_evolution/synteny
falen Caranx_ignobilis.Genome.fasta|sort -k2,2nr # 24 同样按序列长度大小排序:Decapterus_macrosoma:23; Elagatis_bipinnulata:24; 
# Trachinotus_anak:24(未排序)
# Trachinotus_ovatus:24(已排序，但是得重新命名)
# Echeneis naucrates (live sharksucker:24): GCF_900963305.2_fEcheNa1.2_genomic.fna; 得重新命名
# Remorina albescens (white suckerfish:23): GCA_034769265.1_ASM3476926v1_genomic.fna; 得重新命名

# 染色体序列根据长度全部统一重新命名
# (base) shichuang001@login01 Tue Oct 27 2026 17:35:23 ~/jlkang/Carangoidei_evolution/synteny
mkdir rename/
perl Rename.pl Caranx_ignobilis.Genome.fasta 24 > rename/Caranx_ignobilis.fasta
perl Rename.pl Decapterus_macrosoma_assembly_final.fasta 23 > rename/Decapterus_macrosoma.fasta
perl Rename.pl Elagatis_bipinnulata_assembly_final.fasta 24 > rename/Elagatis_bipinnulata.fasta
perl Rename.pl Trachinotus_anak_GCA_977005145.1.fasta 24 > rename/Trachinotus_anak.fasta
perl Rename.pl Trachinotus_ovatus_GCA_022709315.2.fasta 24 > rename/Trachinotus_ovatus.fasta
perl Rename.pl ../GCF_900963305.2_fEcheNa1.2_genomic.fna 24 > rename/Echeneis_naucrates.fasta
perl Rename.pl ../GCA_034769265.1_ASM3476926v1_genomic.fna 23 > rename/Remorina_albescens.fasta
```

- Rename.pl
```Rename.pl
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
```

```bash
# 更新gff3中的序列信息
# (base) shichuang001@login01 Tue Oct 27 2026 22:48:16 ~/jlkang/Carangoidei_evolution/synteny
perl Rename_gff.pl Trachinotus_anak_GCA_977005145.1.fasta 24 Trachinotus_anak_GCA_977005145.1.gff3 > rename/Trachinotus_anak.gff3
perl Rename_gff.pl Trachinotus_ovatus_GCA_022709315.2.fasta 24 Trachinotus_ovatus_GCA_022709315.2.gff3 > rename/Trachinotus_ovatus.gff3
perl Rename_gff.pl Caranx_ignobilis.Genome.fasta 24 ../Caranx_ignobilis.Genome.V1.gff3 > rename/Caranx_ignobilis.gff3
perl Rename_gff.pl Decapterus_macrosoma_assembly_final.fasta 23 ../Decapterus_macrosoma.Genome.V1.gff3 > rename/Decapterus_macrosoma.gff3
perl Rename_gff.pl Elagatis_bipinnulata_assembly_final.fasta 24 ../Elagatis_bipinnulata.Genome.V1.gff3 > rename/Elagatis_bipinnulata.gff3
perl Rename_gff.pl ../GCF_900963305.2_fEcheNa1.2_genomic.fna 24 ../GCF_900963305.2_fEcheNa1.2_genomic.gff > rename/Echeneis_naucrates.gff3 
```

## Pairwise synteny
**需要把所有的fasta文件都改为60个核苷酸/行**
```bash
# Trachinotus_anak vs. Trachinotus_ovatus
# Trachinotus_ovatus vs. Caranx_ignobilis
# Caranx_ignobilis vs. Decapterus_macrosoma
# Decapterus_macrosoma vs. Elagatis_bipinnulata
# Elagatis_bipinnulata vs. Echeneis_naucrates

# 1. 所有的fasta文件都改为60个核苷酸/行
# (base) shichuang001@login01 Wed Oct 28 2026 00:31:14 ~/jlkang/Carangoidei_evolution/synteny/rename
for i in *.fasta;do seqkit seq -w 60 ${i} > genome_wrapped.fa;mv genome_wrapped.fa ${i};done

# 2. Getting the CDS from the genomes using AGAT
# (base) shichuang001@login01 Tue Oct 27 2026 23:45:02 ~/jlkang/Carangoidei_evolution/synteny/rename
conda activate agat
# (base) shichuang001@login01 Tue Oct 28 2026 12:01:02 ~/jlkang/Carangoidei_evolution/synteny/rename
for fa in *.fasta; do
    prefix="${fa%.fasta}"
    gff="${prefix}.gff3"

    if [[ -f "$gff" ]]; then
        agat_sp_extract_sequences.pl \
            -g "$gff" \
            -f "$fa" \
            -t cds \
            -o "${prefix}_genome.cds.fa"
    else
        echo "Warning: $gff not found, skipping $fa" >&2
    fi
done
```
**Trachinotus_ovatus好像有点问题**


```bash
# 3. Getting BED files using AGAT
# agat_convert_sp_gff2bed.pl --gff /projects/sandbox/andre/Lepidium-synteny/data/L_campestre_rc3_functional.gff -o campestre.bed
# (agat) shichuang001@login01 Wed Oct 28 2026 00:58:06 ~/jlkang/Carangoidei_evolution/synteny/rename
for gff in *.gff3; do
    prefix="${gff%.gff3}"
    agat_convert_sp_gff2bed.pl \
            --gff "$gff" \
            -o "${prefix}.bed"    
done

# 4. reformat the FASTA files
# python -m jcvi.formats.fasta format /projects/sandbox/andre/Lepidium-synteny/data/campestre_genome.cds.fa campestre.cds
conda activate jcvi
for fa in *.fasta; do
	prefix="${fa%.fasta}"
	cdsfa="${prefix}_genome.cds.fa"
	cds="${prefix}.cds"
	python -m jcvi.formats.fasta format \
        $cdsfa \
        $cds
done

# 4. Pairwise synteny search
# Caranx_ignobilis vs. Decapterus_macrosoma
# Decapterus_macrosoma vs. Echeneis_naucrates
# Echeneis_naucrates vs. Elagatis_bipinnulata
# Elagatis_bipinnulata vs. Trachinotus_anak
# Trachinotus_anak vs. Trachinotus_ovatus

#!/bin/bash
#SBATCH -J jcvi # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 20 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.jcvi # 标准输出文件
#SBATCH -e logs/%j.jcvi # 错误输出文件
echo ${SLURM_JOB_NODELIST} # 作业占用节点列表
echo start on $(date) # 开始时间
python -m jcvi.compara.catalog ortholog Caranx_ignobilis Decapterus_macrosoma --no_strip_names
python -m jcvi.compara.catalog ortholog Decapterus_macrosoma Echeneis_naucrates --no_strip_names
python -m jcvi.compara.catalog ortholog Echeneis_naucrates Elagatis_bipinnulata --no_strip_names
python -m jcvi.compara.catalog ortholog Elagatis_bipinnulata Trachinotus_anak --no_strip_names
python -m jcvi.compara.catalog ortholog Trachinotus_anak Trachinotus_ovatus --no_strip_names
echo end on $(date) # 结束时间

# 5. generates a .simple file from the .anchors file
python -m jcvi.compara.synteny screen --minspan=30 --simple Caranx_ignobilis.Decapterus_macrosoma.anchors Caranx_ignobilis.Decapterus_macrosoma.anchors.new
python -m jcvi.compara.synteny screen --minspan=30 --simple Decapterus_macrosoma.Echeneis_naucrates.anchors Decapterus_macrosoma.Echeneis_naucrates.anchors.new
python -m jcvi.compara.synteny screen --minspan=30 --simple Echeneis_naucrates.Elagatis_bipinnulata.anchors Echeneis_naucrates.Elagatis_bipinnulata.anchors.new
python -m jcvi.compara.synteny screen --minspan=30 --simple Elagatis_bipinnulata.Trachinotus_anak.anchors Elagatis_bipinnulata.Trachinotus_anak.anchors.new
python -m jcvi.compara.synteny screen --minspan=30 --simple Trachinotus_anak.Trachinotus_ovatus.anchors Trachinotus_anak.Trachinotus_ovatus.anchors.new
```

**6.seqids**
```seqids
Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8,Chr9,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24
Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8,Chr9,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23
Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8,Chr9,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24
Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8,Chr9,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24
Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8,Chr9,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24
Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8,Chr9,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24
```

**7.layout**
```layout
# y, xstart, xend, rotation, color, label, va,  bed
 .85,     .1,    .8,       0,      , Caranx_ignobilis, top, Caranx_ignobilis.bed
 .7,     .1,    .8,       0,      , Decapterus_macrosoma, top, Decapterus_macrosoma.bed
 .55,     .1,    .8,       0,      , Echeneis_naucrates, top, Echeneis_naucrates.bed
 .4,     .1,    .8,       0,      , Elagatis_bipinnulata, top, Elagatis_bipinnulata.bed
 .25,     .1,    .8,       0,      , Trachinotus_anak, top, Trachinotus_anak.bed
 .1,     .1,    .8,       0,      , Trachinotus_ovatus, top, Trachinotus_ovatus.bed
# edges
e, 0, 1, Caranx_ignobilis.Decapterus_macrosoma.anchors.simple
e, 1, 2, Decapterus_macrosoma.Echeneis_naucrates.anchors.simple
e, 2, 3, Echeneis_naucrates.Elagatis_bipinnulata.anchors.simple
e, 3, 4, Elagatis_bipinnulata.Trachinotus_anak.anchors.simple
e, 4, 5, Trachinotus_anak.Trachinotus_ovatus.anchors.simple
```

**8.Plot**
```bash
python -m jcvi.graphics.karyotype seqids layout
```

**9.标记出重要的synteny**
- marker_synteny.pl
```marker_synteny.pl
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

my $bed=$ARGV[0]; my (%hash1, %hash2);
open BED, $bed or die "can not open $bed\n";
while (<BED>) {
	chomp;
	my @a=split;
	$hash1{$a[3]}=$a[0];
}

my $simp=$ARGV[1];
my @chrs=split /\,/, $ARGV[2];
foreach my $chr (@chrs) {
	$hash2{$chr}++;
}

open SIMP, $simp or die "can not open $simp\n";
while (<SIMP>) {
	chomp;
	my @a=split;
	my $chr=$hash1{$a[0]};
	if ($hash2{$chr}) {
		print "m*$_\n";
	} else {
		print "$_\n";
	}
}
```

```bash
perl marker_synteny.pl Caranx_ignobilis.bed Caranx_ignobilis.Decapterus_macrosoma.anchors.simple Chr7,Chr24 > Caranx_ignobilis.Decapterus_macrosoma.anchors.simple.1; mv Caranx_ignobilis.Decapterus_macrosoma.anchors.simple.1 Caranx_ignobilis.Decapterus_macrosoma.anchors.simple
perl marker_synteny.pl Decapterus_macrosoma.bed Decapterus_macrosoma.Echeneis_naucrates.anchors.simple Chr1 > Decapterus_macrosoma.Echeneis_naucrates.anchors.simple.1; mv Decapterus_macrosoma.Echeneis_naucrates.anchors.simple.1 Decapterus_macrosoma.Echeneis_naucrates.anchors.simple
perl marker_synteny.pl Echeneis_naucrates.bed Echeneis_naucrates.Elagatis_bipinnulata.anchors.simple Chr9,Chr17 > Echeneis_naucrates.Elagatis_bipinnulata.anchors.simple.1; mv Echeneis_naucrates.Elagatis_bipinnulata.anchors.simple.1 Echeneis_naucrates.Elagatis_bipinnulata.anchors.simple
perl marker_synteny.pl Elagatis_bipinnulata.bed Elagatis_bipinnulata.Trachinotus_anak.anchors.simple Chr6,Chr24 > Elagatis_bipinnulata.Trachinotus_anak.anchors.simple.1; mv Elagatis_bipinnulata.Trachinotus_anak.anchors.simple.1 Elagatis_bipinnulata.Trachinotus_anak.anchors.simple
perl marker_synteny.pl Trachinotus_anak.bed Trachinotus_anak.Trachinotus_ovatus.anchors.simple Chr5,Chr24 > Trachinotus_anak.Trachinotus_ovatus.anchors.simple.1; mv Trachinotus_anak.Trachinotus_ovatus.anchors.simple.1 Trachinotus_anak.Trachinotus_ovatus.anchors.simple

# plot
python -m jcvi.graphics.karyotype seqids layout
```

```bash
# fission
perl marker_synteny_fission.pl Caranx_ignobilis.bed Caranx_ignobilis.Decapterus_macrosoma.anchors.simple Chr20 > Caranx_ignobilis.Decapterus_macrosoma.anchors.simple.1; mv Caranx_ignobilis.Decapterus_macrosoma.anchors.simple.1 Caranx_ignobilis.Decapterus_macrosoma.anchors.simple
perl marker_synteny_fission.pl Decapterus_macrosoma.bed Decapterus_macrosoma.Echeneis_naucrates.anchors.simple Chr20 > Decapterus_macrosoma.Echeneis_naucrates.anchors.simple.1; mv Decapterus_macrosoma.Echeneis_naucrates.anchors.simple.1 Decapterus_macrosoma.Echeneis_naucrates.anchors.simple
perl marker_synteny_fission.pl Echeneis_naucrates.bed Echeneis_naucrates.Elagatis_bipinnulata.anchors.simple Chr20 > Echeneis_naucrates.Elagatis_bipinnulata.anchors.simple.1; mv Echeneis_naucrates.Elagatis_bipinnulata.anchors.simple.1 Echeneis_naucrates.Elagatis_bipinnulata.anchors.simple
perl marker_synteny_fission.pl Elagatis_bipinnulata.bed Elagatis_bipinnulata.Trachinotus_anak.anchors.simple Chr19 > Elagatis_bipinnulata.Trachinotus_anak.anchors.simple.1; mv Elagatis_bipinnulata.Trachinotus_anak.anchors.simple.1 Elagatis_bipinnulata.Trachinotus_anak.anchors.simple
perl marker_synteny_fission.pl Trachinotus_anak.bed Trachinotus_anak.Trachinotus_ovatus.anchors.simple Chr18 > Trachinotus_anak.Trachinotus_ovatus.anchors.simple.1; mv Trachinotus_anak.Trachinotus_ovatus.anchors.simple.1 Trachinotus_anak.Trachinotus_ovatus.anchors.simple

# fussion
perl marker_synteny_fussion.1.pl Caranx_ignobilis.bed Caranx_ignobilis.Decapterus_macrosoma.anchors.simple Chr19 > Caranx_ignobilis.Decapterus_macrosoma.anchors.simple.1; mv Caranx_ignobilis.Decapterus_macrosoma.anchors.simple.1 Caranx_ignobilis.Decapterus_macrosoma.anchors.simple
perl marker_synteny_fussion.1.pl Decapterus_macrosoma.bed Decapterus_macrosoma.Echeneis_naucrates.anchors.simple Chr19 > Decapterus_macrosoma.Echeneis_naucrates.anchors.simple.1; mv Decapterus_macrosoma.Echeneis_naucrates.anchors.simple.1 Decapterus_macrosoma.Echeneis_naucrates.anchors.simple
perl marker_synteny_fussion.1.pl Echeneis_naucrates.bed Echeneis_naucrates.Elagatis_bipinnulata.anchors.simple Chr18 > Echeneis_naucrates.Elagatis_bipinnulata.anchors.simple.1; mv Echeneis_naucrates.Elagatis_bipinnulata.anchors.simple.1 Echeneis_naucrates.Elagatis_bipinnulata.anchors.simple
perl marker_synteny_fussion.1.pl Elagatis_bipinnulata.bed Elagatis_bipinnulata.Trachinotus_anak.anchors.simple Chr18 > Elagatis_bipinnulata.Trachinotus_anak.anchors.simple.1; mv Elagatis_bipinnulata.Trachinotus_anak.anchors.simple.1 Elagatis_bipinnulata.Trachinotus_anak.anchors.simple
perl marker_synteny_fussion.1.pl Trachinotus_anak.bed Trachinotus_anak.Trachinotus_ovatus.anchors.simple Chr19 > Trachinotus_anak.Trachinotus_ovatus.anchors.simple.1; mv Trachinotus_anak.Trachinotus_ovatus.anchors.simple.1 Trachinotus_anak.Trachinotus_ovatus.anchors.simple

# fussion
perl marker_synteny_fussion.2.pl Caranx_ignobilis.bed Caranx_ignobilis.Decapterus_macrosoma.anchors.simple Chr13 > Caranx_ignobilis.Decapterus_macrosoma.anchors.simple.1; mv Caranx_ignobilis.Decapterus_macrosoma.anchors.simple.1 Caranx_ignobilis.Decapterus_macrosoma.anchors.simple
perl marker_synteny_fussion.2.pl Decapterus_macrosoma.bed Decapterus_macrosoma.Echeneis_naucrates.anchors.simple Chr13 > Decapterus_macrosoma.Echeneis_naucrates.anchors.simple.1; mv Decapterus_macrosoma.Echeneis_naucrates.anchors.simple.1 Decapterus_macrosoma.Echeneis_naucrates.anchors.simple
perl marker_synteny_fussion.2.pl Echeneis_naucrates.bed Echeneis_naucrates.Elagatis_bipinnulata.anchors.simple Chr14 > Echeneis_naucrates.Elagatis_bipinnulata.anchors.simple.1; mv Echeneis_naucrates.Elagatis_bipinnulata.anchors.simple.1 Echeneis_naucrates.Elagatis_bipinnulata.anchors.simple
perl marker_synteny_fussion.2.pl Elagatis_bipinnulata.bed Elagatis_bipinnulata.Trachinotus_anak.anchors.simple Chr14 > Elagatis_bipinnulata.Trachinotus_anak.anchors.simple.1; mv Elagatis_bipinnulata.Trachinotus_anak.anchors.simple.1 Elagatis_bipinnulata.Trachinotus_anak.anchors.simple
perl marker_synteny_fussion.2.pl Trachinotus_anak.bed Trachinotus_anak.Trachinotus_ovatus.anchors.simple Chr12 > Trachinotus_anak.Trachinotus_ovatus.anchors.simple.1; mv Trachinotus_anak.Trachinotus_ovatus.anchors.simple.1 Trachinotus_anak.Trachinotus_ovatus.anchors.simple

# plot
python -m jcvi.graphics.karyotype seqids layout
```

## 补充数据：下载所有Caranginae，Naucratinae，Echeneidae，Trachinotidae的基因组数据，以Toxotes，Xiphiidae，Istiophoridae里面的物种作为外类群（只要有基因组数据都下载下来）
- species_download.txt
```species_download.txt
Caranx heberi
Caranx bucculentus
Caranx crysos
Alepes apercna
Alepes kleinii
Gnathanodon speciosus
Selene dorsalis
Decapterus maruadsi
Trachurus japonicus
Trachurus trachurus
Pseudocaranx dentex
Seriola dumerili
Seriolina nigrofasciata
Seriola lalandi
Seriola quinqueradiata
Trachinotus baillonii
Toxotes jaculatrix
Toxotes lorentzi
Xiphias gladius
Istiompax indica
```

```bash
# (jcvi) shichuang001@login01 Thu Oct 29 2026 13:05:05 ~/jlkang/Carangoidei_evolution
conda activate ncbi_datasets
datasets summary genome taxon --inputfile species_download.txt \
--reference \
--as-json-lines | \
dataformat tsv genome \
--fields accession,assminfo-name,assminfo-level,assmstats-total-number-of-chromosomes,assmstats-number-of-scaffolds,assmstats-total-sequence-len,assmstats-contig-n50,assmstats-scaffold-n50,annotinfo-status,organism-name \
> species_download_genomes.tsv
```

**筛选出染色体组装级别的基因组**
Caranx crysos,Caranx heberi,Decapterus maruadsi,Toxotes jaculatrix,Trachurus japonicus,Trachurus trachurus,Xiphias gladius
其中: Toxotes jaculatrix,Trachurus japonicus,Trachurus trachurus,Xiphias gladius都已经有注释信息了

- Taget_annotated_genomes.txt
```Taget_annotated_genomes.txt
GCF_017976425.1
GCF_045865225.1
GCF_905171665.1
GCF_016859285.1
```

- Taget_noannotated_genomes.txt
```Taget_noannotated_genomes.txt
GCA_984588375.1
GCA_054643105.1
GCA_030347415.2
```

根据Assembly Accession的号进行下载
```bash
# 1. Taget_annotated_genomes.txt
# (ncbi_datasets) kangjingliang@KangdeMBP-2 四 10 29 2026 13:18:01 ~/Documents/2026/鲳鲹进化
conda activate ncbi_datasets
datasets download genome accession --inputfile Taget_annotated_genomes.txt \
--include genome,gff3 \
--dehydrated --filename Taget_annotated_genomes.zip

unzip Taget_annotated_genomes.zip -d Taget_annotated_genomes

datasets rehydrate --directory Taget_annotated_genomes/

# 2. Taget_noannotated_genomes.txt
datasets download genome accession --inputfile Taget_noannotated_genomes.txt \
--include genome,gff3 \
--dehydrated --filename Taget_noannotated_genomes.zip

unzip Taget_noannotated_genomes.zip -d Taget_noannotated_genomes

datasets rehydrate --directory Taget_noannotated_genomes/
```

## Genome annotation
**Decapterus maruadsi:GCA_030347415.2**
```bash
# Decapterus maruadsi: 已经是soft mask
# (base) shichuang001@login01 Thu Oct 29 2026 20:43:59 ~/jlkang/Carangoidei_evolution/Taget_noannotated_genomes/ncbi_dataset/data/GCA_030347415.2
conda activate base
sbatch braker.sh
# Submitted batch job 906151
```

- braker.sh
```braker.sh
#!/bin/bash
#SBATCH -J braker_Decapterus_maruadsi # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 30 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.braker_Decapterus_maruadsi # 标准输出文件
#SBATCH -e logs/%j.braker_Decapterus_maruadsi # 错误输出文件

braker.pl --softmasking --AUGUSTUS_ab_initio --gff3 \
--species Decapterus_maruadsi \
--threads 30 \
--genome=GCA_030347415.2_ASM3034741v2_genomic.fna \
--prot_seq=/jinxianstor/home/shichuang001/chentongyu/databases/orthodb/Vertebrata.fa \
--workingdir ./
```

**Caranx heberi:GCA_054643105.1**
```bash
# Caranx heberi: GCA_054643105.1
# (base) shichuang001@login01 Thu Oct 29 2026 21:51:27 ~/jlkang/Carangoidei_evolution/Taget_noannotated_genomes/ncbi_dataset/data/GCA_054643105.1
sbatch braker.sh
# Submitted batch job 906155
```

- braker.sh
```braker.sh
#!/bin/bash
#SBATCH -J braker_Caranx_heberi # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 30 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.braker_Caranx_heberi # 标准输出文件
#SBATCH -e logs/%j.braker_Caranx_heberi # 错误输出文件

braker.pl --softmasking --AUGUSTUS_ab_initio --gff3 \
--species Caranx_heberi \
--threads 30 \
--genome=GCA_054643105.1_CIBA_Cheb_1.1_genomic.fna \
--prot_seq=/jinxianstor/home/shichuang001/chentongyu/databases/orthodb/Vertebrata.fa \
--workingdir ./
```

**Caranx_crysos:GCA_984588375.1**
```bash
# Caranx_crysos:GCA_984588375.1
# (base) shichuang001@login01 Thu Oct 29 2026 21:56:47 ~/jlkang/Carangoidei_evolution/Taget_noannotated_genomes/ncbi_dataset/data/GCA_984588375.1
sbatch braker.sh
# Submitted batch job 906158
```

- braker.sh
```braker.sh
#!/bin/bash
#SBATCH -J braker_Caranx_crysos # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 30 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.braker_Caranx_crysos # 标准输出文件
#SBATCH -e logs/%j.braker_Caranx_crysos # 错误输出文件

braker.pl --softmasking --AUGUSTUS_ab_initio --gff3 \
--species Caranx_crysos \
--threads 30 \
--genome=GCA_984588375.1_fCarCry1.1_genomic.fna \
--prot_seq=/jinxianstor/home/shichuang001/chentongyu/databases/orthodb/Vertebrata.fa \
--workingdir ./
```

**Remorina_albescens:GCA_034769265.1**
```bash
# (base) shichuang001@login01 Thu Oct 29 2026 22:05:33 ~/jlkang/Carangoidei_evolution/Remorina_albescens
sbatch braker.sh
# Submitted batch job 906164
```

- braker.sh
```braker.sh
#!/bin/bash
#SBATCH -J braker_Remorina_albescens # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 30 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.braker_Remorina_albescens # 标准输出文件
#SBATCH -e logs/%j.braker_Remorina_albescens # 错误输出文件

braker.pl --softmasking --AUGUSTUS_ab_initio --gff3 \
--species Remorina_albescens \
--threads 30 \
--genome=GCA_034769265.1_ASM3476926v1_genomic.fna \
--prot_seq=/jinxianstor/home/shichuang001/chentongyu/databases/orthodb/Vertebrata.fa \
--workingdir ./
```

**Trachinotus_ovatus的gff不包括scaffold，重新再注释一下**
```bash
# (base) shichuang001@login01 Thu Oct 29 2026 22:12:20 ~/jlkang/Carangoidei_evolution/synteny/rename/Trachinotus_ovatus
sbatch braker.sh
# Submitted batch job 906165
```

- braker.sh
```bash
#!/bin/bash
#SBATCH -J braker_Trachinotus_ovatus # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 30 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.braker_Trachinotus_ovatus # 标准输出文件
#SBATCH -e logs/%j.braker_Trachinotus_ovatus # 错误输出文件

braker.pl --softmasking --AUGUSTUS_ab_initio --gff3 \
--species Trachinotus_ovatus \
--threads 30 \
--genome=Trachinotus_ovatus.fasta \
--prot_seq=/jinxianstor/home/shichuang001/chentongyu/databases/orthodb/Vertebrata.fa \
--workingdir ./
```

**再增加四个物种:Lates calcarifer (GCF_001640805.2, 有gff), Sphyraena barracuda (GCA_988282095.1), Sphyraena novaehollandiae (GCA_050436865.1), Polydactylus sextarius (GCA_016801845.1)**
**现在鲹形目所有具有染色体级别基因组的鱼类都已经凑齐**

**Polydactylus sextarius (GCA_016801845.1)**
```bash
# (base) shichuang001@login01 Fri Oct 30 2026 17:42:52 ~/jlkang/Carangoidei_evolution
mkdir extra_genomes
sbatch braker.sh
# Submitted batch job 906702
```

- braker.sh
```bash
#!/bin/bash
#SBATCH -J Polydactylus_sextarius # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 30 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.braker_Polydactylus_sextarius # 标准输出文件
#SBATCH -e logs/%j.braker_Polydactylus_sextarius # 错误输出文件

braker.pl --softmasking --AUGUSTUS_ab_initio --gff3 \
--species Polydactylus_sextarius \
--threads 30 \
--genome=GCA_016801845.1_ASM1680184v1_genomic.fna \
--prot_seq=/jinxianstor/home/shichuang001/chentongyu/databases/orthodb/Vertebrata.fa \
--workingdir ./
```

**Sphyraena novaehollandiae (GCA_050436865.1)**
```bash
# (base) shichuang001@login01 Fri Oct 30 2026 18:05:20 ~/jlkang/Carangoidei_evolution/extra_genomes/GCA_050436865.1
sbatch braker.sh
# Submitted batch job 906952
```

- braker.sh
```bash
#!/bin/bash
#SBATCH -J Sphyraena_novaehollandiae # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 30 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.braker_Sphyraena_novaehollandiae # 标准输出文件
#SBATCH -e logs/%j.braker_Sphyraena_novaehollandiae # 错误输出文件

braker.pl --softmasking --AUGUSTUS_ab_initio --gff3 \
--species Sphyraena_novaehollandiae \
--threads 30 \
--genome=GCA_050436865.1_fSphNov1.hap2.cur.20230825_genomic.fna \
--prot_seq=/jinxianstor/home/shichuang001/chentongyu/databases/orthodb/Vertebrata.fa \
--workingdir ./
```

**Sphyraena barracuda (GCA_988282095.1)**
```bash
# (base) shichuang001@login01 Fri Oct 30 2026 18:08:35 ~/jlkang/Carangoidei_evolution/extra_genomes/GCA_988282095.1
sbatch braker.sh
# Submitted batch job 906946
```

- braker.sh
```bash
#!/bin/bash
#SBATCH -J Sphyraena_barracuda # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 30 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.braker_Sphyraena_barracuda # 标准输出文件
#SBATCH -e logs/%j.braker_Sphyraena_barracuda # 错误输出文件

braker.pl --softmasking --AUGUSTUS_ab_initio --gff3 \
--species Sphyraena_barracuda \
--threads 30 \
--genome=GCA_988282095.1_Sphyraena_barracuda_haplotype_1_assembly_genomic.fna \
--prot_seq=/jinxianstor/home/shichuang001/chentongyu/databases/orthodb/Vertebrata.fa \
--workingdir ./
```

**重新下载卵形鲳鲹基因组(T2T:GCA_061193975.1_ASM6119397v1)**
```bash
# h2076@h2076 Fri Oct 30 2026 20:49:49 ~/Carangoidei_evolution
nohup wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/061/193/975/GCA_061193975.1_ASM6119397v1/GCA_061193975.1_ASM6119397v1_genomic.fna.gz > GCA_061193975_download.txt 2>&1 &
# [1] 484994

# (base) shichuang001@login01 Sat Oct 31 2026 10:17:22 ~/jlkang/Carangoidei_evolution/Trachinotus_ovatus_T2T
sbatch braker.sh
# Submitted batch job 907438
```

- braker.sh
```bash
#!/bin/bash
#SBATCH -J Trachinotus_ovatus # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 30 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.braker_Trachinotus_ovatus # 标准输出文件
#SBATCH -e logs/%j.braker_Trachinotus_ovatus # 错误输出文件

braker.pl --softmasking --AUGUSTUS_ab_initio --gff3 \
--species Trachinotus_ovatus_T2T \
--threads 30 \
--genome=GCA_061193975.1_ASM6119397v1_genomic.fna \
--prot_seq=/jinxianstor/home/shichuang001/chentongyu/databases/orthodb/Vertebrata.fa \
--workingdir ./
```

## braker注释的名字是把fasta的header全部连接到一起了，重新改一下
```bash
# (base) shichuang001@login01 Sun Nov 01 2026 12:50:00 ~/jlkang/Carangoidei_evolution
mkdir rename

# Remorina_albescens: 23 chromosomes
# (base) shichuang001@login01 Sun Nov 01 2026 11:37:05 ~/jlkang/Carangoidei_evolution/Remorina_albescens
falen Remorina_albescens.fasta|sort -k2,2nr
Rename.pl Remorina_albescens.fasta 23 > ../rename/Remorina_albescens.fasta
Rename_gff_braker.pl Remorina_albescens.fasta 23 Remorina_albescens.gff3 > ../rename/Remorina_albescens.gff3

# Xiphias_gladius: 24 chromosomes
# (base) shichuang001@login01 Sun Nov 01 2026 14:42:15 ~/jlkang/Carangoidei_evolution/Taget_annotated_genomes/ncbi_dataset/data/GCF_016859285.1
Rename.pl Xiphias_gladius.fasta 24 > ~/jlkang/Carangoidei_evolution/rename/Xiphias_gladius.fasta
Rename_gff.pl Xiphias_gladius.fasta 24 Xiphias_gladius.gff3 > ~/jlkang/Carangoidei_evolution/rename/Xiphias_gladius.gff3

# Toxotes_jaculatrix: 24 chromosomes
# (base) shichuang001@login01 Sun Nov 01 2026 17:26:32 ~/jlkang/Carangoidei_evolution/Taget_annotated_genomes/ncbi_dataset/data/GCF_017976425.1
Rename.pl Toxotes_jaculatrix.fasta 24 > ~/jlkang/Carangoidei_evolution/rename/Toxotes_jaculatrix.fasta
Rename_gff.pl Toxotes_jaculatrix.fasta 24 Toxotes_jaculatrix.gff3 > ~/jlkang/Carangoidei_evolution/rename/Toxotes_jaculatrix.gff3

# Trachurus_japonicus: 24 chromosomes
# (base) shichuang001@login01 Sun Nov 01 2026 17:29:26 ~/jlkang/Carangoidei_evolution/Taget_annotated_genomes/ncbi_dataset/data/GCF_045865225.1
Rename.pl Trachurus_japonicus.fasta 24 > ~/jlkang/Carangoidei_evolution/rename/Trachurus_japonicus.fasta
Rename_gff.pl Trachurus_japonicus.fasta 24 Trachurus_japonicus.gff3 > ~/jlkang/Carangoidei_evolution/rename/Trachurus_japonicus.gff3

# Trachurus_trachurus: 24 chromosomes
# (base) shichuang001@login01 Sun Nov 01 2026 17:36:04 ~/jlkang/Carangoidei_evolution/Taget_annotated_genomes/ncbi_dataset/data/GCF_905171665.1
Rename.pl Trachurus_trachurus.fasta 24 > ~/jlkang/Carangoidei_evolution/rename/Trachurus_trachurus.fasta
Rename_gff.pl Trachurus_trachurus.fasta 24 Trachurus_trachurus.gff3 > ~/jlkang/Carangoidei_evolution/rename/Trachurus_trachurus.gff3
```

**braker注释**
```bash
# Decapterus_maruadsi: 23 chromosomes
# (base) shichuang001@login01 Sun Nov 01 2026 17:43:26 ~/jlkang/Carangoidei_evolution/Taget_noannotated_genomes/ncbi_dataset/data/GCA_030347415.2
Rename.pl Decapterus_maruadsi.fasta 23 > ~/jlkang/Carangoidei_evolution/rename/Decapterus_maruadsi.fasta
Rename_gff_braker.pl Decapterus_maruadsi.fasta 23 Decapterus_maruadsi.gff3 > ~/jlkang/Carangoidei_evolution/rename/Decapterus_maruadsi.gff3

# Caranx_heberi: 24 chromosomes
# (base) shichuang001@login01 Sun Nov 01 2026 17:49:49 ~/jlkang/Carangoidei_evolution/Taget_noannotated_genomes/ncbi_dataset/data/GCA_054643105.1
Rename.pl Caranx_heberi.fasta 24 > ~/jlkang/Carangoidei_evolution/rename/Caranx_heberi.fasta
Rename_gff_braker.pl Caranx_heberi.fasta 24 Caranx_heberi.gff3 > ~/jlkang/Carangoidei_evolution/rename/Caranx_heberi.gff3

# Caranx_crysos: 24 chromosomes
# (base) shichuang001@login01 Sun Nov 01 2026 17:53:03 ~/jlkang/Carangoidei_evolution/Taget_noannotated_genomes/ncbi_dataset/data/GCA_984588375.1
Rename.pl Caranx_crysos.fasta 24 > ~/jlkang/Carangoidei_evolution/rename/Caranx_crysos.fasta
Rename_gff_braker.pl Caranx_crysos.fasta 24 Caranx_crysos.gff3 > ~/jlkang/Carangoidei_evolution/rename/Caranx_crysos.gff3

# Trachinotus_ovatus_T2T: 24 chromosomes
# (base) shichuang001@login01 Sun Nov 01 2026 17:56:28 ~/jlkang/Carangoidei_evolution/Trachinotus_ovatus_T2T
Rename.pl GCA_061193975.1_ASM6119397v1_genomic.fna 24 > ~/jlkang/Carangoidei_evolution/rename/Trachinotus_ovatus.fasta
Rename_gff_braker.pl GCA_061193975.1_ASM6119397v1_genomic.fna 24 braker.gff3 > ~/jlkang/Carangoidei_evolution/rename/Trachinotus_ovatus.gff3

# Polydactylus sextarius: 26 chromosomes
# (base) shichuang001@login01 Sun Nov 01 2026 17:58:54 ~/jlkang/Carangoidei_evolution/extra_genomes/GCA_016801845.1
Rename.pl GCA_016801845.1_ASM1680184v1_genomic.fna 26 > ~/jlkang/Carangoidei_evolution/rename/Polydactylus_sextarius.fasta
Rename_gff_braker.pl GCA_016801845.1_ASM1680184v1_genomic.fna 26 braker.gff3 > ~/jlkang/Carangoidei_evolution/rename/Polydactylus_sextarius.gff3

# Sphyraena novaehollandiae: 24 chromosomes
# (base) shichuang001@login01 Sun Nov 01 2026 18:11:50 ~/jlkang/Carangoidei_evolution/extra_genomes/GCA_050436865.1
Rename.pl GCA_050436865.1_fSphNov1.hap2.cur.20230825_genomic.fna 24 > ~/jlkang/Carangoidei_evolution/rename/Sphyraena_novaehollandiae.fasta
Rename_gff_braker.pl GCA_050436865.1_fSphNov1.hap2.cur.20230825_genomic.fna 24 braker.gff3 > ~/jlkang/Carangoidei_evolution/rename/Sphyraena_novaehollandiae.gff3

# Sphyraena barracuda: 24 chromosomes
# (base) shichuang001@login01 Sun Nov 01 2026 18:13:16 ~/jlkang/Carangoidei_evolution/extra_genomes/GCA_988282095.1
Rename.pl GCA_988282095.1_Sphyraena_barracuda_haplotype_1_assembly_genomic.fna 24 > ~/jlkang/Carangoidei_evolution/rename/Sphyraena_barracuda.fasta
Rename_gff_braker.pl GCA_988282095.1_Sphyraena_barracuda_haplotype_1_assembly_genomic.fna 24 braker.gff3 > ~/jlkang/Carangoidei_evolution/rename/Sphyraena_barracuda.gff3

# Lates_calcarifer: 24 chromosomes
# (base) shichuang001@login01 Sun Nov 01 2026 18:15:50 ~/jlkang/Carangoidei_evolution/extra_genomes/GCF_001640805.2
Rename.pl Lates_calcarifer.fasta 24 > ~/jlkang/Carangoidei_evolution/rename/Lates_calcarifer.fasta
Rename_gff.pl Lates_calcarifer.fasta 24 Lates_calcarifer.gff3 > ~/jlkang/Carangoidei_evolution/rename/Lates_calcarifer.gff3
```

**自己测序的基因组**
```bash
# Alectis_indicus: 24 chromosomes, header里面多了"AI_"
# (base) shichuang001@login01 Sun Nov 01 2026 19:04:57 ~/jlkang/Carangoidei_evolution
Rename.pl Alectis_indicus.Genome.V1.genome.fasta 24 > rename/Alectis_indicus.fasta
./Rename_gff.pl Alectis_indicus.Genome.V1.genome.fasta 24 Alectis_indicus.Genome.V1.gff3 > rename/Alectis_indicus.gff3

# Alepes_djedaba: 26 chromosomes
cp Alepes_djedaba_assembly_final.fasta rename/Alepes_djedaba.fasta
cp Alepes_djedaba.Genome.V1.gff3 rename/Alepes_djedaba.gff3

# Caranx_equula: 24 chromosomes
Rename.pl Caranx_equula_assembly_final.fasta 24 > rename/Caranx_equula.fasta
Rename_gff.pl Caranx_equula_assembly_final.fasta 24 Caranx_equula.Genome.V1.gff3 > rename/Caranx_equula.gff3

# Caranx_ignobilis: 24 chromosomes
cp Caranx_ignobilis.Genome.fasta rename/Caranx_ignobilis.fasta
cp Caranx_ignobilis.Genome.V1.gff3 rename/Caranx_ignobilis.gff3

# Caranx_melampygus: 不是chromosome level

# Caranx_sexfasciatus: 24 chromosomes
Rename.pl Caranx_sexfasciatus.Genome.fasta 24 > rename/Caranx_sexfasciatus.fasta
Rename_gff.pl Caranx_sexfasciatus.Genome.fasta 24 Caranx_sexfasciatus.Genome.V1.gff3 > rename/Caranx_sexfasciatus.gff3

# Decapterus_macrosoma: 23 chromosomes
cp Decapterus_macrosoma_assembly_final.fasta rename/Decapterus_macrosoma.fasta
cp Decapterus_macrosoma.Genome.V1.gff3 rename/Decapterus_macrosoma.gff3

# Echeneis_naucrates: 24 chromosomes
Rename.pl Echeneis_naucrates.fasta 24 > rename/Echeneis_naucrates.fasta
Rename_gff.pl Echeneis_naucrates.fasta 24 Echeneis_naucrates.gff3 > rename/Echeneis_naucrates.gff3

# Elagatis_bipinnulata: 24 chromosomes
cp Elagatis_bipinnulata_assembly_final.fasta rename/Elagatis_bipinnulata.fasta
cp Elagatis_bipinnulata.Genome.V1.gff3 rename/Elagatis_bipinnulata.gff3

# Megalaspis_cordyla: 24 chromosomes
Rename.pl Megalaspis_cordyla_assembly_final.fasta 24 > rename/Megalaspis_cordyla.fasta
Rename_gff.pl Megalaspis_cordyla_assembly_final.fasta 24 Megalaspis_cordyla.Genome.V1.gff3 > rename/Megalaspis_cordyla.gff3

# Parastromateus_niger: 24 chromosomes
Rename.pl Parastromateus_niger.Genome.fasta 24 > rename/Parastromateus_niger.fasta
Rename_gff.pl Parastromateus_niger.Genome.fasta 24 Parastromateus_niger.Genome.V1.gff3 > rename/Parastromateus_niger.gff3

# Scomberoides_lysan: 24 chromosomes
Rename.pl Scomberoides_lysan.Genome.fasta 24 > rename/Scomberoides_lysan.fasta
Rename_gff.pl Scomberoides_lysan.Genome.fasta 24 Scomberoides_lysan.Genome.V1.gff3 > rename/Scomberoides_lysan.gff3

# Selaroides_leptolepis: 30 chromosomes
cp Selaroides_leptolepis_assembly_final.fasta rename/Selaroides_leptolepis.fasta
cp Selaroides_leptolepis.Genome.V1.gff3 rename/Selaroides_leptolepis.gff3

# Trachinotus_anak: T2T, 24 chromosomes
Rename.pl Trachinotus_anak_GCA_977005145.1.fasta 24 > rename/Trachinotus_anak.fasta
Rename_gff.pl Trachinotus_anak_GCA_977005145.1.fasta 24 Trachinotus_anak_GCA_977005145.1.gff3 > rename/Trachinotus_anak.gff3

# Trachinotus_blochii: 24 chromosomes
Rename.pl Trachinotus_blochii_our.fasta 24 > rename/Trachinotus_blochii.fasta
Rename_gff.pl Trachinotus_blochii_our.fasta 24 Trachinotus_blochii_our.gff3 > rename/Trachinotus_blochii.gff3

# Trachinotus_carolinus: 24 chromosomes
Rename.pl Trachinotus_carolinus_GCA_040938265.1_ragtag.fasta 24 > rename/Trachinotus_carolinus.fasta
Rename_gff.pl Trachinotus_carolinus_GCA_040938265.1_ragtag.fasta 24 Trachinotus_carolinus_GCA_040938265.1_ragtag.gff3 > rename/Trachinotus_carolinus.gff3
```
**总共28个鲹形目物种（均有染色体级别基因组）**
We set two fossil calibrations: O. latipes–T. nigroviridis (~96.9–150.9 Ma), D. rerio–G. aculeatus (~149.85–165.2 Ma) and a time for the root (< 700 Ma)
选用四个物种：
Japanese Medaka, Oryzias latipes (GCF_053564925.1);
Fugu, Takifugu rubripes (GCF_901000725.3);
Zebrafish, Danio rerio (GCF_052040795.2);
Stickleback, Gasterosteus aculeatus (GCF_964276395.1) 为outgroup
```bash
# Zebrafish, Danio rerio (GCF_052040795.2)
# (base) shichuang001@login01 Sun Nov 01 2026 22:42:37 ~/jlkang/genome_and_gff/others/GCF_052040795.2
cp GCF_052040795.2_GRCz13ab_genomic.fna GCF_052040795.2_GRCz13ab_genomic.gff ~/jlkang/Carangoidei_evolution/
# (base) shichuang001@login01 Sun Nov 01 2026 22:43:52 ~/jlkang/Carangoidei_evolution
cp GCF_052040795.2_GRCz13ab_genomic.fna rename/Danio_rerio.fasta
cp GCF_052040795.2_GRCz13ab_genomic.gff rename/Danio_rerio.gff3

# Stickleback, Gasterosteus aculeatus (GCF_964276395.1)
# (base) shichuang001@login01 Sun Nov 01 2026 22:28:13 ~/jlkang/genome_and_gff/others/GCF_053564925.1
cp GCF_053564925.1_Olatipes_Hd-rR_3.1_genomic.fna GCF_053564925.1_Olatipes_Hd-rR_3.1_genomic.gff ~/jlkang/Carangoidei_evolution/
# (base) shichuang001@login01 Sun Nov 01 2026 22:50:34 ~/jlkang/Carangoidei_evolution
cp GCF_053564925.1_Olatipes_Hd-rR_3.1_genomic.fna rename/Gasterosteus_aculeatus.fasta
cp GCF_053564925.1_Olatipes_Hd-rR_3.1_genomic.gff rename/Gasterosteus_aculeatus.gff3

# Fugu, Takifugu rubripes (GCF_901000725.3)
# (base) shichuang001@login01 Sun Nov 01 2026 22:53:44 ~/jlkang/genome_and_gff/others/GCF_901000725.3
cp GCF_901000725.3_fTakRub1.3_genomic.fna ~/jlkang/Carangoidei_evolution/rename/Takifugu_rubripes.fasta
cp GCF_901000725.3_fTakRub1.3_genomic.gff ~/jlkang/Carangoidei_evolution/rename/Takifugu_rubripes.gff3

# Japanese Medaka, Oryzias latipes (GCF_053564925.1)
# (base) shichuang001@login01 Sun Nov 01 2026 22:58:44 ~/jlkang/genome_and_gff/others/GCF_053564925.1
cp GCF_053564925.1_Olatipes_Hd-rR_3.1_genomic.fna ~/jlkang/Carangoidei_evolution/rename/Oryzias_latipes.fasta
cp GCF_053564925.1_Olatipes_Hd-rR_3.1_genomic.gff ~/jlkang/Carangoidei_evolution/rename/Oryzias_latipes.gff3
```
