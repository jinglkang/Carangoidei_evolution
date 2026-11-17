# 基因组共线性
## FG1
```bash
# (base) shichuang001@login01 Sun Nov 15 2026 10:08:08 ~/jlkang/Carangoidei_evolution/rename
mkdir FG1_synteny FG2_synteny FG3_synteny
cp Xiphias_gladius* Toxotes_jaculatrix* Lates_calcarifer* Sphyraena_barracuda* Sphyraena_novaehollandiae* FG1_synteny/
cd FG1_synteny

# 1. 所有的fasta文件都改为60个核苷酸/行
# (base) shichuang001@login01 Sun Nov 15 2026 10:46:04 ~/jlkang/Carangoidei_evolution/rename/FG1_synteny
for i in *.fasta;do seqkit seq -w 60 ${i} > genome_wrapped.fa;mv genome_wrapped.fa ${i};done

# 2. AGAT获得cds序列
# (base) shichuang001@login01 Sun Nov 15 2026 10:49:19 ~/jlkang/Carangoidei_evolution/rename/FG1_synteny
conda activate agat
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

# 3. AGAT获得bed文件
for gff in *.gff3; do
    prefix="${gff%.gff3}"
    agat_convert_sp_gff2bed.pl \
            --gff "$gff" \
            -o "${prefix}.bed"    
done

# 4. reformat fasta文件
conda activate jcvi
for fa in *.fasta; do
	prefix="${fa%.fasta}"
	cdsfa="${prefix}_genome.cds.fa"
	cds="${prefix}.cds"
	python -m jcvi.formats.fasta format \
        $cdsfa \
        $cds
done

# 5. 相互共线性搜索Pairwise synteny search
# Xiphias_gladius vs. Toxotes_jaculatrix
# Toxotes_jaculatrix vs. Lates_calcarifer
# Lates_calcarifer vs. Sphyraena_barracuda
# Sphyraena_barracuda vs. Sphyraena_novaehollandiae
sbatch FG1_synteny.sh
# Submitted batch job 926650

# 6. 准备画图
# 6.1 准备seqids
# (jcvi) shichuang001@login01 Sun Nov 15 2026 13:56:32 ~/jlkang/Carangoidei_evolution/rename/FG1_synteny
perl generate_seqids.pl > seqids
# 6.2 准备layout
vi layout

# 7. 画图
python -m jcvi.graphics.karyotype seqids layout
```

**layout(FG1)**
```
# y, xstart, xend, rotation, color, label, va,  bed
 .9,     .1,    .8,       0,      , Xiphias_gladius, top, Xiphias_gladius.bed
 .7,     .1,    .8,       0,      , Toxotes_jaculatrix, top, Toxotes_jaculatrix.bed
 .5,     .1,    .8,       0,      , Lates_calcarifer, top, Lates_calcarifer.bed
 .3,     .1,    .8,       0,      , Sphyraena_barracuda, top, Sphyraena_barracuda.bed
 .1,     .1,    .8,       0,      , Sphyraena_novaehollandiae, top, Sphyraena_novaehollandiae.bed


# edges
e, 0, 1, Xiphias_gladius.Toxotes_jaculatrix.anchors.simple
e, 1, 2, Toxotes_jaculatrix.Lates_calcarifer.anchors.simple
e, 2, 3, Lates_calcarifer.Sphyraena_barracuda.anchors.simple
e, 3, 4, Sphyraena_barracuda.Sphyraena_novaehollandiae.anchors.simple
```

- **generate_seqids.pl**
```perl
#!/usr/bin/perl
use strict;
use warnings;

my $spe="spe.txt";
open SPE, $spe or die "can not open $spe:$!\n";
while (<SPE>) {
	chomp;
	my @a=split;
	my $fa=$a[0].".fasta";
	my $info;
	open FA, $fa or die "can not open $fa:$!\n";
	while (<FA>) {
		chomp;
		if (/>chr/i) {
			s/>//;
			$info.=$_.",";
		}
	}
	$info=~s/\,$//;
	print "$info\n";
}
```


```FG1_synteny.sh
#!/bin/bash
#SBATCH -J FG1_jcvi # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 20 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.FG1_jcvi # 标准输出文件
#SBATCH -e logs/%j.FG1_jcvi_error # 错误输出文件
echo ${SLURM_JOB_NODELIST} # 作业占用节点列表
echo start on $(date) # 开始时间

# Xiphias_gladius vs. Toxotes_jaculatrix
# Toxotes_jaculatrix vs. Lates_calcarifer
# Lates_calcarifer vs. Sphyraena_barracuda
# Sphyraena_barracuda vs. Sphyraena_novaehollandiae

python -m jcvi.compara.catalog ortholog Xiphias_gladius Toxotes_jaculatrix --no_strip_names
python -m jcvi.compara.catalog ortholog Toxotes_jaculatrix Lates_calcarifer --no_strip_names
python -m jcvi.compara.catalog ortholog Lates_calcarifer Sphyraena_barracuda --no_strip_names
python -m jcvi.compara.catalog ortholog Sphyraena_barracuda Sphyraena_novaehollandiae --no_strip_names
echo end on $(date) # 结束时间

# 5. generates a .simple file from the .anchors file
python -m jcvi.compara.synteny screen --minspan=30 --simple Xiphias_gladius.Toxotes_jaculatrix.anchors Xiphias_gladius.Toxotes_jaculatrix.anchors.new
python -m jcvi.compara.synteny screen --minspan=30 --simple Toxotes_jaculatrix.Lates_calcarifer.anchors Toxotes_jaculatrix.Lates_calcarifer.anchors.new
python -m jcvi.compara.synteny screen --minspan=30 --simple Lates_calcarifer.Sphyraena_barracuda.anchors Lates_calcarifer.Sphyraena_barracuda.anchors.new
python -m jcvi.compara.synteny screen --minspan=30 --simple Sphyraena_barracuda.Sphyraena_novaehollandiae.anchors Sphyraena_barracuda.Sphyraena_novaehollandiae.anchors.new
```

## FG2
```bash
# (base) shichuang001@login01 Sun Nov 15 2026 10:54:32 ~/jlkang/Carangoidei_evolution/rename
cp Echeneis_naucrates* Remorina_albescens* Scomberoides_lysan* Trachinotus* FG2_synteny/
cd FG2_synteny/

for i in *.fasta;do seqkit seq -w 60 ${i} > genome_wrapped.fa;mv genome_wrapped.fa ${i};done

conda activate agat
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

for gff in *.gff3; do
    prefix="${gff%.gff3}"
    agat_convert_sp_gff2bed.pl \
            --gff "$gff" \
            -o "${prefix}.bed"    
done

conda activate jcvi
for fa in *.fasta; do
	prefix="${fa%.fasta}"
	cdsfa="${prefix}_genome.cds.fa"
	cds="${prefix}.cds"
	python -m jcvi.formats.fasta format \
        $cdsfa \
        $cds
done

vi pairwise.txt
perl generate_pairwise.pl # 把输出的内容复制到FG2_synteny.sh
sbatch FG2_synteny.sh
# Submitted batch job 926701

perl generate_seqids.pl > seqids
vi layout
python -m jcvi.graphics.karyotype seqids layout
```

- **FG2_synteny.sh**
```FG2_synteny.sh
#!/bin/bash
#SBATCH -J FG2_jcvi # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 20 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.FG2_jcvi # 标准输出文件
#SBATCH -e logs/%j.FG2_jcvi_error # 错误输出文件
echo ${SLURM_JOB_NODELIST} # 作业占用节点列表
echo start on $(date) # 开始时间

python -m jcvi.compara.catalog ortholog Echeneis_naucrates Remorina_albescens --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple Echeneis_naucrates.Remorina_albescens.anchors Echeneis_naucrates.Remorina_albescens.anchors.new
python -m jcvi.compara.catalog ortholog Remorina_albescens Scomberoides_lysan --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple Remorina_albescens.Scomberoides_lysan.anchors Remorina_albescens.Scomberoides_lysan.anchors.new
python -m jcvi.compara.catalog ortholog Scomberoides_lysan Trachinotus_ovatus --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple Scomberoides_lysan.Trachinotus_ovatus.anchors Scomberoides_lysan.Trachinotus_ovatus.anchors.new
python -m jcvi.compara.catalog ortholog Trachinotus_ovatus Trachinotus_anak --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple Trachinotus_ovatus.Trachinotus_anak.anchors Trachinotus_ovatus.Trachinotus_anak.anchors.new
python -m jcvi.compara.catalog ortholog Trachinotus_anak Trachinotus_carolinus --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple Trachinotus_anak.Trachinotus_carolinus.anchors Trachinotus_anak.Trachinotus_carolinus.anchors.new
python -m jcvi.compara.catalog ortholog Trachinotus_carolinus Trachinotus_blochii --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple Trachinotus_carolinus.Trachinotus_blochii.anchors Trachinotus_carolinus.Trachinotus_blochii.anchors.new
```

- **generate_pairwise.pl**
```perl
#!/usr/bin/perl
use strict;
use warnings;

# python -m jcvi.compara.catalog ortholog Xiphias_gladius Toxotes_jaculatrix --no_strip_names
# python -m jcvi.compara.synteny screen --minspan=30 --simple Xiphias_gladius.Toxotes_jaculatrix.anchors Xiphias_gladius.Toxotes_jaculatrix.anchors.new

my @pairs;
my $pair="pairwise.txt";
open PAIR, $pair or die "can not open $pair:$!\n";
while (<PAIR>) {
	chomp;
	my @a=split;
	my $cmd1="python -m jcvi.compara.catalog ortholog $a[0] $a[1] --no_strip_names";
	my $cmd2="python -m jcvi.compara.synteny screen --minspan=30 --simple $a[0].$a[1].anchors $a[0].$a[1].anchors.new";
	print "$cmd1\n$cmd2\n";
}
```

- **pairwise.txt**
```txt
Echeneis_naucrates Remorina_albescens
Remorina_albescens Scomberoides_lysan
Scomberoides_lysan Trachinotus_ovatus
Trachinotus_ovatus Trachinotus_anak
Trachinotus_anak Trachinotus_carolinus
Trachinotus_carolinus Trachinotus_blochii
```

## FG3
```bash
(base) shichuang001@login01 Sun Nov 15 2026 11:00:09 ~/jlkang/Carangoidei_evolution/rename
cp Elagatis_bipinnulata* Decapterus* Trachurus* Parastromateus_niger* Caranx* Alectis_indicus* Selaroides_leptolepis* Alepes_djedaba* Megalaspis_cordyla* FG3_synteny/
cd FG3_synteny/

for i in *.fasta;do seqkit seq -w 60 ${i} > genome_wrapped.fa;mv genome_wrapped.fa ${i};done

conda activate agat
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

for gff in *.gff3; do
    prefix="${gff%.gff3}"
    agat_convert_sp_gff2bed.pl \
            --gff "$gff" \
            -o "${prefix}.bed"    
done

conda activate jcvi
for fa in *.fasta; do
	prefix="${fa%.fasta}"
	cdsfa="${prefix}_genome.cds.fa"
	cds="${prefix}.cds"
	python -m jcvi.formats.fasta format \
        $cdsfa \
        $cds
done

vi pairwise.txt
perl generate_pairwise.pl
vi FG3_synteny.sh
sbatch FG3_synteny.sh
# Submitted batch job 926702

perl generate_seqids.pl > seqids
vi layout
python -m jcvi.graphics.karyotype seqids layout
```

- **FG3_synteny.sh**
```bash
#!/bin/bash
#SBATCH -J FG3_jcvi # 指定作业名
#SBATCH -p base # 指定队列
#SBATCH -N 1 # 请求节点数
#SBATCH -n 20 # 请求核心数
#SBATCH --gres=gpu:0 # 请求gpu数
#SBATCH -o logs/%j.FG3_jcvi # 标准输出文件
#SBATCH -e logs/%j.FG3_jcvi_error # 错误输出文件
echo ${SLURM_JOB_NODELIST} # 作业占用节点列表
echo start on $(date) # 开始时间

python -m jcvi.compara.catalog ortholog Elagatis_bipinnulata Decapterus_macrosoma --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple Elagatis_bipinnulata.Decapterus_macrosoma.anchors Elagatis_bipinnulata.Decapterus_macrosoma.anchors.new
python -m jcvi.compara.catalog ortholog Decapterus_macrosoma Decapterus_maruadsi --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple Decapterus_macrosoma.Decapterus_maruadsi.anchors Decapterus_macrosoma.Decapterus_maruadsi.anchors.new
python -m jcvi.compara.catalog ortholog Decapterus_maruadsi Trachurus_trachurus --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple Decapterus_maruadsi.Trachurus_trachurus.anchors Decapterus_maruadsi.Trachurus_trachurus.anchors.new
python -m jcvi.compara.catalog ortholog Trachurus_trachurus Trachurus_japonicus --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple Trachurus_trachurus.Trachurus_japonicus.anchors Trachurus_trachurus.Trachurus_japonicus.anchors.new
python -m jcvi.compara.catalog ortholog Trachurus_japonicus Parastromateus_niger --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple Trachurus_japonicus.Parastromateus_niger.anchors Trachurus_japonicus.Parastromateus_niger.anchors.new
python -m jcvi.compara.catalog ortholog Parastromateus_niger Caranx_equula --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple Parastromateus_niger.Caranx_equula.anchors Parastromateus_niger.Caranx_equula.anchors.new
python -m jcvi.compara.catalog ortholog Caranx_equula Alectis_indicus --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple Caranx_equula.Alectis_indicus.anchors Caranx_equula.Alectis_indicus.anchors.new
python -m jcvi.compara.catalog ortholog Alectis_indicus Selaroides_leptolepis --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple Alectis_indicus.Selaroides_leptolepis.anchors Alectis_indicus.Selaroides_leptolepis.anchors.new
python -m jcvi.compara.catalog ortholog Selaroides_leptolepis Alepes_djedaba --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple Selaroides_leptolepis.Alepes_djedaba.anchors Selaroides_leptolepis.Alepes_djedaba.anchors.new
python -m jcvi.compara.catalog ortholog Alepes_djedaba Megalaspis_cordyla --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple Alepes_djedaba.Megalaspis_cordyla.anchors Alepes_djedaba.Megalaspis_cordyla.anchors.new
python -m jcvi.compara.catalog ortholog Megalaspis_cordyla Caranx_crysos --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple Megalaspis_cordyla.Caranx_crysos.anchors Megalaspis_cordyla.Caranx_crysos.anchors.new
python -m jcvi.compara.catalog ortholog Caranx_crysos Caranx_sexfasciatus --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple Caranx_crysos.Caranx_sexfasciatus.anchors Caranx_crysos.Caranx_sexfasciatus.anchors.new
python -m jcvi.compara.catalog ortholog Caranx_sexfasciatus Caranx_ignobilis --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple Caranx_sexfasciatus.Caranx_ignobilis.anchors Caranx_sexfasciatus.Caranx_ignobilis.anchors.new
python -m jcvi.compara.catalog ortholog Caranx_ignobilis Caranx_heberi --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple Caranx_ignobilis.Caranx_heberi.anchors Caranx_ignobilis.Caranx_heberi.anchors.new
```

- **pairwise.txt**
```txt
Elagatis_bipinnulata Decapterus_macrosoma
Decapterus_macrosoma Decapterus_maruadsi
Decapterus_maruadsi Trachurus_trachurus
Trachurus_trachurus Trachurus_japonicus
Trachurus_japonicus Parastromateus_niger
Parastromateus_niger Caranx_equula
Caranx_equula Alectis_indicus
Alectis_indicus Selaroides_leptolepis
Selaroides_leptolepis Alepes_djedaba
Alepes_djedaba Megalaspis_cordyla
Megalaspis_cordyla Caranx_crysos
Caranx_crysos Caranx_sexfasciatus
Caranx_sexfasciatus Caranx_ignobilis
Caranx_ignobilis Caranx_heberi
```
