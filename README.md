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

**共线性jcvi安装**
```bash
# h2076@h2076 Mon Oct 26 2026 16:03:38 ~/sex_change
sudo apt-get install -y texlive texlive-latex-extra texlive-latex-recommended dvipng cm-super msttcorefonts
conda create -n jcvi python=3.9
conda activate jcvi
conda install pip
conda install -c conda-forge biopython
conda install -c conda-forge matplotlib
conda install -c anaconda scipy
conda install -c bioconda last
pip install jcvi

# (jcvi) h2076@h2076 Mon Oct 26 2026 16:03:38 ~/sex_change
```

