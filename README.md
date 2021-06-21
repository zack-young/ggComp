# ggComp
ggComp (Genomic-based Germplasm Compare) is a novel strategy fitting complex genomes to evaluate germplasm resources by identifying shared genomic regions and excluding pervasive copy number variations for pairwise accessions. 
Performing copy number analysis from targeted capture high-throughput sequencing has been a challenging task. This involves binning the targeted region, calculating the log ratio of the read depths between the sample and the reference, and then stitching together thousands of these data points into numerous segments (especially in the context of cancer) to derive the copy number state of genomic regions. Recently, several tools have been developed to adequately detect both somatic as well as germline CNVs. However, review and interpretation of these variants in a clinical as well as research setting is a daunting task. This can involve frequent switches back and forth from a static image to numerous tabular files resulting in an exasperated reviewer.  
 ReconCNV has been developed to overcome this challenge by providing interactive dashboard for hunting copy number variations (CNVs) from high-throughput sequencing data. The tool has been tested for targeted gene panels (including exome data). Python3's powerful visualization and data manipulation modules, namely Bokeh and Pandas, are utilized to create these dynamic visualizations. ReconCNV can be readily applied to most CNV calling algorithms with simple modifications to the configuration file. 
## Installation
The easiest way to get started with reconCNV is via `conda`. Using `conda` ensures you are running Python 3.6 (using which reconCNV was coded) and all dependencies are installed within an virtual environment. This avoids dependency conflicts with existing programs. If `conda` is not available on your system, a minimal installer can be installed using [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

light weight
1. Clone the reconCNV repository
    ```
    git clone https://github.com/zack-young/ggComp.git
    ```
2. check environment bedtools bcftools

3. You are now ready to use reconCNV!

Usage example:


### CNV_detector            
detect CNV regions

```
sh src/ggComp.sh CNV_detector --chr_lis chr1A --single_CNV test/single_CNV.config
```

```
sh src/ggComp.sh CNV_detector --chr_lis chr1A --pair_CNV test/pair_CNV.config
```

### SNP_extractor
extract SNP information from VCF

```
sh src/ggComp.sh SNP_extractor --config test/SNP_extractor.config
```

### DSR_counter
compute the DSR

```
sh src/ggComp.sh DSR_counter --config test/DSR_counter.config
```

### SGR_PHR_definer
identify SGR and PHR

```
sh src/ggComp.sh SGR_PHR_definer --plus_CNV test/SGR_PHR_plus_CNV.config
```

```
sh src/ggComp.sh SGR_PHR_definer --no_CNV test/SGR_PHR_noCNV.config
```

### HMM_smoother
smoothing SGR and PHR result using HMM
#### Smooth only
```
sh WheatComp.sh HMM_smoother \
    -i /data3/user3/wangwx/projs/HMM_for_yzz_comp/210328-modify/test/data \
    --folder_lis /data3/user3/wangwx/projs/HMM_for_yzz_comp/210328-modify/test/folders.txt \
    -o /data3/user3/wangwx/projs/HMM_for_yzz_comp/210328-modify/test/out\
    --processes 21
```
#### Train & Smooth
```
sh WheatComp.sh HMM_smoother \
    -i /data3/user3/wangwx/projs/HMM_for_yzz_comp/210328-modify/test/data \
    --folder_lis /data3/user3/wangwx/projs/HMM_for_yzz_comp/210328-modify/test/folders.txt \
    -o /data3/user3/wangwx/projs/HMM_for_yzz_comp/210328-modify/test/out \
    --processes 21 \
    --train \
    --niter 30
```
-i: 数据所处文件夹路径
-f：样本比较文件夹名称列表文件
-o：输出目录路径
-p：[Optional] 最大并发数，默认2
-t：[Optional] 置位则重新fit模型
-n：[Optional] 训练模型时最大迭代次数。如果不输入此项则默认最多迭代60次
