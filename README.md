# ggComp
ggComp (Genomic-based Germplasm Compare) is a novel strategy fitting complex genomes to evaluate germplasm resources by identifying shared genomic regions and excluding pervasive copy number variations for pairwise accessions. 

## Installation
ggComp is a light weight software that only need to make sure Python 3.8.5, samtools 1.4 and bcftools 1.9 are well installed and set in environment path.

Clone the ggComp repository
```
git clone https://github.com/zack-young/ggComp.git
```


## Usage

```

Program: ggComp (A pairwised comparison method to identify similar genetic regions and shared CNV regions between accessions.)
Version: 1.0

Usage:   ggComp [-v|--version] [-h|--help] <command> <argument>

Commands:

    CNV_detector            detect CNV regions

    SNP_extractor           extract SNP information from VCF

    DSR_counter             compute the DSR

    SGR_PHR_definer         identify SGR and PHR

    HMM_smoother            smoothing SGR and PHR result using HMM
```

### CNV_detector
`detect CNV regions`

```
Usage: ggComp CNV_detector <--chr_lis <STRING>>  <--single_CNV <FILE> | --pair_CNV <FILE>>

    --chr_lis <STRING>      chromosomes that listed in behind file 
                            e.g. 'chr1A chr1B chr1D'

Options:

    --single_CNV            detect CNV region of samples ::: input file contains chromosme
                            list, path of BED file, path of BAM file and output file
                            e.g. chromosome BED BAM output_dev
                            produce three files: *.mask_CNV
                                                 *.mask_CNV_deletion
                                                 *.mask_CNV_duplication
 
    --pair_CNV              identify pairwise sample specific CNV regions and shared CNV
                            regions::: input file contains sample names, path of CNV file
                            and path of output file directory
                            e.g. sample1_name sample1_directory sample2_name sample2_directory pair_dircetory

chromosome length in bed file shall not exceed the limit of bedtools (400Mb)
see test/single_CNV.config and test/pair_CNV.config for more details
(colunms separated by tab)
```

### SNP_extractor
`extract SNP information from VCF for DSR calculation`

```
Usage: ggComp SNP_extractor <--config <FILE>>

Options:"

    --config <FILE>         file contains path of vcf files, sample ID list in vcf for
                            extracting(separated by ','), path of output file
                            e.g. VCF_file    vcf_ID   output_file

see test/SNP_extractor.config for more details
(colunms separated by tab)
```

### DSR_counter
`compute the DSR (Different SNP ratio) by bin`

```
Usage: ggComp DSR_counter <--config <file>>  [Options]

    --config <FILE>         file contains：path of files contain SNP inforamtion;
                            path of output files; end position of chromosomes
                            (use 'defalut' to call built-in end position of each chromosomes)
                            e.g. SNP_file output_file END\default

Options:

    --DP_low <INT>          exclude SNP with DP <= <int>. default 3

    --DP_high <INT>         exclude SNP with DP >= <int>. default 99

    --GQ_sample <INT>       exclude SNP with GQ <= <int>. default 8

    --bin_size <INT>        bln size. default 1000000

see test/DSR_counter.config for more details (colunms separated by tab)
```

### SGR_PHR_definer
`detect SGR (Similar Genetic Regions) and PHR (Polymorphism Hotspot Regions) between sample pairs`

```
Usage: ggComp SGR_PHR_definer <--noCNV <file> | --plus_CNV <file>> [Options]

    --no_CNV <FILE>         ignore CNV
                            file contains path of files contain DSR information
                            and path of output files.
                            e.g. DSR_file  output

    --plus_CNV <FILE>       take CNV into consideration"  
                            file contains path of files contain DSR information
                            path of CNV file and path of output files.
                            e.g. DSR_file  CNV_file output

Options

    --LEVEL <INT>           threshold divide SGR and PHR. default 10

see test/SGR_PHR_definer.config for more details (colunms separated by tab)
```

### HMM_smoother
`smoothing SGR (Similar Genetic Regions) and PHR (Polymorphism Hotspot Regions) result using HMM`

```
Usage: ggComp HMM_smoother [Options]

Options:

    --input <FOLDER>        folder path of SGR and PHR phasing results

    --output <FOLDER>       output path

    --folder_lis <FILE>     containing a list of folder names of SGR and PHR phasing
                            results that are goining to be smoothed
                            e.g. C2_C10
                                C2_C10

    --processes <INT>       (optional) maximum worker processes"
                            default: 2"

    --train                 (optional) train the model using input data (otherwise using"
                            the model published in the article)"

    --niter <INT>           (optional) maximum number of iterations to perform in training"
                            default: 60"
```


## Quickstart with an example
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

