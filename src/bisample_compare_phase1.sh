#!/bin/bash
while getopts 'c:l:' opt
do
    case "$opt" in
    c) config=$OPTARG ;;
    l) chr_lis=$OPTARG ;;
    esac
done

path=$(cd "$(dirname "$0")";pwd)

if ! [ -x "$(command -v bedtools)" ]; then
    echo 'Error: bedtools is not installed or put in PATH.' >&2
    exit 1
fi
while read lines;do
    CHR=`echo $lines|awk '{print $1}'`
    BED=`echo $lines|awk '{print $2}'`
    BAM=`echo $lines|awk '{print $3}'`
    DIR=`echo $lines|awk '{print $4}'`
    if [[ ! -d "${DIR}" ]]; then
        mkdir ${DIR}
    fi
    bedtools coverage -a ${BED} -b ${BAM} -counts -sorted > ${DIR}/${CHR}.DP
    gawk '{print $1"\t"$2"\t"$3"\t"$4*1000/($3-$2)}' ${DIR}/${CHR}.DP > ${DIR}/${CHR}.tmp 
done < ${config}

for DIR in `cut -f4 ${config}|uniq`;do
    cat ${DIR}/*.tmp > ${DIR}/combine_1M_DP
    ave=`python ${path}/count_CNV.py -i ${DIR}/combine_1M_DP`
    for CHR in `echo $chr_lis`;do
        gawk -vOFS="\t" -vave=$ave '{print $1"\t"$2"\t"$3"\t"$4/ave}' ${DIR}/${CHR}.tmp > ${DIR}/${CHR}.norm
        python ${path}/muti_func_snp_compare.py -i ${DIR}/${CHR}.norm --mask_cnv on -o ${DIR}/${CHR}.mask_CNV
        for item in {deletion,duplication};do
            grep "${item}" ${DIR}/${CHR}.mask_CNV > ${DIR}/${CHR}.mask_CNV_${item}
        done
    done
done
echo "Single CNV detection Complete"
#-----------------------------------
#-----------------------------------