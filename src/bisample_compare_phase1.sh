#!/bin/bash
while getopts 'c:l:' opt
do
    case "$opt" in
    c) config=$OPTARG ;;
    l) chr_lis=$OPTARG ;;
    esac
done


while read lines;do
    CHR=`echo $lines|awk '{print $1}'`
    BED=`echo $lines|awk '{print $2}'`
    BAM=`echo $lines|awk '{print $3}'`
    DIR=`echo $lines|awk '{print $4}'`
    bedtools coverage -a ${BED} -b ${BAM} -counts -sorted > ${DIR}/${CHR}.DP
    gawk '{print $1"\t"$2"\t"$3"\t"$4*1000/($3-$2)}' ${DIR}/${CHR}.DP > ${DIR}/${CHR}.tmp 

done < ${config}
cat ${DIR}/*.tmp > ${DIR}/combine_1M_DP
sh script_1M_un1k.sh ${DIR}
ave=`./count_CNV.py -i combine_1M_DP`

for CHR in `echo $chr_lis`;do
    gawk -vOFS="\t" -vave=$ave '{print $1"\t"$2"\t"$3"\t$4/ave}' ${DIR}/${CHR}.tmp > ${DIR}/${CHR}.norm
    muti_func_snp_compare.py -i ${DIR}/${CHR}.norm --mask_cnv on  -o ${DIR}/${CHR}.mask_CNV
done

rm ${DIR}/*tmp
rm ${DIR}/*DP
rm ${DIR}/*norm
#-----------------------------------
#-----------------------------------

for CHR in `echo $chr_lis`;do
    for item in {deletion,duplication};do
        grep "${item}" ${DIR}/${CHR}.mask_CNV > ${DIR}/${CHR}.mask_CNV_${item}
    done
done
#----------------------------------
#----------------------------------
