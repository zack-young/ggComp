#!/bin/bash
#set -euxo pipefail
line=$1
chr=$2
meta_data=$3
SAMPLE1=`echo $line|awk -F ',' '{print $1}'`
SAMPLE2=`echo $line|awk -F ',' '{print $2}'`
declare -A arr_3
while read line;do
arr_3[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $3}'`
done < ${meta_data}
    #./muti_func_snp_compare.py --two_diff_level on -i ${SAMPLE1}_${SAMPLE2}/${chr}.1M.density --sample1 ${SAMPLE1} --sample2 ${SAMPLE2} --chromosome ${chr}
    #line=`cat dev/${SAMPLE1}_${SAMPLE2}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV|wc -l`
    #if test $line = 0;then
    #cat ${SAMPLE1}_${SAMPLE2}/${chr}.only_homo_snp_undefined_level  > ${SAMPLE1}_${SAMPLE2}/${chr}.homo_undefined_snp_level
    #else
    #type=`cut -f4 dev/${SAMPLE1}_${SAMPLE2}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV`
    #awk '{print $1"\t"$2"\t"$3"\t"type}' type="$type" dev/${SAMPLE1}_${SAMPLE2}/${chr}.only_homo_snp_undefined_level > dev/${SAMPLE1}_${SAMPLE2}/${chr}.homo_undefined_snp_level
    #awk -F '\t' 'NR==FNR{a[$1]=$0;next}NR>FNR{if ($1 in a==1) {print a[$1]} else {print $0}}' dev/${SAMPLE1}_${SAMPLE2}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV dev/${SAMPLE1}_${SAMPLE2}/${chr}.only_homo_snp_undefined_level > dev/${SAMPLE1}_${SAMPLE2}/${chr}.homo_undefined_snp_level
    #fi
         #wait
    ./ptf_maker.py -p "/data//user/yangzz/mapping/fieldergenomecompare/BS_300/${SAMPLE1}_${SAMPLE2}" -s ".homo_undefined_snp_level"> plotfile/plotfile_${SAMPLE1}_${SAMPLE2}
    Rscript ~/R/graph_diff_sim.R -d "~/mapping/fieldergenomecompare/BS_300" --sample1 ${SAMPLE1} --sample1_name ${arr_3[${SAMPLE1}]} --sample2 ${SAMPLE2} --sample2_name ${arr_3[${SAMPLE2}]} -s "" 
         #sleep 0.00000000001s
#wait_all
