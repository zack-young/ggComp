#!/bin/bash
set -euxo pipefail

while getopts 'l:b:' opt
do
    case "$opt" in
    l) meta_data=$OPTARG ;;
    b) SAMPLE1_BAM=$OPTARG ;;
    esac
done
if [[ -z "${meta_data}" ]]; then
    echo "No meta_data provided"; exit 1
fi
if [[ -z "${SAMPLE1_BAM}" ]]; then
    echo "No BAM file provided"; exit 1
fi
arr=(`awk '{print $1}' ${meta_data}|tr '\n' ' '|sed s/' '$//g`)
declare -A arr_3
while read line;do
arr_3[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $3}'`
done < ${meta_data}


for SAMPLE1 in ${arr[@]};do
    (if [[ ! -d "${SAMPLE1}" ]]; then
         mkdir ${SAMPLE1}
     fi
         for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do 
         bedtools coverage -a bin/bed/${chr}.1.1M.bed -b ${SAMPLE1_BAM}/${chr}.1.*.bam -counts -sorted > ${SAMPLE1}/${chr}.1.1M.DP
         bedtools coverage -a bin/bed/${chr}.2.1M.bed -b ${SAMPLE1_BAM}/${chr}.2.*.bam -counts -sorted > ${SAMPLE1}/${chr}.2.1M.DP
         gawk '{print $4*1000/($3-$2)}' ${SAMPLE1}/${chr}.1.1M.DP > ${SAMPLE1}/${chr}.1.1M.tmp 
         gawk '{print $4*1000/($3-$2)}' ${SAMPLE1}/${chr}.2.1M.DP > ${SAMPLE1}/${chr}.2.1M.tmp 
         done 
         cat ${SAMPLE1}/chr*.1M.tmp > ${SAMPLE1}/combine_1M_DP
         sh script_1M_un1k.sh ${SAMPLE1}
    ) &
wait_all
done 
wait
#-----------------------------------
#-----------------------------------
for SAMPLE1 in ${arr[@]};do
      for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
     #1
      muti_func_snp_compare.py -b 1000000 -i ${SAMPLE1}/${chr}.1M.norm_un1k --mask_cnv on --sample1 ${arr_3[${SAMPLE1}]}  -o ${SAMPLE1}/${chr}.mask_CNV &
     #2
      done
      wait 
      for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
      for item in {deletion,duplication};do
         grep "${item}" ${SAMPLE1}/${chr}.mask_CNV > ${SAMPLE1}/${chr}.mask_CNV_${item} &
      done
      done
wait_all
done
wait
#----------------------------------
#----------------------------------
parallel -j procfile sh bisample_compare_phase1_2.sh ::: $(eval cat sample_path.txt) ::: 2020-04-17162314
