#!/bin/bash
#set -euxo pipefail
sample_line=$1
meta_data=$2
DEV_PATH=$3
out_file=$4
if [[ -z "${sample_line}" ]]; then
    echo "No first sample provided"; exit 1
fi
if [[ -z "${meta_data}" ]]; then
    echo "No sample data provided."; exit 1
fi

#SAMPLE_list=${SAMPLE1}
#arr=(`echo ${SAMPLE_list} | tr ',' ' '`)

declare -A arr_3
while read line;do
arr_3[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $3}'`
done < ${meta_data}

line_num=1
SAMPLE1=`echo $sample_line | awk -F ',' '{print $1}'`
SAMPLE2=`echo $sample_line | awk -F ',' '{print $2}'`
#if [ "${SAMPLE1}" != "${SAMPLE2}" ]; then
numerator1=`cat ${DEV_PATH}/${SAMPLE1}_${SAMPLE2}/*.homo_undefined_snp_level|grep 'low'|wc -l`
numerator2=`cat ${DEV_PATH}/${SAMPLE1}_${SAMPLE2}/*.homo_undefined_snp_level|grep 'both_CNV'|wc -l`
denominator=`cat ${DEV_PATH}/${SAMPLE1}_${SAMPLE2}/*.homo_undefined_snp_level | wc -l `
decimal=`awk 'BEGIN{printf "%.4f\n",('$numerator1'+'$numerator2')/'$denominator'}'`
#final_ratio=`echo $decimal|awk '{if ($1>=0.5){ratio=1} else if ($1 <0.5&&$1>=0.25){ratio=0.5} else if ($1<0.25&&$1>=0.125){ratio=0.25} else if ($1>=0.0625&&$1<0.125){ratio=0.125} else {ratio=0.0625};print ratio}'`
#final_ratio=`echo $decimal|awk '{if ($1>=0.5){ratio=0.5} else if ($1 <0.5&&$1>=0.4){ratio=0.4} else if ($1<0.4&&$1>=0.2){ratio=0.2};print ratio}'`
#echo -e "${arr_3[${SAMPLE1}]}\t${arr_3[${SAMPLE2}]}\t${final_ratio}" >> All_ratio_torestdis
echo -e "${arr_3[${SAMPLE1}]}\t${arr_3[${SAMPLE2}]}\t${decimal}" >> All_torestdis
#(cp /data/user/yangzz/mapping/fieldergenomecompare/sample_compare/${SAMPLE1}_${SAMPLE2}/${chr}.only_homo_snp_level > dev/${chr}_${SAMPLE1}_${SAMPLE2}_homo_level) &
             #done
#echo $line_num
#line_num=$(( $line_num + 1 ))
wait
#for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#cat  dev/${chr}_*_homo_level |grep 'low'> ${chr}_combine_homo_low
#done
#cat `sed '1d' metadata_cultivar.txt | awk '{print $1"/combine_mask_CNV_deletion_split"}' |tr '\n' ' '|sed s/,$//g`> combine_mask_CNV_deletion_split
#for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
        #num=${chr:3:2}
        #./count_CNV.py --count on --pop_num 7021 -i ${chr}_combine_homo_low |sort -nk1,1 >  ${chr}_combine_homo_low_count
        #./count_CNV.py --compensent on --chrom ${num} -i ${chr}_combine_homo_low_count -o ${chr}.combine_homo_low_compensent
#done

