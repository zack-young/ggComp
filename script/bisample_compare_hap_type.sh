#!/bin/bash
#set -euxo pipefail
chr=$1
n=$2
#dev_path=$6
#out_dir=$7

START=$(( ($n-1)*1000000+1 ))

#$declare -A arr_3
#while read line;do
#arr_3[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $3}'`
#done < ${meta_data}
#SAMPLE_list=${SAMPLE1}
#arr=(`echo ${SAMPLE_list} | tr ',' ' '`)
#arr_3=(`awk '{print $1}' ${meta_data}|tr '\n' ' '|sed s/' '$//g`)
#declare -A arr_2
#while read line;do
#arr_2[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $1}'`
#done < ${meta_data}

num_tmp=0
line_num=1
#dir="/data/user/yangzz/mapping/fieldergenomecompare/statistic/Rht8"
if true;then
    a=`du -h 10+_mcl/${chr}.${n}_mcl_undefined|awk '{print $1}'`
    text=`./hap_assign.py -i 10+_mcl/${chr}.${n}_mcl_undefined  -f $a -s ../1.diff_dev/metadata_cultivar_10+.txt |tr '\n' ' '|sed s/' '$//g`
    echo -e "${chr}\t$START\t$text" >> ${chr}_hap_data.txt
           #if [[ $myjob > 70 ]];then
           #   wait_all
           #fi
fi
wait
#paste `awk '{print "'$dir'""/""'$chr'""_"$1"_""'$suffix2'""_dist"}'  ../metadata_cultivar_final.txt |tr '\n' ' '|sed s/,$//g` > ${dir}/${chr}_combine_${suffix2}_dist
##
#for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
#cat  dev/${chr}_*_homo_level |grep 'low'> ${chr}_combine_homo_low
#done
#cat `sed '1d' metadata_cultivar.txt | awk '{print $1"/combine_mask_CNV_deletion_split"}' |tr '\n' ' '|sed s/,$//g`> combine_mask_CNV_deletion_split
#for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
        #num=${chr:3:2}
        #./count_CNV.py --count on --pop_num 7021 -i ${chr}_combine_homo_low |sort -nk1,1 >  ${chr}_combine_homo_low_count
        #./count_CNV.py --compensent on --chrom ${num} -i ${chr}_combine_homo_low_count -o ${chr}.combine_homo_low_compensent
#done

