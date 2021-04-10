#!/bin/bash
set -euxo pipefail
sample=$1
su=$2
chr=$3
START=$4
suffix=$5
dev_path=$6
out_dir=$7
meta_data=$8

suffix2=${suffix}_${su}

declare -A arr_3
while read line;do
arr_3[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $3}'`
done < ${meta_data}
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
while read line;do
           #myjob=`jobs|wc -l`
           SAMPLE1=`echo $line | awk -F ',' '{print $1}'`
           SAMPLE2=`echo $line | awk -F ',' '{print $2}'`
           num_low=0
           num_diff=0
           num_cnv=10
           #echo $num_tmp
           item=`awk '{if($2==num) print $4}' num="$START" ${dev_path}/${SAMPLE1}_${SAMPLE2}/${chr}.${suffix}`
           if [[ "$item" = "low" ]];then
             num_low=$(($num_low + 1))
           elif [[ "$item" =~ (.*)both_CNV ]];then
             num_cnv=$(($num_cnv + 10))
           else
             num_diff=$(($num_diff + 1))
           fi
           if [ "$num_diff" = 0 ]&&[ "$num_low" = 0 ] ;then
             final_num=$num_cnv
           else
             final_num=`awk 'BEGIN{printf "%.2f\n",'$num_cnv'+'$num_low'/('$num_low'+'$num_diff')}'`
           fi
           if [ "$final_num" != '10.00' ]; then
             echo -e "${arr_3[${SAMPLE1}]}\t${arr_3[${SAMPLE2}]}\t${final_num}" >> ${out_dir}/${chr}_${suffix2}
           fi
           #if [[ $myjob > 70 ]];then
           #   wait_all
           #fi
line_num=$(( $line_num + 1 ))
done < ${sample}
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

