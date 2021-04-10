#!/bin/bash
#set -euxo pipefail
line=$1
time_date=$2
BED=$3
if [[ -z "${line}" ]]; then
    echo "No line provided"; exit 1
fi
if [[ -z "${time_date}" ]]; then
    echo "No time date provided."; exit 1
fi

#SAMPLE_list=${SAMPLE1}
#arr=(`echo ${SAMPLE_list} | tr ',' ' '`)
path="/data/user/yangzz/mapping/fieldergenomecompare/20200418_201_CNV"
#declare -A arr_2
#while read line;do
#arr_2[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $1}'`
#done < ${meta_data}

#declare -A arr_3
#while read line;do
#arr_3[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $4}'`
#done < ${meta_data}

#for sample1 in ${arr[@]};do
#for sample2 in ${arr_2[@]};do
#    unset arr_2[${sample1}]
#    if [ "${sample1}" != "${sample2}" ]; then
SAMPLE1=`echo $line|awk -F ',' '{print $1}'`
SAMPLE2=`echo $line|awk -F ',' '{print $2}'`

if [[ ! -d "dev/${SAMPLE1}_${SAMPLE2}" ]]; then
     mkdir dev/${SAMPLE1}_${SAMPLE2}
fi 
for chr in `cut  -f1 ${BED}|sort | uniq`; do
    for item in {deletion,duplication};do
      a=`du -h ${SAMPLE1}_${time_date}/${chr}.mask_CNV_${item}`
      b=`echo $a |awk '{print $1}'`
      if [ "${b}" == 0 ];then
        #echo 'file zero' 
        awk '{print $1"\t"$2"\t"$3"\t"SAMPLE2"'"$item"'""_CNV"}' SAMPLE2="$SAMPLE2" ${SAMPLE2}_${time_date}/${chr}.mask_CNV_${item}  > dev/${SAMPLE1}_${SAMPLE2}/${chr}_${SAMPLE2}to${SAMPLE1}_own_CNV_${item}
      else
        #echo 'file not zero' 
        awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==0) print $1"\t"$2"\t"$3"\t"SAMPLE2"'"$item"'""_CNV"}' SAMPLE2="$SAMPLE2" ${SAMPLE1}_${time_date}/${chr}.mask_CNV_${item} ${SAMPLE2}_${time_date}/${chr}.mask_CNV_${item} > dev/${SAMPLE1}_${SAMPLE2}/${chr}_${SAMPLE2}to${SAMPLE1}_own_CNV_${item}
      fi
      a=`du -h ${SAMPLE2}_${time_date}/${chr}.mask_CNV_${item}`
      b=`echo $a |awk '{print $1}'`
      #echo ${b}
      if [ "${b}" == 0 ];then
         #echo 'file zero'
         awk '{print $1"\t"$2"\t"$3"\t"SAMPLE1"'"$item"'""_CNV"}' SAMPLE1="$SAMPLE1" ${SAMPLE1}_${time_date}/${chr}.mask_CNV_${item}  > dev/${SAMPLE1}_${SAMPLE2}/${chr}_${SAMPLE1}to${SAMPLE2}_own_CNV_${item}
      else
         #echo 'file not zero'
         awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==0) print $1"\t"$2"\t"$3"\t"SAMPLE1"'"$item"'""_CNV"}' SAMPLE1="$SAMPLE1" ${SAMPLE2}_${time_date}/${chr}.mask_CNV_${item} ${SAMPLE1}_${time_date}/${chr}.mask_CNV_${item} > dev/${SAMPLE1}_${SAMPLE2}/${chr}_${SAMPLE1}to${SAMPLE2}_own_CNV_${item}
      fi
      awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==1) print $1"\t"$2"\t"$3"\t""'"$item"'""_both_CNV"}' ${SAMPLE1}_${time_date}/${chr}.mask_CNV_${item} ${SAMPLE2}_${time_date}/${chr}.mask_CNV_${item} > dev/${SAMPLE1}_${SAMPLE2}/${chr}_${SAMPLE1}to${SAMPLE2}both_CNV_${item}
 #3  
      cat dev/${SAMPLE1}_${SAMPLE2}/${chr}_${SAMPLE1}to${SAMPLE2}_own_CNV_${item} dev/${SAMPLE1}_${SAMPLE2}/${chr}_${SAMPLE2}to${SAMPLE1}_own_CNV_${item} dev/${SAMPLE1}_${SAMPLE2}/${chr}_${SAMPLE1}to${SAMPLE2}both_CNV_${item} > dev/${SAMPLE1}_${SAMPLE2}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV_${item}
    done
    cat dev/${SAMPLE1}_${SAMPLE2}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV_deletion dev/${SAMPLE1}_${SAMPLE2}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV_duplication | sort -u -nk2,2 > dev/${SAMPLE1}_${SAMPLE2}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV
done
#
#wait_all
wait
#
#for sample1 in ${arr[@]};do
#for sample2 in ${arr_2[@]};do
#    unset arr_2[${sample1}]
#    if [ "${sample1}" != "${sample2}" ]; then
#    for chr in chr1A chr2A chr3A chr4A chr5A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do

    #cat tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.homo_snp_density_* > tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.homo_snp_density
    #./muti_func_snp_compare.py --two_diff_level on -i tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.homo_snp_density --sample1 ${sample1} --sample2 ${sample2} --time_data ${time_date} --total ${total} &
#    a1=`du -h tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.${sample1}to${sample2}all_CNV`
#    b1=`echo ${a1} |awk '{print $1}'`
    #cat tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.only_homo_snp_level tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.${sample1}to${sample2}all_CNV > tmp/${total}/${sample1}_${sample2}_${time_date}/${chr}.homo_snp_level
    #if [ "${b}" != 0 ];then
    #awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==0) print $1"\t"$2"\t"$3"\t"SAMPLE1"'"$item"'""_CNV"}' SAMPLE1="$SAMPLE1" tmp/${SAMPLE2}_${time_date}/chr${num}${i}.mask_CNV_${item} tmp/${SAMPLE1}_${time_date}/chr${num}${i}.mask_CNV_${item} > tmp/${SAMPLE1}_${SAMPLE2}_${time_date}/chr${num}${i}_${SAMPLE1}to${SAMPLE2}_own_CNV_${item}
    #fi
   #awk '{print $4}' tmp/${sample}_${SAMPLE2}_${time_date}/combine_homo_snp_density > tmp/${sample}_${SAMPLE2}_${time_date}/combine_homodiff_snp_density
   #sed -i '1i '${sample}_${SAMPLE2}'' tmp/${sample}_${SAMPLE2}_${time_date}/combine_homodiff_snp_density
#    done
#    ./ptf_maker.py -p "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/${total}/${sample1}_${sample2}_${time_date}" -s ".homo_snp_level" > plotfile/plotfile_${sample1}_${sample2}
#    Rscript ~/R/graph_diff_sim.R --sample1 ${sample1} --sample1_name ${arr_3[${sample1}]} --sample2 ${sample2} --sample2_name ${arr_3[${sample2}]} -s ${time_date}
#    fi
#paste `awk '{if ($1!="'"${sample}"'")print "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/""'"${sample}"'""_"$1"_2019-12-15201714/combine_homodiff_snp_density"}' field_cultivar/metadata_cultivar_noheader.txt |tr '\n' ' '|sed s/,$//g` > field_cultivar/combine_${sample}_${SAMPLE2}_homodiff_snp_density
#done
#cat `awk '{if ($1!="PH09")print "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/PH09_"$1"_2019-12-15201714/combine_homo_snp_density"}' field_cultivar/metadata_cultivar_noheader.txt |tr '\n' ' '|sed s/,$//g` > field_cultivar/combine_PH09_homo_snp_density
#done




#for sample in ${arr[@]};do
#while read line;do
#    SAMPLE2=`echo $line | awk '{print $1}'`
#    if [ "${sample}" != "${SAMPLE2}" ]; then
#    if [[ -s "tmp/${sample}_${SAMPLE2}_${time_date}/chr1A.homo_snp_level_20" ]]; then
#    cat tmp/${sample}_${SAMPLE2}_${time_date}/chr*.homo_snp_level_* > tmp/${sample}_${SAMPLE2}_${time_date}/combine_homo_snp_level
#    awk '{print $5}' tmp/${sample}_${SAMPLE2}_${time_date}/combine_homo_snp_level > tmp/${sample}_${SAMPLE2}_${time_date}/combine_homodiff_snp_level
#    sed -i '1i '${sample}_${SAMPLE2}'' tmp/${sample}_${SAMPLE2}_${time_date}/combine_homodiff_snp_level
#    cp tmp/${sample}_${SAMPLE2}_${time_date}/combine_homo_snp_level field_cultivar/tmp_snp_level/${sample}_${SAMPLE2}_homo_snp_level
#    cp tmp/${sample}_${SAMPLE2}_${time_date}/combine_homodiff_snp_level field_cultivar/tmp_snp_level/${sample}_${SAMPLE2}_homodiff_snp_level
#    fi
#    fi
#done < ${meta_data}
#paste `awk '{if ($1!="'"${sample}"'")print "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/""'"${sample}"'""_"$1"_2019-12-15201714/combine_homodiff_snp_density"}' field_cultivar/metadata_cultivar_noheader.txt |tr '\n' ' '|sed s/,$//g` > field_cultivar/combine_${sample}_${SAMPLE2}_homodiff_snp_density
#done
#cat `awk '{if ($1!="PH09")print "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/PH09_"$1"_2019-12-15201714/combine_homo_snp_density"}' field_cultivar/    metadata_cultivar_noheader.txt |tr '\n' ' '|sed s/,$//g` > field_cultivar/combine_PH09_homo_snp_density

