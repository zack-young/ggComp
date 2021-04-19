#!/bin/bash
#set -euxo pipefail
while getopts 'c:l:' opt
do
    case "$opt" in
    c) config=$OPTARG ;;
    l) chr_lis=$OPTARG ;;
    esac
done

while read lines;do
    SAMPLE1=`echo $lines|awk '{print $1}'`
    SAMPLE1_DEV=`echo $lines|awk '{print $2}'`
    SAMPLE2=`echo $lines|awk '{print $3}'`
    SAMPLE2_DEV=`echo $lines|awk '{print $4}'`
    pair_DEV=`echo $lines|awk '{print $5}'`

    for chr in `echo $chr_lis`;do
        for item in {deletion,duplication};do
        b=`du -h ${SAMPLE1_DEV}/${chr}.mask_CNV_${item} | awk '{print $1}'`
        if [ "${b}" == 0 ];then
            #echo 'file zero' 
            awk '{print $1"\t"$2"\t"$3"\t"SAMPLE2"'"$item"'""_CNV"}' SAMPLE2="$SAMPLE2" ${SAMPLE2_DEV}/${chr}.mask_CNV_${item}  > ${pair_DEV}/${chr}_${SAMPLE2}to${SAMPLE1}_own_CNV_${item}
        else
            #echo 'file not zero' 
            awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==0) print $1"\t"$2"\t"$3"\t"SAMPLE2"'"$item"'""_CNV"}' SAMPLE2="$SAMPLE2" ${SAMPLE1_DEV}/${chr}.mask_CNV_${item} ${SAMPLE2_DEV}/${chr}.mask_CNV_${item} > ${pair_DEV}/${chr}_${SAMPLE2}to${SAMPLE1}_own_CNV_${item}
        fi
        b=`du -h ${SAMPLE2_DEV}/${chr}.mask_CNV_${item} | awk '{print $1}'`
        #echo ${b}
        if [ "${b}" == 0 ];then
            #echo 'file zero'
            awk '{print $1"\t"$2"\t"$3"\t"SAMPLE1"'"$item"'""_CNV"}' SAMPLE1="$SAMPLE1" ${SAMPLE1_DEV}/${chr}.mask_CNV_${item}  > ${pair_DEV}/${chr}_${SAMPLE1}to${SAMPLE2}_own_CNV_${item}
        else
            #echo 'file not zero'
            awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==0) print $1"\t"$2"\t"$3"\t"SAMPLE1"'"$item"'""_CNV"}' SAMPLE1="$SAMPLE1" ${SAMPLE2_DEV}/${chr}.mask_CNV_${item} ${SAMPLE1_DEV}/${chr}.mask_CNV_${item} > ${pair_DEV}/${chr}_${SAMPLE1}to${SAMPLE2}_own_CNV_${item}
        fi
        awk 'NR==FNR{a[$2]=$0;next}NR>FNR{if($2 in a ==1) print $1"\t"$2"\t"$3"\t""'"$item"'""_both_CNV"}' ${SAMPLE1_DEV}}/${chr}.mask_CNV_${item} ${SAMPLE2_DEV}}/${chr}.mask_CNV_${item} > ${pair_DEV}/${chr}_${SAMPLE1}to${SAMPLE2}both_CNV_${item}
    #3  
        cat ${pair_DEV}/${chr}_${SAMPLE1}to${SAMPLE2}_own_CNV_${item} ${pair_DEV}/${chr}_${SAMPLE2}to${SAMPLE1}_own_CNV_${item} ${pair_DEV}/${chr}_${SAMPLE1}to${SAMPLE2}both_CNV_${item} > ${pair_DEV}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV_${item}
        done
        cat ${pair_DEV}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV_deletion ${pair_DEV}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV_duplication | sort -u -nk2,2 > ${SAMPLE1}_${SAMPLE2}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV
    done
done < $config
