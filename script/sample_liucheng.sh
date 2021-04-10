#!/bin/bash
a=$1
time_date="`date +%Y-%m-%d`"
mkdir /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/${a}_${time_date}
echo ${a} > /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/samples_${time_date}.txt
while read i;do
  while read j;do
    sh /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/tmp/script.sh $i $j ${time_date} &
    wait_all
  done < /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/chrlst
done < /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/samples_${time_date}.txt
wait
# norm       
while read i;do
  sh /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/normDPbySample.sh $i "1k.DP" "norm" 4 ${time_date} &
  wait_all
done < /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/samples_${time_date}.txt
wait
while read line;do
  sh /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/script_1M.sh $line ${time_date} &
done < /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/samples_${time_date}.txt

