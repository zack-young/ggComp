#!/bin/bash

while read line;do
  if [[ ! -d "${line}" ]]; then
    mkdir ${line}
  fi
  path_SM="/data2/rawdata2/readDepth_ori_1k/${line}"
  for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
      gawk '{sum+=$4} (NR%1000)==0{print sum/1000; sum=0} END{n=(NR%1000);if(n!=0){print sum/n}}' $path_SM/${chr}.1.1k.DP > ${line}/${chr}.1.1M.tmp
      gawk '{sum+=$4} (NR%1000)==0{print sum/1000; sum=0} END{n=(NR%1000);if(n!=0){print sum/n}}' $path_SM/${chr}.2.1k.DP > ${line}/${chr}.2.1M.tmp
  done
  cat ${line}/chr*.1M.tmp > ${line}/combine_1M_DP
  sh script_1M.sh $line &
  wait_all
  #cat ${line}/chr*1M.norm.chr | grep -v 'chr' > ${line}/combine_1M.norm
  #paste ${line}/chr*.1M.norm.chr > ${line}/${line}.1M.norm.chr
done < rest_cultivar

#paste `awk '{print "/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/"$1"/combine_1M.norm"}' samples.txt |tr '\n' ' '|sed s/,$//g` > cultivar_combine_norm
