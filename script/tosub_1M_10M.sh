#!/bin/bash
set -euxo pipefail
while read line;do
  #if [[ ! -d "${line}" ]]; then
  #  mkdir ${line}
  #fi
  #path_SM="/data2/rawdata2/readDepth_ori_1k/${line}"
#  for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
  #    gawk '{sum+=$1} (NR%10)==0{print sum/10; sum=0} END{n=(NR%10);if(n!=0){print sum/n}}' ${line}/${chr}.1.1M.tmp > ${line}/${chr}.1.10M.tmp
  #    gawk '{sum+=$1} (NR%10)==0{print sum/10; sum=0} END{n=(NR%10);if(n!=0){print sum/n}}' ${line}/${chr}.2.1M.tmp > ${line}/${chr}.2.10M.tmp
  #done
#  sed '1d' ${line}/${chr}.1M.norm | gawk '{sum+=$1} (NR%10)==0{print sum/10; sum=0} END{n=(NR%10);if(n!=0){print sum/n}}' > ${line}/${chr}.10M.n.norm
#  sed -i '1i '${line}'' ${line}/${chr}.10M.n.norm
#  done
  cat `awk '{print "'"${line}"'""/"$1".10M.n.norm"}' chr_list|tr '\n'     ' '|sed s/,$//g` | grep -v ${line} > ${line}/combine_10M.norm
  sed -i '1i '$line'' ${line}/combine_10M.norm
  #sh script_10M.sh $line &
  #wait_all
  #cat ${line}/chr*1M.norm.chr | grep -v 'chr' > ${line}/combine_1M.norm
  #paste ${line}/chr*.1M.norm.chr > ${line}/${line}.1M.norm.chr
done < cultivar_list

#paste `awk '{print "/data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/"$1"/combine_1M.norm"}' samples.txt |tr '\n' ' '|sed s/,$//g` > cultivar_combine_norm
