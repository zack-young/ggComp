#!/usr/bin/env usr/sh
set -euxo pipefail

SM=$1

ave=`./count_CNV.py --normalize on -i ${1}/combine_1M_DP`

for CHR in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
  gawk -vOFS="\t" -vave=$ave '{print $4/ave}' ${1}/${CHR}.1.1M.tmp > ${1}/${CHR}.1M.norm_un1k
  gawk -vOFS="\t" -vave=$ave '{print $4/ave}' ${1}/${CHR}.2.1M.tmp >> ${1}/${CHR}.1M.norm_un1k
  sed -i '1i '$SM'' ${1}/${CHR}.1M.norm_un1k
  sed '1 s/^.*$/'$CHR'/' ${1}/${CHR}.1M.norm_un1k > ${1}/${CHR}.1M.norm_un1k.chr
done

