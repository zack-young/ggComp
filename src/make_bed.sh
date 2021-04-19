#!/bin/bash

bedtools makewindows -g /data/genome/wheat/CS_IWGSC/v1/161010_Chinese_Spring_v1.0_pseudomolecules.fasta.fai -w 1000000 > CS_win1M.bed
gawk '{print > $1".1M.bed"}' CS_win1M.bed

for CHR in `cut -f1 CS_win1M.bed|uniq`;do
  declare -i i=0
  total=$(($(wc -l < ${CHR}.1M.bed)/400))
  for ((i=0;i<=$total;i++));do
    # newname=${CHR}.$(gawk -vi=$i 'BEGIN{print tolower(sprintf("%c",i+65))}').bed
    newname=${CHR}.$(printf "%01d" $((${i}+1))).bed
    tail -n +$((i*400+1)) ${CHR}.1M.bed| head -n 400 > ${newname}
  done
done

