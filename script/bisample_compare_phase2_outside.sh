#!/bin/bash
while getopts 'c:m:t:' opt
do
    case "$opt" in
    c) SAMPLE_list=$OPTARG ;;
    esac
done
if [[ -z "${SAMPLE_list}" ]]; then
    echo "No first sample provided"; exit 1
fi

#SAMPLE_list=${SAMPLE1}
#arr=(`echo ${SAMPLE_list} | tr ',' ' '`)
arritemidx(){
  local tmp
  local count=0
  local array=`echo $1`
  for tmp in ${array[@]};do
    if test $2 = $tmp;then
      echo $count
      return
    fi
    count=$(( $count + 1 ))
  done
  echo -1
}

arrslice(){
  array=($1)
  if [ $2 == -1 ];then
    echo ${array[@]}
  elif [ $2 == 0 ];then
    echo ${array[@]:1}
  else
    #echo "${array[@]:0:$2} ${array[@]:$(( $2 + 1 ))}"
    echo "${array[*]:0:$2} ${array[*]:$(( $2 + 1 ))}"
  fi
}
#for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
for CHR in chr5B;do
parallel -j procfile sh bisample_compare_phase2_inside.sh ::: $(eval cat ${SAMPLE_list}) ::: ${CHR}
done

          #cat ${SAMPLE1}_${SAMPLE2}/chr*.1M.delcnv_density > statistic/combine.${SAMPLE1}_${SAMPLE2}_1M.delcnv_density
#
#
