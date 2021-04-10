#!/bin/bash
set -x
#===============================================================================
# FILE: shell_pack.sh
# 
# USAGE: ./shell_pack.sh 
# 
# DESCRIPTION: 
# 
# OPTIONS: ---
# REQUIREMENTS: ---
# BUGS: ---
# NOTES: ---
# AUTHOR: lwq (28120), scue@vip.qq.com
# ORGANIZATION: 
# CREATED: 04/22/2015 02:38:01 PM CST
# REVISION: ---
#===============================================================================

#=== FUNCTION ================================================================
# NAME: usage
# DESCRIPTION: Display usage information.
#===============================================================================
function usage ()
{
 cat <<- EOT

Usage : $0 -p package -s script file1 file2 file3 ..

Options:
 -h|help Display this message
 -p|package The output package name
 -s|script The script will run when unpack package
 Other The all files what you want to pack

EOT
} # ---------- end of function usage ----------

#-----------------------------------------------------------------------
# Handle command line arguments
#-----------------------------------------------------------------------

DP_low=3
DP_high=99
GQ_sample=8
bin_size=1000000
LEVEL="1.5,3.0"

while getopts 'hc:v:a:s:m:b:d:f:n:' opt
do
    case "$opt" in
    h|help ) usage; exit 0 ;;
    d) total=$OPTARG ;;
    f) first=$OPTARG ;;
    n) second=$OPTARG ;;
    \? ) echo -e "\n Option does not exist : $OPTARG\n"
    usage; exit 1 ;;
    esac
done


if [[ -z "${total}" ]]; then
    echo "No total sample bcf name provided."; exit 1
fi


sh bisample_compare_list_generate.sh ${meta_data} ${meta_data} > sample_list.txt

# call CNV
sh bisample_compare_phase1.sh -l ${meta_data} -s sample_list.txt

declare -A arr_3
while read line;do
arr_3[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $2}'`
done < ${meta_data}

while read lines;do
    SAMPLE1=`echo $lines|awk -F ',' '{print $1}'`
    SAMPLE2=`echo $lines|awk -F ',' '{print $2}'`
    for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
        VCF=${BCF_dir}"/"${chr}".bcf.gz"
        echo "${VCF}"    
        bcftools view -v snps -e 'MAF<0.01' --min-ac=1 -M2 -m2 ${snp_path}/${VCF} -Ou| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%DP\t%GQ\t]\n'  -s ${arr_3[${SAMPLE1}]},${arr_3[${SAMPLE2}]} > ${SAMPLE1}_${SAMPLE2}/${CHR}_snp.gt
    done
done < sample_list.txt
wait


## level "1.5,2.85"
while read lines;do
    SAMPLE1=`echo $lines|awk -F ',' '{print $1}'`
    SAMPLE2=`echo $lines|awk -F ',' '{print $2}'`
    for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
        ./muti_func_snp_compare.py --DSR_count on -i ${SAMPLE1}_${SAMPLE2}/${chr}_snp.gt  --DP_low ${DP_low} --DP_high ${DP_high} --GQ_sample ${GQ_sample} -b ${bin_size} -o ${SAMPLE1}_${SAMPLE2}/${chr}.diff_density
    done
done < sample_list.txt
wait

while read lines;do
    SAMPLE1=`echo $lines|awk -F ',' '{print $1}'`
    SAMPLE2=`echo $lines|awk -F ',' '{print $2}'`
    for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
        ./muti_func_snp_compare.py --two_diff_level on -i ${SAMPLE1}_${SAMPLE2}/${chr}.diff_density  --LEVEL ${LEVEL} -o ${SAMPLE1}_${SAMPLE2}/${chr}.homo_undefined_snp_level  &
    done
done < sample_list.txt
wait

while read lines;do
    SAMPLE1=`echo $lines|awk -F ',' '{print $1}'`
    SAMPLE2=`echo $lines|awk -F ',' '{print $2}'`
    line=`cat dev/${SAMPLE1}_${SAMPLE2}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV|wc -l`
    if test $line = 0;then
        cat ${SAMPLE1}_${SAMPLE2}/${chr}.only_homo_snp_undefined_level  > ${SAMPLE1}_${SAMPLE2}/${chr}.homo_undefined_snp_level
    else
        type=`cut -f4 dev/${SAMPLE1}_${SAMPLE2}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV`
        awk -F '\t' 'NR==FNR{a[$1]=$0;next}NR>FNR{if ($1 in a==1) {print a[$1]} else {print $0}}' ${SAMPLE1}_${SAMPLE2}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV ${SAMPLE1}_${SAMPLE2}/${chr}.only_homo_snp_undefined_level > ${SAMPLE1}_${SAMPLE2}/${chr}.homo_undefined_snp_level
    fi
done < sample_list.txt
#./ptf_maker.py -p "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/${SAMPLE1}_${SAMPLE2}_${time_date}" -s ".${SAMPLE1}to${SAMPLE2}all_CNV_deletion" > plotfile/plotfile_${SAMPLE1}_${SAMPLE2}_deletion
#Rscript ~/R/graph_diff_sim.R --sample1 ${SAMPLE1} --sample1_name ${SAMPLE1_name} --sample2 ${SAMPLE2} --sample2_name ${SAMPLE2_name} -s ${time_date}
#Rscript ~/R/graph_CNV.R --sample1 ${SAMPLE1} --sample1_name ${SAMPLE1_name} --sample2 ${SAMPLE2} --sample2_name ${SAMPLE2_name} -s ${time_date}