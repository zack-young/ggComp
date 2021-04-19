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
chr_lis='chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D'
pair_sample_lis="off"
CNV_detection="off"
SNP_extract="off"
DSR_count="off"
SGR_PHR_define="off"
while getopts 'hc:v:a:s:m:b:d:f:n:' opt
do
    case "$opt" in
    h|help ) usage; exit 0 ;;
    p|generate_pair_lis)

    d) meta_data=$OPTARG ;;
    f) lis1=$OPTARG ;;
    n) lis2=$OPTARG ;;
    n) pairlis=$OPTARG ;;
    n) chrlis=$OPTARG ;;
    BCF_dir




    \? ) echo -e "\n Option does not exist : $OPTARG\n"
    usage; exit 1 ;;
    esac
done


if test $generate_pair_lis = "on" ;then  # generate sample list
    if [[ -z "${sample_list1}" ]]; then
        echo "No sample_list1 file provided."; exit 1
    fi
    if [[ -z "${sample_list2}" ]]; then
        echo "No sample_list2 file provided."; exit 1
    fi
    if [[ -z "${sample_list_pairwise}" ]]; then
        echo "No sample_list_pairwise file provided."; exit 1
    fi
    sh bisample_compare_list_generate.sh ${sample_list1} ${sample_list2} > ${sample_list_pairwise}
fi


if test $CNV_detection = "on" ;then    # call CNV
    if [[ -z "${meta_data}" ]]; then
        echo "No meta_data file provided."; exit 1
    fi
    if [[ -z "${sample_list_pairwise}" ]]; then
        echo "No sample_list_pairwise file provided."; exit 1
    fi
    sh bisample_compare_phase1.sh -l ${meta_data} -s $sample_list_pairwise
fi

if test $SNP_extract = "on" ;then    # call CNV
    if [[ -z "${meta_data}" ]]; then
        echo "No meta_data file provided."; exit 1
    fi
    if [[ -z "${BCF_dir}" ]]; then
        echo "No meta_data file provided."; exit 1
    fi
    if [[ -z "$sample_list_pairwise" ]]; then
        echo "No pairwise sample list file provided."; exit 1
    fi
    if [[ -z "$chr_lis" ]]; then
        echo "No chromosome list provided."; exit 1
    fi

    declare -A arr_3
    while read line;do
    arr_3[`echo $line | awk '{print $1}'`]=`echo $line | awk '{print $2}'`
    done < ${meta_data}

    while read lines;do
        SAMPLE1=`echo $lines|awk -F ',' '{print $1}'`
        SAMPLE2=`echo $lines|awk -F ',' '{print $2}'`
        if [[ ! -d "${SAMPLE1}_${SAMPLE2}" ]]; then
            mkdir ${SAMPLE1}_${SAMPLE2}
        fi
        for chr in $chr_lis;do
            VCF=${BCF_dir}"/"${chr}".bcf.gz"
            echo "${VCF}"    
            bcftools view -v snps -e 'MAF<0.01' --min-ac=1 -M2 -m2 ${VCF} -Ou| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%DP\t%GQ\t]\n'  -s ${arr_3[${SAMPLE1}]},${arr_3[${SAMPLE2}]} > ${SAMPLE1}_${SAMPLE2}/${CHR}_snp.gt
        done
    done < $sample_list_pairwise
fi

if test $DSR_count = "on" ;then
    if [[ -z "$sample_list_pairwise" ]]; then
        echo "No pairwise sample list file provided."; exit 1
    fi
## level "1.5,2.85"
    while read lines;do
        SAMPLE1=`echo $lines|awk -F ',' '{print $1}'`
        SAMPLE2=`echo $lines|awk -F ',' '{print $2}'`
        ./muti_func_snp_compare.py --DSR_count on -i ${SAMPLE1}_${SAMPLE2}/${chr}_snp.gt  --DP_low ${DP_low} --DP_high ${DP_high} --GQ_sample ${GQ_sample} -b ${bin_size} -o ${SAMPLE1}_${SAMPLE2}/${chr}.diff_density
    done < $sample_list_pairwise
fi

if test $SGR_PHR_define = "on" ;then
    if [[ -z "$sample_list_pairwise" ]]; then
        echo "No pairwise sample list file provided."; exit 1
    fi
    while read lines;do
        SAMPLE1=`echo $lines|awk -F ',' '{print $1}'`
        SAMPLE2=`echo $lines|awk -F ',' '{print $2}'`
        for chr in  $chr_lis;do
            ./muti_func_snp_compare.py --two_diff_level on -i ${SAMPLE1}_${SAMPLE2}/${chr}.diff_density  --LEVEL ${LEVEL} -o ${SAMPLE1}_${SAMPLE2}/${chr}.homo_undefined_snp_level  &
        done
    done < $sample_list_pairwise
    if test $combine_CNV = "on" ;then
        while read lines;do
            SAMPLE1=`echo $lines|awk -F ',' '{print $1}'`
            SAMPLE2=`echo $lines|awk -F ',' '{print $2}'`
            line=`cat ${SAMPLE1}_${SAMPLE2}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV|wc -l`
            if test $line = 0;then
                :
            else
                type=`cut -f4 ${SAMPLE1}_${SAMPLE2}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV`
                awk -F '\t' 'NR==FNR{a[$1]=$0;next}NR>FNR{if ($1 in a==1) {print a[$1]} else {print $0}}' ${SAMPLE1}_${SAMPLE2}/${chr}.${SAMPLE1}to${SAMPLE2}all_CNV ${SAMPLE1}_${SAMPLE2}/${chr}.homo_undefined_snp_level > ${SAMPLE1}_${SAMPLE2}/${chr}.homo_undefined_snp_level
            fi
        done < $sample_list_pairwise
    fi

fi


#./ptf_maker.py -p "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/${SAMPLE1}_${SAMPLE2}_${time_date}" -s ".${SAMPLE1}to${SAMPLE2}all_CNV_deletion" > plotfile/plotfile_${SAMPLE1}_${SAMPLE2}_deletion
#Rscript ~/R/graph_diff_sim.R --sample1 ${SAMPLE1} --sample1_name ${SAMPLE1_name} --sample2 ${SAMPLE2} --sample2_name ${SAMPLE2_name} -s ${time_date}
#Rscript ~/R/graph_CNV.R --sample1 ${SAMPLE1} --sample1_name ${SAMPLE1_name} --sample2 ${SAMPLE2} --sample2_name ${SAMPLE2_name} -s ${time_date}