#!/bin/bash
set -x
#for i in {1..7};do
#    for a in {A,B,D};do
#        paste `sed '1d' metadata_cultivar_final.txt | awk '{print $1"/chr""'"$i"'""'"$a"'"".10M.norm"}' |tr '\n' ' '|sed s/,$//g`|sed '1d' | awk '{print (NR-1)*10000000"\t"$0}' > chr$i$a.10M.join_norm
#    done
#done
#
#sed '1d' metadata_cultivar_final.txt |awk '{print $1}' |tr '\n' '\t'|sed s/,$//g |awk '{print "loc""\t"$0}' > fielder_header
#
#for i in {1..7};do     
#    for a in {A,B,D};do
#        cat fielder_header chr$i$a.1M.join_norm | sponge chr$i$a.1M.join_norm
#    done
#done
#
#for i in {1..7};do
#    for a in {A,B,D};do
#        Rscript /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/CNV_heatmap.R -i chr$i$a.1M.join_norm -M metadata_cultivar_final.txt
#    done
#done

#while read line;do 
#    #a=`echo $line | awk '{print $1}'`; cd $a ; paste chr*.chr > ${a}.chr;cd ..
#    a=`echo $line | awk '{print $1}'`
#    cat $a/chr*.mask_CNV_deletion_split > $a/combine_mask_CNV_deletion_split &
#    cat $a/chr*.mask_CNV_duplication_split > $a/combine_mask_CNV_duplication_split &
#    cat $a/chr*.mask_CNV_deletion_merged > $a/combine_mask_CNV_deletion_merged &
#    cat $a/chr*.mask_CNV_duplication_merged > $a/combine_mask_CNV_duplication_merged &
#    wait_all
#done < cultivar_list
#cat `sed '1d' metadata_cultivar.txt | awk '{print $1"/combine_mask_CNV_deletion_split"}' |tr '\n' ' '|sed s/,$//g`> combine_mask_CNV_deletion_split
#cat `sed '1d' metadata_cultivar.txt | awk '{print $1"/combine_mask_CNV_duplication_split"}' |tr '\n' ' '|sed s/,$//g`> combine_mask_CNV_duplication_split
#cat `sed '1d' metadata_cultivar.txt | awk '{print $1"/combine_mask_CNV_deletion_merged"}' |tr '\n' ' '|sed s/,$//g`> combine_mask_CNV_deletion_merged
#cat `sed '1d' metadata_cultivar.txt | awk '{print $1"/combine_mask_CNV_duplication_merged"}' |tr '\n' ' '|sed s/,$//g`> combine_mask_CNV_duplication_merged

#

#
#while read line;do
#    a=`echo $line | awk '{print $1}'`
#    Rscript /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/bin/CNV_heat_map_bySample.R -i ${a}/${a}.chr -k /data/user/yangzz/mapping/09.filterVCF/GT_AD/CNV_heatmap/CS_karyotype.txt -p pdf/${a}_DPbySample.pdf
#done < rest_cultivar

#while read line;do
#    a=`echo $line | awk '{print $1}'`
#    for chr in chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D;do
#        mkdir ${a} ; sed '1 s/^.*$/'$chr'/' /data2/rawdata2/readDepth/${a}/${chr}.1M.norm > ${a}/${chr}.1M.norm.chr
#    done
#done < metadata_cultivar.txt

#for num in {1..7}; do
#    for i in {A,B,D}; do
#      #1
#        while read line;do
#            a=`echo $line | awk '{print $1}'`
#            /data/user/yangzz/mapping/09.filterVCF/GT_AD/muti_func_snp_compare.py -b 1000000 -i ${a}/chr${num}${i}.1M.norm --mask_cnv on -o ${a}/chr${num}${i}.mask_CNV &
#        done < rest_cultivar
#        wait_all
#    done
#done

#wait
#for num in {1..7}; do
#    for i in {A,B,D}; do
#        while read line;do
#            a=`echo $line | awk '{print $1}'`
#            for item in {deletion,duplication};do       
#                grep "${item}" ${a}/chr${num}${i}.mask_CNV > ${a}/chr${num}${i}.mask_CNV_${item} 
#            done
#        done < rest_cultivar             
#    done
#done
#while read line;do
#    for num in {1..7}; do
#        for i in {A,B,D}; do
#            for item in {deletion,duplication};do
#                a=`echo $line | awk '{print $1}'`
#                #bedtools merge -i ${a}/chr${num}${i}.mask_CNV_${item} -d 1000000 -c 4 -o distinct  > ${a}/chr${num}${i}.mask_CNV_${item}_merged &
#                ./count_CNV.py -i ${a}/chr${num}${i}.mask_CNV_${item}_merged --split on >  ${a}/chr${num}${i}.mask_CNV_${item}_split &
#            done
#        wait_all
#        done
#    done
#done < cultivar_list

for num in {1..7}; do
    for i in {A,B,D}; do
        for item in {deletion,duplication};do
            #cat `awk -vnum=$num -vi=$i -vitem=$item '{print $1"/chr"num""i".mask_CNV_"item"_split"}' cultivar_list |tr '\n' ' '|sed s/,$//g`|sort -nk1,1 - > field_CNV_count/chr${num}${i}.combine_mask_CNV_${item}
            ./count_CNV.py --count on --pop_num 119 -i field_CNV_count/chr${num}${i}.combine_mask_CNV_${item} |sort -nk1,1 >  field_CNV_count/chr${num}${i}.combine_mask_CNV_${item}_count 
            ./count_CNV.py --compensent_${item} on --chrom ${num}${i} -i field_CNV_count/chr${num}${i}.combine_mask_CNV_${item}_count -o field_CNV_count/chr${num}${i}.combine_mask_CNV_${item}_compensent 

        done
    done
done




