#!/data/user/yangzz/worktools/anaconda3/bin/python
###  #!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""
@author:youngzhengzhao
@file:finalgeno.py
@time:2018/7/2417:35
"""

from __future__ import division  # very important !!!!
from collections import defaultdict
import os
import gzip
import sys
import re
import math
import pandas as pd
import numpy as np
def error(msg):
    print(sys.stderr, 'ERROR: %s' % msg)
    exit(1)
#
def IsGeno(geno):
    if len(geno)== 3 and geno[1] =='/' and geno!= "./.":
        return 1
    else:
        return 0
    #
#
def IsALT(geno):
    if geno[0] != '0' or geno[2] != "0":
        return 1
    else:
        return 0
    #
#    
def IsHete(geno):
    if geno[0] == geno[2]:
        return 0
    else:
        return 1
    #
#
# TODO: 增加更改参数功能
def  IsHomo(geno):  #判断纯合
    if geno[0] == geno[2]:
        return 1
    else:
        return 0
    #
#
def  Isunconsistent(g1, g2):   #判断父母本差异密度
    if g1 == g2: #and g2 == gs :
        return 0
    # if IsGeno(g1) and IsGeno(g2) and IsGeno(gs) :
    #     return 1  
    else :
        return 1
    #
#

def OutputLinehomohigh(chr_cur, start, end, count, count1, count2) :
    if count1 == count2 :
        state = "NA"
    elif count1 > count2 :
        state = "lx987"
    else :
        state = "3097"
    #
    print("\t".join([chr_cur, "%s" % start, "%s" % end,
                     "%s" % count, "%s" % count1, "%s" % count2, state]) )
    #
#

def OutputLinehomomid(chr_cur, start, end, count, count1, count2) :
    if count1 == count2:
        state = "NA"
    elif count1 > count2:
        state = "lx987"
    else:
        state = "3097"
    #
    print("\t".join([chr_cur, "%s" % start, "%s" % end,
                     "%s" % count, "%s" % count1, "%s" % count2, state]) )
    #
#
def DefineAD(gt, ad):
    global lis1
    lis1 = []
    item_ad = ad.strip().split(",")
    if IsHete(gt):
        if '1' in item_ad:
            index = [i for i,x in enumerate(item_ad) if x != '1' and x != '0'] # find all not '1' and '0' item index
            if len(index) == 1:
                gt = str(index[0])+'/'+str(index[0])
            elif len(index) == 0:
                gt = '0/0'
            else:  
                gt = str(index[0])+'/'+str(index[1])
    if ad != '.':
        if gt != './.':
            lis1.append(int(item_ad[int(gt[0])])+int(item_ad[int(gt[2])]))
        else:
            lis1.append(sum(map(int, item_ad)))
    else:
        lis1.append(0)
    lis1.append(gt)
    return(lis1)

def ChangeGQ(gt, gq):
    global lis2
    lis2 = []
    if IsGeno(gt):
        lis2.append(int(gq))
    else:
        lis2.append(0)
    return(lis2)

def ChangeGQ_dot(gq):
    global lis3
    lis3 = []
    if gq != '.':
        lis3.append(int(gq))
    else:
        lis3.append(0)
    return(lis3)

def Makelis(chromsome):
    global list1
    list1 = []
    #f = open("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/MASK/chr" + chromsome + "_combine_total", "r")
    f = open(chromsome, "r")
    mask = f.readlines()
    for line in mask:
        item = line.strip("\n").split("\t")
        list1.append([item[1],item[2], item[3]])
    f.close
    return(list1)

def Makedic(sample):
    dic1 = defaultdict(list)
    for j in range(0, 7):
        for i in range(0, 3):
            lett = ["A", "B", "D"]
            dic1[sample].append([])
            f = open(sample + "/chr" + str(j+1) + lett[i] + ".1M.norm", "r")
            CNVT = f.readlines()
            f.close()
            #lis = []
            for line in CNVT:
            #     lis.append(line)
            # for line in lis:
                item = line.strip()
                dic1[sample][i+j*3].append(item)
    return(dic1)

def GenoPhase(infile, bin_size, gene_snp_filter, chrom_graphy, density_graph_raw_vcf_two,
              two_diff_level, sample1, sample2, time_data, total,
              density, snp_graphy, mask_DP, long_short_count, density_graph_raw_vcf_one,
              two_diff, grapht, per_snp_filter, DP_GQ, filter_and_count_snp, mask_cnv,
              high_threshold, DP_high, DP_low, GQ_sample, LEVEL, centrosome, test, filter_and_count_indel,
              sample, chromosome, cnv_region, transition_transversion_count, indel_count, parents_define,
              chr_col=1, pos_col=2, ref_col=3, alt_col=4, s1_col=5, s2_col=4, s3_col=5):
    # open the infile
    try:
        if infile:
            if infile.endswith(".gz"):
                IN = gzip.open(infile, 'rb')
            else:
                IN = open(infile, 'r')
            #
        else:
            IN = sys.stdin
        #
    except IOError:
        print("\n[Error]:\n\t File cannot be open: %s" % infile)
        exit(-1)
    #
    # init
    count1 = 0
    count2 = 0
    count1_1 = 0
    count2_1 = 0
    count3 = 0
    count4 = 0
    count5 = 0
    count6 = 0
    count7 = 0
    count8 = 0
    count9 = 0
    count5_1 = 0
    count6_1 = 0
    count7_1 = 0
    count8_1 = 0
    DP_count = 0
    length = {'chr1A': '594102056', 'chr1B': '689851870', 'chr1D': '495453186',
              'chr2A': '780798557', 'chr2B': '801256715', 'chr2D': '651852609',
              'chr3A': '750843639', 'chr3B': '830829764', 'chr3D': '615552423', 
              'chr4A': '744588157', 'chr4B': '673617499', 'chr4D': '509857067',
              'chr5A': '709773743', 'chr5B': '713149757', 'chr5D': '566080677',
              'chr6A': '618079260', 'chr6B': '720988478', 'chr6D': '473592718',
              'chr7A': '736706236', 'chr7B': '750620385', 'chr7D': '638686055',
              'chrUn': '480980714'}
#    Makedic(sample)
    #

    if gene_snp_filter == 'on':   
        '''
        find unmatch bin 
        构建好DP区间，ownCNV区间、unmatchCNV区间后输出含
        有不同基因型的目标区域
        '''                             
        firstline = IN.readline().strip()
        item = firstline.strip().split("\t")
        firstpos = int(item[pos_col - 1])
        if int(firstpos) <= int(bin_size):
            start = 1
            end = bin_size
        else:
            num = int(firstpos/bin_size)-1
            start = (num+1) * bin_size
            end = start + bin_size
        IN.seek(0, os.SEEK_SET)
        PART = infile.split("_")
        part = PART[0]
        art = part.lstrip('s')
        f = open("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/MASK/lx99_jm22_mask/chr" + chromsome + ".mask_filtered_combine", "r")
        mask = f.readlines()
        f.close
        a = 0
        for line in IN:
            tokens = line.strip().split("\t")
            chr = tokens[chr_col - 1]
            pos = int(tokens[pos_col - 1])
            g1 = tokens[s1_col - 1]
            g2 = tokens[s2_col - 1]
            DP1 = int(tokens[s1_col])
            DP2 = int(tokens[s2_col])
            chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
            b = chr_cur[0]
            START = str(start)
            judge_num = 0
            while pos > end:
                # #
                start += bin_size
                end += bin_size
            for ther in mask:
                item = ther.strip("").split("\t")               
                if item[2] <= START <= item[3]:
                    judge_num += 1      
            if judge_num == 0:
                if IsGeno(g1) and IsGeno(g2) and IsHomo(g1) and IsHomo(g2):
                    if int(DP_low) <= DP1 <= int(DP_high) and int(DP_low) <= DP2 <= int(DP_high):
                        #if IsALT(g1) or IsALT(g2):
                            #count4 += 1
                        if Isunconsistent(g1, g2):
                            count3 += 1
                            #print(line)
    #
    if chrom_graphy == 'on':   # find discard bin
        firstline = IN.readline().strip()
        item = firstline.strip().split("\t")
        firstpos = int(item[pos_col - 1])
        if int(firstpos) <= int(bin_size):
            start = 1
            end = bin_size
        else:
            num = int(firstpos/bin_size)-1
            start = (num+1) * bin_size
            end = start + bin_size
        IN.seek(0, os.SEEK_SET)
        PART = infile.split("_")
        part = PART[0]
        art = part.lstrip('s')
        for line in IN:
            tokens = line.strip().split("\t")
            chr = tokens[chr_col - 1]
            pos = int(tokens[pos_col - 1])
            ref = tokens[ref_col - 1]
            alt = tokens[alt_col - 1]
            GT = tokens[gt_col - 1]
            DP = int(tokens[5])
            chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
            b = chr_cur[0]
            dic2 = {"A":"0","B":"1","D":"2"}
            num = int(b[3])-1
            ch = int(dic2[b[4]])
            count = int(pos/bin_size)
            while pos > end:
                ratio = (count2+1)/(count1+1)
                DP_ratio = DP_count/count3
                if float(dic1[sample][(num)*3+ch][count]) >= 0.6:
                    if lowun <= DP_ratio <= highun:
                        if ratio >= maxun:
                            print("\t".join([b, art, "%s" % start, "%s" % end, "hete"]))
                        else:
                            print("\t".join([b, art, "%s" % start, "%s" % end, "homo"]))
                    else:
                        print("\t".join([b, art, "%s" % start, "%s" % end, "DP"]))
                else:
                    print("\t".join([b, art, "%s" % start, "%s" % end, "CNV"]))
                count1 = 0
                count2 = 0
                count3 = 0
                DP_count = 0
                start += bin_size
                end += bin_size
                #
                # compare genotype
                DP_count += DP
                count3 += 1
                if IsHomo(GT):                 #continue 会跳出整个循环！！！！
                    count1 += 1                #pass 不起作用 仅是起到站位作用
                if IsHete(GT):
                    count2 += 1
        #
            #
        count += 1
        ratio = (count2+1)/(count1+1)
        DP_ratio = DP_count/count3
        if float(dic1[sample][(num)*3+ch][count]) >= 0.6:
            if lowun <= DP_ratio <= highun:
                if ratio >= maxun:
                    print("\t".join([b, art, "%s" % start, "%s" % end, "hete"]))
                else:
                    print("\t".join([b, art, "%s" % start, "%s" % end, "homo"]))
            else:
                print("\t".join([b, art, "%s" % start, "%s" % end, "DP"]))
        else:
            print("\t".join([b, art, "%s" % start, "%s" % end, "CNV"]))
        #
    #
    if density == 'on':   # counting DP per bin
        firstline = IN.readline().strip()
        item = firstline.strip().split("\t")
        firstpos = int(item[pos_col - 1])
        if int(firstpos) <= int(bin_size):
            start = 1
            end = bin_size
        else:
            num = int(firstpos/bin_size)-1
            start = (num+1) * bin_size
            end = start + bin_size
        IN.seek(0, os.SEEK_SET)
        PART = infile.split("_")
        part = PART[0]
        art = part.lstrip('s')
        Makelis(chromosome)
        #Makedic(sample)
        for line in IN:
            tokens = line.strip().split("\t")
            chr = tokens[chr_col - 1]
            pos = int(tokens[pos_col - 1])
            ref = tokens[ref_col - 1]
            alt = tokens[alt_col - 1]
            GT = tokens[s1_col - 1]
            DP = int(tokens[s1_col])
            GQ = int(tokens[s1_col+1])
            chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
            b = chr_cur[0]
            dic2 = {"A":"0","B":"1","D":"2"}
            num = int(b[3])-1
            ch = int(dic2[b[4]])
            count = int(pos/bin_size)
            while pos > end:                   
                #if count3 == 0:
                #    print("\t".join([b,"%s" % start, "%s" % end,"%s" % count3]))
                #    DP_ratio = 0.01
                #else:
                DP_ratio = DP_count/count3
                START = str(start)
                if START in list1:
                    pass
                else:
                    print("\t".join([b,"%s" % start, "%s" % end, "%.3f" % DP_ratio,"%s" % count3]))
                #if float(dic1[sample][(num)*3+ch][count]) >= 0.6:
                #    print("\t".join([b,"%s" % start, "%s" % end, "%.3f" % DP_ratio,"%s" % count3]))
                count1 = 0
                count2 = 0
                count3 = 0
                DP_count = 0
                start += bin_size
                end += bin_size
                #
            if IsGeno(GT) and IsHomo(GT):
                if int(DP_low) <= DP <= int(DP_high):
                    count3 += 1
                    DP_count += DP
            if IsHomo(GT):                 #continue 会跳出整个循环！！！！
                count1 += 1                #pass 不起作用 仅是起到站位作用
            if IsHete(GT):                 
                count2 += 1                
        DP_ratio = DP_count/count3
        START = str(start)
        if START in list1:
            pass
        else:
            print("\t".join([b,"%s" % start, "%s" % end, "%.3f" % DP_ratio,"%s" % count3]))
        #if float(dic1[sample][(num)*3+ch][count]) >= 0.6:
        #    print("\t".join([b,"%s" % start, "%s" % end, "%.3f" % DP_ratio,"%s" % count3]))
    #
    if mask_DP == 'on':   # mask abnormal DP bin 
        firstline = IN.readline().strip()
        item = firstline.strip().split("\t")
        firstpos = int(item[pos_col - 1])
        if int(firstpos) <= int(bin_size):
            start = 1
            end = bin_size
        else:
            num = int(firstpos/bin_size)-1
            start = (num+1) * bin_size
            end = start + bin_size
        IN.seek(0, os.SEEK_SET)
        PART = infile.split("_")
        part = PART[0]
        art = part.lstrip('s')
        Makelis(chromsome)
        #Makedic(sample)
        for line in IN:
            tokens = line.strip().split("\t")
            chr = tokens[chr_col - 1]
            pos = int(tokens[pos_col - 1])
            ref = tokens[ref_col - 1]
            alt = tokens[alt_col - 1]
            GT = tokens[s1_col - 1]
            #if tokens[s1_col] == '.':
            #    continue
            DP = int(tokens[s1_col])
            # if dp != ".":
            #     dp_lis =dp.strip().split(",")
            # DP = 0
            # for i in dp_lis:
            #     DP += int(i)
            chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
            b = chr_cur[0]
            dic2 = {"A":"0","B":"1","D":"2"}
            num = int(b[3])-1
            ch = int(dic2[b[4]])
            count = int(pos/bin_size)
    
            while pos > end:
                START = str(start)
                if START in list1:
                    pass
                #if float(dic1[sample][(num)*3+ch][count]) < 0.6:
                #    pass
                #    print("\t".join([b, art,"%s" % start, "%s" % end, "CNV"]))
                else:
                    if count3 == 0:
                        print("NA")
                    else:
                        DP_ratio = DP_count/count3
                        if DP_ratio > high_threshold:#lx99
                    #if DP_ratio > 6: #jm22
                            print("\t".join([b, art,"%s" % start, "%s" % end, "DP"]))
                    
                #    
                count1 = 0
                count2 = 0
                count3 = 0
                DP_count = 0
                count3 = 0             
                start += bin_size
                end += bin_size            
                # compare genotype182         
            if IsGeno(GT) and IsHomo(GT):
                if int(DP_low) <= DP <= int(DP_high):
                    count3 += 1
                    DP_count += DP
            if IsHomo(GT):                 #continue 会跳出整个循环！！！！
                count1 += 1                #pass 不起作用 仅是起到站位作用
            if IsHete(GT):                 
                count2 += 1      
            #
        count += 1
        START = str(start)
        if START in list1:
            pass
        #if float(dic1[sample][(num)*3+ch][count]) < 0.6:
        #    print("\t".join([b, art,"%s" % start, "%s" % end, "CNV"]))
        else:
            if count3 == 0:
                print("NA")
            else:
                DP_ratio = DP_count/count3
                if DP_ratio > high_threshold:#lx99
            #if DP_ratio > 6: #jm22
                    print("\t".join([b, art,"%s" % start, "%s" % end, "DP"]))
    #
    if long_short_count == 'on':     # counting long and short arm respectively
        firstline = IN.readline().strip()
        item = firstline.strip().split("\t")
        firstpos = int(item[pos_col - 1])
        count_same=0
        count_diff=0
        count_same_1=0
        count_diff_1=0
        count_same_num=0
        count_diff_num=0        
        if int(firstpos) <= int(bin_size):
            start = 1
            end = bin_size
        else:
            num = int(firstpos/bin_size)-1
            start = (num+1) * bin_size
            end = start + bin_size
        IN.seek(0, os.SEEK_SET)
        PART = infile.split("_")
        part = PART[0]
        art = part.lstrip('s')
        list1 = Makelis(chromosome)
        #print(list1)
        DP_high = DP_high.strip().split(",")
        DP_low = DP_low.strip().split(",")
        for line in IN:
            tokens = line.strip().split("\t")
            chr = tokens[chr_col - 1]
            pos = int(tokens[pos_col - 1])
            ref = tokens[ref_col - 1]
            alt = tokens[alt_col - 1]
            chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
            b = chr_cur[0]
            dic2 = {"A":"0","B":"1","D":"2"}
            num = int(b[3])-1
            ch = int(dic2[b[4]])
            count = int(pos/bin_size)
            #
            len_dic = defaultdict(list)
            #
            GT1 = tokens[s1_col - 1]
            AD1 = tokens[s1_col]
            GQ1 = tokens[s1_col + 1]
            lis3 = ChangeGQ_dot(GQ1)
            GQ1 = lis3[0]
            lis1 = DefineAD(GT1, AD1)
            DP1 = lis1[0]
            GT1 = lis1[1]
            #
            GT2 = tokens[s2_col - 1]
            AD2 = tokens[s2_col]
            GQ2 = tokens[s2_col + 1]
            lis3_1 = ChangeGQ_dot(GQ2)
            GQ2 = lis3_1[0]
            lis1_1 = DefineAD(GT2, AD2)
            DP2 = lis1_1[0]
            GT2 = lis1_1[1]
            #
            ind = (int(b[3])-1)*3 + int(dic2[b[4]])
            #print(ind,int(chromsome[0]),int(dic2[chromsome[1]]))
            centro = centrosome.strip().split(",")
            while pos > end:
                START = str(start)
                #deli = abs(count2-count1)
                #ratio = count3/deli
                if START in list1:
                    pass
                else:
                    if end < int(centro[ind])*1000000:
                        count5 += count1
                        count6 += count2
                        count7 += count3
                        count8 += count4
                        count_diff += count_diff_num
                        count_same += count_same_num

                    else:
                        count5_1 += count1
                        count6_1 += count2
                        count7_1 += count3
                        count8_1 += count4
                        count_diff_1 += count_diff_num
                        count_same_1 += count_same_num
                # #
                start += bin_size
                end += bin_size
                count1 = 0
                count2 = 0
                count3 = 0
                count4 = 0
                count_diff_num = 0
                count_same_num = 0
                #
            if int(DP_low[0]) <= DP1 <= int(DP_high[0]) and GQ1 >= int(GQ_sample):
                if IsGeno(GT1) and IsALT(GT1):
                    if IsHomo(GT1):
                        count1 += 1
                    if IsHete(GT1):
                        count2 += 1

            if int(DP_low[0]) <= DP2 <= int(DP_high[0]) and GQ2 >= int(GQ_sample):
                if IsGeno(GT2) and IsALT(GT2):
                    if IsHomo(GT2):
                        count3 += 1
                    if IsHete(GT2):
                        count4 += 1                
            if int(DP_low[0]) <= DP1 <= int(DP_high[0]) and int(DP_low[0]) <= DP2 <= int(DP_high[0]) and GQ1 >= int(GQ_sample) and GQ2 >= int(GQ_sample):    
                if IsGeno(GT1) and IsGeno(GT2) and IsHomo(GT1) and IsHomo(GT2):
                    if IsALT(GT1) or IsALT(GT2):
                        if Isunconsistent(GT1, GT2):
                            count_diff_num += 1
                        else:
                            count_same_num += 1
        #
        #deli = abs(count2-count1)
        #ratio = count3/deli
        START = str(start)
        if START in list1:
            pass
        else:
            if end < int(centro[ind])*1000000:
                count5 += count1
                count6 += count2
                count7 += count3
                count8 += count4
                count_diff += count_diff_num
                count_same += count_same_num
            else:
                count5_1 += count1
                count6_1 += count2
                count7_1 += count3
                count8_1 += count4
                count_diff_1 += count_diff_num
                count_same_1 += count_same_num
        print("Chr","Length","jm22_homo","jm22_hete","lx99_homo","lx99_hete","same_homo","diff_homo", sep='\t')    
        if 2*int(centro[ind])*1000000 > int(length[b]):
            print(b+"L", int(centro[ind])*1000000, count5, count6, count7, count8, count_same, count_diff, sep='\t')
            print(b+"S", int(length[b])-int(centro[ind])*1000000, count5_1, count6_1, count7_1, count8_1, count_same_1, count_diff_1, sep='\t')
        else:
            print(b+"L", int(length[b])-int(centro[ind])*1000000,count5_1, count6_1, count7_1, count8_1, count_same_1, count_diff_1,sep='\t')
            print(b+"S", int(centro[ind])*1000000,count5, count6, count7, count8, count_same, count_diff, sep='\t')  
    #
    if two_diff == 'on':           # print differnt snp between two samples in snpEFF vcf
        # firstline = IN.readline().strip()
        # item = firstline.strip().split("\t")
        # firstpos = int(item[pos_col - 1])
        # if int(firstpos) <= int(bin_size):
        #     start = 1
        #     end = bin_size
        # else:
        #     num = int(firstpos/bin_size)-1
        #     start = (num+1) * bin_size
        #     end = start + bin_size
        # IN.seek(0, os.SEEK_SET)
        #PART = infile.split("_")
        #part = PART[0]
        #art = part.lstrip('s')
        list1 = []
        list2 = []
        list3 = []
        f = open(chromosome, "r")
        mask = f.readlines()
        for line in mask:
            item = line.strip().split("\t")
            if str(item[3]) == 'low_com':
                list3.append([item[1],item[2]])            
            if str(item[3]) == 'mid_com':
                list1.append([item[1],item[2]])
            if str(item[3]) == 'high_com':
                list2.append([item[1],item[2]])

        f.close
        DP_high = DP_high.strip().split(",")
        DP_low = DP_low.strip().split(",")
        print('CHR','POS','REF','ALT','JiMai22','LiangXing99','ANN',sep='\t')
        for line in IN:
            line_use = line.strip()
            tokens = line.strip().split("\t")
            chr = tokens[chr_col - 1]
            pos = int(tokens[pos_col - 1])
            ref = tokens[ref_col - 1]
            alt = tokens[alt_col - 1]
            chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
            b = chr_cur[0]
            # dic2 = {"A":"0","B":"1","D":"2"}
            # num = int(b[3])-1
            # ch = int(dic2[b[4]])
            # count = int(pos/bin_size)
            # len_dic = defaultdict(list)
            #
            judge = 0
            for i in range(0,len(list1)):
                if pos >= int(list1[i][0]) and pos <= int(list1[i][1]):
                    judge=1
            for i in range(0,len(list2)):    
                if pos >= int(list2[i][0]) and pos <= int(list2[i][1]):
                    judge=2 
            for i in range(0,len(list3)):
                if pos >= int(list3[i][0]) and pos <= int(list3[i][1]):
                    judge=3  
            judge = 4
            ANN = tokens[4].split(",")
            GT1 = tokens[s1_col - 1]
            AD1 = tokens[s1_col]
            GQ1 = tokens[s1_col + 1]
            #print(GQ1)
            lis3 = ChangeGQ_dot(GQ1)
            GQ1 = lis3[0]
            lis1 = DefineAD(GT1, AD1)
            DP1 = lis1[0]
            GT1 = lis1[1]
            #print(tokens[s1_col - 1],lis1[1],pos,sep='\t')
            #
            GT2 = tokens[s2_col - 1]
            AD2 = tokens[s2_col]
            GQ2 = tokens[s2_col + 1]
            lis3_1 = ChangeGQ_dot(GQ2)
            GQ2 = lis3_1[0]
            lis1_1 = DefineAD(GT2, AD2)
            DP2 = lis1_1[0]
            GT2 = lis1_1[1]
            #
            # while pos > end:
            #     START = str(start)
            #     if START in list1:
            #         pass
            #     else:
            #         print("\t".join([b, "%s" % start, "%s" % end, "%s" % count1, "%s" % count2]))
            #     # #
            #     start += bin_size
            #     end += bin_size
            #     count1 = 0
                #
            
            if int(DP_low[0]) <= DP1 <= int(DP_high[0]) and int(DP_low[0]) <= DP2 <= int(DP_high[0]) and GQ1 >= int(GQ_sample) and GQ2 >= int(GQ_sample):    
                if IsGeno(GT1) and IsGeno(GT2) and IsHomo(GT1) and IsHomo(GT2):
                    if IsALT(GT1) or IsALT(GT2):
                        if Isunconsistent(GT1, GT2):
                            lst1=[]
                            lst2=[]
                            #print(pos)
                            count1 += 1
                            if tokens[4] != '.':
                                for ann in ANN:
                                    item = ann.split("|")
                                    FUN = item[1].strip() #.split("&")
                                    lst1.append([FUN,item[2],item[4]])
                                [lst2.append(i) for i in lst1 if not i in lst2]
                                    #lst2.sort(key=FUN.index)
                                for fun in lst2:
                                    #print(b,pos,ref,alt,GT1,GT2,"\t".join(item[0:3]),item[4],sep='\t')
                                    if judge ==1:
                                        print(b,pos,ref,alt,GT1,GT2,fun[0],fun[1],fun[2],'mid',sep='\t')
                                    elif judge ==2:
                                        print(b,pos,ref,alt,GT1,GT2,fun[0],fun[1],fun[2],'high',sep='\t')
                                    elif judge ==3:
                                        print(b,pos,ref,alt,GT1,GT2,fun[0],fun[1],fun[2],'low',sep='\t')
                                    elif judge ==4:
                                        print(b,pos,ref,alt,GT1,GT2,fun[0],fun[1],fun[2],sep='\t')
    #
    if grapht == 'on':       # not used now 2019/4/29
        firstline = IN.readline().strip()
        item = firstline.strip().split("\t")
        firstpos = int(item[pos_col - 1])
        if int(firstpos) <= int(bin_size):
            start = 1
            end = bin_size
        else:
            num = int(firstpos/bin_size)-1
            start = (num+1) * bin_size
            end = start + bin_size
        IN.seek(0, os.SEEK_SET)
        PART = infile.split("_")
        part = PART[0]
        art = part.lstrip('s')
        for line in IN:
            tokens = line.strip().split("\t")
            chr = tokens[chr_col - 1]
            pos = int(tokens[pos_col - 1])
            ref = tokens[ref_col - 1]
            alt = tokens[alt_col - 1]
            GT = tokens[gt_col - 1]
            DP = int(tokens[5])
            chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
            b = chr_cur[0]
            dic2 = {"A":"0","B":"1","D":"2"}
            num = int(b[3])-1
            ch = int(dic2[b[4]])
            count = int(pos/bin_size)
            #
            while pos > end:                        
                ratio = (count2+1)/(count1+1)
                DP_ratio = DP_count/count3
                if float(dic1[sample][(num)*3+ch][count]) >= 0.6:
                    print("\t".join([b, art + ".DP", "%s" % start, "%s" % end, "%s" % count1, "%s" % count2,
                                     "%.2f" % ratio, "RATIO", "%.2f" % DP_ratio, "%s" % count3,"normal"+ "%s" %count4]))  
                    print("\t".join([b, art, "%s" % start, "%s" % end, "%s" % count1, "%s" % count2,
                                     "%.2f" % DP_ratio, "DP", "%.2f" % ratio, "%s" % count3, "normal"+ "%s" %count4]))
                else:
                    count4 += 1
                    print("\t".join([b, art + ".DP", "%s" % start, "%s" % end, "%s" % count1, "%s" % count2,
                                     "%.2f" % ratio, "RATIO", "%.2f" % DP_ratio, "%s" % count3,"CNV"]))
                    print("\t".join([b, art, "%s" % start, "%s" % end, "%s" % count1, "%s" % count2,
                                     "%.2f" % DP_ratio, "DP", "%.2f" % ratio, "%s" % count3, "CNV"]))
                count1 = 0
                count2 = 0
                count3 = 0
                DP_count = 0
                start += bin_size
                end += bin_size
                #
                # compare genotype
            DP_count += DP
            count3 += 1 
            if IsHomo(GT):                 #continue 会跳出整个循环！！！！
                count1 += 1                #pass 不起作用 仅是起到站位作用
            if IsHete(GT):
                count2 += 1
    
        #
            #
        count += 1
        ratio = (count2+1)/(count1+1)
        DP_ratio = DP_count/count3
        if float(dic1[sample][(num)*3+ch][count]) >= 0.6:
            print("\t".join([b, art + ".DP", "%s" % start, "%s" % end, "%s" % count1, "%s" % count2,
                             "%.2f" % ratio, "RATIO", "%.2f" % DP_ratio, "%s" % count3,"normal"+ "%s" %count4]))
            print("\t".join([b, art, "%s" % start, "%s" % end, "%s" % count1, "%s" % count2,
                             "%.2f" % DP_ratio, "DP", "%.2f" % ratio, "%s" % count3, "normal"+ "%s" %count4]))
        else:
            print("\t".join([b, art + ".DP", "%s" % start, "%s" % end, "%s" % count1, "%s" % count2,
                             "%.2f" % ratio, "RATIO", "%.2f" % DP_ratio, "%s" % count3,"CNV"]))
            print("\t".join([b, art, "%s" % start, "%s" % end, "%s" % count1, "%s" % count2,
                             "%.2f" % DP_ratio, "DP", "%.2f" % ratio, "%s" % count3, "CNV"]))
    #
    if test == 'on':
        print(DP_high[1],"\t",DP_low[1])
    #
    if snp_graphy == 'on': #counting differnt snp number between two samples with alt genotype
        print("\t".join(['CHR', 'start', 'end', 'diff_homosnp_counts','same_homosnp_counts']))
        firstline = IN.readline().strip()
        item = firstline.strip().split("\t")
        firstpos = int(item[pos_col - 1])
        count_same=0
        count_diff=0
        count_same_num=0
        count_diff_num=0 
        if int(firstpos) <= int(bin_size):
            start = 1
            end = bin_size
        else:
            num = int(firstpos/bin_size)-1
            start = (num+1) * bin_size
            end = start + bin_size
        IN.seek(0, os.SEEK_SET)
        # PART = infile.split("_")
        # part = PART[0]
        # art = part.lstrip('s')
        list1 = Makelis(chromosome)
        #level = LEVEL.strip().split(",")
        DP_high = DP_high.strip().split(",")
        DP_low = DP_low.strip().split(",")
        for line in IN:
            tokens = line.strip().split("\t")
            chr = tokens[chr_col - 1]
            pos = int(tokens[pos_col - 1])
            ref = tokens[ref_col - 1]
            alt = tokens[alt_col - 1]
            chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
            b = chr_cur[0]
            dic2 = {"A":"0","B":"1","D":"2"}
            num = int(b[3])-1
            ch = int(dic2[b[4]])
            count = int(pos/bin_size)
            #
            len_dic = defaultdict(list)
            #
            GT1 = tokens[s1_col - 1]
            AD1 = tokens[s1_col]
            GQ1 = tokens[s1_col + 1]
            lis3 = ChangeGQ_dot(GQ1)
            GQ1 = lis3[0]
            lis1 = DefineAD(GT1, AD1)
            DP1 = lis1[0]
            GT1 = lis1[1]
            #
            GT2 = tokens[s2_col - 1]
            AD2 = tokens[s2_col]
            GQ2 = tokens[s2_col + 1]
            lis3_1 = ChangeGQ_dot(GQ2)
            GQ2 = lis3_1[0]
            lis1_1 = DefineAD(GT2, AD2)
            DP2 = lis1_1[0]
            GT2 = lis1_1[1]
            while pos > end:
                judge=0
                START = str(start)
                #deli = abs(count2-count1)
                for i in range(0,len(list1)):
                    if START==list1[i][0]:
                        #print("\t".join([b, "%s" % start, "%s" % end, "%s" % list1[i][2],"%s" % count_diff_num]))
                        #print("break1")
                        judge=1
                if judge==0:
                    print("\t".join([b, "%s" % start, "%s" % end,"%s" % count_diff_num,"%s" % count_same_num]))
                # #
                start += bin_size
                end += bin_size
                count1 = 0
                count1_1 = 0
                count2 = 0
                count2_1 = 0
                count3 = 0
                count4 = 0
                count_diff_num = 0
                count_same_num = 0

                #
            # if IsGeno(g1):
            #     count1 += 1
            # if IsGeno(g2):
            #     count2 += 1
            if int(DP_low[0]) <= DP1 <= int(DP_high[0]) and int(DP_low[0]) <= DP2 <= int(DP_high[0]) and GQ1 >= int(GQ_sample) and GQ2 >= int(GQ_sample):
                if IsGeno(GT1) and IsGeno(GT2) and IsHomo(GT1) and IsHomo(GT2):
                    if IsALT(GT1) or IsALT(GT2):
                        if Isunconsistent(GT1, GT2):
                            count_diff_num += 1
                        else:
                            count_same_num += 1
        #
        #deli = abs(count2-count1)
        #ratio = count3/deli
        START = str(start)
        judge=0
        #deli = abs(count2-count1)
        for i in range(0,len(list1)):
            if START==list1[i][0]:
                #print("break1")
                judge=1
        if judge==0:
            print("\t".join([b, "%s" % start, "%s" % end,"%s" % count_diff_num, "%s" % count_same_num]))
    #
    if density_graph_raw_vcf_two == 'on': # level differnt snp number between two samples with alt genotype
        print("\t".join(['CHR', 'start', 'end', 'diff_homosnp_ratio_level','diff_homosnp_count']))
        firstline = IN.readline().strip()
        item = firstline.strip().split("\t")
        firstpos = int(item[pos_col - 1])
        count_same=0
        count_diff=0
        count_same_num=0
        count_diff_num=0 
        if int(firstpos) <= int(bin_size):
            start = 1
            end = bin_size
        else:
            num = int(firstpos/bin_size)-1
            start = (num+1) * bin_size
            end = start + bin_size
        IN.seek(0, os.SEEK_SET)
        # PART = infile.split("_")
        # part = PART[0]
        # art = part.lstrip('s')
        list1 = Makelis(chromosome)
        level = LEVEL.strip().split(",")
        DP_high = DP_high.strip().split(",")
        DP_low = DP_low.strip().split(",")
        for line in IN:
            tokens = line.strip().split("\t")
            chr = tokens[chr_col - 1]
            pos = int(tokens[pos_col - 1])
            ref = tokens[ref_col - 1]
            alt = tokens[alt_col - 1]
            chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
            b = chr_cur[0]
            dic2 = {"A":"0","B":"1","D":"2"}
            num = int(b[3])-1
            #num = int(b[2])-1
            ch = int(dic2[b[4]])
            count = int(pos/bin_size)
            #
            len_dic = defaultdict(list)
            #
            GT1 = tokens[s1_col - 1]
            AD1 = tokens[s1_col]
            GQ1 = tokens[s1_col + 1]
            lis3 = ChangeGQ_dot(GQ1)
            GQ1 = lis3[0]
            lis1 = DefineAD(GT1, AD1)
            DP1 = lis1[0]
            GT1 = lis1[1]
            #
            GT2 = tokens[s2_col - 1]
            AD2 = tokens[s2_col]
            GQ2 = tokens[s2_col + 1]
            lis3_1 = ChangeGQ_dot(GQ2)
            GQ2 = lis3_1[0]
            lis1_1 = DefineAD(GT2, AD2)
            DP2 = lis1_1[0]
            GT2 = lis1_1[1]
            while pos > end:
                START = str(start)
                judge=0
                #deli = abs(count2-count1)
                for i in range(0,len(list1)):
                    if START==list1[i][0]:
                        print("\t".join([b, "%s" % start, "%s" % end, "%s" % list1[i][2],"%s" % count_diff_num]))
                        #print("break1")
                        judge=1
                if judge==0:
                    count_diff += count_diff_num
                    if count_diff_num == 0:
                        ratio =0
                    else:
                        ratio =  math.log(count_diff_num,10)
                    if ratio <= float(level[0]):
                        print("\t".join([b, "%s" % start, "%s" % end, "low_com","%s" % count_diff_num]))
                    if float(level[0]) < ratio <= float(level[1]):
                        print("\t".join([b, "%s" % start, "%s" % end, "mid_com","%s" % count_diff_num]))
                    if float(level[1]) < ratio :
                        print("\t".join([b, "%s" % start, "%s" % end, "high_com","%s" % count_diff_num]))
                # #
                start += bin_size
                end += bin_size
                count1 = 0
                count1_1 = 0
                count2 = 0
                count2_1 = 0
                count3 = 0
                count4 = 0
                count_diff_num = 0
                count_same_num = 0

                #
            # if IsGeno(g1):
            #     count1 += 1
            # if IsGeno(g2):
            #     count2 += 1
            if int(DP_low[0]) <= DP1 <= int(DP_high[0]) and int(DP_low[0]) <= DP2 <= int(DP_high[0]) and GQ1 >= int(GQ_sample) and GQ2 >= int(GQ_sample):
                if IsGeno(GT1) and IsGeno(GT2) and IsHomo(GT1) and IsHomo(GT2):
                    if IsALT(GT1) or IsALT(GT2):
                        if Isunconsistent(GT1, GT2):
                            count_diff_num += 1
                        else:
                            count_same_num += 1
        #
        #deli = abs(count2-count1)
        #ratio = count3/deli
        START = str(start)
        judge = 0
        for i in range(0,len(list1)):
            if start==list1[i][0]:
                print("\t".join([b, "%s" % start, "%s" % end, "%s" % list1[i][2],"%s" % count_diff_num]))
                judge=1
        if judge==0:
            count_diff += count_diff_num
            if count_diff_num == 0:
                ratio =0
            else:
                ratio =  math.log(count_diff_num,10)
            if ratio <= float(level[0]):
                print("\t".join([b, "%s" % start, "%s" % end, "low_com","%s" % count_diff_num]))
            if float(level[0]) < ratio <= float(level[1]):
                print("\t".join([b, "%s" % start, "%s" % end, "mid_com","%s" % count_diff_num]))
            if float(level[1]) < ratio :
                print("\t".join([b, "%s" % start, "%s" % end, "high_com","%s" % count_diff_num]))
    #
    if two_diff_level == 'on':
        firstline = IN.readline().strip()
        lis=[]
        for line in IN:
            tokens = line.strip().split("\t")
            chr = tokens[0]
            start = tokens[1]
            end = tokens[2]
            diff = tokens[3]
            same = tokens[4]
            lis.append([chr,start,end,diff,same])
        frame = pd.DataFrame(lis)
        frame[3] =  frame[3].astype(int)
        frame[1] =  frame[1].astype(int)
        dic={16:[0.0005181483,0.003151449],17:[0.0007397899,0.00249276],18:[0.001025806,0.001988525],19:[0.001386473,0.001598718],20:[0.001832182,0.001294632],21:[ 0.002373268,0.001055425],22:[0.003019867,0.0008657901],23:[0.003781773,0.0007143668],24:[0.004668331, 0.0005926393]}
        frame.loc[(20 < frame[3]) ,5] = 'diff'
        frame.loc[(20 >= frame[3]) ,5] = 'same'
        for num in range(16,25):
            for i in frame.loc[(num == frame[3]) ,1]:
                bayes_frame = frame.loc[((i-5000000) <= frame[1])&(frame[1] <= (i+6000000))&(frame[1] != i),]
                if bayes_frame[(bayes_frame[5]=='diff')].empty :
                    len_diff = 0
                else:
                    len_diff = len(bayes_frame[(bayes_frame[5]=='diff')])
                if bayes_frame[(bayes_frame[5]=='same')].empty :
                    len_same = 0
                else:
                    len_same = len(bayes_frame[(bayes_frame[5]=='same')])
                if len_diff*dic[17][0] > len_same*dic[17][1]:
                    frame.loc[(frame[1] == i),5]='diff'
                else:
                    frame.loc[(frame[1] == i),5]='same'
        #diff_file = "tmp/" + total + '/' +sample1 + '_' +sample2 + "_" + time_data + "/" + chr + ".homo_snp_level"
        #frame.to_csv(diff_file,sep='\t',na_rep='Na',float_format='%.1f',header=['CHR', 'start', 'end', 'diff_homosnp_count', 'same_homosnp_count','diff_homosnp_ratio_level'],index=False)
        diff_file = "tmp/" + total + '/' +sample1 + '_' +sample2 + "_" + time_data + "/" + chr + ".only_homo_snp_level"
        frame.to_csv(diff_file,sep='\t',na_rep='Na',float_format='%.1f', columns=[0,1,2,5], header=False,index=False)
        

    if per_snp_filter == 'on':  # print snp information with alt (first step)
        firstline = IN.readline().strip()
        item = firstline.strip().split("\t")
        firstpos = int(item[pos_col - 1])
        if int(firstpos) <= int(bin_size):
            start = 1
            end = bin_size
        else:
            num = int(firstpos/bin_size)-1
            start = (num+1) * bin_size
            end = start + bin_size
        IN.seek(0, os.SEEK_SET)
        PART = infile.split("_")
        part = PART[0]
        art = part.lstrip('s')
        DP_high = DP_high.strip().split(",")
        DP_low = DP_low.strip().split(",")
        list1 = Makelis(chromosome)
        for line in IN:
            tokens = line.strip().split("\t")
            chr = tokens[chr_col - 1]
            pos = int(tokens[pos_col - 1])
            chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
            b = chr_cur[0]
            judge = 0
            for i in range(0,len(list1)):
                if pos >= int(list1[i][0]) and pos <= int(list1[i][1]):
                    judge=1         
            if judge == 0:
                GT1 = tokens[s1_col - 1]
                AD1 = tokens[s1_col]
                GQ1 = tokens[s1_col + 1]
                lis3 = ChangeGQ_dot(GQ1)
                GQ1 = lis3[0]
                lis1 = DefineAD(GT1, AD1)
                DP1 = lis1[0]
                GT1 = lis1[1]
                #
                GT2 = tokens[s2_col - 1]
                AD2 = tokens[s2_col]
                GQ2 = tokens[s2_col + 1]
                lis3_1 = ChangeGQ_dot(GQ2)
                GQ2 = lis3_1[0]
                lis1_1 = DefineAD(GT2, AD2)
                DP2 = lis1_1[0]
                GT2 = lis1_1[1]

                if IsGeno(GT1) and IsGeno(GT2) and IsHomo(GT1) and IsHomo(GT2):
                    if int(DP_low[0]) <= DP1 <= int(DP_high[0]) and int(DP_low[0]) <= DP2 <= int(DP_high[0]) and GQ1 >= int(GQ_sample) and GQ2 >= int(GQ_sample):
                        if IsALT(GT1) or IsALT(GT2):
                            print(line)
    #
    if filter_and_count_snp == 'on':  # filter on sample snp and count homo_alt snp 
        dic1 = Makedic(sample)
        firstline = IN.readline().strip()
        item = firstline.strip().split("\t")
        firstpos = int(item[pos_col - 1])
        if int(firstpos) <= int(bin_size):
            start = 1
            end = bin_size
        else:
            num = int(firstpos/bin_size)-1
            start = (num+1) * bin_size
            end = start + bin_size
        IN.seek(0, os.SEEK_SET)
        PART = infile.split("_")
        part = PART[0]
        art = part.lstrip('s')
        DP_high = DP_high.strip().split(",")
        DP_low = DP_low.strip().split(",")
        print("\t".join(['CHR', 'start', 'end', 'ALT_HOMO', 'HETE']))
        for line in IN:
            tokens = line.strip().split("\t")
            chr = tokens[chr_col - 1]
            pos = int(tokens[pos_col - 1])
            ref = tokens[ref_col - 1]
            alt = tokens[alt_col - 1]
            chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
            b = chr_cur[0]
            dic2 = {"A":"0","B":"1","D":"2"}
            num = int(b[3])-1
            ch = int(dic2[b[4]])
            count = int(pos/bin_size)
            #
            len_dic = defaultdict(list)
            a = 0
            GT = tokens[s1_col - 1]
            AD = tokens[s1_col]
            GQ = tokens[s1_col + 1]
            lis3 = ChangeGQ_dot(GQ)
            GQ = lis3[0]
            lis1 = DefineAD(GT, AD)
            DP = lis1[0]
            GT = lis1[1]
            judge_num = 0
            while pos > end:

                if float(dic1[sample][(num)*3+ch][count]) < 0.5:
                    pass
                else:
                    print("\t".join([b, "%s" % start, "%s" % end, "%s" % count1, "%s" % count2]))
                # #
                start += bin_size
                end += bin_size
                count1 = 0
                count2 = 0

                judge_num = 1
            if int(DP_low[0]) <= DP <= int(DP_high[0]) and GQ >= int(GQ_sample):
                if IsGeno(GT):
                    if IsHomo(GT) :
                        if IsALT(GT):
                            count1 += 1
                    else:
                        count2 += 1
            
        if float(dic1[sample][(num)*3+ch][count]) < 0.5:
                    pass
        else:
            print("\t".join([b, "%s" % start, "%s" % end, "%s" % count1, "%s" % count2]))
    #
    if density_graph_raw_vcf_one == 'on': # level differnt snp number between two samples with alt genotype
        firstline = IN.readline().strip()
        item = firstline.strip().split("\t")
        firstpos = int(item[pos_col - 1])
        if int(firstpos) <= int(bin_size):
            start = 1
            end = bin_size
        else:
            num = int(firstpos/bin_size)-1
            start = (num+1) * bin_size
            end = start + bin_size
        IN.seek(0, os.SEEK_SET)
        PART = infile.split("_")
        part = PART[0]
        art = part.lstrip('s')
        DP_high = DP_high.strip().split(",")
        DP_low = DP_low.strip().split(",")
        level = LEVEL.strip().split(",")
        list1 = Makelis(chromosome)
        print("\t".join(['CHR', 'start', 'end', 'alt_ratio_level','alt_counts']))
        for line in IN:
            tokens = line.strip().split("\t")
            chr = tokens[chr_col - 1]
            pos = int(tokens[pos_col - 1])
            ref = tokens[ref_col - 1]
            alt = tokens[alt_col - 1]
            chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
            b = chr_cur[0]
            dic2 = {"A":"0","B":"1","D":"2"}
            num = int(b[3])-1
            ch = int(dic2[b[4]])
            count = int(pos/bin_size)
            #
            a = 0
            GT = tokens[s1_col - 1]
            AD = tokens[s1_col]
            GQ = tokens[s1_col + 1]
            lis3 = ChangeGQ_dot(GQ)
            GQ = lis3[0]
            lis1 = DefineAD(GT, AD)
            GT = lis1[1]
            DP = lis1[0]
            judge_num = 0
            while pos > end:
                START = str(start)
                judge=0
                #deli = abs(count2-count1)
                for i in range(0,len(list1)):
                    if START==list1[i][0]:
                        #print("\t".join([b, "%s" % start, "%s" % end, "%s" % list1[i][2],"%s" % count1]))
                        #print("break1")
                        judge=1
                if judge==0:
                    if count1 == 0:
                        print("\t".join([b, "%s" % start, "%s" % end, 'low',"%s" % count1]))
                    else:
                        ratio = math.log(count1,10)
                        if ratio >= float(level[0]):
                            print("\t".join([b, "%s" % start, "%s" % end, 'high',"%s" % count1]))
                        else:
                            print("\t".join([b, "%s" % start, "%s" % end, 'low',"%s" % count1]))
                # #
                start += bin_size
                end += bin_size
                count1 = 0
                count2 = 0

                judge_num = 1
            if IsGeno(GT) and IsHomo(GT) :
                if int(DP_low[0]) <= DP <= int(DP_high[0]) and GQ >= int(GQ_sample):
                    if IsALT(GT):
                        count1 += 1
                    count2 += 1
        START = str(start)
        judge=0
        #deli = abs(count2-count1)
        for i in range(0,len(list1)):
            if START==list1[i][0]:
                #print("\t".join([b, "%s" % start, "%s" % end, "%s" % list1[i][2],"%s" % count1]))
                #print("break1")
                judge=1
        if judge==0:
            if count1 == 0:
                print("\t".join([b, "%s" % start, "%s" % end, 'low',"%s" % count1]))
            else:
                ratio = math.log(count1,10)
                if ratio >= float(level[0]):
                    print("\t".join([b, "%s" % start, "%s" % end, 'high',"%s" % count1]))
                else:
                    print("\t".join([b, "%s" % start, "%s" % end, 'low',"%s" % count1]))
    #
    if DP_GQ == 'on':    # count every DP number and GQ number; see in muti_func_snp_filter.py
        f = open(chromosome, "r")
        mask = f.readlines()
        f.close
        a = 0
        dic1 = {}
        for line in IN:
            tokens = line.strip().split("\t")
            chr = tokens[chr_col - 1]
            pos = int(tokens[pos_col - 1])
            g1 = tokens[s1_col - 1]
            useitem = int(tokens[s2_col-1])
            chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
            b = chr_cur[0]
            judge_num = 0
            for ther in mask:
                item = ther.strip("").split("\t")
                if item[1] <= str(pos) <= item[2]:
                    judge_num += 1       
            if judge_num == 0:
                if IsGeno(g1) and IsHomo(g1):
                    dic1[useitem] = dic1.setdefault(useitem, 0)+1
        for key,value in dic1.items():
            print('{key}\t{value}'.format(key=key, value=value))
        #
    #
    if mask_cnv == 'on':
        readline = IN.readlines()[1:]
        start = 1
        end = 1000001
        for line in readline:
            tokens = line.strip().split("\t")
            RD = float(tokens[0])
            filename1 = infile.strip().split("/")
            filename2 = filename1[-1].strip().split(".")
            filename3 = filename2[0]
            if RD <= 0.5:                
                print(filename3,start,end,"ownCNV_deletion",sample1,sep="\t")
            if RD >= 1.5:                
                print(filename3,start,end,"ownCNV_duplication",sample1,sep="\t")
            start += bin_size
            end += bin_size
    #
    if cnv_region == 'on':
        #print("chromosome","CNV_length")
        filename1 = infile.strip().split("/")
        filename2 = filename1[-1].strip().split(".")
        filename3 = filename2[0]
        centro = centrosome.strip().split(",")
        dic2 = {"A":"0","B":"1","D":"2"}
        ind = (int(filename3[3])-1)*3 + int(dic2[filename3[4]])
        for line in IN:
            tokens = line.strip().split("\t")
            end = int(tokens[2])            
            if int(end) <= int(centro[ind])*1000000:                
                count1 += 1
            else:
                count2 += 1
        #print(length)
        if 2*int(centro[ind])*1000000 > int(length[filename3]):
            print(filename3+"L",count1,sep='\t')
            print(filename3+"S",count2,sep='\t')
        else:
            print(filename3+"L",count2,sep='\t')
            print(filename3+"S",count1,sep='\t')                
    #
    if transition_transversion_count == 'on':
        dic1 = {}
        dic2 = {}
        dic3 = {}
        dic4 = {}
        list1 = []
        list2 = []
        list3 = []
        f = open(chromosome, "r")
        mask = f.readlines()
        for line in mask:
            item = line.strip().split("\t")
            if str(item[3]) == 'mid_com':
                list1.append(item[1])
            if str(item[3]) == 'high_com':
                list2.append(item[1])
            if str(item[3]) == 'low_com':
                list3.append(item[1])            
        f.close
        firstline = IN.readline().strip()
        item = firstline.strip().split("\t")
        firstpos = int(item[pos_col - 1])
        if int(firstpos) <= int(bin_size):
            start = 1
            end = bin_size
        else:
            num = int(firstpos/bin_size)-1
            start = (num+1) * bin_size
            end = start + bin_size
        IN.seek(0, os.SEEK_SET)
        # PART = infile.split("_")
        # part = PART[0]
        # art = part.lstrip('s')
        DP_high = DP_high.strip().split(",")
        DP_low = DP_low.strip().split(",")
        for line in IN:
            tokens = line.strip().split("\t")
            chr = tokens[chr_col - 1]
            pos = int(tokens[pos_col - 1])
            ref = tokens[ref_col - 1]
            alt = tokens[alt_col - 1].strip().split(",")
            chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
            b = chr_cur[0]
            #
            GT1 = tokens[s1_col - 1]
            AD1 = tokens[s1_col]
            GQ1 = tokens[s1_col + 1]
            lis3 = ChangeGQ_dot(GQ1)
            GQ1 = lis3[0]
            lis1 = DefineAD(GT1, AD1)
            DP1 = lis1[0]
            GT1 = lis1[1]
            #
            GT2 = tokens[s2_col - 1]
            AD2 = tokens[s2_col]
            GQ2 = tokens[s2_col + 1]
            lis3_1 = ChangeGQ_dot(GQ2)
            GQ2 = lis3_1[0]
            lis1_1 = DefineAD(GT2, AD2)
            DP2 = lis1_1[0]
            GT2 = lis1_1[1]
            while pos > end:
                START = str(start)
                if START in list1:
                    for k, v in dic1.items():
                        if k in dic2.keys():
                            dic2[k] += v
                        else:
                            dic2[k] = v  
                if START in list2:   
                    for k, v in dic1.items():
                        if k in dic3.keys():
                            dic3[k] += v
                        else:
                            dic3[k] = v 
                if START in list3:
                    count1_1 += count1
                    for k, v in dic1.items():
                        if k in dic4.keys():
                            dic4[k] += v
                        else:
                            dic4[k] = v      
                dic1 = {}  
                count1 = 0
                start += bin_size
                end += bin_size                         
            if int(DP_low[0]) <= DP1 <= int(DP_high[0]) and int(DP_low[0]) <= DP2 <= int(DP_high[0]) and GQ1 >= int(GQ_sample) and GQ2 >= int(GQ_sample):
                if IsGeno(GT1) and IsGeno(GT2) and IsHomo(GT1) and IsHomo(GT2):
                    if IsALT(GT1) or IsALT(GT2):
                            if Isunconsistent(GT1,GT2):
                                count1 += 1
                                if '2' in GT1 or '2' in GT2:
                                    alt_use = alt[1]
                                    ref_alt = ref + alt_use
                                    dic1["SNP_num"] = count1
                                    dic1[ref_alt] = dic1.setdefault(ref_alt, 0)+1 
                                    
                                    
                                else:
                                    alt_use = alt[0]
                                    ref_alt = ref + alt_use
                                    dic1["SNP_num"] = count1
                                    dic1[ref_alt] = dic1.setdefault(ref_alt, 0)+1 
        START = str(start)
        if START in list1:
            for k, v in dic1.items():
                if k in dic2.keys():
                    dic2[k] += v
                else:
                    dic2[k] = v  
        if START in list2:   
            for k, v in dic1.items():
                if k in dic3.keys():
                    dic3[k] += v
                else:
                    dic3[k] = v 
        if START in list3:
            count1_1 += count1
            for k, v in dic1.items():
                if k in dic4.keys():
                    dic4[k] += v
                else:
                    dic4[k] = v                 
                            #print(dic1)
        for key,value in dic4.items():
          print('{key}\t{value}'.format(key=key, value=value), 'low', sep='\t')
        for key,value in dic2.items():
          print('{key}\t{value}'.format(key=key, value=value), 'mid', sep='\t')
        for key,value in dic3.items():
          print('{key}\t{value}'.format(key=key, value=value), 'high', sep='\t')
        
        #temporary
        # for key,value in dic1.items():
        #     print('{key}\t{value}'.format(key=key, value=value), 'high', sep='\t')
        
        #
    #
    if indel_count == 'on':
        dic1 = {}
        dic2 = {}
        dic3 = {}
        dic4 = {}
        list1 = []
        list2 = []
        list3 = []
        f = open(chromosome, "r")
        mask = f.readlines()
        for line in mask:
            item = line.strip().split("\t")
            if str(item[3]) == 'mid_com':
                list1.append(item[1])
            if str(item[3]) == 'high_com':
                list2.append(item[1])
            if str(item[3]) == 'low_com':
                list3.append(item[1]) 
        f.close
        firstline = IN.readline().strip()
        item = firstline.strip().split("\t")
        firstpos = int(item[pos_col - 1])
        if int(firstpos) <= int(bin_size):
            start = 1
            end = bin_size
        else:
            num = int(firstpos/bin_size)-1
            start = (num+1) * bin_size
            end = start + bin_size
        IN.seek(0, os.SEEK_SET)
        # PART = infile.split("_")
        # part = PART[0]
        # art = part.lstrip('s')
        DP_high = DP_high.strip().split(",")
        DP_low = DP_low.strip().split(",")
        for line in IN:
            tokens = line.strip().split("\t")
            chr = tokens[chr_col - 1]
            pos = int(tokens[pos_col - 1])
            ref = tokens[ref_col - 1]
            alt = tokens[alt_col - 1].strip().split(",")
            list_alt_ref=[ref]+alt
            #print(list_alt_ref,list_alt_ref[0],list_alt_ref[1])
            chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
            b = chr_cur[0]
            #
            GT1 = tokens[s1_col - 1]
            AD1 = tokens[s1_col]
            GQ1 = tokens[s1_col + 1]
            lis3 = ChangeGQ_dot(GQ1)
            GQ1 = lis3[0]
            lis1 = DefineAD(GT1, AD1)
            DP1 = lis1[0]
            GT1 = lis1[1]
            #
            GT2 = tokens[s2_col - 1]
            AD2 = tokens[s2_col]
            GQ2 = tokens[s2_col + 1]
            lis3_1 = ChangeGQ_dot(GQ2)
            GQ2 = lis3_1[0]
            lis1_1 = DefineAD(GT2, AD2)
            DP2 = lis1_1[0]
            GT2 = lis1_1[1]
            while pos > end:
                START = str(start)
                if START in list1:
                    for k, v in dic1.items():
                        if k in dic2.keys():
                            dic2[k] += v
                        else:
                            dic2[k] = v  
                if START in list2:   
                    for k, v in dic1.items():
                        if k in dic3.keys():
                            dic3[k] += v
                        else:
                            dic3[k] = v 
                if START in list3:
                    count1_1 += count1
                for k, v in dic1.items():
                    if k in dic4.keys():
                        dic4[k] += v
                    else:
                        dic4[k] = v      
                dic1 = {}  
                count1 = 0
                start += bin_size
                end += bin_size                         
            if int(DP_low[0]) <= DP1 <= int(DP_high[0]) and int(DP_low[0]) <= DP2 <= int(DP_high[0]) and GQ1 >= int(GQ_sample) and GQ2 >= int(GQ_sample):
                if IsGeno(GT1) and IsGeno(GT2) and IsHomo(GT1) and IsHomo(GT2):
                    if IsALT(GT1) or IsALT(GT2):
                            count1 += 1
                            if Isunconsistent(GT1,GT2) and list_alt_ref[int(GT1[0])] != '*' and list_alt_ref[int(GT2[0])] != '*':
                                indel_len=len(list_alt_ref[int(GT1[0])])-len(list_alt_ref[int(GT2[0])])
                                dic1[indel_len] = dic1.setdefault(indel_len, 0)+1
                                # if indel_len == 0:
                                #     print(line)
                                #     print(GT1,GT2,list_alt_ref,list_alt_ref[int(GT1[0])],list_alt_ref[int(GT2[0])])
                                # if '2' in GT1 or '2' in GT2:
                                #     alt_use = alt[1]
                                #     indel_len=len(alt_use)-len(ref)
                                #     dic1[indel_len] = dic1.setdefault(indel_len, 0)+1 
                                # else:
                                #     alt_use = alt[0]
                                #     indel_len=len(alt_use)-len(ref)
                                #     dic1[indel_len] = dic1.setdefault(indel_len, 0)+1 

                            #print(dic1)
        START = str(start)
        if START in list1:
            for k, v in dic1.items():
                if k in dic2.keys():
                    dic2[k] += v
                else:
                    dic2[k] = v  
        if START in list2:   
            for k, v in dic1.items():
                if k in dic3.keys():
                    dic3[k] += v
                else:
                    dic3[k] = v 
        if START in list3:
            count1_1 += count1
            for k, v in dic1.items():
                if k in dic4.keys():
                    dic4[k] += v
                else:
                    dic4[k] = v 
        for key,value in dic4.items():
            print('{key}\t{value}'.format(key=key, value=value), 'low', sep='\t')
        for key,value in dic2.items():
            print('{key}\t{value}'.format(key=key, value=value), 'mid', sep='\t')
        for key,value in dic3.items():
            print('{key}\t{value}'.format(key=key, value=value), 'high', sep='\t')
        
        #temporary
        #for key,value in dic1.items():
        #    print('{key}\t{value}'.format(key=key, value=value), 'high', sep='\t')    
    #
    if filter_and_count_indel == 'on':
        firstline = IN.readline().strip()
        item = firstline.strip().split("\t")
        firstpos = int(item[pos_col - 1])
        count_indel_length=0      
        print("\t".join(["CHR", "start", "end", "indel_count", "indel_length"]))
        if int(firstpos) <= int(bin_size):
            start = 1
            end = bin_size
        else:
            num = int(firstpos/bin_size)-1
            start = (num+1) * bin_size
            end = start + bin_size
        IN.seek(0, os.SEEK_SET)
        PART = infile.split("_")
        part = PART[0]
        art = part.lstrip('s')
        dic1 = Makedic(sample)
        #print(list1)
        DP_high = DP_high.strip().split(",")
        DP_low = DP_low.strip().split(",")
        for line in IN:
            #
            tokens = line.strip().split("\t")
            chr = tokens[chr_col - 1]
            pos = int(tokens[pos_col - 1])
            ref = tokens[ref_col - 1]
            alt = tokens[alt_col - 1]
            chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
            b = chr_cur[0]
            dic2 = {"A":"0","B":"1","D":"2"}
            num = int(b[3])-1
            ch = int(dic2[b[4]])
            count = int(pos/bin_size)
            list_alt_ref=[ref,tokens[alt_col - 1]]
            #
            len_dic = defaultdict(list)
            a = 0
            GT = tokens[s1_col - 1]
            AD = tokens[s1_col]
            GQ = tokens[s1_col + 1]
            lis3 = ChangeGQ_dot(GQ)
            GQ = lis3[0]
            lis1 = DefineAD(GT, AD)
            DP = lis1[0]
            GT = lis1[1]
            judge_num = 0
            while pos > end:
                START = str(start)
                #deli = abs(count2-count1)
                #ratio = count3/deli
                if float(dic1[sample][(num)*3+ch][count]) < 0.5:
                    pass
                else:
                    print("\t".join([b, "%s" % start, "%s" % end, "%s" % count1,"%s" % count_indel_length]))

                # #
                start += bin_size
                end += bin_size
                count1 = 0
                count_indel_length = 0
                #
            if int(DP_low[0]) <= DP <= int(DP_high[0]) and  GQ >= int(GQ_sample):
                if IsGeno(GT) and IsHomo(GT):
                    if IsALT(GT):
                        count1 += 1
                        indel_len=abs(len(ref)-len(list_alt_ref[int(GT[0])]))
                        count_indel_length += indel_len

        #
        #deli = abs(count2-count1)
        #ratio = count3/deli
        START = str(start)
        if float(dic1[sample][(num)*3+ch][count]) < 0.5:
            pass
        else:
            print("\t".join([b, "%s" % start, "%s" % end, "%s" % count1,"%s" % count_indel_length]))
    #
    if parents_define == 'on':
        print("\t".join(['CHR', 'start', 'end', 'diff_homosnp_ratio_level','diff_homosnp_count']))
        firstline = IN.readline().strip()
        item = firstline.strip().split("\t")
        firstpos = int(item[pos_col - 1])
        count_same=0
        count_diff=0
        count_same_num=0
        count_diff_num=0 
        if int(firstpos) <= int(bin_size):
            start = 1
            end = bin_size
        else:
            num = int(firstpos/bin_size)-1
            start = (num+1) * bin_size
            end = start + bin_size
        IN.seek(0, os.SEEK_SET)
        # PART = infile.split("_")
        # part = PART[0]
        # art = part.lstrip('s')
        list1 = Makelis(chromosome)
        level = LEVEL.strip().split(",")
        DP_high = DP_high.strip().split(",")
        DP_low = DP_low.strip().split(",")
        for line in IN:
            tokens = line.strip().split("\t")
            chr = tokens[chr_col - 1]
            pos = int(tokens[pos_col - 1])
            ref = tokens[ref_col - 1]
            alt = tokens[alt_col - 1]
            chr_cur = chr.strip().split(".")  # 不用特意去掉.1 .2 
            b = chr_cur[0]
            dic2 = {"A":"0","B":"1","D":"2"}
            num = int(b[3])-1
            ch = int(dic2[b[4]])
            count = int(pos/bin_size)
            #
            len_dic = defaultdict(list)
            #
            GT1 = tokens[s1_col - 1]
            AD1 = tokens[s1_col]
            GQ1 = tokens[s1_col + 1]
            lis3 = ChangeGQ_dot(GQ1)
            GQ1 = lis3[0]
            lis1 = DefineAD(GT1, AD1)
            DP1 = lis1[0]
            GT1 = lis1[1]
            #
            GT2 = tokens[s2_col - 1]
            AD2 = tokens[s2_col]
            GQ2 = tokens[s2_col + 1]
            lis3_1 = ChangeGQ_dot(GQ2)
            GQ2 = lis3_1[0]
            lis1_1 = DefineAD(GT2, AD2)
            DP2 = lis1_1[0]
            GT2 = lis1_1[1]
            #
            GTs = tokens[s3_col - 1]
            ADs = tokens[s3_col]
            GQs = tokens[s3_col + 1]
            lis3_2 = ChangeGQ_dot(GQs)
            GQs = lis3_2[0]
            lis1_2 = DefineAD(GTs, ADs)
            DPs = lis1_2[0]
            GTs = lis1_2[1]

            while pos > end:
                START = str(start)
                judge=0
                #deli = abs(count2-count1)
                for i in range(0,len(list1)):
                    if START==list1[i][0]:
                        #print("\t".join([b, "%s" % start, "%s" % end, "%s" % list1[i][2],"%s" % count_diff_num]))
                        #print("break1")
                        judge=1
                if judge==0:
                    count_diff_gt1_gt2 += gt1_gts_diff_num 
                    count_diff_gt2_gt3 += gt2_gts_diff_num 
                    print("\t".join([b, "%s" % start, "%s" % end, "%s" % gt1_gts_diff_num, "%s" % gt2_gts_diff_num]))
                    #if count_diff_num == 0:
                    #    ratio =0
                    #else:
                    #    ratio =  math.log(count_diff_num,10)
                    #if ratio <= float(level[0]):
                    #    print("\t".join([b, "%s" % start, "%s" % end, "low_com","%s" % count_diff_num]))
                    #if float(level[0]) < ratio <= float(level[1]):
                    #    print("\t".join([b, "%s" % start, "%s" % end, "mid_com","%s" % count_diff_num]))
                    #if float(level[1]) < ratio :
                    #    print("\t".join([b, "%s" % start, "%s" % end, "high_com","%s" % count_diff_num]))
                # #
                start += bin_size
                end += bin_size
                count1 = 0
                count1_1 = 0
                count2 = 0
                count2_1 = 0
                count3 = 0
                count4 = 0
                gt1_gts_diff_num  = 0
                gt2_gts_diff_num  = 0

                #
            # if IsGeno(g1):
            #     count1 += 1
            # if IsGeno(g2):
            #     count2 += 1
            if (int(DP_low[0]) <= DP1 <= int(DP_high[0]) 
            and int(DP_low[0]) <= DP2 <= int(DP_high[0]) 
            and int(DP_low[0]) <= DPs <= int(DP_high[0]) 
            and GQ1 >= int(GQ_sample) 
            and GQ2 >= int(GQ_sample)
            and GQs >= int(GQ_sample)):
                if (IsGeno(GT1) and IsGeno(GT2) and IsGeno(GTs) 
                and IsHomo(GT1) and IsHomo(GT2) and IsHomo(GTs)):
                    if GT1 != GTs and GT2 != GTs:
                        pass
                    else:
                        if Isunconsistent(GT1, GTs):
                            gt1_gts_diff_num += 1
                        if Isunconsistent(GT2, GTs):
                            gt2_gts_diff_num += 1
        #
        #deli = abs(count2-count1)
        #ratio = count3/deli
        START = str(start)
        for i in range(0,len(list1)):
            if START==list1[i][0]:
                #print("\t".join([b, "%s" % start, "%s" % end, "%s" % list1[i][2],"%s" % count_diff_num]))
                #print("break1")
                judge=1
        if judge==0:
            count_diff_gt1_gt2 += gt1_gts_diff_num 
            count_diff_gt2_gt3 += gt2_gts_diff_num 
            print("\t".join([b, "%s" % start, "%s" % end, "%s" % gt1_gts_diff_num, "%s" % gt2_gts_diff_num]))



    if infile:
        IN.close()
    #
#



from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: %prog [-i <input>] [-b 1000] [-o <output>]\n" \
            "Author : Yang, zhengzhao; yangzhengzhao@cau.edu.cn; 2018-06-07\n" \
            "Description: Identify genotype bins orientation from parents.\n" \
            "Input format:\n" \
            "   chr pos geno_p1 geno_p2 geno_son\n" \
            "Output format:\n" \
            "   chr start end total eqp1 eqp2 orient\n" \
            "   chr1A   10001 20000 20  18  2   1\n"
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="infile",
                  help="Input file, use STDIN if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
    parser.add_option("-b", dest="bin_size",
                  help="size of bin (bp) [default: %default]", metavar="INT",
                      default = 10000)                  
    parser.add_option("--gene_snp_filter", dest="gene_snp_filter",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--chrom_graphy", dest="chrom_graphy",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--density_graph_raw_vcf_two", dest="density_graph_raw_vcf_two",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--density", dest="density",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--two_diff_level", dest="two_diff_level",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--sample1", dest="sample1",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'A')
    parser.add_option("--sample2", dest="sample2",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'A')
    parser.add_option("--time_data", dest="time_data",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'time_data')
    parser.add_option("--total", dest="total",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'total')
    parser.add_option("--snp_graphy", dest="snp_graphy",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--cnv_region", dest="cnv_region",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--mask_DP", dest="mask_DP",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--parents_define", dest="parents_define",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')              
    parser.add_option("--filter_and_count_indel", dest="filter_and_count_indel",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--long_short_count", dest="long_short_count",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--two_diff", dest="two_diff",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--grapht", dest="grapht",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off') 
    parser.add_option("--per_snp_filter", dest="per_snp_filter",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--mask_cnv", dest="mask_cnv",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--transition_transversion_count", dest="transition_transversion_count",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--indel_count", dest="indel_count",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--high_threshold", dest="high_threshold",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = 1)
    parser.add_option("--DP_high", dest="DP_high",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = 1)
    parser.add_option("--DP_low", dest="DP_low",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = 1)
    parser.add_option("--GQ_sample", dest="GQ_sample",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = 1)
    parser.add_option("--LEVEL", dest="LEVEL",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = 1)
    parser.add_option("--centrosome", dest="centrosome",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = "215,239,173,337,346,268,300,346,241,300,317,184,254,206,189,286,327,211,357,287,340")                
    parser.add_option("--test", dest="test",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')   
    parser.add_option("--density_graph_raw_vcf_one", dest="density_graph_raw_vcf_one",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--filter_and_count_snp", dest="filter_and_count_snp",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')                 
    parser.add_option("--DP_GQ", dest="DP_GQ",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')                                                
    parser.add_option("-c", dest="chr_col",
                  help="column id for chromosome [default: %default]", metavar="INT",
                      default = 1)
    parser.add_option("-p", dest="pos_col",
                  help="column id for position [default: %default]", metavar="INT",
                      default = 2)
    parser.add_option("-r", dest="ref_col",
                  help="column id for position [default: %default]", metavar="INT",
                      default = 3)
    parser.add_option("-a", dest="alt_col",
                  help="column id for position [default: %default]", metavar="INT",
                      default = 4)
    parser.add_option("-1", dest="s1_col",
                  help="column id for 1st parent [default: %default]", metavar="INT",
                      default = 5)
    parser.add_option("-2", dest="s2_col",
                  help="column id for 2nd parent [default: %default]", metavar="INT",
                      default = 8)
    parser.add_option("-s", dest="s3_col",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = 6)
    parser.add_option("--sample", dest="sample",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'S1')
    parser.add_option("--chromosome", dest="chromosome",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = '1B')
    parser.add_option("-o", dest="outfile",
                  help="Output file, use STDOUT if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
    #
    (options, args) = parser.parse_args()
    #
    if (options.outfile is not None) :  #??????
        if options.outfile.endswith('.gz') :
            sys.stdout = gzip.open(options.outfile, 'wb')
        else:
            sys.stdout = open(options.outfile, 'w')
        #
    #
    GenoPhase(infile=options.infile, bin_size=int(options.bin_size), gene_snp_filter=options.gene_snp_filter, 
              chrom_graphy=options.chrom_graphy, density_graph_raw_vcf_one=options.density_graph_raw_vcf_one, density=options.density, snp_graphy=options.snp_graphy, 
              mask_DP=options.mask_DP, long_short_count=options.long_short_count, mask_cnv=options.mask_cnv,two_diff=options.two_diff, grapht=options.grapht, 
              per_snp_filter=options.per_snp_filter, DP_GQ=options.DP_GQ, filter_and_count_snp=options.filter_and_count_snp, indel_count=options.indel_count,
              high_threshold=int(options.high_threshold), DP_high=options.DP_high, GQ_sample=options.GQ_sample, two_diff_level=options.two_diff_level,
              sample1=options.sample1, sample2=options.sample2, time_data=options.time_data, total=options.total,
              DP_low=options.DP_low, LEVEL=options.LEVEL, centrosome=options.centrosome, test=options.test, parents_define=options.parents_define,
              sample=options.sample, chromosome=options.chromosome, density_graph_raw_vcf_two=options.density_graph_raw_vcf_two, 
              chr_col=int(options.chr_col), pos_col=int(options.pos_col), cnv_region=options.cnv_region, filter_and_count_indel=options.filter_and_count_indel,
              ref_col=int(options.ref_col), alt_col=int(options.alt_col), transition_transversion_count=options.transition_transversion_count,
              s1_col=int(options.s1_col), s2_col=int(options.s2_col), s3_col=int(options.s3_col))
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#
