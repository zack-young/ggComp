#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""
@author:Young zhengzhao
@file:finalgeno.py
@time:2021/4/1
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

def Makelis(chromosome):
    global list1
    list1 = []
    #f = open("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/MASK/chr" + chromosome + "_combine_total", "r")
    f = open(chromosome, "r")
    mask = f.readlines()
    for line in mask:
        item = line.strip("\n").split("\t")
        list1.append([item[1],item[2], item[3]])
    f.close
    return(list1)

def Makelis_meta_data(chromosome):
    global list_meta1
    list_meta1 = []
    #f = open("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/MASK/chr" + chromosome + "_combine_total", "r")
    f = open(chromosome, "r")
    mask = f.readlines()
    for line in mask:
        item = line.strip("\n").split("\t")
        list_meta1.append([item[0],item[1]])
    f.close
    return(list_meta1)


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

def GenoPhase(infile, bin_size, density_graph_raw_vcf_two, 
              snp_graphy, single_DP_GQ, sample, meta_data, time_data, total,
              high_threshold, DP_high, DP_low, GQ_sample, LEVEL, centrosome,
              chromosome, transition_transversion_count, indel_count, suffix,
              chr_col=1, pos_col=2, ref_col=3, alt_col=4):
    # open the infile
    #print(sys.getsizeof(ob)) show memory
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
    #  muti_func_snp_filter.py single_DP_GQ  
    # snp_graphy
    # density_graph_raw_vcf_two
    # transition_transversion_count
    # indel_count
    with open(meta_data, 'r') as IN:     
        meta_dic = {'CHR':'# [1]CHROM','POS':'[2]POS','REF':'[3]REF','ALT':'[4]ALT'}
        list_use = ['# [1]CHROM','[2]POS','[3]REF','[4]ALT']
        list_matrix = [[],[]]
        for line in IN: 
            tokens = line.strip().split("\t") 
            list_matrix[0].append(tokens[0])
            list_matrix[1].append(tokens[1])
            list_use.append(tokens[1]) 
            meta_dic[tokens[0]]=tokens[1]
    #
    with open(infile, 'r') as IN:
        list_vcf = []
        firstline = IN.readline().strip() 
        header = firstline.strip().split("\t")
        IN.seek(0, os.SEEK_SET)
        list_index=[]
        for item in list_use:
            [list_index.append(i) for i,x in enumerate(header) if x.find(item)!=-1]
        for line in IN: 
            tokens = line.strip().split("\t")
            item_list=[tokens[i] for i in list_index]
            list_vcf.append(item_list)
    header = list_vcf[0]
    del list_vcf[0]
    #
    meta_list = [[],[]]
    for i in range(len(header)):
        meta_list[0].append([x[0] for x in meta_dic.items() if x[1] in header[i]][0])
        meta_list[1].append([header[i] for x in meta_dic.items() if x[1] in header[i]][0])
    tuples = list(zip(*meta_list))
    index = pd.MultiIndex.from_tuples(tuples, names=['sample_index', 'sample_vcf'])
    frame = pd.DataFrame(list_vcf, columns=index)
    #
    tuples = list(zip(*list_matrix))
    index = pd.MultiIndex.from_tuples(tuples, names=['sample_index', 'sample_vcf'])
    s = pd.DataFrame(np.ones((len(frame),len(list_matrix[0])), dtype = np.int),columns=index) 
    s=s.sort_index(axis=1)
    #
    tmp_1=pd.concat([frame[['POS']],s],axis=1)  ## the  main define matrix
    tmp_1[['POS']] =  tmp_1[['POS']].astype(int)
    for item in list_matrix[0]:
        list1 = Makelis("tmp/" + total +'/' + item + "_" + time_data + "/" + chromosome + "." + "mask_CNV") 
        tmp_1.loc[(frame[item].iloc[:,2].isin(['.']))|(frame[item].iloc[:,1].isin(['.']))|(frame[item].iloc[:,0].isin(['./.'])) ,item] = '0'
        GT=frame[item].iloc[:,0].str.split('/',expand=True)
        tmp_1.loc[(GT[0]!=GT[1]) ,item] = '0'    
        for num in range(len(list1)):
            tmp_1.loc[(int(list1[num][0]) <= tmp_1['POS']['[2]POS'])&(tmp_1['POS']['[2]POS'] <= int(list1[num][1])) ,item] = '0'

    # tmp_1 is a CNV matrix
    #sample_list = sample_2.strip().split(",")
    #s1_col = [i for i,x in enumerate(header) if x.find(sample_list[0]+':')!=-1][0]        
    #tmp_ma.iat[1,1] #choose item in matrix based on index


    #
    if snp_graphy == 'on': #counting differnt snp number between two samples with alt genotype
        # for num in range(len(sample_list)):
        #     s2_col = [i for i,x in enumerate(header) if x.find(sample_list[1]+':')!=-1][num]
        #     diff_file = open("tmp/"+ meta_list[num_1][0] + "_" + time_data + "/" + chromosome + suffix, 'w+')
        DP_high = DP_high.strip().split(",")
        DP_low = DP_low.strip().split(",")
        sample = sample.strip().split(",")

        b = chromosome
        for item in list_matrix[0]:
            frame.loc[(frame[item].iloc[:,1].isin(['.'])) ,(item,frame[item].columns[1])] = 0
            frame.loc[(frame[item].iloc[:,2].isin(['.'])) ,(item,frame[item].columns[2])] = 0
            tmp_1.loc[(frame[item].iloc[:,1].astype(int) >=int(DP_high[0]))|(frame[item].iloc[:,1].astype(int) <= int(DP_low[0]))|(frame[item].iloc[:,2].astype(int) < int(GQ_sample)) ,item] = '0'

        sample_list = list_matrix[0]

        for sample_item in sample:
            sample_list.remove(sample_item)
            if len(sample_list) > 0 :
                for item in sample_list:
                    diff_file = open("tmp/" + sample_item + '_' +item + "_" + time_data + "/" + chromosome + ".homo_snp_density" + suffix, 'w+')
       
                    use = pd.concat([frame[['POS']], frame[[sample_item]], tmp_1[[sample_item]], frame[[item]], tmp_1[[item]]], axis=1)
                    use_filter= use[(~ (use[sample_item].iloc[:,3].isin(['0'])))&(~(use[item].iloc[:,3].isin(['0']))) ]
                    if use_filter.empty:
                        diff_file.close()
                        continue
                    else:
                        count_same=0
                        count_diff=0
                        count_same_num=0
                        count_diff_num=0 
                        firstpos = int(use_filter['POS'].iloc[0,0])
                        if int(firstpos) <= int(bin_size):
                            start = 1
                            end = bin_size+1
                        else:
                            num = int(firstpos/bin_size)-1
                            start = (num+1) * bin_size + 1
                            end = start + bin_size
                        for line in range(len(use_filter)):
                            tokens = use_filter.iloc[line]
                            pos = int(tokens[0])
                            GT1 = tokens[1]
                            GT2 = tokens[5]
                            if IsALT(GT1) or IsALT(GT2):
                                if Isunconsistent(GT1, GT2):
                                    count_diff_num += 1
                                else:
                                    count_same_num += 1
                            while pos > end:
                                diff_file.write("\t".join([b, "%s" % start, "%s" % end, "%s" % count_diff_num, "%s" % count_same_num+'\n']))
                                # #
                                start += bin_size
                                end += bin_size
                                count_diff_num = 0
                                count_same_num = 0
                        diff_file.write("\t".join([b, "%s" % start, "%s" % end, "%s" % count_diff_num, "%s" % count_same_num+'\n']))             
                        diff_file.close()
    #
    if density_graph_raw_vcf_two == 'on': # level differnt snp number between two samples with alt genotype
        DP_high = DP_high.strip().split(",")
        DP_low = DP_low.strip().split(",")
        level = LEVEL.strip().split(",")
        sample = sample.strip().split(",")
        b = chromosome
        # snp density level
        frame.replace('.', '0' ,inplace=True) 
        for item in list_matrix[0]:
            tmp_1.loc[(frame[item].iloc[:,1].astype(int) >=int(DP_high[0]))| \
                    (frame[item].iloc[:,1].astype(int) <= int(DP_low[0]))| \
                    (frame[item].iloc[:,2].astype(int) < int(GQ_sample)) ,item] = '0'
        sample_list = list_matrix[0]

        for sample_item in sample:
            sample_list.remove(sample_item) # list_matrix will also remove item
            if len(sample_list) > 0 :
                for item in sample_list:
                    diff_file = open("tmp/" + total + '/'+sample_item + '_' +item + "_" + time_data + "/" + chromosome + ".homo_snp_density"+ suffix, 'w+')
                    for num in range(int(int(tmp_1.loc[0,'POS'])/1000000), int(int(tmp_1.loc[len(tmp_1)-1,'POS'])/1000000)+1):
                        dt_tmp_1 = tmp_1.loc[(frame['POS'].iloc[:,0].astype(int)<1000000*(num+1))&(frame['POS'].iloc[:,0].astype(int)>1000000*num),]
                        dt_frame = frame.loc[(frame['POS'].iloc[:,0].astype(int)<1000000*(num+1))&(frame['POS'].iloc[:,0].astype(int)>1000000*num),]
                        use_filter= dt_frame[(~ (dt_tmp_1[sample_item].iloc[:,0].isin(['0'])))&(~(dt_tmp_1[item].iloc[:,0].isin(['0']))) ]
                        if use_filter.empty:
                            continue
                        else:
                            count_diff_num = len(use_filter[use_filter[sample_item].iloc[:,0] != use_filter[item].iloc[:,0]])
                            count_same_num = len(use_filter[use_filter[sample_item].iloc[:,0] == use_filter[item].iloc[:,0]])
                            start = 1000000*num+1
                            end = 1000000*(num+1)
                            diff_file.write("\t".join([b, "%s" % start, "%s" % end, "%s" % count_diff_num, "%s" % count_same_num + '\n']))
                    diff_file.close()

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
    parser.add_option("--density_graph_raw_vcf_two", dest="density_graph_raw_vcf_two",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--snp_graphy", dest="snp_graphy",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--single_DP_GQ", dest="single_DP_GQ",
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
    parser.add_option("--total", dest="total",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = 'total')
    parser.add_option("--DP_high", dest="DP_high",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = '1')
    parser.add_option("--DP_low", dest="DP_low",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = '1')
    parser.add_option("--time_data", dest="time_data",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = 1)
    parser.add_option("--GQ_sample", dest="GQ_sample",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = 1)
    parser.add_option("--LEVEL", dest="LEVEL",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = 1)
    parser.add_option("--sample", dest="sample",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = "cs,c46")                      
    parser.add_option("--centrosome", dest="centrosome",
                  help="column id for offSpring [default: %default]", metavar="INT",
                      default = "215,239,173,337,346,268,300,346,241,300,317,184,254,206,189,286,327,211,357,287,340")                
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
    parser.add_option("--chromosome", dest="chromosome",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = '1B')
    parser.add_option("--suffix", dest="suffix",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = '1B')
    parser.add_option("--meta_data", dest="meta_data",
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
    GenoPhase(infile=options.infile, bin_size=int(options.bin_size),
              snp_graphy=options.snp_graphy, single_DP_GQ=options.single_DP_GQ,
              indel_count=options.indel_count, suffix=options.suffix, total=options.total,
              high_threshold=int(options.high_threshold), DP_high=options.DP_high, GQ_sample=options.GQ_sample,
              DP_low=str(options.DP_low), LEVEL=options.LEVEL, centrosome=options.centrosome, sample=options.sample, 
              chromosome=options.chromosome, meta_data=options.meta_data, density_graph_raw_vcf_two=options.density_graph_raw_vcf_two, 
              chr_col=int(options.chr_col), pos_col=int(options.pos_col), time_data=options.time_data,
              ref_col=int(options.ref_col), alt_col=int(options.alt_col), transition_transversion_count=options.transition_transversion_count)
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#
