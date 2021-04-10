#!/data/user/yangzz/worktools/anaconda3/bin/python
###  #!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""
@author:Yang Zhengzhao
@file:finalgeno.py
@time:2018/7/2417:35
"""

from __future__ import division  # very important !!!!
import os
import gzip
import sys
import re
import math
import pandas as pd

def error(msg):
    print(sys.stderr, 'ERROR: %s' % msg)
    exit(1)
#


def IsGeno(geno):
    if len(geno)== 3 and geno[1] == '/' and geno != "./.":
        return 1
    else:
        return 0
#


def IsALT(geno):
    if geno[0] != '0' or geno[2] != "0":
        return 1
    else:
        return 0
#


def IsHomo(geno):  #判断纯合
    if geno[0] == geno[2]:
        return 1
    else:
        return 0
#


def Isunconsistent(g1, g2):   #判断父母本差异密度
    if g1 == g2:
        return 0
    else:
        return 1
    #

def Change_dot(gq):
    global lis3
    lis3 = []
    if gq != '.':
        lis3.append(int(gq))
    else:
        lis3.append(0)
    return(lis3)

def ChangeGQ(gt, gq):
    global lis2
    lis2 = []
    if IsGeno(gt):
        lis2.append(int(gq))
    else:
        lis2.append(0)
    return(lis2)
#


def GenoPhase(infile, bin_size, two_diff_level,
              sample, DSR_count, mask_cnv, sample1, sample2, chromosome,
              DP_high, DP_low, GQ_sample, LEVEL, outfile,
              chr_col=1, pos_col=2, ref_col=3, alt_col=4, s1_col=5, s2_col=6):
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
    length = {'chr1A': '594102056', 'chr1B': '689851870', 'chr1D': '495453186',
              'chr2A': '780798557', 'chr2B': '801256715', 'chr2D': '651852609',
              'chr3A': '750843639', 'chr3B': '830829764', 'chr3D': '615552423', 
              'chr4A': '744588157', 'chr4B': '673617499', 'chr4D': '509857067',
              'chr5A': '709773743', 'chr5B': '713149757', 'chr5D': '566080677',
              'chr6A': '618079260', 'chr6B': '720988478', 'chr6D': '473592718',
              'chr7A': '736706236', 'chr7B': '750620385', 'chr7D': '638686055',
              'chrUn': '480980714'}
    #
    if DSR_count == 'on':  # counting differnt snp number between two samples with alt genotype
        print("\t".join(['CHR', 'start', 'end', 'DSR']))
        firstline = IN.readline().strip()
        item = firstline.strip().split("\t")
        firstpos = int(item[pos_col - 1])
        count_miss_num = 0
        count_diff_num = 0
        if int(firstpos) <= int(bin_size):
            start = 1
            end = bin_size
        else:
            num = int(firstpos/bin_size)-1
            start = (num+1) * bin_size
            end = start + bin_size
        IN.seek(0, os.SEEK_SET)
        DP_high = DP_high.strip().split(",")
        DP_low = DP_low.strip().split(",")
        for line in IN:
            tokens = line.strip().split("\t")
            chr = tokens[chr_col - 1]
            pos = int(tokens[pos_col - 1])
            #
            GT1 = tokens[s1_col - 1]
            DP1 = Change_dot(tokens[s1_col])[0]
            GQ1 = tokens[s1_col + 1]
            lis3 = ChangeGQ(GT1, GQ1)
            GQ1 = lis3[0]
            #
            GT2 = tokens[s2_col - 1]
            DP2 = Change_dot(tokens[s2_col])[0]
            GQ2 = tokens[s2_col + 1]
            lis3_1 = ChangeGQ(GT2, GQ2)
            GQ2 = lis3_1[0]
            #
            while pos > end:
                DSR = count_diff_num*bin_size/(bin_size-count_miss_num)
                print("\t".join(["%s" % chr, "%s" % start, "%s" % end, "%s" % DSR]))
                # #
                start += bin_size
                end += bin_size
                count_diff_num = 0
                count_miss_num = 0
            if int(DP_low[0]) <= DP1 <= int(DP_high[0]) and int(DP_low[0]) <= DP2 <= int(DP_high[0]) and GQ1 >= int(GQ_sample) and GQ2 >= int(GQ_sample):
                if IsGeno(GT1) and IsGeno(GT2):
                    if IsHomo(GT1) and IsHomo(GT2):
                        if Isunconsistent(GT1, GT2):
                            count_diff_num += 1
                else:
                    count_miss_num += 1
            else:
                count_miss_num += 1

        #
        final_bin_size = int(length[chr])
        bin_final = final_bin_size - start
        DSR = count_diff_num*bin_final/(bin_final-count_miss_num)
        print("\t".join(["%s" % chr, "%s" % start, "%s" % final_bin_size, "%s" % DSR]))
    #
    if two_diff_level == 'on':
        count = 0
        level = LEVEL.strip().split(",")
        low_level = float(level[0])
        high_level = float(level[1])
        thefile = open(infile)
        while True:
            buffer = thefile.read(1024*8192)
            if not buffer:
                break
            count += buffer.count('\n')
        thefile.close()

        chr = chromosome
        if count > 1:
            firstline = IN.readline().strip()
            lis = []
            for line in IN:
                tokens = line.strip().split("\t")
                start = tokens[1]
                end = tokens[2]
                diff = tokens[3]
                lis.append([chr, start, end, diff])
            frame = pd.DataFrame(lis)
            frame[3] = frame[3].astype(float)
            frame[3] = frame[3].astype(int)
            frame[1] = frame[1].astype(int)
            frame[2] = frame[2].astype(int)
            frame.loc[(high_level < frame[3]), 5] = 'high'
            frame.loc[(low_level >= frame[3]), 5] = 'low'
            frame.loc[(high_level >= frame[3]) & (low_level < frame[3]), 5] = 'mid'
            frame.to_csv(outfile, sep='\t', na_rep='Na', float_format='%.1f', columns=[0, 1, 2, 5], header=False, index=False)
        else:
            diff_file = open(outfile,'w')
            diff_file.close()
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
                print(filename3, start, end, "ownCNV_deletion", sample, sep="\t")
            if RD >= 1.5:
                print(filename3, start, end, "ownCNV_duplication", sample, sep="\t")
            start += bin_size
            end += bin_size
    #
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
                      default=10000)
    parser.add_option("--DSR_count", dest="DSR_count",
                      help="column id for offSpring [default: %default]",
                      metavar="STRING",
                      default='off')                   
    parser.add_option("--two_diff_level", dest="two_diff_level",
                      help="column id for offSpring [default: %default]",
                      metavar="STRING",
                      default='off')
    parser.add_option("--mask_cnv", dest="mask_cnv",
                      help="column id for offSpring [default: %default]",
                      metavar="STRING",
                      default='off')
    parser.add_option("--DP_high", dest="DP_high",
                      help="column id for offSpring [default: %default]",
                      metavar="INT",
                      default=1)
    parser.add_option("--DP_low", dest="DP_low",
                      help="column id for offSpring [default: %default]",
                      metavar="INT",
                      default=1)
    parser.add_option("--GQ_sample", dest="GQ_sample",
                      help="column id for offSpring [default: %default]",
                      metavar="INT",
                      default=1)
    parser.add_option("--sample1", dest="sample1",
                      help="column id for offSpring [default: %default]",
                      metavar="STRING",
                      default='a')
    parser.add_option("--sample2", dest="sample2",
                      help="column id for offSpring [default: %default]",
                      metavar="STRING",
                      default='b')
    parser.add_option("--chromosome", dest="chromosome",
                      help="column id for offSpring [default: %default]",
                      metavar="STRING",
                      default='chr1A')
    parser.add_option("--sample", dest="sample",
                      help="column id for offSpring [default: %default]",
                      metavar="INT",
                      default=1)
    parser.add_option("--LEVEL", dest="LEVEL",
                      help="column id for offSpring [default: %default]",
                      metavar="INT",
                      default="10, 1000")                                                        
    parser.add_option("-c", dest="chr_col",
                      help="column id for chromosome [default: %default]",
                      metavar="INT",
                      default=1)
    parser.add_option("-p", dest="pos_col",
                      help="column id for position [default: %default]",
                      metavar="INT",
                      default=2)
    parser.add_option("-r", dest="ref_col",
                      help="column id for position [default: %default]",
                      metavar="INT",
                      default=3)
    parser.add_option("-a", dest="alt_col",
                      help="column id for position [default: %default]",
                      metavar="INT",
                      default=4)
    parser.add_option("-1", dest="s1_col",
                      help="column id for 1st parent [default: %default]",
                      metavar="INT",
                      default=5)
    parser.add_option("-2", dest="s2_col",
                      help="column id for 2nd parent [default: %default]",
                      metavar="INT",
                      default=8)
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
    GenoPhase(infile=options.infile, bin_size=int(options.bin_size), outfile=options.outfile,
              mask_cnv=options.mask_cnv, DP_high=options.DP_high, DSR_count=options.DSR_count,
              GQ_sample=options.GQ_sample, DP_low=options.DP_low, LEVEL=options.LEVEL,
              sample=options.sample, two_diff_level=options.two_diff_level,
              sample1=options.sample1, sample2=options.sample2, chromosome=options.chromosome,
              chr_col=int(options.chr_col), pos_col=int(options.pos_col),
              ref_col=int(options.ref_col), alt_col=int(options.alt_col),
              s1_col=int(options.s1_col), s2_col=int(options.s2_col))
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#
