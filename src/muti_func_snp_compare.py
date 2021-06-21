
#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""
@author:Yang Zhengzhao
@time:2021/04/19
"""

from __future__ import division
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
    if len(geno) == 3 and geno[1] == '/' and geno != "./.":
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


def IsHomo(geno): 
    if geno[0] == geno[2]:
        return 1
    else:
        return 0
#


def Isunconsistent(g1, g2):  
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


def GenoPhase(infile, bin_size, two_diff_level, length, DSR_count, mask_cnv, DP_high, DP_low,
              GQ_sample, LEVEL, outfile, chr_col=1, pos_col=2, s1_col=5, s2_col=6):
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

    #
    if DSR_count == 'on':  # counting differnt snp number between two samples with alt genotype
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
            #print("\t".join(["%s" % DP_low[0], "%s" % DP_high[0], "%s" % GQ_sample]))
            #
            while pos > end:
                DSR = count_diff_num*bin_size/(bin_size-count_miss_num)
                print("\t".join(["%s" % chr, "%s" % start, "%s" % end, "%.2f" % DSR]))
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
        if length == "default":
            chr_length = {'chr1A': '594102056', 'chr1B': '689851870', 'chr1D': '495453186',
                    'chr2A': '780798557', 'chr2B': '801256715', 'chr2D': '651852609',
                    'chr3A': '750843639', 'chr3B': '830829764', 'chr3D': '615552423', 
                    'chr4A': '744588157', 'chr4B': '673617499', 'chr4D': '509857067',
                    'chr5A': '709773743', 'chr5B': '713149757', 'chr5D': '566080677',
                    'chr6A': '618079260', 'chr6B': '720988478', 'chr6D': '473592718',
                    'chr7A': '736706236', 'chr7B': '750620385', 'chr7D': '638686055',
                    'chrUn': '480980714'}
            final_bin_size = int(chr_length[chr])
        else:
            final_bin_size = int(length)
        bin_final = final_bin_size - start
        DSR = count_diff_num*bin_final/(bin_final-count_miss_num)
        print("\t".join(["%s" % chr, "%s" % start, "%s" % final_bin_size, "%.2f" % DSR]))
    #
    if two_diff_level == 'on':
        count = 0
        level = LEVEL.strip().split(",")
        low_level = float(level[0])
        #high_level = float(level[1])
        thefile = open(infile)
        while True:
            buffer = thefile.read(1024*8192)
            if not buffer:
                break
            count += buffer.count('\n')
        thefile.close()
        print(level,low_level,count)
        if count >= 1:
            lis = []
            for line in IN:
                tokens = line.strip().split("\t")
                chr = tokens[0]
                start = tokens[1]
                end = tokens[2]
                diff = tokens[3]
                lis.append([chr, start, end, diff])
            frame = pd.DataFrame(lis)
            frame[3] = frame[3].astype(float)
            frame[3] = frame[3].astype(int)
            frame[1] = frame[1].astype(int)
            frame[2] = frame[2].astype(int)
            frame.loc[(low_level < frame[3]), 5] = 'PHR'
            frame.loc[(low_level >= frame[3]), 5] = 'SGR'
            #frame.loc[(high_level >= frame[3]) & (low_level < frame[3]), 5] = 'mid'
            frame.to_csv(outfile, sep='\t', na_rep='Na', float_format='%.1f', columns=[0, 1, 2, 5], header=False, index=False)
        else:
            diff_file = open(outfile, 'w')
            diff_file.close()

    if mask_cnv == 'on':
        for line in IN:
            tokens = line.strip().split("\t")
            CHR = tokens[0]
            start = tokens[1]
            end = tokens[2]
            RD = float(tokens[3])
            if RD <= 0.5:
                print(CHR, start, end, "ownCNV_deletion", sep="\t")
            if RD >= 1.5:
                print(CHR, start, end, "ownCNV_duplication", sep="\t")
    #
    if infile:
        IN.close()
    #
#



from optparse import OptionParser

# ===========================================
def main():
    usage = "Author : Yang, zhengzhao; yangzhengzhao@cau.edu.cn; 2021-04\n"
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
    parser.add_option("--length", dest="length",
                      help="column id for offSpring [default: %default]",
                      metavar="INT",
                      default="default")
    parser.add_option("--LEVEL", dest="LEVEL",
                      help="column id for offSpring [default: %default]",
                      metavar="INT",
                      default="10,1000")                                                        
    parser.add_option("-c", dest="chr_col",
                      help="column id for chromosome [default: %default]",
                      metavar="INT",
                      default=1)
    parser.add_option("-p", dest="pos_col",
                      help="column id for position [default: %default]",
                      metavar="INT",
                      default=2)
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
              two_diff_level=options.two_diff_level,
              chr_col=int(options.chr_col), pos_col=int(options.pos_col), length=options.length,
              s1_col=int(options.s1_col), s2_col=int(options.s2_col))
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#
