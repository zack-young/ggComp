#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from __future__ import division  # very important !!!!
from collections import defaultdict
import os
import gzip
import sys
import re
import math
from scipy import stats

def Count(infile, combine_data, count, compensent_deletion, normalize, split, chrom, pop_num, compensent_duplication):
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

    dic1={}
    if count == 'on' :
        for line in IN:
            line = line.strip().split('\t')
            dic1[line[1]] = dic1.setdefault(line[1], 0)+1
        for key,value in dic1.items():
            print('{key}\t{key2}\t{value}'.format(key=key, key2=int(key)+1000000, value=value/int(pop_num)), sep='\t')
#
    if normalize == 'on':
        lis=[]
        for line in IN:
            item = line.strip("\n").split("\t")
            if item[0] != '0':
                lis.append(item[0])
        mode_num=float(stats.mode(lis)[0][0])
        print("%.2f" % mode_num)

    if split == 'on':
        for line in IN:
            item = line.strip("\n").split("\t")
            for i in range(int(item[1]),int(item[2]),1000000):
                print(item[0],i,i+1000000,item[3],sep="\t")
    if compensent_deletion == 'on' :
        dic1 = defaultdict(list)
        for line in IN:
            item = line.strip("\n").split("\t")
            dic1[item[0]].append(item[1])
            dic1[item[0]].append(item[2])
        wholen = (594102056,
                    689851870,
                    495453186,
                    780798557,
                    801256715,
                    651852609,
                    750843639,
                    830829764,
                    615552423,
                    744588157,
                    673617499,
                    509857067,
                    709773743,
                    713149757,
                    566080677,
                    618079260,
                    720988478,
                    473592718,
                    736706236,
                    750620385,
                    638686055)

        b = chrom
        dic2 = {"A":"0","B":"1","D":"2"}
        num = (int(b[0])-1)*3
        ch = int(dic2[b[1]])
        loc = num + ch
        if '1' in dic1:
            use_num = str(1.0-float(dic1["1"][1]))
            print("\t".join(['1', '1000001', "%s" % use_num]))
        else:
            print("\t".join(['1', '1000001', "1.0"]))
        for item in range(1, int(wholen[loc]/1000000)+1):            
            item = str(item*1000000+1)
            if item in dic1:
                use_num = 1.0-float(dic1[item][1])
                print("\t".join(["%s" % item, "%s" % dic1[item][0], "%s" % use_num]))
            else:
                print("\t".join(["%s" % item, "%s" % str(int(item)+1000000), "1.0"]))


    if compensent_duplication == 'on' :
        dic1 = defaultdict(list)
        for line in IN:
            item = line.strip("\n").split("\t")
            dic1[item[0]].append(item[1])
            dic1[item[0]].append(item[2])
        wholen = (594102056,
                    689851870,
                    495453186,
                    780798557,
                    801256715,
                    651852609,
                    750843639,
                    830829764,
                    615552423,
                    744588157,
                    673617499,
                    509857067,
                    709773743,
                    713149757,
                    566080677,
                    618079260,
                    720988478,
                    473592718,
                    736706236,
                    750620385,
                    638686055)

        b = chrom
        dic2 = {"A":"0","B":"1","D":"2"}
        num = (int(b[0])-1)*3
        ch = int(dic2[b[1]])
        loc = num + ch
        if '1' in dic1:
            use_num = str(1.0+float(dic1["1"][1]))
            print("\t".join(['1', '1000001', "%s" % use_num]))
        else:
            print("\t".join(['1', '1000001', "1.0"]))
        for item in range(1, int(wholen[loc]/1000000)+2):            
            item = str(item*1000000+1)
            if item in dic1:
                use_num = 1.0+float(dic1[item][1])
                print("\t".join(["%s" % item, "%s" % dic1[item][0], "%s" % use_num]))
            else:
                print("\t".join(["%s" % item, "%s" % str(int(item)+1000000), "1.0"]))
    if infile:
        IN.close()
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
    parser.add_option("-o", dest="outfile",
                  help="Output file, use STDOUT if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
    parser.add_option("--combine_data", dest="combine_data",
                  help="Input file, use STDIN if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
    parser.add_option("--count", dest="count",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--split", dest="split",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--compensent_deletion", dest="compensent_deletion",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--chrom", dest="chrom",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = '1A')
    parser.add_option("--pop_num", dest="pop_num",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = '81')
    parser.add_option("--normalize", dest="normalize",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
    parser.add_option("--compensent_duplication", dest="compensent_duplication",
                  help="column id for offSpring [default: %default]", metavar="STRING",
                      default = 'off')
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
    Count(infile=options.infile, combine_data=options.combine_data, count=options.count, chrom=options.chrom,
          normalize=options.normalize, compensent_deletion=options.compensent_deletion, split=options.split,
          pop_num=options.pop_num, compensent_duplication=options.compensent_duplication)
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#
