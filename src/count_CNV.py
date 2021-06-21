#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from __future__ import division  # very important !!!!

import gzip
import sys
from scipy import stats
from optparse import OptionParser


def Count(infile):
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

    lis = []
    for line in IN:
        item = line.strip("\n").split("\t")
        if item[3] != '0':
            lis.append(item[3])
    mode_num = float(stats.mode(lis)[0][0])
    print("%.2f" % mode_num)


def main():
    usage = "Usage: %prog [-i <input>] \n" \
            "Author : Yang, zhengzhao; yangzhengzhao@cau.edu.cn; 2021-04\n" \
            "Description: Counting the mode.\n" \
            "Input format:\n" \
            "   chr start end read_depth\n" \
            "Output format:\n" \
            "   mode\n"
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="infile",
                  help="Input file, use STDIN if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
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
    Count(infile=options.infile)
    #
#


if __name__ == "__main__":
    main()
#
