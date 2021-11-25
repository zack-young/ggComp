#!/usr/bin/env python
import os
import gzip
import sys
import re
import math

def combine(l, n=2):
    combine_sample = []
    one = [0] * n
    def next_c(li = 0, ni = 0):
        if ni == n:
            combine_sample.append(copy.copy(one))
            return
        for lj in range(li, len(l)):
            one[ni] = l[lj]
            next_c(lj + 1, ni + 1)
    next_c()
    return combine_sample


def Graphmaker(pathway, sample, num):
    num = int(num)
    pathway = pathway.strip().split(",")
    sample = sample.strip().split(",")
    for i in range(6, -1, -1):
        lett = ["A", "B", "D"]
        for j in range(2, -1, -1):
            if num == 1: 
                print(str(pathway[0]) + "/" + 'chr' + str(i+1) + lett[j] + str(sample[0]) + "\tchr" + str(i+1) + lett[j])
            else:
                print(str(pathway[0]) + "/" + 'chr' + str(i+1) + lett[j] + str(sample[0]) + "\tchr" + str(i+1) + lett[j])
                for item in range(1,num):
                    print(str(pathway[item]) + "/" + 'chr' + str(i+1) + lett[j] + str(sample[item]) + "\tNA")
            print("NA\tNA")

from optparse import OptionParser

# ===========================================
def main():
    usage = '''Usage: %prog [-i <input>] [-b 1000] [-o <output>] 
       [-p "jm22_snp_indel_analysis/snp,lx99_snp_indel_analysis/snp"] [-s "_snp_altLEVEL,_snp_altLEVEL"] [-n 2]'''
    #
    parser = OptionParser(usage)
    parser.add_option("-p", dest="pathway",
                      help="Input file, use STDIN if omit; "
                      "gzip file end with \".gz\".", metavar="STRING")
    parser.add_option("-o", dest="outfile",
                      help="Output file, use STDOUT if omit; "
                      "gzip file end with \".gz\".", metavar="FILE")
    parser.add_option("-s", dest="sample",
                      help="Output file, use STDOUT if omit; "
                      "gzip file end with \".gz\".", metavar="STRING")
    parser.add_option("-n", dest="num",
                      help="column id for offSpring [default: %default]", metavar="INT",
                      default=1)
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
    Graphmaker(pathway=options.pathway, sample=options.sample, num=options.num)
    #
#
if __name__ == "__main__":
    main()
                                       
