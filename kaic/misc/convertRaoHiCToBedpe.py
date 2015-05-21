#!/usr/bin/env python

import argparse;

'''
Created on May 20, 2015

@author: kkruse1
'''



if __name__ == '__main__':
    parser = argparse.ArgumentParser();
    
    parser.add_argument(
        'input',
        help='''Rao 2012 et al. supplementary file of Hi-C coordinates (+ additional information)'''
    );
    
    parser.add_argument(
        'chromosome',
        help='''Chromosome name of the file (e.g. chr1)'''
    );
    
    parser.add_argument(
        'resolution',
        help='''Resolution of the file (e.g. 100000)'''
    );
    
    parser.add_argument(
        'output',
        help='''Output file (BEDPE format)'''
    );
    
    
    args = parser.parse_args()
    in_file = args.input
    out_file = args.output
    
    resolution = int(args.resolution)
    chrom = args.chromosome
    
    with open(in_file, 'r') as f:
        with open(out_file, 'w') as o:
            header = []
            header.append("chrom1")
            header.append("start1")
            header.append("end1")
            header.append("chrom2")
            header.append("start2")
            header.append("end2")
            header.append("score")
            o.write("\t".join(header) + "\n")
            
            for line in f:
                start1, start2, score = line.rstrip().split("\t")
                start1 = int(start1)
                start2 = int(start2)
                score = float(score)
                end1 = start1 + resolution
                end2 = start2 + resolution
                
                o.write("%s\t%d\t%d\t%s\t%d\t%d\t%.6E\n" % (chrom, start1, end1, chrom, start2, end2, score))
                