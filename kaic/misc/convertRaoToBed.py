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
        help='''Rao 2012 et al. supplementary file of genomic location pairs (+ additional information)'''
    );
    
    parser.add_argument(
        'output',
        help='''Output file (BEDPE format)'''
    );
    
    
    args = parser.parse_args()
    in_file = args.input
    out_file = args.output
    
    with open(in_file, 'r') as f:
        with open(out_file, 'w') as o:
            header = f.readline().rstrip().split("\t")
            header[0] = "chrom"
            header[1] = "start"
            header[2] = "end"
            o.write("\t".join(header) + "\n")
            
            line = f.readline()
            while line != '':
                o.write(line)

                line = f.readline()
    