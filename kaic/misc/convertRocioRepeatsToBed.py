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
        help='''Rocio's list of repeat element locations (+ additional information)'''
    );
    
    parser.add_argument(
        'output',
        help='''Output file (BED format)'''
    );
    
    
    args = parser.parse_args()
    in_file = args.input
    out_file = args.output
    
    with open(in_file, 'r') as f:
        with open(out_file, 'w') as o:
            header = "chrom\tstart\tend\tname\tscore/tstrand"
           
            o.write(header + "\n")
            
            i = 0
            for line in f:
                line = line.rstrip()
                chrom, source, t, start, end, score, strand, x, name = line.split("\t")
                
                chrom = "chr" + chrom
                strand = "+1" if strand == '+' else "-1"
                name = name + "_" + str(i)
                i += 1
                
                o.write(chrom + "\t" + start + "\t" + end + "\t" + name + "\t" + score + "\t" + strand + "\n")
    