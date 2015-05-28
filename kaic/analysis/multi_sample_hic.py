#!/usr/bin/env python

import argparse
import gridmap
import logging
import os.path
from kaic.tools.files import get_number_of_lines
logging.basicConfig(level=logging.DEBUG)


def intList(thisList):
    return [int(x) for x in thisList.split(",")]
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser();
    
    parser.add_argument(
        'input',
        nargs='+',
        help='''FASTQ input files basenames (without _*.fastq)'''
    );
    
    parser.add_argument(
        '-g', '--genome', dest='genome',
        help='''Genome file''',
        required=True
    );
    
    parser.add_argument(
        '-s', '--sample-sizes', dest='sample_sizes',
        type=intList,
        help='''Comma-separated list of sample sizes''',
        required=True
    );
    
    parser.add_argument(
        '-o', '--output-folder', dest='output_folder',
        help='''Output folder''',
        required=True
    );
    
    args = parser.parse_args()
    
    
    
    # step one:
    # reconstruct paired FASTQ file names
    pairs = []
    n_lines = []
    n_sum = 0
    for basename in args.input:
        file_name1 = basename + "_1.fastq"
        file_name2 = basename + "_2.fastq"
        
        if not os.path.isfile(file_name1):
            raise IOError("File " + file_name1 + " not found")
        if not os.path.isfile(file_name2):
            raise IOError("File " + file_name2 + " not found")
        
        pairs.append([file_name1, file_name2])
        n = get_number_of_lines(file_name1)
        n_sum += n
        n_lines.append(n)
    
    n_lines_ratios = []
    for i in range(0,len(n_lines)):
        n_lines_ratios[i] = n_lines[i]/n_sum
        
    print n_lines_ratios