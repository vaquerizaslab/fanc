#!/usr/bin/env python

import argparse
import gridmap
import logging
logging.basicConfig(level=logging.DEBUG)


def intList(thisList):
    return [int(x) for x in thisList.split(",")]
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser();
    
    parser.add_argument(
        'input',
        nargs='+',
        help='''FASTQ input files'''
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
    
    print args