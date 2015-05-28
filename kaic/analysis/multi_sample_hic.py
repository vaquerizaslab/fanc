#!/usr/bin/env python

from __future__ import division

import argparse
import gridmap
import logging
import os.path
import os
from kaic.tools.files import get_number_of_lines
from kaic.hrandom.sample_fastq import sample_fastq
logging.basicConfig(level=logging.DEBUG)


def make_dir(dir_name):
    try: 
        os.makedirs(dir_name)
    except OSError:
        if not os.path.isdir(sample_folder):
            raise

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
    
    
    
    
    # step 1:
    # reconstruct paired FASTQ file names
    logging.info("Getting FASTQ file sizes...")
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
        logging.info("\t%s\t%d" % (basename, n))
    
    n_lines_ratios = []
    for n in n_lines:
        n_lines_ratios.append(n/n_sum)
    
    print n_lines_ratios
    
    # step 2:
    # extract samples from FASTQ files
    for sample_size in args.sample_sizes:
        
        # make sure we are not sampling more than we have
        sample_size = min(n_sum, sample_size)
        
        # create new folders
        sample_folder = "%s/sample_%d/" % (args.output_folder,sample_size) 
        make_dir(sample_folder + "fastq")
        make_dir(sample_folder + "sam")
        make_dir(sample_folder + "sam/filtered")
        make_dir(sample_folder + "hic")
        make_dir(sample_folder + "hic/filtered")
        make_dir(sample_folder + "binned")
        make_dir(sample_folder + "binned/filtered")
        make_dir(sample_folder + "corrected")
        
    
        file_sample_sizes = []
        fss_sum = 0
        for ratio in n_lines_ratios:
            fss = ratio*sample_size
            fss_sum += fss
            file_sample_sizes.append(fss)
        
        if fss_sum != sample_size:
            file_sample_sizes[0] += n_sum-fss_sum
        
        print file_sample_sizes
        print n_sum-fss_sum
        
        
    
    