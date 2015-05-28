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
    n_entries = []
    n_sum = 0
    for basename in args.input:
        
        file_name1 = basename + "_1.fastq"
        file_name2 = basename + "_2.fastq"
        
        if not os.path.isfile(file_name1):
            raise IOError("File " + file_name1 + " not found")
        if not os.path.isfile(file_name2):
            raise IOError("File " + file_name2 + " not found")
        
        pairs.append([file_name1, file_name2])
        n = int(get_number_of_lines(file_name1)/4)
        n_sum += n
        n_entries.append(n)
        logging.info("\t%s\t%d" % (basename, n))
    
    n_entries_ratios = []
    for n in n_entries:
        n_entries_ratios.append(n/n_sum)
    
    print n_entries_ratios
    
    # step 2:
    # extract samples from FASTQ files
    for sample_size in args.sample_sizes:
        
        # make sure we are not sampling more than we have
        sample_size = min(n_sum, sample_size)
        
        # create new folders
        
        sample_folder = "%s/sample_%d/" % (args.output_folder,sample_size)
        fastq_folder = sample_folder + "fastq"
        make_dir(fastq_folder)
        sam_folder = sample_folder + "sam"
        make_dir(sam_folder)
        sam_filtered_folder = sample_folder + "sam/filtered"
        make_dir(sam_filtered_folder)
        hic_folder = sample_folder + "hic"
        make_dir(hic_folder)
        hic_filtered_folder = sample_folder + "hic/filtered"
        make_dir(hic_filtered_folder)
        binned_folder = sample_folder + "binned"
        make_dir(binned_folder)
        binned_filtered_folder = sample_folder + "binned/filtered"
        make_dir(binned_filtered_folder)
        corrected_folder = sample_folder + "corrected"
        make_dir(corrected_folder)
        
    
        file_sample_sizes = []
        fss_sum = 0
        for ratio in n_entries_ratios:
            fss = int(ratio*sample_size)
            fss_sum += fss
            file_sample_sizes.append(fss)
        
        # correct for integer conversion of sample sizes
        if fss_sum != sample_size:
            file_sample_sizes[0] += sample_size-fss_sum
        
        print file_sample_sizes
        print sum(file_sample_sizes)
        
        for i in range(0, len(pairs)):
            file1 = pairs[i][0]
            file2 = pairs[i][1]
            fss = file_sample_sizes[i]
            out1 = "%s/%s.%d_1.fastq" % (fastq_folder,args.input[i],fss)
            out2 = "%s/%s.%d_2.fastq" % (fastq_folder,args.input[i],fss)
            
            
            sample_fastq(file1,file2,out1,out2,sample_size=fss)
            
        
        
    
    