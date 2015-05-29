#!/usr/bin/env python

from __future__ import division

import argparse
from gridmap import Job, process_jobs
import logging
import os.path
#import os
from kaic.correcting.filterUnwantedLigations import removeUnwantedLigations
from kaic.genome.genomeTools import loadGenomeObject
from kaic.tools.files import get_number_of_lines
from kaic.hrandom.sample_fastq import sample_fastq
from kaic.mapping.iterativeMapping import iterative_mapping
from kaic.mapping.samToHiCDataSet import sam_to_hic
from kaic.merging.mergeHiCDataSets import merge_hic
from kaic.binning.hicHeatMap import bin_hic
#from kaic.correcting.iterativeCorrection import ice

logging.basicConfig(level=logging.INFO)


def make_dir(dir_name):
    try: 
        os.makedirs(dir_name)
    except OSError:
        if not os.path.isdir(dir_name):
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
        '-s', '--sample-size', dest='sample_size',
        type=int,
        help='''FASTQ sample size (number of reads)''',
        required=True
    );
    
    parser.add_argument(
        '-o', '--output-folder', dest='output_folder',
        help='''Output folder''',
        required=True
    );
    
    parser.add_argument(
        '-mm', '--mapping-min', dest='mapping_min',
        help='''Minimum read length for iterative mapping''',
        type=int,
        default=28
    );
    
    parser.add_argument(
        '-ms', '--mapping-step', dest='mapping_step',
        help='''Step size for iterative mapping''',
        type=int,
        default=2
    );
    
    parser.add_argument(
        '-mi', '--mapping-index', dest='mapping_index',
        help='''Bowtie index for iterative mapping''',
        required=True
    );
    
    parser.add_argument(
        '-r', '--resolutions', dest='resolutions',
        type=intList,
        default=[10000,20000,50000,100000,200000,500000,1000000],
        help='''List of binning resolutions (comma-separated)''',
    );
    
    args = parser.parse_args()
    
    genome = loadGenomeObject(args.genome)
    
    
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
    # create folder structure

    # make sure we are not sampling more than we have
    sample_size = min(n_sum, args.sample_size)
    
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
    
    
    # step 3:
    # calculate file sample sizes
    file_sample_sizes = []
    fss_sum = 0
    for ratio in n_entries_ratios:
        fss = int(ratio*sample_size)
        fss_sum += fss
        file_sample_sizes.append(fss)
    
    # correct for integer conversion of sample sizes
    if fss_sum != sample_size:
        file_sample_sizes[0] += sample_size-fss_sum
    
    
    # step 4:
    # write new FASTQ files
    sample_pairs = []
    sample_jobs = []
    for i in range(0, len(pairs)):
        file1 = pairs[i][0]
        file2 = pairs[i][1]
        fss = file_sample_sizes[i]
        out1 = "%s/%s.%d_1.fastq" % (fastq_folder,args.input[i],fss)
        out2 = "%s/%s.%d_2.fastq" % (fastq_folder,args.input[i],fss)
        
        # a. sample new files
        largs = [file1,file2,out1,out2]
        kwargs = {'sample_size': fss}
        job = Job(sample_fastq,largs,kwlist=kwargs,queue='all.q')
        sample_jobs.append(job)
        sample_pairs.append([out1,out2])
    
    # do the actual sampling
    process_jobs(sample_jobs,max_processes=4)
    
    
    
    # step 5:
    # map new files
    sam_pairs = []
    sam_jobs = []
    for i in range(0, len(sample_pairs)):
        file1 = sample_pairs[i][0]
        file2 = sample_pairs[i][1]
        out1 = "%s/%s.%d_1.sam" % (sam_folder, args.input[i], file_sample_sizes[i])
        out2 = "%s/%s.%d_2.sam" % (sam_folder, args.input[i], file_sample_sizes[i])
        
        job1 = Job(iterative_mapping,[[file1],[out1],args.mapping_index],kwlist={ 'min_length':args.mapping_min, 'step_size': args.mapping_step },queue='all.q',num_slots=8)
        job2 = Job(iterative_mapping,[[file2],[out2],args.mapping_index],kwlist={ 'min_length':args.mapping_min, 'step_size': args.mapping_step },queue='all.q',num_slots=8)
        sam_jobs.append(job1)
        sam_jobs.append(job2)
        sam_pairs.append([out1,out2])
        
    # do the actual mapping
    process_jobs(sam_jobs,max_processes=4)
    
    # step 6:
    # filter SAM files
    filtered_sam_pairs = []
    filtered_sam_jobs = []
    for i in range(0, len(sam_pairs)):
        file1 = sam_pairs[i][0]
        file2 = sam_pairs[i][1]
        out1 = "%s/%s.%d_1.filtered.sam" % (sam_filtered_folder, args.input[i], file_sample_sizes[i])
        out2 = "%s/%s.%d_2.filtered.sam" % (sam_filtered_folder, args.input[i], file_sample_sizes[i])
        
        kwargs = {
            'outputSam1': out1,
            'outputSam2': out2,
            'inwardCutoff': 25000,
            'outwardCutoff': 10000,
            'removeSingle': True,
            'removeSelf': True,
            'removeDuplicates': True
        }
        job = Job(removeUnwantedLigations,[file1,file2,args.genome],kwlist=kwargs,queue='all.q')
        filtered_sam_jobs.append(job)
        filtered_sam_pairs.append([out1,out2])
    
    # do the actual filtering
    process_jobs(filtered_sam_jobs,max_processes=4)
    
    # step 7.a:
    # create Hi-C objects
    hic_files = []
    hic_jobs = []
    for i in range(0, len(filtered_sam_pairs)):
        file1 = filtered_sam_pairs[i][0]
        file2 = filtered_sam_pairs[i][1]
        out = "%s/%s.hic" % (hic_folder, args.input[i])
        
        job = Job(sam_to_hic, [file1, file2, args.genome, out],queue='all.q')
        hic_jobs.append(job)
        hic_files.append(out)
    
    # do the actual hic creation
    process_jobs(hic_jobs,max_processes=4)
    
    # step 7.b:
    # merge Hi-C objects
    hic_file = "%s/all.hic" % hic_folder
    merge_hic(hic_files, genome, hic_file)
    
    # step 8:
    # bin Hi-C object
    binned_files = []
    bin_jobs = []
    for resolution in args.resolutions:
        binned_file = "%s/all.%d.hm" % (binned_folder, resolution)
        
        
        job = Job(bin_hic, [hic_file, args.genome, resolution, binned_file],queue='all.q')
        bin_jobs.append(job)
        binned_files.append(binned_file)
    
    # do the actual binning
    process_jobs(bin_jobs,max_processes=4)
    
    # step 9:
    # correct binned Hi-C maps
    #corrected_file = "%s/all.ice.hm" % corrected_folder
    
    
    
    
    