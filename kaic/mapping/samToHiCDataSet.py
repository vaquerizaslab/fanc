import argparse;
from hiclib import mapping;
from hiclib.fragmentHiC import HiCdataset
from mirnylib import h5dict, genome;
import time;
import os.path;

import string
import random

def main(args):
    print args;
    
    print("Using the following settings");
    for arg, val in args.__dict__.iteritems():
        print arg, " = ", val;
        
    time.sleep(5);
    
    # read in genome object
    genome_db = genome.Genome(args.genomeFolder, readChrms=args.readChrms)
    
    # generate hicrandom temporary filename
    rs = ''.join(random.SystemRandom().choice(string.uppercase + string.digits) for _ in xrange(6))  # @UndefinedVariable
    tmpFilename = args.output + '.' + rs + '.tmp'
    
    # merge SAM files into Hi-C map
    mapped_reads = h5dict.h5dict(tmpFilename)
    
    sam_base1 = os.path.splitext(args.sam1)[0];
    sam_base2 = os.path.splitext(args.sam2)[0];
    
    print sam_base1;
    print sam_base2;
    
    time.sleep(5);
    
    mapping.parse_sam(
        sam_basename1=sam_base1,
        sam_basename2=sam_base2,
        out_dict=mapped_reads,
        genome_db=genome_db,
	    save_seqs=False
        )
    
    TR = HiCdataset(args.output, genome=genome_db)
    TR.parseInputData(dictLike=mapped_reads)
    
    os.unlink(tmpFilename)

def splitList(thisList):
    return thisList.split(",");

if __name__ == '__main__':
    parser = argparse.ArgumentParser();
    
    parser.add_argument(
        'sam1',
        help='''SAM file with first side of Hi-C reads'''
    );
    
    parser.add_argument(
        'sam2',
        help='''SAM file with second side of Hi-C reads'''
    );
    
    parser.add_argument(
        'genomeFolder',
        help='''Path to the folder with the FASTA genome files'''
    );
    
    parser.add_argument(
        'output',
        help='''Output file for Hi-C object'''
    );
    
    parser.add_argument(
        '-r', '--read-chromosomes', dest='readChrms',
        type=splitList,
        default=["#","X"],
        help='''Comma-separated list of chromosomes to read (options: #=numbered, X, Y, M). Default: #,X'''
    );
    
    main(parser.parse_args());