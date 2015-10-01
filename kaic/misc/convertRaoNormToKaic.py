#!/usr/bin/env python

import argparse;
import os.path
from kaic.data.genomic import Chromosome, Genome, HicBasic

'''
Created on May 20, 2015

@author: kkruse1
'''



if __name__ == '__main__':
    parser = argparse.ArgumentParser();
    
    parser.add_argument(
        'raw',
        help='''Rao 2012 et al. raw contact counts'''
    )
    
    parser.add_argument(
        'norm',
        help='''Rao 2012 et al. KR normalized contact counts'''
    )
    
    parser.add_argument(
        'chromosome',
        help='''FASTQ file of the corresponding chromosome'''
    );
    
    parser.add_argument(
        'resolution',
        type=int,
        help='''Resolution of the file (e.g. 100000)'''
    );
    
    parser.add_argument(
        'output',
        help='''Output file'''
    );
    
    
    args = parser.parse_args()
    raw_file = os.path.expanduser(args.raw)
    norm_file = os.path.expanduser(args.norm)
    out_file = os.path.expanduser(args.output)
    chrom_file = os.path.expanduser(args.chromosome)
    
    # read chromosome
    chromosome = Chromosome.from_fasta(chrom_file)
    genome = Genome(chromosomes=[chromosome])
    nodes = genome.get_regions(args.resolution)
    node_dict = {}
    for ix, node in enumerate(nodes):
        print node.start-1
        node_dict[node.start-1] = ix
        
    print node_dict
    
    # build Hi-C data set
    hic = HicBasic(file_name=out_file)
    hic.add_nodes(nodes)
    
    
    # read norm vector
    norm_dict = {}
    with open(norm_file, 'r') as norm:
        current = 0
        for line in norm:
            line = line.rstrip()
            v = float(line)
            
            norm_dict[current] = v
            
            current += args.resolution
    
    with open(raw_file, 'r') as raw:
        for line in raw:
            start1, start2, score = line.rstrip().split("\t")
            ix1 = node_dict[int(start1)]
            ix2 = node_dict[int(start2)]
            norm_score = float(score)/norm_dict[int(start1)]/norm_dict[int(start2)]
            
            hic.add_edge([ix1,ix2,norm_score], flush=False)
        hic.flush()
            
    hic.close()
