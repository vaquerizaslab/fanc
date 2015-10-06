#!/usr/bin/env python

import argparse
import os.path
from kaic.data.genomic import Chromosome, Genome, HicBasic

if __name__ == '__main__':
    parser = argparse.ArgumentParser();
    
    def split_list(my_list):
        return my_list.split(",")
    
    parser.add_argument(
        'output',
        help='''File name of the output Hi-C file'''
    )
    
    parser.add_argument(
        'chromosome_fasta',
        type=split_list,
        help='''Comma-separated list of FASTA files with chromosome data'''
    )
    
    parser.add_argument(
        'resolution',
        type=int,
        help='''Resolution of the heatmap (int)'''
    )
    
    parser.add_argument(
        'intra',
        help='''Root folder of the intra-chromosomal data'''
    )
    
    parser.add_argument(
        'inter',
        nargs="?",
        help='''Root folder of the inter-chromosomal data'''
    )
    
    args = parser.parse_args()
    
    output_path = os.path.expanduser(args.output)
    intra_path = os.path.expanduser(args.intra)
    inter_path = None
    if args.inter:
        inter_path = os.path.expanduser(args.inter)
    
    # convert resolution to string
    resolution_dict = {
        1000000: '1mb',
        100000: '100kb',
        10000: '10kb',
        250000: '250kb',
        25000: '25kb',
        500000: '500kb',
        50000: '50kb',
        5000: '5kb'
    }
    resolution_string = resolution_dict[args.resolution]
    
    # build genome
    genome = Genome()
    chromosome_names = []
    for chromosome_file in args.chromosome_fasta:
        chromosome_file = os.path.expanduser(chromosome_file)
        
        chromosome = Chromosome.from_fasta(chromosome_file)
        genome.add_chromosome(chromosome)
        chromosome_names.append(chromosome.name)
    
    # get nodes and build reference dictionary
    nodes = genome.get_regions(args.resolution)
    node_dict = {}
    for ix, node in enumerate(nodes):
        if not node_dict.has_key(node.chromosome):
            node_dict[node.chromosome] = {}
        node_dict[node.chromosome][node.start-1] = ix
        
    # build Hi-C data set
    hic = HicBasic(file_name=output_path)
    hic.add_nodes(nodes)
    
    # intra-chromosomal data
    data_sets = [intra_path]
    if inter_path:
        data_sets.append(inter_path)
    
    for path in data_sets:
        for chromosome_name in chromosome_names:
            norm_file = "%s/%s/MAPQGE30/%s_%s.KRnorm" % (path,chromosome_name,chromosome_name,resolution_string)
            rawo_file = "%s/%s/MAPQGE30/%s_%s.RAWobserved" % (path,chromosome_name,chromosome_name,resolution_string)
            
            # read norm vector
            norm_dict = {}
            with open(norm_file, 'r') as norm:
                current = 0
                for line in norm:
                    line = line.rstrip()
                    v = float(line)
                    
                    norm_dict[current] = v
                    
                    current += args.resolution
            
            with open(rawo_file, 'r') as raw:
                for line in raw:
                    start1, start2, score = line.rstrip().split("\t")
                    ix1 = node_dict[int(start1)]
                    ix2 = node_dict[int(start2)]
                    norm_score = float(score)/norm_dict[int(start1)]/norm_dict[int(start2)]
                    
                    hic.add_edge([ix1,ix2,norm_score], flush=False)
    
    hic.flush()
    