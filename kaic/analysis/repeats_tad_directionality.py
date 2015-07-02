#!/usr/bin/env python


'''
Created on Jun 5, 2015

@author: kkruse1
'''

import argparse
import kaic.data.genomic as gd
import kaic.plotting.plot_genomic_data as pgd
import os.path
from kaic.genome.genomeTools import loadGenomeObject
from kaic.hrandom.genomic import random_bed

def splitList(thisList):
    return thisList.split(",");

if __name__ == '__main__':
    parser = argparse.ArgumentParser();
    
    parser.add_argument(
        'input',
        nargs='+',
        help='''Bed files'''
    )
    
    parser.add_argument(
        '-res', '--resolution', dest='resolution',
        type=int,
        help='''Resolution of the Hi-C data set''',
        required=True
    )

    parser.add_argument(
        '-w', '--window', dest='window',
        type=int,
        help='''Window size for enrichment plots''',
        required=True
    )
    
    parser.add_argument(
        '-b', '--bins', dest='bins',
        type=int,
        help='''Number of bins per plot''',
        default=2000
    )
    
    parser.add_argument(
        '-o', '--output-file', dest='output',
        help='''Output file''',
        required=True
    )
    
    parser.add_argument(
        '-hic', '--hic-file', dest='hic_file',
        help='''TADs BED file''',
        required=True
    )
    
    parser.add_argument(
        '-s', '--scale', dest='scale',
        action='store_false',
        help='''DO NOT rescale enrichment scores to [0,1]'''
    );
    parser.set_defaults(scale=True)

    args = parser.parse_args()
    
    print "Hi-C file: %s" % args.hic_file 
    hic = gd.Hic(args.hic_file)
    
    print "Calculating directionality index"
    di = hic.directionality(args.resolution)
    
    
    
    if args.genome is not None:
        genome = loadGenomeObject(args.genome)
    
    beds = []
    for bed_file in args.input:
        base = os.path.splitext(os.path.basename(bed_file))[0]
        bed = gd.BedImproved.from_bed_file(bed_file,name=base)
        beds.append(bed)
        
        if args.genome is not None:
            bed_rand = random_bed(genome, bed, name=base + "_random")
            beds.append(bed_rand)
        
    #bd = pgd.BedDistribution(tads, beds, window_size=args.window, n_bins=args.bins, rescale=args.scale)
    #bd.show(args.output)
        
    
        
    