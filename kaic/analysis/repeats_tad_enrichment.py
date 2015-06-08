#!/usr/bin/env python


'''
Created on Jun 5, 2015

@author: kkruse1
'''

import argparse
import kaic.data.genomic as gd
import kaic.plotting.plot_genomic_data as pgd
from kaic.tools.files import make_dir
import os.path

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
        '-t', '--tad-file', dest='tad_file',
        help='''TADs BED file''',
        required=True
    )

    args = parser.parse_args()
    
    print "TAD file: %s" % args.tad_file 
    tads = gd.Bed(args.tad_file)
    
    beds = []
    for bed_file in args.input:
        
        bed = gd.Bed(bed_file)
        beds.append(bed)
        
    bd = pgd.BedDistribution(tads, beds, window_size=args.window, n_bins=args.bins)
    bd.show(args.output)
        

    