from hiclib.fragmentHiC import HiCdataset;
import time;
import argparse;
from hiclib import binnedData
import matplotlib.pyplot as plt
import kaic.genome.genomeTools as gt
import numpy as np

def main(args):
    print("Using the following settings");
    for arg, val in args.__dict__.iteritems():
        print arg, " = ", val;
        
    time.sleep(5);
    
    # read in genome object
    genome_db = gt.loadGenomeObject(args.genome)
    
    BD = binnedData.binnedData(args.resolution, genome_db)
    BD.simpleLoad(args.input, 'hm')
    
    hm = BD.dataDict['hm']
    
    if args.absolute == False:
        nrows = hm.shape[0]
        ncols = hm.shape[1]
        ex = np.sum(hm)/(nrows*ncols)
        print "Expected: ", ex
        hm = np.log2(hm/ex)
    
    fig, ax = plt.subplots()
    hm = ax.imshow(hm, interpolation='none',aspect=1,vmin=args.min,vmax=args.max)
    plt.show()
    

def splitList(thisList):
    return thisList.split(",");

if __name__ == '__main__':
    parser = argparse.ArgumentParser();
    
    parser.add_argument(
        'input',
        help='''Input Hi-C Map file (hdf5 dict)'''
    );
    
    parser.add_argument(
        'genome',
        help='''Genome object file'''
    );
    
    parser.add_argument(
        '-res', '--resolution', dest='resolution',
        type=int,
        help='''Cutoff for filtering very large fragments''',
        required=True
    );
    
    parser.add_argument(
        '-min', dest='min',
        type=int,
        default=-3,
        help='''Lower plotting boundary'''
    );
    
    parser.add_argument(
        '-max', dest='max',
        type=int,
        default=3,
        help='''Upper plotting boundary'''
    );
    
    parser.add_argument(
        '-a', '--absolute', dest='absolute',
        action='store_true',
        help='''Plot absolute values instead of log2-fold enrichment over expectation'''
    );
    parser.set_defaults(absolute=False);
    
    main(parser.parse_args());