from hiclib.fragmentHiC import HiCdataset;
import time;
import argparse;
from hiclib import binnedData
from hiclib import highResBinnedData
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
    
    if args.lr == True:
        BD = binnedData.binnedData(args.resolution, genome_db)
        BD.simpleLoad(args.input, 'hm', chromosomeOrder=args.order)
        #BD.simpleLoad(args.input, 'hm', chromosomeOrder=[11])
        hm = BD.dataDict['hm']
    else:
        BD = highResBinnedData.HiResHiC(args.resolution, genome_db)
        BD.loadData(args.input)
        hm = BD.getCombinedMatrix();
        
    
    
    
    if args.absolute == False:
        #BD.removeZeros()
        #hm = BD.dataDict['hm']
        nrows = hm.shape[0]
        ncols = hm.shape[1]
        ex = np.sum(hm)/(nrows*ncols)
        print "Expected: ", ex
        hm = np.log2(hm/ex)
        #BD.dataDict['hm'] = hm
        #BD.restoreZeros()
        #hm = BD.dataDict['hm']
    
    fig, ax = plt.subplots()
    myPlot = ax.imshow(hm, interpolation='none',aspect=1,vmin=args.min,vmax=args.max)
    myPlot.set_cmap(args.colormap)
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
        '-o', '--order', dest='order',
        type=splitList,
        default=None,
        help='''comma-separated list of chromosome names for heatmap order'''
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
    
    parser.add_argument(
        '-c', '--color-map', dest='colormap',
        default='YlGnBu',
        help='''Matplotlib color map name to use'''
    );
    
    parser.add_argument(
        '-l', '--low-res', dest='lr',
        action='store_true',
        help='''Use low-resolution analysis'''
    );
    parser.set_defaults(lr=False);
    
    parser.set_defaults(absolute=False);
    
    main(parser.parse_args());