from hiclib.fragmentHiC import HiCdataset;
import time;
import argparse;
from hiclib import binnedData
import matplotlib.pyplot as plt
import kaic.genome.genomeTools as gt


def main(args):
    print("Using the following settings");
    for arg, val in args.__dict__.iteritems():
        print arg, " = ", val;
        
    time.sleep(5);
    
    # read in genome object
    genome_db = gt.loadGenomeObject(args.genome)
    
    if args.restrictionEnzyme != '':
        genome_db.setEnzyme(args.restrictionEnzyme);
    
    BD = binnedData.binnedData(args.resolution, genome_db)
    BD.simpleLoad(args.input, 'hm')
    
    fig, ax = plt.subplots()
    hm = ax.imshow(BD.dataDict['hm'], interpolation='none',aspect=1,vmin=args.min,vmax=args.max)
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
        'genomeFolder',
        help='''Genome folder with FASTA files'''
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
        default=0,
        help='''Lower plotting boundary'''
    );
    
    parser.add_argument(
        '-max', dest='max',
        type=int,
        default=3,
        help='''Upper plotting boundary'''
    );
    
    main(parser.parse_args());