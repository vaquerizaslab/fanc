import time;
import argparse;
import kaic.genome.genomeTools as gt
from hiclib import binnedData

def main(args):
    print("Using the following settings");
    for arg, val in args.__dict__.iteritems():
        print arg, " = ", val;

    time.sleep(5);

    # read in genome object
    genome_db = gt.loadGenomeObject(args.genome)


    BD = binnedData.binnedData(args.resolution, genome_db)
    BD.simpleLoad(args.input, 'hm')
    
    
    if args.diagonal == True:
        BD.removeDiagonal()   #we never ever use diagonal
    if args.sf > 0:
        BD.removeBySequencedCount(args.sf)
    if args.poor > 0:
        BD.removePoorRegions(cutoff=args.poor)
    if args.trunc > 0:
        BD.truncTrans(high=args.trunc)
    
    BD.export("hm",args.output,byChromosome=True)
    
    


def splitList(thisList):
    return thisList.split(",");

if __name__ == '__main__':
    parser = argparse.ArgumentParser();

    parser.add_argument(
        'input',
        help='''Input Hi-C Map file (binned data, hdf5 dict)'''
    );

    parser.add_argument(
        'genome',
        help='''Genome object file'''
    );
    
    parser.add_argument(
        'output',
        help='''Output file'''
    );

    parser.add_argument(
        '-res', '--resolution', dest='resolution',
        type=int,
        help='''Cutoff for filtering very large fragments''',
        required=True
    );

    parser.add_argument(
        '-d', '--diagonal', dest='diagonal',
        action='store_true',
        help='''Remove diagonal'''
    );
    parser.set_defaults(diagonal=False);
    
    parser.add_argument(
        '-c', '--sequenced-count', dest='sq',
        action='store_true',
        help='''Remove bins that have less than sequencedFraction*resolution sequenced counts.'''
    );
    parser.set_defaults(sq=False);
    
    parser.add_argument(
        '-sf', '--sequenced-fraction', dest='sf',
        type=float,
        default=0.0,
        help='''Remove bins that have less than sequencedFraction*resolution sequenced counts.'''
    );
    
    parser.add_argument(
        '-p', '--poor', dest='poor',
        type=float,
        default=0.0,
        help='''Remove this percent of bins with least counts'''
    );
    
    parser.add_argument(
        '-t', '--truncate-trans', dest='trunc',
        type=float,
        default=0.0,
        help='''Truncates this fraction of trans contacts to remove blowouts'''
    );
    

    main(parser.parse_args());
