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
    

    BD.iterativeCorrectWithoutSS(force=True)
    
    BD.export("hm",args.output)
    
    


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
    
    main(parser.parse_args());
