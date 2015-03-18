import argparse;
import time;
from hiclib.fragmentHiC import HiCdataset
import kaic.genome.genomeTools as gt

def main(args):
    print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    print("Using the following settings");
    for arg, val in args.__dict__.iteritems():
        print arg, " = ", val;
    print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
    time.sleep(5);
    
    genome_db = gt.loadGenomeObject(args.genome)
    
    merged = HiCdataset(args.output,
                        genome=genome_db,
                        mode="w")
    
    merged.merge(args.datasets)

def splitList(thisList):
    return thisList.split(",");

if __name__ == '__main__':
    parser = argparse.ArgumentParser();
    
    parser.add_argument(
        'datasets',
        type=splitList,
        help='''Mirny HiCdataset structure(s)'''
    );
    
    parser.add_argument(
        'genome',
        help='''Genome object file'''
    );
    
    parser.add_argument(
        'output',
        help='''Output filename'''
    );
    
   
    
    main(parser.parse_args());