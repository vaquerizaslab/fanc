import argparse;
import time;
from hiclib.fragmentHiC import HiCdataset
from mirnylib import genome;

def main(args):
    print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    print("Using the following settings");
    for arg, val in args.__dict__.iteritems():
        print arg, " = ", val;
    print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
    time.sleep(5);
    
    genome_db = genome.Genome(args.genomeFolder, readChrms=args.readChrms)
    
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
        'genomeFolder',
        help='''Path to the folder with the FASTA genome files'''
    );
    
    parser.add_argument(
        'output',
        help='''Output filename'''
    );
    
    parser.add_argument(
        '-r', '--read-chromosomes', dest='readChrms',
        type=splitList,
        default=["#","X"],
        help='''Comma-separated list of chromosomes to read (options: #=numbered, X, Y, M). Default: #,X'''
    );
    
    # NUMERIC OPTION EXAMPLE
    #parser.add_argument(
    #    '-m', '--min-length', dest='minLength',
    #    type=int,
    #    default=25,
    #    help='''Minimum sequence length to attempt the mapping'''
    #);
    
    # CUSTOM OPTION EXAMPLE
    #parser.add_argument(
    #    '-r', '--read-chromosomes', dest='readChrms',
    #    type=splitList,
    #    default=["#","X"],
    #    help='''Comma-separated list of chromosomes to read (options: #=numbered, X, Y, M). Default: #,X'''
    #);
    
    # BOOLEAN OPTION EXAMPLE
    #parser.add_argument(
    #    '-nf', '--no-filter', dest='filter',
    #    action='store_false',
    #    help='''Do not filter output SAM file by mappability, mapping quality, and uniqueness'''
    #);
    #parser.set_defaults(filter=True);
    
    main(parser.parse_args());