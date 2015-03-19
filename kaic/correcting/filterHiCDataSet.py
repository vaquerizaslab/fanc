from hiclib.fragmentHiC import HiCdataset;
import time;
import argparse;
import kaic.genome.genomeTools as gt


def main(args):
    print("Using the following settings");
    for arg, val in args.__dict__.iteritems():
        print arg, " = ", val;
        
    time.sleep(5);
    
    # read in genome object
    genome_db = gt.loadGenomeObject(args.genome)
    

    if args.memory == True:
        print "Mapping in memory";
        TR = HiCdataset("memory", inMemory=True, genome=genome_db);
    else:
        print "Mapping directly to file";
        TR = HiCdataset(args.output, genome=genome_db);
        
    TR.load(args.input);
    
    # filter duplicates
    if args.duplicates:
        print "Filtering duplicates"
        TR.filterDuplicates();
    
    # filter large and small fragments
    if args.large > -1 and args.small > -1:
        print "Filtering very small (<", args.small, ") and very large (>", args.large, ") fragments"
        TR.filterLarge(cutlarge=args.large,cutsmall=args.small);
    else:
        if args.large > -1:
            print "Filtering only very large (>", args.large, ") fragments"
            TR.filterLarge(cutlarge=args.large,cutsmall=0);
        if args.small > -1:
            print "Filtering only very small (<", args.small, ") fragments"
            TR.filterLarge(cutlarge=9999999999,cutsmall=args.small);
    
    # filter extremely mapped fragments
    if args.extremelow > 0 and args.extremehigh > 0:
        print "Filtering fraction of fragments with very low (", args.extremelow, ") and very high number of contacts (", args.extremehigh, ")"
        TR.filterExtreme(cutH=args.extremehigh,cutL=args.extremelow);
    else:
        if args.extremelow > 0:
            print "Filtering fraction of fragments with very low (", args.extremelow, ")  number of contacts"
            TR.filterExtreme(cutH=0,cutL=args.extremelow);
        if args.extremehigh > 0:
            print "Filtering fraction of fragments with very high number of contacts (", args.extremehigh, ")"
            TR.filterExtreme(cutH=args.extremehigh,cutL=0);
    
    # filter 
    if args.restrictionOffset > 0:
        print "Removing reads that start within ", args.restrictionOffset, "bp from the RE site"
        TR.filterRsiteStart(offset=args.restrictionOffset);
    
    # save now if we hadn't saved earlier
    if args.memory == True:
        print "Saving to file";
        TR.save(args.output);

def splitList(thisList):
    return thisList.split(",");

if __name__ == '__main__':
    parser = argparse.ArgumentParser();
    
    parser.add_argument(
        'input',
        help='''Input Hi-C data set (not binned)'''
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
        '-m', '--memory', dest='memory',
        action='store_true',
        help='''Perform operations in memory, only then save'''
    );
    
    parser.add_argument(
        '-d', '--duplicates', dest='duplicates',
        action='store_true',
        help='''Enable duplicate read filtering'''
    );
    parser.set_defaults(duplicates=False);
    
    parser.add_argument(
        '-l', '--large', dest='large',
        type=int,
        default=-1,
        help='''Cutoff for filtering very large fragments'''
    );
    
    parser.add_argument(
        '-s', '--small', dest='small',
        type=int,
        default=-1,
        help='''Cutoff for filtering very small fragments'''
    );
    
    parser.add_argument(
        '-el', '--extreme-low', dest='extremelow',
        type=float,
        default=0.0,
        help='''Cutoff (fraction of fragments) for filtering fragments with a very small number of reads'''
    );
    
    parser.add_argument(
        '-eh', '--extreme-high', dest='extremehigh',
        type=float,
        default=0.0,
        help='''Cutoff (fraction of fragments) for filtering fragments with a very high number of reads'''
    );
    
    parser.add_argument(
        '-ro', '--restriction-offset', dest='restrictionOffset',
        type=int,
        default=0,
        help='''Cutoff (in base-pairs) for filtering reads located less than this value from an the RE site'''
    );
    

    
    main(parser.parse_args());