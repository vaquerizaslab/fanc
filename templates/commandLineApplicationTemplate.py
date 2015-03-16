import argparse;
import time;


#
# TODO FUNCTIONS HERE
#

def main(args):
    print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    print("Using the following settings");
    for arg, val in args.__dict__.iteritems():
        print arg, " = ", val;
    print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
    time.sleep(5);
    
    #
    # TODO MAIN CODE HERE TODO
    #


def splitList(thisList):
    return thisList.split(",");

if __name__ == '__main__':
    parser = argparse.ArgumentParser();
    
    # INPUT EXAMPLE
    #parser.add_argument(
    #    'sam1',
    #    help='''SAM file with first side of Hi-C reads'''
    #);
    
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