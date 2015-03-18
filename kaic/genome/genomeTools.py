import argparse;
import time;
import cPickle as pickle;
from mirnylib import genome;

#
# TODO FUNCTIONS HERE
#

def createGenomeObject(folder, re, readChrms):
    # read in genome object
    genome_db = genome.Genome(folder, readChrms=readChrms)
    
    if args.restrictionEnzyme != '':
        genome_db.setEnzyme(re);
        
    return genome_db;

def saveGenomeObject(genome_db, output):
    time.sleep(5);
    
    with open(output, 'wb') as o:
        pickle.dump(genome_db, o, pickle.HIGHEST_PROTOCOL)
    
def loadGenomeObject(inFile):
    with open(inFile, 'rb') as i:
        genome_db = pickle.load(i)
    return genome_db


def splitList(thisList):
    return thisList.split(",");

if __name__ == '__main__':
    parser = argparse.ArgumentParser();
    
    parser.add_argument(
        'input',
        help='''Genome folder with FASTA files'''
    );
    
    parser.add_argument(
        'output',
        help='''Output file for genome object'''
    );
    

    parser.add_argument(
        '-re', '--restriction-enzyme', dest='restrictionEnzyme',
        default='',
        help='''Restriction enzyme name (e.g. HindIII)''',
        required=True
    );
    
    
    parser.add_argument(
        '-r', '--read-chromosomes', dest='readChrms',
        type=splitList,
        default=["#","X"],
        help='''Comma-separated list of chromosomes to read (options: #=numbered, X, Y, M). Default: #,X'''
    );
    
    
    args = parser.parse_args()
    
    print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    print("Using the following settings");
    for arg, val in args.__dict__.iteritems():
        print arg, " = ", val;
    print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
    
    genome = createGenomeObject(args.input, args.restrictionEnzyme, args.readChrms)
    saveGenomeObject(genome, args.output);