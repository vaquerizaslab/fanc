import time;
import argparse;
import kaic.genome.genomeTools as gt
from hiclib import binnedData
import matplotlib.pyplot as plt
import numpy as np
from bisect import bisect_right
import random
from matplotlib.backends.backend_pdf import PdfPages

def findBins(start, end, chromosome, genome):
    if isinstance( chromosome, (int, long) ):
        chrIdx = chromosome
    else:
        chrIdx = genome.label2idx[chromosome]
    
    mask = genome.chrmIdxBinCont == chrIdx
    chrBinStart = genome.chrmStartsBinCont[chrIdx]
    binStarts = genome.posBinCont[mask]
    
    startBin = bisect_right(binStarts, start) + chrBinStart;
    endBin = bisect_right(binStarts, end) + chrBinStart;
    
    return startBin, endBin

def randomPosition(length, chromosome, genome):
    if isinstance( chromosome, (int, long) ):
        chrIdx = chromosome
    else:
        chrIdx = genome.label2idx[chromosome]
        
    chrLen = genome.chrmLens[chrIdx]-length
    
    start = random.randint(0,chrLen);
    end = start+length
    
    return start, end, chromosome;
    

def compare(genome, resolution, hicMap, positions1, positions2 = None, output = None):
    print("Using the following settings");
    print "genome = ", genome
    print "hicMap = ", hicMap
    print "resolution = ", resolution
    print "positions1 = ", positions1
    print "position2 = ", positions2 
        
    time.sleep(5);
    
    # read in genome object
    genome_db = gt.loadGenomeObject(genome)
    
    # read in Hi-C map
    BD = binnedData.binnedData(resolution, genome_db)
    BD.simpleLoad(hicMap, 'hm')
    hm = np.array(BD.dataDict['hm'])
    
    nrows = hm.shape[0]
    ncols = hm.shape[1]
    ex = np.sum(hm)/(nrows*ncols)
    
    x = range(-10,11)
    ix = range(0,21)
    
    # read in element file
    enrichment1 = []
    enrichment2 = []
    randomEnrichment1 = []
    randomEnrichment2 = []
    for i in ix:
        enrichment1.append( [] )
        enrichment2.append( [] )
        randomEnrichment1.append( [] )
        randomEnrichment2.append( [] )
        
    with open(positions1) as f:
        for line in f:
            d = line.rstrip().split('\t');
            
            start = int(d[1])
            end = int(d[2])
            chrm = d[0]
            
            startBin, endBin = findBins(start, end, chrm, genome_db)
            randomPos = randomPosition(end-start,chrm,genome_db)
            randStartBin, randEndBin = findBins(randomPos[0], randomPos[1], randomPos[2], genome_db)
            
            total = np.sum(hm[startBin,])
            
            for i in ix:
                enrichment1[i].append( hm[startBin,startBin+x[i]]/ex )
                randomEnrichment1[i].append( hm[randStartBin,randStartBin+x[i]]/ex )
    
    if positions2 != None:
        with open(positions2) as f:
            for line in f:
                d = line.rstrip().split('\t');
                
                start = int(d[1])
                end = int(d[2])
                chrm = d[0]
                
                startBin, endBin = findBins(start, end, chrm, genome_db)
                randomPos = randomPosition(end-start,chrm,genome_db)
                randStartBin, randEndBin = findBins(randomPos[0], randomPos[1], randomPos[2], genome_db)
                
                total = np.sum(hm[startBin,])
                
                for i in ix:
                    enrichment2[i].append( hm[startBin,startBin+x[i]]/ex )
                    randomEnrichment2[i].append( hm[randStartBin,randStartBin+x[i]]/ex )
        
        enrichment2[10] = enrichment2[9]
        randomEnrichment2[10] = randomEnrichment2[9]
    
    enrichment1[10] = enrichment1[9]
    randomEnrichment1[10] = randomEnrichment1[9]
    
    
    
#
#     fig, ax = plt.subplots()
#     hm = ax.imshow(BD.dataDict['hm'], interpolation='none',aspect=1,vmin=args.min,vmax=args.max)
#     plt.show()
    
    fig, ax = plt.subplots(nrows=1,ncols=1,sharex=True)
    ax.errorbar(x, np.median(enrichment1,axis=1), yerr=np.std(enrichment1,axis=1), label='pos1')
    ax.errorbar(x, np.median(randomEnrichment1,axis=1), yerr=np.std(randomEnrichment1,axis=1), label='rand1')
    if positions2 != None:
        ax.errorbar(x, np.median(enrichment2,axis=1), yerr=np.std(enrichment2,axis=1), label='pos2')
        ax.errorbar(x, np.median(randomEnrichment2,axis=1), yerr=np.std(randomEnrichment2,axis=1), label='rand2')
    plt.legend(loc='upper left')
    
    if output == None:
        plt.show();
    else:
        pp = PdfPages(output)
        pp.saveFig();
    

def splitList(thisList):
    return thisList.split(",");

if __name__ == '__main__':
    parser = argparse.ArgumentParser();
    
    parser.add_argument(
        'hicMap',
        help='''Input Hi-C Map file (hdf5 dict)'''
    );
    
    parser.add_argument(
        'genome',
        help='''Genome object file'''
    );
    
    parser.add_argument(
        'elementPositions',
        help='''Tab-separated '''
    );
    
    parser.add_argument(
        '-res', '--resolution', dest='resolution',
        type=int,
        help='''Datset resolution''',
        required=True
    );
    
    parser.add_argument(
        '-cmp', '--compare', dest='elementPositions2',
        default= None,
        help='''Comparison data set (element positions)'''
    );
    
    parser.add_argument(
        '-o', '--output', dest='output',
        default= None,
        help='''Output file (pdf) - suppresses plotting window'''
    );
    
    args = parser.parse_args()
    
    if args.output == None:
        if args.elementPositions2 != None:
            compare(args.genome, args.resolution, args.hicMap, args.elementPositions, positions2 = args.elementPositions2)
        else:
            compare(args.genome, args.resolution, args.hicMap, args.elementPositions)
    else:
        if args.elementPositions2 != None:
            compare(args.genome, args.resolution, args.hicMap, args.elementPositions, positions2 = args.elementPositions2, output = args.output)
        else:
            compare(args.genome, args.resolution, args.hicMap, args.elementPositions, output = args.output)