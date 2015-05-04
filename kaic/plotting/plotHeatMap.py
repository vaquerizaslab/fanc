from hiclib.fragmentHiC import HiCdataset;
import time;
import argparse;
from hiclib import binnedData
from hiclib import highResBinnedData
import matplotlib
import matplotlib.pyplot as plt
import kaic.genome.genomeTools as gt
import numpy as np



def main(args):
    print("Using the following settings");
    for arg, val in args.__dict__.iteritems():
        print arg, " = ", val;
        
    time.sleep(5);
    
    plotHeatmap(args.input, args.genome, args.resolution, args.chromosomeOrder, args.absolute, args.colormap, args.min, args.max, args.lr)
    
def plotHeatmap(binnedData, 
                genome, 
                resolution, 
                chromosomeOrder=None, 
                absolute=False, 
                colormap=False, 
                vmin=-3, 
                vmax=3, 
                lr=False,
                iStartIndex=None,
                iEndIndex=None,
                jStartIndex=None,
                jEndIndex=None):
    
    
    
    # read in genome object
    genome_db = gt.loadGenomeObject(genome)
    
    if iStartIndex == None:
        iStartIndex = 0
    if iEndIndex == None:
        iEndIndex = len(genome_db.posBinCont)-1
    if jStartIndex == None:
        jStartIndex = 0
    if jEndIndex == None:
        jEndIndex = len(genome_db.posBinCont)-1
    
    if lr == True:
        BD = binnedData.binnedData(resolution, genome_db)
        BD.simpleLoad(binnedData, 'hm', chromosomeOrder=chromosomeOrder)
        hm = BD.dataDict['hm']
        hm = hm[iStartIndex:iEndIndex,jStartIndex:jEndIndex]
    else:
        BD = highResBinnedData.HiResHiC(genome_db, resolution)
        BD.loadData(binnedData)
        
        hm = np.zeros((iEndIndex-iStartIndex, jEndIndex-jStartIndex), float)
        for chr1, chr2 in BD:
            M = BD.data[(chr1, chr2)].getData()
            
            # to convert chromosome into global (genome-wide) indices
            chr1StartBin = BD.genome.chrmStartsBinCont[chr1]
            chr2StartBin = BD.genome.chrmStartsBinCont[chr2]
            chr1EndBin = BD.genome.chrmEndsBinCont[chr1]
            chr2EndBin = BD.genome.chrmEndsBinCont[chr2]
            
            if chr1StartBin > iEndIndex or chr1EndBin < iStartIndex:
                continue
            if chr2StartBin > jEndIndex or chr2EndBin < jStartIndex:
                continue
            
            iStart = max(0, iStartIndex-chr1StartBin)
            iEnd = min(M.shape[0]-1,chr1EndBin-iStartIndex)
            jStart = max(0, jStartIndex-chr2StartBin)
            jEnd = min(M.shape[1]-1,chr2EndBin-jStartIndex)
            
            hm[iStart+chr1StartBin:iEnd+chr1StartBin, jStart+chr2StartBin:jEnd+chr2StartBin] = M[iStart:iEnd,jStart:jEnd]
            
            
        
        # hm = BD.getCombinedMatrix();
        
    
    
    
    if absolute == False:
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
    
    cdict = {'red': ((0.0, 1.0, 1.0),
                    (0.28, 0.18, 0.18),
                    (0.72, 0.78, 0.78),
                    (1.0, 1.0, 1.0)),
            'green': ((0.0, 1.0, 1.0),
                    (0.36, 0.05, 0.05),
                    (0.49, 0.12, 0.12),
                    (1.0, 1.0, 1.0)),
            'blue': ((0.0, 1.0, 1.0),
                    (0.26, 0.62, 0.62),
                    (0.37, 0.5, 0.5),
                    (0.77, 0.2, 0.2),
                    (0.92, 0.64, 0.64),
                    (1.0, 1.0, 1.0))
            }
    cmap = matplotlib.colors.LinearSegmentedColormap("Sexton colormap", cdict, 256)
    
    fig, ax = plt.subplots()
    myPlot = ax.imshow(hm, interpolation='none',aspect=1,norm=matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax))
    if colormap == None:
        myPlot.set_cmap(cmap)
    else:
        myPlot.set_cmap(colormap)
    plt.show()
