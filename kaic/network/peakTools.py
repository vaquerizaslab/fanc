'''
Created on May 4, 2015

@author: kkruse1
'''

from __future__ import division
from kaic.plotting.plotHeatMap import plotHeatmap
from scipy.sparse import lil_matrix

def readPeakList(inFile):
    peaks = []
    with open(inFile, 'r') as peakFile:
        line = peakFile.readline().rstrip()
        print line, "\n"
        
        while line != '':
            line = peakFile.readline().rstrip()
            if line =='':
                continue
            fields = line.split("\t")
            i = int(fields[0])
            j = int(fields[1])
            inter = bool(int(fields[2]))
            obs = int(fields[3])
            fdr = float(fields[4])
            d = float(fields[5])
            ll = float(fields[6])
            v = float(fields[7])
            h = float(fields[8])
            d_fdr = float(fields[9])
            ll_fdr = float(fields[10])
            v_fdr = float(fields[11])
            h_fdr = float(fields[12])
            
            peak = {
                    'i': i,
                    'j': j,
                    'inter': inter,
                    'obs': obs,
                    'fdr': fdr,
                    'd': d,
                    'll': ll,
                    'v': v,
                    'h': h,
                    'd_fdr': d_fdr,
                    'll_fdr': ll_fdr,
                    'v_fdr': v_fdr,
                    'h_fdr': h_fdr
                    }
            peaks.append(peak)
    
    return peaks


def filterPeakList(peaks, 
                   fdr=0.1, 
                   d_fc=1.75, 
                   ll_fc=1.75, 
                   h_fc=1.5, 
                   v_fc=1.5, 
                   d_or_ll_fc=2,
                   inter_fdr=None, 
                   inter_d_fc=None, 
                   inter_ll_fc=None, 
                   inter_h_fc=None, 
                   inter_v_fc=None, 
                   inter_d_or_ll_fc=None,
                   filterFc=True,
                   filterFdr=True):
    
    if inter_fdr == None:
        inter_fdr = fdr
    if inter_d_fc == None:
        inter_d_fc = d_fc
    if inter_ll_fc == None:
        inter_ll_fc = ll_fc
    if inter_h_fc == None:
        inter_h_fc = h_fc
    if inter_v_fc == None:
        inter_v_fc = v_fc
    if inter_d_or_ll_fc == None:
        inter_d_or_ll_fc = d_or_ll_fc
    
    filteredPeaks = []
    inter = 0
    interAfter = 0
    for peak in peaks:
        if peak['inter']:
            inter += 1
            
            if filterFdr:
                # check FDRs
                if peak['d_fdr'] > inter_fdr:
                    continue
                if peak['ll_fdr'] > inter_fdr:
                    continue
                if peak['h_fdr'] > inter_fdr:
                    continue
                if peak['v_fdr'] > inter_fdr:
                    continue
            
            if filterFc:
                # check fold-changes
                if peak['d'] == 0 or peak['obs']/peak['d'] < inter_d_fc:
                    continue
                if peak['ll'] == 0 or peak['obs']/peak['ll'] < inter_ll_fc:
                    continue
                if peak['v'] == 0 or peak['obs']/peak['v'] < inter_v_fc:
                    continue
                if peak['h'] == 0 or peak['obs']/peak['h'] < inter_h_fc:
                    continue
                
                # check either/or fold changes
                if peak['obs']/peak['d'] < inter_d_or_ll_fc and peak['obs']/peak['ll'] < inter_d_or_ll_fc:
                    continue
            
            interAfter += 1
        else:
            if filterFdr:
                # check FDRs
                if peak['d_fdr'] > fdr:
                    continue
                if peak['ll_fdr'] > fdr:
                    continue
                if peak['h_fdr'] > fdr:
                    continue
                if peak['v_fdr'] > fdr:
                    continue
            
            if filterFc:
                # check fold-changes
                if peak['d'] == 0 or peak['obs']/peak['d'] < d_fc:
                    continue
                if peak['ll'] == 0 or peak['obs']/peak['ll'] < ll_fc:
                    continue
                if peak['v'] == 0 or peak['obs']/peak['v'] < v_fc:
                    continue
                if peak['h'] == 0 or peak['obs']/peak['h'] < h_fc:
                    continue
                
                # check either/or fold changes
                if peak['obs']/peak['d'] < d_or_ll_fc and peak['obs']/peak['ll'] < d_or_ll_fc:
                    continue
                    
        filteredPeaks.append(peak)
    
    print "Before filtering: %d/%d/%d peaks/intra/inter" % (len(peaks),len(peaks)-inter,inter)
    print "After filtering:  %d/%d/%d peaks/intra/inter" % (len(filteredPeaks),len(filteredPeaks)-interAfter,interAfter)
        
    return filteredPeaks

def mergePeaks(peaks, genome, dist=1, maxDiam=5):        
    M = lil_matrix((len(genome.posBinCont),len(genome.posBinCont)),dtype=bool)
    s = M.shape
    
    print "Sorting"
    peaksSort = sorted(peaks, key=lambda x: x['obs'], reverse=True)
    
    print "Building matrix"
    peakTuples = {}
    for peak in peaks:
        i = peak['i']
        j = peak['j']
        M[i,j] = True
        peakTuples[(i,j)] = peak
        
    
    print "Merging"
    
    mergedPeaks = []
    while len(peakTuples) > 0:
        print len(peakTuples)
        
        # find new maximum
        while len(peaksSort) > 0 and not M[peaksSort[0]['i'],peaksSort[0]['j']]:
            peaksSort.pop(0)
        
        maxPeak = peaksSort.pop(0)
        mergedPeak = [maxPeak]
        M[maxPeak['i'],maxPeak['j']] = False
        del peakTuples[(maxPeak['i'],maxPeak['j'])]
        
        diam = 1
        noPeakDist = 0
        while len(peakTuples) > 0 and diam < maxDiam and noPeakDist<=dist:
            foundPeak = False
            iMin = max(0,maxPeak['i']-diam)
            iMax = min(s[0]-1,maxPeak['i']+diam)
            jMin = max(0,maxPeak['j']-diam)
            jMax = min(s[0]-1,maxPeak['j']+diam)
            
            # top row
            for j in range(jMin,jMax+1):
                if M[iMin,j]:
                    mergedPeak.append(peakTuples[(iMin,j)])
                    del peakTuples[(iMin,j)]
                    M[iMin,j] = False
                    foundPeak = True
            # bottom row
            for j in range(jMin,jMax+1):
                if M[iMax,j]:
                    mergedPeak.append(peakTuples[(iMax,j)])
                    del peakTuples[(iMax,j)]
                    M[iMax,j] = False
                    foundPeak = True
            # left col
            for i in range(iMin,iMax+1):
                if M[i,jMin]:
                    mergedPeak.append(peakTuples[(i,jMin)])
                    del peakTuples[(i,jMin)]
                    M[i,jMin] = False
                    foundPeak = True
            # right col
            for i in range(iMin,iMax+1):
                if M[i,jMax]:
                    mergedPeak.append(peakTuples[(i,jMax)])
                    del peakTuples[(i,jMax)]
                    M[i,jMax] = False
                    foundPeak = True
            
            if not foundPeak:
                noPeakDist += 1
            
            diam += 1
        
        newPeak = {
                   'i': maxPeak['i'],
                   'j': maxPeak['j'],
                   'obs': maxPeak['obs'],
                   'fdr': maxPeak['fdr'],
                   'size': len(mergedPeak),
                   'diam': diam
                   }
        mergedPeaks.append(newPeak)
    return mergedPeaks
        
def plotPeaks(hic, genomeFile, resolution, peaks, minIx=None, maxIx=None, vmin=-3, vmax=3):
    peakTuples = []
    for peak in peaks:
        peakTuples.append((peak['i'],peak['j']))
    
    plotHeatmap(hic, genomeFile, resolution, iStartIndex=minIx, iEndIndex=maxIx, jStartIndex=minIx, jEndIndex=maxIx, highlightPixels=peakTuples, vmin=vmin, vmax=vmax)
    
    
    
        
        
        
        
        
    
    
    