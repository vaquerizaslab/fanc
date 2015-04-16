'''
Created on Apr 14, 2015

@author: kkruse1
'''
from __future__ import division
import os.path
from bisect import bisect_right
from matplotlib import pyplot as plt
from warnings import warn

class Read(object):
    def __init__(self, chromosome, position, reverse=False):
        self.chromosome = chromosome
        self.position = position
        self.reverse = reverse
        self.fragmentStart = None
        self.fragmentEnd = None
        self.fragmentMid = None
    
    def setFragmentPosition(self, start, end):
        self.fragmentStart = start
        self.fragmentEnd = end
        self.fragmentMid = self.fragmentStart + int((self.fragmentEnd-self.fragmentStart)/2)
        
    def getRestrictionSiteDistance(self):
        return min(self.position-self.fragmentStart,self.fragmentEnd-self.position);

class ReadPair(object):
    def __init__(self, name):
        self.name = name
        
        self.left = None
        self.right = None
        
    def setLeftRead(self, read):
        self.left = read
        self._swapReads()
    
    def setRightRead(self, read):
        self.right = read
        self._swapReads()
        
    def _swapReads(self):
        if self.isPair() \
           and self.left.chromosome == self.right.chromosome \
           and self.left.position > self.right.position:
            tmp = self.left
            self.left = self.right
            self.right = tmp
    
    def hasLeftRead(self):
        return self.left != None
    
    def hasRightRead(self):
        return self.right != None
    
    def isPair(self):
        return self.hasRightRead() and self.hasLeftRead()
    
    def isSameChromosome(self):
        return self.isPair() and self.left.chromosome == self.right.chromosome
    
    def isSameFragment(self):
        return self.isSameChromosome() and self.left.fragmentStart == self.right.fragmentStart
    
    def getDistance(self):
        if not self.isSameChromosome():
            return None
        
        return self.right.position - self.left.position
    
    def getGapSize(self):
        if not self.isSameChromosome():
            return None
        
        if self.isSameFragment():
            return 0
        
        gap = self.right.fragmentStart - self.left.fragmentEnd
        
        if gap == 1: # neighboring fragments
            return 0
        
        return gap
        
    def isInwardPair(self):
        if self.isSameChromosome() and not self.left.reverse and self.right.reverse:
            return True
        return False
    
    def isOutwardPair(self):
        if self.isSameChromosome() and self.left.reverse and not self.right.reverse:
            return True
        return False
    
    def isSameStrandPair(self):
        if self.isSameChromosome() and self.left.reverse == self.right.reverse:
            return True
        return False



class ReadPairs(object):
    
    
    def __init__(self, genome=None):
        
        self.pairs = {}
        self.samFiles = []
        self.genome = genome
    

    def _extractReadInformation(self, line):
        x = line.rstrip()
        fields = x.split("\t")
        readId = fields[0]
        flag = int(fields[1])
        chromosome = fields[2]
        position = int(fields[3])
        
        return readId, flag, chromosome, position
    
    def loadSamFiles(self, sam1, sam2, warnings=False):
        self.loadSamFile(sam1, warnings=warnings)
        self.loadSamFile(sam2, warnings=warnings)
    
    def loadSamFile(self, sam, warnings=False):
        sam = os.path.abspath(os.path.expanduser(sam))
        self.samFiles.append(sam)
        
        leftReadMemory = set()
        
        # read reads
        with open(sam, 'r') as s:
            for x in s:
                if x.startswith("@") or len(x) == 0:
                    continue;
                readId, flag, chromosome, position = self._extractReadInformation(x)
                
                if flag == 16:
                    reverse = True
                elif flag == 0:
                    reverse = False
                else:
                    continue
                
                chrLabel = self.genome._extractChrmLabel(chromosome)
                if not chrLabel in self.genome.label2idx:
                    continue
                
                
                read = Read(chromosome, position, reverse)
                
                if not readId in self.pairs:
                    self._updateFragmentPosition(read)
                    self.pairs[readId] = ReadPair(readId)
                    self.pairs[readId].setLeftRead(read)
                    leftReadMemory.add(readId)
                else:
                    if readId in leftReadMemory:
                        if warnings == True:
                            warn("Already have left with this ID: " + readId)
                    else:
                        if self.pairs[readId].hasRightRead():
                            if warnings == True:
                                warn("Already have right read matching this pair: " + readId)
                        else:
                            self._updateFragmentPosition(read)
                            self.pairs[readId].setRightRead(read)
    
    def getInwardPairs(self, distanceCutoff=None):
        inwardPairs = {}
        for name in self.pairs:
            pair = self.pairs[name]
            if pair.isInwardPair():
                if distanceCutoff == None or pair.getGapSize() >= distanceCutoff:
                    inwardPairs[name] = pair
        return inwardPairs
    
    def getOutwardPairs(self, distanceCutoff=None):
        outwardPairs = {}
        for name in self.pairs:
            pair = self.pairs[name]
            if pair.isOutwardPair():
                if distanceCutoff == None or pair.getGapSize() >= distanceCutoff:
                    outwardPairs[name] = pair
        return outwardPairs
    
    def getSameStrandPairs(self, distanceCutoff=None):
        sameStrandPairs = {}
        for name in self.pairs:
            pair = self.pairs[name]
            if pair.isSameStrandPair():
                if distanceCutoff == None or pair.getGapSize() >= distanceCutoff:
                    sameStrandPairs[name] = pair
        return sameStrandPairs
    
    def getInterChromosomalPairs(self):
        interChromosomalPairs = {}
        for name in self.pairs:
            if not self.pairs[name].isSameChromosome():
                interChromosomalPairs[name] = self.pairs[name]
        return interChromosomalPairs
    
    def getIntraChromosomalPairs(self, distanceCutoff=None):
        intraChromosomalPairs = []
        for name in self.pairs:
            pair = self.pairs[name]
            if pair.isSameChromosome():
                if distanceCutoff == None or pair.getGapSize() >= distanceCutoff:
                    intraChromosomalPairs.append(pair)
        return intraChromosomalPairs
    
    def removeLargeDistanceToRestrictionSite(self, cutoffDistance=500):
        for name in self.pairs.keys():
            pair = self.pairs[name]
            if pair.isPair():
                if pair.left.getRestrictionSiteDistance() > cutoffDistance:
                    del self.pairs[name]
                elif pair.right.getRestrictionSiteDistance() > cutoffDistance:
                    del self.pairs[name]
    
    def removeSelfLigations(self):
        for name in self.pairs.keys():
            pair = self.pairs[name]
            if pair.isSameFragment():
                del self.pairs[name]
                
    def removeInwardReads(self, cutoffDistance=1000):
        for name in self.pairs.keys():
            pair = self.pairs[name]
            if pair.isInwardPair() and pair.getGapSize() < cutoffDistance:
                del self.pairs[name]
                
    def removeOutwardReads(self, cutoffDistance=25000):
        for name in self.pairs.keys():
            pair = self.pairs[name]
            if pair.isOutwardPair() and pair.getGapSize() < cutoffDistance:
                del self.pairs[name]
                
    def removeUnwantedLigations(self, inwardCutoff=1000, outwardCutoff=25000, reDistCutoff=500, removeSingle=True, removeSelf=True):
        inwardRemoved = 0
        outwardRemoved = 0
        reRemoved = 0
        singleRemoved = 0
        selfRemoved = 0
        inwardTotal = 0
        outwardTotal = 0
        selfTotal = 0
        total = 0
        removed = 0
        for name in self.pairs.keys():
            total += 1
            pair = self.pairs[name]
            if pair.isPair():
                if pair.isSameFragment():
                    selfTotal += 1
                    if removeSelf == True:
                        del self.pairs[name]
                        removed += 1
                        selfRemoved += 1
                elif reDistCutoff != None and (pair.left.getRestrictionSiteDistance() > reDistCutoff or
                      pair.right.getRestrictionSiteDistance() > reDistCutoff):
                    del self.pairs[name]
                    removed += 1
                    reRemoved += 1
                elif pair.isOutwardPair():
                    outwardTotal += 1
                    if outwardCutoff != None and pair.getGapSize() < outwardCutoff:
                        del self.pairs[name]
                        removed += 1
                        outwardRemoved += 1
                elif pair.isInwardPair():
                    inwardTotal += 1
                    if inwardCutoff != None and pair.getGapSize() < inwardCutoff:
                        del self.pairs[name]
                        removed += 1
                        inwardRemoved += 1
            elif removeSingle == True:
                del self.pairs[name]
                singleRemoved += 1
            elif reDistCutoff != None and (  (pair.hasLeftRead() and
                     pair.left.getRestrictionSiteDistance() > reDistCutoff) or
                    (pair.hasRightRead() and 
                     pair.right.getRestrictionSiteDistance() > reDistCutoff)  ):
                del self.pairs[name]
                removed += 1
                reRemoved += 1
        stat = ("Statistics (removed total):\n" +
                "Total\t%d\t%d\n" % (removed, total) + 
                "Single\t%d\t%d\n" % (singleRemoved, singleRemoved) +
                "RE site\t%d\t%d\n" % (reRemoved, total) +
                "Self\t%d\t%d\n" % (selfRemoved, selfTotal) +
                "Inward\t%d\t%d\n" % (inwardRemoved, inwardTotal) +
                "Outward\t%d\t%d" % (outwardRemoved, outwardTotal)
                )
        print stat
    
    def removeSingleReads(self):
        for name in self.pairs.keys():
            pair = self.pairs[name]
            if not pair.isPair():
                del self.pairs[name]
        
    def getReadPairs(self):
        truePairs = {}
        for name in self.pairs:
            if self.pairs[name].isPair():
                truePairs[name] = self.pairs[name]
        return truePairs
        
    def getIsolatedReads(self):
        isolated = {}
        for name in self.pairs:
            if not self.pairs[name].isPair():
                if self.pairs[name].hasLeftRead():
                    isolated[name] = self.pairs[name].left
                else:
                    isolated[name] = self.pairs[name].right
        return isolated
        
    def _updateFragmentPosition(self, read):
        chrLabel = self.genome._extractChrmLabel(read.chromosome)
        chrIdx = self.genome.label2idx[chrLabel]
        
        fragmentBin = bisect_right(self.genome.rsites[chrIdx], read.position)
        if fragmentBin == len(self.genome.rsites[chrIdx]):
            fragmentStart = self.genome.rsites[chrIdx][fragmentBin-1]
            fragmentEnd = self.genome.chrmLens[chrIdx]
        elif fragmentBin == 0:
            fragmentStart = 0
            fragmentEnd = self.genome.rsites[chrIdx][fragmentBin]
        else:
            fragmentStart = self.genome.rsites[chrIdx][fragmentBin-1]
            fragmentEnd = self.genome.rsites[chrIdx][fragmentBin]
        read.setFragmentPosition(fragmentStart, fragmentEnd)

        
    def plotLigationProductRatios(self, output = None, cutoffs=[0,10,20,60,100,200,600,1000,2000,6000,10000,20000,60000,100000,200000,600000,1000000]):
        inwards = []
        outwards = []
        same = []
        
        for c in cutoffs:
            total = len(self.getIntraChromosomalPairs(c))
            inwards.append(len(self.getInwardPairs(c))/total)
            outwards.append(len(self.getOutwardPairs(c))/total)
            same.append(len(self.getSameStrandPairs(c))/total)
        
        if output != None:
            plt.ioff()
        
        fig = plt.figure()
        plt.plot(cutoffs, inwards, color='r')
        plt.plot(cutoffs, outwards, color='g')
        plt.plot(cutoffs, same, color='b')
        plt.xscale('log')
        
        if output == None:
            plt.show();
        else:
            fig.savefig(output)
            plt.close(fig)
            plt.ion()
        
    def plotGapSizes(self, output=None, step=500):
        intra = self.getIntraChromosomalPairs()
        gapSizes = []
        for pair in intra:
            gapSizes.append(pair.getGapSize())
        
        if output != None:
            plt.ioff()
        
        fig = plt.figure()
        plt.hist(gapSizes, bins=range(0,max(gapSizes),step))
        
        if output == None:
            plt.show();
        else:
            fig.savefig(output)
            plt.close(fig)
            plt.ion()
        
    def plotErrorStructure(self, output=None, dataPoints=100):
        x = []
        inwardRatios = []
        outwardRatios = []

        data = []
        intra = self.getIntraChromosomalPairs()
        intra.sort(key=lambda x: x.getGapSize())
        for pair in intra:
            if len(data) < dataPoints:
                data.append(pair)
            else:
                inwards = 0
                outwards = 0
                same = 0
                mids = 0
                for dataPair in data:
                    if dataPair.isSameStrandPair():
                        same = same + 1
                    elif dataPair.isInwardPair():
                        inwards = inwards + 1
                    else:
                        outwards = outwards + 1
                    mids = mids + dataPair.getGapSize()
                x.append(mids/len(data))
                inwardRatios.append(inwards/same)
                outwardRatios.append(outwards/same)
                
                print mids/len(data), ", ", inwards/same, ", ", outwards/same
                
                data = []
        
        if output != None:
            plt.ioff()
        
        fig = plt.figure()
        fig.suptitle("Error structure by distance")
        plt.plot(x,inwardRatios, 'b', label="inward/same strand")
        plt.plot(x,outwardRatios, 'r', label="outward/same strand")
        plt.xscale('log')
        plt.axhline(y=0.5,color='black',ls='dashed')
        plt.ylim(0,3)
        plt.xlabel('gap size between fragments')
        plt.ylabel('ratio of number of reads')
        plt.legend(loc='upper right')

        if output == None:
            plt.show();
        else:
            fig.savefig(output)
            plt.close(fig)
            plt.ion()
        
    def plotRestrictionSiteDistances(self, output=None, bins=200, cutoff=500):
        distances = []
        for name in self.pairs:
            pair = self.pairs[name]
            if pair.hasLeftRead():
                distances.append(pair.left.getRestrictionSiteDistance())
            if pair.hasRightRead():
                distances.append(pair.right.getRestrictionSiteDistance())
        
        if output != None:
            plt.ioff()
        
        fig = plt.figure()
        fig.suptitle('Restriction site distances')
        plt.hist(distances, bins, log=True)
        plt.xscale('log')
        plt.axvline(x=cutoff,color='red',ls='dashed')
        plt.xlabel('Distance of read to RE site')
        plt.ylabel('Occurrence')
        
        if output == None:
            plt.show();
        else:
            fig.savefig(output)
            plt.close(fig)
            plt.ion()
            
            
            
    def filterSamFiles(self, output=None):
        assert ( output == None or
                 ( isinstance(output, list) and
                   len(output) == len(self.samFiles) ) )
            
        for i in range(0,len(self.samFiles)):
        
            path = self.samFiles[i]
            
            if output != None:
                assert len(output)-1 >= i
                outputSam = os.path.abspath(os.path.expanduser(output[i]))
            else:
                folder, fileName = os.path.split(path)
                base = os.path.splitext(fileName)[0]
                outputSam = folder + '/' + base + '.filtered.sam'
            
            print outputSam
            
            with open(path, 'r') as s:
                with open(outputSam, 'w') as o:
                    for x in s:
                        if x.startswith("@") or len(x) == 0:
                            o.write(x)
                        else:
                            readInfo = self._extractReadInformation(x)
                            if readInfo[0] in self.pairs:
                                o.write(x)
                            
