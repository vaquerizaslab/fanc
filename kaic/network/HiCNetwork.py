'''
Created on Apr 23, 2015

@author: kkruse1
'''

from __future__ import division
from hiclib import highResBinnedData
import kaic.genome.genomeTools
import graph_tool as gt
import tables as t
from scipy.stats import poisson
import numpy as np
from scikits.statsmodels.sandbox.stats.multicomp import multipletests
from __builtin__ import classmethod
from collections import Counter

class Edge(t.IsDescription):
    idNumber = t.Int64Col()   # @UndefinedVariable
    source   = t.Int32Col()   # @UndefinedVariable
    sink     = t.Int32Col()   # @UndefinedVariable
    weight   = t.Float32Col() # @UndefinedVariable
    inter    = t.BoolCol()    # @UndefinedVariable
    p        = t.Float64Col() # @UndefinedVariable
    fdr      = t.Float64Col() # @UndefinedVariable

class Locus(t.IsDescription):
    idNumber     = t.Int64Col() # @UndefinedVariable
    chromosomeId = t.Int64Col() # @UndefinedVariable
    start        = t.Int64Col() # @UndefinedVariable
    end          = t.Int64Col() # @UndefinedVariable
    mappable     = t.BoolCol() # @UndefinedVariable
    
class Chromosome(t.IsDescription):
    idNumber = t.Int64Col()  # @UndefinedVariable
    name     = t.StringCol(16) # @UndefinedVariable
    length   = t.Int64Col()  # @UndefinedVariable



class HiCNetwork(object):
    def __init__(self, data=None, fdrCutoff=0.01):
        self.data = data
#         self.file = data
#         self.network = self.file.get_node("/", 'network')
#         self.nodes = self.file.get_node(self.network, 'nodes')
#         self.edges = self.file.get_node(self.network, 'edges')
#         self.genome = self.file.get_node("/", 'genome')
#         self.chromosomes = self.file.get_node(self.genome, 'chromosomes')
#         
#         # turn into graph
#         self.g = gt.Graph(directed=False)
#         for node in self.nodes:
#             n = self.g.add_vertex()
        #for edge in self.edges:
            
        
        
    
    def __del__(self):
        print "Closing hdf5 file"
        try:
            self.file.close()
        except AttributeError:
            print "Nothing to close"
        
    
    def _initializeData(self, fileName=None):
        if not hasattr(self, 'file'):
            self.file = t.open_file(fileName, mode = "w", title = "HiCNetwork")
        
        # network
        if not hasattr(self, 'network'):
            self.network = self.file.create_group("/", 'network', 'Hi-C network')
        if not hasattr(self, 'nodes'):
            self.nodes = self.file.create_table(self.network, 'nodes', Locus, 'Network nodes')
        if not hasattr(self, 'edges'):
            self.edges = self.file.create_table(self.network, 'edges', Edge, 'Network edges')
        
        # genome
        if not hasattr(self, 'genome'):
            self.genome = self.file.create_group("/", 'genome', 'Genome information')
        if not hasattr(self, 'chromosomes'):
            self.file.create_table(self.genome, 'chromosomes', Chromosome, 'Chromosome information')
    
    @classmethod
    def liebermannBinnedHiC(cls, fileName, hic, genome=None, resolution=None):
        if type(genome) == str:
            genome = kaic.genome.genomeTools.loadGenomeObject(genome)
        
        if type(hic) == str:
            if type(resolution) == int and genome != None:
                tmpHic = highResBinnedData.HiResHiC(genome, resolution)
                tmpHic.loadData(hic)
                hic = tmpHic
        
        # calculating reads at all distances
        nDistances = max(Counter(genome.chrmIdxBinCont).values())
        pixelsAtDistances = np.zeros(nDistances)
        readsAtDistances = np.zeros(nDistances)
        for chr1, chr2 in hic.data:
            if chr1 == chr2:
                data = hic.data[(chr1, chr2)].getData()
                
                for i in range(0,data.shape[0]):
                    for j in range(i,data.shape[1]):
                        distance = j-i   
                        readsAtDistances[distance] += data[i,j]                 
                        pixelsAtDistances[distance] += 1
        
        # smoothing distance vector
        smoothedPixelsAtDistances = np.zeros(len(pixelsAtDistances))
        smoothedReadsAtDistances = np.zeros(len(readsAtDistances))
        for i in range(0,len(readsAtDistances)):
            newReads = readsAtDistances[i]
            newPixels = pixelsAtDistances[i]
            windowSize = 0
            canExtend = True
            while newReads < 400 and canExtend: # smoothing
                windowSize += 1
                canExtend = False
                if i-windowSize >= 0:
                    newReads += readsAtDistances[i-windowSize]
                    newPixels += pixelsAtDistances[i-windowSize]
                    canExtend = True
                if i+windowSize < len(readsAtDistances):
                    newReads += readsAtDistances[i+windowSize]
                    newPixels += pixelsAtDistances[i+windowSize]
                    canExtend = True
            smoothedReadsAtDistances[i] = newReads
            smoothedPixelsAtDistances[i] = newPixels
        
        OE = smoothedReadsAtDistances/smoothedPixelsAtDistances
        
        def E(i,j):
            d = abs(i-j)
            return OE[d]
        
        # lower-left neighborhood
        def E_ll(M,i,j,w=1,p=0):
            sum1 = 0
            for a in range(i+1, i+w+1):
                for b in range(j-w, j):
                    sum1 += M[a,b]
            
            sum2 = 0
            for a in range(i+1, i+p+1):
                for b in range (j-p, j):
                    sum2 += M[a,b]
            
            sum3 = 0
            for a in range(i+1, i+w+1):
                for b in range(j-w, j):
                    sum3 += E(a,b)
            
            sum4 = 0
            for a in range(i+1, i+p+1):
                for b in range (j-p, j):
                    sum4 += E(a,b)
            
            return (sum1-sum2)/(sum3-sum4)*E(i,j)
        
        # sum of reads in lower-left neighborhood
        def ll_sum(M,i,j,w=1,p=0):
            sum1 = 0
            for a in range(i+1, i+w+1):
                for b in range(j-w, j):
                    sum1 += M[a,b]
            
            sum2 = 0
            for a in range(i+1, i+p+1):
                for b in range (j-p, j):
                    sum2 += M[a,b]
            
            return (sum1-sum2)
            
        # vertical neighborhood
        def E_v(M,i,j,w=1,p=0):
            sum1 = 0
            for a in range(i-w, i-p):
                for b in range(j-1, j+2):
                    sum1 += M[a,b]
            
            sum2 = 0
            for a in range(i+p+1, i+w+1):
                for b in range (j-1, j+2):
                    sum2 += M[a,b]
            
            sum3 = 0
            for a in range(i-w, i-p):
                for b in range(j-1, j+2):
                    sum3 += E(a,b)
            
            sum4 = 0
            for a in range(i+p+1, i+w+1):
                for b in range (j-1, j+2):
                    sum4 += E(a,b)
            
            return (sum1-sum2)/(sum3-sum4)*E(i,j)
        
        #horizontal neighborhood
        def E_h(M,i,j,w=1,p=0):
            sum1 = 0
            for a in range(i-1, i+2):
                for b in range(j-w, j-p):
                    sum1 += M[a,b]
            
            sum2 = 0
            for a in range(i-1, i+2):
                for b in range (j+p+1, j+w+1):
                    sum2 += M[a,b]
            
            sum3 = 0
            for a in range(i-1, i+2):
                for b in range(j-w, j-p):
                    sum3 += E(a,b)
            
            sum4 = 0
            for a in range(i-1, i+2):
                for b in range (j+p+1, j+w+1):
                    sum4 += E(a,b)
            
            return (sum1-sum2)/(sum3-sum4)*E(i,j)
        
        # donut neighborhood
        def E_d(M,i,j,w=1,p=0):
            topSum1 = 0
            for a in range(i-w, i+w+1):
                for b in range(j-w, j+w+1):
                    topSum1 += M[a,b]
        
            topSum2 = 0
            for a in range(i-p, i+p+1):
                for b in range (j-p, j+p+1):
                    topSum2 += M[a,b]
                    
            topSum3 = 0
            for a in range(i-w, i-p):
                topSum3 += M[a,j]
            
            topSum4 = 0
            for a in range(i+p+1, i+w+1):
                topSum4 += M[a,j]
                
            topSum5 = 0
            for b in range(j-w,j-p):
                topSum5 += M[i,b]
            
            topSum6 = 0
            for b in range(j+p+1,j+w+1):
                topSum6 += M[i,b]
            
            bottomSum1 = 0
            for a in range(i-w, i+w+1):
                for b in range(j-w, j+w+1):
                    bottomSum1 += E(a,b)
        
            bottomSum2 = 0
            for a in range(i-p, i+p+1):
                for b in range (j-p, j+p+1):
                    bottomSum2 += E(a,b)
                    
            bottomSum3 = 0
            for a in range(i-w, i-p):
                bottomSum3 += E(a,j)
            
            bottomSum4 = 0
            for a in range(i+p+1, i+w+1):
                bottomSum4 += E(a,j)
                
            bottomSum5 = 0
            for b in range(j-w,j-p):
                bottomSum5 += E(i,b)
            
            bottomSum6 = 0
            for b in range(j+p+1,j+w+1):
                bottomSum6 += E(i,b)
                
            return (topSum1-topSum2-topSum3-topSum4-topSum5-topSum6) / \
                (bottomSum1-bottomSum2-bottomSum3-bottomSum4-bottomSum5-bottomSum6) * \
                E(i,j)
        
        
        # determine p according to resolution
        if resolution > 25000:
            p = 1
            wInit = 3
        else:
            p = int(24999/resolution)
            wInit = round(25000/resolution) + 2
        
        # for every pixel calculate neighborhoods
        for chr1, chr2 in hic.data:
            if chr1 == chr2:
                M = hic.data[(chr1, chr2)].getData()
                
                for i in range(0,data.shape[0]):
                    for j in range(i,data.shape[1]):
                        # do not examine loci closer than p+3
                        if abs(i-j) > p+2:
                            w = wInit
                            while ll_sum(M,i,j,w=w,p=p) < 16 or w > 19:
                                w += 1
                            
                            # neighborhood expected values
                            ll = E_ll(M,i,j,w=w,p=p)
                            #p_ll = 1-poisson.cdf(obs,expectedInterReads)
                            h = E_h(M,i,j,w=w,p=p)
                            v = E_v(M,i,j,w=w,p=p)
                            d = E_d(M,i,j,w=w,p=p)
                                
        
        return cls([OE,smoothedReadsAtDistances,smoothedPixelsAtDistances])
        
        
    @classmethod
    def fromBinnedHiC(cls, fileName, hic, genome=None, resolution=None, intra=True, inter=True, total=False, fdrCutoff=0.05):
        
        if type(genome) == str:
            genome = kaic.genome.genomeTools.loadGenomeObject(genome)
        
        if type(hic) == str:
            if type(resolution) == int and genome != None:
                tmpHic = highResBinnedData.HiResHiC(genome, resolution)
                tmpHic.loadData(hic)
                hic = tmpHic
        
        
        
        container = getEmptyNetworkContainer(fileName)
        
        genomeContainer = container.get_node("/", 'genome')
        chromosomes = container.get_node(genomeContainer, 'chromosomes')
        
        network = container.get_node("/", 'network')
        nodes = container.get_node(network, 'nodes')
        edges = container.get_node(network, 'edges')
        
        
        chromosome = chromosomes.row
        mappableChrm = []
        for i in hic.genome.idx2label:
            chromosome['idNumber'] = i
            chromosome['name'] = hic.genome.idx2label[i]
            chromosome['length'] = hic.genome.chrmLens[i]
            chromosome.append()
            mappableChrm.append(0)
        chromosomes.flush()
        
        
        # nodes
        node = nodes.row
        nNodes = 0
        for i in range(0,len(hic.genome.chrmIdxBinCont)):
            chrm = hic.genome.chrmIdxBinCont[i]
            start = hic.genome.posBinCont[i]
            if i < len(hic.genome.chrmIdxBinCont)-1:
                end = hic.genome.posBinCont[i+1]-1
                if end < 0:
                    end = hic.genome.chrmLens[chrm]
            else:
                end = hic.genome.chrmLens[chrm]
            node['idNumber'] = i
            node['chromosomeId'] = chrm
            node['start'] = start
            node['end'] = end
            node['mappable'] = False
            node.append()
            nNodes += 1
        nodes.flush()
        
        # edges
        edge = edges.row
        edgeCounter = 0
        mappability = [False] * nNodes
        totalReads = 0
        intraReads = 0
        interReads = 0
        intraCounter = 0
        interCounter = 0
        for chr1, chr2 in hic.data:
            data = hic.data[(chr1, chr2)].getData()
            chr1StartBin = hic.genome.chrmStartsBinCont[chr1]
            chr2StartBin = hic.genome.chrmStartsBinCont[chr2]
            
            for i in range(0,data.shape[0]):
                iNode = i+chr1StartBin
                for j in range(i,data.shape[1]):
                    jNode = j+chr2StartBin
                    
                    if data[i,j] != 0:
                        edge['idNumber'] = edgeCounter
                        edge['source'] = iNode
                        edge['sink'] = jNode
                        edge['weight'] = data[i,j]
                        totalReads += data[i,j]
                        if chr1 == chr2:
                            intraReads += data[i,j]
                            intraCounter += 1
                        else:
                            interReads += data[i,j]
                            edge['inter'] = True
                            interCounter += 1
                        edgeCounter += 1
                        edge.append()
                        mappability[iNode] = True
                        mappability[jNode] = True
        edges.flush()
        
        print "Added %d edges" % edgeCounter
        print "Added %d intra edges" % intraCounter
        print "Added %d inter edges" % interCounter
        
        nMappable = 0
        for i in range(0,len(mappability)):
            chrm = nodes.cols.chromosomeId[i]
            
            nodes.cols.mappable[i] = mappability[i]
            if mappability[i] == True:
                nMappable += 1
                mappableChrm[chrm] += 1
        nodes.flush()
        
        
        # calculate possible intra- and inter-chromosomal contacts
        totalPossible = (nMappable**2 - nMappable) / 2
        intraPossible = 0
        interPossible = 0
        for i in range(0,len(mappableChrm)):
            for j in range(i, len(mappableChrm)):
                if i == j:
                    # diagonal is removed by default
                    intraPossible += mappableChrm[i]**2/2 - mappableChrm[i]/2
                else:
                    interPossible += mappableChrm[i] * mappableChrm[j]
        
        print "Possible contacts: %d" % totalPossible
        print "Possible intra-chromosomal contacts: %d" % intraPossible
        print "Possible inter-chromosomal contacts: %d" % interPossible
        
        expectedTotalReads = totalReads/totalPossible
        expectedIntraReads = intraReads/intraPossible
        expectedInterReads = interReads/interPossible
        pTotal = expectedTotalReads/totalPossible
        pIntra = expectedIntraReads/intraPossible
        pInter = expectedInterReads/interPossible
        
        print "Expected reads per contact: %.2f (p = %.2E)" % (expectedTotalReads, pTotal)
        print "Expected intra-chromosomal reads per contact: %.2f (p = %.2E)" % (expectedIntraReads, pIntra)
        print "Expected inter-chromosomal reads per contact: %.2f (p = %.2E)" % (expectedInterReads, pInter)
        
        # calculate probabilities
        psIntra = np.ones(intraPossible)
        psInter = np.ones(interPossible)
        intraCounter = 0
        interCounter = 0
        for i in range(0,edgeCounter):
            if i % int(edgeCounter/20) == 0:
                percent = int(i/int(edgeCounter/20))
                print "%d%% done" % (percent*5)
            obs = edges.cols.weight[i]

            if edges.cols.inter[i] == True:
                psInter[interCounter] = 1-poisson.cdf(obs,expectedInterReads)
                interCounter += 1
            else:
                psIntra[intraCounter] = 1-poisson.cdf(obs,expectedIntraReads)
                intraCounter += 1

        
        # correct for multiple testing
        fdrInter = multipletests(psInter,method='fdr_bh')[1]
        fdrIntra = multipletests(psIntra,method='fdr_bh')[1]
        interCounter = 0
        intraCounter = 0
        for i in range(0,edgeCounter):
            if edges.cols.inter[i] == True:
                edges.cols.p[i] = psInter[interCounter]
                edges.cols.fdr[i] = fdrInter[interCounter]
                interCounter += 1
            else:
                edges.cols.p[i] = psIntra[intraCounter]
                edges.cols.fdr[i] = fdrIntra[intraCounter]
                intraCounter += 1
        
        
        
        return cls(container)
    


def getEmptyNetworkContainer(fileName):
        h5file = t.open_file(fileName, mode = "w", title = "HiCNetwork")
        
        # network
        network = h5file.create_group("/", 'network', 'Hi-C network')
        h5file.create_table(network, 'nodes', Locus, 'Network nodes')
        h5file.create_table(network, 'edges', Edge, 'Network edges')
            
        # genome
        genome = h5file.create_group("/", 'genome', 'Genome information')
        h5file.create_table(genome, 'chromosomes', Chromosome, 'Chromosome information')
        
        return h5file
            
