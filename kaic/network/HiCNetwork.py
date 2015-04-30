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
from bisect import bisect_right

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
    def liebermannBinnedHiC(cls, hic, genome=None, resolution=None, fdrIntra=None, fdrInter=None, findInterPeaks=True):
        peaks = get_peaks_from_hic(hic,genome,resolution,fdrIntra,fdrInter,findInterPeaks)
#         print "Merging peaks"
#         mergedPeaks = {}
#         while len(peaks) > 0: 
#             # find maximum peak
#             max_peak_height = 0
#             max_peak = None
#             for peak in peaks:
#                 if max_peak_height < peaks[peak]['obs']:
#                     max_peak_height = peaks[peak]['obs']
#                     max_peak = peak
#             
#             mergedPeak = { 'obs': max_peak_height,
#                            'list': [max_peak],
#                            'centroid': max_peak,
#                            'radius': 0
#                           }
#             
#             stepSize = 20000/resolution + 1 # +1 because neighboring bins have dist 0
#             foundClusterMember = True
#             while len(peaks) > 0 and foundClusterMember:
#                 foundClusterMember = False
#                 maxDist = mergedPeak['radius'] + stepSize
#                 for peak in peaks.keys():
#                     d = np.sqrt((mergedPeak['centroid'][0]-peak[0])**2+(mergedPeak['centroid'][1]-peak[1])**2)
#                     if d < maxDist:
#                         print "found peak to merge ", mergedPeak['list']
#                         foundClusterMember = True
#                         mergedPeak['list'].append(peak)
#                         
#                         #calculate new centroid
#                         x = 0
#                         y = 0
#                         r = 0
#                         for mpeak in mergedPeak['list']:
#                             x += mpeak[0]
#                             y += mpeak[1]
#                             r = max(r,np.sqrt((mergedPeak['centroid'][0]-mpeak[0])**2+(mergedPeak['centroid'][1]-mpeak[1])**2))
#                         mergedPeak['centroid'] = (x/len(mergedPeak['list']), y/len(mergedPeak['list']))
#                         mergedPeak['radius'] = r
#                         
#                         del peaks[peak]
#                         
#                         break
#             
#             print "Adding new merged peak"
#             mergedPeaks[mergedPeak['centroid']] = mergedPeak
#             print mergedPeak['list']
            

        return cls(peaks)
        
        # lambda-chunking
        
        
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
            



def get_peaks_from_hic(hic, genome=None, resolution=None, fdrIntra=None, fdrInter=None, findInterPeaks=True, output=None):
    if type(genome) == str:
        genome = kaic.genome.genomeTools.loadGenomeObject(genome)
    
    if type(hic) == str:
        if type(resolution) == int and genome != None:
            tmpHic = highResBinnedData.HiResHiC(genome, resolution)
            tmpHic.loadData(hic)
            hic = tmpHic
    
    print "Checking bin mappability"
    mappableChrm = {}
    for i in hic.genome.idx2label:
        mappableChrm[i] = []
        
    interObserved = 0
    for i in range(0,len(hic.genome.chrmIdxBinCont)):
        chrm = hic.genome.chrmIdxBinCont[i]
        mappableChrm[chrm].append(False)
    for chr1, chr2 in hic.data:
        M = hic.data[(chr1, chr2)].getData()
        
        for i in range(0,M.shape[0]):
            for j in range(i,M.shape[1]):
                
                if M[i,j] > 0:
                    mappableChrm[chr1][i] = True
                    mappableChrm[chr2][j] = True
                    if chr1 != chr2:
                        interObserved += M[i,j]
    
    print "Calculating possible contacts"
    interPossible = 0
    for i in range(0,len(mappableChrm)):
        i_mappable = sum(mappableChrm[i])
        for j in range(i, len(mappableChrm)):
            j_mappable = sum(mappableChrm[j])
            if i != j:
                interPossible += i_mappable * j_mappable
    
    interExpected = interPossible/interObserved
    
    print "Creating 1D model"
    oeDiff = {}
    # calculating reads at all distances
    nDistances = max(Counter(genome.chrmIdxBinCont).values())
    pixelsAtDistances = np.zeros(nDistances)
    readsAtDistances = np.zeros(nDistances)
    for chr1, chr2 in hic.data:
        if chr1 == chr2:
            oeDiff[chr1] = []
            
            M = hic.data[(chr1, chr2)].getData()
            
            for i in range(0,M.shape[0]):
                for j in range(i,M.shape[1]):
                    # only count mappable bins
                    if mappableChrm[chr1][i] and mappableChrm[chr2][j]:
                        distance = j-i
                        readsAtDistances[distance] += M[i,j]
                        pixelsAtDistances[distance] += 1

                    
    print "Smoothing 1D model"
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
    print "1D model length is %d" %len(OE)
    
    def E(i,j,inter=False):
        if inter:
            return interExpected
        d = abs(i-j)
        return OE[d]
    
    # sum of reads in lower-left neighborhood
    def ll_sum(M,i,j,w=1,p=0):
        i_max, j_max = M.shape
        
        sum1 = 0
        for a in range(max(0,i+1), min(i_max,i+w+1)):
            for b in range(max(0,j-w), min(j_max,j)):
                sum1 += M[a,b]
        
        sum2 = 0
        for a in range(max(0,i+1), min(i_max,i+p+1)):
            for b in range (max(0,j-p), min(j_max,j)):
                sum2 += M[a,b]
        
        return (sum1-sum2)
    
    # lower-left neighborhood
    def E_ll(M,i,j,w=1,p=0,inter=False):
        i_max, j_max = M.shape
        
        sum1 = 0
        for a in range(max(0,i+1), min(i+w+1,i_max)):
            for b in range(max(0,j-w), min(j_max,j)):
                sum1 += M[a,b]
        
        sum2 = 0
        for a in range(max(0,i+1), min(i_max,i+p+1)):
            for b in range (max(0,j-p), min(j_max,j)):
                sum2 += M[a,b]
        
        sum3 = 0
        for a in range(max(0,i+1), min(i_max,i+w+1)):
            for b in range(max(0,j-w), min(j_max,j)):
                sum3 += E(a,b,inter)
        
        sum4 = 0
        for a in range(max(0,i+1), min(i_max,i+p+1)):
            for b in range (max(0,j-p), min(j_max,j)):
                sum4 += E(a,b,inter)
        
        if sum3-sum4 == 0:
            print "sum1: %.2f, sum2: %.2f, sum3: %.2f, sum4: %.2f" % (float(sum1),float(sum2),float(sum3),float(sum4))
            print "start sum1/3a: %d, end sum1/3a: %d" % (max(0,i+1),min(i+w+1,i_max))
            print "start sum1/3b: %d, end sum1/3b: %d" % (max(0,j-w),min(j_max,j))
            print "start sum2/4a: %d, end sum2/4a: %d" % (max(0,i+1),min(i_max,i+p+1))
            print "start sum2/4b: %d, end sum2/4b: %d" % (max(0,j-p),min(j_max,j))
        return (sum1-sum2)/(sum3-sum4)*E(i,j,inter)
                    
    # vertical neighborhood
    def E_v(M,i,j,w=1,p=0,inter=False):
        i_max, j_max = M.shape
        
        sum1 = 0
        for a in range(max(0,i-w), min(i_max,i+w)):
            for b in range(max(0,j-1), min(j_max,j+2)):
                sum1 += M[a,b]
        
        sum2 = 0
        for a in range(max(0,i-p), min(i_max,i+p+1)):
            for b in range (max(0,j-1), min(j_max,j+2)):
                sum2 += M[a,b]
        
        sum3 = 0
        for a in range(max(0,i-w), min(i_max,i+w)):
            for b in range(max(0,j-1), min(j_max,j+2)):
                sum3 += E(a,b,inter)
        
        sum4 = 0
        for a in range(max(0,i-p), min(i_max,i+p+1)):
            for b in range (max(0,j-1), min(j_max,j+2)):
                sum4 += E(a,b,inter)
        
        return (sum1-sum2)/(sum3-sum4)*E(i,j,inter)
    
    #horizontal neighborhood
    def E_h(M,i,j,w=1,p=0,inter=False):
        i_max, j_max = M.shape
        
        sum1 = 0
        for a in range(max(0,i-1), min(i_max,i+2)):
            for b in range(max(0,j-w), min(j_max,j+w)):
                sum1 += M[a,b]
        
        sum2 = 0
        for a in range(max(0,i-1), min(i_max,i+2)):
            for b in range (max(0,j-p), min(j_max,j+p+1)):
                sum2 += M[a,b]
        
        sum3 = 0
        for a in range(max(0,i-1), min(i_max,i+2)):
            for b in range(max(0,j-w), min(j_max,j+w)):
                sum3 += E(a,b,inter)
        
        sum4 = 0
        for a in range(max(0,i-1), min(i_max,i+2)):
            for b in range (max(0,j-p), min(j_max,j+p+1)):
                sum4 += E(a,b,inter)
        
        return (sum1-sum2)/(sum3-sum4)*E(i,j,inter)
    
    # donut neighborhood
    def E_d(M,i,j,w=1,p=0,inter=False):
        i_max, j_max = M.shape
        
        topSum1 = 0
        for a in range(max(0,i-w), min(i_max,i+w+1)):
            for b in range(max(0,j-w), min(j_max,j+w+1)):
                topSum1 += M[a,b]
    
        topSum2 = 0
        for a in range(max(0,i-p), min(i_max,i+p+1)):
            for b in range (max(0,j-p), min(j_max,j+p+1)):
                topSum2 += M[a,b]
                
        topSum3 = 0
        for a in range(max(0,i-w), min(i_max,i-p)):
            topSum3 += M[a,j]
        
        topSum4 = 0
        for a in range(max(0,i+p+1), min(i_max,i+w+1)):
            topSum4 += M[a,j]
            
        topSum5 = 0
        for b in range(max(0,j-w),min(j_max,j-p)):
            topSum5 += M[i,b]
        
        topSum6 = 0
        for b in range(max(0,j+p+1),min(j_max,j+w+1)):
            topSum6 += M[i,b]
        
        bottomSum1 = 0
        for a in range(max(0,i-w), min(i_max,i+w+1)):
            for b in range(max(0,j-w), min(j_max,j+w+1)):
                bottomSum1 += E(a,b,inter)
    
        bottomSum2 = 0
        for a in range(max(0,i-p), min(i_max,i+p+1)):
            for b in range (max(0,j-p), min(j_max,j+p+1)):
                bottomSum2 += E(a,b,inter)
                
        bottomSum3 = 0
        for a in range(max(0,i-w), min(i_max,i-p)):
            bottomSum3 += E(a,j,inter)
        
        bottomSum4 = 0
        for a in range(max(0,i+p+1), min(i_max,i+w+1)):
            bottomSum4 += E(a,j,inter)
            
        bottomSum5 = 0
        for b in range(max(0,j-w),min(j_max,j-p)):
            bottomSum5 += E(i,b,inter)
        
        bottomSum6 = 0
        for b in range(max(0,j+p+1),min(j_max,j+w+1)):
            bottomSum6 += E(i,b,inter)
            
        return (topSum1-topSum2-topSum3-topSum4-topSum5-topSum6) / \
            (bottomSum1-bottomSum2-bottomSum3-bottomSum4-bottomSum5-bottomSum6) * \
            E(i,j,inter)
    
    
    # determine p according to resolution
    if resolution > 25000:
        p = 1
        wInit = 3
    else:
        p = int(24999/resolution)
        wInit = int(round(25000/resolution) + 2)
    C = hic.biases
    
    print "Initial values:\np=%d, w=%d" % (p,wInit)
    
    
    def calculate_neighborhood(E_local=E_d):
        print "Calculating neighborhood"
        
        ij = []
        exp = []
        obs = []
        ij_inter = []
        exp_inter = []
        obs_inter = []
        for chr1, chr2 in hic.data:
            if chr1 == chr2:
                print "Checking chr %d (%s)" % (chr1, genome.idx2label[chr1])
                
                M = hic.data[(chr1, chr2)].getData()
                
                # to convert chromosome into global (genome-wide) indices
                chr1StartBin = hic.genome.chrmStartsBinCont[chr1]
                chr2StartBin = hic.genome.chrmStartsBinCont[chr2]
                
                for i in range(0,M.shape[0]):
                    iNode = i+chr1StartBin
                    for j in range(i,M.shape[1]):
                        jNode = j+chr2StartBin
                                                
                        # do not examine loci closer than p+3
                        if abs(i-j) > p+2:
                            w = wInit
                            while ll_sum(M,i,j,w=w,p=p) < 16 and w < 20:
                                w += 1
                            
                            if w >= 20:
                                continue
                            # reproduce original value
                            o = int(M[i,j]*C[iNode]*C[jNode])
                            
                            # neighborhood expected values
                            ij.append([iNode,jNode])
                            e = E_local(M,i,j,w=w,p=p)
                            exp.append(e)
                            obs.append(o)
                            
                            oeDiff[chr1].append(o-e)
            elif findInterPeaks:
                print "Checking chr %d and chr %d " % (chr1, chr2)
                M = hic.data[(chr1, chr2)].getData()
                
                # to convert chromosome into global (genome-wide) indices
                chr1StartBin = hic.genome.chrmStartsBinCont[chr1]
                chr2StartBin = hic.genome.chrmStartsBinCont[chr2]
                
                for i in range(0,M.shape[0]):
                    iNode = i+chr1StartBin
                    for j in range(0,M.shape[1]):
                        jNode = j+chr2StartBin
                        
                        w = wInit
                        while ll_sum(M,i,j,w=w,p=p) < 16 and w < 20:
                            w += 1
                        
                        if w >= 20:
                            continue
                            
                        # reproduce original value
                        o = int(M[i,j]*C[iNode]*C[jNode])
                        
                        # neighborhood expected values
                        ij_inter.append([iNode,jNode])
                        e = E_local(M,i,j,w=w,p=p,inter=True)
                        exp_inter.append(e)
                        obs_inter.append(o)
                        
                
        return ij, exp, obs, ij_inter, exp_inter, obs_inter
    
    def get_chunks(ij,exp,obs):
        print "Chunking results"
        
        # lambda-chunking
        e_max = 1
        e_exp = 0
        chunks_max = []
        chunks = []
        max_exp = max(exp)
        while e_max < max_exp:
            chunks_max.append(e_max)
            
            chunks.append({'ij':[], 'exp':[], 'obs':[], 'max': e_max})
            e_exp += 1
            e_max = 2**(e_exp/3)
        chunks_max.append(e_max)
        chunks.append({'ij':[], 'exp':[], 'obs':[], 'max': e_max})
        
        
        
        for i in range(0,len(exp)):
            chunk = bisect_right(chunks_max, exp[i])
            chunks[chunk]['exp'].append(exp[i])
            chunks[chunk]['obs'].append(obs[i])
            chunks[chunk]['ij'].append(ij[i])
            #print ij[i]
        
        return chunks
    
    def get_peaks(chunks):
        print "Calculating peaks"
        
        peaks = {}
        total = 0
        for chunk in chunks:
            obs_counts = Counter(chunk['obs'])
            obs_sum = 0
            for count in obs_counts:
                obs_sum += obs_counts[count]
            x = []
            y = []
            for count in obs_counts:
                x.append(count)
                y.append(obs_counts[count]/obs_sum)
            
            # calculate x where FDR=fdr
            obs_integral_left = 0
            x_fdrs = {}
            for i in range(0,len(x)):
                obs_integral_left += y[i]
                obs_integral = 1-obs_integral_left
                poisson_integral = 1-poisson(chunk['max']).cdf(x[i])
                
                if poisson_integral > obs_integral or obs_integral == 0:
                    currentFdr = 1
                else:
                    currentFdr = poisson_integral/obs_integral
                
                x_fdrs[x[i]] = currentFdr
                #if poisson_integral < obs_integral*fdr:
                #    x_cutoff = x[i]
                #    break
            
            
            
            #if x_cutoff != None:
            for i in range(0,len(chunk['ij'])):
                total += 1
                if chunk['obs'][i] == 0 or (fdrIntra != None and x_fdrs[chunk['obs'][i]] >= fdrIntra):
                    continue
                peaks[(chunk['ij'][i][0], chunk['ij'][i][1])] = [chunk['obs'][i], chunk['exp'][i], x_fdrs[chunk['obs'][i]]]
                
        if isinstance(fdrIntra,float):
            print "\tFound %d/%d peaks with FDR<%.2f after lambda chunking correction" % (len(peaks),total,fdrIntra)
        else:
            print "\tFound %d/%d peaks with o > 0" % (len(peaks),total,fdrIntra)
                #if chunk['obs'][i] > x_cutoff:
                #    peaks[(chunk['ij'][i][0], chunk['ij'][i][1])] = [chunk['obs'][i], chunk['exp'][i]]
        return peaks
    
    def get_inter_peaks(ij, exp, obs):
        print "Calculating inter-chromosomal peaks"
        print "\t p-values..."
        
        l = len(ij)
        ps = [1] * l
        for i in range(0,l):
            if i % int(l/20) == 0:
                percent = int(i/int(l/20))
                print "\t%d%% done" % (percent*5)
                
            o = obs[i]
            e = exp[i]
            if o > 0:
                ps[i] = 1-poisson.cdf(o,e)
                if not (ps[i] >=0 and ps[i] < 1):
                    ps[i] = 1
        
        print "\t fdr values..."
        
        fdrs = multipletests(ps,method='fdr_bh')[1]
        
        print "\t calling significant neighborhood peaks..."
        
        peaks = {}
        for i in range(0,len(fdrs)):
            if obs[i] == 0 or (fdrInter != None and fdrs[i] > fdrInter):
                continue
            peaks[(ij[i][0],ij[i][1])] = [obs[i],exp[i],fdrs[i]]
        
        if isinstance(fdrInter,float):
            print "\tFound %d/%d peaks with FDR<%.2f after BH correction" % (len(peaks),l, fdrInter)
        else:
            print "\tFound %d/%d peaks with o > 0" % (len(peaks),l, fdrInter)
        return peaks
            
    
    print "Processing h neighborhood"
    h_ij, h_e, h_o, h_ij_inter, h_e_inter, h_o_inter = calculate_neighborhood(E_local=E_h)
    h_chunks = get_chunks(h_ij, h_e, h_o)
    h_peaks_intra = get_peaks(h_chunks)
    h_peaks_inter = get_inter_peaks(h_ij_inter, h_e_inter, h_o_inter)
    
    print "Processing ll neighborhood"
    ll_ij, ll_e, ll_o, ll_ij_inter, ll_e_inter, ll_o_inter = calculate_neighborhood(E_local=E_ll)
    ll_chunks = get_chunks(ll_ij, ll_e, ll_o)
    ll_peaks_intra = get_peaks(ll_chunks)
    ll_peaks_inter = get_inter_peaks(ll_ij_inter, ll_e_inter, ll_o_inter)
    
    print "Processing d neighborhood"
    d_ij, d_e, d_o, d_ij_inter, d_e_inter, d_o_inter = calculate_neighborhood(E_local=E_d)
    d_chunks = get_chunks(d_ij, d_e, d_o)
    d_peaks_intra = get_peaks(d_chunks)
    d_peaks_inter = get_inter_peaks(d_ij_inter, d_e_inter, d_o_inter)
    
    print "Processing v neighborhood"
    v_ij, v_e, v_o, v_ij_inter, v_e_inter, v_o_inter = calculate_neighborhood(E_local=E_v)
    v_chunks = get_chunks(v_ij, v_e, v_o)
    v_peaks_intra = get_peaks(v_chunks)
    v_peaks_inter = get_inter_peaks(v_ij_inter, v_e_inter, v_o_inter)
    
    peaks = {}
    print "Calling significant peaks"
    for peak in ll_peaks_intra:
        if (peak in d_peaks_intra and
            peak in h_peaks_intra and
            peak in v_peaks_intra):
                
            peaks[peak] = {'ll': ll_peaks_intra[peak][1],
                           'd': d_peaks_intra[peak][1],
                           'h': h_peaks_intra[peak][1],
                           'v': v_peaks_intra[peak][1],
                           'll_fdr': ll_peaks_intra[peak][2],
                           'd_fdr': d_peaks_intra[peak][2],
                           'h_fdr': h_peaks_intra[peak][2],
                           'v_fdr': v_peaks_intra[peak][2],
                           'obs': v_peaks_intra[peak][0],
                           'inter': 0,
                           'fdr': max(ll_peaks_intra[peak][2],
                                      h_peaks_intra[peak][2],
                                      v_peaks_intra[peak][2],
                                      d_peaks_intra[peak][2])}
            
    for peak in ll_peaks_inter:
        if (peak in d_peaks_inter and
            peak in h_peaks_inter and
            peak in v_peaks_inter):

            peaks[peak] = {'ll': ll_peaks_inter[peak][1],
                           'd': d_peaks_inter[peak][1],
                           'h': h_peaks_inter[peak][1],
                           'v': v_peaks_inter[peak][1],
                           'll_fdr': ll_peaks_inter[peak][2],
                           'd_fdr': d_peaks_inter[peak][2],
                           'h_fdr': h_peaks_inter[peak][2],
                           'v_fdr': v_peaks_inter[peak][2],
                           'obs': v_peaks_inter[peak][0],
                           'inter': 1,
                           'fdr': max(ll_peaks_inter[peak][2],
                                      h_peaks_inter[peak][2],
                                      v_peaks_inter[peak][2],
                                      d_peaks_inter[peak][2])}
            

#         print "Merging peaks"
#         mergedPeaks = {}
#         while len(peaks) > 0: 
#             # find maximum peak
#             max_peak_height = 0
#             max_peak = None
#             for peak in peaks:
#                 if max_peak_height < peaks[peak]['obs']:
#                     max_peak_height = peaks[peak]['obs']
#                     max_peak = peak
#             
#             mergedPeak = { 'obs': max_peak_height,
#                            'list': [max_peak],
#                            'centroid': max_peak,
#                            'radius': 0
#                           }
#             
#             stepSize = 20000/resolution + 1 # +1 because neighboring bins have dist 0
#             foundClusterMember = True
#             while len(peaks) > 0 and foundClusterMember:
#                 foundClusterMember = False
#                 maxDist = mergedPeak['radius'] + stepSize
#                 for peak in peaks.keys():
#                     d = np.sqrt((mergedPeak['centroid'][0]-peak[0])**2+(mergedPeak['centroid'][1]-peak[1])**2)
#                     if d < maxDist:
#                         print "found peak to merge ", mergedPeak['list']
#                         foundClusterMember = True
#                         mergedPeak['list'].append(peak)
#                         
#                         #calculate new centroid
#                         x = 0
#                         y = 0
#                         r = 0
#                         for mpeak in mergedPeak['list']:
#                             x += mpeak[0]
#                             y += mpeak[1]
#                             r = max(r,np.sqrt((mergedPeak['centroid'][0]-mpeak[0])**2+(mergedPeak['centroid'][1]-mpeak[1])**2))
#                         mergedPeak['centroid'] = (x/len(mergedPeak['list']), y/len(mergedPeak['list']))
#                         mergedPeak['radius'] = r
#                         
#                         del peaks[peak]
#                         
#                         break
#             
#             print "Adding new merged peak"
#             mergedPeaks[mergedPeak['centroid']] = mergedPeak
#             print mergedPeak['list']
        
    
    if output != None:
        with open(output, 'w') as out:
            out.write("i\tj\tinter\t\tobs\tfdr\td\tll\tv\th\td_fdr\tll_fdr\tv_fdr\th_fdr\n")
            
            for peak in peaks:
                i = peak[0]
                j = peak[1]
                inter = peaks[peak]['inter']
                o = int(peaks[peak]['obs'])
                fdr = float(peaks[peak]['fdr'])
                d = float(peaks[peak]['d'])
                ll = float(peaks[peak]['ll'])
                v = float(peaks[peak]['v'])
                h = float(peaks[peak]['h'])
                d_fdr = float(peaks[peak]['d_fdr'])
                ll_fdr = float(peaks[peak]['ll_fdr'])
                h_fdr = float(peaks[peak]['h_fdr'])
                v_fdr = float(peaks[peak]['v_fdr'])
                line = "%d\t%d\t%d\t%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n" % (i,j,inter,o,fdr,d,ll,v,h,d_fdr,ll_fdr,h_fdr,v_fdr)
                out.write(line)
        
    return peaks
        
    return peaks