'''
Created on May 13, 2015

@author: kkruse1
'''

from hiclib import highResBinnedData
from kaic.genome.genomeTools import loadGenomeObject

def getChromosomeMatrix(hic,chrm1,genome=None,chrm2=None):
    if genome is None:
        genome = hic.genome
    if type(chrm1) == str:
        chrm1 = genome.label2idx[chrm1]
    if type(chrm2) == str:
        chrm2 = genome.label2idx[chrm2]
    if not chrm2:
        chrm2 = chrm1
    
    if chrm1 > chrm2:
        tmp = chrm1
        chrm1 = chrm2
        chrm2 = tmp
    
    #rows = genome.chrmEndsBinCont[chrm1]-genome.chrmStartsBinCont[chrm1]
    #cols = genome.chrmEndsBinCont[chrm2]-genome.chrmStartsBinCont[chrm2]
    #M = np.zeros((rows,cols),float)
    print chrm1
    print chrm2
    data = hic.data[(chrm1,chrm2)].getData()
    
    return data


def load_mirny_binned_hic(file_name, genome, resolution):
    
    if type(file_name) == highResBinnedData.HiResHiC:
        return file_name
    
    genome = loadGenomeObject(genome)
    
    hic = highResBinnedData.HiResHiC(genome, resolution)
    hic.loadData(file_name)
    
    return hic
    

