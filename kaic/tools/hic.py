'''
Created on May 13, 2015

@author: kkruse1
'''

#import numpy as np

def getChromosomeMatrix(hic,genome,chrm1,chrm2=None):
    if type(chrm1) == str:
        chrm1 = genome.label2idx(chrm1)
    if type(chrm2) == str:
        chrm2 = genome.label2idx(chrm2)
    if not chrm2:
        chrm2 = chrm1
    
    if chrm1 > chrm2:
        tmp = chrm1
        chrm1 = chrm2
        chrm2 = tmp
    
    #rows = genome.chrmEndsBinCont[chrm1]-genome.chrmStartsBinCont[chrm1]
    #cols = genome.chrmEndsBinCont[chrm2]-genome.chrmStartsBinCont[chrm2]
    #M = np.zeros((rows,cols),float)
    data = hic.data[(chrm1,chrm2)].getData()
    
    return data