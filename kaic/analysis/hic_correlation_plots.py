'''
Created on Jun 1, 2015

@author: kkruse1
'''

import pandas as p
from kaic.tools.hic import load_mirny_binned_hic
from rpy2.robjects import pandas2ri as p2r
from rpy2.robjects.packages import importr
import numpy as np
from scipy.stats.stats import pearsonr
from kaic.plotting.plot_genomic_data import open_graphics_file, close_graphics_file

def pairwise_data_frame(hic1, hic2, genome, resolution, include_zeros=False, chromosome1=None, chromosome2=None):
    #genome = loadGenomeObject(genome)
    hic1 = load_mirny_binned_hic(hic1, genome, resolution)
    hic2 = load_mirny_binned_hic(hic2, genome, resolution)
    
    chromosome2 = chromosome1 if chromosome2 is None else chromosome2
    
    l = []
    for chr1, chr2 in hic1.data:
        
        if chromosome1 is not None and chr1 != chromosome1:
            continue
        if chromosome2 is not None and chr2 != chromosome2:
            continue
        
        
        data1 = hic1.data[(chr1, chr2)].getData()
        data2 = hic2.data[(chr1, chr2)].getData()
        
        if not data1.shape == data2.shape:
            raise ValueError("The Hi-C objects do not appear to be comparable (different shapes)")
        
        chr1StartBin = hic1.genome.chrmStartsBinCont[chr1]
        chr2StartBin = hic1.genome.chrmStartsBinCont[chr2]
        
        
        for i in range(0,data1.shape[0]):
            iNode = i+chr1StartBin
            pos1 = hic1.genome.posBinCont[iNode]
            for j in range(i,data1.shape[1]):
                jNode = j+chr2StartBin
                pos2 = hic1.genome.posBinCont[jNode]
                
                if not include_zeros and data1[i,j] == 0 and data2[i,j] == 0:
                    continue
                
                l.append([chr1, pos1, chr2, pos2, data1[i,j], data2[i,j]])
        
        
    df = p.DataFrame(l, columns=["chr1", "pos1", "chr2", "pos2", "val1", "val2"])
    
    return df


def correlation_data_frame(hic1, hic2, genome, resolution, include_zeros=False, order=None):
    #genome = loadGenomeObject(genome)
    hic1 = load_mirny_binned_hic(hic1, genome, resolution)
    hic2 = load_mirny_binned_hic(hic2, genome, resolution)
    
    if order is None:
        order = []
        for label in hic1.genome.chrmLabels:
            order.append(label)
    
    indexes = []
    for label in order:
        indexes.append(hic1.genome.label2idx[label])
    
    
    l = np.zeros((len(indexes),len(indexes)))
    df = p.DataFrame(l, index=order, columns=order)
    
    for chr1_i in range(0,len(indexes)):
        for chr2_j in range(chr1_i,len(indexes)):
            # swap if wrong order
            if chr1_i > chr2_j:
                tmp = chr1_i
                chr1_i = chr2_j
                chr2_j = tmp
            
            chr1 = indexes[chr1_i]
            chr2 = indexes[chr2_j]

            data1 = hic1.data[(chr1, chr2)].getData()
            data2 = hic2.data[(chr1, chr2)].getData()
            
            if not data1.shape == data2.shape:
                raise ValueError("The Hi-C objects do not appear to be comparable (different shapes)")
   
            l1 = []
            l2 = []
            for i in range(0,data1.shape[0]):
                for j in range(i,data1.shape[1]):
                    if not include_zeros and data1[i,j] == 0 and data2[i,j] == 0:
                        continue
                    
                    l1.append(data1[i,j])
                    l2.append(data2[i,j])
            
            df[order[chr1_i]][order[chr2_j]] = pearsonr(l1,l2)[0]
            df[order[chr2_j]][order[chr1_i]] = pearsonr(l1,l2)[0]
                
    return df


def plot_chromosome_correlation(hic1, hic2, genome, resolution, include_zeros=False, order=None, output=None, width=9,height=9):
    df = correlation_data_frame(hic1, hic2, genome, resolution, include_zeros, order)
    df_p = df.copy()
    
    # set upper to -1
    labels = df_p.columns
    for i in range(0,len(labels)):
        il = labels[i]
        for j in range(i+1,len(labels)):
            jl = labels[j]
            
            df_p[jl][il] = -1
            
    
    p2r.activate()
    corrplot = importr("corrplot")
    base = importr("base")
            
    
    if output:
        open_graphics_file(output,width,height)
    
    # plot
    m = base.as_matrix(df)
    m.colnames = m.rownames
    m_p = base.as_matrix(df_p)
    m_p.colnames = m_p.rownames
    
    corrplot.corrplot(m,method="circle",tl_pos="lt",p_mat=m_p, insig="p-value", sig_level=-1)
    
    if output:
        close_graphics_file()
    
    
    