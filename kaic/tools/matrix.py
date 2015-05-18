'''
Created on May 11, 2015

@author: kkruse1
'''

import numpy as np
import matplotlib
import re

from matplotlib import pyplot as plt

def readMatrixFromFile(file_name, delim="\t"):
    with open(file_name, 'r') as f:
        M = [ map(float,line.split(delim)) for line in f ]
    return M



def writeMatrixToFile(M,file_name,delim="\t",row_names=None,col_names=None):
    # check if M is 2D matrix
    try:
        n_rows = len(M)
        n_cols = len(M[0])
    except IndexError:
        raise IndexError("Input must be 2D matrix")

    with open(file_name, 'w') as o:
        if col_names != None:
            for name in col_names:
                o.write(name + delim)
            o.write("\n")
        
        for i in range(0,n_rows):
            if row_names:
                o.write(row_names[i] + delim)
            
            for j in range(0,n_cols):
                o.write("%.6E%s" % (M[i,j],delim))
            o.write("\n")

    print "Done writing to file."

def removeSparseRows(M,cutoff=None):
    s = np.sum(M,0)
    
    if cutoff == None:
        cutoff = min(s)
    
    idxs = np.where(s <= cutoff)[0]
    A = np.delete(M,idxs,0)
    A = np.delete(A,idxs,1)
    
    return A, idxs
    
    
def restoreSparseRows(M,idxs,rows=None):
    idxsn = idxs.copy()
    for i in range(0,len(idxs)):
        idxsn[i] = idxs[i]-i
    
    
    A = np.insert(M,idxsn,0,axis=0)
    if len(M.shape) > 1:
        A = np.insert(A,idxsn,0,axis=1)

    return A

def compare(A,M):
    return sum(abs(M-A),0)

def is_symmetric(M, tol=1e-10):
    for i in range(0,M.shape[0]):
        for j in range(i,M.shape[1]):
            if abs(M[i,j]-M[j,i]) > tol:
                print "(%d,%d) %.6f != %.6f (%d,%d)" % (i,j,M[i,j],M[j,i],j,i)
                return False
    return True
            
            
def fromEdgeListFile(inFile, resolution, output=None):
    l = []
    maxLocus = 0
    with open(inFile,'r') as f:
        for line in f:
            line = line.rstrip()
            l1,l2,v = line.split("\t")
            l1 = int(l1)
            l2 = int(l2)
            l.append([l1,l2,v])
            if maxLocus < l1:
                maxLocus = l1
            if maxLocus < l2:
                maxLocus = l2
    
    print maxLocus
    names = range(0,maxLocus+resolution,resolution)
    locus2idx = {}
    for i in range(0,len(names)):
        locus2idx[names[i]] = i
    
    M = np.zeros((len(names),len(names)))
    for l1, l2, v in l:
        M[locus2idx[l1],locus2idx[l2]] = v
        M[locus2idx[l2],locus2idx[l1]] = v
    
    if output:
        with open(output,'w') as o:
            line = ""
            for name in names:
                line = str(name) + "\t"
                #o.write(name, "\t")
            line = re.sub("\t$", "", line)
            o.write(line + "\n")
            
            for i in range(0,M.shape[0]):
                line = str(names[i]) + "\t"
                #o.write(names[i], "\t")
                for j in range(0,M.shape[1]):
                    line += "%.6e\t" % M[i,j]
                    #o.write("%.6e\t" % M[i,j])
                line = re.sub("\t$", "", line)
                o.write(line + "\n")
    
    return M
    
    
def plot(M, absolute=False, colormap = None, vmin=-3, vmax=3, iStartIndex = None, iEndIndex = None, jStartIndex = None, jEndIndex = None, highlightPixels=None):
    
    if iStartIndex == None:
        iStartIndex = 0
    if iEndIndex == None:
        iEndIndex = M.shape[0]-1
    if jStartIndex == None:
        jStartIndex = 0
    if jEndIndex == None:
        jEndIndex = M.shape[1]-1
    
    hm = M[iStartIndex:iEndIndex+1, jStartIndex:jEndIndex+1]
    
    
    if absolute == False:
        ex = np.sum(hm)/(hm.shape[0]*hm.shape[1])
        print "Expected: ", ex
        hm = np.log2(hm/ex)

    
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
    myPlot = ax.imshow(hm, interpolation='none',aspect=1,vmin=vmin,vmax=vmax)
    if colormap == None:
        myPlot.set_cmap(cmap)
    else:
        myPlot.set_cmap(colormap)
    
    
    if highlightPixels != None:
        for pixel in highlightPixels:
            if pixel[1] < iStartIndex or pixel[1] > iEndIndex:
                continue
            if pixel[0] < jStartIndex or pixel[0] > jEndIndex:
                continue
            i = pixel[1]-iStartIndex
            j = pixel[0]-jStartIndex
            
            c = plt.Circle((i,j),radius=1,fill=False,color='r')
            ax.add_patch(c)
    plt.show()
            
            
            
            
            
            
            
            