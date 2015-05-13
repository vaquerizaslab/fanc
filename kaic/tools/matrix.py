'''
Created on May 11, 2015

@author: kkruse1
'''

import numpy as np

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
        for name in col_names:
            o.write(name + delim)
        o.write("\n")
        
        for i in n_rows:
            if row_names:
                o.write(row_names[i] + delim)
            
            for j in n_cols:
                o.write("%.6E%s" + (M[i,j],delim))

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