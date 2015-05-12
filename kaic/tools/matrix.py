'''
Created on May 11, 2015

@author: kkruse1
'''

def readMatrixFromFile(fileName, delim=","):
    with open(fileName, 'r') as f:
        M = [ map(float,line.split(delim)) for line in f ]
    return M




