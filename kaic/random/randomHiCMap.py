'''
Created on Mar 12, 2015

@author: kkruse1
'''

import numpy as np
import argparse
import time
import matplotlib.pyplot as plt


def kth_diag_indices(a, k):
    rows, cols = np.diag_indices_from(a)
    if k < 0:
        return rows[:k], cols[-k:]
    elif k > 0:
        return rows[k:], cols[:-k]
    else:
        return rows, cols
    
def randomHiCMap(size=1000,upper=500,middle=100,lower=0,nTads=4,zeros=10,exponent=0.008):    
    M = np.random.randint(lower, middle, size=(size, size))
    
    step = size/nTads
    tadStarts = range(0,size,step)
    if tadStarts[len(tadStarts)-1] > size-step:
        tadStarts.pop()
    
    # create TADs
    for i in range(0,len(tadStarts)):
        tadStart = tadStarts[i]
        if i < len(tadStarts)-1:
            tadEnd = tadStarts[i+1]-1
        else:
            tadEnd = size-1
        
        tadSize = tadEnd-tadStart
        
        TAD = np.random.randint( middle, upper, size=(tadSize, tadSize))
        
        TAD[50,0:50-1] = 1000
        
        for tadBreakPoint in range(0,tadSize):

            if tadBreakPoint > 0:
                expSum = np.ones(tadBreakPoint-1)
                
                for j in range(0,len(expSum)):
                    expSum[j] = np.exp(-1*exponent*j)
                
                TAD[0:tadBreakPoint-1,tadBreakPoint] = TAD[0:tadBreakPoint-1,tadBreakPoint]*expSum[::-1]+lower

        M[tadStart:tadEnd,tadStart:tadEnd] = TAD

                
    
    
    
    # add zero bins if requested
    for i in range(0,zeros):
        zeroBin = np.random.randint(0,size-1);
        M[zeroBin,:] = 0
        M[:,zeroBin] = 0
        
    
    # make matrix symmetrical
    for i in range(0,size):
        for j in range(i+1,size):
            M[j,i] = M[i,j]
    
    return M

def main(args):
    print("Using the following settings")
    for arg, val in args.__dict__.iteritems():
        print arg, " = ", val
        
    time.sleep(2);
    
    hicMap = randomHiCMap(args.size,args.upper,args.middle,args.lower,args.tads,args.zeros,args.exponent)
    
    if args.plot == True:
        fig, ax = plt.subplots()
        hm = ax.imshow(hicMap, aspect=1)
        plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser();
    
    parser.add_argument(
        '-s', '--size', dest='size',
        default=1000,
        type=int,
        help='''Size of the heatmap (in bins)'''
    );
    
    parser.add_argument(
        '-u', '--upper-limit', dest='upper',
        default=500,
        type=int,
        help='''Upper limit of reads per contact'''
    );
    
    parser.add_argument(
        '-l', '--lower-limit', dest='lower',
        default=0,
        type=int,
        help='''Lower limit of reads per contact'''
    );
    
    parser.add_argument(
        '-m', '--middle-point', dest='middle',
        default=100,
        type=int,
        help='''Middle point of reads per contact inside and outside of a TAD/chromosome'''
    );
    
    parser.add_argument(
        '-t', '--tads', dest='tads',
        default=4,
        type=int,
        help='''Number of TADs to generate'''
    );
    
    parser.add_argument(
        '-z', '--zeros', dest='zeros',
        default=20,
        type=int,
        help='''Number of bins with zero sum to generate (unmappable regions)'''
    );
    
    parser.add_argument(
        '-e', '--exponent', dest='exponent',
        default=0.008,
        type=float,
        help='''Tweak exponential decay with this exponent for intra-chromosomal contacts/TADs'''
    );
    
    parser.add_argument(
        '-o', '--output', dest='output',
        default = '',
        help='''Output file'''
    );
    
    
    
    parser.add_argument(
        '-p', '--plot', dest='plot',
        action='store_true',
        help='''Plot with matplotlib after generating'''
    );
    

    main(parser.parse_args());