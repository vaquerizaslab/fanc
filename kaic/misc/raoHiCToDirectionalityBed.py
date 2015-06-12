#!/usr/bin/env python


'''
Created on Jun 5, 2015

@author: kkruse1
'''

import argparse
import kaic.data.genomic as gd


def splitList(thisList):
    return thisList.split(",");

if __name__ == '__main__':
    parser = argparse.ArgumentParser();
    
    parser.add_argument(
        'input',
        nargs='+',
        help='''Hi-C files (Bedpe)'''
    )
    
    parser.add_argument(
        '-res', '--resolution', dest='resolution',
        type=int,
        help='''Resolution of the Hi-C data set''',
        required=True
    )
    
    parser.add_argument(
        '-c', '--chromosomes', dest='chromosomes',
        type=splitList,
        help='''Comma-separated list of chromosome names''',
        required=True
    )
    
    parser.add_argument(
        '-w', '--window', dest='window',
        type=int,
        default=2000000,
        help='''Window size for directionality calculation'''
    )
   
    parser.add_argument(
        '-o', '--output-file', dest='output',
        help='''Output file''',
        required=True
    )

    args = parser.parse_args()
    
    with open(args.output, 'w') as o:
        for i in range(0,len(args.input)):
            hic = gd.Hic(args.input[i])
            chrm = args.chromosomes[i]
            
            di = hic.directionality(args.resolution, window_size=args.window)
        
            for j in range(0,len(di)):
                start = j*args.resolution
                end = start + args.resolution-1
                name = "di_%s_%d" %(chrm, j)
                score = di[j]
                
                o.write("%s\t%d\t%d\t%s\t%f\t+\n" % (chrm, start, end, name, score)) 
        
        
    
        
    