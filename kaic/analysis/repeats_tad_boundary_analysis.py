#!/usr/bin/env python


'''
Created on Jun 5, 2015

@author: kkruse1
'''

import argparse
import kaic.data.genomic as gd
import pandas

def splitList(thisList):
    return thisList.split(",");

if __name__ == '__main__':
    parser = argparse.ArgumentParser();
    
    parser.add_argument(
        'input1',
        help='''Bed file 1'''
    )
    
    parser.add_argument(
        'input2',
        help='''Bed file 2'''
    )

    
    parser.add_argument(
        '-t', '--tad-file', dest='tad_file',
        help='''TADs BED file''',
        required=True
    )
    
    parser.add_argument(
        '-w', '--window', dest='window',
        type=int,
        help='''Window around the repeat to look for TAD boundaries''',
        required=True
    )

    args = parser.parse_args()
    
    print "TAD file: %s" % args.tad_file 
    tads = gd.Bed(args.tad_file)
    
    bed1 = gd.Bed(args.input1)
    bed2 = gd.Bed(args.input2)
    
    
    bed1_single = []
    bed1_double = []
    bed2_single = []
    bed2_double = []
    
    bed1_df = bed1.as_data_frame()
    bed2_df = bed2.as_data_frame()
    
    for i in range(0,len(bed1_df)):
        start = bed1_df["start"][i]
        end = bed1_df["end"][i]
        chrom = bed1_df["chrom"][i]
        
        end_tads = tads.as_data_frame(chrom, end=[start-args.window,start+args.window])
        start_tads = tads.as_data_frame(chrom, start=[end-args.window,end+args.window])
        
        if len(end_tads) > 0 or len(start_tads) > 0:
            bed1_single.append(True)
        else:
            bed1_single.append(False)
        
        if len(end_tads) > 0 and len(start_tads) > 0:
            bed1_double.append(True)
        else:
            bed1_double.append(False)
        
        
    for i in range(0,len(bed2_df)):
        start = bed2_df["start"][i]
        end = bed2_df["end"][i]
        chrom = bed2_df["chrom"][i]
        
        
        end_tads = tads.as_data_frame(chrom, end=[start-args.window,start+args.window])
        start_tads = tads.as_data_frame(chrom, start=[end-args.window,end+args.window])
        
        if len(end_tads) > 0 or len(start_tads) > 0:
            bed2_single.append(True)
        else:
            bed2_single.append(False)
        
        if len(end_tads) > 0 and len(start_tads) > 0:
            # also check if the TADs are actually different
            n_pairs = 0
            for j in range(0,len(start_tads)):
                found = False
                for k in range(0,len(end_tads)):
                    if (start_tads["start"][j] == end_tads["start"][k] and
                        start_tads["end"][j] == end_tads["end"][k] and
                        start_tads["chrom"][j] == end_tads["chrom"][k]):
                        found = True
                if found:
                    n_pairs += 1
            
            if not n_pairs == len(start_tads):
                bed2_double.append(True)
            else:
                bed2_double.append(False)
        else:
            bed2_double.append(False)
        
        
    colnames = ["bed1", "bed2"]
    rownames = ["True", "False"]
    
    result_single = [[sum(bed1_single),sum(bed2_single)],[len(bed1_single)-sum(bed1_single),len(bed2_single)-sum(bed2_single)]]
    result_double = [[sum(bed1_double),sum(bed2_double)],[len(bed1_double)-sum(bed1_double),len(bed2_double)-sum(bed2_double)]]
    
    
    print "Single TAD boundary in vicinity:"
    print pandas.DataFrame(result_single, index=rownames, columns=colnames)
    
    print "Dual TAD boundary in vicinity:"
    print pandas.DataFrame(result_double, index=rownames, columns=colnames)
    
    
    
        
    