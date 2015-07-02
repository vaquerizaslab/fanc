'''
Created on Jun 10, 2015

@author: kkruse1
'''

from kaic.genome.genomeTools import loadGenomeObject
from kaic.data.genomic import BedImproved
from random import randint

def random_bed(genome, n, size=100, keep_chromosome_distribution=True, name=None):
    genome = loadGenomeObject(genome)
    
    
    new = []
    if isinstance(n, BedImproved):
        original_bed = n
        
        for i in range(0,len(original_bed)):
            row = original_bed[i]
            chrom = row[:,"chrom"][3:]
            start = row[:,"start"]
            end = row[:,"end"]
            name = row[:,"name"]
            score= row[:,"score"]
            strand = row[:,"strand"]
            size = end-start
            
            chrom_ix = genome.label2idx[chrom]
            
            new_start = randint(1,genome.chrmLens[chrom_ix]-size)
            new_end = new_start + size
            
            
            new.append(['chr' + chrom, new_start, new_end, name, score, strand])
    elif type(n) is int:
        for i in range(0,n):
            chrom_ix = randint(0,len(genome.chrmLabels))
            chrom = genome.chrmLabels[chrom_ix]
            
            new_start = randint(1,genome.chrmLens[chrom_ix]-size)
            new_end = new_start + size
        
            new.append([chrom,new_start,new_end])
    
    bed = BedImproved(col_names=["chrom","start","end","name","score","strand"], col_types=[str,int,int,str,float,str], data=new, name=name)
    
    return bed