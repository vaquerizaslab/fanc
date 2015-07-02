'''
Created on Jun 29, 2015

@author: kkruse1
'''

import numpy as np
from kaic.data.genomic import Chromosome, Genome

class TestChromosome:
    
    @classmethod
    def setup_method(self, method):
        self.chromosome = Chromosome(name='chr1', length=10000, sequence='agcgctgctgaagcttcgatcgtaagcttc')
        
    def test_attributes(self):
        assert self.chromosome.name == 'chr1'
        assert self.chromosome.length == 10000
        assert len(self.chromosome) == 10000
        assert self.chromosome.sequence == 'agcgctgctgaagcttcgatcgtaagcttc'
        
    def test_re_sites(self):
        res = self.chromosome.get_restriction_sites('HindIII')
        assert len(res) == 2
        assert np.array_equal(res, [12,25])


class TestGenome:
    
    @classmethod
    def setup_method(self, method):
        chr1 = Chromosome(name='chr1', length=10000, sequence='agcgctgctgaagcttcgatcgtaagcttc')
        chr2 = Chromosome(name='chr2', length=5000, sequence='gcgctgctgaagcttcgatcgtaagcttc')
        self.genome = Genome(chromosomes=[chr1,chr2])
        
    def test_iter(self):
        i = 0
        for chromosome in self.genome:
            if i == 0:
                assert chromosome.name == 'chr1'
                assert chromosome.length == 10000
                assert chromosome.sequence == 'agcgctgctgaagcttcgatcgtaagcttc'
            if i == 1:
                assert chromosome.name == 'chr2'
                assert chromosome.length == 5000
                assert chromosome.sequence == 'gcgctgctgaagcttcgatcgtaagcttc'
            i+=1
            
    def test_node_list(self):
        nl = self.genome.get_node_list('HindIII')
        
        assert len(nl) == 6
        for i in range(0,len(nl)):
            node = nl[i]
            if i == 0:
                assert node['chromosome'] == 'chr1'
                assert node['start'] == 1
                assert node['end'] == 12
            if i == 5:
                assert node['chromosome'] == 'chr2'
                assert node['start'] == 24
                assert node['end'] == 5000
            i += 1
            
            
        nl = self.genome.get_node_list(4000)
        
        assert len(nl) == 5
        for i in range(0,len(nl)):
            node = nl[i]
            if i == 0:
                assert node['chromosome'] == 'chr1'
                assert node['start'] == 1
                assert node['end'] == 4000
            if i == 5:
                assert node['chromosome'] == 'chr2'
                assert node['start'] == 4001
                assert node['end'] == 5000
            i += 1