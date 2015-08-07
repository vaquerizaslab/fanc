'''
Created on Jun 29, 2015

@author: kkruse1
'''

import numpy as np
from kaic.data.genomic import Chromosome, Genome, HicBasic, HicNode, HicEdge
import os.path
import random

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
        nl = self.genome.get_regions('HindIII')
        
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
            
            
        nl = self.genome.get_regions(4000)
        
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
            
class TestHicBasic:
    
    @classmethod
    def setup_method(self, method):
        self.hic = HicBasic()
        
        # add some nodes (120 to be exact)
        nodes = []
        for i in range(1,50000,1000):
            nodes.append(HicNode(chromosome="chr1",start=i,end=i+1000-1))
        for i in range(1,30000,1000):
            nodes.append(HicNode(chromosome="chr2",start=i,end=i+1000-1))
        for i in range(1,20000,500):
            nodes.append(HicNode(chromosome="chr3",start=i,end=i+1000-1))
        self.hic.add_nodes(nodes)
        
        # add some edges randomly
        edges = []
        for i in range(0,1000):
            ix1 = random.randint(0,len(nodes)-1)
            ix2 = random.randint(0,len(nodes)-1)
            edges.append(HicEdge(source=ix1,sink=ix2,weight=random.random()))
        self.hic.add_edges(edges)
    
    def test_initialize_xml(self):
        current_dir = os.path.dirname(os.path.realpath(__file__))
        
        # from XML
        hic1 = HicBasic(current_dir + "/test_genomic/hic.example.xml")
        nodes1 = hic1.nodes()
        edges1 = hic1.edges()
        assert len(nodes1) == 2
        assert len(edges1) == 1
    
    def test_initialize_empty(self):
        hic = HicBasic()
        nodes = hic.nodes()
        edges = hic.edges()
        assert len(nodes) == 0
        assert len(edges) == 0
    
    def test_save_and_load(self, tmpdir):
        current_dir = os.path.dirname(os.path.realpath(__file__))
        dest_file = str(tmpdir) + "/hic.h5" 
        
        # from XML
        hic1 = HicBasic(current_dir + "/test_genomic/hic.example.xml", file_name=dest_file)
        hic1.close()
        hic2 = HicBasic(dest_file)
        nodes2 = hic2.nodes()
        edges2 = hic2.edges()
        assert len(nodes2) == 2
        assert len(edges2) == 1
        
        # dedicated save function
        # ...
    
    
    def test_nodes(self):
        nodes = self.hic.nodes()
        assert len(nodes) == 120
    
    def test_edges(self):
        edges = self.hic.edges()
        assert len(edges) == 1000
        
    def test_from_mirny(self):
        pass
        
        
        