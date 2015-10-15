'''
Created on Jun 29, 2015

@author: kkruse1
'''

import numpy as np
from kaic.data.genomic import Chromosome, Genome, HicBasic, HicNode, HicEdge,\
    GenomicRegion, GenomicRegions, _get_overlap_map, _edge_overlap_split_rao
import os.path
import pytest
from kaic.construct.seq import Reads, FragmentMappedReadPairs
import logging

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
        regions = self.genome.get_regions('HindIII')
        
        assert len(regions) == 6
        for i in range(0,len(regions)):
            region = regions[i]
            if i == 0:
                assert region['chromosome'] == 'chr1'
                assert region['start'] == 1
                assert region['end'] == 12
            if i == 5:
                assert region['chromosome'] == 'chr2'
                assert region['start'] == 25
                assert region['end'] == 5000
            i += 1
            
            
        nl = self.genome.get_regions(4000)
        
        print nl
        
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
            
class TestGenomicRegions:
    @classmethod
    def setup_method(self, method):
        chromosomes = [
            {'name': 'chr1', 'end': 10000},
            {'name': 'chr2', 'end': 15000},
            {'name': 'chr3', 'end': 7000}
        ]
        
        regions = []
        for chromosome in chromosomes:
            for start in range(1,chromosome["end"]-1000, 1000):
                regions.append(GenomicRegion(start,start+999,chromosome=chromosome["name"]))
        self.regions = GenomicRegions(regions)
    
    def test_get_item(self):
        region = self.regions[0]
        assert isinstance(region, GenomicRegion)
        assert region.chromosome == 'chr1'
        assert region.start == 1
        assert region.end == 1000
        assert region.strand == None
        
    def test_from_string(self):
        region1 = GenomicRegion.from_string('chr1')
        assert region1.chromosome == 'chr1'
        assert region1.start is None
        assert region1.end is None
        assert region1.strand == None
        
        region2 = GenomicRegion.from_string('chr1:0')
        assert region2.chromosome == 'chr1'
        assert region2.start == 0
        assert region2.end is None
        assert region2.strand == None
        
        region3 = GenomicRegion.from_string('chr1:0-4956')
        assert region3.chromosome == 'chr1'
        assert region3.start == 0
        assert region3.end == 4956
        assert region3.strand == None
        
        region4 = GenomicRegion.from_string('chr1:0-4956:-')
        assert region4.chromosome == 'chr1'
        assert region4.start == 0
        assert region4.end == 4956
        assert region4.strand == -1
        
        region5 = GenomicRegion.from_string('chr1:0-4956:+1')
        assert region5.chromosome == 'chr1'
        assert region5.start == 0
        assert region5.end == 4956
        assert region5.strand == 1
        
        with pytest.raises(ValueError):
            # invalid start
            GenomicRegion.from_string('chr1:x-4956:-')
        with pytest.raises(ValueError):
            # too many fields
            GenomicRegion.from_string('chr1:0:4956:-')
        with pytest.raises(ValueError):
            # invalid strand
            GenomicRegion.from_string('chr1:0-4956:0')
    
        
class TestHicBasic:
    
    @classmethod
    def setup_method(self, method):
        self.dir = os.path.dirname(os.path.realpath(__file__))
        
        hic = HicBasic()
        
        # add some nodes (120 to be exact)
        nodes = []
        for i in range(1,5000,1000):
            nodes.append(HicNode(chromosome="chr1",start=i,end=i+1000-1))
        for i in range(1,3000,1000):
            nodes.append(HicNode(chromosome="chr2",start=i,end=i+1000-1))
        for i in range(1,2000,500):
            nodes.append(HicNode(chromosome="chr3",start=i,end=i+1000-1))
        hic.add_nodes(nodes)
        
        # add some edges with increasing weight for testing
        edges = []
        weight = 1
        for i in range(0,len(nodes)):
            for j in range(i,len(nodes)):
                edges.append(HicEdge(source=i,sink=j,weight=weight))
                weight += 1

        hic.add_edges(edges)
        
        self.hic = hic
    
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
        
#         from subprocess import check_output
#         print check_output(["h5dump", dest_file])
#         print class_id_dict
        
        hic2 = HicBasic(dest_file)
        nodes2 = hic2.nodes()
        edges2 = hic2.edges()
        assert len(nodes2) == 2
        assert len(edges2) == 1
        
        # dedicated save function
        # ...
    
    
    def test_nodes(self):
        nodes = self.hic.nodes()
        assert len(nodes) == 12
    
    def test_edges(self):
        edges = self.hic.edges()
        assert len(edges) == 78
    
    def test_get_node_x_by_region(self):
        region1 = GenomicRegion.from_string('chr1')
        nodes1 = self.hic._getitem_nodes(region1)
        assert len(nodes1) == 5
        
        region2 = GenomicRegion.from_string('chr2')
        nodes2 = self.hic._getitem_nodes(region2)
        assert len(nodes2) == 3
        
        region3 = GenomicRegion.from_string('chr3')
        nodes3 = self.hic._getitem_nodes(region3)
        assert len(nodes3) == 4
        
        region4 = GenomicRegion.from_string('chr1:3452-6000')
        nodes4 = self.hic._getitem_nodes(region4)
        assert len(nodes4) == 2
        
        region5 = GenomicRegion.from_string('chr1:1-51000')
        nodes5 = self.hic._getitem_nodes(region5)
        assert len(nodes5) == 5
        
    def test_getitem_nodes(self):
        # all
        node_ix1 = self.hic._getitem_nodes(slice(None,None,None), as_index=True)
        assert np.array_equal(node_ix1, [0,1,2,3,4,5,6,7,8,9,10,11])
        
        # smaller slice
        node_ix2 = self.hic._getitem_nodes(slice(4,10,1), as_index=True)
        assert np.array_equal(node_ix2, [4,5,6,7,8,9])
        
        # single ix
        node_ix3 = self.hic._getitem_nodes(1, as_index=True)
        assert node_ix3 == 1
        
        # single chromosome
        node_ix4 = self.hic._getitem_nodes('chr1', as_index=True)
        assert np.array_equal(node_ix4, [0,1,2,3,4])
        
        # HicNode
        node_ix5 = self.hic._getitem_nodes(HicNode(ix=1), as_index=True)
        assert node_ix5 == 1
        
        # list of items
        node_ix6 = self.hic._getitem_nodes(['chr1','chr3'], as_index=True)
        assert np.array_equal(node_ix6, [0,1,2,3,4,8,9,10,11])
        
        # nested list of items
        node_ix7 = self.hic._getitem_nodes(['chr1',['chr2','chr3']], as_index=True)
        assert np.array_equal(node_ix7, [0,1,2,3,4,5,6,7,8,9,10,11])
        
        # item repetition
        node_ix8 = self.hic._getitem_nodes(['chr3','chr3'], as_index=True)
        assert np.array_equal(node_ix8, [8,9,10,11,8,9,10,11])

    def test_get_matrix(self):
        # whole matrix
        m = self.hic[:,:]
        # spot checks
        assert np.array_equal(m.shape, [12,12])
        assert m[0,0] == 1
        assert m[11,11] == 78
        assert m[11,0] == 12
        assert m[1,10] == 22
        # symmetry check
        for i in range(0,12):
            for j in range(i,12):
                assert m[i,j] == m[j,i]
                
        # only upper half
        m = self.hic[:6,:]
        # spot checks
        assert np.array_equal(m.shape, [6,12])
        assert m[0,0] == 1
        assert m[5,11] == 57
        assert m[5,0] == 6
        assert m[1,10] == 22
        
        # only lower half
        m = self.hic[6:,:]
        # spot checks
        assert np.array_equal(m.shape, [6,12])
        assert m[0,0] == 7
        assert m[5,11] == 78
        assert m[5,0] == 12
        assert m[1,10] == 67
        
        # only left half
        m = self.hic[:,:6]
        # spot checks
        assert np.array_equal(m.shape, [12,6])
        assert m[0,0] == 1
        assert m[11,5] == 57
        assert m[11,0] == 12
        assert m[1,4] == 16
        
        # only right half
        m = self.hic[:,6:]
        # spot checks
        assert np.array_equal(m.shape, [12,6])
        assert m[0,0] == 7
        assert m[11,5] == 78
        assert m[11,0] == 63
        assert m[1,4] == 22
        
        # top-left chunk
        m = self.hic[:6,:6]
        # spot checks
        assert np.array_equal(m.shape, [6,6])
        assert m[0,0] == 1
        assert m[5,5] == 51
        assert m[5,0] == 6 
        assert m[1,4] == 16
        
        # bottom_right chunk
        m = self.hic[6:,6:]
        # spot checks
        assert np.array_equal(m.shape, [6,6])
        assert m[0,0] == 58
        assert m[5,5] == 78
        assert m[5,0] == 63
        assert m[1,4] == 67
        
        # central chunk
        m = self.hic[3:9,3:9]
        # spot checks
        assert np.array_equal(m.shape, [6,6])
        assert m[0,0] == 34
        assert m[5,5] == 69
        assert m[5,0] == 39
        assert m[1,4] == 46
        
        # disjunct pieces
        m = self.hic[[1,9],[4,6]]
        # spot checks
        assert np.array_equal(m.shape, [2,2])
        assert m[0,0] == 16
        assert m[1,1] == 61
        assert m[1,0] == 48
        assert m[0,1] == 18
        
        # single row
        m = self.hic[1,0:3]
        assert np.array_equal(m, [2,13,14])
        
        # single row but only one value
        m = self.hic[1,1:2]
        assert np.array_equal(m, [13])
        
        # single col
        m = self.hic[0:3,1]
        assert np.array_equal(m, [2,13,14])
        
        # single col but only one value
        m = self.hic[1:2,1]
        assert np.array_equal(m, [13])
        
        # single value
        m = self.hic[1,1]
        assert m == 13
        
        # empty array
        m = self.hic[1:1,2:2]
        assert np.array_equal(m.shape, [0,0])
    
    def test_set_matrix(self):
        
        hic = HicBasic(self.hic)
        
        # whole matrix
        old = hic[:,:]
        # set diagonal to zero
        for i in range(0,old.shape[0]):
            old[i,i] = 0
        hic[:,:] = old
        m = hic[:,:]
        assert np.array_equal(m.shape, old.shape)
        for i in range(0,m.shape[0]):
            for j in range(0, m.shape[1]):
                if i == j:
                    assert m[i,j] == 0
                else:
                    assert m[i,j] == old[i,j]
        
        # central matrix
        hic = HicBasic(self.hic)
        old = hic[2:8,2:10]
        # set border elements to zero
        # set checkerboard pattern
        for i in range(0,old.shape[0]):
            for j in range(0,old.shape[1]):
                if i == 0 or j == 0:
                    old[i,j] = 0
                elif i % 2 == 0 and j % 2 == 1:
                    old[i,j] = 0
                elif i % 2 == 1 and j % 2 == 0:
                    old[i,j] = 0
        hic[2:8,2:10] = old
        m = hic[2:8,2:10]
        
        assert np.array_equal(m.shape, old.shape)
        for i in range(0,m.shape[0]):
            for j in range(0, m.shape[1]):
                if i == 0 or j == 0:
                    assert m[i,j] == 0
                elif i % 2 == 0 and j % 2 == 1:
                    assert m[i,j] == 0
                elif i % 2 == 1 and j % 2 == 0:
                    assert m[i,j] == 0
                else:
                    assert m[i,j] == old[i,j]
        
        hic = HicBasic(self.hic)
        # row
        old = hic[1,2:10]
        for i in range(0,8,2):
            old[i] = 0
        hic[1,2:10] = old
        
        assert np.array_equal(hic[1,:], [2,13,0,15,0,17,0,19,0,21,22,23])
        assert np.array_equal(hic[:,1], [2,13,0,15,0,17,0,19,0,21,22,23])
        
        hic = HicBasic(self.hic)
        # col
        old = hic[2:10,1]
        for i in range(0,8,2):
            old[i] = 0
        hic[2:10,1] = old
        
        assert np.array_equal(hic[1,:], [2,13,0,15,0,17,0,19,0,21,22,23])
        assert np.array_equal(hic[:,1], [2,13,0,15,0,17,0,19,0,21,22,23])
        
        # individual
        hic = HicBasic(self.hic)
        hic[2,1] = 0
        assert hic[2,1] == 0
        assert hic[1,2] == 0
    
    def test_as_data_frame(self):
        df = self.hic.as_data_frame(('chr1','chr1'))
        assert np.array_equal(df.shape, [5,5])
        assert np.array_equal(df.index, [1,1001,2001,3001,4001])
        assert np.array_equal(df.columns, [1,1001,2001,3001,4001])
        
    def test_from_mirny(self):
        # TODO
        pass
    
    def test_merge(self):
        hic = HicBasic()
        
        # add some nodes (120 to be exact)
        nodes = []
        for i in range(1,5000,1000):
            nodes.append(HicNode(chromosome="chr1",start=i,end=i+1000-1))
        for i in range(1,3000,1000):
            nodes.append(HicNode(chromosome="chr2",start=i,end=i+1000-1))
        for i in range(1,2000,400):
            nodes.append(HicNode(chromosome="chr4",start=i,end=i+100-1))
        hic.add_nodes(nodes)
        
        # add some edges with increasing weight for testing
        edges = []
        weight = 1
        for i in range(0,len(nodes)):
            for j in range(i,len(nodes)):
                edges.append(HicEdge(source=i,sink=j,weight=weight))
                weight += 1

        hic.add_edges(edges)
        
        left = self.hic[:,:]
        right = hic[:,:]
        
        # check length
        original_length = len(self.hic.nodes())
        self.hic.merge(hic)
        assert len(self.hic.nodes()) == original_length + 5
            
        merged = self.hic[:,:]
        double = [0,1,2,3,4,5,6,7]
        for i in double:
            for j in double:
                assert merged[i,j] == left[i,j] + right[i,j]
        
        three = [8,9,10,11]
        for i in double:
            for j in three:
                assert merged[i,j] == left[i,j]
        
        four = [12,13,14,15,16]
        for i in three:
            for j in four:
                assert merged[i,j] == 0
        
        for i in double:
            for j in four:
                assert merged[i,j] == right[i,j-4]

    def test_from_pairs(self):
        reads1 = Reads(self.dir + "/test_genomic/yeast.sample.chrI.1.sam")
        reads2 = Reads(self.dir + "/test_genomic/yeast.sample.chrI.2.sam")
        chrI = Chromosome.from_fasta(self.dir + "/test_genomic/chrI.fa")
        genome = Genome(chromosomes=[chrI])
        pairs = FragmentMappedReadPairs()
        pairs.load(reads1,reads2,genome.get_regions('HindIII'))
        
        pl = len(pairs)
        
        hic = HicBasic()
        hic._from_read_fragment_pairs(pairs, _max_buffer_size=1000)
        
        assert len(hic._regions) == len(pairs._regions)
        
        reads = 0
        edge_dict = {}
        for edge in hic.edges():
            edge_string = "%d-%d" % (edge.source, edge.sink) 
            assert edge_string not in edge_dict
            edge_dict[edge_string] = 1
            reads += edge.weight
        
        assert reads == pl
        
    def test_overlap_map(self):
        # ----|----|----|----|---|-----|-| new
        # -------|-------|-------|-------| old
        old_regions = []
        old_regions.append(HicNode(chromosome='chr1', start=1, end=8))
        old_regions.append(HicNode(chromosome='chr1', start=9, end=16))
        old_regions.append(HicNode(chromosome='chr1', start=17, end=24))
        old_regions.append(HicNode(chromosome='chr1', start=25, end=32))
        
        new_regions = []
        new_regions.append(HicNode(chromosome='chr1', start=1, end=5))
        new_regions.append(HicNode(chromosome='chr1', start=6, end=10))
        new_regions.append(HicNode(chromosome='chr1', start=11, end=15))
        new_regions.append(HicNode(chromosome='chr1', start=16, end=20))
        new_regions.append(HicNode(chromosome='chr1', start=21, end=24))
        new_regions.append(HicNode(chromosome='chr1', start=25, end=30))
        new_regions.append(HicNode(chromosome='chr1', start=31, end=32))
        
        overlap_map = _get_overlap_map(old_regions, new_regions)
        assert len(overlap_map[0]) == 2
        assert np.array_equal(overlap_map[0][0], [0,1.0])
        assert np.array_equal(overlap_map[0][1], [1,0.6])
        assert len(overlap_map[1]) == 3
        assert np.array_equal(overlap_map[1][0], [1,0.4])
        assert np.array_equal(overlap_map[1][1], [2,1.0])
        assert np.array_equal(overlap_map[1][2], [3,0.2])
        assert len(overlap_map[2]) == 2
        assert np.array_equal(overlap_map[2][0], [3,0.8])
        assert np.array_equal(overlap_map[2][1], [4,1.0])
        assert len(overlap_map[3]) == 2
        assert np.array_equal(overlap_map[3][0], [5,1.0])
        assert np.array_equal(overlap_map[3][1], [6,1.0])
        
        # ----|----|-| new
        # --|--|--|--| old
        old_regions = []
        old_regions.append(HicNode(chromosome='chr1', start=1, end=3))
        old_regions.append(HicNode(chromosome='chr1', start=4, end=6))
        old_regions.append(HicNode(chromosome='chr1', start=7, end=9))
        old_regions.append(HicNode(chromosome='chr1', start=10, end=12))
        
        new_regions = []
        new_regions.append(HicNode(chromosome='chr1', start=1, end=5))
        new_regions.append(HicNode(chromosome='chr1', start=6, end=10))
        new_regions.append(HicNode(chromosome='chr1', start=11, end=12))
        
        overlap_map = _get_overlap_map(old_regions, new_regions)
        assert len(overlap_map[0]) == 1
        assert np.array_equal(overlap_map[0][0], [0,0.6])
        assert len(overlap_map[1]) == 2
        assert np.array_equal(overlap_map[1][0], [0,0.4])
        assert np.array_equal(overlap_map[1][1], [1,0.2])
        assert len(overlap_map[2]) == 1
        assert np.array_equal(overlap_map[2][0], [1,0.6])
        assert len(overlap_map[3]) == 2
        assert np.array_equal(overlap_map[3][0], [1,0.2])
        assert np.array_equal(overlap_map[3][1], [2,1.0])
        
    def test_edge_splitting_rao(self):
        #     1         2         3         4
        # ---------|---------|---------|---------| old
        # ----|----|----|----|---------|---|--|--| new
        #  1    2    3    4       5      6  7  8
        overlap_map = {}
        overlap_map[1] = [[1,1.0], [2,1.0]]
        overlap_map[2] = [[3,1.0], [4,1.0]]
        overlap_map[3] = [[5,1.0]]
        overlap_map[4] = [[6,1.0],[7,1.0],[8,1.0]]
        
        original_edge = [1,2,12.0]
        new_edges = _edge_overlap_split_rao(original_edge,overlap_map)
        assert len(new_edges) == 4
        weight_sum = 0
        for new_edge in new_edges:
            weight_sum += new_edge[2]
        assert weight_sum == original_edge[2]
        
        original_edge = [1,1,12.0]
        new_edges = _edge_overlap_split_rao(original_edge,overlap_map)
        assert len(new_edges) == 3
        weight_sum = 0
        for new_edge in new_edges:
            weight_sum += new_edge[2]
        assert weight_sum == original_edge[2]
        
        original_edge = [1,3,9.0]
        new_edges = _edge_overlap_split_rao(original_edge,overlap_map)
        assert len(new_edges) == 2
        weight_sum = 0
        for new_edge in new_edges:
            weight_sum += new_edge[2]
        assert weight_sum == original_edge[2]
        
        original_edge = [3,3,9.0]
        new_edges = _edge_overlap_split_rao(original_edge,overlap_map)
        assert len(new_edges) == 1
        assert new_edges[0][2] == original_edge[2]
        
        original_edge = [1,4,9.0]
        new_edges = _edge_overlap_split_rao(original_edge,overlap_map)
        assert len(new_edges) == 4
        weight_sum = 0
        for new_edge in new_edges:
            weight_sum += new_edge[2]
        assert weight_sum == original_edge[2]
        
    def test_from_hic(self):
        reads1 = Reads(self.dir + "/test_genomic/yeast.sample.chrI.1.sam")
        reads2 = Reads(self.dir + "/test_genomic/yeast.sample.chrI.2.sam")
        chrI = Chromosome.from_fasta(self.dir + "/test_genomic/chrI.fa")
        genome = Genome(chromosomes=[chrI])
        pairs = FragmentMappedReadPairs()
        pairs.load(reads1,reads2,genome.get_regions('HindIII'))
                
        hic = HicBasic()
        hic._from_read_fragment_pairs(pairs, _max_buffer_size=1000)
        original_reads = 0
        for edge in hic.edges():
            original_reads += edge.weight
        
        
        def assert_binning(hic, bin_size, buffer_size):
            binned = HicBasic()
            assert len(binned.nodes()) == 0
            binned.add_regions(genome.get_regions(bin_size))
            binned._from_hic(hic, _edge_buffer_size=buffer_size)
            
            new_reads = 0
            for edge in binned.edges():
                new_reads += edge.weight
            
            # search for duplicated edges
            edge_dict = {}
            for edge in binned.edges():
                assert (edge.source, edge.sink) not in edge_dict
                
            # make sure that the total number
            # of reads stays the same
            assert original_reads == new_reads
            print len(binned.regions())
        
        bin_sizes = [500,1000,5000,10000,20000]
        buffer_sizes = [10,100,500,1000,10000,50000]
        for bin_size in bin_sizes:
            for buffer_size in buffer_sizes:
                assert_binning(hic, bin_size, buffer_size)
        
    def test_from_hic_sample(self):
        hic = HicBasic()
        hic.add_region(GenomicRegion(chromosome='chr1',start=1,end=100))
        hic.add_region(GenomicRegion(chromosome='chr1',start=101,end=200))
        hic.add_edge([0,0,12])
        hic.add_edge([0,1,36])
        hic.add_edge([1,1,24])
        
        binned = HicBasic()
        binned.add_region(GenomicRegion(chromosome='chr1', start=1, end=50))
        binned.add_region(GenomicRegion(chromosome='chr1', start=51, end=100))
        binned.add_region(GenomicRegion(chromosome='chr1', start=101, end=150))
        binned.add_region(GenomicRegion(chromosome='chr1', start=151, end=200))
        
        binned._from_hic(hic, _edge_buffer_size=1000)
        
        original_reads = 0
        for edge in hic.edges():
            original_reads += edge.weight
        
        new_reads = 0
        for edge in binned.edges():
            new_reads += edge.weight
        
        # search for duplicated edges
        edge_dict = {}
        for edge in binned.edges():
            assert (edge.source, edge.sink) not in edge_dict
        
        # make sure that the total number
        # of reads stays the same
        assert original_reads == new_reads
        
        