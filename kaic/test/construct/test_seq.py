'''
Created on Jul 15, 2015

@author: kkruse1
'''

import os.path
import numpy as np
from kaic.construct.seq import Reads, FragmentMappedReadPairs,\
    FragmentRead, InwardPairsFilter, UnmappedFilter, OutwardPairsFilter,\
    ReDistanceFilter, FragmentReadPair
from kaic.data.genomic import Genome, GenomicRegion, Chromosome
import numpy as np


class TestReads:
    
    @classmethod
    def setup_method(self, method):
        self.dir = os.path.dirname(os.path.realpath(__file__))
        self.sam1_file = self.dir + "/test_seq/test1.sam"
        self.sam2_file = self.dir + "/test_seq/test2.sam"
        self.lambda_sam1_file = self.dir + "/test_seq/lambda_reads1.sam"
        
    def test_load(self):
        def compare(read, values):
            assert read.qname == values[0]
            assert read.flag == values[1]
            assert read.ref == values[2]
            assert read.pos == values[3]
            assert read.mapq == values[4]
            assert np.array_equal(read.cigar, values[5])
            assert read.rnext == values[6]
            assert read.pnext == values[7]
            assert read.tlen == values[8]
            assert read.seq == values[9]
            assert read.qual == values[10]
            assert len(read.tags) == values[11]
            assert read.strand == values[12]
        
        reads = Reads()
        reads.load(self.sam1_file, is_sorted=False)
        
        compare(reads[0], ['SRR038105.1',0,'chrXI',128390,35,[(0,15)],-1,-1,0,'GATATGATGGATTTG','FDFFFFFFFFFFFCF',9,1])
        
        # SRR038105.1000167    0    chrXIV    703158    42    15M    *    0    0    TACGGTATTGGTCGG    FFFFCFFFFFFFFCF    AS:i:0    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:15    YT:Z:UU
        res = reads.get_read_by_qname('SRR038105.1000167')
        compare(res, ['SRR038105.1000167',0,'chrXIV',703158,42,[(0,15)],-1,-1,0,'TACGGTATTGGTCGG','FFFFCFFFFFFFFCF',8,1])
        
        # SRR038105.1000320    0    chrXVI    577162    35    15M    *    0    0    TTGATAAAATAGTCC    <<@FF<FFFFAFAFA    AS:i:0    XS:i:-5    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:15    YT:Z:UU
        res = reads.get_read_by_qname('SRR038105.1000320')
        compare(res, ['SRR038105.1000320',0,'chrXVI',577162,35,[(0,15)],-1,-1,0,'TTGATAAAATAGTCC','<<@FF<FFFFAFAFA',9,1])

        # check unpaired right
        # SRR038105.1000002    16    chrIV    203242    42    16M    *    0    0    ACCCATTATTTCTCGA    IIIIIFIICIFIIIII    AS:i:0    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:16    YT:Z:UU
        res = reads.get_read_by_qname('SRR038105.1000002')
        assert res is None
        
        # check unpaired left
        # SRR038105.1000011    16    chrIV    526796    42    16M    *    0    0    GGTGAATTAGAAGATA    FFFFFFFFFFFFFFFF    AS:i:0    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:16    YT:Z:UU
        res = reads.get_read_by_qname('SRR038105.1000011')
        compare(res, ['SRR038105.1000011',16,'chrIV',526796,42,[(0,16)],-1,-1,0,'GGTGAATTAGAAGATA','FFFFFFFFFFFFFFFF',8,-1])
        
    def test_ix(self):
        pairs = Reads(self.sam1_file)
        i = 0
        for pair in pairs._reads:
            assert pair['ix'] == i
            i += 1
    
    def test_strand(self):
        reads = Reads(self.sam1_file)
        for read in reads:
            if read.flag == 16:
                assert read.strand == -1
            if read.flag == 4:
                assert read.strand == 1
    
    def test_iter(self):
        pairs = Reads(self.sam1_file)
        counter = 0
        for _ in pairs:
            counter += 1
        
        assert counter == 271
        
        pairs.filter_non_unique()
        after_counter = 0
        for _ in pairs:
            after_counter += 1
            
        assert after_counter < counter
    
    def test_select(self):
        reads = Reads(self.sam1_file)
        
        assert reads[0].qname == 'SRR038105.1'
        reads.filter_non_unique()
        assert reads[0].qname == 'SRR038105.10'
        
    def test_build_from_scratch(self):
        reads = Reads()
        reads.load(self.sam1_file)
        
        assert len(reads) == 271

    def test_read_alen(self):
        reads = Reads(self.sam1_file)
        read = reads[0]
        read.cigar = [(4,22), (0,42), (2,1), (0,18), (1,4), (0,5)]
        assert read.alen == 42 + 18 + 5

    def test_read_get_tag(self):
        reads = Reads(self.sam1_file)
        read = reads[0]
        assert read.get_tag('AS') == 0
        assert read.get_tag('MD') == '15'
        assert read.get_tag('X0') == None
    
    def test_quality_filter(self):
        reads = Reads(self.sam1_file)
        
        reads.filter_quality(30, queue=False)
        for row in reads._reads.all():
            if row['mapq'] < 30:
                assert row[reads._reads._mask_field] == 2
            else:
                assert row[reads._reads._mask_field] == 0
        
    def test_uniqueness_filter(self):
        reads = Reads(self.sam1_file)
        
        reads.filter_non_unique(strict=True)
        for row in reads._reads.all():
            if row['pos'] > 0:
                tags = reads._tags[row['ix']]
                has_xs = False
                for tag in tags:
                    if tag[0] == 'XS':
                        has_xs = True
                        break
                if has_xs:
                    assert row[reads._reads._mask_field] == 2
            else:
                assert row[reads._reads._mask_field] == 0
    
    def test_unmapped_filter(self):
        reads = Reads(self.lambda_sam1_file)
        
        l = len(reads)
        
        unmapped_filter = UnmappedFilter(reads.add_mask_description("unmapped", "Filter unmapped reads"))
        reads.filter(unmapped_filter)
        
        assert len(reads) < l

    def test_queue_filters(self):
        reads = Reads(self.sam1_file)
        
        l = len(reads)
        
        reads.filter_quality(30, queue=True)
        reads.filter_non_unique(strict=True, queue=True)
        
        assert len(reads) == l
        
        reads.run_queued_filters()
        
        assert len(reads) < l

    def test_iter_qname_sorted(self):
        reads = Reads(self.sam1_file)
        previous = 0
        for read in reads.reads(sort_by_qname_ix=True):
            assert read.qname_ix > previous
            previous = read.qname_ix
        assert previous != 0


class TestBWAReads:
    @classmethod
    def setup_method(self, method):
        self.dir = os.path.dirname(os.path.realpath(__file__))
        self.bwamem_sam1_file = self.dir + "/test_seq/test_bwa1.sam"
        self.bwamem_sam2_file = self.dir + "/test_seq/test_bwa2.sam"

    def test_infer_mapper(self):
        reads = Reads(self.bwamem_sam1_file)
        assert reads.mapper == 'bwa'
        reads = Reads(self.bwamem_sam1_file, mapper='bowtie2')
        assert reads.mapper == 'bowtie2'

    def test_bwamem_quality_filter(self):
        reads = Reads(self.bwamem_sam1_file)
        assert len(reads) == 992
        reads.filter_quality(cutoff=0.90, queue=False)
        assert len(reads) == 924
        for read in reads:
            assert float(read.get_tag('AS')) / read.alen >= 0.90

    def test_bwamem_uniqueness_filter(self):
        reads = Reads(self.bwamem_sam1_file)
        assert len(reads) == 992
        reads.filter_non_unique(cutoff=3, queue=False)
        assert len(reads) == 623
        for read in reads:
            assert read.mapq > 3


class TestFragmentMappedReads:
    @classmethod
    def setup_method(self, method):
        self.dir = os.path.dirname(os.path.realpath(__file__))
        sam1_file = self.dir + "/test_seq/lambda_reads1.sam"
        sam2_file = self.dir + "/test_seq/lambda_reads2.sam"
        self.reads1 = Reads(sam1_file)
        self.reads2 = Reads(sam2_file)
        self.reads1.filter_unmapped()
        self.reads2.filter_unmapped()
        self.genome = Genome.from_folder(self.dir + "/test_seq/lambda_genome/")
        
        self.pairs = FragmentMappedReadPairs()
        self.pairs.load(self.reads1, self.reads2, regions=self.genome.get_regions(1000))
        
    def test_select(self):
        pair = self.pairs[0]
        assert isinstance(pair, FragmentReadPair)
        assert pair.left.position == 40884
        assert pair.right.position == 41039
        assert pair.left.strand == 1
        assert pair.right.strand == -1
        assert isinstance(pair.left.fragment, GenomicRegion)
        assert isinstance(pair.right.fragment, GenomicRegion)
        assert pair.left.fragment.start == 40001
        assert pair.left.fragment.end == 41000
        assert pair.left.fragment.chromosome == 'gi|9626243|ref|NC_001416.1|'
        assert pair.right.fragment.start == 41001
        assert pair.right.fragment.end == 42000
        assert pair.right.fragment.chromosome == 'gi|9626243|ref|NC_001416.1|'
        
        pair = self.pairs[1]
        assert isinstance(pair, FragmentReadPair)
        assert pair.left.position == 19617
        assert pair.right.position == 19736
        assert pair.left.strand == 1
        assert pair.right.strand == -1
        assert isinstance(pair.left.fragment, GenomicRegion)
        assert isinstance(pair.right.fragment, GenomicRegion)
        assert pair.left.fragment.start == 19001
        assert pair.left.fragment.end == 20000
        assert pair.left.fragment.chromosome == 'gi|9626243|ref|NC_001416.1|'
        assert pair.right.fragment.start == 19001
        assert pair.right.fragment.end == 20000
        assert pair.right.fragment.chromosome == 'gi|9626243|ref|NC_001416.1|'
        
        pair = self.pairs[-1]
        assert isinstance(pair, FragmentReadPair)
        assert pair.left.position == 25765
        assert pair.right.position == 25622
        assert pair.left.strand == -1
        assert pair.right.strand == 1
        assert isinstance(pair.left.fragment, GenomicRegion)
        assert isinstance(pair.right.fragment, GenomicRegion)
        assert pair.left.fragment.start == 25001
        assert pair.left.fragment.end == 26000
        assert pair.left.fragment.chromosome == 'gi|9626243|ref|NC_001416.1|'
        assert pair.right.fragment.start == 25001
        assert pair.right.fragment.end == 26000
        assert pair.right.fragment.chromosome == 'gi|9626243|ref|NC_001416.1|'
    
    def test_iter(self):
        for pair in self.pairs:
            assert isinstance(pair, FragmentReadPair)
            assert isinstance(pair.left, FragmentRead)
            assert isinstance(pair.right, FragmentRead)
            assert isinstance(pair.left.fragment, GenomicRegion)
            assert isinstance(pair.right.fragment, GenomicRegion)
            assert pair.left.position > 0 or pair.left.position == -1
            assert pair.right.position > 0 or pair.right.position == -1
            assert pair.left.strand == -1 or pair.left.strand == 1
            assert pair.right.strand == -1 or pair.right.strand == 1
            
            assert 0 < pair.left.fragment.start <= pair.left.fragment.end
            assert 0 < pair.right.fragment.start <= pair.right.fragment.end
            if pair.left.position > 0:
                assert pair.left.fragment.start <= pair.left.position <= pair.left.fragment.end
            if pair.right.position > 0:
                assert pair.right.fragment.start <= pair.right.position <= pair.right.fragment.end
            
    def test_len(self):
        assert len(self.pairs) == 44
            
    def test_single(self):
        assert len(self.pairs._single) == 6

    def test_auto_mindist(self):
        pairs = FragmentMappedReadPairs()
        dists = range(10)
        inward_ratios = [4.00, 2.00, 1.00, 0.50, 1.00, 0.40, 0.60, 0.55, 0.51, 0.50]
        outward_ratios = [0.10, 0.20, 0.40, 0.30, 0.55, 0.50, 0.45, 0.55, 0.49, 0.50]
        sane_inward = pairs._auto_dist(dists, inward_ratios, 0.2, 0.2, 3)
        sane_outward = pairs._auto_dist(dists, outward_ratios, 0.2, 0.2, 3)
        insane_inward = pairs._auto_dist(dists, inward_ratios, 0.001, 0.001, 3)
        insane_outward = pairs._auto_dist(dists, outward_ratios, 0.001, 0.001, 3)
        assert sane_inward == dists[8]
        assert sane_outward == dists[6]
        assert insane_inward is None
        assert insane_outward is None

    def test_filter_inward(self):
        mask = self.pairs.add_mask_description('inwards', 'Mask read pairs that inward facing and closer than 100bp')
        in_filter = InwardPairsFilter(minimum_distance=100, mask=mask)
        
        assert len(self.pairs) == 44
        self.pairs.filter(in_filter)
        
#         print "Valid pairs:"
#         for pair in self.pairs:
#             print pair[0]
#             print pair[1]
        assert len(self.pairs) == 18

    def test_filter_outward(self):
        mask = self.pairs.add_mask_description('outwards', 'Mask read pairs that outward facing and closer than 100bp')
        out_filter = OutwardPairsFilter(minimum_distance=100, mask=mask)
        
        assert len(self.pairs) == 44
        self.pairs.filter(out_filter)
        assert len(self.pairs) == 28
        
    def test_filter_redist(self):
        mask = self.pairs.add_mask_description('re-dist', 'Mask read pairs where one half maps more than 100bp away from both RE sites')
        re_filter = ReDistanceFilter(maximum_distance=300, mask=mask)
        
        assert len(self.pairs) == 44
        self.pairs.filter(re_filter)
        assert len(self.pairs) == 13

    def test_get_error_structure(self):
        reads1 = Reads(self.dir + "/../data/test_genomic/yeast.sample.chrI.1.sam")
        reads2 = Reads(self.dir + "/../data/test_genomic/yeast.sample.chrI.2.sam")
        chrI = Chromosome.from_fasta(self.dir + "/../data/test_genomic/chrI.fa")
        genome = Genome(chromosomes=[chrI])
        pairs = FragmentMappedReadPairs()
        pairs.load(reads1, reads2, genome.get_regions('HindIII'))
        x, i, o = pairs.get_error_structure(data_points=200, skip_self_ligations=False)
        assert len(x) == len(i)
        assert len(i) == len(o)
        assert len(o) == 3
        assert x == [494.03856041131104, 4487.5800970873788, 19399.908018867925]
        assert i == [2.616915422885572, 0.8059701492537313, 0.6417910447761194]
        assert o == [0.2537313432835821, 0.24378109452736318, 0.46766169154228854]

    def test_re_dist(self):
        read1 = FragmentRead(GenomicRegion(chromosome='chr1', start=1, end=1000), position=200, strand=-1)
        assert read1.re_distance() == 199
        read2 = FragmentRead(GenomicRegion(chromosome='chr1', start=1, end=1000), position=990, strand=-1)
        assert read2.re_distance() == 10


class TestBWAFragmentMappedReads:
    @classmethod
    def setup_method(self, method):
        self.dir = os.path.dirname(os.path.realpath(__file__))
        sam1_file = self.dir + "/test_seq/test_bwa1.sam"
        sam2_file = self.dir + "/test_seq/test_bwa2.sam"
        self.reads1 = Reads()
        self.reads2 = Reads()
        self.reads1.load(sam1_file)
        self.reads2.load(sam2_file)
        self.reads1.filter_unmapped()
        self.reads2.filter_unmapped()
        self.genome = Genome.from_folder(self.dir + "/test_seq/dmel_genome/")
        self.pairs = FragmentMappedReadPairs()
        self.pairs.load(self.reads1, self.reads2, regions=self.genome.get_regions('MboI'))
        
    def test_loaded_bwamem_pairs(self):
        assert self.pairs._single_count == 896
        assert self.pairs._pair_count == 512
