import os.path
import pytest
from kaic.construct.seq import Reads, FragmentMappedReadPairs,\
    FragmentRead, InwardPairsFilter, UnmappedFilter, OutwardPairsFilter,\
    ReDistanceFilter, FragmentReadPair, SelfLigationFilter, PCRDuplicateFilter,\
    LazyFragmentRead, ContaminantFilter
from kaic.data.genomic import Genome, GenomicRegion, Chromosome
import numpy as np


class TestReads:
    
    def setup_method(self, method):
        self.dir = os.path.dirname(os.path.realpath(__file__))
        self.sam1_file = self.dir + "/test_seq/test1.sam"
        self.sam2_file = self.dir + "/test_seq/test2.sam"
        self.lambda_sam1_file = self.dir + "/test_seq/lambda_reads1.sam"
        self.lambda_sam1_contaminant_file = self.dir + "/test_seq/lambda_reads1_contaminant.sam"
        
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
        reads.load(self.sam1_file)

        compare(reads[0],
                ['SRR038105.1', 0, 'chrXI', 128390, 35, [(0, 15)], -1, -1, 0, 'GATATGATGGATTTG', 'FDFFFFFFFFFFFCF', 9,
                 1])

        # SRR038105.1000167    0    chrXIV    703158    42    15M    *    0    0    TACGGTATTGGTCGG    FFFFCFFFFFFFFCF    AS:i:0    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:15    YT:Z:UU
        res = reads.get_read_by_qname('SRR038105.1000167')
        compare(res, ['SRR038105.1000167', 0, 'chrXIV', 703158, 42, [(0, 15)], -1, -1, 0, 'TACGGTATTGGTCGG',
                      'FFFFCFFFFFFFFCF', 8, 1])

        # SRR038105.1000320    0    chrXVI    577162    35    15M    *    0    0    TTGATAAAATAGTCC    <<@FF<FFFFAFAFA    AS:i:0    XS:i:-5    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:15    YT:Z:UU
        res = reads.get_read_by_qname('SRR038105.1000320')
        compare(res, ['SRR038105.1000320', 0, 'chrXVI', 577162, 35, [(0, 15)], -1, -1, 0, 'TTGATAAAATAGTCC',
                      '<<@FF<FFFFAFAFA', 9, 1])

        # check unpaired right
        # SRR038105.1000002    16    chrIV    203242    42    16M    *    0    0    ACCCATTATTTCTCGA    IIIIIFIICIFIIIII    AS:i:0    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:16    YT:Z:UU
        res = reads.get_read_by_qname('SRR038105.1000002')
        assert res is None

        # check unpaired left
        # SRR038105.1000011    16    chrIV    526796    42    16M    *    0    0    GGTGAATTAGAAGATA    FFFFFFFFFFFFFFFF    AS:i:0    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:16    YT:Z:UU
        res = reads.get_read_by_qname('SRR038105.1000011')
        compare(res, ['SRR038105.1000011', 16, 'chrIV', 526796, 42, [(0, 16)], -1, -1, 0, 'GGTGAATTAGAAGATA',
                      'FFFFFFFFFFFFFFFF', 8, -1])

        reads.close()
        
    def test_ix(self):
        pairs = Reads(self.sam1_file)
        i = 0
        for pair in pairs._reads:
            assert pair['ix'] == i
            i += 1
        pairs.close()
    
    def test_strand(self):
        reads = Reads(self.sam1_file)
        for read in reads:
            if read.flag == 16:
                assert read.strand == -1
            if read.flag == 4:
                assert read.strand == 1
        reads.close()
    
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
        pairs.close()
    
    def test_select(self):
        reads = Reads(self.sam1_file)
        
        assert reads[0].qname == 'SRR038105.1'
        reads.filter_non_unique()
        assert reads[0].qname == 'SRR038105.10'
        reads.close()
        
    def test_build_from_scratch(self):
        reads = Reads()
        reads.load(self.sam1_file)
        
        assert len(reads) == 271
        reads.close()

    def test_read_alen(self):
        reads = Reads(self.sam1_file)
        read = reads[0]
        read.cigar = [(4,22), (0,42), (2,1), (0,18), (1,4), (0,5)]
        assert read.alen == 42 + 18 + 5
        reads.close()

    def test_read_get_tag(self):
        reads = Reads(self.sam1_file)
        read = reads[0]
        assert read.get_tag('AS') == 0
        assert read.get_tag('MD') == '15'
        assert read.get_tag('X0') is None
        reads.close()
    
    def test_quality_filter(self):
        reads = Reads(self.sam1_file)
        
        reads.filter_quality(30, queue=False)
        for row in reads._reads._iter_visible_and_masked():
            if row['mapq'] < 30:
                assert row[reads._reads._mask_field] == 2
            else:
                assert row[reads._reads._mask_field] == 0
        reads.close()
        
    def test_uniqueness_filter(self):
        reads = Reads(self.sam1_file)
        
        reads.filter_non_unique(strict=True)
        for row in reads._reads._iter_visible_and_masked():
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
        reads.close()
    
    def test_unmapped_filter(self):
        reads = Reads(self.lambda_sam1_file)
        
        l = len(reads)
        
        unmapped_filter = UnmappedFilter(reads.add_mask_description("unmapped", "Filter unmapped reads"))
        reads.filter(unmapped_filter)
        
        assert len(reads) < l
        reads.close()

    def test_contaminant_filter(self):
        reads = Reads(self.lambda_sam1_file)

        l = len(reads)

        contaminant = Reads(self.lambda_sam1_contaminant_file)
        contaminant_filter = ContaminantFilter(contaminant)

        reads.filter(contaminant_filter)

        assert len(reads) == l-9
        reads.close()
        contaminant.close()

    def test_queue_filters(self):
        reads = Reads(self.sam1_file)
        
        l = len(reads)
        
        reads.filter_quality(30, queue=True)
        reads.filter_non_unique(strict=True, queue=True)
        
        assert len(reads) == l
        
        reads.run_queued_filters()
        
        assert len(reads) < l
        reads.close()

    def test_iter_qname_sorted(self):
        reads = Reads(self.sam1_file)
        previous = 0
        for read in reads.reads(sort_by_qname_ix=True):
            assert read.qname_ix > previous
            previous = read.qname_ix
        assert previous != 0
        reads.close()

    def test_iterate_exclude_filters(self):
        reads = Reads(self.sam1_file)
        reads.filter_unmapped(queue=True)
        reads.filter_quality(35, queue=True)
        reads.filter_non_unique(strict=True, queue=True)
        reads.run_queued_filters()
        assert len(list(reads.reads(excluded_filters=['unmapped', 'mapq', 'uniqueness']))) == 271
        assert len(list(reads.reads(excluded_filters=['unmapped', 'uniqueness']))) == 246
        assert len(list(reads.reads(excluded_filters=['unmapped', 'mapq']))) == 153
        assert len(list(reads.reads(excluded_filters=['mapq', 'uniqueness']))) == 271
        assert len(list(reads.reads(excluded_filters=['unmapped']))) == 153
        assert len(list(reads.reads(excluded_filters=['mapq']))) == 153
        assert len(list(reads.reads(excluded_filters=['uniqueness']))) == 246
        reads.close()


class TestFileOpsReads:
    def test_tmp_with(self, tmpdir):
        filename = str(tmpdir) + "/test.file"
        with Reads(file_name=filename, mode='a', tmpdir='/tmp') as f:
            assert os.path.isfile(filename) == False
            assert os.path.isfile(f.tmp_file_name) == True
        assert os.path.isfile(filename) == True
        assert os.path.isfile(f.tmp_file_name) == False

    def test_tmp_with_exception(self, tmpdir):
        filename = str(tmpdir) + "/test.file"
        with pytest.raises(Exception):
            with Reads(file_name=filename, mode='a', tmpdir='/tmp') as f:
                assert os.path.isfile(filename) == False
                assert os.path.isfile(f.tmp_file_name) == True
                try:
                    raise Exception
                except:
                    assert os.path.isfile(filename) == False
                    assert os.path.isfile(f.tmp_file_name) == False


class TestBWAReads:
    @classmethod
    def setup_method(self, method):
        self.dir = os.path.dirname(os.path.realpath(__file__))
        self.bwamem_sam1_file = self.dir + "/test_seq/test_bwa1.sam"
        self.bwamem_sam2_file = self.dir + "/test_seq/test_bwa2.sam"

    def test_infer_mapper(self):
        with Reads(self.bwamem_sam1_file) as reads:
            assert reads.mapper == 'bwa'
        with Reads(self.bwamem_sam1_file, mapper='bowtie2') as reads:
            assert reads.mapper == 'bowtie2'

    def test_bwamem_quality_filter(self):
        with Reads(self.bwamem_sam1_file) as reads:
            assert len(reads) == 995
            reads.filter_quality(cutoff=0.90, queue=False)
            assert len(reads) == 927
            for read in reads:
                assert float(read.get_tag('AS')) / read.alen >= 0.90

    def test_bwamem_uniqueness_filter(self):
        with Reads(self.bwamem_sam1_file) as reads:
            assert len(reads) == 995
            reads.filter_non_unique(cutoff=3, queue=False)
            assert len(reads) == 626
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
        regions = self.genome.get_regions(1000)
        self.pairs.load(self.reads1, self.reads2, regions=regions)
        regions.close()

    def teardown_method(self, method):
        self.reads1.close()
        self.reads2.close()
        self.genome.close()
        self.pairs.close()
        
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
        ad = FragmentMappedReadPairs._auto_dist
        np.random.seed(101)
        x = np.linspace(1, 100, 10)
        i = [4.0, 3.0, 1.0, 1.2, 0.4, 0.8, 0.5, 0.4, 0.6, 0.5, 0.6, 0.5, 0.5, 0.5]
        o = [0.2, 0.4, 0.3, 0.6, 0.5, 0.4, 0.5, 0.5, 0.6, 0.5, 0.5, 0.4, 0.5, 0.5]
        b = np.random.normal(500, 200, 14).astype(int)
        assert ad(x, i, b, 0.05) == 67
        assert ad(x, o, b, 0.05) == 45

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

    def test_filter_self_ligated(self):
        mask = self.pairs.add_mask_description('self_ligated', 'Mask read pairs that represent self-ligated fragments')
        self_ligation_filter = SelfLigationFilter(mask=mask)

        assert len(self.pairs) == 44
        self.pairs.filter(self_ligation_filter)
        assert len(self.pairs) == 7

    def test_get_ligation_structure_biases(self):
        reads1 = Reads(self.dir + "/../data/test_genomic/yeast.sample.chrI.1.sam")
        reads2 = Reads(self.dir + "/../data/test_genomic/yeast.sample.chrI.2.sam")
        chrI = Chromosome.from_fasta(self.dir + "/../data/test_genomic/chrI.fa")
        genome = Genome(chromosomes=[chrI])
        pairs = FragmentMappedReadPairs()
        regions = genome.get_regions('HindIII')
        pairs.load(reads1, reads2, regions)
        reads1.close()
        reads2.close()
        genome.close()
        regions.close()
        x, i, o, b = pairs.get_ligation_structure_biases(sampling=200, skip_self_ligations=False)
        assert len(x) == len(i) == len(o) == len(b) == 3
        assert x.tolist() == [494, 4487, 19399]
        assert i.tolist() == [2.616915422885572, 0.8059701492537313, 0.6417910447761194]
        assert o.tolist() == [0.2537313432835821, 0.24378109452736318, 0.46766169154228854]
        assert b.tolist() == [778, 412, 424]
        pairs.close()

    def test_re_dist(self):
        read1 = FragmentRead(GenomicRegion(chromosome='chr1', start=1, end=1000), position=200, strand=-1)
        assert read1.re_distance() == 199
        read2 = FragmentRead(GenomicRegion(chromosome='chr1', start=1, end=1000), position=990, strand=-1)
        assert read2.re_distance() == 10

    def test_iterate_exclude_filters(self):
        mask = self.pairs.add_mask_description('inwards', 'Mask read pairs that inward facing and closer than 100bp')
        in_filter = InwardPairsFilter(minimum_distance=100, mask=mask)
        self.pairs.filter(in_filter)
        mask = self.pairs.add_mask_description('outwards', 'Mask read pairs that outward facing and closer than 100bp')
        out_filter = OutwardPairsFilter(minimum_distance=100, mask=mask)
        self.pairs.filter(out_filter)
        mask = self.pairs.add_mask_description('re-dist', 'Mask read pairs where one half maps more than 100bp away from both RE sites')
        re_filter = ReDistanceFilter(maximum_distance=300, mask=mask)
        self.pairs.filter(re_filter)
        assert len(self.pairs) == 0
        assert len(list(self.pairs.pairs(excluded_filters=['inwards', 'outwards', 're-dist']))) == 44
        assert len(list(self.pairs.pairs(excluded_filters=['inwards', 'outwards']))) == 13
        assert len(list(self.pairs.pairs(excluded_filters=['inwards', 're-dist']))) == 28
        assert len(list(self.pairs.pairs(excluded_filters=['outwards', 're-dist']))) == 18
        assert len(list(self.pairs.pairs(excluded_filters=['inwards']))) == 11
        assert len(list(self.pairs.pairs(excluded_filters=['outwards']))) == 2
        assert len(list(self.pairs.pairs(excluded_filters=['re-dist']))) == 2
        assert len(list(self.pairs.pairs(excluded_filters=[in_filter, re_filter]))) == 28
        assert len(list(self.pairs.pairs(excluded_filters=[in_filter, mask]))) == 28
        assert len(list(self.pairs.pairs(excluded_filters=[in_filter, 3]))) == 28


class TestFragmentRead:
    def setup_method(self, method):
        fragment1 = GenomicRegion(1, 1000, chromosome='chr1', strand=1, ix=0)
        fragment2 = GenomicRegion(1001, 2000, chromosome='chr2', strand=-1, ix=1)
        self.read1 = FragmentRead(fragment1, position=500, strand=1, qname_ix=1)
        self.read2 = FragmentRead(fragment2, position=1200, strand=1, qname_ix=2)

        class DummyPairs(object):
            def __init__(self):
                self._ix_to_chromosome = {
                    0: 'chr1',
                    1: 'chr2'
                }

        row = dict()
        row['ix'] = 0
        row['left_read_qname_ix'] = 1
        row['left_read_position'] = 500
        row['left_read_strand'] = 1
        row['left_fragment'] = 0
        row['left_fragment_start'] = 1
        row['left_fragment_end'] = 1000
        row['left_fragment_chromosome'] = 0
        row['right_read_qname_ix'] = 2
        row['right_read_position'] = 1200
        row['right_read_strand'] = -1
        row['right_fragment'] = 1
        row['right_fragment_start'] = 1001
        row['right_fragment_end'] = 2000
        row['right_fragment_chromosome'] = 1
        dummy_pairs = DummyPairs()
        self.lazy_read1 = LazyFragmentRead(row, dummy_pairs, side='left')
        self.lazy_read2 = LazyFragmentRead(row, dummy_pairs, side='right')

    def test_attributes(self):
        assert isinstance(self.read1.fragment, GenomicRegion)
        assert self.read1.fragment.chromosome == 'chr1'
        assert self.read1.fragment.start == 1
        assert self.read1.fragment.end == 1000
        assert self.read1.fragment.strand == 1
        assert self.read1.position == 500
        assert self.read1.strand == 1
        assert self.read1.qname_ix == 1

    def test_lazy_attributes(self):
        assert isinstance(self.lazy_read1.fragment, GenomicRegion)
        assert self.lazy_read1.fragment.chromosome == 'chr1'
        assert self.lazy_read1.fragment.start == 1
        assert self.lazy_read1.fragment.end == 1000
        assert self.lazy_read1.fragment.strand == 1
        assert self.lazy_read1.position == 500
        assert self.lazy_read1.strand == 1
        assert self.lazy_read1.qname_ix == 1

        assert isinstance(self.lazy_read2.fragment, GenomicRegion)
        assert self.lazy_read2.fragment.chromosome == 'chr2'
        assert self.lazy_read2.fragment.start == 1001
        assert self.lazy_read2.fragment.end == 2000
        assert self.lazy_read2.fragment.strand == 1
        assert self.lazy_read2.position == 1200
        assert self.lazy_read2.strand == -1
        assert self.lazy_read2.qname_ix == 2


class TestFragmentReadPair:

    def test_convenience_functions(self):
        pair = FragmentReadPair(
            FragmentRead(
                GenomicRegion(1, 1000, chromosome='chr1'),
                position=500, strand=1
            ),
            FragmentRead(
                GenomicRegion(10001, 11000, chromosome='chr1'),
                position=10500, strand=-1
            )
        )
        assert pair.is_same_chromosome()
        assert pair.get_gap_size() == 9001
        assert pair.is_inward_pair()
        assert not pair.is_outward_pair()
        assert not pair.is_same_fragment()
        assert not pair.is_same_pair()

        pair = FragmentReadPair(
            FragmentRead(
                GenomicRegion(1, 1000, chromosome='chr1'),
                position=500, strand=1
            ),
            FragmentRead(
                GenomicRegion(1, 1000, chromosome='chr1'),
                position=600, strand=1
            )
        )
        assert pair.is_same_chromosome()
        assert pair.get_gap_size() == 0
        assert not pair.is_inward_pair()
        assert not pair.is_outward_pair()
        assert pair.is_same_fragment()
        assert pair.is_same_pair()

        pair = FragmentReadPair(
            FragmentRead(
                GenomicRegion(1, 1000, chromosome='chr1'),
                position=500, strand=-1
            ),
            FragmentRead(
                GenomicRegion(1, 1000, chromosome='chr2'),
                position=600, strand=1
            )
        )
        assert not pair.is_same_chromosome()
        assert pair.get_gap_size() is None
        assert not pair.is_inward_pair()
        assert not pair.is_outward_pair()
        assert not pair.is_same_fragment()
        assert not pair.is_same_pair()

        pair = FragmentReadPair(
            FragmentRead(
                GenomicRegion(1, 1000, chromosome='chr1'),
                position=500, strand=-1
            ),
            FragmentRead(
                GenomicRegion(1001, 2000, chromosome='chr1'),
                position=1200, strand=1
            )
        )
        assert pair.is_same_chromosome()
        assert pair.get_gap_size() == 0
        assert not pair.is_inward_pair()
        assert pair.is_outward_pair()
        assert not pair.is_same_fragment()
        assert not pair.is_same_pair()


class TestFileOpsFragmentMappedReadPairs:
    def test_tmp_with(self, tmpdir):
        filename = str(tmpdir) + "/test.file"
        with FragmentMappedReadPairs(file_name=filename, mode='a', tmpdir='/tmp') as f:
            assert os.path.isfile(filename) == False
            assert os.path.isfile(f.tmp_file_name) == True
        assert os.path.isfile(filename) == True
        assert os.path.isfile(f.tmp_file_name) == False

    def test_tmp_with_exception(self, tmpdir):
        filename = str(tmpdir) + "/test.file"
        with pytest.raises(Exception):
            with FragmentMappedReadPairs(file_name=filename, mode='a', tmpdir='/tmp') as f:
                assert os.path.isfile(filename) == False
                assert os.path.isfile(f.tmp_file_name) == True
                try:
                    raise Exception
                except:
                    assert os.path.isfile(filename) == False
                    assert os.path.isfile(f.tmp_file_name) == False


class TestBWAFragmentMappedReads:
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
        regions = self.genome.get_regions('MboI')
        self.pairs.load(self.reads1, self.reads2, regions=regions)
        regions.close()

    def teardown_method(self, method):
        self.reads1.close()
        self.reads2.close()
        self.genome.close()
        self.pairs.close()
        
    def test_loaded_bwamem_pairs(self):
        assert self.pairs._single_count == 896
        assert self.pairs._pair_count == 515

    def test_pcr_duplicate_filter(self):
        mask = self.pairs.add_mask_description('pcr_duplicate', 'Mask suspected PCR duplicates.')
        pcr_duplicate_filter = PCRDuplicateFilter(pairs=self.pairs, threshold=3)

        assert len(self.pairs) == 515
        self.pairs.filter(pcr_duplicate_filter)
        assert len(self.pairs) == 512
