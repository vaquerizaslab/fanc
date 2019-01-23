import os.path
import pytest
from kaic.pairs import SamBamReadPairGenerator, ReadPairs, UnmappedFilter, FragmentReadPair, \
    FragmentRead, InwardPairsFilter, OutwardPairsFilter, ContaminantFilter, QualityFilter, \
    BwaMemQualityFilter, ReDistanceFilter, SelfLigationFilter, LazyFragment, LazyFragmentRead, \
    PCRDuplicateFilter
from genomic_regions import GenomicRegion
from kaic.regions import Genome, Chromosome
from kaic.general import Mask
import numpy as np


class TestReadPairs:
    @classmethod
    def setup_method(self, method):
        self.dir = os.path.dirname(os.path.realpath(__file__))
        sam1_file = os.path.join(self.dir, "test_pairs", "lambda_reads1_sort.sam")
        sam2_file = os.path.join(self.dir, "test_pairs", "lambda_reads2_sort.sam")
        pair_generator = SamBamReadPairGenerator(sam1_file, sam2_file)
        self.pairs = ReadPairs()
        f = UnmappedFilter(mask=Mask(ix=0, name='unmapped'))
        pair_generator.add_filter(f)
        self.genome = Genome.from_folder(os.path.join(self.dir, "test_pairs", "lambda_genome"))
        regions = self.genome.get_regions(1000)
        self.pairs.add_regions(regions.regions)
        regions.close()
        self.pairs.add_read_pairs(pair_generator)
        self.pairs_class = ReadPairs

    def teardown_method(self, method):
        self.genome.close()
        self.pairs.close()

    def test_select(self):
        pair = self.pairs[11]
        assert isinstance(pair, FragmentReadPair)
        assert pair.left.position == 40883
        assert pair.right.position == 41038
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

        pair = self.pairs[21]
        assert isinstance(pair, FragmentReadPair)
        assert pair.left.position == 19616
        assert pair.right.position == 19735
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
        assert pair.left.position == 5066
        assert pair.right.position == 5199
        assert pair.left.strand == 1
        assert pair.right.strand == -1
        assert isinstance(pair.left.fragment, GenomicRegion)
        assert isinstance(pair.right.fragment, GenomicRegion)
        assert pair.left.fragment.start == 5001
        assert pair.left.fragment.end == 6000
        assert pair.left.fragment.chromosome == 'gi|9626243|ref|NC_001416.1|'
        assert pair.right.fragment.start == 5001
        assert pair.right.fragment.end == 6000
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

    def test_auto_mindist(self):
        ad = self.pairs_class._auto_dist
        np.random.seed(101)
        x = np.linspace(1, 100, 10)
        i = [4.0, 3.0, 1.0, 1.2, 0.4, 0.8, 0.5, 0.4, 0.6, 0.5, 0.6, 0.5, 0.5, 0.5]
        o = [0.2, 0.4, 0.3, 0.6, 0.5, 0.4, 0.5, 0.5, 0.6, 0.5, 0.5, 0.4, 0.5, 0.5]
        b = np.random.normal(500, 200, 14).astype(int)
        assert ad(x, i, b, 0.05) == 67
        assert ad(x, o, b, 0.05) == 45

    def test_filter_inward(self):
        mask = self.pairs.add_mask_description('inwards', 'Mask read pairs that are inward '
                                                          'facing and closer than 100bp')
        in_filter = InwardPairsFilter(minimum_distance=100, mask=mask)

        assert len(self.pairs) == 44
        self.pairs.filter(in_filter)

        assert len(self.pairs) == 18

    def test_filter_outward(self):
        mask = self.pairs.add_mask_description('outwards', 'Mask read pairs that are outward '
                                                           'facing and closer than 100bp')
        out_filter = OutwardPairsFilter(minimum_distance=100, mask=mask)

        assert len(self.pairs) == 44
        self.pairs.filter(out_filter)
        assert len(self.pairs) == 28

    def test_filter_redist(self):
        mask = self.pairs.add_mask_description('re-dist',
                                               'Mask read pairs where one half maps more '
                                               'than 100bp away from both RE sites')
        re_filter = ReDistanceFilter(maximum_distance=300, mask=mask)

        assert len(self.pairs) == 44
        self.pairs.filter(re_filter)
        assert len(self.pairs) == 10

    def test_filter_self_ligated(self):
        mask = self.pairs.add_mask_description('self_ligated', 'Mask read pairs that represent self-ligated fragments')
        self_ligation_filter = SelfLigationFilter(mask=mask)

        assert len(self.pairs) == 44
        self.pairs.filter(self_ligation_filter)
        assert len(self.pairs) == 7

    def test_get_ligation_structure_biases(self):
        sam_file1 = os.path.join(self.dir, "test_matrix", "yeast.sample.chrI.1_sorted.sam")
        sam_file2 = os.path.join(self.dir, "test_matrix", "yeast.sample.chrI.2_sorted.sam")

        chrI = Chromosome.from_fasta(os.path.join(self.dir, "test_matrix", "chrI.fa"))
        genome = Genome(chromosomes=[chrI])
        pairs = self.pairs_class()
        regions = genome.get_regions('HindIII')
        pairs.add_regions(regions.regions)
        pair_generator = SamBamReadPairGenerator(sam_file1, sam_file2)
        pairs.add_read_pairs(pair_generator)
        genome.close()
        regions.close()

        x, i, o, b = pairs.get_ligation_structure_biases(sampling=200, skip_self_ligations=False)
        assert len(x) == len(i) == len(o) == len(b) == 3
        assert x.tolist() == [470, 4489, 19259]
        assert i.tolist() == [2.8756218905472637, 0.8059701492537313, 0.6368159203980099]
        assert o.tolist() == [0.2537313432835821, 0.24875621890547264, 0.46766169154228854]
        assert b.tolist() == [830, 413, 423]
        pairs.close()

    def test_re_dist(self):
        read1 = FragmentRead(GenomicRegion(chromosome='chr1', start=1, end=1000), position=200, strand=-1)
        assert read1.re_distance() == 199
        read2 = FragmentRead(GenomicRegion(chromosome='chr1', start=1, end=1000), position=990, strand=-1)
        assert read2.re_distance() == 10

    def test_iterate_exclude_filters(self):
        mask = self.pairs.add_mask_description('inwards', 'Mask read pairs that are inward '
                                                          'facing and closer than 100bp')
        in_filter = InwardPairsFilter(minimum_distance=100, mask=mask)
        self.pairs.filter(in_filter)
        mask = self.pairs.add_mask_description('outwards', 'Mask read pairs that are outward '
                                                           'facing and closer than 100bp')
        out_filter = OutwardPairsFilter(minimum_distance=100, mask=mask)
        self.pairs.filter(out_filter)
        mask = self.pairs.add_mask_description('re-dist',
                                               'Mask read pairs where one half maps more '
                                               'than 300bp away from both RE sites')
        re_filter = ReDistanceFilter(maximum_distance=300, mask=mask)
        self.pairs.filter(re_filter)
        assert len(self.pairs) == 0
        assert len(list(self.pairs.pairs(excluded_filters=['inwards', 'outwards', 're-dist']))) == 44
        assert len(list(self.pairs.pairs(excluded_filters=['inwards', 'outwards']))) == 10
        assert len(list(self.pairs.pairs(excluded_filters=['inwards', 're-dist']))) == 28
        assert len(list(self.pairs.pairs(excluded_filters=['outwards', 're-dist']))) == 18
        assert len(list(self.pairs.pairs(excluded_filters=['inwards']))) == 8
        assert len(list(self.pairs.pairs(excluded_filters=['outwards']))) == 2
        assert len(list(self.pairs.pairs(excluded_filters=['re-dist']))) == 2
        assert len(list(self.pairs.pairs(excluded_filters=[in_filter, re_filter]))) == 28
        assert len(list(self.pairs.pairs(excluded_filters=[in_filter, mask]))) == 28
        assert len(list(self.pairs.pairs(excluded_filters=[in_filter, 3]))) == 28


class TestFragmentRead:
    def setup_method(self, method):
        fragment1 = GenomicRegion(start=1, end=1000, chromosome='chr1', strand=1, ix=0)
        fragment2 = GenomicRegion(start=1001, end=2000, chromosome='chr2', strand=-1, ix=1)
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
                GenomicRegion(start=1, end=1000, chromosome='chr1'),
                position=500, strand=1
            ),
            FragmentRead(
                GenomicRegion(start=10001, end=11000, chromosome='chr1'),
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
                GenomicRegion(start=1, end=1000, chromosome='chr1'),
                position=500, strand=1
            ),
            FragmentRead(
                GenomicRegion(start=1, end=1000, chromosome='chr1'),
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
                GenomicRegion(start=1, end=1000, chromosome='chr1'),
                position=500, strand=-1
            ),
            FragmentRead(
                GenomicRegion(start=1, end=1000, chromosome='chr2'),
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
                GenomicRegion(start=1, end=1000, chromosome='chr1'),
                position=500, strand=-1
            ),
            FragmentRead(
                GenomicRegion(start=1001, end=2000, chromosome='chr1'),
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
        with ReadPairs(file_name=filename, mode='a', tmpdir='/tmp') as f:
            assert not os.path.isfile(filename)
            assert os.path.isfile(f.tmp_file_name)
        assert os.path.isfile(filename)
        assert not os.path.isfile(f.tmp_file_name)

    def test_tmp_with_exception(self, tmpdir):
        filename = str(tmpdir) + "/test.file"
        with pytest.raises(Exception):
            with ReadPairs(file_name=filename, mode='a', tmpdir='/tmp') as f:
                assert not os.path.isfile(filename)
                assert os.path.isfile(f.tmp_file_name)
                try:
                    raise Exception
                except Exception:
                    assert not os.path.isfile(filename)
                    assert not os.path.isfile(f.tmp_file_name)


