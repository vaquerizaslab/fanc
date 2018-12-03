from __future__ import division
import kaic
from kaic.data.network import RaoPeakCaller, RaoPeakInfo
from kaic.data.genomic import Hic, RegionMatrix, GenomicRegion
from kaic.tools.general import pairwise
import numpy as np
import math
import pickle
import os.path


class TestRaoPeakCaller:
    def setup_method(self,  method):
        l = [
            [0, 1, 2, 3, 4, 5, 6],
            [0, 1, 2, 3, 4, 5, 6],
            [0, 1, 2, 3, 4, 5, 6],
            [0, 1, 2, 9, 4, 5, 6],
            [0, 1, 2, 3, 4, 5, 6],
            [0, 1, 2, 3, 4, 5, 6],
            [0, 1, 2, 3, 4, 5, 6]
        ]

        regions = [GenomicRegion(start=1, end=1000, chromosome='chr1', ix=0),
                   GenomicRegion(start=1001, end=2000, chromosome='chr1', ix=1),
                   GenomicRegion(start=2001, end=3000, chromosome='chr1', ix=2),
                   GenomicRegion(start=3001, end=4000, chromosome='chr1', ix=3),
                   GenomicRegion(start=4001, end=5000, chromosome='chr1', ix=4),
                   GenomicRegion(start=5001, end=6000, chromosome='chr1', ix=5),
                   GenomicRegion(start=6001, end=7000, chromosome='chr1', ix=6)]
        self.regions_dict = {i: region for i, region in enumerate(regions)}
        self.m = RegionMatrix(np.array(l), col_regions=regions, row_regions=regions)
        self.e = np.ones(self.m.shape)

    def test_ll_sum(self):
        # change p
        assert RaoPeakCaller.ll_sum(self.m, 3, 3, w=3, p=0) == 9
        assert RaoPeakCaller.ll_sum(self.m, 3, 3, w=3, p=1) == 7
        assert RaoPeakCaller.ll_sum(self.m, 3, 3, w=3, p=2) == 3
        assert RaoPeakCaller.ll_sum(self.m, 3, 3, w=3, p=3) == 0
        # change w
        assert RaoPeakCaller.ll_sum(self.m, 3, 3, w=4, p=1) == 7
        assert RaoPeakCaller.ll_sum(self.m, 3, 3, w=2, p=1) == 4
        assert RaoPeakCaller.ll_sum(self.m, 3, 3, w=1, p=1) == 0
        # change p and w
        assert RaoPeakCaller.ll_sum(self.m, 3, 3, w=1, p=0) == 2

        # edge case
        assert RaoPeakCaller.ll_sum(self.m, 2, 2, w=3, p=1) == 2

    def test_e_ll(self):
        assert RaoPeakCaller.e_ll(self.m, 3, 3, self.e, w=3, p=0) == 1
        assert RaoPeakCaller.e_ll(self.m, 3, 3, self.e, w=3, p=1) == 7/8
        assert RaoPeakCaller.e_ll(self.m, 3, 3, self.e, w=3, p=2) == 3/5
        assert np.isnan(RaoPeakCaller.e_ll(self.m, 3, 3, self.e, w=3, p=3))

        # change w
        assert RaoPeakCaller.e_ll(self.m, 3, 3, self.e, w=2, p=1) == 4/3
        assert np.isnan(RaoPeakCaller.e_ll(self.m, 3, 3, self.e, w=1, p=1))
        # change p and w
        assert RaoPeakCaller.e_ll(self.m, 3, 3, self.e, w=1, p=0) == 2

        # edge case
        assert RaoPeakCaller.e_ll(self.m, 2, 2, self.e, w=3, p=1) == 2/5

    def test_e_h(self):
        assert RaoPeakCaller.e_h(self.m, 3, 3, self.e, w=3, p=0) == (9+45)/(9+9)
        assert RaoPeakCaller.e_h(self.m, 3, 3, self.e, w=3, p=1) == (3+33)/(6+6)
        assert RaoPeakCaller.e_h(self.m, 3, 3, self.e, w=3, p=2) == (0+18)/(3+3)
        assert np.isnan(RaoPeakCaller.e_h(self.m, 3, 3, self.e, w=3, p=3))

        # change w
        assert RaoPeakCaller.e_h(self.m, 3, 3, self.e, w=2, p=1) == (3+15)/6
        assert np.isnan(RaoPeakCaller.e_h(self.m, 3, 3, self.e, w=1, p=1))
        # change p and w
        assert RaoPeakCaller.e_h(self.m, 3, 3, self.e, w=1, p=0) == (6+12)/6

    def test_e_v(self):
        assert RaoPeakCaller.e_v(self.m, 3, 3, self.e, w=3, p=0) == (27+27)/(9+9)
        assert RaoPeakCaller.e_v(self.m, 3, 3, self.e, w=3, p=1) == (18+18)/(6+6)
        assert RaoPeakCaller.e_v(self.m, 3, 3, self.e, w=3, p=2) == (9+9)/(3+3)
        assert np.isnan(RaoPeakCaller.e_v(self.m, 3, 3, self.e, w=3, p=3))

        # change w
        assert RaoPeakCaller.e_v(self.m, 3, 3, self.e, w=2, p=1) == (18+18)/12
        assert np.isnan(RaoPeakCaller.e_v(self.m, 3, 3, self.e, w=1, p=1))
        # change p and w
        assert RaoPeakCaller.e_v(self.m, 3, 3, self.e, w=1, p=0) == (9+9)/6

    def test_e_d(self):
        assert RaoPeakCaller.e_d(self.m, 3, 3, self.e, w=3, p=0) == (9+9+45+45)/(9*4)
        assert RaoPeakCaller.e_d(self.m, 3, 3, self.e, w=3, p=1) == (7+7+41+41)/(8*4)
        assert RaoPeakCaller.e_d(self.m, 3, 3, self.e, w=3, p=2) == (3+3+27+27)/(5*4)
        assert np.isnan(RaoPeakCaller.e_d(self.m, 3, 3, self.e, w=3, p=3))

        # change w
        assert RaoPeakCaller.e_d(self.m, 3, 3, self.e, w=2, p=1) == (4+4+14+14)/(3*4)
        assert np.isnan(RaoPeakCaller.e_d(self.m, 3, 3, self.e, w=1, p=1))
        # change p and w
        assert RaoPeakCaller.e_d(self.m, 3, 3, self.e, w=1, p=0) == (2+2+4+4)/4

    def test_find_chunk(self):
        assert RaoPeakCaller.find_chunk(None) is None
        assert RaoPeakCaller.find_chunk(0) == 0
        assert RaoPeakCaller.find_chunk(1) == 1
        assert RaoPeakCaller.find_chunk(1.001) == 1
        assert RaoPeakCaller.find_chunk(1.5) == 2
        assert RaoPeakCaller.find_chunk(1.7) == 3
        assert RaoPeakCaller.find_chunk(30) == 15
        assert RaoPeakCaller.find_chunk(1024) == 31

    def test_call_peaks(self):
        dir = os.path.dirname(os.path.realpath(__file__))
        hic_10kb = Hic(dir + "/test_network/rao2014.chr11_77400000_78600000.hic", mode='r')

        peak_caller = RaoPeakCaller()
        peaks = peak_caller.call_peaks(hic_10kb)

        assert len(peaks) == 6408

        valid_peaks = []

        has_43_57 = False
        for peak in peaks.edges:
            if peak.fdr_ll < 0.1 and peak.fdr_v < 0.1 and peak.fdr_h < 0.1 and peak.fdr_d < 0.1:
                valid_peaks.append(peak)
            if peak.source == 43 and peak.sink == 57:
                has_43_57 = True

        assert len(valid_peaks) == 124
        assert has_43_57
        hic_10kb.close()
        peaks.close()

    def test_merge_peaks(self):
        directory = os.path.dirname(os.path.realpath(__file__))
        peaks = RaoPeakInfo(directory + "/test_network/rao2014.chr11_77400000_78600000.peaks_filtered", mode='r')

        merged_peaks = peaks.merged_peaks()
        assert len(merged_peaks) < len(peaks)
        assert len(merged_peaks) == 24
        peaks.close()
        merged_peaks.close()

    def test_pickle_matrix(self):
        assert hasattr(self.m, 'row_regions')
        assert hasattr(self.m, 'col_regions')
        m_dump = pickle.dumps(self.m)
        m_load = pickle.loads(m_dump)
        assert hasattr(m_load, 'row_regions')
        assert hasattr(m_load, 'col_regions')

        # matrix identical
        for i in range(m_load.shape[0]):
            for j in range(m_load.shape[1]):
                assert self.m[i, j] == m_load[i, j]

        # row regions identical
        for region1, region2 in zip(self.m.row_regions, m_load.row_regions):
            assert region1.start == region2.start
            assert region1.end == region2.end

        # col regions identical
        for region1, region2 in zip(self.m.col_regions, m_load.col_regions):
            assert region1.start == region2.start
            assert region1.end == region2.end


class TestOverlapPeaks:
    def setup_method(self, method):
        peaks = [
            [(10, 12), (1, 4), (15, 8)],
            [(11, 12), (5, 8), (15.2, 7.8)],
            [(10.5, 11.8), (5.5, 8.1)],
        ]
        regions = [GenomicRegion('chr1', a + 1, b, ix=i) for i, (a, b) in enumerate(pairwise(np.arange(0, 10001, 100)))]
        self.peaks = {}
        for i in range(3):
            p = kaic.data.network.PeakInfo()
            p.add_regions(regions)
            p.add_edges(kaic.data.network.Peak(x=x, y=y, source=math.floor(x), sink=math.floor(y)) for x, y in (tuple(sorted(xy)) for xy in peaks[i]))
            self.peaks[i] = p

    def test_overlap(self):
        stats, merged = kaic.data.network.overlap_peaks(self.peaks, max_distance=100)
        expected = {
            (0, 1, 2): 1,
            (0,): 1,
            (1, 2): 1,
            (0, 1): 1,
        }
        assert all(len(merged[frozenset(x)]) == n for x, n in expected.items())
        assert all(stats["n"] == 1)
