from __future__ import division
from kaic.data.network import RaoPeakCaller, process_matrix_range, RaoPeakInfo
from kaic.data.genomic import Hic, RegionMatrix, GenomicRegion
import numpy as np
import tables as t
import pickle
import os.path
import msgpack
from scipy.stats import poisson


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

        regions = [GenomicRegion(1, 1000, 'chr1', ix=0), GenomicRegion(1001, 2000, 'chr1', ix=1),
                   GenomicRegion(2001, 3000, 'chr1', ix=2), GenomicRegion(3001, 4000, 'chr1', ix=3),
                   GenomicRegion(4001, 5000, 'chr1', ix=4), GenomicRegion(5001, 6000, 'chr1', ix=5),
                   GenomicRegion(6001, 7000, 'chr1', ix=6)]
        self.regions_dict = {i: region for i, region in enumerate(regions)}
        self.m = RegionMatrix(np.array(l), col_regions=regions, row_regions=regions)

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
        def e(i, j):
            return 1
        assert RaoPeakCaller.e_ll(self.m, 3, 3, e, w=3, p=0) == 1
        assert RaoPeakCaller.e_ll(self.m, 3, 3, e, w=3, p=1) == 7/8
        assert RaoPeakCaller.e_ll(self.m, 3, 3, e, w=3, p=2) == 3/5
        assert RaoPeakCaller.e_ll(self.m, 3, 3, e, w=3, p=3) == 0

        # change w
        assert RaoPeakCaller.e_ll(self.m, 3, 3, e, w=2, p=1) == 4/3
        assert RaoPeakCaller.e_ll(self.m, 3, 3, e, w=1, p=1) == 0
        # change p and w
        assert RaoPeakCaller.e_ll(self.m, 3, 3, e, w=1, p=0) == 2

        # edge case
        assert RaoPeakCaller.e_ll(self.m, 2, 2, e, w=3, p=1) == 2/5

    def test_e_h(self):
        def e(i, j):
            return 1

        assert RaoPeakCaller.e_h(self.m, 3, 3, e, w=3, p=0) == (9+45)/(9+9)
        assert RaoPeakCaller.e_h(self.m, 3, 3, e, w=3, p=1) == (3+33)/(6+6)
        assert RaoPeakCaller.e_h(self.m, 3, 3, e, w=3, p=2) == (0+18)/(3+3)
        assert RaoPeakCaller.e_h(self.m, 3, 3, e, w=3, p=3) == 0

        # change w
        assert RaoPeakCaller.e_h(self.m, 3, 3, e, w=2, p=1) == (3+15)/6
        assert RaoPeakCaller.e_h(self.m, 3, 3, e, w=1, p=1) == 0
        # change p and w
        assert RaoPeakCaller.e_h(self.m, 3, 3, e, w=1, p=0) == (6+12)/6

    def test_e_v(self):
        def e(i, j):
            return 1

        assert RaoPeakCaller.e_v(self.m, 3, 3, e, w=3, p=0) == (27+27)/(9+9)
        assert RaoPeakCaller.e_v(self.m, 3, 3, e, w=3, p=1) == (18+18)/(6+6)
        assert RaoPeakCaller.e_v(self.m, 3, 3, e, w=3, p=2) == (9+9)/(3+3)
        assert RaoPeakCaller.e_v(self.m, 3, 3, e, w=3, p=3) == 0

        # change w
        assert RaoPeakCaller.e_v(self.m, 3, 3, e, w=2, p=1) == (18+18)/12
        assert RaoPeakCaller.e_v(self.m, 3, 3, e, w=1, p=1) == 0
        # change p and w
        assert RaoPeakCaller.e_v(self.m, 3, 3, e, w=1, p=0) == (9+9)/6

    def test_e_d(self):
        def e(i, j):
            return 1

        assert RaoPeakCaller.e_d(self.m, 3, 3, e, w=3, p=0) == (9+9+45+45)/(9*4)
        assert RaoPeakCaller.e_d(self.m, 3, 3, e, w=3, p=1) == (7+7+41+41)/(8*4)
        assert RaoPeakCaller.e_d(self.m, 3, 3, e, w=3, p=2) == (3+3+27+27)/(5*4)
        assert RaoPeakCaller.e_d(self.m, 3, 3, e, w=3, p=3) == 0

        # change w
        assert RaoPeakCaller.e_d(self.m, 3, 3, e, w=2, p=1) == (4+4+14+14)/(3*4)
        assert RaoPeakCaller.e_d(self.m, 3, 3, e, w=1, p=1) == 0
        # change p and w
        assert RaoPeakCaller.e_d(self.m, 3, 3, e, w=1, p=0) == (2+2+4+4)/4

    def test_lambda_chunks(self):
        assert np.array_equal(RaoPeakCaller._lambda_chunks(1000), [2**(0/3), 2**(1/3), 2**(2/3), 2**(3/3), 2**(4/3),
                                                                   2**(5/3), 2**(6/3), 2**(7/3), 2**(8/3), 2**(9/3),
                                                                   2**(10/3), 2**(11/3), 2**(12/3), 2**(13/3),
                                                                   2**(14/3), 2**(15/3), 2**(16/3), 2**(17/3),
                                                                   2**(18/3), 2**(19/3), 2**(20/3), 2**(21/3),
                                                                   2**(22/3), 2**(23/3), 2**(24/3), 2**(25/3),
                                                                   2**(26/3), 2**(27/3), 2**(28/3), 2**(29/3),
                                                                   2**(30/3)])

    def test_submatrix_indices(self):
        # top-left
        ij_pairs = [
            (0, 0),
            (0, 1),
            (0, 2),
            (1, 0),
            (1, 1)
        ]

        m_sub, ij_update = RaoPeakCaller._submatrix_indices(self.m, ij_pairs, w_max=3)

        assert np.array_equal(m_sub, [
            [0, 1, 2, 3, 4, 5],
            [0, 1, 2, 3, 4, 5],
            [0, 1, 2, 3, 4, 5],
            [0, 1, 2, 9, 4, 5],
            [0, 1, 2, 3, 4, 5]
        ])
        assert np.array_equal(ij_pairs, ij_update)

        # center
        ij_pairs = [
            (3, 3)
        ]

        m_sub, ij_update = RaoPeakCaller._submatrix_indices(self.m, ij_pairs, w_max=2)

        assert np.array_equal(m_sub, [
            [1, 2, 3, 4, 5],
            [1, 2, 3, 4, 5],
            [1, 2, 9, 4, 5],
            [1, 2, 3, 4, 5],
            [1, 2, 3, 4, 5],
        ])
        assert np.array_equal(ij_update, [
            (2,2)
        ])

        # bottom-right
        ij_pairs = [
            (5, 5),
            (5, 6),
            (6, 5),
            (6, 6)
        ]

        m_sub, ij_update = RaoPeakCaller._submatrix_indices(self.m, ij_pairs, w_max=2)

        assert np.array_equal(m_sub, [
            [9, 4, 5, 6],
            [3, 4, 5, 6],
            [3, 4, 5, 6],
            [3, 4, 5, 6]
        ])
        assert np.array_equal(ij_update, [
            (2, 2),
            (2, 3),
            (3, 2),
            (3, 3)
        ])

        # ij_pairs = [(3, 6), (4, 4), (4, 5), (4, 6), (5, 5), (5, 6), (6, 6)]
        #
        # m_sub, ij_update = RaoPeakCaller._submatrix_indices(self.m, ij_pairs, w_max=2)
        #
        # assert np.array_equal(m_sub, [
        #     [2 3 4 5 6],
        #     [2 3 4 5 6],
        #     [2 9 4 5 6],
        #     [2 3 4 5 6],
        #     [2 3 4 5 6],
        #     [2 3 4 5 6]
        # ])
        # assert np.array_equal(ij_update, [
        #     (2, 2),
        #     (2, 3),
        #     (3, 2),
        #     (3, 3)
        # ])


    def test_find_chunk(self):
        chunks = RaoPeakCaller._lambda_chunks(1000)
        assert RaoPeakCaller._find_chunk(chunks, None) is None
        assert RaoPeakCaller._find_chunk(chunks, 0) == 0
        assert RaoPeakCaller._find_chunk(chunks, 1) == 0
        assert RaoPeakCaller._find_chunk(chunks, 1.001) == 1
        assert RaoPeakCaller._find_chunk(chunks, 1.5) == 2
        assert RaoPeakCaller._find_chunk(chunks, 2) == 3
        assert RaoPeakCaller._find_chunk(chunks, 30) == 15
        assert RaoPeakCaller._find_chunk(chunks, 1024) == 30

    def test_process_matrix_range(self):
        ij_pairs = [
            (0, 0), (0, 1), (0, 2),
            (0, 3), (0, 4), (0, 5),
            (0, 6), (1, 0), (1, 1),
            (1, 2), (1, 3), (1, 4),
            (1, 5), (1, 6), (3, 0),
            (3, 1), (3, 2), (3, 3),
            (3, 4), (3, 5)
        ]

        def e(i, j):
            return 1

        min_locus_dist = 1
        #c = [2, 2, 2, 2, 2, 2, 2]
        c = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
        chunks = RaoPeakCaller._lambda_chunks(36)
        rv = process_matrix_range(self.m, msgpack.dumps(ij_pairs), msgpack.dumps(ij_pairs), 1, c,
                                  chunks, w=2, p=0, min_locus_dist=min_locus_dist,
                                  min_ll_reads=2, max_w=2)
        (region_list, observed_list,
         e_ll_list, e_h_list,
         e_v_list, e_d_list,
         observed_counts) = msgpack.loads(rv)

        print ij_pairs
        print region_list
        print observed_list

        for ix, pair in enumerate(region_list):
            print ix
            i = pair[0]
            j = pair[1]
            observed = self.m[i, j]

            if observed == 0:
                continue

            assert observed/(c[i]*c[j]) == observed_list[ix][0]
            e_ll, e_h, e_v, e_d = RaoPeakCaller.e_all(self.m, i, j, 1, w=2, p=0)
            if j-i <= min_locus_dist:
                assert e_ll_list[ix][0] is None
                assert e_h_list[ix][0] is None
                assert e_v_list[ix][0] is None
                assert e_d_list[ix][0] is None
            else:
                assert e_ll_list[ix][0] == RaoPeakCaller.e_ll(self.m, i, j, e, w=2, p=0)/(c[i]*c[j]) == e_ll/(c[i]*c[j])
                assert e_h_list[ix][0] == RaoPeakCaller.e_h(self.m, i, j, e, w=2, p=0)/(c[i]*c[j]) == e_h/(c[i]*c[j])
                assert e_v_list[ix][0] == RaoPeakCaller.e_v(self.m, i, j, e, w=2, p=0)/(c[i]*c[j]) == e_v/(c[i]*c[j])
                assert e_d_list[ix][0] == RaoPeakCaller.e_d(self.m, i, j, e, w=2, p=0)/(c[i]*c[j]) == e_d/(c[i]*c[j])

        # lower part of matrix
        ij_region_pairs = [
            (3, 3),
            (3, 4), (3, 5),
            (4, 4), (4, 5)
        ]

        m_sub, ij_updated_pairs = RaoPeakCaller._submatrix_indices(self.m, ij_region_pairs, w_max=2)

        def e(i, j):
            return 1

        min_locus_dist = 1
        c = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
        chunks = RaoPeakCaller._lambda_chunks(36)
        rv = process_matrix_range(m_sub, msgpack.dumps(ij_updated_pairs), msgpack.dumps(ij_region_pairs),
                                  [1, 1, 1, 1, 1, 1, 1], c,
                                  chunks, w=2, p=0, min_locus_dist=min_locus_dist,
                                  min_ll_reads=2, max_w=2)

        (region_list, observed_list,
         e_ll_list, e_h_list,
         e_v_list, e_d_list,
         observed_counts) = msgpack.loads(rv)

        for ix, pair in enumerate(region_list):
            i = pair[0]
            j = pair[1]
            observed = self.m[i, j]

            if observed == 0:
                continue

            assert observed/(c[i]*c[j]) == observed_list[ix][0]
            e_ll, e_h, e_v, e_d = RaoPeakCaller.e_all(self.m, i, j, 1, w=2, p=0)
            if j-i <= min_locus_dist:
                assert e_ll_list[ix][0] is None
                assert e_h_list[ix][0] is None
                assert e_v_list[ix][0] is None
                assert e_d_list[ix][0] is None
            else:
                assert e_ll_list[ix][0] == RaoPeakCaller.e_ll(self.m, i, j, e, w=2, p=0)/(c[i]*c[j]) == e_ll/(c[i]*c[j])
                assert e_h_list[ix][0] == RaoPeakCaller.e_h(self.m, i, j, e, w=2, p=0)/(c[i]*c[j]) == e_h/(c[i]*c[j])
                assert e_v_list[ix][0] == RaoPeakCaller.e_v(self.m, i, j, e, w=2, p=0)/(c[i]*c[j]) == e_v/(c[i]*c[j])
                assert e_d_list[ix][0] == RaoPeakCaller.e_d(self.m, i, j, e, w=2, p=0)/(c[i]*c[j]) == e_d/(c[i]*c[j])
            ix += 1

    def test_find_peaks_in_matrix(self):
        def e(i,j):
            return 1
        #c = [2, 2, 2, 2, 2, 2, 2]
        c = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
        mappable = [True, True, True, True, True, True, True]

        peaks = RaoPeakInfo(regions=self.m.row_regions)

        lambda_chunks = RaoPeakCaller._lambda_chunks(36)
        observed_chunk_distribution = RaoPeakCaller._get_chunk_distribution_container(lambda_chunks)

        peak_caller = RaoPeakCaller(max_w=2, min_ll_reads=2, min_locus_dist=1, batch_size=2, e_ll_cutoff=None,
                                    e_v_cutoff=None, e_d_cutoff=None, e_h_cutoff=None)
        peak_caller._find_peaks_in_matrix(self.m, 1, c, mappable, peaks,
                                          observed_chunk_distribution, lambda_chunks, w=2, p=0)

        assert len(peaks) == 5+4+3+2+1

        for peak in peaks.peaks():
            print peak
            i = peak.source
            j = peak.sink
            observed = self.m[i, j]
            assert observed/(c[i]*c[j]) == peak.observed
            e_ll, e_h, e_v, e_d = RaoPeakCaller.e_all(self.m, i, j, 1, w=2, p=0)
            if j-i <= 1:
                assert peak.e_ll is None
                assert peak.e_h is None
                assert peak.e_v is None
                assert peak.e_d is None
            else:
                assert abs(peak.e_h-RaoPeakCaller.e_h(self.m, i, j, e, w=2, p=0)/(c[i]*c[j])) < 0.01
                assert abs(peak.e_h-e_h/(c[i]*c[j])) < 0.01
                assert abs(peak.e_ll-RaoPeakCaller.e_ll(self.m, i, j, e, w=2, p=0)/(c[i]*c[j])) < 0.01
                assert abs(peak.e_ll-e_ll/(c[i]*c[j])) < 0.01
                assert abs(peak.e_d-RaoPeakCaller.e_d(self.m, i, j, e, w=2, p=0)/(c[i]*c[j])) < 0.01
                assert abs(peak.e_d-e_d/(c[i]*c[j])) < 0.01
                assert abs(peak.e_v-RaoPeakCaller.e_v(self.m, i, j, e, w=2, p=0)/(c[i]*c[j])) < 0.01
                assert abs(peak.e_v-e_v/(c[i]*c[j])) < 0.01

        def cmp_batches(peaks1, peaks2):
            peaks1_d = {(peak.source, peak.sink): peak for peak in peaks1}
            peaks2_d = {(peak.source, peak.sink): peak for peak in peaks2}

            for source, sink in peaks1_d.iterkeys():
                peak1 = peaks1_d[(source, sink)]
                peak2 = peaks2_d[(source, sink)]

                assert peak1.source == peak2.source
                assert peak1.sink == peak2.sink
                assert peak1.observed == peak2.observed
                assert peak1.e_ll == peak2.e_ll
                assert peak1.e_h == peak2.e_h
                assert peak1.e_v == peak2.e_v
                assert peak1.e_d == peak2.e_d
                assert peak1.e_ll_chunk == peak2.e_ll_chunk
                assert peak1.e_h_chunk == peak2.e_h_chunk
                assert peak1.e_v_chunk == peak2.e_v_chunk
                assert peak1.e_d_chunk == peak2.e_d_chunk

        peaks2 = RaoPeakInfo(regions=self.m.row_regions)
        peaks3 = RaoPeakInfo(regions=self.m.row_regions)

        peak_caller = RaoPeakCaller(max_w=2, min_ll_reads=2, min_locus_dist=1, batch_size=1, e_ll_cutoff=None,
                                    e_v_cutoff=None, e_d_cutoff=None, e_h_cutoff=None)
        lambda_chunks = RaoPeakCaller._lambda_chunks(36)
        observed_chunk_distribution = RaoPeakCaller._get_chunk_distribution_container(lambda_chunks)
        peak_caller._find_peaks_in_matrix(self.m, 1, c, mappable, peaks2,
                                          observed_chunk_distribution, lambda_chunks, w=2, p=0)
        cmp_batches(peaks, peaks2)

        peak_caller = RaoPeakCaller(max_w=2, min_ll_reads=2, min_locus_dist=1, batch_size=100, e_ll_cutoff=None,
                                    e_v_cutoff=None, e_d_cutoff=None, e_h_cutoff=None)
        lambda_chunks = RaoPeakCaller._lambda_chunks(36)
        observed_chunk_distribution = RaoPeakCaller._get_chunk_distribution_container(lambda_chunks)
        peak_caller._find_peaks_in_matrix(self.m, 1, c, mappable, peaks3,
                                          observed_chunk_distribution, lambda_chunks, w=2, p=0)
        cmp_batches(peaks, peaks3)
        peaks.close()
        peaks2.close()
        peaks3.close()

    def test_find_peaks_in_inter_matrix(self):
        def e(i, j):
            return 1

        c = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
        mappable = [True, True, True, True, True, True, True]

        peaks = RaoPeakInfo(regions=self.m.row_regions)
        peak_caller = RaoPeakCaller(max_w=2, min_ll_reads=2, min_locus_dist=1, batch_size=2, e_ll_cutoff=None,
                                    e_v_cutoff=None, e_d_cutoff=None, e_h_cutoff=None)
        peak_caller._find_peaks_in_matrix(self.m, 1, c, mappable, peaks,
                                          None, None, w=2, p=0)

        assert len(peaks) == 18

        seen = set()
        for peak in peaks.peaks():
            print peak
            i = peak.source
            j = peak.sink

            if (i, j) in seen:
                continue

            seen.add((i, j))

            observed = self.m[i, j]
            assert observed/(c[i]*c[j]) == peak.observed
            e_ll, e_h, e_v, e_d = RaoPeakCaller.e_all(self.m, i, j, 1, w=2, p=0)
            if j-i <= 1:
                assert peak.e_ll is None
                assert peak.e_h is None
                assert peak.e_v is None
                assert peak.e_d is None
            else:
                assert abs(peak.e_h-RaoPeakCaller.e_h(self.m, i, j, e, w=2, p=0)/(c[i]*c[j])) < 0.01
                assert abs(peak.e_h-e_h/(c[i]*c[j])) < 0.01
                assert abs(peak.e_ll-RaoPeakCaller.e_ll(self.m, i, j, e, w=2, p=0)/(c[i]*c[j])) < 0.01
                assert abs(peak.e_ll-e_ll/(c[i]*c[j])) < 0.01
                assert abs(peak.e_d-RaoPeakCaller.e_d(self.m, i, j, e, w=2, p=0)/(c[i]*c[j])) < 0.01
                assert abs(peak.e_d-e_d/(c[i]*c[j])) < 0.01
                assert abs(peak.e_v-RaoPeakCaller.e_v(self.m, i, j, e, w=2, p=0)/(c[i]*c[j])) < 0.01
                assert abs(peak.e_v-e_v/(c[i]*c[j])) < 0.01

        def cmp_batches(peaks1, peaks2):
            peaks1_d = {(peak.source, peak.sink): peak for peak in peaks1}
            peaks2_d = {(peak.source, peak.sink): peak for peak in peaks2}

            for source, sink in peaks1_d.iterkeys():
                peak1 = peaks1_d[(source, sink)]
                peak2 = peaks2_d[(source, sink)]

                assert peak1.source == peak2.source
                assert peak1.sink == peak2.sink
                assert peak1.observed == peak2.observed
                assert peak1.e_ll == peak2.e_ll
                assert peak1.e_h == peak2.e_h
                assert peak1.e_v == peak2.e_v
                assert peak1.e_d == peak2.e_d
                assert peak1.e_ll_chunk == peak2.e_ll_chunk
                assert peak1.e_h_chunk == peak2.e_h_chunk
                assert peak1.e_v_chunk == peak2.e_v_chunk
                assert peak1.e_d_chunk == peak2.e_d_chunk

        peaks2 = RaoPeakInfo(regions=self.m.row_regions)
        peaks3 = RaoPeakInfo(regions=self.m.row_regions)

        peak_caller = RaoPeakCaller(max_w=2, min_ll_reads=2, min_locus_dist=1, batch_size=1, e_ll_cutoff=None,
                                    e_v_cutoff=None, e_d_cutoff=None, e_h_cutoff=None)
        peak_caller._find_peaks_in_matrix(self.m, 1, c, mappable, peaks2,
                                          None, None, w=2, p=0)
        cmp_batches(peaks, peaks2)

        peak_caller = RaoPeakCaller(max_w=2, min_ll_reads=2, min_locus_dist=1, batch_size=100, e_ll_cutoff=None,
                                    e_v_cutoff=None, e_d_cutoff=None, e_h_cutoff=None)
        peak_caller._find_peaks_in_matrix(self.m, 1, c, mappable, peaks3,
                                          None, None, w=2, p=0)
        cmp_batches(peaks, peaks3)
        peaks.close()
        peaks2.close()
        peaks3.close()

    def test_fdr_cutoffs(self):
        lambda_chunks = RaoPeakCaller._lambda_chunks(4)
        observed_chunk_distribution = RaoPeakCaller._get_chunk_distribution_container(lambda_chunks)

        xy = [(0, 10), (1, 10), (2, 5), (3, 2), (4, 1), (5, 1), (6, 1), (7, 1), (8, 1), (9, 1)]
        for x, y in xy:
            observed_chunk_distribution['ll'][0][x] = y

        fdr_cutoffs = RaoPeakCaller._get_fdr_cutoffs(lambda_chunks, observed_chunk_distribution)
        for fdr in fdr_cutoffs['ll'][0].itervalues():
            assert 0 <= fdr <= 1

    def test_chromosome_map(self):
        dir = os.path.dirname(os.path.realpath(__file__))
        hic_10kb = Hic(dir + "/test_network/rao2014.chr11_77400000_78600000.hic", mode='r')
        chromosome_map = RaoPeakCaller.chromosome_map(hic_10kb)
        for i, region in enumerate(hic_10kb.regions(lazy=True)):
            if region.chromosome == 'chr11':
                assert chromosome_map[i] == 0
            else:
                print region
                assert chromosome_map[i] == 1
        hic_10kb.close()

    def test_call_peaks(self):
        dir = os.path.dirname(os.path.realpath(__file__))
        hic_10kb = Hic(dir + "/test_network/rao2014.chr11_77400000_78600000.hic", mode='r')

        peak_caller = RaoPeakCaller(process_inter=False, e_ll_cutoff=1.75,
                                    e_d_cutoff=1.75, e_h_cutoff=1.5, e_v_cutoff=1.5)
        peaks = peak_caller.call_peaks(hic_10kb)
        peak_info = peaks.peak_table

        assert len(peak_info) == 219

        valid_peaks = []

        has_43_57 = False
        for peak in peak_info:
            if peak['fdr_ll'] < 0.1 and peak['fdr_v'] < 0.1 and peak['fdr_h'] < 0.1 and peak['fdr_d'] < 0.1:
                valid_peaks.append(peak.fetch_all_fields())
            if peak['source'] == 43 and peak['sink'] == 57:
                has_43_57 = True

        assert len(valid_peaks) == 6
        assert has_43_57
        hic_10kb.close()
        peaks.close()

    def test_call_peaks_multiprocessing(self):
        dir = os.path.dirname(os.path.realpath(__file__))
        hic_10kb = Hic(dir + "/test_network/rao2014.chr11_77400000_78600000.hic", mode='r')

        peak_caller = RaoPeakCaller(process_inter=False, e_ll_cutoff=1.75,
                                    e_d_cutoff=1.75, e_h_cutoff=1.5, e_v_cutoff=1.5,
                                    cluster=False)
        peaks = peak_caller.call_peaks(hic_10kb)
        peak_info = peaks.peak_table

        assert len(peak_info) == 219

        valid_peaks = []

        has_43_57 = False
        for peak in peak_info:
            if peak['fdr_ll'] < 0.1 and peak['fdr_v'] < 0.1 and peak['fdr_h'] < 0.1 and peak['fdr_d'] < 0.1:
                valid_peaks.append(peak.fetch_all_fields())
            if peak['source'] == 43 and peak['sink'] == 57:
                has_43_57 = True

        assert len(valid_peaks) == 6
        assert has_43_57
        hic_10kb.close()
        peaks.close()

    def test_merge_peaks(self):
        directory = os.path.dirname(os.path.realpath(__file__))
        peaks = RaoPeakInfo(directory + "/test_network/rao2014.chr11_77400000_78600000.peaks", mode='r')

        merged_peaks = peaks.merged_peaks()
        assert len(merged_peaks) < len(peaks)
        assert len(merged_peaks) == 4
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
        for i in xrange(m_load.shape[0]):
            for j in xrange(m_load.shape[1]):
                assert self.m[i, j] == m_load[i, j]

        # row regions identical
        for region1, region2 in zip(self.m.row_regions, m_load.row_regions):
            assert region1.start == region2.start
            assert region1.end == region2.end

        # col regions identical
        for region1, region2 in zip(self.m.col_regions, m_load.col_regions):
            assert region1.start == region2.start
            assert region1.end == region2.end
