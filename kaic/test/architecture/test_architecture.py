from kaic.architecture.architecture import VectorArchitecturalRegionFeature
from kaic.data.genomic import GenomicRegion
import tables as t
import types
import numpy as np


class TestVectorArchitecturalRegionFeature:
    def setup_method(self, method):
        self.vaf = VectorArchitecturalRegionFeature(data_fields={'a': t.Int32Col(), 'b': t.StringCol(10)})
        self.vaf.add_regions([GenomicRegion(1,1000,'chr1', a=1, b='a'),
                              GenomicRegion(1001, 2000, 'chr1', a=2, b='b'),
                              GenomicRegion(2001, 3000, 'chr1', a=3, b='c'),
                              GenomicRegion(1, 1000, 'chr2', a=4, b='d')])

    def test_get_rows(self):
        assert isinstance(self.vaf[0], GenomicRegion)
        assert self.vaf[0].a == 1
        assert self.vaf[0].b == 'a'

        regions = self.vaf[1:3]
        assert isinstance(regions, types.GeneratorType)
        for i, r in enumerate(regions):
            assert r.chromosome == 'chr1'
            assert r.a == i+2
            assert r.b == 'abcd'[i+1]

        regions = self.vaf['chr1']
        assert isinstance(regions, types.GeneratorType)
        for i, r in enumerate(regions):
            assert r.chromosome == 'chr1'
            assert r.a == i+1
            assert r.b == 'abcd'[i]

        regions = self.vaf['chr1:1-2000']
        assert isinstance(regions, types.GeneratorType)
        for i, r in enumerate(regions):
            assert r.chromosome == 'chr1'
            assert r.a == i+1
            assert r.b == 'abcd'[i]
        regions = self.vaf['chr1:1-2000']
        assert len(list(regions)) == 2

        regions = self.vaf['chr1:1-2001']
        assert isinstance(regions, types.GeneratorType)
        for i, r in enumerate(regions):
            assert r.chromosome == 'chr1'
            assert r.a == i+1
            assert r.b == 'abcd'[i]
        regions = self.vaf['chr1:1-2001']
        assert len(list(regions)) == 3

        regions = self.vaf[GenomicRegion(1, 2000, None)]
        assert isinstance(regions, types.GeneratorType)
        for i, r in enumerate(regions):
            if i < 2:
                assert r.chromosome == 'chr1'
                assert r.a == i+1
                assert r.b == 'abcd'[i]
            else:
                assert r.chromosome == 'chr2'
                assert r.a == 4
                assert r.b == 'd'

        regions = self.vaf[GenomicRegion(1, 2000, None)]
        assert len(list(regions)) == 3

        regions = self.vaf[GenomicRegion(None, None, None)]
        assert len(list(regions)) == 4

    def test_get_columns(self):
        # let's test single ones
        results = self.vaf[0, 'a']
        assert results == 1
        results = self.vaf[1, 'b']
        assert results == 'b'
        results = self.vaf[2, 'chromosome']
        assert results == 'chr1'

        # int
        results = self.vaf['chr1', 1]  # chromosome
        assert isinstance(results, list)
        assert np.array_equal(['chr1', 'chr1', 'chr1'], results)

        results = self.vaf['chr1', 6]  # b
        assert isinstance(results, list)
        assert np.array_equal(['a', 'b', 'c'], results)

        # str
        results = self.vaf['chr1', 'chromosome']  # chromosome
        assert isinstance(results, list)
        assert np.array_equal(['chr1', 'chr1', 'chr1'], results)

        results = self.vaf['chr1', 'b']  # b
        assert isinstance(results, list)
        assert np.array_equal(['a', 'b', 'c'], results)

        # slice
        results = self.vaf['chr1', 5:7]  # a, b
        assert isinstance(results, dict)
        assert 'a' in results
        assert 'b' in results
        assert np.array_equal([1, 2, 3], results['a'])
        assert np.array_equal(['a', 'b', 'c'], results['b'])

        # list
        results = self.vaf['chr1', ['a', 'b']]  # a, b
        assert isinstance(results, dict)
        assert 'a' in results
        assert 'b' in results
        assert np.array_equal([1, 2, 3], results['a'])
        assert np.array_equal(['a', 'b', 'c'], results['b'])




