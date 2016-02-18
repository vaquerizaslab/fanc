from kaic.architecture.architecture import TableArchitecturalFeature, calculateondemand
from kaic.data.genomic import GenomicRegion, VectorArchitecturalRegionFeature
import tables as t
import types
import numpy as np


class VAF(VectorArchitecturalRegionFeature):
    """
    This only exists so we can instantiate a VectorArchitecturalRegionFeature
    for testing.
    """
    def __init__(self, file_name=None, mode='a', data_fields=None,
                 regions=None, data=None, _table_name_data='region_data',
                 tmpdir=None):
        VectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, data_fields=data_fields,
                                                  regions=regions, data=data, _table_name_data=_table_name_data,
                                                  tmpdir=tmpdir)

    def _calculate(self, *args, **kwargs):
        self.add_regions([GenomicRegion(1, 1000, 'chr1', a=1, b='a'),
                          GenomicRegion(1001, 2000, 'chr1', a=2, b='b'),
                          GenomicRegion(2001, 3000, 'chr1', a=3, b='c'),
                          GenomicRegion(1, 1000, 'chr2', a=4, b='d')])


class TestVectorArchitecturalRegionFeature:
    def setup_method(self, method):
        self.vaf = VAF(data_fields={'a': t.Int32Col(), 'b': t.StringCol(10)})

    def teardown_method(self, method):
        self.vaf.close()

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


class TAF(TableArchitecturalFeature):
    """
    This only exists so we can instantiate a VectorArchitecturalRegionFeature
    for testing.
    """
    def __init__(self, group, fields, file_name=None, mode='a', tmpdir=None):
        TableArchitecturalFeature.__init__(self, group, fields,
                                           file_name=file_name, mode=mode, tmpdir=tmpdir)

    def _calculate(self, *args, **kwargs):
        self.data('a', [1, 2, 3, 4, 5])
        self.data('b', ['a', 'b', 'c', 'd', 'e'])

    @calculateondemand
    def a(self, key=None):
        if key is None:
            return self[:, 'a']
        return self[key, 'a']


class TestTableArchitecturalFeature:
    def setup_method(self, method):
        self.taf = TAF('test', {'a': t.Int32Col(), 'b': t.StringCol(10)})

    def teardown_method(self, method):
        self.taf.close()

    def test_get_rows(self):
        assert isinstance(self.taf[0], dict)
        assert self.taf[0]['a'] == 1
        assert self.taf[0]['b'] == 'a'

        results = self.taf[1:3]
        assert isinstance(results, types.GeneratorType)
        for i, r in enumerate(results):
            assert r['a'] == i+2
            assert r['b'] == 'abcd'[i+1]

    def test_get_columns(self):
        # let's test single ones
        results = self.taf[0, 'a']
        assert results == 1
        results = self.taf[1, 'b']
        assert results == 'b'

        # int
        results = self.taf[0, 0]
        assert results == 1
        results = self.taf[1, 1]
        assert results == 'b'

        # slice
        results = self.taf[0, 0:2]
        assert results == {'a': 1, 'b': 'a'}
        results = self.taf[1, 0:2]
        assert results == {'a': 2, 'b': 'b'}

    def test_calculateondemand_decorator(self):
        assert self.taf.a() == [1, 2, 3, 4, 5]
        assert self.taf.a(1) == 2
