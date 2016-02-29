from __future__ import division
from kaic.data.genomic import Hic, Node, Edge
from kaic.architecture.architecture import calculateondemand
from kaic.architecture.hic_architecture import PossibleContacts, ExpectedContacts
from kaic.architecture.genome_architecture import MatrixArchitecturalRegionFeature, VectorArchitecturalRegionFeature
from kaic.data.genomic import GenomicRegion, Node, Edge
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

    def test_setitem(self):
        assert self.vaf[0, 'a'] == 1
        self.vaf[0, 'a'] = 9
        assert self.vaf[0, 'a'] == 9

        assert self.vaf[2, 'chromosome'] == 'chr1'
        self.vaf[2, 'chromosome'] = 'chr2'
        assert self.vaf[2, 'chromosome'] == 'chr2'
        self.vaf[2, 'chromosome'] = 'chr1'

        assert np.array_equal(['a', 'b', 'c'], self.vaf['chr1', 6])
        self.vaf['chr1', 6] = ['d', 'e', 'f']
        assert np.array_equal(['d', 'e', 'f'], self.vaf['chr1', 6])


class MAF(MatrixArchitecturalRegionFeature):
    """
    This only exists so we can instantiate a MatrixArchitecturalRegionFeature
    for testing.
    """
    def __init__(self, file_name=None, mode='a', data_fields=None,
                 regions=None, edges=None, tmpdir=None):
        MatrixArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, data_fields=data_fields,
                                                  regions=regions, edges=edges, tmpdir=tmpdir)

    def _calculate(self, *args, **kwargs):
        for i in xrange(10):
            if i < 5:
                chromosome = 'chr1'
                start = i*1000
                end = (i+1)*1000
            elif i < 8:
                chromosome = 'chr2'
                start = (i-5)*1000
                end = (i+1-5)*1000
            else:
                chromosome = 'chr3'
                start = (i-8)*1000
                end = (i+1-8)*1000
            node = Node(chromosome=chromosome, start=start, end=end)
            self.add_region(node, flush=False)
        self.flush()

        for i in xrange(10):
            for j in xrange(i, 10):
                edge = Edge(source=i, sink=j, weight=i*j, foo=i, bar=j, baz='x' + str(i*j))
                self.add_edge(edge, flush=False)
        self.flush()

    @calculateondemand
    def foo(self, key=None):
        return self.as_matrix(key, values_from='foo')


class TestMatrixArchitecturalRegionFeature:
    def setup_method(self, method):
        self.maf = MAF(data_fields={'foo': t.Int32Col(pos=0),
                                    'bar': t.Float32Col(pos=1),
                                    'baz': t.StringCol(50, pos=2)})

    def teardown_method(self, method):
        self.maf.close()

    def test_edges(self):
        assert len(self.maf.edges()) == 55


class TestPossbibleContacts:
    def setup_method(self, method):
        hic = Hic()

        # add some nodes (120 to be exact)
        nodes = []
        for i in range(1, 5000, 1000):
            nodes.append(Node(chromosome="chr1", start=i, end=i+1000-1))
        for i in range(1, 3000, 1000):
            nodes.append(Node(chromosome="chr2", start=i, end=i+1000-1))
        for i in range(1, 2000, 500):
            nodes.append(Node(chromosome="chr3", start=i, end=i+1000-1))
        hic.add_nodes(nodes)

        # add some edges with increasing weight for testing
        edges = []
        weight = 1
        for i in range(0, len(nodes)):
            for j in range(i, len(nodes)):
                edges.append(Edge(source=i, sink=j, weight=weight))
                weight += 1

        hic.add_edges(edges)

        self.hic = hic

    def teardown_method(self, method):
        self.hic.close()

    def test_no_region(self):
        print self.hic[:]
        with PossibleContacts(self.hic) as pc:

            assert pc.intra_possible() == 15+6+10
            assert pc.inter_possible() == 47

    def test_with_region(self):
        with PossibleContacts(self.hic, regions='chr1') as pc:

            assert pc.intra_possible() == 15
            assert pc.inter_possible() == 0

        with PossibleContacts(self.hic, regions=['chr2', 'chr3']) as pc:

            assert pc.intra_possible() == 16
            assert pc.inter_possible() == 12


class TestExpectedContacts:
    def setup_method(self, method):
        hic = Hic()

        # add some nodes (120 to be exact)
        nodes = []
        for i in range(1, 5000, 1000):
            nodes.append(Node(chromosome="chr1", start=i, end=i+1000-1))
        for i in range(1, 3000, 1000):
            nodes.append(Node(chromosome="chr2", start=i, end=i+1000-1))
        for i in range(1, 2000, 500):
            nodes.append(Node(chromosome="chr3", start=i, end=i+1000-1))
        hic.add_nodes(nodes)

        # add some edges with increasing weight for testing
        edges = []
        weight = 1
        for i in range(0, len(nodes)):
            for j in range(i, len(nodes)):
                edges.append(Edge(source=i, sink=j, weight=weight))
                weight += 1

        hic.add_edges(edges)

        self.hic = hic

    def teardown_method(self, method):
        self.hic.close()

    def test_intra_expected_no_smooting(self):
        with ExpectedContacts(self.hic, smooth=False) as ec:
            assert abs(ec.inter_expected() - (405+610+734)/47) < 0.001
            assert len(ec.intra_expected()) == 5
            assert abs(ec.intra_expected()[0] - (584/12)) < 0.001
            assert abs(ec.intra_expected()[1] - (408/9)) < 0.001
            assert abs(ec.intra_expected()[2] - (243/6)) < 0.001
            assert abs(ec.intra_expected()[3] - (92/3)) < 0.001
            assert abs(ec.intra_expected()[4] - 5) < 0.001

        with ExpectedContacts(self.hic, smooth=False, regions=['chr2', 'chr3']) as ec:
            assert abs(ec.inter_expected() - 734/12) < 0.001
            assert len(ec.intra_expected()) == 4
            assert abs(ec.intra_expected()[0] - (469/7)) < 0.001
            assert abs(ec.intra_expected()[1] - (332/5)) < 0.001
            assert abs(ec.intra_expected()[2] - (199/3)) < 0.001
            assert abs(ec.intra_expected()[3] - (72/1)) < 0.001

    def test_intra_expected_with_smooting(self):
        with ExpectedContacts(self.hic, smooth=True) as ec:
            assert abs(ec.inter_expected() - (405+610+734)/47) < 0.001
            assert len(ec.intra_expected()) == 5
            assert abs(ec.intra_expected()[0] - (584/12)) < 0.001
            assert abs(ec.intra_expected()[1] - (408/9)) < 0.001
            assert abs(ec.intra_expected()[2] - (243+408+92)/(9+6+3)) < 0.001
            assert abs(ec.intra_expected()[3] - (243++408+92+5)/(9+6+3+1)) < 0.001
            assert abs(ec.intra_expected()[4] - (243++408+92+5)/(9+6+3+1)) < 0.001

        with ExpectedContacts(self.hic, smooth=True, regions=['chr2', 'chr3']) as ec:
            assert abs(ec.inter_expected() - 734/12) < 0.001
            assert len(ec.intra_expected()) == 4
            assert abs(ec.intra_expected()[0] - (469/7)) < 0.001
            assert abs(ec.intra_expected()[1] - (332+469+199)/(7+5+3)) < 0.001
            assert abs(ec.intra_expected()[2] - (199+72+332)/(3+5+1)) < 0.001
            assert abs(ec.intra_expected()[3] - (199+72+332)/(3+5+1)) < 0.001
