import os
import numpy as np
from kaic.compatibility.cooler import to_cooler
from genomic_regions import GenomicRegion
from kaic.matrix import Edge, RegionPairsTable, RegionMatrixTable, RegionMatrix
from kaic.hic import Hic, _get_overlap_map, _edge_overlap_split_rao, kr_balancing, ice_balancing
from kaic.regions import Chromosome, Genome
from kaic.pairs import ReadPairs, SamBamReadPairGenerator
from kaic.tools.matrix import is_symmetric
import tables
import pytest

test_dir = os.path.dirname(os.path.realpath(__file__))


def _get_test_regions(folder=test_dir, with_bias=False):
    biases_file = os.path.join(folder, 'test_matrix', 'test_biases.txt')

    biases = []
    with open(biases_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line == '':
                continue

            b = 1/float(line)

            biases.append(b)

    if not with_bias:
        biases = [1.0] * len(biases)

    regions_file = os.path.join(test_dir, 'test_matrix', 'test_regions.bed')

    regions = []
    with open(regions_file, 'r') as f:
        for i, line in enumerate(f):
            line = line.rstrip()
            if line == '':
                continue

            fields = line.split("\t")

            chromosome = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            region = GenomicRegion(chromosome=chromosome,
                                   start=start,
                                   end=end,
                                   bias=biases[i],
                                   ix=i)
            regions.append(region)
    return regions


def _get_test_edges(folder=test_dir, norm=False):
    if norm:
        edges_file = os.path.join(folder, 'test_matrix', 'test_edges_kr.txt')
    else:
        edges_file = os.path.join(folder, 'test_matrix', 'test_edges_uncorrected.txt')

    edges = []
    with open(edges_file, 'r') as f:
        for i, line in enumerate(f):
            line = line.rstrip()
            if line == '':
                continue

            fields = line.split("\t")
            source = int(fields[0])
            sink = int(fields[1])
            weight = float(fields[2])

            if not np.isnan(weight):
                edge = Edge(source=source, sink=sink, weight=weight)
                edges.append(edge)
    return edges


class RegionMatrixContainerTestFactory:
    def setup_method(self, method):
        self.matrix = None

    def test_edges_iter(self):
        pass

    @pytest.mark.parametrize("lazy", [True, False])
    def test_get_edges_uncorrected(self, lazy):
        edges_dict = {(e.source, e.sink): e.weight for e in _get_test_edges(norm=False)}

        for edge in self.matrix.edges(norm=False):
            assert np.isclose(edge.weight,
                              edges_dict[(edge.source, edge.sink)],
                              rtol=1e-03)

    @pytest.mark.parametrize("lazy", [True, False])
    def test_get_edges(self, lazy):
        edges_dict = {(e.source, e.sink): e.weight for e in _get_test_edges(norm=True)}

        for edge in self.matrix.edges(norm=True, lazy=lazy):
            if (edge.source, edge.sink) not in edges_dict:
                continue

            assert np.isclose(edge.weight,
                              edges_dict[(edge.source, edge.sink)],
                              rtol=1e-03)


class TestHic(RegionMatrixContainerTestFactory):
    def setup_method(self, method):
        hic_file = os.path.join(test_dir, 'test_matrix', 'test_kaic.hic')
        self.matrix = Hic(hic_file, mode='r')

    def teardown_method(self, method):
        self.matrix.close()


class TestRegionPairs:
    def setup_method(self, method):
        self.rmt = RegionPairsTable(additional_edge_fields={'weight': tables.Int32Col(pos=0),
                                                            'foo': tables.Int32Col(pos=1),
                                                            'bar': tables.Float32Col(pos=2),
                                                            'baz': tables.StringCol(50, pos=3)})

        for i in range(10):
            if i < 5:
                chromosome = 'chr1'
                start = i * 1000
                end = (i + 1) * 1000
            elif i < 8:
                chromosome = 'chr2'
                start = (i - 5) * 1000
                end = (i + 1 - 5) * 1000
            else:
                chromosome = 'chr3'
                start = (i - 8) * 1000
                end = (i + 1 - 8) * 1000
            node = GenomicRegion(chromosome=chromosome, start=start, end=end)
            self.rmt.add_region(node)
        self.rmt.flush()

        for i in range(10):
            for j in range(i, 10):
                edge = Edge(source=i, sink=j, weight=i * j, foo=i, bar=j, baz='x' + str(i * j))
                self.rmt.add_edge(edge)
        self.rmt.flush()

        self.rp_class = RegionPairsTable

    def teardown_method(self, method):
        self.rmt.close()

    def test_create(self):
        rmt1 = RegionPairsTable()
        assert rmt1.field_names == ['source', 'sink']
        assert len(list(rmt1.edges())) == 0
        rmt1.close()

        rmt2 = RegionPairsTable(additional_edge_fields={'foo': tables.Int32Col(pos=0),
                                                        'bar': tables.Float32Col(pos=1)})
        assert rmt2.field_names == ['source', 'sink', 'foo', 'bar']
        assert len(list(rmt2.edges())) == 0
        rmt2.close()

        additional_edge_fields = {
            'foo': tables.Int32Col(pos=0),
            'bar': tables.Float32Col(pos=1)
        }

        rmt3 = RegionPairsTable(additional_edge_fields=additional_edge_fields)
        assert rmt3.field_names == ['source', 'sink', 'foo', 'bar']
        assert len(list(rmt3.edges())) == 0
        rmt3.close()

        assert len(list(self.rmt.edges())) == 55
        assert self.rmt.field_names == ['source', 'sink', 'weight', 'foo', 'bar', 'baz']

    def test_edges(self):

        for edge in self.rmt.edges():
            i = edge.source
            j = edge.sink
            assert edge.weight == i * j
            assert edge.foo == i
            assert edge.bar == j
            assert edge.baz == 'x' + str(i * j)

            with pytest.raises(AttributeError):
                assert edge.qux is None

    def test_edges_nodup(self):
        covered = set()
        for edge in self.rmt.edges((slice(0, 2), slice(1, 3))):
            pair = (edge.source, edge.sink)
            if pair in covered:
                assert 0
            covered.add(pair)

        covered = set()
        for edge in self.rmt.edges((slice(1, 3), slice(0, 2))):
            pair = (edge.source, edge.sink)
            if pair in covered:
                assert 0
            covered.add(pair)

        covered = set()
        for edge in self.rmt.edges((slice(0, 3), slice(1, 2))):
            pair = (edge.source, edge.sink)
            if pair in covered:
                assert 0
            covered.add(pair)

    def test_lazy_edges(self):
        for edge in self.rmt.edges(lazy=True):
            i = edge.source
            j = edge.sink
            assert edge.weight == i * j
            assert edge.foo == i
            assert edge.bar == j
            assert edge.baz.decode('utf-8') == 'x' + str(i * j)

            with pytest.raises(AttributeError):
                assert edge.qux is None

    def test_edges_set_attribute(self):
        for edge in self.rmt.edges():
            edge.foo = 999

        for edge in self.rmt.edges():
            assert edge.foo != 999

    def test_lazy_edges_set_attribute(self):
        for edge in self.rmt.edges(lazy=True):
            edge.foo = 999
            edge.update()

        for edge in self.rmt.edges():
            assert edge.foo == 999

    def test_edge_subset(self):
        edges = self.rmt.edge_subset(key=('chr2', 'chr2'))
        for edge in edges:
            assert edge.bar == max(edge.sink, edge.source)

        edges = self.rmt.edge_subset(key=('chr2', 'chr3'))
        for edge in edges:
            assert edge.bar == max(edge.sink, edge.source)

        edges = self.rmt.edge_subset(key=slice(0, None, None))
        s = 0
        for edge in edges:
            s += 1
            assert edge.bar == max(edge.sink, edge.source)
        assert s == 55

    def test_add_edge(self):
        rmt = self.rp_class(additional_edge_fields={'weight': tables.Float64Col()})
        rmt.add_region(GenomicRegion(chromosome='1', start=1, end=1000))
        rmt.flush()
        rmt.add_edge(Edge(0, 0, weight=100))
        rmt.flush()
        edge = rmt.edges[0]
        assert edge.source == 0
        assert edge.sink == 0
        assert edge.weight == 100
        rmt.close()

        rmt = self.rp_class(additional_edge_fields={'weight': tables.Float64Col()})
        rmt.add_region(GenomicRegion(chromosome='1', start=1, end=1000))
        rmt.flush()
        rmt.add_edge([0, 0, 100])
        rmt.flush()
        edge = rmt.edges[0]
        assert edge.source == 0
        assert edge.sink == 0
        assert edge.weight == 100
        rmt.close()

        rmt = self.rp_class(additional_edge_fields={'weight': tables.Float64Col()})
        rmt.add_region(GenomicRegion(chromosome='1', start=1, end=1000))
        rmt.flush()
        rmt.add_edge({'source': 0, 'sink': 0, 'weight': 100})
        rmt.flush()
        edge = rmt.edges[0]
        assert edge.source == 0
        assert edge.sink == 0
        assert edge.weight == 100
        rmt.close()


class TestRegionMatrixTable:
    def setup_method(self, method):
        self.rmt = RegionMatrixTable(additional_edge_fields={'weight': tables.Int32Col(pos=0),
                                                             'foo': tables.Int32Col(pos=1),
                                                             'bar': tables.Float32Col(pos=2),
                                                             'baz': tables.StringCol(50, pos=3)})

        for i in range(10):
            if i < 5:
                chromosome = 'chr1'
                start = i * 1000
                end = (i + 1) * 1000
            elif i < 8:
                chromosome = 'chr2'
                start = (i - 5) * 1000
                end = (i + 1 - 5) * 1000
            else:
                chromosome = 'chr3'
                start = (i - 8) * 1000
                end = (i + 1 - 8) * 1000
            node = GenomicRegion(chromosome=chromosome, start=start, end=end)
            self.rmt.add_region(node)
        self.rmt.flush()

        for i in range(10):
            for j in range(i, 10):
                edge = Edge(source=i, sink=j, weight=i * j, foo=i, bar=j, baz='x' + str(i * j))
                self.rmt.add_edge(edge)
        self.rmt.flush()

    def teardown_method(self, method):
        self.rmt.close()

    def test_matrix(self):
        m = self.rmt.matrix()
        for row_region in m.row_regions:
            i = row_region.ix
            for col_region in m.col_regions:
                j = col_region.ix
                if np.ma.is_masked(m[i, j]):
                    assert np.ma.is_masked(m[j, i])
                else:
                    assert m[i, j] == m[j, i] == i * j

        m = self.rmt.matrix(score_field='foo')
        for row_region in m.row_regions:
            i = row_region.ix
            for col_region in m.col_regions:
                j = col_region.ix
                if np.ma.is_masked(m[i, j]):
                    assert np.ma.is_masked(m[j, i])
                else:
                    assert m[i, j] == m[j, i] == min(i, j)

        m = self.rmt.matrix(score_field='bar')
        for row_region in m.row_regions:
            i = row_region.ix
            for col_region in m.col_regions:
                j = col_region.ix
                if np.ma.is_masked(m[i, j]):
                    assert np.ma.is_masked(m[j, i])
                else:
                    assert m[i, j] == m[j, i] == max(i, j)

    def test_matrix_subset(self):
        m = self.rmt.matrix(key=('chr2', 'chr2'), score_field='bar')
        for i, row_region in enumerate(m.row_regions):
            for j, col_region in enumerate(m.col_regions):
                assert (np.ma.is_masked(m[i, j]) and np.ma.is_masked(m[j, i])) or \
                       (m[i, j] == m[j, i] == max(row_region.ix, col_region.ix))

        m = self.rmt.matrix(key=('chr2', 'chr3'), score_field='bar', mask=False)

        for i, row_region in enumerate(m.row_regions):
            for j, col_region in enumerate(m.col_regions):
                print(i, j, row_region.ix, col_region.ix, m[i, j])
                assert m[i, j] == max(row_region.ix, col_region.ix)

        m = self.rmt.matrix(key=('chr3', 'chr2'), score_field='bar', mask=False)
        for i, row_region in enumerate(m.row_regions):
            for j, col_region in enumerate(m.col_regions):
                assert m[i, j] == max(row_region.ix, col_region.ix)


class TestHicBasic:
    def setup_method(self, method):
        self.dir = os.path.dirname(os.path.realpath(__file__))

        hic = Hic()

        # add some nodes (120 to be exact)
        nodes = []
        for i in range(1, 5000, 1000):
            nodes.append(GenomicRegion(chromosome="chr1", start=i, end=i + 1000 - 1))
        for i in range(1, 3000, 1000):
            nodes.append(GenomicRegion(chromosome="chr2", start=i, end=i + 1000 - 1))
        for i in range(1, 2000, 500):
            nodes.append(GenomicRegion(chromosome="chr3", start=i, end=i + 1000 - 1))
        hic.add_regions(nodes)

        # add some edges with increasing weight for testing
        edges = []
        weight = 1
        for i in range(0, len(nodes)):
            for j in range(i, len(nodes)):
                edges.append(Edge(source=i, sink=j, weight=weight))
                weight += 1

        hic.add_edges(edges)

        self.hic = hic
        self.hic_cerevisiae = Hic(self.dir + "/test_matrix/cerevisiae.chrI.HindIII_upgrade.hic")
        self.hic_class = Hic

    def teardown_method(self, method):
        self.hic_cerevisiae.close()
        self.hic.close()

    def test_initialize_empty(self):
        hic = self.hic_class()
        nodes = list(hic.regions())
        edges = list(hic.edges())
        assert len(nodes) == 0
        assert len(edges) == 0
        hic.close()

    def test_save_and_load(self, tmpdir):
        dest_file = str(tmpdir) + "/hic.h5"

        hic1 = self.hic_class(file_name=dest_file, mode='w')
        hic1.add_region(GenomicRegion('chr1', 1, 1000))
        hic1.add_region(GenomicRegion('chr2', 1, 1000))
        hic1.flush()
        hic1.add_edge([0, 1])
        hic1.flush()
        hic1.close()

        hic2 = self.hic_class(dest_file, mode='r')
        nodes2 = list(hic2.regions())
        edges2 = list(hic2.edges())
        assert len(nodes2) == 2
        assert len(edges2) == 1

        hic2.close()

    def test_nodes(self):
        nodes = list(self.hic.regions())
        assert len(nodes) == 12

    def test_edges(self):
        edges = list(self.hic.edges())
        assert len(edges) == 78

    def test_intra_edges(self):
        edges = list(self.hic.edges(inter_chromosomal=False))
        assert sum(1 for _ in edges) == 31

    def test_get_node_x_by_region(self):
        region1 = GenomicRegion.from_string('chr1')
        nodes1 = list(self.hic.regions(region1))
        assert len(nodes1) == 5

        region2 = GenomicRegion.from_string('chr2')
        nodes2 = list(self.hic.regions(region2))
        assert len(nodes2) == 3

        region3 = GenomicRegion.from_string('chr3')
        nodes3 = list(self.hic.regions(region3))
        assert len(nodes3) == 4

        region4 = GenomicRegion.from_string('chr1:3452-6000')
        nodes4 = list(self.hic.regions(region4))
        assert len(nodes4) == 2

        region5 = GenomicRegion.from_string('chr1:1-51000')
        nodes5 = list(self.hic.regions(region5))
        assert len(nodes5) == 5

    def test_getitem_nodes(self):
        # all
        node_ix1 = [r.ix for r in self.hic.regions(slice(None, None, None))]
        assert np.array_equal(node_ix1, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])

        # smaller slice
        node_ix2 = [r.ix for r in self.hic.regions(slice(4, 10, 1))]
        assert np.array_equal(node_ix2, [4, 5, 6, 7, 8, 9])

        # single ix
        node_ix3 = self.hic.regions(1)
        assert node_ix3.ix == 1

        # single chromosome
        node_ix4 = [r.ix for r in self.hic.regions('chr1')]
        assert np.array_equal(node_ix4, [0, 1, 2, 3, 4])

    def test_get_matrix(self):
        # whole matrix
        m = self.hic[:, :]
        # spot checks
        assert np.array_equal(m.shape, [12, 12])
        assert m[0, 0] == 1
        assert m[11, 11] == 78
        assert m[11, 0] == 12
        assert m[1, 10] == 22
        # symmetry check
        for i in range(0, 12):
            for j in range(i, 12):
                assert m[i, j] == m[j, i]

        # only upper half
        m = self.hic[:6, :]
        # spot checks
        assert np.array_equal(m.shape, [6, 12])
        assert m[0, 0] == 1
        assert m[5, 11] == 57
        assert m[5, 0] == 6
        assert m[1, 10] == 22

        # only lower half
        m = self.hic[6:, :]
        # spot checks
        assert np.array_equal(m.shape, [6, 12])
        assert m[0, 0] == 7
        assert m[5, 11] == 78
        assert m[5, 0] == 12
        assert m[1, 10] == 67

        # only left half
        m = self.hic[:, :6]
        # spot checks
        assert np.array_equal(m.shape, [12, 6])
        assert m[0, 0] == 1
        assert m[11, 5] == 57
        assert m[11, 0] == 12
        assert m[1, 4] == 16

        # only right half
        m = self.hic[:, 6:]
        # spot checks
        assert np.array_equal(m.shape, [12, 6])
        assert m[0, 0] == 7
        assert m[11, 5] == 78
        assert m[11, 0] == 63
        assert m[1, 4] == 22

        # top-left chunk
        m = self.hic[:6, :6]
        # spot checks
        assert np.array_equal(m.shape, [6, 6])
        assert m[0, 0] == 1
        assert m[5, 5] == 51
        assert m[5, 0] == 6
        assert m[1, 4] == 16

        # bottom_right chunk
        m = self.hic[6:, 6:]
        # spot checks
        assert np.array_equal(m.shape, [6, 6])
        assert m[0, 0] == 58
        assert m[5, 5] == 78
        assert m[5, 0] == 63
        assert m[1, 4] == 67

        # central chunk
        m = self.hic[3:9, 3:9]
        # spot checks
        assert np.array_equal(m.shape, [6, 6])
        assert m[0, 0] == 34
        assert m[5, 5] == 69
        assert m[5, 0] == 39
        assert m[1, 4] == 46

        # single row
        m = self.hic[1, 0:3]
        assert np.array_equal(m, [2., 13., 14.])

        # single row but only one value
        m = self.hic[1, 1:2]
        assert np.array_equal(m, [13.])

        # single col
        m = self.hic[0:3, 1]
        assert np.array_equal(m, [2., 13., 14.])

        # single col but only one value
        m = self.hic[1:2, 1]
        assert np.array_equal(m, [13.])

        # single value
        m = self.hic[1, 1]
        assert m == 13

        # empty array
        m = self.hic[1:1, 2:2]
        assert np.array_equal(m.shape, [0, 0])

    def test_merge(self):
        hic = self.hic_class()

        # add some nodes (120 to be exact)
        nodes = []
        for i in range(1, 5000, 1000):
            nodes.append(GenomicRegion(chromosome="chr1", start=i, end=i + 1000 - 1))
        for i in range(1, 3000, 1000):
            nodes.append(GenomicRegion(chromosome="chr2", start=i, end=i + 1000 - 1))
        for i in range(1, 2000, 500):
            nodes.append(GenomicRegion(chromosome="chr3", start=i, end=i + 1000 - 1))
        hic.add_regions(nodes, preserve_attributes=False)

        # add some edges with increasing weight for testing
        edges = []
        weight = 1
        for i in range(0, len(nodes)):
            for j in range(i, len(nodes)):
                edges.append(Edge(source=i, sink=j, weight=weight))
                weight += 1

        hic.add_edges(edges)

        # check length
        merged_hic_2x = Hic.merge([self.hic, hic])
        merged_hic_3x = Hic.merge([self.hic, hic, hic])
        hic.close()

        m = self.hic[:, :]
        m_merged_2x = merged_hic_2x[:, :]
        m_merged_3x = merged_hic_3x[:, :]

        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                assert m[i, j] == 0 or m[i, j] == m_merged_2x[i, j] / 2
                assert m[i, j] == 0 or m[i, j] == m_merged_3x[i, j] / 3
        merged_hic_2x.close()
        merged_hic_3x.close()

    def test_from_pairs(self):
        sam_file1 = os.path.join(self.dir, "test_matrix", "yeast.sample.chrI.1_sorted.sam")
        sam_file2 = os.path.join(self.dir, "test_matrix", "yeast.sample.chrI.2_sorted.sam")
        chrI = Chromosome.from_fasta(os.path.join(self.dir, "test_matrix", "chrI.fa"))
        genome = Genome(chromosomes=[chrI])

        pairs = ReadPairs()
        regions = genome.get_regions('HindIII')
        pairs.add_regions(regions)
        g = SamBamReadPairGenerator(sam_file1, sam_file2)
        pairs.add_read_pairs(g)
        genome.close()
        regions.close()

        pl = len(pairs)

        hic = pairs.to_hic(_hic_class=self.hic_class)

        assert len(hic.regions) == len(pairs.regions)
        pairs.close()

        reads = 0
        edge_set = set()
        for edge in hic.edges():
            key = (edge.source, edge.sink)
            assert key not in edge_set
            edge_set.add(key)
            reads += edge.weight

        assert reads == pl
        hic.close()

    def test_overlap_map(self):
        # ----|----|----|----|---|-----|-| new
        # -------|-------|-------|-------| old
        old_regions = [
            GenomicRegion(chromosome='chr1', start=1, end=8),
            GenomicRegion(chromosome='chr1', start=9, end=16),
            GenomicRegion(chromosome='chr1', start=17, end=24),
            GenomicRegion(chromosome='chr1', start=25, end=32)
        ]

        new_regions = [
            GenomicRegion(chromosome='chr1', start=1, end=5),
            GenomicRegion(chromosome='chr1', start=6, end=10),
            GenomicRegion(chromosome='chr1', start=11, end=15),
            GenomicRegion(chromosome='chr1', start=16, end=20),
            GenomicRegion(chromosome='chr1', start=21, end=24),
            GenomicRegion(chromosome='chr1', start=25, end=30),
            GenomicRegion(chromosome='chr1', start=31, end=32)
        ]

        overlap_map = _get_overlap_map(old_regions, new_regions)
        assert len(overlap_map[0]) == 2
        assert np.array_equal(overlap_map[0][0], [0, 1.0])
        assert np.array_equal(overlap_map[0][1], [1, 0.6])
        assert len(overlap_map[1]) == 3
        assert np.array_equal(overlap_map[1][0], [1, 0.4])
        assert np.array_equal(overlap_map[1][1], [2, 1.0])
        assert np.array_equal(overlap_map[1][2], [3, 0.2])
        assert len(overlap_map[2]) == 2
        assert np.array_equal(overlap_map[2][0], [3, 0.8])
        assert np.array_equal(overlap_map[2][1], [4, 1.0])
        assert len(overlap_map[3]) == 2
        assert np.array_equal(overlap_map[3][0], [5, 1.0])
        assert np.array_equal(overlap_map[3][1], [6, 1.0])

        # ----|----|-| new
        # --|--|--|--| old
        old_regions = list()
        old_regions.append(GenomicRegion(chromosome='chr1', start=1, end=3))
        old_regions.append(GenomicRegion(chromosome='chr1', start=4, end=6))
        old_regions.append(GenomicRegion(chromosome='chr1', start=7, end=9))
        old_regions.append(GenomicRegion(chromosome='chr1', start=10, end=12))

        new_regions = list()
        new_regions.append(GenomicRegion(chromosome='chr1', start=1, end=5))
        new_regions.append(GenomicRegion(chromosome='chr1', start=6, end=10))
        new_regions.append(GenomicRegion(chromosome='chr1', start=11, end=12))

        overlap_map = _get_overlap_map(old_regions, new_regions)
        assert len(overlap_map[0]) == 1
        assert np.array_equal(overlap_map[0][0], [0, 0.6])
        assert len(overlap_map[1]) == 2
        assert np.array_equal(overlap_map[1][0], [0, 0.4])
        assert np.array_equal(overlap_map[1][1], [1, 0.2])
        assert len(overlap_map[2]) == 1
        assert np.array_equal(overlap_map[2][0], [1, 0.6])
        assert len(overlap_map[3]) == 2
        assert np.array_equal(overlap_map[3][0], [1, 0.2])
        assert np.array_equal(overlap_map[3][1], [2, 1.0])

    def test_edge_splitting_rao(self):
        #     1         2         3         4
        # ---------|---------|---------|---------| old
        # ----|----|----|----|---------|---|--|--| new
        #  1    2    3    4       5      6  7  8
        overlap_map = {}
        overlap_map[1] = [[1, 1.0], [2, 1.0]]
        overlap_map[2] = [[3, 1.0], [4, 1.0]]
        overlap_map[3] = [[5, 1.0]]
        overlap_map[4] = [[6, 1.0], [7, 1.0], [8, 1.0]]

        original_edge = [1, 2, 12.0]
        new_edges = _edge_overlap_split_rao(original_edge, overlap_map)
        assert len(new_edges) == 4
        weight_sum = 0
        for new_edge in new_edges:
            weight_sum += new_edge[2]
        assert weight_sum == original_edge[2]

        original_edge = [1, 1, 12.0]
        new_edges = _edge_overlap_split_rao(original_edge, overlap_map)
        assert len(new_edges) == 3
        weight_sum = 0
        for new_edge in new_edges:
            weight_sum += new_edge[2]
        assert weight_sum == original_edge[2]

        original_edge = [1, 3, 9.0]
        new_edges = _edge_overlap_split_rao(original_edge, overlap_map)
        assert len(new_edges) == 2
        weight_sum = 0
        for new_edge in new_edges:
            weight_sum += new_edge[2]
        assert weight_sum == original_edge[2]

        original_edge = [3, 3, 9.0]
        new_edges = _edge_overlap_split_rao(original_edge, overlap_map)
        assert len(new_edges) == 1
        assert new_edges[0][2] == original_edge[2]

        original_edge = [1, 4, 9.0]
        new_edges = _edge_overlap_split_rao(original_edge, overlap_map)
        assert len(new_edges) == 4
        weight_sum = 0
        for new_edge in new_edges:
            weight_sum += new_edge[2]
        assert weight_sum == original_edge[2]

    def test_bin(self):
        original_reads = 0
        for edge in self.hic_cerevisiae.edges():
            original_reads += edge.weight

        def assert_binning(bin_size):
            binned = self.hic_cerevisiae.bin(bin_size)

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
            binned.close()

        bin_sizes = [500, 1000, 5000, 10000, 20000]
        for bin_size in bin_sizes:
            assert_binning(bin_size)

    def test_from_hic_sample(self):
        hic = self.hic_class()
        hic.add_region(GenomicRegion(chromosome='chr1', start=1, end=100))
        hic.add_region(GenomicRegion(chromosome='chr1', start=101, end=200))
        hic.flush()
        hic.add_edge([0, 0, 12])
        hic.add_edge([0, 1, 36])
        hic.add_edge([1, 1, 24])
        hic.flush()

        binned = self.hic_class()
        binned.add_region(GenomicRegion(chromosome='chr1', start=1, end=50))
        binned.add_region(GenomicRegion(chromosome='chr1', start=51, end=100))
        binned.add_region(GenomicRegion(chromosome='chr1', start=101, end=150))
        binned.add_region(GenomicRegion(chromosome='chr1', start=151, end=200))
        binned.flush()

        binned.load_from_hic(hic)

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
        hic.close()
        binned.close()

    def test_builtin_bin(self):
        hic = self.hic_class()
        hic.add_region(GenomicRegion(chromosome='chr1', start=1, end=100))
        hic.add_region(GenomicRegion(chromosome='chr1', start=101, end=200))
        hic.flush()
        hic.add_edge([0, 0, 12])
        hic.add_edge([0, 1, 36])
        hic.add_edge([1, 1, 24])
        hic.flush()

        binned = hic.bin(50)

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
        hic.close()
        binned.close()

    def test_knight_matrix_balancing(self):
        chrI = Chromosome.from_fasta(self.dir + "/test_matrix/chrI.fa")
        genome = Genome(chromosomes=[chrI])

        hic = self.hic_class()
        regions = genome.get_regions(10000)
        genome.close()
        hic.add_regions(regions)
        regions.close()
        hic.load_from_hic(self.hic_cerevisiae)

        m = hic[:, :]
        assert is_symmetric(m)

        kr_balancing(hic)
        m_corr = hic[:, :]
        assert is_symmetric(m_corr)

        for n in sum(m_corr):
            if np.ma.is_masked(n):
                continue
            assert abs(1.0 - n) < 1e-5 or n == 0
        hic.close()

    def test_knight_matrix_balancing_per_chromosome(self):
        chrI = Chromosome.from_fasta(self.dir + "/test_matrix/chrI.fa")
        genome = Genome(chromosomes=[chrI])

        hic = self.hic_class()
        regions = genome.get_regions(10000)
        genome.close()
        hic.add_regions(regions)
        regions.close()
        hic.load_from_hic(self.hic_cerevisiae)

        m = hic[:, :]
        assert is_symmetric(m)

        kr_balancing(hic, whole_matrix=False)
        m_corr = hic[:, :]
        assert is_symmetric(m_corr)

        for n in sum(m_corr):
            if np.ma.is_masked(n):
                continue
            assert abs(1.0 - n) < 1e-5 or n == 0
        hic.close()

    def test_ice_matrix_balancing(self):
        chrI = Chromosome.from_fasta(self.dir + "/test_matrix/chrI.fa")
        genome = Genome(chromosomes=[chrI])

        hic = self.hic_class()
        regions = genome.get_regions(10000)
        genome.close()
        hic.add_regions(regions)
        regions.close()
        hic.load_from_hic(self.hic_cerevisiae)

        m = hic[:, :]
        assert is_symmetric(m)

        ice_balancing(hic)
        m_corr = hic[:, :]
        assert is_symmetric(m_corr)

        sum_m_corr = sum(m_corr)
        for n in sum_m_corr:
            if np.ma.is_masked(n):
                continue
            assert (sum_m_corr[0] - 5 < n < sum_m_corr[0] + 5) or n == 0
        hic.close()

    def test_diagonal_filter(self):
        hic = self.hic

        m = hic[:]
        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                if i == j:
                    assert m[i, j] != 0

        hic.filter_diagonal(distance=1)
        m = hic[:]
        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                if abs(i - j) <= 1:
                    assert m[i, j] == 0

    def test_low_coverage_filter(self):
        hic = self.hic

        hic.filter_low_coverage_regions(cutoff=201)
        m = hic[:]
        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                if np.ma.is_masked(m[i, j]):
                    continue
                if i == 0 or j == 0 or i == 1 or j == 1:
                    assert m[i, j] == 0
                else:
                    assert m[i, j] != 0

    def test_to_cooler(self, tmpdir):
        cooler = pytest.importorskip("cooler")
        out = str(tmpdir.join("test_to_cooler.cool"))
        to_cooler(self.hic, out)
        c = cooler.Cooler(out)
        assert np.all(np.isclose(c.matrix(balance=False)[:], self.hic[:]))


class TestRegionMatrix:
    def setup_method(self, method):
        hic = Hic()

        m = np.zeros((12, 12))
        row_regions = []
        col_regions = []
        # add some nodes (120 to be exact)
        nodes = []
        for i in range(1, 5000, 1000):
            node = GenomicRegion(chromosome="chr1", start=i, end=i + 1000 - 1)
            nodes.append(node)
            row_regions.append(node)
            col_regions.append(node)
        for i in range(1, 3000, 1000):
            node = GenomicRegion(chromosome="chr2", start=i, end=i + 1000 - 1)
            nodes.append(node)
            row_regions.append(node)
            col_regions.append(node)
        for i in range(1, 2000, 500):
            node = GenomicRegion(chromosome="chr3", start=i, end=i + 1000 - 1)
            nodes.append(node)
            row_regions.append(node)
            col_regions.append(node)
        hic.add_regions(nodes)

        # add some edges with increasing weight for testing
        edges = []
        weight = 1
        for i in range(0, len(nodes)):
            for j in range(i, len(nodes)):
                edges.append(Edge(source=i, sink=j, weight=weight))
                m[i, j] = weight
                m[j, i] = weight
                weight += 1

        hic.add_edges(edges)

        self.hic = hic
        self.m = RegionMatrix(m, row_regions=row_regions, col_regions=col_regions)

    def teardown_method(self, method):
        self.hic.close()

    def test_create(self):
        hm = self.hic[:, :]
        m = self.m

        assert np.array_equal(hm, m)

        assert hasattr(m, 'row_regions')
        assert hasattr(m, 'col_regions')

    def test_repr(self):
        repr(self.m)

    def test_convert_key(self):
        key = self.m._convert_key('chr1:2001-5000', self.m._row_region_trees)
        assert key.start == 2
        assert key.stop == 5

        key = self.m._convert_key('chr1', self.m._row_region_trees)
        assert key.start == 0
        assert key.stop == 5

        key = self.m._convert_key('chr2', self.m._row_region_trees)
        assert key.start == 5
        assert key.stop == 8

    def test_select(self):
        def _equal(regions1, regions2):
            for i, region in enumerate(regions1):
                if region != regions2[i]:
                    return False
            return True

        res_all = self.m[:, :]
        assert _equal(res_all.shape, (12, 12))
        assert _equal(res_all.row_regions, self.m.row_regions[0:12])
        assert _equal(res_all.col_regions, self.m.col_regions[0:12])
        res_square = self.m[2:6, 2:6]
        assert _equal(res_square.shape, (4, 4))
        assert _equal(res_square.row_regions, self.m.row_regions[2:6])
        assert _equal(res_square.col_regions, self.m.col_regions[2:6])
        res_rect = self.m[2:6, 5:7]
        assert _equal(res_rect.shape, (4, 2))
        assert _equal(res_rect.row_regions, self.m.row_regions[2:6])
        assert _equal(res_rect.col_regions, self.m.col_regions[5:7])
        res_row = self.m[1]
        assert _equal(res_row.shape, (12,))
        assert _equal(res_row.row_regions, [self.m.row_regions[1]])
        assert _equal(res_row.col_regions, self.m.col_regions[:])
        res_col = self.m[:, 1]
        assert _equal(res_col.shape, (12,))
        assert _equal(res_col.row_regions, self.m.row_regions[:])
        assert _equal(res_col.col_regions, [self.m.col_regions[1]])
        res_row_sub = self.m[1, 2:6]
        assert _equal(res_row_sub.shape, (4,))
        assert _equal(res_row_sub.row_regions, [self.m.row_regions[1]])
        assert _equal(res_row_sub.col_regions, self.m.col_regions[2:6])
        res_col_sub = self.m[2:6, 1]
        assert _equal(res_col_sub.shape, (4,))
        assert _equal(res_col_sub.row_regions, self.m.row_regions[2:6])
        assert _equal(res_col_sub.col_regions, [self.m.col_regions[1]])
        res_single = self.m[0, 0]
        assert isinstance(res_single, float)

        hm = self.hic[:, :]
        for i, row_region in enumerate(hm.row_regions):
            assert row_region.start == self.m.row_regions[i].start
            assert row_region.end == self.m.row_regions[i].end
            assert row_region.chromosome == self.m.row_regions[i].chromosome

        for i, col_region in enumerate(hm.col_regions):
            assert col_region.start == self.m.col_regions[i].start
            assert col_region.end == self.m.col_regions[i].end
            assert col_region.chromosome == self.m.col_regions[i].chromosome

    def test_masked_matrix(self):
        hic = Hic()

        m = np.zeros((12, 12))
        row_regions = []
        col_regions = []
        # add some nodes (120 to be exact)
        nodes = []
        for i in range(1, 5000, 1000):
            node = GenomicRegion(chromosome="chr1", start=i, end=i + 1000 - 1)
            nodes.append(node)
            row_regions.append(node)
            col_regions.append(node)
        for i in range(1, 3000, 1000):
            node = GenomicRegion(chromosome="chr2", start=i, end=i + 1000 - 1)
            nodes.append(node)
            row_regions.append(node)
            col_regions.append(node)
        for i in range(1, 2000, 500):
            node = GenomicRegion(chromosome="chr3", start=i, end=i + 1000 - 1)
            nodes.append(node)
            row_regions.append(node)
            col_regions.append(node)
        hic.add_regions(nodes)

        # add some edges with increasing weight for testing
        edges = []
        weight = 1
        for i in range(0, len(nodes)):
            for j in range(i, len(nodes)):
                if i != 1 and j != 1 and i != 5 and j != 5:
                    edges.append(Edge(source=i, sink=j, weight=weight))
                    m[i, j] = weight
                    m[j, i] = weight
                weight += 1

        hic.add_edges(edges)

        m = hic.matrix()
        hic.close()

        # check masking
        for i in range(m.shape[0]):
            assert np.ma.is_masked(m[1, i])
            assert np.ma.is_masked(m[i, 1])
            assert np.ma.is_masked(m[5, i])
            assert np.ma.is_masked(m[i, 5])

        # check not masked
        not_masked = {0, 2, 3, 4, 6, 7, 8, 9, 10, 11}
        masked = {1, 5}

        for j in not_masked:
            for i in range(m.shape[0]):
                if i not in masked:
                    assert not np.ma.is_masked(m[i, j])
                    assert not np.ma.is_masked(m[j, i])
                else:
                    assert np.ma.is_masked(m[i, j])
                    assert np.ma.is_masked(m[j, i])
