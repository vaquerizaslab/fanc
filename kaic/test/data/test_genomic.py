from __future__ import division
import numpy as np
from kaic.data.genomic import Chromosome, Genome, Hic, Node, Edge, \
    GenomicRegion, GenomicRegions, _get_overlap_map, _edge_overlap_split_rao, \
    RegionMatrix, RegionsTable, RegionMatrixTable, RegionPairs, AccessOptimisedRegionPairs, \
    AccessOptimisedRegionMatrixTable, AccessOptimisedHic, Bed, BigWig, Tabix
from kaic.tools.files import write_bed, write_gff, write_bigwig
from kaic.tools.general import which
from kaic.architecture.hic_architecture import BackgroundLigationFilter, ExpectedObservedEnrichmentFilter
import os.path
import pytest
from kaic.construct.seq import Reads, FragmentMappedReadPairs, AccessOptimisedReadPairs
from kaic.tools.matrix import is_symmetric
import kaic.correcting.knight_matrix_balancing as knight
import kaic.correcting.ice_matrix_balancing as ice
import tables as t


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
        assert np.array_equal(res, [12, 25])


class TestGenome:
    def setup_method(self, method):
        chr1 = Chromosome(name='chr1', length=10000, sequence='agcgctgctgaagcttcgatcgtaagcttc')
        chr2 = Chromosome(name='chr2', length=5000, sequence='gcgctgctgaagcttcgatcgtaagcttc')
        self.genome = Genome(chromosomes=[chr1, chr2])

    def teardown_method(self, method):
        self.genome.close()

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
            i += 1

    def test_node_list(self):
        regions = self.genome.get_regions('HindIII')

        assert len(regions) == 6
        for i in range(0, len(regions)):
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

        regions.close()

        nl = self.genome.get_regions(4000)

        assert len(nl) == 5
        for i in range(0, len(nl)):
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

        nl.close()

    def test_from_string(self):
        dir = os.path.dirname(os.path.realpath(__file__))
        genome = Genome.from_string(dir + '/test_genomic/chromosomes.fa')
        chr1 = genome[0]
        assert len(chr1) == 5
        chr2 = genome[1]
        assert len(chr2) == 3
        chr3 = genome[2]
        assert len(chr3) == 4
        chr4 = genome[3]
        assert len(chr4) == 2
        genome.close()


class TestGenomicRegion:
    def test_from_string(self):
        region1 = GenomicRegion.from_string('chr1')
        assert region1.chromosome == 'chr1'
        assert region1.start is None
        assert region1.end is None
        assert region1.strand is None

        region2 = GenomicRegion.from_string('chr1:0')
        assert region2.chromosome == 'chr1'
        assert region2.start == 0
        assert region2.end == 0
        assert region2.strand is None

        region3 = GenomicRegion.from_string('chr1:0-4956')
        assert region3.chromosome == 'chr1'
        assert region3.start == 0
        assert region3.end == 4956
        assert region3.strand is None

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


class RegionBasedTestFactory:
    def setup_method(self, method):
        self.regions = None
        self.empty_regions = None

    def test_get_item(self):
        region = self.regions.regions[0]
        assert isinstance(region, GenomicRegion)
        assert region.chromosome == 'chr1'
        assert region.start == 1
        assert region.end == 1000
        assert region.strand is None

    def test_len(self):
        assert len(self.regions.regions) == 29

    def test_iter(self):
        region_iter = self.regions.regions

        for i, region in enumerate(region_iter):
            start = 1 + i * 1000
            chromosome = 'chr1'
            if i > 22:
                start -= 23000
                chromosome = 'chr3'
            elif i > 8:
                start -= 9000
                chromosome = 'chr2'

            assert region.chromosome == chromosome
            assert region.start == start

    def test_region_bins(self):
        bins = self.regions.region_bins(GenomicRegion(chromosome='chr1', start=3400, end=8100))
        assert bins.start == 3
        assert bins.stop == 9

        bins = self.regions.region_bins('chr2:1-5000')
        assert bins.start == 9
        assert bins.stop == 14

        bins = self.regions.region_bins('chr2:1-5001')
        assert bins.start == 9
        assert bins.stop == 15

    def test_subset(self):
        # this is essentially the same as region_bins
        intersect = self.regions.subset(GenomicRegion(chromosome='chr1', start=3400, end=8100))
        assert len(list(intersect)) == 6


class TestGenomicRegions(RegionBasedTestFactory):
    def setup_method(self, method):
        chromosomes = [
            {'name': 'chr1', 'end': 10000},
            {'name': 'chr2', 'end': 15000},
            {'name': 'chr3', 'end': 7000}
        ]

        regions = []
        for chromosome in chromosomes:
            for start in range(1, chromosome["end"] - 1000, 1000):
                regions.append(GenomicRegion(start, start + 999, chromosome=chromosome["name"]))
        self.regions = GenomicRegions(regions)
        self.empty_regions = GenomicRegions()

    def test_add_region(self):
        # GenomicRegion
        self.empty_regions.add_region(GenomicRegion(start=1, end=1000, chromosome='chr1'))
        assert self.empty_regions[0].start == 1
        assert self.empty_regions[0].end == 1000
        assert self.empty_regions[0].chromosome == 'chr1'

        # dict
        self.empty_regions.add_region({'start': 1001, 'end': 2000, 'chromosome': 'chr1'})
        assert self.empty_regions[1].start == 1001
        assert self.empty_regions[1].end == 2000
        assert self.empty_regions[1].chromosome == 'chr1'

        # list
        self.empty_regions.add_region(['chr1', 2001, 3000])
        assert self.empty_regions[2].start == 2001
        assert self.empty_regions[2].end == 3000
        assert self.empty_regions[2].chromosome == 'chr1'


class TestRegionsTable(TestGenomicRegions):
    def setup_method(self, method):
        chromosomes = [
            {'name': 'chr1', 'end': 10000},
            {'name': 'chr2', 'end': 15000},
            {'name': 'chr3', 'end': 7000}
        ]

        regions = []
        for chromosome in chromosomes:
            for start in range(1, chromosome["end"] - 1000, 1000):
                regions.append(GenomicRegion(start, start + 999, chromosome=chromosome["name"]))
        self.regions = RegionsTable(regions)
        self.empty_regions = RegionsTable(additional_fields={'a': t.Int32Col(), 'b': t.StringCol(10)})

    def teardown_method(self, method):
        self.regions.close()
        self.empty_regions.close()

    def test_add_additional_fields(self):
        # GenomicRegion
        self.empty_regions.add_region(GenomicRegion(start=1, end=1000, chromosome='chr1', a=10, b='ten'))
        assert self.empty_regions[0].start == 1
        assert self.empty_regions[0].end == 1000
        assert self.empty_regions[0].chromosome == 'chr1'
        assert self.empty_regions[0].a == 10
        assert self.empty_regions[0].b == 'ten'

        # dict
        self.empty_regions.add_region({'start': 1001, 'end': 2000, 'chromosome': 'chr1', 'a': 11, 'b': 'eleven'})
        assert self.empty_regions[1].start == 1001
        assert self.empty_regions[1].end == 2000
        assert self.empty_regions[1].chromosome == 'chr1'
        assert self.empty_regions[1].a == 11
        assert self.empty_regions[1].b == 'eleven'

        # list
        self.empty_regions.add_region(['chr1', 2001, 3000])
        assert self.empty_regions[2].start == 2001
        assert self.empty_regions[2].end == 3000
        assert self.empty_regions[2].chromosome == 'chr1'
        assert self.empty_regions[2].a == 0
        assert self.empty_regions[2].b == ''


class TestBed(RegionBasedTestFactory):
    @pytest.fixture(autouse=True)
    def setup_method(self, tmpdir):
        chromosomes = [
            {'name': 'chr1', 'end': 10000},
            {'name': 'chr2', 'end': 15000},
            {'name': 'chr3', 'end': 7000}
        ]

        regions = []
        for chromosome in chromosomes:
            for start in range(1, chromosome["end"] - 1000, 1000):
                regions.append(GenomicRegion(start, start + 999, chromosome=chromosome["name"]))

        bed_file = os.path.join(str(tmpdir), 'test.bed')
        write_bed(bed_file, regions)

        self.regions = Bed(bed_file)


class TestBigWig(RegionBasedTestFactory):
    @pytest.fixture(autouse=True)
    def setup_method(self, tmpdir):
        chromosomes = [
            {'name': 'chr1', 'end': 10000},
            {'name': 'chr2', 'end': 15000},
            {'name': 'chr3', 'end': 7000}
        ]

        regions = []
        for chromosome in chromosomes:
            for start in range(1, chromosome["end"] - 1000, 1000):
                regions.append(GenomicRegion(start, start + 999, chromosome=chromosome["name"]))

        bw_file = os.path.join(str(tmpdir), 'test.bw')
        write_bigwig(bw_file, regions)

        self.regions = BigWig(bw_file)

    def test_get_item(self):
        pass


@pytest.mark.skipif(which('tabix') is None, reason="Cannot find tabix in PATH")
class TestTabix(RegionBasedTestFactory):
    @pytest.fixture(autouse=True)
    def setup_method(self, tmpdir):
        chromosomes = [
            {'name': 'chr1', 'end': 10000},
            {'name': 'chr2', 'end': 15000},
            {'name': 'chr3', 'end': 7000}
        ]

        regions = []
        for chromosome in chromosomes:
            for start in range(1, chromosome["end"] - 1000, 1000):
                regions.append(GenomicRegion(start, start + 999, chromosome=chromosome["name"]))

        bed_file = os.path.join(str(tmpdir), 'test.bed')
        write_bed(bed_file, regions)
        bed = Bed(bed_file)
        bed.tabix()

        tabix_file = os.path.join(str(tmpdir), 'test.bed.gz')

        self.regions = Tabix(tabix_file)

    def test_get_item(self):
        pass


class TestRegionPairs:
    def setup_method(self, method):
        self.rmt = RegionPairs(additional_fields={'weight': t.Int32Col(pos=0),
                                                  'foo': t.Int32Col(pos=1),
                                                  'bar': t.Float32Col(pos=2),
                                                  'baz': t.StringCol(50, pos=3)})

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
            node = Node(chromosome=chromosome, start=start, end=end)
            self.rmt.add_region(node, flush=False)
        self.rmt.flush()

        for i in range(10):
            for j in range(i, 10):
                edge = Edge(source=i, sink=j, weight=i * j, foo=i, bar=j, baz='x' + str(i * j))
                self.rmt.add_edge(edge, flush=False)
        self.rmt.flush()

        self.rp_class = RegionPairs

    def teardown_method(self, method):
        self.rmt.close()

    def test_create(self):
        rmt1 = RegionPairs()
        assert rmt1.field_names == ['source', 'sink']
        assert len(rmt1.edges()) == 0
        rmt1.close()

        rmt2 = RegionPairs(additional_fields={'foo': t.Int32Col(pos=0), 'bar': t.Float32Col(pos=1)})
        assert rmt2.field_names == ['source', 'sink', 'foo', 'bar']
        assert len(rmt2.edges()) == 0
        rmt2.close()

        class AdditionalFields(t.IsDescription):
            foo = t.Int32Col(pos=0)
            bar = t.Float32Col(pos=1)

        rmt3 = RegionPairs(additional_fields=AdditionalFields)
        assert rmt3.field_names == ['source', 'sink', 'foo', 'bar']
        assert len(rmt3.edges()) == 0
        rmt3.close()

        assert len(self.rmt.edges()) == 55
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
        for edge in self.rmt.edge_subset((slice(0, 2), slice(1, 3))):
            pair = (edge.source, edge.sink)
            if pair in covered:
                assert 0
            covered.add(pair)

        covered = set()
        for edge in self.rmt.edge_subset((slice(1, 3), slice(0, 2))):
            pair = (edge.source, edge.sink)
            if pair in covered:
                assert 0
            covered.add(pair)

        covered = set()
        for edge in self.rmt.edge_subset((slice(0, 3), slice(1, 2))):
            pair = (edge.source, edge.sink)
            if pair in covered:
                assert 0
            covered.add(pair)

    def test_edges_sorted(self):

        previous_weight = -1
        for edge in self.rmt.edges_sorted('weight'):
            assert edge.weight == edge.source * edge.sink
            assert edge.weight >= previous_weight
            previous_weight = edge.weight
            assert edge.foo == edge.source
            assert edge.bar == edge.sink
            assert edge.baz == 'x' + str(edge.source * edge.sink)

            with pytest.raises(AttributeError):
                assert edge.qux is None

    def test_lazy_edges(self):
        for edge in self.rmt.edges(lazy=True):
            i = edge.source
            j = edge.sink
            assert edge.weight == i * j
            assert edge.foo == i
            assert edge.bar == j
            assert edge.baz == 'x' + str(i * j)

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
        rmt = self.rp_class(additional_fields={'weight': t.Float64Col()})
        rmt.add_node(Node(chromosome='1', start=1, end=1000))
        rmt.add_edge(Edge(0, 0, weight=100))

        edge = rmt.edges[0]
        assert edge.source == 0
        assert edge.sink == 0
        assert edge.weight == 100
        rmt.close()

        rmt = self.rp_class(additional_fields={'weight': t.Float64Col()})
        rmt.add_node(Node(chromosome='1', start=1, end=1000))
        rmt.add_edge([0, 0, 100])

        edge = rmt.edges[0]
        assert edge.source == 0
        assert edge.sink == 0
        assert edge.weight == 100
        rmt.close()

        rmt = self.rp_class(additional_fields={'weight': t.Float64Col()})
        rmt.add_node(Node(chromosome='1', start=1, end=1000))
        rmt.add_edge({'source': 0, 'sink': 0, 'weight': 100})

        edge = rmt.edges[0]
        assert edge.source == 0
        assert edge.sink == 0
        assert edge.weight == 100
        rmt.close()


class TestAccessOptimisedRegionPairs(TestRegionPairs):
    def setup_method(self, method):
        self.rmt = AccessOptimisedRegionPairs(additional_fields={'weight': t.Int32Col(pos=0),
                                                                 'foo': t.Int32Col(pos=1),
                                                                 'bar': t.Float32Col(pos=2),
                                                                 'baz': t.StringCol(50, pos=3)})

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
            node = Node(chromosome=chromosome, start=start, end=end)
            self.rmt.add_region(node, flush=False)
        self.rmt.flush()

        for i in range(10):
            for j in range(i, 10):
                edge = Edge(source=i, sink=j, weight=i * j, foo=i, bar=j, baz='x' + str(i * j))
                self.rmt.add_edge(edge, flush=False)
        self.rmt.flush()

        self.rp_class = AccessOptimisedRegionPairs


class TestRegionMatrixTable:
    def setup_method(self, method):
        self.rmt = RegionMatrixTable(additional_fields={'weight': t.Int32Col(pos=0),
                                                        'foo': t.Int32Col(pos=1),
                                                        'bar': t.Float32Col(pos=2),
                                                        'baz': t.StringCol(50, pos=3)})

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
            node = Node(chromosome=chromosome, start=start, end=end)
            self.rmt.add_region(node, flush=False)
        self.rmt.flush()

        for i in range(10):
            for j in range(i, 10):
                edge = Edge(source=i, sink=j, weight=i * j, foo=i, bar=j, baz='x' + str(i * j))
                self.rmt.add_edge(edge, flush=False)
        self.rmt.flush()

    def teardown_method(self, method):
        self.rmt.close()

    def test_matrix(self):
        m = self.rmt.as_matrix()
        for row_region in m.row_regions:
            i = row_region.ix
            for col_region in m.col_regions:
                j = col_region.ix
                assert m[i, j] == m[j, i] == i * j

        m = self.rmt.as_matrix(values_from='foo')
        for row_region in m.row_regions:
            i = row_region.ix
            for col_region in m.col_regions:
                j = col_region.ix
                assert m[i, j] == m[j, i] == min(i, j)

        m = self.rmt.as_matrix(values_from='bar')
        for row_region in m.row_regions:
            i = row_region.ix
            for col_region in m.col_regions:
                j = col_region.ix
                assert m[i, j] == m[j, i] == max(i, j)

    def test_matrix_subset(self):
        m = self.rmt.as_matrix(key=('chr2', 'chr2'), values_from='bar')
        for i, row_region in enumerate(m.row_regions):
            for j, col_region in enumerate(m.col_regions):
                assert m[i, j] == m[j, i] == max(row_region.ix, col_region.ix)

        m = self.rmt.as_matrix(key=('chr2', 'chr3'), values_from='bar')
        for i, row_region in enumerate(m.row_regions):
            for j, col_region in enumerate(m.col_regions):
                assert m[i, j] == max(row_region.ix, col_region.ix)


class TestAccessOptimisedRegionMatrixTable(TestRegionMatrixTable):
    def setup_method(self, method):
        self.rmt = AccessOptimisedRegionMatrixTable(additional_fields={'weight': t.Int32Col(pos=0),
                                                                       'foo': t.Int32Col(pos=1),
                                                                       'bar': t.Float32Col(pos=2),
                                                                       'baz': t.StringCol(50, pos=3)})

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
            node = Node(chromosome=chromosome, start=start, end=end)
            self.rmt.add_region(node, flush=False)
        self.rmt.flush()

        for i in range(10):
            for j in range(i, 10):
                edge = Edge(source=i, sink=j, weight=i * j, foo=i, bar=j, baz='x' + str(i * j))
                self.rmt.add_edge(edge, flush=False)
        self.rmt.flush()


class TestHicBasic:
    def setup_method(self, method):
        self.dir = os.path.dirname(os.path.realpath(__file__))

        hic = Hic()

        # add some nodes (120 to be exact)
        nodes = []
        for i in range(1, 5000, 1000):
            nodes.append(Node(chromosome="chr1", start=i, end=i + 1000 - 1))
        for i in range(1, 3000, 1000):
            nodes.append(Node(chromosome="chr2", start=i, end=i + 1000 - 1))
        for i in range(1, 2000, 500):
            nodes.append(Node(chromosome="chr3", start=i, end=i + 1000 - 1))
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
        self.hic_cerevisiae = Hic(self.dir + "/test_genomic/cerevisiae.chrI.HindIII.hic")
        self.hic_class = Hic

    def teardown_method(self, method):
        self.hic_cerevisiae.close()
        self.hic.close()

    def test_initialize_empty(self):
        hic = self.hic_class()
        nodes = hic.nodes()
        edges = hic.edges()
        assert len(nodes) == 0
        assert len(edges) == 0
        hic.close()

    def test_save_and_load(self, tmpdir):
        dest_file = str(tmpdir) + "/hic.h5"

        hic1 = self.hic_class(file_name=dest_file, mode='w')
        hic1.add_node(GenomicRegion(1, 1000, 'chr1'))
        hic1.add_node(GenomicRegion(1, 1000, 'chr2'))
        hic1.add_edge([0, 1])
        hic1.close()

        hic2 = self.hic_class(dest_file, mode='r')
        nodes2 = hic2.nodes()
        edges2 = hic2.edges()
        assert len(nodes2) == 2
        assert len(edges2) == 1

        hic2.close()

    def test_nodes(self):
        nodes = self.hic.nodes()
        assert len(nodes) == 12

    def test_edges(self):
        edges = self.hic.edges()
        assert len(edges) == 78

    def test_intra_edges(self):
        edges = self.hic.edges(only_intrachromosomal=True)
        assert sum(1 for _ in edges) == 31

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
        node_ix1 = self.hic._getitem_nodes(slice(None, None, None), as_index=True)
        assert np.array_equal(node_ix1, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])

        # smaller slice
        node_ix2 = self.hic._getitem_nodes(slice(4, 10, 1), as_index=True)
        assert np.array_equal(node_ix2, [4, 5, 6, 7, 8, 9])

        # single ix
        node_ix3 = self.hic._getitem_nodes(1, as_index=True)
        assert node_ix3 == 1

        # single chromosome
        node_ix4 = self.hic._getitem_nodes('chr1', as_index=True)
        assert np.array_equal(node_ix4, [0, 1, 2, 3, 4])

        # HicNode
        node_ix5 = self.hic._getitem_nodes(Node(ix=1), as_index=True)
        assert node_ix5 == 1

        # list of items
        node_ix6 = self.hic._getitem_nodes(['chr1', 'chr3'], as_index=True)
        assert np.array_equal(node_ix6, [0, 1, 2, 3, 4, 8, 9, 10, 11])

        # nested list of items
        node_ix7 = self.hic._getitem_nodes(['chr1', ['chr2', 'chr3']], as_index=True)
        assert np.array_equal(node_ix7, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])

        # item repetition
        node_ix8 = self.hic._getitem_nodes(['chr3', 'chr3'], as_index=True)
        assert np.array_equal(node_ix8, [8, 9, 10, 11, 8, 9, 10, 11])

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

        # disjunct pieces
        m = self.hic[[1, 9], [4, 6]]
        # spot checks
        assert np.array_equal(m.shape, [2, 2])
        assert m[0, 0] == 16
        assert m[1, 1] == 61
        assert m[1, 0] == 48
        assert m[0, 1] == 18

        # single row
        m = self.hic[1, 0:3]
        assert np.array_equal(m, [2, 13, 14])

        # single row but only one value
        m = self.hic[1, 1:2]
        assert np.array_equal(m, [13])

        # single col
        m = self.hic[0:3, 1]
        assert np.array_equal(m, [2, 13, 14])

        # single col but only one value
        m = self.hic[1:2, 1]
        assert np.array_equal(m, [13])

        # single value
        m = self.hic[1, 1]
        assert m == 13

        # empty array
        m = self.hic[1:1, 2:2]
        assert np.array_equal(m.shape, [0, 0])

    def test_set_matrix(self):

        hic = self.hic_class(self.hic)

        n_edges = len(hic.edges)
        # for edge_table in hic._edge_table_iter():
        #    print edge_table

        # whole matrix
        old = hic[:, :]
        # set diagonal to zero
        for i in range(0, old.shape[0]):
            old[i, i] = 0
        hic[:, :] = old
        m = hic[:, :]

        assert np.array_equal(m.shape, old.shape)
        for i in range(0, m.shape[0]):
            for j in range(0, m.shape[1]):
                if i == j:
                    assert m[i, j] == 0
                else:
                    assert m[i, j] == old[i, j]

        assert len(hic.edges) < n_edges
        hic.close()

        # central matrix
        hic = self.hic_class(self.hic)
        old = hic[2:8, 2:10]
        # set border elements to zero
        # set checkerboard pattern
        for i in range(0, old.shape[0]):
            for j in range(0, old.shape[1]):
                if i == 0 or j == 0:
                    old[i, j] = 0
                elif i % 2 == 0 and j % 2 == 1:
                    old[i, j] = 0
                elif i % 2 == 1 and j % 2 == 0:
                    old[i, j] = 0
        hic[2:8, 2:10] = old
        m = hic[2:8, 2:10]
        hic.close()

        assert np.array_equal(m.shape, old.shape)
        for i in range(0, m.shape[0]):
            for j in range(0, m.shape[1]):
                if i == 0 or j == 0:
                    assert m[i, j] == 0
                elif i % 2 == 0 and j % 2 == 1:
                    assert m[i, j] == 0
                elif i % 2 == 1 and j % 2 == 0:
                    assert m[i, j] == 0
                else:
                    assert m[i, j] == old[i, j]

        hic = self.hic_class(self.hic)
        # row
        old = hic[1, 2:10]
        for i in range(0, 8, 2):
            old[i] = 0
        hic[1, 2:10] = old

        assert np.array_equal(hic[1, :], [2, 13, 0, 15, 0, 17, 0, 19, 0, 21, 22, 23])
        assert np.array_equal(hic[:, 1], [2, 13, 0, 15, 0, 17, 0, 19, 0, 21, 22, 23])
        hic.close()

        hic = self.hic_class(self.hic)
        # col
        old = hic[2:10, 1]
        for i in range(0, 8, 2):
            old[i] = 0
        hic[2:10, 1] = old

        assert np.array_equal(hic[1, :], [2, 13, 0, 15, 0, 17, 0, 19, 0, 21, 22, 23])
        assert np.array_equal(hic[:, 1], [2, 13, 0, 15, 0, 17, 0, 19, 0, 21, 22, 23])
        hic.close()

        # individual
        hic = self.hic_class(self.hic)
        hic[2, 1] = 0
        assert hic[2, 1] == 0
        assert hic[1, 2] == 0
        hic.close()

    def test_as_data_frame(self):
        df = self.hic.as_data_frame(('chr1', 'chr1'))
        assert np.array_equal(df.shape, [5, 5])
        assert np.array_equal(df.index, [1, 1001, 2001, 3001, 4001])
        assert np.array_equal(df.columns, [1, 1001, 2001, 3001, 4001])

    def test_merge(self):
        hic = self.hic_class()

        # add some nodes (120 to be exact)
        nodes = []
        for i in range(1, 5000, 1000):
            nodes.append(Node(chromosome="chr1", start=i, end=i + 1000 - 1))
        for i in range(1, 3000, 1000):
            nodes.append(Node(chromosome="chr2", start=i, end=i + 1000 - 1))
        for i in range(1, 2000, 400):
            nodes.append(Node(chromosome="chr4", start=i, end=i + 100 - 1))
        hic.add_nodes(nodes)

        # add some edges with increasing weight for testing
        edges = []
        weight = 1
        for i in range(0, len(nodes)):
            for j in range(i, len(nodes)):
                edges.append(Edge(source=i, sink=j, weight=weight))
                weight += 1

        hic.add_edges(edges)

        left = self.hic[:, :]
        right = hic[:, :]

        # check length
        original_length = len(self.hic.nodes())
        self.hic.merge(hic, _edge_buffer_size=5)

        assert len(self.hic.nodes()) == original_length + 5
        hic.close()

        merged = self.hic[:, :]
        double = [0, 1, 2, 3, 4, 5, 6, 7]
        for i in double:
            for j in double:
                assert merged[i, j] == left[i, j] + right[i, j]

        three = [8, 9, 10, 11]
        for i in double:
            for j in three:
                assert merged[i, j] == left[i, j]

        four = [12, 13, 14, 15, 16]
        for i in three:
            for j in four:
                assert merged[i, j] == 0

        for i in double:
            for j in four:
                assert merged[i, j] == right[i, j - 4]

    def test_multi_merge(self):
        def populate_hic(hic, s_s):
            # add some nodes (13 to be exact)
            nodes = []
            for i in range(1, 5000, 1000):
                nodes.append(Node(chromosome="chr1", start=i, end=i + 1000 - 1))
            for i in range(1, 3000, 1000):
                nodes.append(Node(chromosome="chr2", start=i, end=i + 1000 - 1))
            for i in range(1, 2000, 400):
                nodes.append(Node(chromosome="chr4", start=i, end=i + 100 - 1))
            hic.add_nodes(nodes)

            # add half as many random edges
            edges = []
            weight = 1

            for i, j in s_s:
                e = Edge(source=i, sink=j, weight=weight)
                edges.append(e)
                weight += 1
            hic.add_edges(edges)

        hic1 = self.hic_class()
        populate_hic(hic1, s_s=[(7, 4), (8, 3), (7, 3), (5, 5), (6, 4), (3, 0), (7, 6), (7, 0), (6, 0), (2, 1),
                                (6, 3), (2, 0), (6, 2), (5, 1), (8, 6), (8, 5), (6, 5), (3, 3)])
        assert hic1[:, :].sum() == 320  # 411
        hic2 = self.hic_class()
        populate_hic(hic2, s_s=[(7, 3), (3, 2), (0, 0), (3, 3), (8, 1), (6, 1), (6, 0), (8, 8), (7, 4), (2, 0),
                                (6, 2), (8, 7), (5, 1), (5, 4), (8, 3), (1, 0), (4, 1), (7, 6), (6, 5), (5, 5)])
        assert hic2[:, :].sum() == 385  # 443

        hic3 = self.hic_class()
        populate_hic(hic3, s_s=[(6, 4), (7, 5), (0, 0), (7, 0), (8, 1), (5, 2), (7, 6), (4, 4), (7, 4), (2, 0),
                                (6, 2), (7, 2), (4, 3), (4, 2), (8, 6), (5, 1), (4, 0), (8, 5), (6, 5), (8, 4)])
        assert hic3[:, :].sum() == 409  # 331

        hic_sum = hic1[:, :] + hic2[:, :] + hic3[:, :]
        hic1.merge([hic2, hic3], _edge_buffer_size=5)
        assert (hic1[:, :] == hic_sum).all()
        hic1.close()
        hic2.close()
        hic3.close()

    def test_merge_into(self):
        hic1 = self.hic_class()
        hic2 = self.hic_class()

        # add some nodes (120 to be exact)
        nodes = []
        for i in range(1, 5000, 1000):
            nodes.append(Node(chromosome="chr1", start=i, end=i + 1000 - 1))
        for i in range(1, 3000, 1000):
            nodes.append(Node(chromosome="chr2", start=i, end=i + 1000 - 1))
        for i in range(1, 2000, 500):
            nodes.append(Node(chromosome="chr3", start=i, end=i + 1000 - 1))
        hic1.add_nodes(nodes)
        hic2.add_nodes(nodes)

        # add some edges with increasing weight for testing
        edges = []
        weight = 1
        for i in range(0, len(nodes)):
            for j in range(i, len(nodes)):
                edges.append(Edge(source=i, sink=j, weight=weight))
                weight += 1

        hic1.add_edges(edges)
        hic2.add_edges(edges)

        merged_hic = self.hic_class.from_hic([hic1, hic2])

        left = hic1[:, :]
        right = hic2[:, :]

        # check length
        original_length = len(hic1.regions)
        assert len(merged_hic.regions) == original_length
        hic1.close()
        hic2.close()

        merged = merged_hic[:, :]
        merged_hic.close()

        for i in range(merged.shape[0]):
            for j in range(merged.shape[1]):
                assert merged[i, j] == left[i, j] + right[i, j]

        double = [0, 1, 2, 3, 4, 5, 6, 7]
        for i in double:
            for j in double:
                assert merged[i, j] == left[i, j] + right[i, j]

    def test_from_pairs(self):
        reads1 = Reads(self.dir + "/test_genomic/yeast.sample.chrI.1.sam")
        reads2 = Reads(self.dir + "/test_genomic/yeast.sample.chrI.2.sam")
        chrI = Chromosome.from_fasta(self.dir + "/test_genomic/chrI.fa")
        genome = Genome(chromosomes=[chrI])
        pairs = FragmentMappedReadPairs()
        regions = genome.get_regions('HindIII')
        pairs.load(reads1, reads2, regions)
        reads1.close()
        reads2.close()
        genome.close()
        regions.close()

        pl = len(pairs)

        hic = self.hic_class()
        hic.load_read_fragment_pairs(pairs)

        assert len(hic._regions) == len(pairs._regions)
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

    def test_from_ao_pairs(self):
        reads1 = Reads(self.dir + "/test_genomic/yeast.sample.chrI.1.sam")
        reads2 = Reads(self.dir + "/test_genomic/yeast.sample.chrI.2.sam")
        chrI = Chromosome.from_fasta(self.dir + "/test_genomic/chrI.fa")
        genome = Genome(chromosomes=[chrI])
        pairs_old = FragmentMappedReadPairs()
        regions = genome.get_regions('HindIII')
        pairs_old.load(reads1, reads2, regions)
        reads1.close()
        reads2.close()
        genome.close()
        regions.close()

        reads1 = Reads(self.dir + "/test_genomic/yeast.sample.chrI.1.sam")
        reads2 = Reads(self.dir + "/test_genomic/yeast.sample.chrI.2.sam")
        chrI = Chromosome.from_fasta(self.dir + "/test_genomic/chrI.fa")
        genome = Genome(chromosomes=[chrI])
        pairs = AccessOptimisedReadPairs()
        regions = genome.get_regions('HindIII')
        pairs.load(reads1, reads2, regions)
        reads1.close()
        reads2.close()
        genome.close()
        regions.close()

        pl = len(pairs)
        pl_old = len(pairs_old)

        assert pl == pl_old

        hic = self.hic_class()
        hic.load_read_fragment_pairs(pairs)

        assert len(hic._regions) == len(pairs._regions)
        pairs.close()
        pairs_old.close()

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
            Node(chromosome='chr1', start=1, end=8),
            Node(chromosome='chr1', start=9, end=16),
            Node(chromosome='chr1', start=17, end=24),
            Node(chromosome='chr1', start=25, end=32)
        ]

        new_regions = [
            Node(chromosome='chr1', start=1, end=5),
            Node(chromosome='chr1', start=6, end=10),
            Node(chromosome='chr1', start=11, end=15),
            Node(chromosome='chr1', start=16, end=20),
            Node(chromosome='chr1', start=21, end=24),
            Node(chromosome='chr1', start=25, end=30),
            Node(chromosome='chr1', start=31, end=32)
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
        old_regions = []
        old_regions.append(Node(chromosome='chr1', start=1, end=3))
        old_regions.append(Node(chromosome='chr1', start=4, end=6))
        old_regions.append(Node(chromosome='chr1', start=7, end=9))
        old_regions.append(Node(chromosome='chr1', start=10, end=12))

        new_regions = []
        new_regions.append(Node(chromosome='chr1', start=1, end=5))
        new_regions.append(Node(chromosome='chr1', start=6, end=10))
        new_regions.append(Node(chromosome='chr1', start=11, end=12))

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

    def test_from_hic(self):
        chrI = Chromosome.from_fasta(self.dir + "/test_genomic/chrI.fa")
        genome = Genome(chromosomes=[chrI])

        original_reads = 0
        for edge in self.hic_cerevisiae.edges():
            original_reads += edge.weight

        def assert_binning(hic, bin_size):
            binned = self.hic_class()
            assert len(binned.nodes()) == 0
            regions = genome.get_regions(bin_size)
            binned.add_regions(regions)
            binned.load_from_hic(hic)

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
            regions.close()
            binned.close()

        bin_sizes = [500, 1000, 5000, 10000, 20000]
        for bin_size in bin_sizes:
            assert_binning(self.hic_cerevisiae, bin_size)
        genome.close()

    def test_from_hic_sample(self):
        hic = self.hic_class()
        hic.add_region(GenomicRegion(chromosome='chr1', start=1, end=100))
        hic.add_region(GenomicRegion(chromosome='chr1', start=101, end=200))
        hic.add_edge([0, 0, 12])
        hic.add_edge([0, 1, 36])
        hic.add_edge([1, 1, 24])

        binned = self.hic_class()
        binned.add_region(GenomicRegion(chromosome='chr1', start=1, end=50))
        binned.add_region(GenomicRegion(chromosome='chr1', start=51, end=100))
        binned.add_region(GenomicRegion(chromosome='chr1', start=101, end=150))
        binned.add_region(GenomicRegion(chromosome='chr1', start=151, end=200))

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
        hic.add_edge([0, 0, 12])
        hic.add_edge([0, 1, 36])
        hic.add_edge([1, 1, 24])

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
        chrI = Chromosome.from_fasta(self.dir + "/test_genomic/chrI.fa")
        genome = Genome(chromosomes=[chrI])

        hic = self.hic_class()
        regions = genome.get_regions(10000)
        genome.close()
        hic.add_regions(regions)
        regions.close()
        hic.load_from_hic(self.hic_cerevisiae)

        m = hic[:, :]
        assert is_symmetric(m)

        knight.correct(hic)
        m_corr = hic[:, :]
        assert is_symmetric(m_corr)

        for n in sum(m_corr):
            assert abs(1.0 - n) < 1e-5 or n == 0
        hic.close()

    def test_knight_matrix_balancing_copy(self):
        chrI = Chromosome.from_fasta(self.dir + "/test_genomic/chrI.fa")
        genome = Genome(chromosomes=[chrI])

        hic = self.hic_class()
        regions = genome.get_regions(10000)
        genome.close()
        hic.add_regions(regions)
        regions.close()
        hic.load_from_hic(self.hic_cerevisiae)

        m = hic[:, :]
        assert is_symmetric(m)

        hic_new = knight.correct(hic, copy=True)
        m_corr = hic_new[:, :]
        assert is_symmetric(m_corr)
        assert m_corr.shape == m.shape

        assert hic is not hic_new

        for n in sum(m_corr):
            assert abs(1.0 - n) < 1e-5 or n == 0

        hic_new2 = knight.correct(hic, copy=True, whole_matrix=False)
        m_corr_pc = hic_new2[:, :]
        assert is_symmetric(m_corr_pc)
        assert m_corr_pc.shape == m.shape

        assert hic is not hic_new2

        for i in range(m_corr.shape[0]):
            for j in range(m_corr.shape[1]):
                assert abs(m_corr[i, j] - m_corr_pc[i, j]) < 0.0001

        hic.close()
        hic_new.close()
        hic_new2.close()

    def test_ice_matrix_balancing(self):
        chrI = Chromosome.from_fasta(self.dir + "/test_genomic/chrI.fa")
        genome = Genome(chromosomes=[chrI])

        hic = self.hic_class()
        regions = genome.get_regions(10000)
        genome.close()
        hic.add_regions(regions)
        regions.close()
        hic.load_from_hic(self.hic_cerevisiae)

        m = hic[:, :]
        assert is_symmetric(m)

        ice.correct(hic)
        m_corr = hic[:, :]
        assert is_symmetric(m_corr)

        sum_m_corr = sum(m_corr)
        for n in sum_m_corr:
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
                if i == 0 or j == 0 or i == 1 or j == 1:
                    assert m[i, j] == 0
                else:
                    assert m[i, j] != 0

    def test_filter_background_ligation(self):
        blf = BackgroundLigationFilter(self.hic, fold_change=2)
        assert blf.cutoff - 2 * (610 + 405 + 734) / 47 < 0.001
        for edge in self.hic.edges(lazy=True):
            if edge.weight < blf.cutoff:
                assert not blf.valid_edge(edge)
            else:
                assert blf.valid_edge(edge)

    def test_filter_expected_observed_enrichment(self):
        eof = ExpectedObservedEnrichmentFilter(self.hic, fold_change=1)
        previous = len(self.hic.edges)
        self.hic.filter(eof)
        assert len(self.hic.edges) == previous - 14 - 23  # 15 intra, 23 inter filtered

    def test_to_cooler(self, tmpdir):
        cooler = pytest.importorskip("cooler")
        out = str(tmpdir.join("test_to_cooler.cool"))
        self.hic.to_cooler(out)
        c = cooler.Cooler(out)
        assert np.all(np.isclose(c.matrix(balance=False)[:], self.hic[:]))


class TestAccessOptimisedHic(TestHicBasic):
    def setup_method(self, method):
        self.dir = os.path.dirname(os.path.realpath(__file__))
        self.hic_class = AccessOptimisedHic

        hic = self.hic_class()

        # add some nodes (120 to be exact)
        nodes = []
        for i in range(1, 5000, 1000):
            nodes.append(Node(chromosome="chr1", start=i, end=i + 1000 - 1))
        for i in range(1, 3000, 1000):
            nodes.append(Node(chromosome="chr2", start=i, end=i + 1000 - 1))
        for i in range(1, 2000, 500):
            nodes.append(Node(chromosome="chr3", start=i, end=i + 1000 - 1))
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
        self.hic_cerevisiae = AccessOptimisedHic(self.dir + "/test_genomic/cerevisiae.chrI.HindIII.fa.hic")

    def teardown_method(self, method):
        self.hic.close()
        self.hic_cerevisiae.close()


class TestRegionMatrix:
    def setup_method(self, method):
        hic = Hic()

        m = np.zeros((12, 12))
        row_regions = []
        col_regions = []
        # add some nodes (120 to be exact)
        nodes = []
        for i in range(1, 5000, 1000):
            node = Node(chromosome="chr1", start=i, end=i + 1000 - 1)
            nodes.append(node)
            row_regions.append(node)
            col_regions.append(node)
        for i in range(1, 3000, 1000):
            node = Node(chromosome="chr2", start=i, end=i + 1000 - 1)
            nodes.append(node)
            row_regions.append(node)
            col_regions.append(node)
        for i in range(1, 2000, 500):
            node = Node(chromosome="chr3", start=i, end=i + 1000 - 1)
            nodes.append(node)
            row_regions.append(node)
            col_regions.append(node)
        hic.add_nodes(nodes)

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
        key = self.m._convert_key('chr1:2001-5000', self.m.row_regions)
        assert key.start == 2
        assert key.stop == 5

        key = self.m._convert_key('chr1', self.m.row_regions)
        assert key.start == 0
        assert key.stop == 5

        key = self.m._convert_key('chr2', self.m.row_regions)
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
        assert res_row.row_regions == self.m.row_regions[1]
        assert _equal(res_row.col_regions, self.m.col_regions[:])
        res_col = self.m[:, 1]
        assert _equal(res_col.shape, (12,))
        assert _equal(res_col.row_regions, self.m.row_regions[:])
        assert res_col.col_regions == self.m.col_regions[1]
        res_row_sub = self.m[1, 2:6]
        assert _equal(res_row_sub.shape, (4,))
        assert res_row_sub.row_regions == self.m.row_regions[1]
        assert _equal(res_row_sub.col_regions, self.m.col_regions[2:6])
        res_col_sub = self.m[2:6, 1]
        assert _equal(res_col_sub.shape, (4,))
        assert _equal(res_col_sub.row_regions, self.m.row_regions[2:6])
        assert res_col_sub.col_regions == self.m.col_regions[1]
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
        self.hic[1, :] = np.zeros(12)
        self.hic[5, :] = np.zeros(12)

        m = self.hic.as_matrix(mask_missing=True)

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

    def test_imputed_matrix(self):
        self.hic[1, :] = np.zeros(12)
        self.hic[5, :] = np.zeros(12)

        m = self.hic.as_matrix(impute_missing=True)

        assert is_symmetric(m)

        assert m[1, 1] == 52.0
