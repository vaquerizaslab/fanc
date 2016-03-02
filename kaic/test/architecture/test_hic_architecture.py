from __future__ import division
from kaic.data.genomic import Hic, Node, Edge
from kaic.architecture.hic_architecture import PossibleContacts, ExpectedContacts, DirectionalityIndex, InsulationIndex
import pytest
import numpy as np


class TestHicArchitecture:
    def setup_method(self, method):
        # make TAD-like structures for testing
        hic = Hic()

        nodes = []
        for i in range(1, 12000, 1000):
            node = Node(chromosome="chr1", start=i, end=i+1000-1)
            nodes.append(node)
        for i in range(1, 4000, 500):
            node = Node(chromosome="chr2", start=i, end=i+500-1)
            nodes.append(node)
        hic.add_nodes(nodes)

        edges = []
        for i in xrange(0, 5):
            for j in xrange(i, 5):
                edges.append(Edge(source=i, sink=j, weight=50))
        for i in xrange(6, 12):
            for j in xrange(i, 12):
                edges.append(Edge(source=i, sink=j, weight=75))
        for i in xrange(13, 18):
            for j in xrange(i, 18):
                edges.append(Edge(source=i, sink=j, weight=30))
        for i in xrange(18, 20):
            for j in xrange(i, 20):
                edges.append(Edge(source=i, sink=j, weight=50))

        hic.add_edges(edges)
        self.hic = hic

    def teardown_method(self, method):
        self.hic.close()

    def test_call_architecture(self):
        with self.hic.architecture.possible_contacts:
            pass

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


class TestDirectionalityIndex:
    def setup_method(self, method):
        # make TAD-like structures for testing
        hic = Hic()

        nodes = []
        for i in range(1, 12000, 1000):
            node = Node(chromosome="chr1", start=i, end=i+1000-1)
            nodes.append(node)
        for i in range(1, 4000, 500):
            node = Node(chromosome="chr2", start=i, end=i+500-1)
            nodes.append(node)
        hic.add_nodes(nodes)

        edges = []
        for i in xrange(0, 5):
            for j in xrange(i, 5):
                edges.append(Edge(source=i, sink=j, weight=50))
        for i in xrange(6, 12):
            for j in xrange(i, 12):
                edges.append(Edge(source=i, sink=j, weight=75))
        for i in xrange(13, 18):
            for j in xrange(i, 18):
                edges.append(Edge(source=i, sink=j, weight=30))
        for i in xrange(18, 20):
            for j in xrange(i, 20):
                edges.append(Edge(source=i, sink=j, weight=50))

        hic.add_edges(edges)
        self.hic = hic

    def teardown_method(self, method):
        self.hic.close()

    def test_directionality_index(self):

        with DirectionalityIndex(self.hic, window_sizes=(3000, 5000)) as dip:
            print len(dip.regions)
            d = dip.directionality_index(window_size=3000)

            assert sum(d) > 0  # juuuuust making sure...

            # beginning of second TAD in chr1
            assert d[6] == max(d)
            # declining over TAD till end of chr1
            assert d[6] >= d[7] >= d[8] >= d[9] >= d[10] >= d[11]

            # beginning of second TAD in chr2
            assert d[18] == max(d[13:20])

            d5000 = dip.directionality_index(5000)
            assert len(d5000) == len(self.hic.regions)

            with pytest.raises(AttributeError):
                dip.directionality_index(10000)

    def test_boundary_distances(self):

        with DirectionalityIndex(self.hic, window_sizes=(3000, 5000)) as dip:
            boundary_dist = dip._get_boundary_distances()
            assert len(boundary_dist) == len(self.hic.regions)
            assert np.array_equal(boundary_dist, [0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0, 0, 1, 2, 3, 3, 2, 1, 0])


class TestInsulationIndex:
    def setup_method(self, method):
        # make TAD-like structures for testing
        hic = Hic()

        nodes = []
        for i in range(1, 12000, 1000):
            node = Node(chromosome="chr1", start=i, end=i+1000-1)
            nodes.append(node)
        for i in range(1, 4000, 500):
            node = Node(chromosome="chr2", start=i, end=i+500-1)
            nodes.append(node)
        hic.add_nodes(nodes)

        edges = []
        for i in xrange(0, 5):
            for j in xrange(i, 5):
                edges.append(Edge(source=i, sink=j, weight=50))
        for i in xrange(6, 12):
            for j in xrange(i, 12):
                edges.append(Edge(source=i, sink=j, weight=75))
        for i in xrange(13, 18):
            for j in xrange(i, 18):
                edges.append(Edge(source=i, sink=j, weight=30))
        for i in xrange(18, 20):
            for j in xrange(i, 20):
                edges.append(Edge(source=i, sink=j, weight=50))

        hic.add_edges(edges)
        self.hic = hic

    def teardown_method(self, method):
        self.hic.close()

    def test_insulation_index(self):

        with InsulationIndex(self.hic, window_sizes=(2000, 3000)) as ins:
            print len(ins.regions)
            d = ins.insulation_index(window_size=2000)

            assert np.isnan(d[0])
            assert np.isnan(d[1])
            assert d[2] == 50.0
            assert d[3] - 38.77083206176758 < 0.00001

            with pytest.raises(AttributeError):
                ins.directionality_index(10000)
