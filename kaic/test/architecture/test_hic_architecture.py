from __future__ import division
from kaic.data.genomic import Hic, Node, Edge
from kaic.architecture.hic_architecture import PossibleContacts, ExpectedContacts, DirectionalityIndex, \
    InsulationIndex, ObservedExpectedRatio
import pytest
import numpy as np
from kaic.tools import dummy
import os.path


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
        for i in range(0, 5):
            for j in range(i, 5):
                edges.append(Edge(source=i, sink=j, weight=50))
        for i in range(6, 12):
            for j in range(i, 12):
                edges.append(Edge(source=i, sink=j, weight=75))
        for i in range(13, 18):
            for j in range(i, 18):
                edges.append(Edge(source=i, sink=j, weight=30))
        for i in range(18, 20):
            for j in range(i, 20):
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

            # find largest distance
            max_d = 0
            hic_matrix = self.hic.as_matrix(mask_missing=True)
            for i in range(hic_matrix.shape[0]):
                row_region = hic_matrix.row_regions[i]
                for j in range(hic_matrix.shape[1]):
                    col_region = hic_matrix.col_regions[j]

                    if row_region.chromosome == col_region.chromosome:
                        d = abs(row_region.ix - col_region.ix)
                        max_d = max(d, max_d)
            assert len(ec.intra_expected()) == max_d + 1

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


class TestObservedExpectedRatio:
    def setup_method(self, method):
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

    def teardown_method(self, method):
        self.hic.close()

    def test_observed_expected_ratio(self):
        obs = self.hic[:]
        with ExpectedContacts(self.hic) as ex:
            intra_expected = ex.intra_expected()
            inter_expected = ex.inter_expected()

            ex = np.empty(obs.shape)
            rd = self.hic.regions_dict
            for i in range(obs.shape[0]):
                for j in range(obs.shape[1]):
                    if rd[i].chromosome == rd[j].chromosome:
                        ex[i, j] = intra_expected[abs(i - j)]
                    else:
                        ex[i, j] = inter_expected

        with ObservedExpectedRatio(self.hic, per_chromosome=False) as oer:
            oer_m = oer[:]
            assert oer_m.shape == obs.shape

            for i in range(obs.shape[0]):
                for j in range(obs.shape[1]):
                    assert oer_m[i, j] - (obs[i, j] / ex[i, j]) < 0.001

    def test_observed_expected_ratio_region(self):
        obs = self.hic['chr1', 'chr1']
        with ExpectedContacts(self.hic, regions='chr1') as ex:
            intra_expected = ex.intra_expected()
            inter_expected = ex.inter_expected()

            ex = np.empty(obs.shape)
            rd = self.hic.regions_dict
            for i in range(obs.shape[0]):
                for j in range(obs.shape[1]):
                    if rd[i].chromosome == rd[j].chromosome:
                        ex[i, j] = intra_expected[abs(i - j)]
                    else:
                        ex[i, j] = inter_expected

        with ObservedExpectedRatio(self.hic, regions='chr1') as oer:
            oer_m = oer[:]
            assert oer_m.shape == obs.shape

            for i in range(obs.shape[0]):
                for j in range(obs.shape[1]):
                    assert oer_m[i, j] - (obs[i, j] / ex[i, j]) < 0.001


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
        for i in range(0, 5):
            for j in range(i, 5):
                edges.append(Edge(source=i, sink=j, weight=50))
        for i in range(6, 12):
            for j in range(i, 12):
                edges.append(Edge(source=i, sink=j, weight=75))
        for i in range(13, 18):
            for j in range(i, 18):
                edges.append(Edge(source=i, sink=j, weight=30))
        for i in range(18, 20):
            for j in range(i, 20):
                edges.append(Edge(source=i, sink=j, weight=50))

        hic.add_edges(edges)
        self.hic = hic

    def teardown_method(self, method):
        self.hic.close()

    def test_directionality_index(self):

        with DirectionalityIndex(self.hic, window_sizes=(3000, 5000)) as dip:
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

    def test_directionality_index_regions(self):
        with DirectionalityIndex(self.hic, window_sizes=(3000, 5000)) as dip:
            do = dip.directionality_index(window_size=3000)

        with DirectionalityIndex(self.hic, window_sizes=(3000, 5000), regions='chr1') as dip:
            d = dip.directionality_index(window_size=3000)

            assert d == do[:12]

            assert sum(d) > 0  # juuuuust making sure...

            # beginning of second TAD in chr1
            assert d[6] == max(d)
            # declining over TAD till end of chr1
            assert d[6] >= d[7] >= d[8] >= d[9] >= d[10] >= d[11]

            d5000 = dip.directionality_index(5000)
            assert len(d5000) == 12

            with pytest.raises(AttributeError):
                dip.directionality_index(10000)

        with DirectionalityIndex(self.hic, window_sizes=(3000, 5000), regions='chr2') as dip:
            d = dip.directionality_index(window_size=3000)

            assert d == do[12:]

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
        for i in range(0, 5):
            for j in range(i, 5):
                edges.append(Edge(source=i, sink=j, weight=50))
        for i in range(6, 12):
            for j in range(i, 12):
                edges.append(Edge(source=i, sink=j, weight=75))
        for i in range(13, 18):
            for j in range(i, 18):
                edges.append(Edge(source=i, sink=j, weight=30))
        for i in range(18, 20):
            for j in range(i, 20):
                edges.append(Edge(source=i, sink=j, weight=50))

        hic.add_edges(edges)
        self.hic = hic

    def teardown_method(self, method):
        self.hic.close()

    def test_insulation_index(self):
        with InsulationIndex(self.hic, window_sizes=(2000, 3000), impute_missing=True) as ins:
            d = ins.insulation_index(window_size=2000)

            assert np.isnan(d[0])
            assert d[2] == 50.0
            assert d[3] - 40.13888931274414 < 0.00001
            with pytest.raises(AttributeError):
                ins.directionality_index(10000)

    def test_insulation_index_regions(self):
        with InsulationIndex(self.hic, window_sizes=(2000, 3000)) as ins:
            do = ins.insulation_index(window_size=2000)

        with InsulationIndex(self.hic, window_sizes=(2000, 3000), regions='chr1') as ins:
            d = ins.insulation_index(window_size=2000)

            assert len(d) == len(do[:12])
            for i in range(len(d)):
                if np.isnan(d[i]):
                    assert np.isnan(do[i])
                else:
                    assert d[i] == do[i]

        with InsulationIndex(self.hic, window_sizes=(2000, 3000), regions='chr2') as ins:
            d = ins.insulation_index(window_size=2000)

            assert len(d) == len(do[12:])
            for i in range(len(d)):
                if np.isnan(d[i]):
                    assert np.isnan(do[i+12])
                else:
                    assert d[i] == do[i+12]

            with pytest.raises(AttributeError):
                ins.directionality_index(10000)

    def test_sparse_insulation_index(self):
        with dummy.sample_hic_matrix2() as hic:
            with InsulationIndex(hic, window_sizes=(1000, 2000, 3000, 4000, 5000)) as ii:
                ii_1000 = ii.insulation_index(1000)
                assert np.isnan(ii_1000[0])
                assert ii_1000[1] == 0
                assert ii_1000[2] == 7
                assert np.isnan(ii_1000[3])
                assert ii_1000[4] == 8
                assert np.isnan(ii_1000[5])
                assert np.isnan(ii_1000[6])
                assert ii_1000[7] == 0
                assert np.isnan(ii_1000[8])
                assert np.isnan(ii_1000[9])

                ii_2000 = ii.insulation_index(2000)
                assert np.isnan(ii_2000[0])
                assert np.isnan(ii_2000[1])
                assert ii_2000[2] == (7 + 2) / 2
                assert ii_2000[3] == (8 + 0) / 2
                assert ii_2000[4] == (8 + 5 + 0 + 9) / 4
                assert np.isnan(ii_2000[5])
                assert np.isnan(ii_2000[6])
                assert ii_2000[7] == (8 + 9 + 0 + 0) / 4
                assert np.isnan(ii_2000[8])
                assert np.isnan(ii_2000[9])

                ii_3000 = ii.insulation_index(3000)
                assert np.isnan(ii_3000[0])
                assert np.isnan(ii_3000[1])
                assert np.isnan(ii_3000[2])
                assert ii_3000[3] < (3+8+5)/6 + 0.0001
                assert ii_3000[4] < (8 + 8 + 5 + 9) / 6 + 0.0001
                assert np.isnan(ii_3000[5])
                assert np.isnan(ii_3000[6])
                assert np.isnan(ii_3000[7])
                assert np.isnan(ii_3000[8])
                assert np.isnan(ii_3000[9])


class TestBoundaryCalling:
    expected_values = {
        True: ("chr11:77633337-77643336", 6),
        False: ("chr11:77640001-77650000", 6)
    }

    def setup_method(self, method):
        self.dir = os.path.dirname(os.path.realpath(__file__))
        self.hic = Hic(os.path.join(self.dir, "../data/test_network/rao2014.chr11_77400000_78600000.hic"), mode="r")
        self.ins = InsulationIndex(self.hic, window_sizes=(50000,))
        self.ii = self.ins.insulation_index(50000)

    def teardown_method(self, method):
        self.hic.close()
        self.ins.close()

    @pytest.mark.parametrize("sub_bin_precision", [True, False])
    def test_boundaries(self, sub_bin_precision):
        boundaries = self.ins.boundaries(50000, delta_window=3, sub_bin_precision=sub_bin_precision)
        assert self.expected_values[sub_bin_precision][0] == str(boundaries[0])
        assert self.expected_values[sub_bin_precision][1] == len(boundaries)
