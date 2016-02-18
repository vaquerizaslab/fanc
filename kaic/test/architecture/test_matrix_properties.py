from __future__ import division
from kaic.data.genomic import Hic, Node, Edge, PossibleContacts, ExpectedContacts


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
