from kaic.data.genomic import Hic, Node, Edge


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
        with self.hic.architecture.possible_contacts() as pc:
            pass
