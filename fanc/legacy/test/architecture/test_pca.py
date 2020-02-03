from __future__ import division
from fanc.architecture.pca import do_pca
from fanc.data.genomic import Hic, Node, Edge
import math
import random
random.seed(1)


class TestPCA:
    def setup_method(self, method):

        def get_hic(max_random_offset=2.0, amplify=False):
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
            weight = 1.0
            for i in range(0, len(nodes)):
                for j in range(i, len(nodes)):
                    if amplify and i > 7:
                        a = 10.0
                    else:
                        a = 1.0
                    edges.append(Edge(source=i, sink=j, weight=a+weight*max_random_offset*random.random()))
                    weight += 1.0

            hic.add_edges(edges)
            return hic

        self.hic1 = get_hic()
        self.hic2 = get_hic()
        self.hic3 = get_hic(amplify=True)
        self.hic4 = get_hic(amplify=True)

    def teardown_method(self, method):
        self.hic1.close()
        self.hic2.close()
        self.hic3.close()
        self.hic4.close()

    def test_pca(self):
        hics = [self.hic1, self.hic2, self.hic3, self.hic4]

        pca, res = do_pca(hics, sample_size=25)

        # euclidian distances
        def euclid(xy1, xy2):
            return math.sqrt((xy1[0] - xy2[0]) ** 2 + (xy1[1] - xy2[1]) ** 2)

        related_pair_dist = euclid((res[2, 0], res[2, 1]), (res[3, 0], res[3, 1]))
        pairs = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3))
        for pair in pairs:
            i = pair[0]
            j = pair[1]
            assert related_pair_dist < euclid((res[i, 0], res[i, 1]), (res[j, 0], res[j, 1]))

    def test_pca_chromosome(self):
        hics = [self.hic1, self.hic2, self.hic3, self.hic4]

        pca, res = do_pca(hics, sample_size=25, regions='chr3')

        # euclidian distances
        def euclid(xy1, xy2):
            return math.sqrt((xy1[0] - xy2[0]) ** 2 + (xy1[1] - xy2[1]) ** 2)

        related_pair_dist = euclid((res[2, 0], res[2, 1]), (res[3, 0], res[3, 1]))
        pairs = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3))
        for pair in pairs:
            i = pair[0]
            j = pair[1]
            assert related_pair_dist < euclid((res[i, 0], res[i, 1]), (res[j, 0], res[j, 1]))

