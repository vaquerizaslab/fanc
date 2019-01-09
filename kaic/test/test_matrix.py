import os
import numpy as np
from genomic_regions import GenomicRegion
from kaic.matrix import Edge
from kaic.hic import Hic


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

    def test_get_edges_uncorrected(self):
        edges_dict = {(e.source, e.sink): e.weight for e in _get_test_edges(norm=False)}

        for edge in self.matrix.edges(norm=False):
            assert np.isclose(edge.weight,
                              edges_dict[(edge.source, edge.sink)],
                              rtol=1e-03)

    def test_get_edges(self):
        edges_dict = {(e.source, e.sink): e.weight for e in _get_test_edges(norm=True)}

        for edge in self.matrix.edges(norm=True):
            if np.isnan(edges_dict[(edge.source, edge.sink)]):
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
