import os
import gzip
import numpy as np
from genomic_regions import GenomicRegion
from kaic.matrix import Edge, RegionMatrixTable


test_dir = os.path.dirname(os.path.realpath(__file__))


def _get_test_regions(folder=test_dir, with_bias=False):
    chromosome_lengths = [('chr1', 228252215), ('chr2', 189746636)]

    region_ix = 0
    regions = []
    for chromosome, length in chromosome_lengths:
        if with_bias:
            bias_file = os.path.join(folder, 'test_matrix', 'rhesus_{}_biases.txt'.format(chromosome))
            bias = np.loadtxt(bias_file)
        else:
            bins = int(length / 2500000) + 1
            bias = np.repeat(1.0, bins)

        for i, start in enumerate(range(1, length, 2500000)):
            end = min(length, start + 2500000 - 1)
            region = GenomicRegion(chromosome=chromosome, start=start, end=end,
                                   ix=region_ix, bias=bias[i])
            regions.append(region)
            region_ix += 1

    return regions


def _get_test_edges(folder=test_dir):
    edges = []
    chromosomes = ['chr1', 'chr2']
    for i in range(len(chromosomes)):
        chromosome1 = chromosomes[i]
        for j in range(i, len(chromosomes)):
            chromosome2 = chromosomes[j]
            edges_file = os.path.join(folder, 'test_matrix',
                                      'rhesus_{}_{}_edges_uncorrected.txt.gz'.format(chromosome1,
                                                                                     chromosome2))

            with gzip.open(edges_file, 'r') as f:
                for line in f:
                    line = line.decode('utf-8')
                    line = line.rstrip()
                    if line == '':
                        continue

                    fields = line.split("\t")
                    source, sink, weight = int(fields[0]), int(fields[1]), float(fields[2])

                    edge = Edge(source, sink, weight=weight)
                    edges.append(edge)
    return edges


class RegionMatrixContainerTestFactory:

    def test_edges_iter(self):
        pass

    def test_get_edges(self):
        edges_dict = {(e.source, e.sink): e.weight for e in _get_test_edges()}

        for edge in self.matrix.edges(correct=False):
            assert np.isclose(edge.weight, edges_dict[(edge.source, edge.sink)])


class TestRegionMatrixTable(RegionMatrixContainerTestFactory):
    def setup_method(self, method):
        self.matrix = RegionMatrixTable()
        self.matrix.add_regions(_get_test_regions(with_bias=True))
        self.matrix.add_edges(_get_test_edges())

    def teardown_method(self, method):
        self.matrix.close()