from __future__ import division
import numpy as np
from kaic.regions import Chromosome, Genome, RegionsTable
from genomic_regions import GenomicRegion
import pytest

import os.path

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
                assert region.chromosome == 'chr1'
                assert region.start == 1
                assert region.end == 12
            if i == 5:
                assert region.chromosome == 'chr2'
                assert region.start == 25
                assert region.end == 5000
            i += 1

        regions.close()

        nl = self.genome.get_regions(4000)

        assert len(nl) == 5
        for i in range(0, len(nl)):
            node = nl[i]
            if i == 0:
                assert node.chromosome == 'chr1'
                assert node.start == 1
                assert node.end == 4000
            if i == 5:
                assert node.chromosome == 'chr2'
                assert node.start == 4001
                assert node.end == 5000
            i += 1

        nl.close()

    def test_from_string(self):
        dir = os.path.dirname(os.path.realpath(__file__))
        genome = Genome.from_string(dir + '/test_regions/chromosomes.fa')
        chr1 = genome[0]
        assert len(chr1) == 5
        chr2 = genome[1]
        assert len(chr2) == 3
        chr3 = genome[2]
        assert len(chr3) == 4
        chr4 = genome[3]
        assert len(chr4) == 2
        genome.close()


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

    @pytest.mark.parametrize("lazy", [True, False])
    def test_iter(self, lazy):
        region_iter = self.regions.regions(lazy=lazy)

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

    @pytest.mark.parametrize("lazy", [True, False])
    def test_region_subset(self, lazy):
        region_iter = self.regions.regions('chr1', lazy=lazy)

        for i, region in enumerate(region_iter):
            start = 1 + i * 1000

            assert region.chromosome == 'chr1'
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


class TestRegionsTable(RegionBasedTestFactory):
    def setup_method(self, method):
        chromosomes = [
            {'name': 'chr1', 'end': 10000},
            {'name': 'chr2', 'end': 15000},
            {'name': 'chr3', 'end': 7000}
        ]

        regions = []
        for chromosome in chromosomes:
            for start in range(1, chromosome["end"] - 1000, 1000):
                regions.append(GenomicRegion(start=start, end=start + 999, chromosome=chromosome["name"]))
        self.regions = RegionsTable()
        self.regions.add_regions(regions)
        self.empty_regions = RegionsTable(additional_fields={'a': t.Int32Col(), 'b': t.StringCol(10)})

    def teardown_method(self, method):
        self.regions.close()
        self.empty_regions.close()

    def test_add_additional_fields(self):
        # GenomicRegion
        self.empty_regions.add_region(GenomicRegion(start=1, end=1000, chromosome='chr1', a=10, b='ten'))
        self.empty_regions.flush()
        assert self.empty_regions[0].start == 1
        assert self.empty_regions[0].end == 1000
        assert self.empty_regions[0].chromosome == 'chr1'
        assert self.empty_regions[0].a == 10
        assert self.empty_regions[0].b == 'ten'

        # dict
        self.empty_regions.add_region({'start': 1001, 'end': 2000, 'chromosome': 'chr1', 'a': 11, 'b': 'eleven'})
        self.empty_regions.flush()
        assert self.empty_regions[1].start == 1001
        assert self.empty_regions[1].end == 2000
        assert self.empty_regions[1].chromosome == 'chr1'
        assert self.empty_regions[1].a == 11
        assert self.empty_regions[1].b == 'eleven'

        # list
        self.empty_regions.add_region(['chr1', 2001, 3000])
        self.empty_regions.flush()
        assert self.empty_regions[2].start == 2001
        assert self.empty_regions[2].end == 3000
        assert self.empty_regions[2].chromosome == 'chr1'
        assert self.empty_regions[2].a == 0
        assert self.empty_regions[2].b == ''

    def test_add_region(self):
        # GenomicRegion
        self.empty_regions.add_region(GenomicRegion(start=1, end=1000, chromosome='chr1'))
        self.empty_regions.flush()
        assert self.empty_regions[0].start == 1
        assert self.empty_regions[0].end == 1000
        assert self.empty_regions[0].chromosome == 'chr1'

        # dict
        self.empty_regions.add_region({'start': 1001, 'end': 2000, 'chromosome': 'chr1'})
        self.empty_regions.flush()
        assert self.empty_regions[1].start == 1001
        assert self.empty_regions[1].end == 2000
        assert self.empty_regions[1].chromosome == 'chr1'

        # list
        self.empty_regions.add_region(['chr1', 2001, 3000])
        self.empty_regions.flush()
        assert self.empty_regions[2].start == 2001
        assert self.empty_regions[2].end == 3000
        assert self.empty_regions[2].chromosome == 'chr1'

