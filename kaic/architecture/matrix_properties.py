from __future__ import division
from kaic.architecture.architecture import TableArchitecturalFeature, calculateondemand
from kaic.data.genomic import GenomicRegion
import tables as t
from collections import defaultdict
import numpy as np


class ExpectedContacts(TableArchitecturalFeature):
    def __init__(self, hic, file_name=None, mode='a', tmpdir=None, smooth=True, min_reads=400,
                 regions=None, _table_name='expected_contacts'):
        if isinstance(hic, str):
            file_name = hic
            hic = None

        TableArchitecturalFeature.__init__(self, _table_name,
                                           {'distance': t.Int64Col(), 'intra': t.Float32Col(),
                                            'contacts': t.Float32Col(), 'pixels': t.Float32Col()},
                                           file_name=file_name, mode=mode, tmpdir=tmpdir)

        self.hic = hic
        self.smooth = smooth
        self.min_reads = min_reads
        self.regions = regions

    def _calculate(self):
        """
        Get intra- and inter-chromosomal expected contact counts.
        """

        # extract mappable regions from Hic object
        marginals = self.hic.marginals()
        regions = []
        region_ix = set()
        if self.regions is None:
            for region in self.hic.regions(lazy=False):
                if marginals[region.ix] > 0:
                    regions.append(region)
                    region_ix.add(region.ix)
        else:
            if isinstance(self.regions, str) or isinstance(self.regions, GenomicRegion):
                self.regions = [self.regions]

            for region in self.regions:
                if isinstance(region, str):
                    region = GenomicRegion.from_string(region)

                for r in self.hic.subset(region, lazy=False):
                    if marginals[r.ix] > 0:
                        regions.append(r)
                        region_ix.add(r.ix)

        # count the number of regions per chromosome
        # and find the first and last region in each chromosome
        regions_by_chromosome = defaultdict(int)
        min_region_by_chromosome = dict()
        max_region_by_chromosome = dict()
        for region in regions:
            if (region.chromosome not in max_region_by_chromosome or
                    max_region_by_chromosome[region.chromosome] < region.ix):
                max_region_by_chromosome[region.chromosome] = region.ix

            if (region.chromosome not in min_region_by_chromosome or
                    min_region_by_chromosome[region.chromosome] > region.ix):
                min_region_by_chromosome[region.chromosome] = region.ix

            regions_by_chromosome[region.chromosome] += 1

        # find the largest distance between two regions
        # in the entire intra-chromosomal genome
        max_distance = 0
        for chromosome in min_region_by_chromosome:
            max_distance = max(max_distance, max_region_by_chromosome[chromosome]-min_region_by_chromosome[chromosome])

        # get the number of pixels at a given bin distance
        pixels_by_distance = np.zeros(max_distance + 1)
        for chromosome, n in regions_by_chromosome.iteritems():
            for distance in xrange(0, n):
                pixels_by_distance[distance] += n-distance

        # build a reverse-lookup chromosome map to quickly
        # determine if an edge is intra-chromosomal
        chromosome_map = dict()
        for i, chromosome in enumerate(self.hic.chromosomes()):
            chromosome_map[chromosome] = i

        chromosomes = np.zeros(len(self.hic.regions), dtype=int)
        for i, region in enumerate(self.hic.regions(lazy=True)):
            chromosomes[i] = chromosome_map[region.chromosome]

        # get the number of reads at a given bin distance
        reads_by_distance = np.zeros(max_distance + 1)
        inter_observed = 0
        for edge in self.hic.edges(lazy=True):
            source = edge.source
            sink = edge.sink
            # skip excluded regions
            if source not in region_ix or sink not in region_ix:
                continue
            # only intra-chromosomal distances
            if chromosomes[edge.source] == chromosomes[edge.sink]:
                reads_by_distance[edge.sink - edge.source] += edge.weight
            else:
                inter_observed += edge.weight

        with PossibleContacts(self.hic, regions=self.regions) as pc:
            intra_possible, inter_possible = pc.intra_possible(), pc.inter_possible()

        try:
            inter_expected = inter_observed/inter_possible
        except ZeroDivisionError:
            inter_expected = 0

        self.data('distance', np.arange(len(reads_by_distance))*self.hic.bin_size())

        # return here if smoothing not requested
        if not self.smooth:
            intra_expected = reads_by_distance/pixels_by_distance

            self.data('intra', intra_expected)
            self.data('contacts', reads_by_distance)
            self.data('pixels', pixels_by_distance)
            self._table.attrs['inter'] = inter_expected
            return

        # smoothing
        smoothed_reads_by_distance = np.zeros(max_distance + 1)
        smoothed_pixels_by_distance = np.zeros(max_distance + 1)
        for i in xrange(len(reads_by_distance)):
            smoothed_reads = reads_by_distance[i]
            smoothed_pixels = pixels_by_distance[i]
            window_size = 0
            can_extend = True
            # smooth to a minimum number of reads per distance
            while smoothed_reads < self.min_reads and can_extend:
                window_size += 1
                can_extend = False
                # check if we can increase the window to the left
                if i - window_size >= 0:
                    smoothed_reads += reads_by_distance[i-window_size]
                    smoothed_pixels += pixels_by_distance[i-window_size]
                    can_extend = True
                # check if we can increase the window to the right
                if i + window_size < len(reads_by_distance):
                    smoothed_reads += reads_by_distance[i+window_size]
                    smoothed_pixels += pixels_by_distance[i+window_size]
                    can_extend = True
            smoothed_reads_by_distance[i] = smoothed_reads
            smoothed_pixels_by_distance[i] = smoothed_pixels

        intra_expected = smoothed_reads_by_distance/smoothed_pixels_by_distance

        self.data('intra', intra_expected)
        self.data('contacts', smoothed_reads_by_distance)
        self.data('pixels', smoothed_pixels_by_distance)
        self._table.attrs['inter'] = inter_expected

    @calculateondemand
    def intra_expected(self):
        return self[:, 'intra']

    @calculateondemand
    def inter_expected(self):
        return self._table.attrs['inter']

    @calculateondemand
    def distance(self):
        return self[:, 'distance']

    @calculateondemand
    def intra_contacts(self):
        return self[:, 'contacts']

    @calculateondemand
    def intra_pixels(self):
        return self[:, 'pixels']


class PossibleContacts(TableArchitecturalFeature):
    def __init__(self, hic, file_name=None, mode='a', tmpdir=None, regions=None,
                 _table_name='expected_contacts'):
        if isinstance(hic, str):
            file_name = hic
            hic = None

        TableArchitecturalFeature.__init__(self, _table_name,
                                           {'intra': t.Int32Col(), 'inter': t.Int32Col()},
                                           file_name=file_name, mode=mode, tmpdir=tmpdir)

        self.hic = hic
        self.regions = regions

    def _calculate(self):
        marginals = self.hic.marginals()

        mappable = defaultdict(int)
        if self.regions is None:
            for r in self.hic.regions(lazy=True):
                if marginals[r.ix] > 0:
                    mappable[r.chromosome] += 1
        else:
            if isinstance(self.regions, str) or isinstance(self.regions, GenomicRegion):
                self.regions = [self.regions]

            for region in self.regions:
                if isinstance(region, str):
                    region = GenomicRegion.from_string(region)

                for r in self.hic.subset(region, lazy=True):
                    if marginals[r.ix] > 0:
                        mappable[r.chromosome] += 1

        # calculate possible combinations
        intra_possible = 0
        inter_possible = 0
        chromosomes = mappable.keys()
        for i in xrange(len(chromosomes)):
            chromosome1 = chromosomes[i]
            n1 = mappable[chromosome1]
            intra_possible += n1**2/2 + n1/2
            for j in xrange(i+1, len(chromosomes)):
                chromosome2 = chromosomes[j]
                n2 = mappable[chromosome2]
                inter_possible += n1*n2

        self.data('intra', [intra_possible])
        self.data('inter', [inter_possible])

    @calculateondemand
    def intra_possible(self):
        return self[0, 'intra']

    @calculateondemand
    def inter_possible(self):
        return self[0, 'inter']
