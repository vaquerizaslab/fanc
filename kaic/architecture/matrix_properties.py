from kaic.architecture.architecture import TableArchitecturalFeature, calculateondemand
from kaic.data.genomic import GenomicRegion
import tables as t
from collections import defaultdict
import numpy as np


class ExpectedContacts(TableArchitecturalFeature):
    def __init__(self, hic, file_name=None, mode='a', tmpdir=None, smooth=True, min_reads=400):
        if isinstance(hic, str):
            file_name = hic
            hic = None

        TableArchitecturalFeature.__init__(self, 'expected_contacts',
                                           {'distance': t.Int64Col(), 'intra': t.Float32Col(),
                                            'contacts': t.Float32Col(), 'pixels': t.Float32Col()},
                                           file_name=file_name, mode=mode, tmpdir=tmpdir)

        self.hic = hic
        self.smooth = smooth
        self.min_reads = min_reads

    def _calculate(self, smooth=True, min_smoothed_reads=400):
        """
        Get intra- and inter-chromosomal expected contact counts.

        :param smooth: Smoothe intra-chromosomal expected counts
        :param min_smoothed_reads: Minimum number of reads/counts per
                                   expected value
        :return: (np.array, float), where the first argument is a numpy
                 array with expected intra-chromosomal counts at any given
                 distance from the diagonal of the Hi-C matrix (loci distance)
                 and float is the average number of inter-chromosomal reads
                 per contact
        """
        regions_by_chromosome = defaultdict(int)
        for region in self.hic.regions(lazy=True):
            regions_by_chromosome[region.chromosome] += 1

        max_distance = max(regions_by_chromosome.values())
        # get the number of pixels at a given bin distance
        pixels_by_distance = np.zeros(max_distance + 1)
        for chromosome, n in regions_by_chromosome.iteritems():
            current_pixels = n + 1
            for distance in xrange(0, n + 1):
                pixels_by_distance[distance] += current_pixels
                current_pixels -= 1

        # get the number of reads at a given bin distance
        chromosome_map = dict()
        for i, chromosome in enumerate(self.hic.chromosomes()):
            chromosome_map[chromosome] = i

        chromosomes = np.zeros(len(self.hic.regions()), dtype=int)
        for i, region in enumerate(self.hic.regions(lazy=True)):
            chromosomes[i] = chromosome_map[region.chromosome]

        reads_by_distance = np.zeros(max_distance + 1)
        inter_observed = 0
        for edge in self.hic.edges(lazy=True):
            # only intra-chromosomal distances
            if chromosomes[edge.source] == chromosomes[edge.sink]:
                reads_by_distance[edge.sink - edge.source] += edge.weight
            else:
                inter_observed += edge.weight

        pc = PossibleContacts(self.hic)
        intra_possible, inter_possible = pc.intra_possible(), pc.inter_possible()
        pc.close()

        try:
            inter_expected = inter_observed/inter_possible
        except ZeroDivisionError:
            inter_expected = 0

        self.data('distance', np.arange(len(reads_by_distance))*self.hic.bin_size())

        # return here if smoothing not requested
        if not smooth:
            intra_expected = reads_by_distance/pixels_by_distance

            self.data('intra', intra_expected)
            self.data('contacts', reads_by_distance)
            self.data('pixels', pixels_by_distance)
            self._table.attrs['inter'] = inter_expected

        # smoothing
        smoothed_reads_by_distance = np.zeros(max_distance + 1)
        smoothed_pixels_by_distance = np.zeros(max_distance + 1)
        for i in xrange(len(reads_by_distance)):
            smoothed_reads = reads_by_distance[i]
            smoothed_pixels = pixels_by_distance[i]
            window_size = 0
            can_extend = True
            # smooth to a minimum number of reads per distance
            while smoothed_reads < min_smoothed_reads and can_extend:
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
    def __init__(self, hic, file_name=None, mode='a', tmpdir=None, regions=None):
        if isinstance(hic, str):
            file_name = hic
            hic = None

        TableArchitecturalFeature.__init__(self, 'possible_region_contacts',
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
