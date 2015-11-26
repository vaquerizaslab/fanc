from abc import abstractmethod, ABCMeta
import numpy as np
import logging
from collections import defaultdict


class PeakCaller(object):

    __metaclass__ = ABCMeta

    def __init__(self):
        pass

    @abstractmethod
    def call_peaks(self, hic):
        pass


class RaoPeakCaller(PeakCaller):
    def __init__(self, hic, fdr_intra=None, fdr_inter=None, p=None, w_init=None, min_locus_dist=3):
        self.fdr_intra = fdr_intra
        self.fdr_inter = fdr_inter
        self.p = p
        self.w_init = w_init
        self.min_locus_dist = min_locus_dist
        super(RaoPeakCaller, self).__init__()

    def chromosome_map(self, hic):
        # make a quick-lookup chromosome map
        chromosome_counter = 0
        chromosome_map = dict()
        for chromosome in hic.chromosomes():
            chromosome_map[chromosome] = chromosome_counter
            chromosome_counter += 1

        chromosomes = np.zeros(len(hic.regions()), dtype=int)
        for i, region in enumerate(hic.regions(lazy=True)):
            chromosomes[i] = chromosome_map[region.chromosome]
        return chromosomes

    def get_expected(self, hic, smooth=True, min_smoothed_reads=400, _mappable=None, _chromosomes=None):
        if _mappable is None:
            _mappable = hic.get_mappable_regions()

        max_distance = max(_mappable.values())
        # get the number of pixels at a given bin distance
        pixels_by_distance = np.zeros(max_distance + 1)
        for chromosome, n in _mappable.iteritems():
            current_pixels = n + 1
            for distance in xrange(0, n + 1):
                pixels_by_distance[distance] += current_pixels
                current_pixels -= 1

        # get the number of reads at a given bin distance
        if _chromosomes is None:
            _chromosomes = self.chromosome_map(hic)
        reads_by_distance = np.zeros(max_distance + 1)
        inter_observed = 0
        for edge in hic.edges(lazy=True):
            # only intra-chromosomal distances
            if _chromosomes[edge.source] == _chromosomes[edge.sink]:
                reads_by_distance[edge.sink - edge.source] += edge.weight
            else:
                inter_observed += edge.weight

        intra_possible, inter_possible = hic.get_possible_contacts(_mappable=_mappable)

        # return here if smoothing not requested
        if not smooth:
            return reads_by_distance/pixels_by_distance, inter_observed/inter_possible

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
        return smoothed_reads_by_distance/smoothed_pixels_by_distance, inter_observed/inter_possible

    # sum of reads in lower-left neighborhood
    @staticmethod
    def ll_sum(m, i, j, w=1, p=0):
        i_max, j_max = m.shape

        sum1 = 0
        for a in xrange(max(0, i+1), min(i_max, i+w+1)):
            for b in xrange(max(0, j-w), min(j_max, j)):
                sum1 += m[a, b]

        sum2 = 0
        for a in xrange(max(0, i+1), min(i_max, i+p+1)):
            for b in xrange(max(0, j-p), min(j_max, j)):
                sum2 += m[a, b]

        return sum1 - sum2

    # lower-left neighborhood
    @staticmethod
    def e_ll(m, i, j, e, w=1, p=0):
        i_max, j_max = m.shape

        sum1 = 0
        for a in xrange(max(0, i+1), min(i+w+1, i_max)):
            for b in xrange(max(0, j-w), min(j_max, j)):
                sum1 += m[a, b]

        sum2 = 0
        for a in xrange(max(0, i+1), min(i_max, i+p+1)):
            for b in xrange(max(0, j-p), min(j_max, j)):
                sum2 += m[a, b]

        sum3 = 0
        for a in xrange(max(0, i+1), min(i_max, i+w+1)):
            for b in xrange(max(0, j-w), min(j_max, j)):
                sum3 += e(a, b)

        sum4 = 0
        for a in xrange(max(0, i+1), min(i_max, i+p+1)):
            for b in xrange(max(0, j-p), min(j_max, j)):
                sum4 += e(a, b)

        return (sum1-sum2)/(sum3-sum4)*e(i, j)

    # vertical neighborhood
    @staticmethod
    def e_v(m, i, j, e, w=1, p=0):
        i_max, j_max = m.shape

        sum1 = 0
        for a in xrange(max(0, i-w), min(i_max, i+w)):
            for b in xrange(max(0, j-1), min(j_max, j+2)):
                sum1 += m[a, b]

        sum2 = 0
        for a in xrange(max(0, i-p), min(i_max, i+p+1)):
            for b in xrange(max(0, j-1), min(j_max, j+2)):
                sum2 += m[a, b]

        sum3 = 0
        for a in xrange(max(0, i-w), min(i_max, i+w)):
            for b in xrange(max(0, j-1), min(j_max, j+2)):
                sum3 += e(a, b)

        sum4 = 0
        for a in xrange(max(0, i-p), min(i_max, i+p+1)):
            for b in xrange(max(0, j-1), min(j_max, j+2)):
                sum4 += e(a, b)

        return (sum1-sum2)/(sum3-sum4)*e(i, j)

    # horizontal neighborhood
    @staticmethod
    def e_h(m, i, j, e, w=1, p=0):
        i_max, j_max = m.shape

        sum1 = 0
        for a in xrange(max(0, i-1), min(i_max, i+2)):
            for b in xrange(max(0, j-w), min(j_max, j+w)):
                sum1 += m[a, b]

        sum2 = 0
        for a in xrange(max(0, i-1), min(i_max, i+2)):
            for b in xrange(max(0, j-p), min(j_max, j+p+1)):
                sum2 += m[a, b]

        sum3 = 0
        for a in xrange(max(0, i-1), min(i_max, i+2)):
            for b in xrange(max(0, j-w), min(j_max, j+w)):
                sum3 += e(a, b)

        sum4 = 0
        for a in xrange(max(0, i-1), min(i_max, i+2)):
            for b in xrange(max(0, j-p), min(j_max, j+p+1)):
                sum4 += e(a, b)

        return (sum1-sum2)/(sum3-sum4)*e(i, j)

    # donut neighborhood
    @staticmethod
    def e_d(m, i, j, e, w=1, p=0):
        i_max, j_max = m.shape

        top_sum1 = 0
        for a in xrange(max(0, i-w), min(i_max, i+w+1)):
            for b in xrange(max(0, j-w), min(j_max, j+w+1)):
                top_sum1 += m[a, b]

        top_sum2 = 0
        for a in xrange(max(0, i-p), min(i_max, i+p+1)):
            for b in xrange(max(0, j-p), min(j_max, j+p+1)):
                top_sum2 += m[a, b]

        top_sum3 = 0
        for a in xrange(max(0, i-w), min(i_max, i-p)):
            top_sum3 += m[a, j]

        top_sum4 = 0
        for a in xrange(max(0, i+p+1), min(i_max, i+w+1)):
            top_sum4 += m[a, j]

        top_sum5 = 0
        for b in xrange(max(0, j-w), min(j_max, j-p)):
            top_sum5 += m[i, b]

        top_sum6 = 0
        for b in xrange(max(0, j+p+1), min(j_max, j+w+1)):
            top_sum6 += m[i, b]

        bottom_sum1 = 0
        for a in xrange(max(0, i-w), min(i_max, i+w+1)):
            for b in xrange(max(0, j-w), min(j_max, j+w+1)):
                bottom_sum1 += e(a, b)

        bottom_sum2 = 0
        for a in xrange(max(0, i-p), min(i_max, i+p+1)):
            for b in xrange(max(0, j-p), min(j_max, j+p+1)):
                bottom_sum2 += e(a, b)

        bottom_sum3 = 0
        for a in xrange(max(0, i-w), min(i_max, i-p)):
            bottom_sum3 += e(a, j)

        bottom_sum4 = 0
        for a in xrange(max(0, i+p+1), min(i_max, i+w+1)):
            bottom_sum4 += e(a, j)

        bottom_sum5 = 0
        for b in xrange(max(0, j-w), min(j_max, j-p)):
            bottom_sum5 += e(i, b)

        bottom_sum6 = 0
        for b in xrange(max(0, j+p+1), min(j_max, j+w+1)):
            bottom_sum6 += e(i, b)

        return (top_sum1-top_sum2-top_sum3-top_sum4-top_sum5-top_sum6) / \
            (bottom_sum1-bottom_sum2-bottom_sum3-bottom_sum4-bottom_sum5-bottom_sum6) * \
            e(i, j)

    @staticmethod
    def e_all(m, i, j, e, w=1, p=0, min_locus_dist=3, min_ll_reads=16, max_w=20):
        # do not examine loci closer than p+3
        if abs(i-j) <= p + min_locus_dist:
            return None, None, None, None

        # assure minimum number of reads
        while RaoPeakCaller.ll_sum(m, i, j, w=w, p=p) < min_ll_reads and w < max_w:
            w += 1

        # fail-safe
        if w >= 20:
            return None, None, None, None

        # neighborhood expected values
        e_ll = RaoPeakCaller.e_ll(m, i, j, e, w=w, p=p)
        e_h = RaoPeakCaller.e_h(m, i, j, e, w=w, p=p)
        e_v = RaoPeakCaller.e_v(m, i, j, e, w=w, p=p)
        e_d = RaoPeakCaller.e_d(m, i, j, e, w=w, p=p)

        return e_ll, e_h, e_v, e_d

    def call_peaks(self, hic):
        logging.info("Calculating expected values...")
        intra_expected, inter_expected = self.get_expected(hic, smooth=True)
        logging.info("Done.")

        # define helper functions for
        # expected number of contacts
        def e_intra(i, j):
            d = abs(i-j)
            return intra_expected[d]

        def e_inter(i, j):
            return inter_expected

        # initialize peak parameters
        p = self.p
        w_init = self.w_init

        # empirically choose p and w_init
        # if not already chosen
        if p is None or w_init is None:
            bin_size = hic.bin_size()
            if bin_size > 25000:
                p = 1 if p is not None else p
                w_init = 3 if w_init is not None else w_init
            else:
                p = int(24999/bin_size) if p is not None else p
                w_init = int(round(25000/bin_size) + 2) if w_init is not None else w_init

        c = hic.bias_vector()
        oe_diff = defaultdict(list)

        logging.info("Initial parameter values: p=%d, w=%d" % (p, w_init))

        # start processing chromosome pairs
        chromosomes = list(hic.chromosomes())
        for i_chr, chromosome1 in enumerate(chromosomes):
            for j_chr in xrange(i_chr, len(chromosomes)):
                chromosome2 = chromosomes[j_chr]

                m = hic[chromosome1, chromosome2]

                if chromosome1 == chromosome2:
                    for i in xrange(m.shape[0]):
                        i_region = m.row_regions[i].ix
                        for j in xrange(i, m.shape[1]):
                            j_region = m.col_regions[j].ix

                            observed = int(m[i, j]*c[i_region]*c[j_region])
                            e_ll, e_h, e_v, e_d = RaoPeakCaller.e_all(m, i, j, e_intra, w=w_init, p=p)


        def calculate_neighborhood(e_local=RaoPeakCaller.e_d, find_inter_peaks=True):
            ij = []
            exp = []
            obs = []
            ij_inter = []
            exp_inter = []
            obs_inter = []

            chromosomes = list(hic.chromosomes())

            # intra-chromosomal
            for chromosome in chromosomes:
                m = hic[chromosome, chromosome]

                for i in xrange(m.shape[0]):
                    i_region = m.row_regions[i].ix
                    for j in xrange(i, m.shape[1]):
                        j_region = m.col_regions[j].ix

                        # do not examine loci closer than p+3
                        if abs(i-j) > p + self.min_locus_dist:
                            w = w_init

                            # assure minimum number of reads
                            while ll_sum(m, i, j, w=w, p=p) < 16 and w < 20:
                                w += 1

                            # fail-safe
                            if w >= 20:
                                continue

                            # reproduce original value
                            observed = int(m[i, j]*c[i_region]*c[j_region])

                            # neighborhood expected values
                            ij.append([i_region, j_region])
                            expected = e_local(m, i, j, w=w, p=p)
                            exp.append(expected)
                            obs.append(observed)

                            oe_diff[chromosome].append(observed-expected)

            # inter-chromosomal
            if find_inter_peaks:
                for chromosome_i in xrange(len(chromosomes)):
                    chromosome1 = chromosomes[chromosome_i]
                    for chromosome_j in xrange(chromosome_i+1, len(chromosomes)):
                        chromosome2 = chromosomes[chromosome_j]

                        m = hic[chromosome1, chromosome2]

                        for i in xrange(m.shape[0]):
                            i_region = m.row_regions[i].ix
                            for j in xrange(i, m.shape[1]):
                                j_region = m.col_regions[j].ix

                                w = w_init

                                # assure minimum number of reads
                                while ll_sum(m, i, j, w=w, p=p) < 16 and w < 20:
                                    w += 1

                                # fail-safe
                                if w >= 20:
                                    continue

                                # reproduce original value
                                observed = int(m[i, j]*c[i_region]*c[j_region])

                                # neighborhood expected values
                                ij_inter.append([i_region, j_region])
                                expected = e_local(m, i, j, w=w, p=p, inter=True)
                                exp_inter.append(expected)
                                obs_inter.append(observed)

            return ij, exp, obs, ij_inter, exp_inter, obs_inter

