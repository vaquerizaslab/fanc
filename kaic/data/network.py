from __future__ import division
from abc import abstractmethod, ABCMeta
import numpy as np
from scipy.stats import poisson
import logging
from collections import defaultdict
import tables as t
from bisect import bisect_left
from functools import partial
from kaic.data.genomic import RegionsTable
from kaic.data.general import FileBased, MaskedTable, MaskFilter
import msgpack
import time

try:
    import gridmap
    _has_gridmap = True
except ImportError:
    import multiprocessing
    _has_gridmap = False


class PeakCaller(object):

    __metaclass__ = ABCMeta

    def __init__(self):
        pass

    @abstractmethod
    def call_peaks(self, hic):
        pass


class RaoPeakInfo(RegionsTable, FileBased):
    class PeakInformation(t.IsDescription):
        source = t.Int32Col(pos=0)
        sink = t.Int32Col(pos=1)
        observed = t.Int32Col(pos=2)
        e_ll = t.Float32Col(pos=3)
        e_h = t.Float32Col(pos=4)
        e_v = t.Float32Col(pos=5)
        e_d = t.Float32Col(pos=6)
        e_ll_chunk = t.Int32Col(pos=7)
        e_h_chunk = t.Int32Col(pos=8)
        e_v_chunk = t.Int32Col(pos=9)
        e_d_chunk = t.Int32Col(pos=10)
        fdr_ll = t.Float32Col(pos=11)
        fdr_h = t.Float32Col(pos=12)
        fdr_v = t.Float32Col(pos=13)
        fdr_d = t.Float32Col(pos=14)

    def __init__(self, file_name, mode='a', _table_name_regions='regions',
                 _table_name_peaks_intra='peak_info_intra',
                 _table_name_peaks_inter='peak_info_inter'):
        FileBased.__init__(self, file_name, mode=mode)
        RegionsTable.__init__(self, file_name=self.file, _table_name_regions=_table_name_regions)

        if _table_name_peaks_intra in self.file.root:
            self._edges = self.file.get_node('/', _table_name_peaks_intra)
        else:
            self.peak_table_intra = MaskedTable(self.file.root, _table_name_peaks_intra, RaoPeakInfo.PeakInformation)

        if _table_name_peaks_inter in self.file.root:
            self._edges = self.file.get_node('/', _table_name_peaks_inter)
        else:
            self.peak_table_inter = MaskedTable(self.file.root, _table_name_peaks_inter, RaoPeakInfo.PeakInformation)

    def _row2peak(self, row, lazy=False, auto_update=True):
        if not lazy:
            peak = RaoPeak(row['source'], row['sink'], row['observed'],
                           e_ll=row['e_ll'], e_h=row['e_h'], e_v=row['e_v'], e_d=row['e_d'],
                           e_ll_chunk=row['e_ll_chunk'], e_h_chunk=row['e_h_chunk'],
                           e_v_chunk=row['e_v_chunk'], e_d_chunk=row['e_d_chunk'],
                           fdr_ll=row['fdr_ll'], fdr_h=row['fdr_h'], fdr_v=row['fdr_v'], fdr_d=row['fdr_d'])
        else:
            peak = LazyRaoPeak(row, auto_update=auto_update)
        return peak

    def intra_peaks(self, lazy=False, auto_update=True):
        it = self.peak_table_intra.iterrows()
        return (self._row2peak(row, lazy=lazy, auto_update=auto_update) for row in it)

    def filter_intra(self, peak_filter, queue=False, log_progress=False):
        """
        Filter edges in this object by using a :class:`~PeakFilter`.

        :param peak_filter: Class implementing :class:`~PeakFilter`.
                            Must override valid_peak method, ideally sets mask parameter
                            during initialization.
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        :param log_progress: If true, process iterating through all edges
                             will be continuously reported.
        """
        if not queue:
            self.peak_table_intra.filter(peak_filter, _logging=log_progress)
        else:
            self.peak_table_intra.queue_filter(peak_filter)


class RaoPeak(object):
    def __init__(self, source, sink, observed=0, e_ll=0, e_h=0, e_v=0, e_d=0,
                 e_ll_chunk=0, e_h_chunk=0, e_v_chunk=0, e_d_chunk=0,
                 fdr_ll=1, fdr_h=1, fdr_v=1, fdr_d=1):
        self.source = source
        self.sink = sink
        self.observed = observed
        self.e_ll = e_ll
        self.e_h = e_h
        self.e_v = e_v
        self.e_d = e_d
        self.e_ll_chunk = e_ll_chunk
        self.e_h_chunk = e_h_chunk
        self.e_v_chunk = e_v_chunk
        self.e_d_chunk = e_d_chunk
        self.fdr_ll = fdr_ll
        self.fdr_h = fdr_h
        self.fdr_v = fdr_v
        self.fdr_d = fdr_d


class LazyRaoPeak(object):
    def __init__(self, row, auto_update=True):
        self.row = row
        self.auto_update = auto_update

    def _set_item(self, item, value):
        self.row[item] = value
        if self.auto_update:
            self.row.update()

    def __getattr__(self, item):
        if item == 'row' or item == 'auto_update':
            return super(LazyRaoPeak, self).__getattr__(item)
        return self.row[item]

    def __setattr__(self, key, value):
        if key == 'row' or key == 'auto_update':
            super(LazyRaoPeak, self).__setattr__(key, value)
        else:
            self.row[key] = value
            if self.auto_update:
                self.row.update()

    def update(self):
        self.row.update()


class PeakFilter(MaskFilter):
    """
    Abstract class that provides filtering functionality for the
    peaks in a :class:`~RaoPeakInfo` object.

    Extends MaskFilter and overrides valid(self, row) to make
    :class:`~RaoPeakInfo` filtering more "natural".

    To create custom filters for the :class:`~RapPeakInfo` object, extend this
    class and override the valid_peak(self, peak) method.
    valid_peak should return False for a specific :class:`~RaoPeak` object
    if the object is supposed to be filtered/masked and True
    otherwise. See :class:`~DiagonalFilter` for an example.

    Pass a custom filter to the :func:`~Hic.filter` method in :class:`~Hic`
    to apply it.
    """

    __metaclass__ = ABCMeta

    def __init__(self, mask=None):
        """
        Initialize PeakFilter.

        :param mask: The Mask object that should be used to mask
                     filtered :class:`~RaoPeak` objects. If None the default
                     Mask will be used.
        """
        super(PeakFilter, self).__init__(mask)

    @abstractmethod
    def valid_peak(self, peak):
        """
        Determine if a :class:`~RaoPeak` object is valid or should
        be filtered.

        When implementing custom PeakFilter this method must be
        overridden. It should return False for :class:`~RaoPeak` objects that
        are to be fitered and True otherwise.

        Internally, the :class:`~RaoPeakInfo` object will iterate over all RaoPeak
        instances to determine their validity on an individual
        basis.

        :param peak: A :class:`~RaoPeak` object
        :return: True if :class:`~PeakFilter` is valid, False otherwise
        """
        pass

    def valid(self, row):
        """
        Map valid_peak to MaskFilter.valid(self, row).

        :param row: A pytables Table row.
        :return: The boolean value returned by valid_edge.
        """
        peak = LazyRaoPeak(row)
        return self.valid_peak(peak)


class RaoPeakCaller(PeakCaller):
    """
    Class that calls peaks the same way Rao et al. (2014) propose.

    Every pixel in a Hi-C matrix is evaluated based on its local
    "neighborhood", i.e. the pixel's observed value is compared
    to the expected value calculated from all pixels in its close
    surroundings.

    If a pixel is significantly enriched with respect to all
    investigated neighborhoods, it is assumed to be a "peak".

    Four neighborhood types are calculated:

    - "donut": Pixels surrounding the investigated pixel in a
      certain distance range
    - "lower-left" Pixels to the "lower-left" of a given pixel
    - "horizontal": Pixels left and right of a given pixel
    - "vertical": Pixels above and below a given pixel

    While the first neighborhood type most generally calculates
    enrichment of local background, the other types of
    neighborhoods serve mostly to exclude false-positive results
    from non-peak structures, such as TAD boundaries.

    The :class:`~RaoPeakCaller` is initialized with the peak
    calling parameters and run using :func:`~RaoPeakCaller.call_peaks`.
    """

    def __init__(self, p=None, w_init=None, min_locus_dist=3, max_w=20, min_ll_reads=16,
                 observed_cutoff=1, e_ll_cutoff=1.0, e_h_cutoff=1.0, e_v_cutoff=1.0,
                 e_d_cutoff=1.0, process_inter=False, n_processes=4,
                 batch_size=500000):
        """
        Initialize the peak caller with parameters.

        :param p:
        :param w_init:
        :param min_locus_dist:
        :param max_w:
        :param min_ll_reads:
        :param observed_cutoff:
        :param e_ll_cutoff:
        :param e_h_cutoff:
        :param e_v_cutoff:
        :param e_d_cutoff:
        :param process_inter:
        :param n_processes:
        :param batch_size:
        :return:
        """
        self.p = p
        self.w_init = w_init
        self.min_locus_dist = min_locus_dist
        self.max_w = max_w
        self.min_ll_reads = min_ll_reads
        self.observed_cutoff = observed_cutoff
        self.e_ll_cutoff = e_ll_cutoff
        self.e_h_cutoff = e_h_cutoff
        self.e_v_cutoff = e_v_cutoff
        self.e_d_cutoff = e_d_cutoff
        self.process_inter = process_inter
        self.n_processes = n_processes
        self.batch_size = batch_size
        super(RaoPeakCaller, self).__init__()

    @staticmethod
    def chromosome_map(hic):
        # make a quick-lookup chromosome map
        chromosome_map = dict()
        for i, chromosome in enumerate(hic.chromosomes()):
            chromosome_map[chromosome] = i

        chromosomes = np.zeros(len(hic.regions()), dtype=int)
        for i, region in enumerate(hic.regions(lazy=True)):
            chromosomes[i] = chromosome_map[region.chromosome]
        return chromosomes

    @staticmethod
    def get_expected(hic, smooth=True, min_smoothed_reads=400, _mappable=None, _chromosomes=None):
        if _mappable is None:
            _mappable = hic.mappable_regions()

        regions_by_chromosome = defaultdict(int)
        for region in hic.regions(lazy=True):
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
        if _chromosomes is None:
            _chromosomes = RaoPeakCaller.chromosome_map(hic)

        reads_by_distance = np.zeros(max_distance + 1)
        inter_observed = 0
        for edge in hic.edges(lazy=True):
            # only intra-chromosomal distances
            if _chromosomes[edge.source] == _chromosomes[edge.sink]:
                reads_by_distance[edge.sink - edge.source] += edge.weight
            else:
                inter_observed += edge.weight

        intra_possible, inter_possible = hic.possible_contacts(_mappable=_mappable)

        try:
            inter_expected = inter_observed/inter_possible
        except ZeroDivisionError:
            inter_expected = 0

        # return here if smoothing not requested
        if not smooth:
            return reads_by_distance/pixels_by_distance, inter_expected

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
        return smoothed_reads_by_distance/smoothed_pixels_by_distance, inter_expected

    # sum of reads in lower-left neighborhood
    @staticmethod
    def ll_sum(m, i, j, w=1, p=0):
        i_max, j_max = m.shape

        sum1 = np.sum(m[max(0, i+1):min(i_max, i+w+1), max(0, j-w):min(j_max, j)])
        sum2 = np.sum(m[max(0, i+1):min(i_max, i+p+1), max(0, j-p):min(j_max, j)])

        return sum1 - sum2

    # lower-left neighborhood
    @staticmethod
    def e_ll(m, i, j, e, w=1, p=0):
        i_max, j_max = m.shape

        # dividend
        sum1 = np.sum(m[max(0, i+1):min(i+w+1, i_max), max(0, j-w):min(j_max, j)])
        sum2 = np.sum(m[max(0, i+1):min(i_max, i+p+1), max(0, j-p):min(j_max, j)])

        if sum1-sum2 == 0:
            return 0

        # divisor
        sum3 = 0
        for a in xrange(max(0, i+1), min(i_max, i+w+1)):
            for b in xrange(max(0, j-w), min(j_max, j)):
                sum3 += e(a, b)

        sum4 = 0
        for a in xrange(max(0, i+1), min(i_max, i+p+1)):
            for b in xrange(max(0, j-p), min(j_max, j)):
                sum4 += e(a, b)

        return (sum1-sum2)/(sum3-sum4)*e(i, j)

    # horizontal neighborhood
    @staticmethod
    def e_h(m, i, j, e, w=1, p=0):
        i_max, j_max = m.shape

        # dividend
        sum1 = np.sum(m[max(0, i-1):min(i_max, i+2), max(0, j-w):min(j_max, j+w+1)])
        sum2 = np.sum(m[max(0, i-1):min(i_max, i+2), max(0, j-p):min(j_max, j+p+1)])

        if sum1-sum2 == 0:
            return 0

        # divisor
        sum3 = 0
        for a in xrange(max(0, i-1), min(i_max, i+2)):
            for b in xrange(max(0, j-w), min(j_max, j+w+1)):
                sum3 += e(a, b)

        sum4 = 0
        for a in xrange(max(0, i-1), min(i_max, i+2)):
            for b in xrange(max(0, j-p), min(j_max, j+p+1)):
                sum4 += e(a, b)

        return (sum1-sum2)/(sum3-sum4)*e(i, j)

    # vertical neighborhood
    @staticmethod
    def e_v(m, i, j, e, w=1, p=0):
        i_max, j_max = m.shape

        # dividend
        sum1 = np.sum(m[max(0, i-w):min(i_max, i+w+1), max(0, j-1):min(j_max, j+2)])
        sum2 = np.sum(m[max(0, i-p):min(i_max, i+p+1), max(0, j-1):min(j_max, j+2)])

        if sum1-sum2 == 0:
            return 0

        # divisor
        sum3 = 0
        for a in xrange(max(0, i-w), min(i_max, i+w+1)):
            for b in xrange(max(0, j-1), min(j_max, j+2)):
                sum3 += e(a, b)

        sum4 = 0
        for a in xrange(max(0, i-p), min(i_max, i+p+1)):
            for b in xrange(max(0, j-1), min(j_max, j+2)):
                sum4 += e(a, b)

        return (sum1-sum2)/(sum3-sum4)*e(i, j)

    # donut neighborhood
    @staticmethod
    def e_d(m, i, j, e, w=1, p=0):
        i_max, j_max = m.shape

        top_sum1 = np.sum(m[max(0, i-w):min(i_max, i+w+1), max(0, j-w):min(j_max, j+w+1)])
        top_sum2 = np.sum(m[max(0, i-p):min(i_max, i+p+1), max(0, j-p):min(j_max, j+p+1)])
        top_sum3 = np.sum(m[max(0, i-w):min(i_max, i-p), j])
        top_sum4 = np.sum(m[max(0, i+p+1):min(i_max, i+w+1), j])
        top_sum5 = np.sum(m[i, max(0, j-w):min(j_max, j-p)])
        top_sum6 = np.sum(m[i, max(0, j+p+1):min(j_max, j+w+1)])

        top = (top_sum1-top_sum2-top_sum3-top_sum4-top_sum5-top_sum6)
        if top == 0:
            return 0

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

        return top / \
            (bottom_sum1-bottom_sum2-bottom_sum3-bottom_sum4-bottom_sum5-bottom_sum6) * \
            e(i, j)

    @staticmethod
    def e_all(m, i, j, e, w=1, p=0):

        if isinstance(e, int) or isinstance(e, float):
            def e_f(ix, jx):
                return e
        else:
            def e_f(ix, jx):
                return e[abs(ix-jx)]

        # neighborhood expected values
        e_ll = RaoPeakCaller.e_ll(m, i, j, e_f, w=w, p=p)
        e_h = RaoPeakCaller.e_h(m, i, j, e_f, w=w, p=p)
        e_v = RaoPeakCaller.e_v(m, i, j, e_f, w=w, p=p)
        e_d = RaoPeakCaller.e_d(m, i, j, e_f, w=w, p=p)

        e_ll = None if e_ll is None else float(e_ll)
        e_h = None if e_h is None else float(e_h)
        e_v = None if e_v is None else float(e_v)
        e_d = None if e_d is None else float(e_d)

        return e_ll, e_h, e_v, e_d

    @staticmethod
    def _lambda_chunks(max_expected, max_chunk_function=lambda x: 2**(x/3)):
        logging.info("Finding maximum expected values")

        e_max_chunk = 1
        e_exp = 0
        chunks_max_list = []
        while e_max_chunk < max_expected:
            chunks_max_list.append(e_max_chunk)
            e_exp += 1
            # increment power in steps of 1/3
            e_max_chunk = max_chunk_function(e_exp)
        # final chunk
        chunks_max_list.append(e_max_chunk)

        return chunks_max_list

    @staticmethod
    def _find_chunk(chunk_list, value):
        if value is None:
            return None
        return bisect_left(chunk_list, value)

    @staticmethod
    def _submatrix_indices(m, ij_pairs, w_max=20):
        if len(ij_pairs) == 0:
            return None, None

        # find boundaries including padding through w
        min_i, min_j, max_i, max_j = ij_pairs[0][0]-w_max, ij_pairs[0][1]-w_max, w_max+1, w_max+1
        for i, j in ij_pairs:
            min_i, min_j, max_i, max_j = (min(min_i, i-w_max), min(min_j, j-w_max),
                                          max(max_i, i+w_max+1), max(max_j, j+w_max+1))

        # ensure that we don't cross matrix boundaries
        min_i, min_j, max_i, max_j = (max(0, min_i), max(0, min_j),
                                      min(m.shape[0], max_i), min(m.shape[1], max_j))

        # extract sub-matrix
        m_sub = m[min_i:max_i, min_j:max_j]

        # convert ij_pairs
        ij_converted = []
        for i, j in ij_pairs:
            ij_converted.append((i-min_i, j-min_j))

        return m_sub, ij_converted

    def _process_jobs(self, jobs, peak_info, observed_chunk_distribution):
        if _has_gridmap:
            t1 = time.time()
            job_outputs = gridmap.process_jobs(jobs, max_processes=self.n_processes)

            t2 = time.time()

            row = peak_info.row
            for output in job_outputs:
                rv = output

                (region_pairs, observed_list, e_ll_list,
                 e_h_list, e_v_list, e_d_list, observed_chunk_distribution_part) = msgpack.loads(rv)

                for ix in xrange(len(region_pairs)):
                    source = region_pairs[ix][0]
                    sink = region_pairs[ix][1]
                    observed, observed_chunk = observed_list[ix]
                    e_ll, e_ll_chunk = e_ll_list[ix]
                    e_h, e_h_chunk = e_h_list[ix]
                    e_v, e_v_chunk = e_v_list[ix]
                    e_d, e_d_chunk = e_d_list[ix]

                    if (e_ll is not None and e_h is not None and
                            e_v is not None and e_d is not None):
                        row['source'] = source
                        row['sink'] = sink
                        row['observed'] = observed
                        row['e_ll'] = e_ll
                        row['e_h'] = e_h
                        row['e_v'] = e_v
                        row['e_d'] = e_d
                        row['e_ll_chunk'] = e_ll_chunk
                        row['e_h_chunk'] = e_h_chunk
                        row['e_v_chunk'] = e_v_chunk
                        row['e_d_chunk'] = e_d_chunk
                        row.append()

                for e_type in observed_chunk_distribution_part.iterkeys():
                    for chunk_ix in xrange(len(observed_chunk_distribution_part[e_type])):
                        for o in observed_chunk_distribution_part[e_type][chunk_ix].iterkeys():
                            observed_chunk_distribution[e_type][chunk_ix][o] += observed_chunk_distribution_part[e_type][chunk_ix][o]

            t3 = time.time()
        else:
            raise RuntimeError("gridmap not installed, multiprocessing not enabled")

    @staticmethod
    def _get_chunk_distribution_container(lambda_chunks):
        observed_chunk_distribution = {
            'll': [],
            'h': [],
            'v': [],
            'd': []
        }
        for _ in xrange(len(lambda_chunks)):
            observed_chunk_distribution['ll'].append(defaultdict(int))
            observed_chunk_distribution['h'].append(defaultdict(int))
            observed_chunk_distribution['v'].append(defaultdict(int))
            observed_chunk_distribution['d'].append(defaultdict(int))
        return observed_chunk_distribution

    @staticmethod
    def _get_fdr_cutoffs(lambda_chunks, observed_chunk_distribution):
        # determine all possible observed values

        fdr_cutoffs = dict()
        for e_type in observed_chunk_distribution:  # ll, h, v, d
            fdr_cutoffs[e_type] = []
            for chunk, max_e in enumerate(lambda_chunks):  # (0, 1), (1, 1.26), (2, 1.59), ...
                fdr_cutoffs[e_type].append(dict())
                poisson_e = poisson(max_e)

                observed_sum = 0
                for observed_count in observed_chunk_distribution[e_type][chunk].itervalues():
                    observed_sum += observed_count

                observed_distribution_integral_left = 0
                for observed in sorted(observed_chunk_distribution[e_type][chunk].keys()):
                    observed_count = observed_chunk_distribution[e_type][chunk][observed]
                    observed_distribution_integral_left += observed_count/observed_sum
                    observed_distribution_integral_right = 1 - observed_distribution_integral_left
                    poisson_integral_right = 1 - poisson_e.cdf(observed)

                    if observed_distribution_integral_right > 0:
                        integral_ratio = poisson_integral_right/observed_distribution_integral_right
                        fdr_cutoff = min(1, integral_ratio)
                        fdr_cutoffs[e_type][chunk][observed] = fdr_cutoff
                    else:
                        fdr_cutoffs[e_type][chunk][observed] = 0
        return fdr_cutoffs

    @staticmethod
    def _create_e_job(m, ij_pairs, ij_region_pairs, e=1, c=None,
                      chunks=None, w=1, p=0, min_locus_dist=3,
                      min_ll_reads=16, max_w=20, observed_cutoff=1,
                      e_ll_cutoff=None, e_h_cutoff=None,
                      e_v_cutoff=None, e_d_cutoff=None):

        if _has_gridmap:
            ij_pairs_compressed = msgpack.dumps(ij_pairs)
            ij_region_pairs_compressed = msgpack.dumps(ij_region_pairs)
            job = gridmap.Job(process_matrix_range, [m, ij_pairs_compressed, ij_region_pairs_compressed, e, c,
                              chunks, w, p, min_locus_dist, min_ll_reads,
                              max_w, observed_cutoff, e_ll_cutoff, e_h_cutoff,
                              e_v_cutoff, e_d_cutoff])
        else:
            # process with multiprocessing
            job = None
        return job

    def _find_peaks_in_matrix(self, m, e, c, inter, mappable, peak_info,
                              observed_chunk_distribution, lambda_chunks, w, p):

        jobs = []
        ij_pairs = []
        ij_region_pairs = []
        for i in xrange(m.shape[0]):
            i_region = m.row_regions[i].ix

            if not mappable[i_region]:
                continue

            for j in xrange(i, m.shape[1]):
                j_region = m.col_regions[j].ix

                if not mappable[j_region]:
                    continue

                ij_pairs.append((i, j))
                ij_region_pairs.append((i_region, j_region))

                if len(ij_pairs) > self.batch_size:
                    m_segment, updated_ij_pairs = RaoPeakCaller._submatrix_indices(m, ij_pairs, w_max=self.max_w)

                    job = RaoPeakCaller._create_e_job(m_segment, updated_ij_pairs,
                                                      ij_region_pairs[:], e, c, lambda_chunks,
                                                      w, p, self.min_locus_dist, self.min_ll_reads,
                                                      self.max_w, observed_cutoff=self.observed_cutoff,
                                                      e_ll_cutoff=self.e_ll_cutoff, e_h_cutoff=self.e_h_cutoff,
                                                      e_v_cutoff=self.e_v_cutoff, e_d_cutoff=self.e_d_cutoff)
                    jobs.append(job)
                    if len(jobs) >= self.n_processes:
                        self._process_jobs(jobs, peak_info, observed_chunk_distribution)
                        jobs = []
                    ij_pairs = []
                    ij_region_pairs = []

        if len(ij_pairs) > 0:
            # one last flush
            m_segment, updated_ij_pairs = RaoPeakCaller._submatrix_indices(m, ij_pairs, w_max=self.max_w)

            job = RaoPeakCaller._create_e_job(m_segment, updated_ij_pairs,
                                              ij_region_pairs[:], e, c, lambda_chunks,
                                              w, p, self.min_locus_dist, self.min_ll_reads,
                                              self.max_w, observed_cutoff=self.observed_cutoff,
                                              e_ll_cutoff=self.e_ll_cutoff, e_h_cutoff=self.e_h_cutoff,
                                              e_v_cutoff=self.e_v_cutoff, e_d_cutoff=self.e_d_cutoff)
            jobs.append(job)

        if len(jobs) > 0:
            self._process_jobs(jobs, peak_info, observed_chunk_distribution)
        peak_info.flush()

    def call_peaks(self, hic, chromosomes=None, file_name=None):
        peaks = RaoPeakInfo(file_name)

        peak_info = peaks.peak_table_intra
        peak_info_inter = peaks.peak_table_inter

        # mappability
        logging.info("Calculating visibility of regions...")
        mappable = hic.marginals() > 0
        logging.info("Done.")

        logging.info("Calculating expected values...")
        intra_expected, inter_expected = RaoPeakCaller.get_expected(hic, smooth=True)
        logging.info("Done.")

        # initialize peak parameters
        p = self.p
        w_init = self.w_init

        # empirically choose p and w_init
        # if not already chosen
        if p is None or w_init is None:
            bin_size = hic.bin_size()
            if bin_size > 25000:
                p = 1 if p is None else p
                w_init = 3 if w_init is None else w_init
            else:
                p = int(24999/bin_size) if p is None else p
                w_init = int(round(25000/bin_size) + 2) if w_init is None else w_init
        logging.info("Initial parameter values: p=%d, w=%d" % (p, w_init))

        logging.info("Obtaining bias vector...")
        c = hic.bias_vector()
        logging.info("Done.")

        logging.info("Finding maximum observed value...")
        max_observed = 0
        for edge in hic.edges(lazy=True):
            max_observed = max(max_observed, edge['weight']/(c[edge['source']]*c[edge['sink']]))
        logging.info("Done.")

        logging.info("Calculating lambda-chunk boundaries...")
        lambda_chunks = RaoPeakCaller._lambda_chunks(max_observed)
        logging.info("Done.")

        observed_chunk_distribution = RaoPeakCaller._get_chunk_distribution_container(lambda_chunks)

        # start processing chromosome pairs
        if chromosomes is None:
            chromosomes = list(hic.chromosomes())

        for i_chr, chromosome1 in enumerate(chromosomes):
            for j_chr in xrange(i_chr, len(chromosomes)):
                chromosome2 = chromosomes[j_chr]

                logging.info("Processing %s-%s" % (chromosome1, chromosome2))

                m = hic[chromosome1, chromosome2]
                if chromosome1 == chromosome2:
                    self._find_peaks_in_matrix(m, intra_expected, c, False, mappable, peak_info,
                                               observed_chunk_distribution, lambda_chunks, w_init, p)
                elif self.process_inter:
                    self._find_peaks_in_matrix(m, inter_expected, c, True, mappable, peak_info_inter,
                                               observed_chunk_distribution, lambda_chunks, w_init, p)
        peak_info.flush()

        # calculate fdrs
        fdr_cutoffs = RaoPeakCaller._get_fdr_cutoffs(lambda_chunks, observed_chunk_distribution)

        # inter_pvalues = []
        for peak in peaks.intra_peaks(lazy=True, auto_update=False):
            peak.fdr_ll = fdr_cutoffs['ll'][peak.e_ll_chunk][peak.observed]
            peak.fdr_h = fdr_cutoffs['h'][peak.e_h_chunk][peak.observed]
            peak.fdr_v = fdr_cutoffs['v'][peak.e_v_chunk][peak.observed]
            peak.fdr_d = fdr_cutoffs['d'][peak.e_d_chunk][peak.observed]
            peak.update()
        peak_info.flush()

        # return peak_info, fdr_cutoffs, observed_chunk_distribution
        return peaks

    def merge_peaks(self, peak_list, hic, euclidian_distance=20000):
        merged_list = []
        peak_list_ixs = []
        current_peak = None




def process_matrix_range(m, ij_pairs, ij_region_pairs, e, c, chunks, w=1, p=0,
                         min_locus_dist=3, min_ll_reads=16, max_w=20,
                         observed_cutoff=0, e_ll_cutoff=None, e_h_cutoff=None,
                         e_v_cutoff=None, e_d_cutoff=None):

    t1 = time.time()
    ij_pairs = msgpack.loads(ij_pairs)
    ij_region_pairs = msgpack.loads(ij_region_pairs)

    observed_chunk_distribution = RaoPeakCaller._get_chunk_distribution_container(chunks)

    e_all = partial(RaoPeakCaller.e_all, e=e, w=w, p=p)

    c_row_sub = c[m.row_regions[0].ix:m.row_regions[-1].ix+1]
    c_col_sub = c[m.col_regions[0].ix:m.col_regions[-1].ix+1]
    m_corr = np.zeros(m.shape)

    for i in xrange(m.shape[0]):
        for j in xrange(m.shape[1]):
            m_corr[i, j] = m[i, j]/c_row_sub[i]/c_col_sub[j]

    observed_list = []
    e_ll_list = []
    e_h_list = []
    e_v_list = []
    e_d_list = []
    region_list = []
    for ij_pair, ij_region_pair in zip(ij_pairs, ij_region_pairs):
        i = ij_pair[0]
        j = ij_pair[1]
        i_region = ij_region_pair[0]
        j_region = ij_region_pair[1]
        cf = c[i_region]*c[j_region]
        observed = m[i, j]
        observed_c = int(round(observed/cf))

        # do not examine loci closer than p+3
        if abs(i_region-j_region) <= p + min_locus_dist:
            e_ll, e_h, e_v, e_d = None, None, None, None
        else:
            # check the lower-left condition
            # assure minimum number of reads
            w_corr = w
            while RaoPeakCaller.ll_sum(m_corr, i, j, w=w_corr, p=p) < min_ll_reads and w_corr <= max_w:
                w_corr += 1

            if w_corr > max_w:
                e_ll, e_h, e_v, e_d = None, None, None, None
            else:
                e_ll, e_h, e_v, e_d = e_all(m, i, j, w=w_corr)
        e_ll_c = None if e_ll is None else e_ll/cf
        e_h_c = None if e_h is None else e_h/cf
        e_v_c = None if e_v is None else e_v/cf
        e_d_c = None if e_d is None else e_d/cf

        e_ll_chunk = RaoPeakCaller._find_chunk(chunks, e_ll_c)
        e_h_chunk = RaoPeakCaller._find_chunk(chunks, e_h_c)
        e_v_chunk = RaoPeakCaller._find_chunk(chunks, e_v_c)
        e_d_chunk = RaoPeakCaller._find_chunk(chunks, e_d_c)

        # update observed distribution
        if e_ll_chunk is not None:
            observed_chunk_distribution['ll'][e_ll_chunk][observed_c] += 1
        if e_h_chunk is not None:
            observed_chunk_distribution['h'][e_h_chunk][observed_c] += 1
        if e_v_chunk is not None:
            observed_chunk_distribution['v'][e_v_chunk][observed_c] += 1
        if e_d_chunk is not None:
            observed_chunk_distribution['d'][e_d_chunk][observed_c] += 1

        # do filtering here to reduce message size
        # check that we have enough observed reads
        if not observed_c > observed_cutoff:
            continue
        # check that we have a high-enough o/e ratio for ll
        if e_ll_cutoff is not None and e_ll_c > 0 and observed_c/e_ll_c < e_ll_cutoff:
            continue
        # check that we have a high-enough o/e ratio for ll
        if e_h_cutoff is not None and e_h_c > 0 and observed_c/e_h_c < e_h_cutoff:
            continue
        # check that we have a high-enough o/e ratio for ll
        if e_v_cutoff is not None and e_v_c > 0 and observed_c/e_v_c < e_v_cutoff:
            continue
        # check that we have a high-enough o/e ratio for ll
        if e_d_cutoff is not None and e_d_c > 0 and observed_c/e_d_c < e_d_cutoff:
            continue

        observed_list.append((observed_c, RaoPeakCaller._find_chunk(chunks, observed_c)))
        e_ll_list.append((e_ll_c, e_ll_chunk))
        e_h_list.append((e_h_c, e_h_chunk))
        e_v_list.append((e_v_c, e_v_chunk))
        e_d_list.append((e_d_c, e_d_chunk))
        region_list.append((i_region, j_region))

    t2 = time.time()

    # speed up information-passing between processes
    rv = msgpack.dumps((region_list, observed_list, e_ll_list,
                        e_h_list, e_v_list, e_d_list, observed_chunk_distribution))
    return rv
