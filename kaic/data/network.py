from __future__ import division
import os
from kaic.config import config
from abc import abstractmethod, ABCMeta
import numpy as np
from scipy.stats import poisson
from collections import defaultdict
import tables as t
from bisect import bisect_left
from functools import partial
from kaic.data.genomic import RegionMatrixTable, Edge, LazyEdge
from kaic.data.general import MaskFilter
import msgpack
import time
import multiprocessing
import math
from kaic.architecture.hic_architecture import ExpectedContacts
from kaic.tools.general import RareUpdateProgressBar
import warnings
from future.utils import with_metaclass
import logging
logger = logging.getLogger(__name__)

try:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        os.environ["CREATE_PLOTS"] = "0"  # prevent gridmap from overriding the matplotlib backend
        gridmap_logger = logging.getLogger('gridmap.conf')
        gridmap_logger.disabled = True  # shut gridmap up about not having drmaa available
        import gridmap
        gridmap_logger.disabled = False
    has_gridmap = True
except ImportError:
    has_gridmap = False


class PeakCaller(with_metaclass(ABCMeta, object)):

    def __init__(self):
        pass

    @abstractmethod
    def call_peaks(self, hic):
        pass


class PeakInfo(RegionMatrixTable):
    """
    General-purpose class for recording peaks in Hic (and similar) data.

    A peak has the following information:
    source, sink: coordinates of the highest peak pixel in the Hi-C matrix
    observed: observed value of the peak in the Hi-C matrix, generally uncorrected
    expected: expected value of the peak at this position in the Hi-C matrix
    p_value: a P-value that reflects how likely it is to observe a peak with
             these properties at random
    x, y: coordinates of the peak centroid, if it is larger than one pixel
    radius: radius of the peak, expressed in bins (can be converted to base pairs)
    """

    _classid = 'PEAKINFO'

    class MergedPeakInformation(t.IsDescription):
        source = t.Int32Col(pos=0)
        sink = t.Int32Col(pos=1)
        observed = t.Int32Col(pos=2)
        expected = t.Float32Col(pos=3)
        p_value = t.Float32Col(pos=4)
        q_value_sum = t.Float32Col(pos=5)
        x = t.Float32Col(pos=6)
        y = t.Float32Col(pos=7)
        radius = t.Float32Col(pos=8)

    def __init__(self, file_name=None, mode='a', tmpdir=None, regions=None, _table_name_regions='regions',
                 _table_name_peaks='edges'):
        """
        Initialize a PeakInfo object.

        :param file_name: If None, will create a working file in memory. If path to an
                          existing peak info file, will load its information. If path
                          to a non-existant file will create the file.
        :param mode: File mode, use 'a' for append, 'r' for read, and 'w' for write
        :param regions: Iterable with :class:`~GenomicRegion` objects to be loaded
        :param _table_name_regions: Internal, controls name of the region PyTables table
        :param _table_name_peaks: Internal, controls name of the peak PyTables table
        """

        RegionMatrixTable.__init__(self, file_name, mode=mode, tmpdir=tmpdir,
                                   additional_fields=PeakInfo.MergedPeakInformation,
                                   _table_name_nodes=_table_name_regions, _table_name_edges=_table_name_peaks)

        self.peak_table = self._edges

        if regions is not None:
            for region in regions:
                self.add_region(region, flush=False)
            self._regions.flush()

    def _row_to_edge(self, row, lazy=False, distances_in_bp=False, auto_update=True):
        if distances_in_bp:
            bin_size = self.bin_size
        else:
            bin_size = 1

        if not lazy:
            source = row["source"]
            sink = row["sink"]
            d = dict()
            for field in self.field_names:
                if field != 'source' and field != 'sink':
                    if field in ('x', 'y', 'radius'):
                        d[field] = row[field]*bin_size
                    else:
                        d[field] = row[field]

            source_node_row = self._regions[source]
            source_node = self._row_to_node(source_node_row)
            sink_node_row = self._regions[sink]
            sink_node = self._row_to_node(sink_node_row)
            return Peak(source_node, sink_node, **d)
        else:
            return LazyPeak(row, self._regions, bin_size=bin_size, auto_update=auto_update)

    def peaks(self, distances_in_bp=False, lazy=False, auto_update=True):
        return self.edges(lazy=lazy, distances_in_bp=distances_in_bp, auto_update=auto_update)

    def __iter__(self):
        return self.peaks()

    def filter(self, peak_filter, queue=False, log_progress=False):
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
            self.peak_table.filter(peak_filter, _logging=log_progress)
        else:
            self.peak_table.queue_filter(peak_filter)

    def filter_rao(self, queue=False):
        """
        Convenience function that applies a :class:`~RaoMergedPeakFilter`.

        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        """
        mask = self.add_mask_description('rao', 'Mask singlet peaks with a q-value sum < .02')
        rao_filter = RaoMergedPeakFilter(mask=mask)
        self.filter(rao_filter, queue)


class RaoPeakInfo(RegionMatrixTable):
    """
    Information about peaks called by :class:`~RaoPeakCaller`.

    A peak has the following information:

    source, sink: coordinates of the highest peak pixel in the Hi-C matrix
    observed: observed value of the peak in the Hi-C matrix, generally uncorrected
    e_ll: expected value of the peak given its lower-left neighborhood
    e_h: expected value of the peak given its horizontal neighborhood
    e_v: expected value of the peak given its vertical neighborhood
    e_d: expected value of the peak given its surrounding (donut) neighborhood
    e_ll_chunk: "lambda-chunk" this peak falls into given its 'll' neighborhood
    e_h_chunk: "lambda-chunk" this peak falls into given its 'h' neighborhood
    e_v_chunk: "lambda-chunk" this peak falls into given its 'v' neighborhood
    e_d_chunk: "lambda-chunk" this peak falls into given its 'd' neighborhood
    fdr_ll: FDR of the peak given its lower-left neighborhood
    fdr_h: FDR of the peak given its horizontal neighborhood
    fdr_v: FDR of the peak given its vertical neighborhood
    fdr_d: FDR of the peak given its surrounding (donut) neighborhood

    For more information about neighborhoods and peak infomration,
    see :class:`~RaoPeakCaller`.
    """

    _classid = 'RAOPEAKINFO'

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

    def __init__(self, file_name=None, mode='a', tmpdir=None, regions=None,
                 _table_name_regions='regions', _table_name_peaks='edges'):
        """
        Initialize a RaoPeakInfo object.

        :param file_name: If None, will create a working file in memory. If path to an
                          existing peak info file, will load its information. If path
                          to a non-existant file will create the file.
        :param mode: File mode, use 'a' for append, 'r' for read, and 'w' for write
        :param regions: Iterable with :class:`~GenomicRegion` objects to be loaded
        :param _table_name_regions: Internal, controls name of the region PyTables table
        :param _table_name_peaks: Internal, controls name of the peak PyTables table
        """

        RegionMatrixTable.__init__(self, file_name, mode=mode, tmpdir=tmpdir,
                                   additional_fields=RaoPeakInfo.PeakInformation,
                                   _table_name_nodes=_table_name_regions, _table_name_edges=_table_name_peaks)

        self.peak_table = self._edges

        if regions is not None:
            for region in regions:
                self.add_region(region, flush=False)
            self._regions.flush()

    def peaks(self, lazy=False, auto_update=True):
        return self.edges(lazy=lazy, auto_update=auto_update)

    def peaks_sorted(self, sortby, lazy=False, auto_update=True):
        return self.edges_sorted(sortby, lazy=lazy, auto_update=auto_update)

    def __iter__(self):
        return self.peaks()

    def filter(self, peak_filter, queue=False, log_progress=False):
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
            self.peak_table.filter(peak_filter, _logging=log_progress)
        else:
            self.peak_table.queue_filter(peak_filter)

    def filter_fdr(self, fdr_cutoff, queue=False):
        """
        Convenience function that applies a :class:`~FdrPeakFilter`.
        The actual algorithm and rationale used for filtering will depend on the
        internal _mapper attribute.

        :param fdr_cutoff: The false-discovery rate of every neighborhood enrichment
                           must be lower or equal to this threshold
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        """
        mask = self.add_mask_description('fdr', 'Mask peaks with an FDR higher than %e' % fdr_cutoff)
        fdr_filter = FdrPeakFilter(fdr_cutoff=fdr_cutoff, mask=mask)
        self.filter(fdr_filter, queue)

    def filter_observed_expected_ratio(self, ll_ratio=1.0, h_ratio=1.0, v_ratio=1.0, d_ratio=1.0, queue=False):
        """
        Convenience function that applies a :class:`~ObservedExpectedRatioPeakFilter`.
        The actual algorithm and rationale used for filtering will depend on the
        internal _mapper attribute.

        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        """
        mask = self.add_mask_description('o/e', 'Mask peaks with a low observed/expected ratio')
        oe_filter = ObservedExpectedRatioPeakFilter(ll_ratio=ll_ratio, h_ratio=h_ratio,
                                                    v_ratio=v_ratio, d_ratio=d_ratio,
                                                    mask=mask)
        self.filter(oe_filter, queue)

    def filter_rao(self, queue=False):
        """
        Convenience function that applies all filters Rao et al. (2014) do.

        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        """
        mask = self.add_mask_description('rao', 'Mask peaks that do not pass the RaoPeakFilter')
        rao_filter = RaoPeakFilter(mask=mask)
        self.filter(rao_filter, queue)

    @staticmethod
    def _euclidian_distance(x1, y1, x2, y2):
        """
        Determine the 2D euclidian distance between two points.
        """
        return math.sqrt((x1-x2)**2+(y1-y2)**2)

    @staticmethod
    def _centroid_and_radius(peak_list):
        """
        Determine the centroid coordinates and radius of list of peaks.
        """
        x = 0
        y = 0
        for peak in peak_list:
            x += peak.source
            y += peak.sink

        if len(peak_list) > 0:
            x /= len(peak_list)
            y /= len(peak_list)

        radius = 0
        for peak in peak_list:
            distance = RaoPeakInfo._euclidian_distance(x, y, peak.source, peak.sink)
            if distance > radius:
                radius = distance
        return x, y, radius

    def merged_peaks(self, file_name=None, euclidian_distance=20000):
        """
        Merge spatially proximal peaks.

        :param file_name: Optional file to save merged peak info to.
        :param euclidian_distance: Maximal distance in base pairs to still
                                   consider two peaks to be the same
        :return: :class:`~PeakInfo`
        """
        merged_peaks = PeakInfo(file_name=file_name, mode='w', regions=self.regions(lazy=True))

        # get region index
        regions_dict = self.regions_dict
        bin_size = self.bin_size

        def _append_merged_peak(peak_list):
            # add merged peak
            hp = peak_list[0]
            q_value_sum = hp.fdr_ll + hp.fdr_d + hp.fdr_h + hp.fdr_v
            merged_peak = Peak(source=hp.source, sink=hp.sink,
                               observed=hp.observed, expected=hp.e_d,
                               p_value=hp.fdr_d, q_value_sum=q_value_sum, x=x, y=y,
                               radius=radius)
            merged_peaks.add_edge(merged_peak, flush=False)

        chromosome_names = self.chromosomes()
        for i, chromosome_name1 in enumerate(chromosome_names):
            for j in range(i, len(chromosome_names)):
                chromosome_name2 = chromosome_names[j]

                logger.info("Merging peaks in %s/%s" % (chromosome_name1, chromosome_name2))
                remaining_peaks_set = set()
                for peak in self.peaks():
                    region1 = regions_dict[peak.source]
                    region2 = regions_dict[peak.sink]
                    if region1.chromosome == chromosome_name1 and region2.chromosome == chromosome_name2:
                        remaining_peaks_set.add(peak)

                last_peak_number = 0
                current_peaks = []
                l = len(remaining_peaks_set)
                with RareUpdateProgressBar(max_value=l, poll_interval=20, silent=config.hide_progressbars) as pb:
                    while len(remaining_peaks_set) > 0:
                        x, y, radius = RaoPeakInfo._centroid_and_radius(current_peaks)

                        if len(current_peaks) == last_peak_number:
                            if len(current_peaks) > 0:
                                _append_merged_peak(current_peaks)
                                current_peaks = []
                                last_peak_number = 0

                            # find highest peak
                            highest_peak = None
                            for peak in remaining_peaks_set:
                                if highest_peak is None:
                                    highest_peak = peak
                                else:
                                    if highest_peak.observed < peak.observed:
                                        highest_peak = peak

                            current_peaks.append(highest_peak)
                            remaining_peaks_set.remove(highest_peak)
                        else:
                            last_peak_number = len(current_peaks)

                            closest_peak = None
                            closest_distance = None
                            for peak in remaining_peaks_set:
                                distance = RaoPeakInfo._euclidian_distance(x, y, peak.source, peak.sink)
                                if closest_peak is None or distance < closest_distance:
                                    closest_peak = peak
                                    closest_distance = distance

                            if closest_distance*bin_size <= euclidian_distance+(radius*bin_size):
                                current_peaks.append(closest_peak)
                                remaining_peaks_set.remove(closest_peak)

                        pb.update(l-len(remaining_peaks_set))

                if len(current_peaks) > 0:
                    _append_merged_peak(current_peaks)
            merged_peaks.flush()
        return merged_peaks


class Peak(Edge):
    """
    Container for a Peak/enriched contact in a Hi-C matrix.
    """
    def __init__(self, source, sink, *args, **kwargs):
        super(Peak, self).__init__(source, sink, *args, **kwargs)


class LazyPeak(LazyEdge):
    """
    Container for a Peak/enriched contact in a Hi-C matrix.

    This class implements :class:`~LazyPeak`, which provides lazy
    loading of attributes from a PyTables table row.
    """
    def __init__(self, row, nodes_table, auto_update=True, bin_size=1):
        super(LazyPeak, self).__init__(row, nodes_table, auto_update=auto_update)
        self.reserved.add('bin_size')
        self.bin_size = bin_size

    def __getattr__(self, item):
        res = super(LazyPeak, self).__getattr__(item)
        if item in ('x', 'y', 'radius'):
            return self.bin_size*res
        return res


class PeakFilter(with_metaclass(ABCMeta, MaskFilter)):
    """
    Abstract class that provides filtering functionality for the
    peaks in a :class:`~RaoPeakInfo` object.

    Extends MaskFilter and overrides valid(self, row) to make
    :class:`~RaoPeakInfo` filtering more "natural".

    To create custom filters for the :class:`~RapPeakInfo` object, extend this
    class and override the valid_peak(self, peak) method.
    valid_peak should return False for a specific :class:`~Edge` object
    if the object is supposed to be filtered/masked and True
    otherwise. See :class:`~DiagonalFilter` for an example.

    Pass a custom filter to the :func:`~RaoPeakInfo.filter` method in :class:`~Hic`
    to apply it.
    """

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
        peak = LazyEdge(row)
        return self.valid_peak(peak)


class FdrPeakFilter(PeakFilter):
    """
    Filter for peaks that do not pass a certain FDR cutoff.
    """
    def __init__(self, mask=None, fdr_cutoff=None, fdr_ll_cutoff=0.1, fdr_v_cutoff=0.1,
                 fdr_h_cutoff=0.1, fdr_d_cutoff=0.1):
        """
        Initialize filter object.

        :param mask: A :class:`~kaic.data.general.Mask` object.
        :param fdr_cutoff: Global FDR cutoff. Is overridden by the
                           neighborhood-specific cutoffs
        :param fdr_ll_cutoff: FDR cutoff for the lower-left neighborhood
        :param fdr_v_cutoff: FDR cutoff for the vertical neighborhood
        :param fdr_h_cutoff: FDR cutoff for the horizontal neighborhood
        :param fdr_d_cutoff: FDR cutoff for the donut neighborhood
        """
        super(FdrPeakFilter, self).__init__(mask=mask)
        if fdr_cutoff is not None:
            fdr_ll_cutoff = fdr_cutoff
            fdr_h_cutoff = fdr_cutoff
            fdr_v_cutoff = fdr_cutoff
            fdr_d_cutoff = fdr_cutoff
        self.fdr_cutoff = fdr_cutoff
        self.fdr_ll_cutoff = fdr_ll_cutoff
        self.fdr_h_cutoff = fdr_h_cutoff
        self.fdr_v_cutoff = fdr_v_cutoff
        self.fdr_d_cutoff = fdr_d_cutoff

    def valid_peak(self, peak):
        """
        Evaluate whether a peak passes FDR cutoffs set in __init__
        :param peak: An :class:`~kaic.data.genomic.Edge` object
        :return: True if peak passes interal FDR cutoffs, False otherwise
        """
        if peak.fdr_ll > self.fdr_ll_cutoff:
            return False
        if peak.fdr_h > self.fdr_h_cutoff:
            return False
        if peak.fdr_v > self.fdr_v_cutoff:
            return False
        if peak.fdr_d > self.fdr_d_cutoff:
            return False
        return True


class ObservedExpectedRatioPeakFilter(PeakFilter):
    """
    Filter peaks that do not have a sufficiently strong observed/expected ratio.
    """
    def __init__(self, ll_ratio=1.0, h_ratio=1.0, v_ratio=1.0, d_ratio=1.0, mask=None):
        """
        Initialize filter object.

        :param ll_ratio: Minimum observed/e_ll ratio
        :param h_ratio: Minimum observed/e_h ratio
        :param v_ratio: Minimum observed/e_v ratio
        :param d_ratio: Minimum observed/e_d ratio
        :param mask: A :class:`~kaic.data.general.Mask` object.
        :return: True if all observed/expected ratios pass the thresholds,
                 False otherwise
        """
        super(ObservedExpectedRatioPeakFilter, self).__init__(mask=mask)
        self.ll_ratio = ll_ratio
        self.h_ratio = h_ratio
        self.v_ratio = v_ratio
        self.d_ratio = d_ratio

    def valid_peak(self, peak):
        # covering my ass for legacy bug
        if peak.e_d == 0 or peak.e_ll == 0 or peak.e_h == 0 or peak.e_v == 0:
            return False

        if self.ll_ratio is not None and peak.observed/peak.e_ll < self.ll_ratio:
            return False
        if self.h_ratio is not None and peak.observed/peak.e_h < self.h_ratio:
            return False
        if self.v_ratio is not None and peak.observed/peak.e_v < self.v_ratio:
            return False
        if self.d_ratio is not None and peak.observed/peak.e_d < self.d_ratio:
            return False
        return True


class RaoPeakFilter(PeakFilter):
    """
    Filter peaks exactly the same way that Rao et al. (2014) do.

    It only retains peaks that

    1. are at least 2-fold enriched over either the donut or lower-left neighborhood
    2. are at least 1.5-fold enriched over the horizontal and vertical neighborhoods
    3. are at least 1.75-fold enriched over both the donut and lower-left neighborhood
    4. have an FDR <= 0.1 in every neighborhood
    """
    def __init__(self, mask=None):
        PeakFilter.__init__(self, mask=mask)

    def valid_peak(self, peak):
        # covering my ass for legacy bug
        if peak.e_d == 0 or peak.e_ll == 0 or peak.e_h == 0 or peak.e_v == 0:
            return False

        # 1.
        if peak.observed/peak.e_d <= peak.observed/peak.e_ll < 2.0:
            return False

        # 2.
        if peak.observed/peak.e_h < 1.5 and peak.observed/peak.e_v < 1.5:
            return False

        # 3.
        if peak.observed/peak.e_d < 1.75 or peak.observed/peak.e_ll < 1.75:
            return False

        # 4.
        if peak.fdr_d > .1:
            return False
        if peak.fdr_ll > .1:
            return False
        if peak.fdr_h > .1:
            return False
        if peak.fdr_v > .1:
            return False

        return True


class RaoMergedPeakFilter(PeakFilter):
    """
    Filter merged peaks exactly the same way that Rao et al. (2014) do.

    It removes peaks that are singlets and have a q-value sum >.02.
    """
    def __init__(self, cutoff=.02, mask=None):
        PeakFilter.__init__(self, mask=mask)
        self.cutoff = cutoff

    def valid_peak(self, peak):
        if peak.radius == 0 and peak.q_value_sum > self.cutoff:
            return False

        return True


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

    FDRs for intra-chromosomal peaks are automatically corrected for multiple
    testing using the "lamda-chunking" methodology introduced in Rao et al. 2014.
    FDRs for inter-chromosomal peaks are corrected by default using the Benjamini
    Hochberg false-discovery rate correction (but 'bonferroni' is also an option)
    """

    def __init__(self, p=None, w_init=None, min_locus_dist=3, max_w=20, min_ll_reads=16,
                 observed_cutoff=1, e_ll_cutoff=1.0, e_h_cutoff=1.0, e_v_cutoff=1.0,
                 e_d_cutoff=1.0, process_inter=False, correct_inter='fdr', n_processes=4,
                 batch_size=500000, cluster=False):
        """
        Initialize RaoPeakCaller with peak calling parameters.

        :param p: (int) "padding" of pixels around invesitgated peak
        :param w_init: initial width of the area around a pixel to investigate
        :param min_locus_dist: Minimal distance in bins between two loci to consider peak
        :param max_w: Maximal width after extending investigated area around peak
        :param min_ll_reads: Threshold for the number of reads in the lower-left
                             neighborhood of a pixel to consider it as a peak
        :param observed_cutoff: Minimum (uncorrected) observed contact count for
                                a pixel to be reported
        :param e_ll_cutoff: Only report peaks with an observed/e_ll ratio >= e_ll_cutoff
        :param e_h_cutoff: Only report peaks with an observed/e_h ratio >= e_h_cutoff
        :param e_v_cutoff: Only report peaks with an observed/e_v ratio >= e_v_cutoff
        :param e_d_cutoff: Only report peaks with an observed/e_d ratio >= e_d_cutoff
        :param process_inter: If False, ignores inter-chromosomal peaks
        :param correct_inter: If None, does not correct inter-chromosomal peaks for
                              multiple testing. Other options are 'fdr' for Benjamini-
                              Hochberg correction of 'bonferroni' for Bonferroni correction
        :param n_processes: Number of processes to use for peak calling.
        :param batch_size: Number of pixels to investigate per batch.
        :param cluster: If True, attempts to call peaks using an SGE cluster. If False,
                        will use multiprocessing.
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
        self.correct_inter = correct_inter
        self.n_processes = n_processes
        self.batch_size = batch_size
        self.mpqueue = None
        self.cluster = cluster
        if self.cluster is True:
            if not has_gridmap:
                logger.warn("Cannot use the cluster because of previous error.")
                self.cluster = False

        super(RaoPeakCaller, self).__init__()

    @staticmethod
    def chromosome_map(hic):
        """
        Make a quick-lookup chromosome map.
        """
        chromosome_map = dict()
        for i, chromosome in enumerate(hic.chromosomes()):
            chromosome_map[chromosome] = i

        chromosomes = np.zeros(len(hic.regions), dtype=int)
        for i, region in enumerate(hic.regions(lazy=True)):
            chromosomes[i] = chromosome_map[region.chromosome]
        return chromosomes

    # sum of reads in lower-left neighborhood
    @staticmethod
    def ll_sum(m, i, j, w=1, p=0):
        """
        Compute the sum of pixels in the lower-left neighborhood of a pixel.
        """
        i_max, j_max = m.shape

        sum1 = np.sum(m[max(0, i+1):min(i_max, i+w+1), max(0, j-w):min(j_max, j)])
        sum2 = np.sum(m[max(0, i+1):min(i_max, i+p+1), max(0, j-p):min(j_max, j)])

        return sum1 - sum2

    # lower-left neighborhood
    @staticmethod
    def e_ll(m, i, j, e, w=1, p=0):
        """
        Compute the average value of pixels in the lower-left neighborhood of a pixel.
        """
        i_max, j_max = m.shape

        # dividend
        sum1 = np.sum(m[max(0, i+1):min(i+w+1, i_max), max(0, j-w):min(j_max, j)])
        sum2 = np.sum(m[max(0, i+1):min(i_max, i+p+1), max(0, j-p):min(j_max, j)])

        if sum1-sum2 == 0:
            return 0

        # divisor
        sum3 = 0
        for a in range(max(0, i+1), min(i_max, i+w+1)):
            for b in range(max(0, j-w), min(j_max, j)):
                sum3 += e(a, b)

        sum4 = 0
        for a in range(max(0, i+1), min(i_max, i+p+1)):
            for b in range(max(0, j-p), min(j_max, j)):
                sum4 += e(a, b)

        return (sum1-sum2)/(sum3-sum4)*e(i, j)

    # horizontal neighborhood
    @staticmethod
    def e_h(m, i, j, e, w=1, p=0):
        """
        Compute the average value of pixels in the horizontal neighborhood of a pixel.
        """
        i_max, j_max = m.shape

        # dividend
        sum1 = np.sum(m[max(0, i-1):min(i_max, i+2), max(0, j-w):min(j_max, j+w+1)])
        sum2 = np.sum(m[max(0, i-1):min(i_max, i+2), max(0, j-p):min(j_max, j+p+1)])

        if sum1-sum2 == 0:
            return 0

        # divisor
        sum3 = 0
        for a in range(max(0, i-1), min(i_max, i+2)):
            for b in range(max(0, j-w), min(j_max, j+w+1)):
                sum3 += e(a, b)

        sum4 = 0
        for a in range(max(0, i-1), min(i_max, i+2)):
            for b in range(max(0, j-p), min(j_max, j+p+1)):
                sum4 += e(a, b)

        return (sum1-sum2)/(sum3-sum4)*e(i, j)

    # vertical neighborhood
    @staticmethod
    def e_v(m, i, j, e, w=1, p=0):
        """
        Compute the average value of pixels in the vertical neighborhood of a pixel.
        """
        i_max, j_max = m.shape

        # dividend
        sum1 = np.sum(m[max(0, i-w):min(i_max, i+w+1), max(0, j-1):min(j_max, j+2)])
        sum2 = np.sum(m[max(0, i-p):min(i_max, i+p+1), max(0, j-1):min(j_max, j+2)])

        if sum1-sum2 == 0:
            return 0

        # divisor
        sum3 = 0
        for a in range(max(0, i-w), min(i_max, i+w+1)):
            for b in range(max(0, j-1), min(j_max, j+2)):
                sum3 += e(a, b)

        sum4 = 0
        for a in range(max(0, i-p), min(i_max, i+p+1)):
            for b in range(max(0, j-1), min(j_max, j+2)):
                sum4 += e(a, b)

        return (sum1-sum2)/(sum3-sum4)*e(i, j)

    # donut neighborhood
    @staticmethod
    def e_d(m, i, j, e, w=1, p=0):
        """
        Compute the average value of pixels in the "donut" neighborhood of a pixel.
        """
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
        for a in range(max(0, i-w), min(i_max, i+w+1)):
            for b in range(max(0, j-w), min(j_max, j+w+1)):
                bottom_sum1 += e(a, b)

        bottom_sum2 = 0
        for a in range(max(0, i-p), min(i_max, i+p+1)):
            for b in range(max(0, j-p), min(j_max, j+p+1)):
                bottom_sum2 += e(a, b)

        bottom_sum3 = 0
        for a in range(max(0, i-w), min(i_max, i-p)):
            bottom_sum3 += e(a, j)

        bottom_sum4 = 0
        for a in range(max(0, i+p+1), min(i_max, i+w+1)):
            bottom_sum4 += e(a, j)

        bottom_sum5 = 0
        for b in range(max(0, j-w), min(j_max, j-p)):
            bottom_sum5 += e(i, b)

        bottom_sum6 = 0
        for b in range(max(0, j+p+1), min(j_max, j+w+1)):
            bottom_sum6 += e(i, b)

        return top / \
            (bottom_sum1-bottom_sum2-bottom_sum3-bottom_sum4-bottom_sum5-bottom_sum6) * \
            e(i, j)

    @staticmethod
    def e_all(m, i, j, e, w=1, p=0):
        """
        Compute the average value of pixels in all neighborhoods of a pixel.
        """
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
        """
        Compute the expected value "lambda chunks" to classify pixels.
        """
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
        """
        Use bisection to find a matching lambda chunk for a given expected value.
        """
        if value is None:
            return None
        return bisect_left(chunk_list, value)

    @staticmethod
    def _submatrix_indices(m, ij_pairs, w_max=20):
        """
        Given a matrix and a list of index pairs, return a submatrix (and matching converted indices)
        that accommodates all pairs.
        """
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

    def _process_jobs(self, jobs, peaks, observed_chunk_distribution):
        """
        Process the output from :func:`~process_matrix_range` and save in peak table.
        """
        if self.cluster:
            # if the grid does not work for some reason, this will fall back on
            # multiprocessing itself
            job_outputs = gridmap.process_jobs(jobs, max_processes=self.n_processes)
        else:
            # use multiprocessing
            for p in jobs:
                p.start()

            # Exit the completed processes
            for p in jobs:
                p.join()

            # Get process results from the output queue
            job_outputs = [self.mpqueue.get() for _ in jobs]
            self.mpqueue = None

        for output in job_outputs:
            rv = output

            (region_pairs, observed_list, e_ll_list,
             e_h_list, e_v_list, e_d_list,
             observed_chunk_distribution_part) = msgpack.loads(rv)
            logger.info("Got output")

            # intra-chromosomal
            if observed_chunk_distribution_part is not None:
                for ix in range(len(region_pairs)):
                    source = region_pairs[ix][0]
                    sink = region_pairs[ix][1]

                    observed, observed_chunk = observed_list[ix]
                    e_ll, e_ll_chunk = e_ll_list[ix]
                    e_h, e_h_chunk = e_h_list[ix]
                    e_v, e_v_chunk = e_v_list[ix]
                    e_d, e_d_chunk = e_d_list[ix]

                    if (e_ll is not None and e_h is not None and
                            e_v is not None and e_d is not None):
                        peak = Edge(source=source, sink=sink, observed=observed, e_ll=e_ll, e_h=e_h, e_v=e_v, e_d=e_d,
                                    e_ll_chunk=e_ll_chunk, e_h_chunk=e_h_chunk, e_v_chunk=e_v_chunk, e_d_chunk=e_d_chunk)
                        peaks.add_edge(peak, flush=False)

                for e_type in observed_chunk_distribution_part.keys():
                    for chunk_ix in range(len(observed_chunk_distribution_part[e_type])):
                        for o in observed_chunk_distribution_part[e_type][chunk_ix].keys():
                            et = e_type.decode() if isinstance(e_type, bytes) else e_type
                            observed_chunk_distribution[et][chunk_ix][o] += observed_chunk_distribution_part[e_type][chunk_ix][o]
            else:
                for ix in range(len(region_pairs)):
                    source = region_pairs[ix][0]
                    sink = region_pairs[ix][1]

                    observed, _ = observed_list[ix]
                    e_ll, fdr_ll = e_ll_list[ix]
                    e_h, fdr_h = e_h_list[ix]
                    e_v, fdr_v = e_v_list[ix]
                    e_d, fdr_d = e_d_list[ix]

                    if (e_ll is not None and e_h is not None and
                            e_v is not None and e_d is not None):
                        peak = Edge(source=source, sink=sink, observed=observed, e_ll=e_ll, e_h=e_h, e_v=e_v, e_d=e_d,
                                    fdr_ll=fdr_ll, fdr_v=fdr_v, fdr_h=fdr_h, fdr_d=fdr_d)
                        peaks.add_edge(peak, flush=False)
        peaks.flush()

    @staticmethod
    def _get_chunk_distribution_container(lambda_chunks):
        observed_chunk_distribution = {
            'll': [],
            'h': [],
            'v': [],
            'd': []
        }
        for _ in range(len(lambda_chunks)):
            observed_chunk_distribution['ll'].append(defaultdict(int))
            observed_chunk_distribution['h'].append(defaultdict(int))
            observed_chunk_distribution['v'].append(defaultdict(int))
            observed_chunk_distribution['d'].append(defaultdict(int))
        return observed_chunk_distribution

    @staticmethod
    def _get_fdr_cutoffs(lambda_chunks, observed_chunk_distribution):
        """
        For all possible observed values in each lambda chunk, determine the
        FDR cutoff that denotes the lower significance bound.
        """
        # determine all possible observed values

        fdr_cutoffs = dict()
        for e_type in observed_chunk_distribution:  # ll, h, v, d
            fdr_cutoffs[e_type] = []
            for chunk, max_e in enumerate(lambda_chunks):  # (0, 1), (1, 1.26), (2, 1.59), ...
                fdr_cutoffs[e_type].append(dict())
                poisson_e = poisson(max_e)

                observed_sum = 0
                for observed_count in observed_chunk_distribution[e_type][chunk].values():
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

    def _create_e_job(self, m, ij_pairs, ij_region_pairs, e=1, c=None,
                      chunks=None, w=1, p=0, min_locus_dist=3,
                      min_ll_reads=16, max_w=20, observed_cutoff=1,
                      e_ll_cutoff=None, e_h_cutoff=None,
                      e_v_cutoff=None, e_d_cutoff=None):

        ij_pairs_compressed = msgpack.dumps(ij_pairs)
        ij_region_pairs_compressed = msgpack.dumps(ij_region_pairs)

        if self.cluster:
            job = gridmap.Job(process_matrix_range, [m, ij_pairs_compressed, ij_region_pairs_compressed, e, c,
                              chunks, w, p, min_locus_dist, min_ll_reads,
                              max_w, observed_cutoff, e_ll_cutoff, e_h_cutoff,
                              e_v_cutoff, e_d_cutoff])
        else:
            if self.mpqueue is None:
                self.mpqueue = multiprocessing.Queue(self.n_processes)
            # process with multiprocessing
            job = multiprocessing.Process(target=multiprocessing_matrix_range,
                                          args=(m, ij_pairs_compressed, ij_region_pairs_compressed, e, c,
                                                chunks, w, p, min_locus_dist, min_ll_reads,
                                                max_w, observed_cutoff, e_ll_cutoff, e_h_cutoff,
                                                e_v_cutoff, e_d_cutoff, self.mpqueue))
        return job

    def _find_peaks_in_matrix(self, m, e, c, mappable, peak_info,
                              observed_chunk_distribution, lambda_chunks, w, p):
        """
        Given a matrix (strictly inter- OR intra-chromosomal), calculate peak
        information for all pixels.
        """
        jobs = []
        ij_pairs = []
        ij_region_pairs = []
        for i in range(m.shape[0]):
            i_region = m.row_regions[i].ix

            if not mappable[i_region]:
                continue

            # make sure we investigate whole inter-chromosomal matrix
            if lambda_chunks is None:
                start_j = 0
            else:
                start_j = i

            for j in range(start_j, m.shape[1]):
                j_region = m.col_regions[j].ix

                if not mappable[j_region]:
                    continue

                # if this is inter-chromosomal data, and we don't
                # have a positive value, skip it
                if lambda_chunks is None and m[i, j] == 0:
                    continue

                ij_pairs.append((i, j))
                ij_region_pairs.append((i_region, j_region))

                if len(ij_pairs) > self.batch_size:
                    m_segment, updated_ij_pairs = RaoPeakCaller._submatrix_indices(m, ij_pairs, w_max=self.max_w)

                    job = self._create_e_job(m_segment, updated_ij_pairs,
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

            job = self._create_e_job(m_segment, updated_ij_pairs,
                                     ij_region_pairs[:], e, c, lambda_chunks,
                                     w, p, self.min_locus_dist, self.min_ll_reads,
                                     self.max_w, observed_cutoff=self.observed_cutoff,
                                     e_ll_cutoff=self.e_ll_cutoff, e_h_cutoff=self.e_h_cutoff,
                                     e_v_cutoff=self.e_v_cutoff, e_d_cutoff=self.e_d_cutoff)
            jobs.append(job)

        if len(jobs) > 0:
            self._process_jobs(jobs, peak_info, observed_chunk_distribution)
        peak_info.flush()

    def call_peaks(self, hic, chromosome_pairs=None, file_name=None, expected=None):
        """
        Call peaks in Hi-C matrix.

        This method will determine each pixel's likelihood to
        be a "true" peak. By default, only pixels with non-zero count and
        an observed/expected ratio > 1.0 for each neighborhood will be
        reported, because these can by defninition not be true peaks.

        The peak calling behavior can be influenced by modifying
        the object attributes set when initializing :class:`~RaoPeakCaller`.

        :param hic: A :class:`~kaic.data.genomic.Hic` object
        :param chromosome_pairs: If None, all chromosome pairs will be
                                 investigated for peaks. Otherwise
                                 specify a list of chromosome name
                                 tuples (e.g. [('chr1', 'chr1'),
                                 ('chr1', 'chr3'), ...])
        :param file_name: An optional filename that backs the returned
                          :class:`~RaoPeakInfo` object
        :return: :class:`~RaoPeakInfo` object
        """
        peaks = RaoPeakInfo(file_name, regions=hic.regions(lazy=True), mode='w')

        # mappability
        calculated_expected = False
        if expected is None:
            logger.info("Expected values...")
            calculated_expected = True
            expected = ExpectedContacts(hic, smooth=True)
        intra_expected = expected.intra_expected()
        inter_expected = expected.inter_expected()

        intra_possible, inter_possible = hic.possible_contacts()

        # mappability
        mappable = expected.marginals() > 0

        logger.info("Done.")
        
        # initialize peak parameters
        p = self.p
        w_init = self.w_init

        # empirically choose p and w_init
        # if not already chosen
        if p is None or w_init is None:
            bin_size = hic.bin_size
            if bin_size > 25000:
                p = 1 if p is None else p
                w_init = 3 if w_init is None else w_init
            else:
                p = int(24999/bin_size) if p is None else p
                w_init = int(25000/bin_size + 0.5) + 2 if w_init is None else w_init
        logger.info("Initial parameter values: p=%d, w=%d" % (p, w_init))

        logger.info("Obtaining bias vector...")
        c = hic.bias_vector()
        logger.info("Done.")

        logger.info("Finding maximum observed value...")
        max_observed = 0
        for edge in hic.edges(lazy=True):
            new_max = edge.weight/(c[edge.source]*c[edge.sink])
            if not math.isinf(new_max):
                max_observed = max(max_observed, new_max)
        logger.info("Done.")

        logger.info("Calculating lambda-chunk boundaries...")
        lambda_chunks = RaoPeakCaller._lambda_chunks(max_observed*2)
        logger.info("Done.")

        observed_chunk_distribution = RaoPeakCaller._get_chunk_distribution_container(lambda_chunks)
        inter_stats = {'total': 0, 'observed': 0}

        # start processing chromosome pairs
        if chromosome_pairs is None:
            chromosome_pairs = []
            chromosomes = list(hic.chromosomes())

            for i_chr, chromosome1 in enumerate(chromosomes):
                for j_chr in range(i_chr, len(chromosomes)):
                    chromosome2 = chromosomes[j_chr]
                    chromosome_pairs.append((chromosome1, chromosome2))

        for chromosome1, chromosome2 in chromosome_pairs:
            logger.info("Processing %s-%s" % (chromosome1, chromosome2))

            m = hic[chromosome1, chromosome2]
            if chromosome1 == chromosome2:
                self._find_peaks_in_matrix(m, intra_expected, c, mappable, peaks,
                                           observed_chunk_distribution, lambda_chunks, w_init, p)
            elif self.process_inter:
                self._find_peaks_in_matrix(m, inter_expected, c, mappable, peaks,
                                           None, None, w_init, p)
        peaks.flush()

        # calculate fdrs
        fdr_cutoffs = RaoPeakCaller._get_fdr_cutoffs(lambda_chunks, observed_chunk_distribution)

        regions_dict = peaks.regions_dict
        for peak in peaks.peaks(lazy=True, auto_update=False):
            region1 = regions_dict[peak.source]
            region2 = regions_dict[peak.sink]
            if region1.chromosome == region2.chromosome:
                try:
                    peak.fdr_ll = fdr_cutoffs['ll'][peak.e_ll_chunk][peak.observed]
                    peak.fdr_h = fdr_cutoffs['h'][peak.e_h_chunk][peak.observed]
                    peak.fdr_v = fdr_cutoffs['v'][peak.e_v_chunk][peak.observed]
                    peak.fdr_d = fdr_cutoffs['d'][peak.e_d_chunk][peak.observed]
                    peak.update()
                except KeyError:
                    peak.fdr_ll = 0
                    peak.fdr_h = 0
                    peak.fdr_v = 0
                    peak.fdr_d = 0
            else:
                # Bonferroni correction
                if self.correct_inter == 'bonferroni':
                    peak.fdr_ll *= inter_possible
                    peak.fdr_h *= inter_possible
                    peak.fdr_v *= inter_possible
                    peak.fdr_d *= inter_possible
                    peak.update()
        peaks.flush()

        if self.process_inter and self.correct_inter == 'fdr':
            # fdr_ll
            for i, peak in enumerate(peaks.peaks_sorted('fdr_ll', lazy=True, auto_update=False)):
                region1 = regions_dict[peak.source]
                region2 = regions_dict[peak.sink]
                if region1.chromosome != region2.chromosome:
                    peak.fdr_ll *= inter_possible/(i+1)
                    peak.update()
            peaks.flush()
            # fdr_h
            for i, peak in enumerate(peaks.peaks_sorted('fdr_h', lazy=True, auto_update=False)):
                region1 = regions_dict[peak.source]
                region2 = regions_dict[peak.sink]
                if region1.chromosome != region2.chromosome:
                    peak.fdr_h *= inter_possible/(i+1)
                    peak.update()
            peaks.flush()
            # fdr_v
            for i, peak in enumerate(peaks.peaks_sorted('fdr_v', lazy=True, auto_update=False)):
                region1 = regions_dict[peak.source]
                region2 = regions_dict[peak.sink]
                if region1.chromosome != region2.chromosome:
                    peak.fdr_v *= inter_possible/(i+1)
                    peak.update()
            peaks.flush()
            # fdr_d
            for i, peak in enumerate(peaks.peaks_sorted('fdr_d', lazy=True, auto_update=False)):
                region1 = regions_dict[peak.source]
                region2 = regions_dict[peak.sink]
                if region1.chromosome != region2.chromosome:
                    peak.fdr_d *= inter_possible/(i+1)
                    peak.update()
            peaks.flush()

        if calculated_expected:
            expected.close()

        return peaks


def multiprocessing_matrix_range(m, ij_pairs, ij_region_pairs, e, c, chunks, w=1, p=0,
                                 min_locus_dist=3, min_ll_reads=16, max_w=20,
                                 observed_cutoff=0, e_ll_cutoff=None, e_h_cutoff=None,
                                 e_v_cutoff=None, e_d_cutoff=None, queue=None):
    output = process_matrix_range(m, ij_pairs, ij_region_pairs, e, c, chunks, w=w, p=p,
                                  min_locus_dist=min_locus_dist, min_ll_reads=min_ll_reads, max_w=max_w,
                                  observed_cutoff=observed_cutoff, e_ll_cutoff=e_ll_cutoff, e_h_cutoff=e_h_cutoff,
                                  e_v_cutoff=e_v_cutoff, e_d_cutoff=e_d_cutoff)
    queue.put(output)


def process_matrix_range(m, ij_pairs, ij_region_pairs, e, c, chunks, w=1, p=0,
                         min_locus_dist=3, min_ll_reads=16, max_w=20,
                         observed_cutoff=0, e_ll_cutoff=None, e_h_cutoff=None,
                         e_v_cutoff=None, e_d_cutoff=None):

    t1 = time.time()
    ij_pairs = msgpack.loads(ij_pairs)
    ij_region_pairs = msgpack.loads(ij_region_pairs)

    if chunks is not None:
        observed_chunk_distribution = RaoPeakCaller._get_chunk_distribution_container(chunks)
    else:
        observed_chunk_distribution = None

    e_all = partial(RaoPeakCaller.e_all, e=e, w=w, p=p)

    c_row_sub = c[m.row_regions[0].ix:m.row_regions[-1].ix+1]
    c_col_sub = c[m.col_regions[0].ix:m.col_regions[-1].ix+1]
    m_corr = np.zeros(m.shape)

    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
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
        observed_c = int(observed/cf + 0.5)

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

        if e_ll is None or e_ll == 0:
            continue
        if e_d is None or e_d == 0:
            continue
        if e_h is None or e_h == 0:
            continue
        if e_v is None or e_v == 0:
            continue

        e_ll_c = e_ll/cf
        e_h_c = e_h/cf
        e_v_c = e_v/cf
        e_d_c = e_d/cf

        if observed_chunk_distribution is not None:
            e_ll_chunk = RaoPeakCaller._find_chunk(chunks, e_ll_c)
            e_h_chunk = RaoPeakCaller._find_chunk(chunks, e_h_c)
            e_v_chunk = RaoPeakCaller._find_chunk(chunks, e_v_c)
            e_d_chunk = RaoPeakCaller._find_chunk(chunks, e_d_c)

            # update observed distribution
            try:
                if e_ll_chunk is not None:
                    observed_chunk_distribution['ll'][e_ll_chunk][observed_c] += 1

                if e_h_chunk is not None:
                    observed_chunk_distribution['h'][e_h_chunk][observed_c] += 1

                if e_v_chunk is not None:
                    observed_chunk_distribution['v'][e_v_chunk][observed_c] += 1

                if e_d_chunk is not None:
                    observed_chunk_distribution['d'][e_d_chunk][observed_c] += 1

            except IndexError:
                logger.error("Chunk distribution index error")
                logger.error("observed_c: %d" % observed_c)
                logger.error("e_ll_chunk: %s" % str(e_ll_chunk))
                logger.error("e_h_chunk: %s" % str(e_h_chunk))
                logger.error("e_v_chunk: %s" % str(e_v_chunk))
                logger.error("e_d_chunk: %s" % str(e_d_chunk))
                #continue

        else:
            # calculate fdrs instead
            e_ll_chunk = 1-poisson(e_ll_c).cdf(observed_c)
            e_h_chunk = 1-poisson(e_h_c).cdf(observed_c)
            e_v_chunk = 1-poisson(e_v_c).cdf(observed_c)
            e_d_chunk = 1-poisson(e_d_c).cdf(observed_c)

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

        if observed_chunk_distribution is not None:
            observed_list.append((observed_c, RaoPeakCaller._find_chunk(chunks, observed_c)))
        else:
            observed_list.append((observed_c, None))
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
