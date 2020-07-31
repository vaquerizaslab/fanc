from __future__ import division
import os
from .config import config
from abc import abstractmethod, ABCMeta
import numpy as np
from scipy.stats import poisson
from collections import defaultdict, OrderedDict
import tables as t
from .matrix import RegionMatrixTable, Edge, LazyEdge
from .general import MaskFilter
import msgpack
import msgpack_numpy
import math
import pandas as pd
from .tools.general import RareUpdateProgressBar, pairwise
import warnings
from future.utils import with_metaclass, viewitems
from itertools import tee
try:
    from itertools import izip as zip
except ImportError:
    pass
import logging
logger = logging.getLogger(__name__)


msgpack_numpy.patch()


try:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        os.environ["CREATE_PLOTS"] = "0"  # prevent gridmap from overriding the matplotlib backend
        gridmap_logger = logging.getLogger('gridmap.conf')
        gridmap_logger.disabled = True  # shut gridmap up about not having drmaa available
        import gridmap
        gridmap_logger.disabled = False
    has_gridmap = True

    # prepare environment for better success rate
    if 'MAX_TIME_BETWEEN_HEARTBEATS' not in os.environ:
        os.environ['MAX_TIME_BETWEEN_HEARTBEATS'] = '600'
    if 'SEND_ERROR_MAIL' not in os.environ:
        os.environ['SEND_ERROR_MAIL'] = '0'
    if 'MAX_IDLE_HEARTBEATS' not in os.environ:
        os.environ['MAX_IDLE_HEARTBEATS'] = '20'
    if 'IDLE_THRESHOLD' not in os.environ:
        os.environ['IDLE_THRESHOLD'] = '0.5'
except ImportError:
    has_gridmap = False


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

    def __init__(self, file_name=None, mode='a', tmpdir=None,
                 _table_name_regions='regions', _table_name_peaks='edges'):
        """
        Initialize a PeakInfo object.

        :param file_name: If None, will create a working file in memory. If path to an
                          existing peak info file, will load its information. If path
                          to a non-existent file will create the file.
        :param mode: File mode, use 'a' for append, 'r' for read, and 'w' for write
        :param _table_name_regions: Internal, controls name of the region PyTables table
        :param _table_name_peaks: Internal, controls name of the peak PyTables table
        """
        additional_fields = {
            'weight': t.Float32Col(pos=2),
            'uncorrected': t.Int32Col(pos=3),
            'expected_local': t.Float32Col(pos=4),
            'oe': t.Float32Col(pos=5),
            'p_value': t.Float32Col(pos=6),
            'q_value_sum': t.Float32Col(pos=7),
            'x': t.Float32Col(pos=8),
            'y': t.Float32Col(pos=9),
            'radius': t.Float32Col(pos=10),
            'merged_pixels': t.Int32Col(pos=11),
        }

        RegionMatrixTable.__init__(self, file_name, mode=mode, tmpdir=tmpdir,
                                   additional_edge_fields=additional_fields,
                                   _table_name_regions=_table_name_regions,
                                   _table_name_edges=_table_name_peaks)

        self.peak_table = self._edges

    def _row_to_edge(self, row, lazy_edge=False, distances_in_bp=False):
        if distances_in_bp:
            bin_size = self.bin_size
        else:
            bin_size = 1

        if lazy_edge is None:
            source = row["source"]
            sink = row["sink"]
            d = dict()
            for field in self.field_names:
                if field != 'source' and field != 'sink':
                    value = row[field]
                    value = value.decode() if isinstance(value, bytes) else value
                    if field in ('x', 'y', 'radius'):
                        d[field] = value * bin_size
                    else:
                        d[field] = value

            source_node_row = self._regions[source]
            source_node = self.regions[source]
            sink_node_row = self._regions[sink]
            sink_node = self.regions[sink]
            return Peak(source_node, sink_node, **d)

        lazy_edge._row = row
        return lazy_edge

    def peaks(self, distances_in_bp=False, lazy=False):
        return self.edges(lazy=lazy, distances_in_bp=distances_in_bp)

    def _edges_subset(self, key=None, row_regions=None, col_regions=None,
                      lazy=False, lazy_edge=None, *args, **kwargs):
        if lazy:
            lazy_edge = LazyPeak(None, self._regions)
        else:
            lazy_edge = None
        return RegionMatrixTable._edges_subset(self, key=key, row_regions=row_regions,
                                               col_regions=col_regions, lazy=lazy,
                                               lazy_edge=lazy_edge, *args, **kwargs)

    def _edges_iter(self, lazy=False, lazy_edge=None, *args, **kwargs):
        if lazy:
            lazy_edge = LazyPeak(None, self._regions)
        else:
            lazy_edge = None
        return RegionMatrixTable._edges_iter(self, lazy=lazy,
                                             lazy_edge=lazy_edge,
                                             *args, **kwargs)

    def __iter__(self):
        return self.peaks()

    def filter_rao(self, queue=False):
        """
        Convenience function that applies a :class:`~RaoMergedPeakFilter`.

        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        """
        mask = self.add_mask_description('rao', 'Mask singlet peaks with a q-value sum > .02')
        rao_filter = RaoMergedPeakFilter(mask=mask)
        self.filter(rao_filter, queue)

    def to_bedpe(self, file_name, anchor_radius=True, score_field='q_value_sum', name_field=None):
        regions = list(self.regions)
        with open(file_name, 'w') as f:
            for peak in self.peaks(lazy=True, distances_in_bp=True):
                r1 = regions[peak.source]
                r2 = regions[peak.sink]
                if score_field is not None:
                    score = getattr(peak, score_field)
                else:
                    score = '.'
                if name_field is not None:
                    name = getattr(peak, name_field)
                else:
                    name = '.'

                if anchor_radius:
                    r = peak.radius
                    start1 = r1.start - r
                    end1 = r1.end + r
                    start2 = r2.start - r
                    end2 = r2.end + r
                else:
                    start1 = r1.start
                    end1 = r1.end
                    start2 = r2.start
                    end2 = r2.end

                f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    r1.chromosome, int(start1), int(end1),
                    r2.chromosome, int(start2), int(end2),
                    name, score
                ))


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

    For more information about neighborhoods and peak information,
    see :class:`~RaoPeakCaller`.
    """

    _classid = 'RAOPEAKINFO'

    def __init__(self, file_name=None, mode='a', tmpdir=None,
                 _table_name_regions='regions', _table_name_peaks='edges'):
        """
        Initialize a RaoPeakInfo object.

        :param file_name: If None, will create a working file in memory. If path to an
                          existing peak info file, will load its information. If path
                          to a non-existent file will create the file.
        :param mode: File mode, use 'a' for append, 'r' for read, and 'w' for write
        :param _table_name_regions: Internal, controls name of the region PyTables table
        :param _table_name_peaks: Internal, controls name of the peak PyTables table
        """

        additional_fields = {
            'source': t.Int32Col(pos=0),
            'sink': t.Int32Col(pos=1),
            'weight': t.Float32Col(pos=2),
            'uncorrected': t.Int32Col(pos=3),
            'w': t.Int32Col(pos=4),
            'p': t.Int32Col(pos=5),
            'e_ll': t.Float32Col(pos=6, dflt=1.0),
            'e_h': t.Float32Col(pos=7, dflt=1.0),
            'e_v': t.Float32Col(pos=8, dflt=1.0),
            'e_d': t.Float32Col(pos=9, dflt=1.0),
            'oe_ll': t.Float32Col(pos=10, dflt=1.0),
            'oe_h': t.Float32Col(pos=11, dflt=1.0),
            'oe_v': t.Float32Col(pos=12, dflt=1.0),
            'oe_d': t.Float32Col(pos=13, dflt=1.0),
            'mappability_ll': t.Float32Col(pos=14, dflt=1.0),
            'mappability_h': t.Float32Col(pos=15, dflt=1.0),
            'mappability_v': t.Float32Col(pos=16, dflt=1.0),
            'mappability_d': t.Float32Col(pos=17, dflt=1.0),
            'fdr_ll': t.Float32Col(pos=18, dflt=1.0),
            'fdr_h': t.Float32Col(pos=19, dflt=1.0),
            'fdr_v': t.Float32Col(pos=20, dflt=1.0),
            'fdr_d': t.Float32Col(pos=21, dflt=1.0),
            'e_ll_chunk': t.Int32Col(pos=22, dflt=1.0),
            'e_h_chunk': t.Int32Col(pos=23, dflt=1.0),
            'e_v_chunk': t.Int32Col(pos=24, dflt=1.0),
            'e_d_chunk': t.Int32Col(pos=25, dflt=1.0),
            'll_sum': t.Int32Col(pos=26, dflt=0),
        }

        RegionMatrixTable.__init__(self, file_name, mode=mode, tmpdir=tmpdir,
                                   additional_edge_fields=additional_fields,
                                   _table_name_regions=_table_name_regions,
                                   _table_name_edges=_table_name_peaks)

        self.peak_table = self._edges

    def peaks(self, lazy=False, auto_update=True, **kwargs):
        return self.edges(lazy=lazy, auto_update=auto_update, **kwargs)

    def peaks_sorted(self, sortby, lazy=False, auto_update=True, **kwargs):
        return self.edges_sorted(sortby, lazy=lazy, auto_update=auto_update, **kwargs)

    def __iter__(self):
        return self.peaks()

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

    def filter_observed(self, cutoff, queue=False):
        """
        Convenience function that applies a :class:`~ObservedPeakFilter`.
        The actual algorithm and rationale used for filtering will depend on the
        internal _mapper attribute.

        :param cutoff: Minimum observed value
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        """
        mask = self.add_mask_description('fdr', 'Mask peaks with an observed value lower than %e' % cutoff)
        observed_filter = ObservedPeakFilter(cutoff=cutoff, mask=mask)
        self.filter(observed_filter, queue)

    def filter_mappability(self, cutoff, queue=False):
        """
        Convenience function that applies a :class:`~MappabilityPeakFilter`.
        The actual algorithm and rationale used for filtering will depend on the
        internal _mapper attribute.

        :param cutoff: Minimum mappability (fraction of 1)
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        """
        mask = self.add_mask_description('mappability',
                                         'Mask peaks with a mappability lower than %e' % cutoff)
        mappability_filter = MappabilityPeakFilter(mappability_cutoff=cutoff, mask=mask)
        self.filter(mappability_filter, queue)

    def filter_enrichment(self, enrichment_ll_cutoff=1.0, enrichment_h_cutoff=1.0,
                          enrichment_v_cutoff=1.0, enrichment_d_cutoff=1.0, queue=False):
        """
        Convenience function that applies a :class:`~ObservedExpectedRatioPeakFilter`.
        The actual algorithm and rationale used for filtering will depend on the
        internal _mapper attribute.

        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        """
        mask = self.add_mask_description('o/e', 'Mask peaks with a low observed/expected ratio')
        enrichment_filter = EnrichmentPeakFilter(enrichment_ll_cutoff=enrichment_ll_cutoff,
                                                 enrichment_h_cutoff=enrichment_h_cutoff,
                                                 enrichment_v_cutoff=enrichment_v_cutoff,
                                                 enrichment_d_cutoff=enrichment_d_cutoff,
                                                 mask=mask)
        self.filter(enrichment_filter, queue)

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
        merged_peaks = PeakInfo(file_name=file_name, mode='w')
        merged_peaks.add_regions(self.regions(lazy=True), preserve_attributes=False)

        bin_size = self.bin_size

        def _append_merged_peak(peak_list):
            # add merged peak
            hp = peak_list[0]
            q_value_sum = hp.fdr_ll + hp.fdr_d + hp.fdr_h + hp.fdr_v
            oe = 1 if hp.e_d == 0 else hp.weight/hp.e_d
            merged_peak = Peak(source=hp.source, sink=hp.sink, weight=hp.weight,
                               uncorrected=hp.uncorrected, expected_local=hp.e_d,
                               p_value=hp.fdr_d, q_value_sum=q_value_sum, x=x, y=y,
                               radius=radius, oe=oe)
            merged_peaks.add_edge(merged_peak)

        merged_peak_counter = 0
        chromosome_names = self.chromosomes()
        for i, chromosome_name1 in enumerate(chromosome_names):
            for j in range(i, len(chromosome_names)):
                chromosome_name2 = chromosome_names[j]

                remaining_peaks_set = set()
                for peak in self.edge_subset((chromosome_name1, chromosome_name2), lazy=False):
                    remaining_peaks_set.add(peak)

                if len(remaining_peaks_set) == 0:
                    continue

                logger.info("Merging peaks in %s/%s" % (chromosome_name1, chromosome_name2))

                last_peak_number = 0
                current_peaks = []
                n_remaining = len(remaining_peaks_set)
                with RareUpdateProgressBar(max_value=n_remaining, silent=config.hide_progressbars,
                                           prefix="Peak merge") as pb:
                    while len(remaining_peaks_set) > 0:
                        x, y, radius = RaoPeakInfo._centroid_and_radius(current_peaks)

                        if len(current_peaks) == last_peak_number:
                            if len(current_peaks) > 0:
                                _append_merged_peak(current_peaks)
                                current_peaks = []
                                last_peak_number = 0
                                merged_peak_counter += 1

                            # find highest peak
                            highest_peak = None
                            for peak in remaining_peaks_set:
                                if highest_peak is None:
                                    highest_peak = peak
                                else:
                                    if highest_peak.weight < peak.weight:
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

                        pb.update(n_remaining-len(remaining_peaks_set))

                if len(current_peaks) > 0:
                    _append_merged_peak(current_peaks)
                    merged_peak_counter += 1

        logger.info("Total merged peaks: {}".format(merged_peak_counter))
        merged_peaks.flush()
        return merged_peaks

    def parameter_tool(self, *region_pairs, **kwargs):
        from .plotting.hic_plotter import PeakParameterPlot
        import matplotlib.pyplot as plt
        parameters_plot = PeakParameterPlot(self, **kwargs)
        parameters_plot.plot(*region_pairs)
        plt.show()

        return parameters_plot.observed_cutoff, \
               parameters_plot.oe_cutoffs, \
               parameters_plot.fdr_cutoffs, \
               parameters_plot.mappability_cutoffs


class Peak(Edge):
    """
    Container for a Peak/enriched contact in a Hi-C matrix.
    """
    def __init__(self, source, sink, *args, **kwargs):
        object.__setattr__(self, 'e_ll', None)
        object.__setattr__(self, 'e_h', None)
        object.__setattr__(self, 'e_d', None)
        object.__setattr__(self, 'e_v', None)
        super(Peak, self).__init__(source, sink, *args, **kwargs)


class LazyPeak(LazyEdge):
    """
    Container for a Peak/enriched contact in a Hi-C matrix.

    This class implements :class:`~LazyPeak`, which provides lazy
    loading of attributes from a PyTables table row.
    """
    def __init__(self, row, nodes_table, bin_size=1):
        LazyEdge.__init__(self, row, nodes_table)
        self._bin_size = bin_size

    def __getattr__(self, item):
        res = LazyEdge.__getattr__(self, item)
        if item in ('x', 'y', 'radius'):
            return self._bin_size * res
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


class ObservedPeakFilter(PeakFilter):
    """
    Filter for peaks that do not pass a certain FDR cutoff.
    """
    def __init__(self, cutoff=1, mask=None):
        """
        Initialize filter object.

        :param mask: A :class:`~fanc.data.general.Mask` object.
        :param cutoff: Minimum observed value to consider peak
        """
        super(ObservedPeakFilter, self).__init__(mask=mask)
        self.cutoff = cutoff

    def valid_peak(self, peak):
        """
        Evaluate whether a peak passes FDR cutoffs set in __init__
        :param peak: An :class:`~fanc.data.genomic.Edge` object
        :return: True if peak passes interal FDR cutoffs, False otherwise
        """
        if peak.uncorrected < self.cutoff:
            return False
        return True


class DistancePeakFilter(PeakFilter):
    """
    Filter for peaks where regions are closer than this cutoff in bins.
    """
    def __init__(self, cutoff=1, mask=None):
        """
        Initialize filter object.

        :param mask: A :class:`~fanc.data.general.Mask` object.
        :param cutoff: Minimum observed value to consider peak
        """
        super(DistancePeakFilter, self).__init__(mask=mask)
        self.cutoff = cutoff

    def valid_peak(self, peak):
        """
        Evaluate whether a peak passes FDR cutoffs set in __init__
        :param peak: An :class:`~fanc.data.genomic.Edge` object
        :return: True if peak passes internal FDR cutoffs, False otherwise
        """
        if abs(peak.source - peak.sink) < self.cutoff:
            return False
        return True


class FdrPeakFilter(PeakFilter):
    """
    Filter for peaks that do not pass a certain FDR cutoff.
    """
    def __init__(self, mask=None, fdr_cutoff=None,
                 fdr_ll_cutoff=None, fdr_v_cutoff=None,
                 fdr_h_cutoff=None, fdr_d_cutoff=None):
        """
        Initialize filter object.

        :param mask: A :class:`~fanc.data.general.Mask` object.
        :param fdr_cutoff: Global FDR cutoff. Is overridden by the
                           neighborhood-specific cutoffs
        :param fdr_ll_cutoff: FDR cutoff for the lower-left neighborhood
        :param fdr_v_cutoff: FDR cutoff for the vertical neighborhood
        :param fdr_h_cutoff: FDR cutoff for the horizontal neighborhood
        :param fdr_d_cutoff: FDR cutoff for the donut neighborhood
        """
        super(FdrPeakFilter, self).__init__(mask=mask)
        self.fdr_ll_cutoff = fdr_cutoff
        self.fdr_h_cutoff = fdr_cutoff
        self.fdr_v_cutoff = fdr_cutoff
        self.fdr_d_cutoff = fdr_cutoff
        if fdr_ll_cutoff is not None:
            self.fdr_ll_cutoff = fdr_ll_cutoff
        if fdr_h_cutoff is not None:
            self.fdr_h_cutoff = fdr_h_cutoff
        if fdr_v_cutoff is not None:
            self.fdr_v_cutoff = fdr_v_cutoff
        if fdr_d_cutoff is not None:
            self.fdr_d_cutoff = fdr_d_cutoff

    def valid_peak(self, peak):
        """
        Evaluate whether a peak passes FDR cutoffs set in __init__
        :param peak: An :class:`~fanc.data.genomic.Edge` object
        :return: True if peak passes interal FDR cutoffs, False otherwise
        """
        if self.fdr_ll_cutoff is not None and peak.fdr_ll > self.fdr_ll_cutoff:
            return False
        if self.fdr_h_cutoff is not None and peak.fdr_h > self.fdr_h_cutoff:
            return False
        if self.fdr_v_cutoff is not None and peak.fdr_v > self.fdr_v_cutoff:
            return False
        if self.fdr_d_cutoff is not None and peak.fdr_d > self.fdr_d_cutoff:
            return False
        return True


class MappabilityPeakFilter(PeakFilter):
    """
    Filter for peaks that do not pass a certain FDR cutoff.
    """
    def __init__(self, mask=None, mappability_cutoff=None,
                 mappability_ll_cutoff=None, mappability_v_cutoff=None,
                 mappability_h_cutoff=None, mappability_d_cutoff=None):
        """
        Initialize filter object.

        :param mask: A :class:`~fanc.data.general.Mask` object.
        :param mappability_cutoff: Global mappability cutoff. Is overridden by the
                                   neighborhood-specific cutoffs
        :param mappability_ll_cutoff: Mappability cutoff for the lower-left neighborhood
        :param mappability_v_cutoff: Mappability cutoff for the vertical neighborhood
        :param mappability_h_cutoff: Mappability cutoff for the horizontal neighborhood
        :param mappability_d_cutoff: Mappability cutoff for the donut neighborhood
        """
        super(MappabilityPeakFilter, self).__init__(mask=mask)
        self.mappability_ll_cutoff = mappability_cutoff
        self.mappability_h_cutoff = mappability_cutoff
        self.mappability_v_cutoff = mappability_cutoff
        self.mappability_d_cutoff = mappability_cutoff
        if mappability_ll_cutoff is not None:
            self.mappability_ll_cutoff = mappability_ll_cutoff
        if mappability_h_cutoff is not None:
            self.mappability_h_cutoff = mappability_h_cutoff
        if mappability_v_cutoff is not None:
            self.mappability_v_cutoff = mappability_v_cutoff
        if mappability_d_cutoff is not None:
            self.mappability_d_cutoff = mappability_d_cutoff

    def valid_peak(self, peak):
        """
        Evaluate whether a peak passes FDR cutoffs set in __init__
        :param peak: An :class:`~fanc.data.genomic.Edge` object
        :return: True if peak passes interal FDR cutoffs, False otherwise
        """
        if self.mappability_ll_cutoff is not None and peak.mappability_ll < self.mappability_ll_cutoff:
            return False
        if self.mappability_h_cutoff is not None and peak.mappability_h < self.mappability_h_cutoff:
            return False
        if self.mappability_v_cutoff is not None and peak.mappability_v < self.mappability_v_cutoff:
            return False
        if self.mappability_d_cutoff is not None and peak.mappability_d < self.mappability_d_cutoff:
            return False
        return True


class EnrichmentPeakFilter(PeakFilter):
    """
    Filter peaks that do not have a sufficiently strong observed/expected ratio.
    """
    def __init__(self, enrichment_cutoff=None, 
                 enrichment_ll_cutoff=None, enrichment_h_cutoff=None,
                 enrichment_v_cutoff=None, enrichment_d_cutoff=None, 
                 mask=None):
        """
        Initialize filter object.

        :param ll_ratio: Minimum observed/e_ll ratio
        :param h_ratio: Minimum observed/e_h ratio
        :param v_ratio: Minimum observed/e_v ratio
        :param d_ratio: Minimum observed/e_d ratio
        :param mask: A :class:`~fanc.data.general.Mask` object.
        :return: True if all observed/expected ratios pass the thresholds,
                 False otherwise
        """
        super(EnrichmentPeakFilter, self).__init__(mask=mask)
        self.enrichment_ll_cutoff = enrichment_cutoff
        self.enrichment_h_cutoff = enrichment_cutoff
        self.enrichment_v_cutoff = enrichment_cutoff
        self.enrichment_d_cutoff = enrichment_cutoff
        if enrichment_ll_cutoff is not None:
            self.enrichment_ll_cutoff = enrichment_ll_cutoff
        if enrichment_h_cutoff is not None:
            self.enrichment_h_cutoff = enrichment_h_cutoff
        if enrichment_v_cutoff is not None:
            self.enrichment_v_cutoff = enrichment_v_cutoff
        if enrichment_d_cutoff is not None:
            self.enrichment_d_cutoff = enrichment_d_cutoff

    def valid_peak(self, peak):
        # covering my ass for legacy_old bug
        if peak.e_d == 0 or peak.e_ll == 0 or peak.e_h == 0 or peak.e_v == 0:
            return False

        if self.enrichment_ll_cutoff is not None and peak.oe_ll < self.enrichment_ll_cutoff:
            return False
        if self.enrichment_h_cutoff is not None and peak.oe_h < self.enrichment_h_cutoff:
            return False
        if self.enrichment_v_cutoff is not None and peak.oe_v < self.enrichment_v_cutoff:
            return False
        if self.enrichment_d_cutoff is not None and peak.oe_d < self.enrichment_d_cutoff:
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
        # covering my ass for legacy_old bug
        if peak.e_d == 0 or peak.e_ll == 0 or peak.e_h == 0 or peak.e_v == 0:
            return False

        # 1.
        if peak.oe_d <= peak.oe_ll < 2.0:
            return False

        # 2.
        if peak.oe_h < 1.5 and peak.oe_v < 1.5:
            return False

        # 3.
        if peak.oe_d < 1.75 or peak.oe_ll < 1.75:
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


class FdrSumFilter(PeakFilter):
    """
    Remove peaks that have a q-value sum > cutoff.
    """
    def __init__(self, cutoff=1.0, mask=None):
        PeakFilter.__init__(self, mask=mask)
        self.cutoff = cutoff

    def valid_peak(self, peak):
        if peak.q_value_sum > self.cutoff:
            return False

        return True


class RaoPeakCaller(object):
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

    def __init__(self, p=None, w_init=None, min_locus_dist=None, max_w=20, min_ll_reads=16,
                 process_inter=False, correct_inter='fdr', n_processes=4,
                 slice_size=2000, min_mappable_fraction=0.7, cluster=False):
        """
        Initialize RaoPeakCaller with peak calling parameters.

        :param p: (int) "padding" of pixels around invesitgated peak
        :param w_init: initial width of the area around a pixel to investigate
        :param min_locus_dist: Minimal distance in bins between two loci to consider peak
        :param max_w: Maximal width after extending investigated area around peak
        :param min_ll_reads: Threshold for the number of reads in the lower-left
                             neighborhood of a pixel to consider it as a peak
        :param process_inter: If False, ignores inter-chromosomal peaks
        :param correct_inter: If None, does not correct inter-chromosomal peaks for
                              multiple testing. Other options are 'fdr' for Benjamini-
                              Hochberg correction of 'bonferroni' for Bonferroni correction
        :param n_processes: Number of processes to use for peak calling.
        :param slice_size: length of the matrix square investigated by each process.
        :param cluster: If True, attempts to call peaks using an SGE cluster. If False,
                        will use multiprocessing.
        """
        self.p = p
        self.w_init = w_init
        if p is not None and w_init is not None:
            if not w_init > p:
                raise ValueError("w_init ({}) must be strictly greater than p ({})!".format(w_init, p))

        self.min_locus_dist = min_locus_dist
        self.max_w = max_w
        self.min_ll_reads = min_ll_reads
        self.process_inter = process_inter
        self.correct_inter = correct_inter
        self.n_processes = n_processes
        self.slice_size = slice_size
        self.min_mappable_fraction = min_mappable_fraction
        self.cluster = cluster
        if self.cluster is True:
            if not has_gridmap:
                logger.warning("Cannot use the cluster because of previous error.")
                self.cluster = False

        super(RaoPeakCaller, self).__init__()

    # sum of reads in lower-left neighborhood
    @staticmethod
    def ll_sum(m, i, j, w=1, p=0):
        """
        Compute the sum of pixels in the lower-left neighborhood of a pixel.
        """
        i_max, j_max = m.shape

        sum1 = np.ma.sum(m[max(0, i+1):min(i_max, i+w+1), max(0, j-w):min(j_max, j)])
        sum2 = np.ma.sum(m[max(0, i+1):min(i_max, i+p+1), max(0, j-p):min(j_max, j)])

        return sum1 - sum2

    @staticmethod
    def enrichment(m, e, i, j, neighborhood_method, w=1, p=0):
        dividend = neighborhood_method(m, i, j, w, p)
        divisor = neighborhood_method(e, i, j, w, p)
        return float(dividend / divisor * e[i, j])

    @staticmethod
    def e_ll_sum(m, i, j, w=1, p=0):
        i_max, j_max = m.shape
        sum1 = np.ma.sum(m[max(0, i + 1):min(i + w + 1, i_max), max(0, j - w):min(j_max, j)])
        sum2 = np.ma.sum(m[max(0, i + 1):min(i_max, i + p + 1), max(0, j - p):min(j_max, j)])
        return sum1 - sum2

    @staticmethod
    def e_h_sum(m, i, j, w=1, p=0):
        i_max, j_max = m.shape
        sum1 = np.ma.sum(m[max(0, i - 1):min(i_max, i + 2), max(0, j - w):min(j_max, j + w + 1)])
        sum2 = np.ma.sum(m[max(0, i - 1):min(i_max, i + 2), max(0, j - p):min(j_max, j + p + 1)])
        return sum1 - sum2

    @staticmethod
    def e_v_sum(m, i, j, w=1, p=0):
        i_max, j_max = m.shape
        sum1 = np.ma.sum(m[max(0, i - w):min(i_max, i + w + 1), max(0, j - 1):min(j_max, j + 2)])
        sum2 = np.ma.sum(m[max(0, i - p):min(i_max, i + p + 1), max(0, j - 1):min(j_max, j + 2)])
        return sum1 - sum2

    @staticmethod
    def e_d_sum(m, i, j, w=1, p=0):
        i_max, j_max = m.shape
        sum1 = np.ma.sum(m[max(0, i - w):min(i_max, i + w + 1), max(0, j - w):min(j_max, j + w + 1)])
        sum2 = np.ma.sum(m[max(0, i - p):min(i_max, i + p + 1), max(0, j - p):min(j_max, j + p + 1)])
        sum3 = np.ma.sum(m[max(0, i - w):max(0, i - p), j])
        sum4 = np.ma.sum(m[max(0, i + p + 1):min(i_max, i + w + 1), j])
        sum5 = np.ma.sum(m[i, max(0, j - w):max(0, j - p)])
        sum6 = np.ma.sum(m[i, max(0, j + p + 1):min(j_max, j + w + 1)])

        sum1 = 0 if np.ma.is_masked(sum1) else sum1
        sum2 = 0 if np.ma.is_masked(sum2) else sum2
        sum3 = 0 if np.ma.is_masked(sum3) else sum3
        sum4 = 0 if np.ma.is_masked(sum4) else sum4
        sum5 = 0 if np.ma.is_masked(sum5) else sum5
        sum6 = 0 if np.ma.is_masked(sum6) else sum6

        return sum1 - sum2 - sum3 - sum4 - sum5 - sum6

    # lower-left neighborhood
    @staticmethod
    def e_ll(m, i, j, e, w=1, p=0):
        """
        Compute the average value of pixels in the lower-left neighborhood of a pixel.
        """
        return RaoPeakCaller.enrichment(m, e, i, j, RaoPeakCaller.e_ll_sum, w=w, p=p)

    # horizontal neighborhood
    @staticmethod
    def e_h(m, i, j, e, w=1, p=0):
        """
        Compute the average value of pixels in the horizontal neighborhood of a pixel.
        """
        return RaoPeakCaller.enrichment(m, e, i, j, RaoPeakCaller.e_h_sum, w=w, p=p)

    # vertical neighborhood
    @staticmethod
    def e_v(m, i, j, e, w=1, p=0):
        """
        Compute the average value of pixels in the vertical neighborhood of a pixel.
        """
        return RaoPeakCaller.enrichment(m, e, i, j, RaoPeakCaller.e_v_sum, w=w, p=p)

    # donut neighborhood
    @staticmethod
    def e_d(m, i, j, e, w=1, p=0):
        """
        Compute the average value of pixels in the "donut" neighborhood of a pixel.
        """
        return RaoPeakCaller.enrichment(m, e, i, j, RaoPeakCaller.e_d_sum, w=w, p=p)

    @staticmethod
    def find_chunk(value, chunk_func=lambda x: 3*np.log2(x)):
        """
        Use bisection to find a matching lambda chunk for a given expected value.
        """
        if value is None or np.isnan(value):
            return None
        if value < 1:
            return 0
        v = chunk_func(value)
        if not np.isfinite(v):
            return None
        return max(0, int(v) + 1)

    def _process_jobs(self, jobs, peaks, observed_chunk_distribution):
        """
        Process the output from :func:`~process_matrix_range` and save in peak table.
        """
        # if the grid does not work for some reason, this will fall back on
        # multiprocessing itself
        logger.debug("Getting gridmap output...")

        job_kwargs = {}
        if config.gridmap_tmpdir is not None:
            job_kwargs['temp_dir'] = config.gridmap_tmpdir
        job_outputs = gridmap.process_jobs(jobs, max_processes=self.n_processes,
                                           local=not self.cluster, **job_kwargs)
        logger.debug("Got gridmap output.")

        for compressed_results in job_outputs:
            results = msgpack.loads(compressed_results, strict_map_key=False)
            for result in results:
                found_none = False
                for value in result:
                    if value is None:
                        found_none = True
                if found_none:
                    continue

                source, sink, weight, w_corr, p, observed, \
                    ll_sum, e_ll, e_v, e_h, e_d, \
                    o_chunk, e_ll_chunk, e_v_chunk, e_h_chunk, e_d_chunk, \
                    e_ll_mappable, e_v_mappable, e_h_mappable, e_d_mappable = result

                # update observed distribution
                observed_chunk_distribution['ll'][e_ll_chunk][observed] += 1
                observed_chunk_distribution['h'][e_h_chunk][observed] += 1
                observed_chunk_distribution['v'][e_v_chunk][observed] += 1
                observed_chunk_distribution['d'][e_d_chunk][observed] += 1

                if weight == 0.0:
                    continue

                oe_ll = 1 if e_ll == 0 else weight / e_ll
                oe_h = 1 if e_h == 0 else weight / e_h
                oe_v = 1 if e_v == 0 else weight / e_v
                oe_d = 1 if e_d == 0 else weight / e_d

                peak = Edge(source=source, sink=sink,
                            weight=weight, uncorrected=observed,
                            w=w_corr, p=p, ll_sum=ll_sum,
                            e_ll=e_ll, e_h=e_h, e_v=e_v, e_d=e_d,
                            oe_ll=oe_ll, oe_h=oe_h, oe_v=oe_v, oe_d=oe_d,
                            e_ll_chunk=e_ll_chunk, e_v_chunk=e_v_chunk,
                            e_h_chunk=e_h_chunk, e_d_chunk=e_d_chunk,
                            mappability_ll=e_ll_mappable, mappability_v=e_v_mappable,
                            mappability_h=e_h_mappable, mappability_d=e_d_mappable)
                peaks.add_edge(peak)

    @staticmethod
    def _get_fdr_cutoffs(observed_chunk_distribution, e_func=lambda x: 2**(x/3)):
        """
        For all possible observed values in each lambda chunk, determine the
        FDR cutoff that denotes the lower significance bound.
        """
        # determine all possible observed values

        fdr_cutoffs = dict()
        for e_type in observed_chunk_distribution:  # ll, h, v, d
            fdr_cutoffs[e_type] = defaultdict(lambda: defaultdict(int))
            for chunk in observed_chunk_distribution[e_type].keys():
                max_e = e_func(chunk)
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

    @staticmethod
    def segment_matrix_intra(m, chunk_size, w_max):
        for i in range(0, m.shape[0], chunk_size):
            i_start = max(0, i - w_max)
            i_end = min(i + chunk_size + w_max, m.shape[0])
            i_range = (i_start, i_end)
            i_inspect = (i, min(i + chunk_size, m.shape[0]))
            for j in range(i, m.shape[1], chunk_size):
                j_start = max(0, j - w_max)
                j_end = min(j + chunk_size + w_max, m.shape[1])
                j_range = (j_start, j_end)
                j_inspect = (j, min(j + chunk_size, m.shape[1]))
                ms = m[i_start:i_end, j_start:j_end]
                yield ms, i_range, i_inspect, j_range, j_inspect

    def _find_peaks_intra_matrix(self, m, e, c, peak_info, mappable, ix_offset,
                                 observed_chunk_distribution, w, p):
        """
        Given a matrix (strictly intra-chromosomal), calculate peak
        information for all pixels.
        """

        jobs = []
        for segment in RaoPeakCaller.segment_matrix_intra(m, self.slice_size, self.max_w):
            ms, i_range, i_inspect, j_range, j_inspect = segment

            args = [ms, e, ix_offset,
                    i_range, i_inspect, mappable[i_range[0]:i_range[1]], c[i_range[0]:i_range[1]],
                    j_range, j_inspect, mappable[j_range[0]:j_range[1]], c[j_range[0]:j_range[1]],
                    w, p, self.min_locus_dist, self.min_ll_reads, self.min_mappable_fraction,
                    self.max_w]

            args = msgpack.dumps(args)
            job = gridmap.Job(process_matrix_segment_intra, [args])
            jobs.append(job)

            # submit intermediate segments if maximum number of jobs reached
            if len(jobs) >= self.n_processes:
                self._process_jobs(jobs, peak_info, observed_chunk_distribution)
                jobs = []

        if len(jobs) > 0:
            self._process_jobs(jobs, peak_info, observed_chunk_distribution)

    def call_peaks(self, hic, chromosome_pairs=None, file_name=None, intra_expected=None, inter_expected=None):
        """
        Call peaks in Hi-C matrix.

        This method will determine each pixel's likelihood to
        be a "true" peak. By default, only pixels with non-zero count and
        an observed/expected ratio > 1.0 for each neighborhood will be
        reported, because these can by definition not be true peaks.

        The peak calling behavior can be influenced by modifying
        the object attributes set when initializing :class:`~RaoPeakCaller`.

        :param hic: A :class:`~fanc.data.genomic.Hic` object
        :param chromosome_pairs: If None, all chromosome pairs will be
                                 investigated for peaks. Otherwise
                                 specify a list of chromosome name
                                 tuples (e.g. [('chr1', 'chr1'),
                                 ('chr1', 'chr3'), ...])
        :param file_name: An optional filename that backs the returned
                          :class:`~RaoPeakInfo` object
        :param intra_expected: A dict of the form
                               <chromosome>:<list of expected values> to override
                               expected value calculation
        :param inter_expected: A float describing the expected value
                               for inter-chromosomal contact matrix entries
        :return: :class:`~RaoPeakInfo` object
        """
        if self.process_inter:
            raise RuntimeError("Inter-chromosomal peak calling not currently supported!")

        peaks = RaoPeakInfo(file_name, mode='w')
        peaks.add_regions(hic.regions, preserve_attributes=False)

        # expected values
        if intra_expected is None:
            logger.info("Calculating intra-chromosomal expected values...")
            _, intra_expected, _ = hic.expected_values()
            logger.info("Done.")
        # if self.process_inter and inter_expected is None:
        #     logger.info("Inter-chromosomal expected values...")
        #     with ExpectedContacts(hic, smooth=False) as expected:
        #         inter_expected = expected.inter_expected()
        #
        # if self.process_inter:
        #     intra_possible, inter_possible = hic.possible_contacts()
        # else:
        #     intra_possible, inter_possible = None, None

        # mappability
        mappable = hic.mappable()
        
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
        if self.min_locus_dist is None:
            self.min_locus_dist = p
        logger.info("Initial parameter values: p=%d, w=%d" % (p, w_init))

        logger.info("Obtaining bias vector...")
        c = hic.bias_vector()
        logger.info("Done.")

        # lambda chunks container
        observed_chunk_distribution = dict()
        for e_type in ('ll', 'h', 'v', 'd'):
            observed_chunk_distribution[e_type] = defaultdict(lambda: defaultdict(int))

        # start processing chromosome pairs
        if chromosome_pairs is None:
            chromosome_pairs = []
            chromosomes = list(hic.chromosomes())

            for i_chr, chromosome1 in enumerate(chromosomes):
                for j_chr in range(i_chr, len(chromosomes)):
                    chromosome2 = chromosomes[j_chr]
                    chromosome_pairs.append((chromosome1, chromosome2))

        chromosome_bins = hic.chromosome_bins
        for chromosome1, chromosome2 in chromosome_pairs:
            logger.info("Processing %s-%s" % (chromosome1, chromosome2))

            start1, end1 = chromosome_bins[chromosome1]
            ix_offset = start1
            start2, end2 = chromosome_bins[chromosome2]
            if chromosome1 == chromosome2:
                m = hic.matrix((chromosome1, chromosome2))
                self._find_peaks_intra_matrix(m, intra_expected[chromosome1], c[start1:end1],
                                              peaks, mappable[start1:end1], ix_offset,
                                              observed_chunk_distribution, w_init, p)
            elif self.process_inter:
                warnings.warn("Inter-chromosomal peak calling not currently supported!")
                # self._find_peaks_inter_matrix(m, inter_expected, c[start1:end1], c[start2:end2],
                #                              peaks, mappable[start1:end1], mappable[start2:end2],
                #                              None, None, w_init, p)
        peaks.flush()

        # calculate fdrs
        logger.info("Finding FDR cutoffs...")
        fdr_cutoffs = RaoPeakCaller._get_fdr_cutoffs(observed_chunk_distribution)

        # regions_dict = peaks.regions_dict
        with RareUpdateProgressBar(max_value=len(peaks.edges), prefix='FDR') as pb:
            for i, peak in enumerate(peaks.peaks(lazy=True, writable=True)):
                # region1 = regions_dict[peak.source]
                # region2 = regions_dict[peak.sink]
                # if region1.chromosome == region2.chromosome:
                try:
                    peak.fdr_ll = fdr_cutoffs['ll'][peak.e_ll_chunk][peak.uncorrected]
                    peak.fdr_h = fdr_cutoffs['h'][peak.e_h_chunk][peak.uncorrected]
                    peak.fdr_v = fdr_cutoffs['v'][peak.e_v_chunk][peak.uncorrected]
                    peak.fdr_d = fdr_cutoffs['d'][peak.e_d_chunk][peak.uncorrected]
                except KeyError:
                    peak.fdr_ll = 1
                    peak.fdr_h = 1
                    peak.fdr_v = 1
                    peak.fdr_d = 1
                peak.update()

                pb.update(i)
                # else:
                #     # Bonferroni correction
                #     if self.correct_inter == 'bonferroni':
                #         peak.fdr_ll *= inter_possible
                #         peak.fdr_h *= inter_possible
                #         peak.fdr_v *= inter_possible
                #         peak.fdr_d *= inter_possible
                #         peak.update()
        peaks.flush()

        # if self.process_inter and self.correct_inter == 'fdr':
        #     # fdr_ll
        #     for i, peak in enumerate(peaks.peaks_sorted('fdr_ll', lazy=True, auto_update=False)):
        #         region1 = regions_dict[peak.source]
        #         region2 = regions_dict[peak.sink]
        #         if region1.chromosome != region2.chromosome:
        #             peak.fdr_ll *= inter_possible/(i+1)
        #             peak.update()
        #     peaks.flush()
        #     # fdr_h
        #     for i, peak in enumerate(peaks.peaks_sorted('fdr_h', lazy=True, auto_update=False)):
        #         region1 = regions_dict[peak.source]
        #         region2 = regions_dict[peak.sink]
        #         if region1.chromosome != region2.chromosome:
        #             peak.fdr_h *= inter_possible/(i+1)
        #             peak.update()
        #     peaks.flush()
        #     # fdr_v
        #     for i, peak in enumerate(peaks.peaks_sorted('fdr_v', lazy=True, auto_update=False)):
        #         region1 = regions_dict[peak.source]
        #         region2 = regions_dict[peak.sink]
        #         if region1.chromosome != region2.chromosome:
        #             peak.fdr_v *= inter_possible/(i+1)
        #             peak.update()
        #     peaks.flush()
        #     # fdr_d
        #     for i, peak in enumerate(peaks.peaks_sorted('fdr_d', lazy=True, auto_update=False)):
        #         region1 = regions_dict[peak.source]
        #         region2 = regions_dict[peak.sink]
        #         if region1.chromosome != region2.chromosome:
        #             peak.fdr_d *= inter_possible/(i+1)
        #             peak.update()
        #     peaks.flush()

        return peaks


def process_matrix_segment_intra(data):
    m_original, e, ix_offset, \
        i_range, i_inspect, mappable_i, c_i, \
        j_range, j_inspect, mappable_j, c_j, \
        w, p, min_locus_dist, min_ll_reads, min_mappable, \
        max_w = msgpack.loads(data, strict_map_key=False)

    m_original = np.array(m_original)
    # construct convenient matrices
    row_ixs = np.arange(i_range[0], i_range[1])
    col_ixs = np.arange(j_range[0], j_range[1])
    with np.errstate(divide='ignore', invalid='ignore'):
        m_uncorrected = np.rint(m_original/c_i[:, None]/c_j)
    m_distance = np.array([abs(col_ixs - i) for i in row_ixs])
    expected_f = np.vectorize(lambda x: e[x])
    m_expected = expected_f(m_distance)

    # mask above matrices by mappability
    mask = np.zeros(m_original.shape, dtype=bool)
    mask[np.logical_not(mappable_i)] = True
    mask[:, np.logical_not(mappable_j)] = True

    m_original = np.ma.masked_where(mask, m_original)
    m_uncorrected = np.ma.masked_where(mask, m_uncorrected)
    m_expected = np.ma.masked_where(mask, m_expected)
    m_ones = np.ones(m_original.shape)

    results = []
    for o_i in range(i_inspect[0], i_inspect[1]):
        i = o_i - i_range[0]
        for o_j in range(j_inspect[0], j_inspect[1]):
            j = o_j - j_range[0]

            # only inspect pixels at a certain distance
            # above the diagonal
            if o_j - o_i < p + min_locus_dist:
                continue

            # only inspect mappable pixels
            if mask[i, j]:
                continue

            # only inspect pixels if they have more than
            # a minimum number of reads
            w_corr = w
            ll_sum = 0
            while ll_sum < min_ll_reads and w_corr <= max_w:
                ll_sum = RaoPeakCaller.ll_sum(m_uncorrected, i, j, w=w_corr, p=p)
                if np.ma.is_masked(ll_sum):
                    ll_sum = 0
                w_corr += 1

            if w_corr > max_w:
                continue

            # calculate mappability and enrichment values
            ll_mappable = 1 - RaoPeakCaller.e_ll_sum(mask, i, j, w, p) / RaoPeakCaller.e_ll_sum(m_ones, i, j, w, p)
            if ll_mappable < min_mappable:
                continue
            e_ll = RaoPeakCaller.e_ll(m_original, i, j, m_expected, w=w_corr, p=p)

            v_mappable = 1 - RaoPeakCaller.e_v_sum(mask, i, j, w, p) / RaoPeakCaller.e_v_sum(m_ones, i, j, w, p)
            if v_mappable < min_mappable:
                continue
            e_v = RaoPeakCaller.e_v(m_original, i, j, m_expected, w=w_corr, p=p)

            h_mappable = 1 - RaoPeakCaller.e_h_sum(mask, i, j, w, p) / RaoPeakCaller.e_h_sum(m_ones, i, j, w, p)
            if h_mappable < min_mappable:
                continue
            e_h = RaoPeakCaller.e_h(m_original, i, j, m_expected, w=w_corr, p=p)

            d_mappable = 1 - RaoPeakCaller.e_d_sum(mask, i, j, w, p) / RaoPeakCaller.e_d_sum(m_ones, i, j, w, p)
            if d_mappable < min_mappable:
                continue
            e_d = RaoPeakCaller.e_d(m_original, i, j, m_expected, w=w_corr, p=p)

            # find chunks
            cf = c_i[i] * c_j[j]
            o_chunk = RaoPeakCaller.find_chunk(m_uncorrected[i, j])
            e_ll_chunk = RaoPeakCaller.find_chunk(e_ll/cf)
            e_h_chunk = RaoPeakCaller.find_chunk(e_h/cf)
            e_v_chunk = RaoPeakCaller.find_chunk(e_v/cf)
            e_d_chunk = RaoPeakCaller.find_chunk(e_d/cf)

            result = [o_i + ix_offset, o_j + ix_offset, float(m_original.data[i, j]), w_corr, p,
                      int(m_uncorrected.data[i, j]),
                      int(ll_sum), e_ll, e_v, e_h, e_d,
                      o_chunk, e_ll_chunk, e_v_chunk, e_h_chunk, e_d_chunk,
                      ll_mappable, v_mappable, h_mappable, d_mappable]

            results.append(result)
    return msgpack.dumps(results)


def overlap_peaks(peaks, max_distance=6000):
    """
    Calculate overlap between different peak calls.

    Useful for comparing peak calls across different samples
    or conditions.

    :param dict peaks: Peaks to overlap. Dictionary of
                       :class:`fanc.data.network.PeakInfo`,
                       keys are dataset names.
    :param int max_distance: Maximum distance between peaks for overlap
    :return: DataFrame of overlap statistics and dictionary
             containing overlapping peaks. Keys are sets
             of dataset names.
    :rtype: (pandas.DataFrame, fanc.data.network.PeakInfo)
    """
    # Algorithm from https://github.com/theaidenlab/juicebox/
    # blob/cb5999cb1e8e430dd29d4114fb208aca4b8d35ac/src/juicebox/
    # tools/utils/juicer/hiccups/HiCCUPSUtils.java#L235

    def key_func(p):
        try:
            return p[1].weight
        except AttributeError:
            return 1

    def hypotenuse(x, y):
        return math.sqrt(x*x + y*y)

    def mean(data):
        # https://stackoverflow.com/a/31062966
        n = 0
        mean = 0.0
        for x in data:
            n += 1
            mean += (x - mean)/n
        if n < 1:
            return float("nan")
        else:
            return mean

    summarize_attrs = [
        (mean, "weight"),
        (mean, "oe"),
        (mean, "uncorrected"),
        (mean, "expected_local"),
        (sum, "p_value"),
        (sum, "q_value_sum"),
    ]

    if not all(a == b for regions in zip(peaks.values()) for a, b in pairwise(regions)):
        raise ValueError("All peak calls must have the same regions.")

    peaks1 = next(iter(peaks.values()))
    bin_size = peaks1.bin_size
    max_distance = max_distance/bin_size
    logger.info("Fetching and sorting peaks...")

    all_peaks = list(sorted(((s, p) for s, pgen in viewitems(peaks) for p in pgen),
                            key=key_func, reverse=True))

    logger.info("Done.")
    logger.info("Finding overlaps...")
    out_peaks = defaultdict(list)
    total_n = len(all_peaks)

    with RareUpdateProgressBar(max_value=total_n, silent=config.hide_progressbars,
                               prefix="Overlap") as pb:
        while len(all_peaks) > 0:
            cur_p = all_peaks.pop(0)
            cur_p_list = [cur_p]
            cur_x = cur_p[1].x
            cur_y = cur_p[1].y
            cluster_radius = max_distance
            for p in all_peaks:
                if hypotenuse(cur_x - p[1].x, cur_y - p[1].y) <= cluster_radius:
                    cur_p_list.append(p)
                    cur_x = mean(_p.x for _s, _p in cur_p_list)
                    cur_y = mean(_p.y for _s, _p in cur_p_list)
                    r = max(hypotenuse(cur_x - _p.x, cur_y - _p.y) for _s, _p in cur_p_list)
                    cluster_radius = max_distance + r

            summed_attrs = {}
            for sum_func, attr in summarize_attrs:
                summed_attrs[attr] = sum_func(getattr(p, attr) for s, p in cur_p_list)

            cons_p = Peak(
                x=cur_x,
                y=cur_y,
                radius=r if len(cur_p_list) > 1 else 0.,
                source=math.floor(min(cur_x, cur_y)),
                sink=math.floor(max(cur_x, cur_y)),
                **summed_attrs
            )
            out_peaks[frozenset(s for s, p in cur_p_list)].append(cons_p)
            for p in cur_p_list[1:]:
                all_peaks.remove(p)
            pb.update(total_n - len(all_peaks))
    logger.info("Done.")
    logger.info("Gathering overlapped peaks.")
    out_dict = {}
    out_stats = []
    for sample_set, p_list in viewitems(out_peaks):
        pi = PeakInfo()
        pi.add_regions(peaks1.regions(), preserve_attributes=False)
        pi.add_edges(p_list)
        out_dict[sample_set] = pi
        stat = OrderedDict((s, s in sample_set) for s in peaks.keys())
        stat["n"] = len(p_list)
        out_stats.append(stat)
    logger.info("Done.")

    return pd.DataFrame(out_stats), out_dict
