"""
Module for working with genomic data.

This module provides classes and functions to work with objects in the context
of the genome.

:class:`~Chromosome`, :class:`~Genome`, and :class:`~GenomicRegion` simplify
working with reference sequence data by providing easy access and many convenience
functions.

Examples:

.. code:: python

    # assemble custom genome
    chr1 = Chromosome.from_fasta("/path/to/chr1_fasta_file")
    chr2 = Chromosome.from_fasta("/path/to/chr2_fasta_file")
    genome = Genome(chromosomes=[chr1,chr2])

    # extract genomic regions binned every 10000bp
    regions = genome.get_regions(10000)

.. code:: python

    # assemble genome from folder with FASTA files
    genome = Genome.from_folder("/path/to/fasta_folder/")

:class:`~Hic` is the central class for working with Hi-C data. It provides
matrix-like selectors and convenient access to specific genomic regions. In the
fanc pipeline, a Hic object is assembled at the fragment level from
:class:`~fanc.construct.seq.FragmentMappedReadPairs`. From there, it can be
binned to equi-distant genomic regions.

.. code:: python

    # use previously existing FragmentMappedReadPairs object 'pairs'
    hic = Hic(file_name="/path/to/save_file")
    hic.load_read_fragment_pairs(pairs)

    # bin Hi-C object
    binned = hic.bin(10000)

    # ... further processing


Alternatively, a Hic object can be assembled from scratch using genomic
regions and edges (contacts) between them.

Example:

.. code:: python

    hic = Hic()
    genome = Genome.from_folder("/path/to/fasta_folder")
    hic.add_regions(genome.get_regions(10000))

    hic.add_edges(list_of_edges)

"""

from __future__ import division, print_function
from fanc.config import config
import tables as t
import pandas as p
import numpy as np
import pysam
import pybedtools
from fanc.tools.files import is_fasta_file, write_bigwig, write_bed, write_gff
from fanc.tools.matrix import apply_sliding_func
from fanc.tools.general import natural_sort, natural_cmp
from Bio import SeqIO, Restriction, Seq
from fanc.legacy.data.general import TableObject, Maskable, MaskedTable, MaskFilter, FileGroup
from abc import abstractmethod, ABCMeta
import os.path
from fanc.tools.general import ranges, distribute_integer, create_col_index, \
    RareUpdateProgressBar, range_overlap, str_to_int
try:
    from itertools import izip as zip
except ImportError:
    pass
import pickle
from collections import defaultdict
import copy
import re
import shlex
import warnings
from bisect import bisect_right, bisect_left
from future.utils import with_metaclass, string_types, viewitems
from builtins import object
from pandas import DataFrame, read_table
import subprocess
import pyBigWig
import intervaltree
import logging
logger = logging.getLogger(__name__)


def _edge_overlap_split_rao(original_edge, overlap_map):
    """
    Resolve the distribution of contacts when binning using
    Rao et al. 2014 approach.
    """
    original_source = original_edge[0]
    original_sink = original_edge[1]
    original_weight = original_edge[2]

    new_source_nodes = overlap_map[original_source]
    new_sink_nodes = overlap_map[original_sink]
    
    if len(new_source_nodes) == 0:
        return []
    elif len(new_source_nodes) == 1:
        new_source_nodes = [new_source_nodes[0][0]]
    else:
        new_source_nodes = [new_source_nodes[0][0], new_source_nodes[-1][0]]
    
    if len(new_sink_nodes) == 0:
        return []
    elif len(new_sink_nodes) == 1:
        new_sink_nodes = [new_sink_nodes[0][0]]
    else:
        new_sink_nodes = [new_sink_nodes[0][0], new_sink_nodes[-1][0]]
    
    edges = {}
    for new_source in new_source_nodes:
        for new_sink in new_sink_nodes:
            if new_source <= new_sink:
                edges[(new_source, new_sink)] = 0
            else:
                edges[(new_sink, new_source)] = 0
    
    weights = distribute_integer(original_weight, len(edges))
    edges_list = []
    for i, key_pair in enumerate(edges):
        edges_list.append([key_pair[0], key_pair[1], weights[i]])

    return edges_list


def _weighted_mean(intervals):
    intervals = np.array(intervals)
    if len(intervals) == 0:
        return np.nan
    mask = np.isfinite(intervals[:, 2])
    valid = intervals[mask]
    if len(valid) == 0:
        return np.nan
    weights = (valid[:, 1] - valid[:, 0])
    weights += 1
    # safety
    weights = [weight if weight > 0 else 1 for weight in weights]
    return np.average(valid[:, 2], weights=weights)


def as_region(region):
    if isinstance(region, string_types):
        return GenomicRegion.from_string(region)
    elif isinstance(region, GenomicRegion):
        return region
    raise ValueError("region parameter cannot be converted to GenomicRegion!")


class RegionBased(object):
    """
    Base class for working with genomic regions.
    
    Guide for inheriting classes which functions to override:
    
    MUST (basic functionality):
        _region_iter
        _get_regions
    
    SHOULD (works if above are implemented, but is highly inefficient):
        _region_subset
        _region_intervals
    
    CAN (override for potential speed benefits or added functionality):
        _region_len
        chromosomes
        chromosome_lens
        region_bins
    """
    def __init__(self):
        pass

    @property
    def _estimate_region_bounds(self):
        return True

    @property
    def file_type(self):
        return 'region'

    def _region_iter(self, *args, **kwargs):
        raise NotImplementedError("Function not implemented")

    def _get_regions(self, item, *args, **kwargs):
        raise NotImplementedError("Function not implemented")

    def _region_subset(self, region, *args, **kwargs):
        return self.regions[self.region_bins(region)]

    def _region_intervals(self, region, *args, **kwargs):
        intervals = []
        for region in self.regions(region, *args, **kwargs):
            intervals.append((region.start, region.end, region.score))
        return intervals

    def _region_len(self):
        return sum(1 for _ in self.regions)

    def __len__(self):
        return self._region_len()

    def __getitem__(self, item):
        return self._get_regions(item)

    def __iter__(self):
        return self.regions()

    @property
    def regions(self):
        """
        Iterate over genomic regions in this object.

        Will return a :class:`~GenomicRegion` object in every iteration.
        Can also be used to get the number of regions by calling
        len() on the object returned by this method.

        :return: RegionIter
        """
        class RegionIter(object):
            def __init__(self, region_based):
                self._region_based = region_based

            def __len__(self):
                return self._region_based._region_len()

            def __iter__(self):
                return self()

            def _fix_chromosome(self, regions):
                for r in regions:
                    if r.chromosome.startswith('chr'):
                        r.chromosome = r.chromosome[3:]
                    else:
                        try:
                            r.chromosome = 'chr' + r.chromosome
                        except AttributeError:
                            r = r.copy()
                            r.chromosome = 'chr' + r.chromosome
                    yield r

            def __call__(self, key=None, *args, **kwargs):
                fix_chromosome = kwargs.pop('fix_chromosome', False)

                if key is None:
                    iterator = self._region_based._region_iter(*args, **kwargs)
                else:
                    iterator = self._region_based.subset(key, *args, **kwargs)

                if fix_chromosome:
                    return self._fix_chromosome(iterator)
                else:
                    return iterator

            def __getitem__(self, item):
                return self._region_based._get_regions(item)

        return RegionIter(self)

    def _convert_region(self, region):
        """
        Take any object that can be interpreted as a region and return a :class:`GenomicRegion`.
        
        :param region: Any object interpretable as genomic region (string, :class:`GenomicRegion`)
        :return: :class:`GenomicRegion`
        """
        if isinstance(region, string_types):
            region = GenomicRegion.from_string(region)

        if isinstance(region, GenomicRegion):
            if region.start is None and self._estimate_region_bounds:
                region.start = 0

            if region.end is None and self._estimate_region_bounds:
                chromosome_lengths = self.chromosome_lengths
                if region.chromosome in chromosome_lengths:
                    region.end = chromosome_lengths[region.chromosome]
        return region

    def chromosomes(self):
        """
        Get a list of chromosome names.
        """
        chromosomes_set = set()
        chromosomes = []
        for region in self.regions:
            if region.chromosome not in chromosomes_set:
                chromosomes_set.add(region.chromosome)
                chromosomes.append(region.chromosome)
        return chromosomes

    @property
    def chromosome_lens(self):
        return self.chromosome_lengths

    @property
    def chromosome_lengths(self):
        """
        Returns a dictionary of chromosomes and their length
        in bp.
        """
        chr_lens = {}
        for r in self.regions:
            if chr_lens.get(r.chromosome) is None:
                chr_lens[r.chromosome] = r.end
                continue
            if r.end > chr_lens[r.chromosome]:
                chr_lens[r.chromosome] = r.end
        return chr_lens

    def region_bins(self, region):
        """
        Takes a genomic region and returns a slice of the bin
        indices that are covered by the region.

        :param region: String or class:`~GenomicRegion`
                       object for which covered bins will
                       be returned.
        :return: slice
        """
        region = self._convert_region(region)

        start_ix = None
        end_ix = None
        for i, r in enumerate(self.regions):
            ix = r.ix if hasattr(r, 'ix') and r.ix is not None else i
            if not (r.chromosome == region.chromosome and r.start <= region.end and r.end >= region.start):
                continue
            if start_ix is None:
                start_ix = ix
                end_ix = ix + 1
                continue
            end_ix = ix + 1
        return slice(start_ix, end_ix)

    def find_region(self, query_regions, _regions_dict=None, _region_ends=None, _chromosomes=None):
        """
        Find the region that is at the center of a region.

        :param query_regions: Region selector string, :class:~GenomicRegion, or
                              list of the former
        :return: index (or list of indexes) of the region at the center of the
                 query region
        """
        is_single = False
        if isinstance(query_regions, string_types):
            is_single = True
            query_regions = [GenomicRegion.from_string(query_regions)]

        if isinstance(query_regions, GenomicRegion):
            is_single = True
            query_regions = [query_regions]

        if _regions_dict is None or _region_ends is None or _chromosomes is None:
            regions_dict = defaultdict(list)
            region_ends = defaultdict(list)
            chromosomes = set()

            for region in self.regions:
                regions_dict[region.chromosome].append(region)
                region_ends[region.chromosome].append(region.end)
                chromosomes.add(region.chromosome)
        else:
            regions_dict = _regions_dict
            region_ends = _region_ends
            chromosomes = _chromosomes

        hit_regions = []
        for query_region in query_regions:
            if isinstance(query_region, string_types):
                query_region = GenomicRegion.from_string(query_region)

            if query_region.chromosome not in chromosomes:
                hit_regions.append(None)
                continue

            center = query_region.start + (query_region.end-query_region.start)/2
            ix = bisect_left(region_ends[query_region.chromosome], center)
            try:
                hit_regions.append(regions_dict[query_region.chromosome][ix])
            except IndexError:
                hit_regions.append(None)
        if is_single:
            return hit_regions[0]
        return hit_regions

    def subset(self, region, *args, **kwargs):
        """
        Takes a class:`~GenomicRegion` and returns all regions that
        overlap with the supplied region.

        :param region: String or class:`~GenomicRegion`
                       object for which covered bins will
                       be returned.
        """
        region = self._convert_region(region)
        return self._region_subset(region, *args, **kwargs)

    def region_intervals(self, region, bins=None, bin_size=None, smoothing_window=None,
                         nan_replacement=None, zero_to_nan=False, *args, **kwargs):
        region = self._convert_region(region)
        if not isinstance(region, GenomicRegion):
            raise ValueError("Region must be a GenomicRegion object or equivalent string!")

        raw_intervals = self._region_intervals(region, *args, **kwargs)

        if raw_intervals is None:
            raw_intervals = []

        if bins is None and bin_size is None:
            return raw_intervals

        if bins is not None:
            return RegionBased.bin_intervals(raw_intervals, bins,
                                             interval_range=[region.start, region.end],
                                             smoothing_window=smoothing_window,
                                             nan_replacement=nan_replacement,
                                             zero_to_nan=zero_to_nan)

        if bin_size is not None:
            return RegionBased.bin_intervals_equidistant(raw_intervals, bin_size,
                                                         interval_range=[region.start, region.end],
                                                         smoothing_window=smoothing_window,
                                                         nan_replacement=nan_replacement,
                                                         zero_to_nan=zero_to_nan)

    def intervals(self, *args, **kwargs):
        return self.region_intervals(*args, **kwargs)

    def binned_regions(self, region=None, bins=None, bin_size=None, smoothing_window=None,
                       nan_replacement=None, zero_to_nan=False, *args, **kwargs):
        region = self._convert_region(region)
        br = []
        if region is None:
            for chromosome in self.chromosomes():
                interval_bins = self.region_intervals(chromosome, bins=bins, bin_size=bin_size,
                                                      smoothing_window=smoothing_window,
                                                      nan_replacement=nan_replacement,
                                                      zero_to_nan=zero_to_nan, *args, **kwargs)
                br += [GenomicRegion(chromosome=chromosome, start=interval_bin[0],
                                     end=interval_bin[1], score=interval_bin[2])
                       for interval_bin in interval_bins]
        else:
            interval_bins = self.region_intervals(region, bins=bins, bin_size=bin_size,
                                                  smoothing_window=smoothing_window,
                                                  nan_replacement=nan_replacement,
                                                  zero_to_nan=zero_to_nan, *args, **kwargs)
            br += [GenomicRegion(chromosome=region.chromosome, start=interval_bin[0],
                                 end=interval_bin[1], score=interval_bin[2])
                   for interval_bin in interval_bins]
        return br

    @staticmethod
    def bin_intervals(intervals, bins, interval_range=None, smoothing_window=None,
                      nan_replacement=None, zero_to_nan=False):
        if intervals is None:
            return []

        intervals = np.array(list(intervals))

        if interval_range is None:
            try:
                interval_range = (min(intervals[:, 0]), max(intervals[:, 1]))
            except (IndexError, TypeError):
                raise ValueError("intervals cannot be None or length 0 if not providing interval_range!")

        bin_size = (interval_range[1] - interval_range[0] + 1) / bins
        logger.debug("Bin size: {}".format(bin_size))

        return RegionBased._bin_intervals_equidist(intervals, bin_size, interval_range, bins=bins,
                                                   smoothing_window=smoothing_window,
                                                   nan_replacement=nan_replacement,
                                                   zero_to_nan=zero_to_nan)

    @staticmethod
    def bin_intervals_equidistant(intervals, bin_size, interval_range=None, smoothing_window=None,
                                  nan_replacement=None, zero_to_nan=False):
        if intervals is None:
            return []

        intervals = np.array(list(intervals))

        if interval_range is None:
            try:
                interval_range = (min(intervals[:, 0]), max(intervals[:, 1]))
            except (IndexError, TypeError):
                raise ValueError("intervals cannot be None or length 0 if not providing interval_range!")

        if isinstance(interval_range, GenomicRegion):
            interval_range = (interval_range.start, interval_range.end)

        return RegionBased._bin_intervals_equidist(intervals, bin_size, interval_range,
                                                   smoothing_window=smoothing_window,
                                                   nan_replacement=nan_replacement,
                                                   zero_to_nan=zero_to_nan)

    @staticmethod
    def _bin_intervals_equidist(intervals, bin_size, interval_range, bins=None, smoothing_window=None,
                                nan_replacement=None, zero_to_nan=False):
        if bins is None:
            bins = int((interval_range[1] - interval_range[0] + 1) / bin_size + .5)

        current_interval = 0
        bin_coordinates = []
        bin_weighted_sum = [0.0] * bins
        bin_weighted_count = [0.0] * bins
        bin_start = interval_range[0]
        for bin_counter in range(bins):
            bin_end = int(interval_range[0] + bin_size + (bin_size * bin_counter) + 0.5) - 1
            bin_coordinates.append((bin_start, bin_end))

            if current_interval < len(intervals):
                interval = intervals[current_interval]
            else:
                interval = None

            # add all successive, fully-contained intervals to bin
            while interval is not None and (interval[0] <= interval[1] <= bin_end and interval[1] >= bin_start):
                value = interval[2]
                if zero_to_nan and value < 10e-8:
                    value = np.nan

                if not np.isfinite(value):
                    if nan_replacement is not None:
                        value = nan_replacement
                    else:
                        value = None

                if value is not None:
                    f = (interval[1] + 1 - interval[0]) / bin_size
                    bin_weighted_sum[bin_counter] += f * value
                    bin_weighted_count[bin_counter] += f

                current_interval += 1
                if current_interval < len(intervals):
                    interval = intervals[current_interval]
                else:
                    interval = None

            # add partially-contained interval to bin
            if interval is not None and (interval[0] <= bin_end and interval[1] >= bin_start):
                value = interval[2]
                if zero_to_nan and value < 10e-8:
                    value = np.nan

                if not np.isfinite(value):
                    if nan_replacement is not None:
                        value = nan_replacement
                    else:
                        value = None

                if value is not None:
                    f = (min(bin_end, interval[1] + 1) - max(bin_start, interval[0])) / bin_size
                    bin_weighted_sum[bin_counter] += f * value
                    bin_weighted_count[bin_counter] += f

            bin_start = bin_end + 1

        with np.errstate(divide='ignore', invalid='ignore'):
            result = np.true_divide(bin_weighted_sum, bin_weighted_count)

            if nan_replacement is not None:
                result[~ np.isfinite(result)] = nan_replacement  # -inf inf NaN
            else:
                result[~ np.isfinite(result)] = np.nan  # -inf inf NaN

        if smoothing_window is not None:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                result = apply_sliding_func(result, smoothing_window)

        return tuple((bin_coordinates[i][0], bin_coordinates[i][1], result[i]) for i in range(len(result)))

    def binned_values(self, region, bins, smoothing_window=None, zero_to_nan=False):
        region = self._convert_region(region)
        return RegionBased.bin_intervals(self.region_intervals(region), bins,
                                         interval_range=(region.start, region.end),
                                         smoothing_window=smoothing_window,
                                         zero_to_nan=zero_to_nan)

    def binned_values_equidistant(self, region, bin_size, smoothing_window=None, zero_to_nan=False):
        region = self._convert_region(region)
        return RegionBased.bin_intervals_equidistant(self.region_intervals(region), bin_size,
                                                     interval_range=(region.start, region.end),
                                                     smoothing_window=smoothing_window,
                                                     zero_to_nan=zero_to_nan)

    def to_bed(self, file_name, subset=None, **kwargs):
        """
        Export regions as BED file
        """
        write_bed(file_name, self.regions(subset), **kwargs)

    def to_gff(self, file_name, subset=None, **kwargs):
        """
        Export regions as GFF file
        """
        write_gff(file_name, self.regions(subset), **kwargs)

    def to_bigwig(self, file_name, subset=None, **kwargs):
        """
        Export regions as BigWig file.
        """
        write_bigwig(file_name, self.regions(subset), mode='w', **kwargs)


class RegionWrapper(RegionBased):
    def __init__(self, regions):
        super(RegionWrapper, self).__init__()
        region_intervals = defaultdict(list)

        self._regions = []
        for region in regions:
            self._regions.append(region)
            interval = intervaltree.Interval(region.start - 1, region.end, data=region)
            region_intervals[region.chromosome].append(interval)

        self.region_trees = {}
        for chromosome, intervals in region_intervals.items():
            self.region_trees[chromosome] = intervaltree.IntervalTree(intervals)

    def _get_regions(self, item, *args, **kwargs):
        return self._regions[item]

    def _region_iter(self, *args, **kwargs):
        for region in self._regions:
            yield region

    def _region_subset(self, region, *args, **kwargs):
        tree = self.region_trees[region.chromosome]
        for interval in tree[region.start:region.end]:
            yield interval.data

    def _region_len(self):
        return len(self._regions)

    def chromosomes(self):
        return list(self.region_trees.keys())


class Bed(pybedtools.BedTool, RegionBased):
    """
    Data type representing a BED file.

    Only exists to support 'with' statements
    """

    def __init__(self, *args, **kwargs):
        pybedtools.BedTool.__init__(self, *args, **kwargs)
    
    def __exit__(self, exec_type, exec_val, exec_tb):
        pass

    def __enter__(self):
        return self

    def _region_iter(self, *args, **kwargs):
        for interval in self.intervals:
            yield self._interval_to_region(interval)

    def _get_regions(self, item, *args, **kwargs):
        if isinstance(item, string_types):
            item = GenomicRegion.from_string(item)

        if not isinstance(item, GenomicRegion):
            intervals = pybedtools.BedTool.__getitem__(self, item)
            if isinstance(intervals, pybedtools.Interval):
                return self._interval_to_region(intervals)
            elif isinstance(intervals, GenomicRegion):
                return intervals
            else:
                regions = []
                for interval in intervals:
                    regions.append(self._interval_to_region(interval))
                return regions

        start = item.start if item.start is not None else 1

        query_interval = pybedtools.cbedtools.Interval(chrom=item.chromosome,
                                                       start=start,
                                                       end=item.end)

        regions = []
        for interval in self.all_hits(query_interval):
            region = self._interval_to_region(interval)
            regions.append(region)
        return regions

    def _region_subset(self, region, *args, **kwargs):
        for interval in self.filter(lambda i: i.chrom == region.chromosome
                                    and i.start <= region.end
                                    and i.end >= region.start):
            yield self._interval_to_region(interval)

    def _region_len(self):
        return sum(1 for _ in self.intervals)

    def _interval_to_region(self, interval):
        try:
            score = float(interval.score)
        except (TypeError, ValueError):
            score = None

        if score is None:
            if len(interval.fields) == 4:  # likely bedGraph!
                try:
                    score = float(interval.fields[3])
                except ValueError:
                    score = np.nan
            else:
                score = np.nan

        try:
            name = interval.name
        except (TypeError, ValueError):
            warnings.warn("Pybedtools could not retrieve interval name. Continuing anyways.")
            name = None

        if self.intervals.file_type == 'gff':
            try:
                attributes = {key: value for key, value in interval.attrs.items()}
            except ValueError:
                attributes = {}

            attributes['chromosome'] = interval.chrom
            attributes['start'] = interval.start
            attributes['end'] = interval.end
            attributes['strand'] = interval.strand
            attributes['score'] = score
            attributes['fields'] = interval.fields
            attributes['source'] = interval.fields[1]
            attributes['feature'] = interval.fields[2]
            attributes['frame'] = interval.fields[7] if len(interval.fields) > 7 else '.'

            region = GenomicRegion(**attributes)
        else:
            region = GenomicRegion(chromosome=interval.chrom, start=interval.start, end=interval.end,
                                   strand=interval.strand, score=score, fields=interval.fields,
                                   name=name)
        return region

    def merge_overlapping(self, stat=_weighted_mean, sort=True):
        if sort:
            bed = self
        else:
            bed = self.sort()

        current_intervals = []
        for interval in bed:
            if len(current_intervals) == 0 or (current_intervals[-1].start < interval.end and
                                               current_intervals[-1].end > interval.start and
                                               current_intervals[-1].chrom == interval.chrom):
                current_intervals.append(interval)
            else:
                # merge
                intervals = np.array([(current.start, current.end,
                                       float(current.score) if current.score != '.' else np.nan)
                                      for current in current_intervals])
                merged_score = "{:0.6f}".format(stat(intervals))
                merged_strand = current_intervals[0].strand
                merged_start = min(intervals[:, 0])
                merged_end = max(intervals[:, 1])
                merged_chrom = current_intervals[0].chrom
                merged_name = current_intervals[0].name
                merged_interval = pybedtools.Interval(merged_chrom, merged_start, merged_end, name=merged_name,
                                                      score=merged_score, strand=merged_strand)
                current_intervals = [interval]
                yield merged_interval


class Bedpe(Bed):
    def __init__(self, *args, **kwargs):
        Bed.__init__(self, *args, **kwargs)

    @property
    def file_type(self):
        return 'bedpe'

    def _interval_to_region(self, interval):
        fields = interval.fields

        if len(fields) < 6:
            raise ValueError("File does not appear to be BEDPE (columns: {})".format(len(fields)))

        try:
            score = float(fields[7])
        except (IndexError, TypeError, ValueError):
            score = np.nan

        try:
            name = fields[6]
        except IndexError:
            name = '.'

        try:
            strand1 = fields[8]
        except IndexError:
            strand1 = '.'

        try:
            strand2 = fields[9]
        except IndexError:
            strand2 = '.'

        region = GenomicRegion(chromosome=fields[0], start=int(fields[1]), end=int(fields[2]),
                               chromosome1=fields[0], start1=int(fields[1]), end1=int(fields[2]),
                               chromosome2=fields[3], start2=int(fields[4]), end2=int(fields[5]),
                               strand=strand1, strand1=strand1, strand2=strand2,
                               score=score, fields=fields,
                               name=name)
        return region


class BigWig(RegionBased):
    def __init__(self, bw):
        RegionBased.__init__(self)
        if isinstance(bw, string_types):
            bw = pyBigWig.open(bw)
        self.bw = bw
        self._intervals = None

    @property
    def file_type(self):
        return 'bw'

    def __exit__(self, exec_type, exec_val, exec_tb):
        pass

    def __enter__(self):
        return self

    def _region_iter(self, *args, **kwargs):
        chromosome_lengths = self.chromosome_lengths
        chromosomes = self.chromosomes()
        for chromosome in chromosomes:
            for start, end, score in self.bw.intervals(chromosome, 1, chromosome_lengths[chromosome]):
                yield GenomicRegion(chromosome=chromosome, start=start+1, end=end, score=score)

    def _get_regions(self, item, *args, **kwargs):
        if isinstance(item, string_types):
            item = GenomicRegion.from_string(item)

        if not isinstance(item, GenomicRegion):
            return self.bw[item]

        return self.subset(item)

    def _region_subset(self, region, *args, **kwargs):
        if isinstance(region, GenomicRegion):
            regions = [region]
        else:
            regions = region

        for r in regions:
            for start, end, score in self.region_intervals(r):
                yield GenomicRegion(chromosome=r.chromosome, start=start, end=end, score=score)

    def _region_intervals(self, region, *args, **kwargs):
        if self._intervals is None:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                try:
                    intervals = self.bw.intervals(region.chromosome, region.start, region.end)
                except RuntimeError:
                    logger.debug("Invalid interval bounds? {}".format(region))
                    raise
        else:
            intervals = self._memory_intervals(region)

        interval_list = []
        if intervals is not None:
            for interval in intervals:
                interval_list.append((interval[0]+1, interval[1], interval[2]))
        return interval_list

    def _region_len(self):
        return sum(1 for _ in self.regions)

    def chromosomes(self):
        return natural_sort(list(self.chromosome_lengths.keys()))

    @property
    def chromosome_lengths(self):
        return self.bw.chroms()

    def __getattr__(self, name):
        try:
            func = getattr(self.__dict__['bw'], name)
            return func
        except AttributeError:
            if name == '__enter__':
                return BigWig.__enter__
            elif name == '__exit__':
                return BigWig.__exit__
            raise

    def __getitem__(self, item):
        return self._get_regions(item)

    def load_intervals_into_memory(self):
        self._intervals = dict()
        for chromosome in self.bw.chroms().keys():
            interval_starts = []
            interval_ends = []
            interval_values = []
            for start, end, score in self.bw.intervals(chromosome):
                interval_starts.append(start)
                interval_ends.append(end)
                interval_values.append(score)
            self._intervals[chromosome] = (interval_starts,
                                           interval_ends,
                                           interval_values)

    def _memory_intervals(self, region):
        all_intervals = self._intervals[region.chromosome]
        start_ix = bisect_right(all_intervals[0], region.start) - 1
        end_ix = bisect_left(all_intervals[1], region.end)
        return [(all_intervals[0][ix], all_intervals[1][ix], all_intervals[2][ix])
                for ix in range(max(0, start_ix), min(len(all_intervals[0]), end_ix+1))]

    def region_stats(self, region, bins=1, stat='mean'):
        if isinstance(region, string_types):
            region = GenomicRegion.from_string(region)

        chroms = self.bw.chroms()
        r_start = region.start - 1 if region.start is not None else 0
        r_end = region.end if region.end is not None else chroms[region.chromosome]

        return self.stats(region.chromosome, r_start, r_end, type=stat, nBins=bins)

    def intervals(self, region, bins=None, bin_size=None, smoothing_window=None,
                  nan_replacement=None, zero_to_nan=False, *args, **kwargs):
        return self.region_intervals(region, bins=bins, bin_size=bin_size, smoothing_window=smoothing_window,
                                     nan_replacement=nan_replacement, zero_to_nan=zero_to_nan,
                                     *args, **kwargs)


class Tabix(RegionBased):
    def __init__(self, file_name, preset=None):
        self._file_name = file_name
        self._file = pysam.TabixFile(file_name, parser=pysam.asTuple())
        RegionBased.__init__(self)

        self._file_type = self._get_file_extension()
        if preset is None:
            preset = self.file_type

        if isinstance(preset, string_types):
            if preset == 'gff' or preset == 'gtf':
                self._region_object = GffRegion
            elif preset == 'bed' or preset == 'bdg':
                self._region_object = BedRegion
            elif preset == 'vcf':
                self._region_object = BedRegion
            else:
                raise ValueError("Preset {} not valid".format(preset))
        else:
            self._region_object = preset

    @property
    def _estimate_region_bounds(self):
        return False

    @property
    def file_type(self):
        return self._file_type

    def _get_file_extension(self):
        fn = self._file_name
        if fn.endswith('.gz') or fn.endswith('.gzip'):
            fn = os.path.splitext(fn)[0]
        extension = os.path.splitext(fn)[1]
        return extension[1:]

    def _region_iter(self, *args, **kwargs):
        for chromosome in self.chromosomes():
            for region in self.subset(chromosome):
                yield region

    def _region_subset(self, region):
        try:
            for fields in self._file.fetch(region.chromosome, region.start, region.end):
                yield self._region_object(fields)
        except ValueError:
            if region.chromosome not in self.chromosomes():
                warnings.warn('{} not in list of contigs'.format(region.chromosome))
            else:
                raise

    def chromosomes(self):
        return self._file.contigs

    @staticmethod
    def to_tabix(file_name, preset=None, _tabix_path='tabix'):
        tabix_command = [_tabix_path]
        if preset is not None:
            tabix_command += ['-p', preset]
        tabix_command += file_name


class GenomicDataFrame(DataFrame):
    @property
    def regions(self):
        class RegionIter(object):
            def __init__(self, df):
                self.df = df

            def __iter__(self):
                for ix, row in self.df.iterrows():
                    yield self.df._row_to_region(row, ix=ix)

            def __call__(self):
                return iter(self)

        return RegionIter(self)

    def subset(self, region):
        for ix, row in self._sub_rows(region):
            yield self._row_to_region(row, ix=ix)

    def _sub_rows(self, region):
        if isinstance(region, string_types):
            region = GenomicRegion.from_string(region)

        if isinstance(region, GenomicRegion):
            regions = [region]
        else:
            regions = region

        for r in regions:
            if isinstance(r, string_types):
                r = GenomicRegion.from_string(r)

            query = ''
            if r.chromosome is not None:
                query += 'chromosome == "{}" and '.format(r.chromosome)
            if r.start is not None:
                query += 'start >= {} and '.format(r.start)
            if r.end is not None:
                query += 'end <= {} and '.format(r.end)
            query = query[:-5]

            sub_df = self.query(query)
            for ix, row in sub_df.iterrows():
                yield ix, row

    def _row_to_region(self, row, ix=None):
        attributes = {'ix': ix}
        for key, value in row.items():
            attributes[key] = value
        return GenomicRegion(**attributes)

    @classmethod
    def read_table(cls, file_name, **kwargs):
        return cls(read_table(file_name, **kwargs))


class Chromosome(object):
    """
    Chromosome data type.

    .. attribute:: name

        Name of the chromosome

    .. attribute:: length

        Length of the chromosome in base-pairs

    .. attribute:: sequence

        Base-pair sequence of DNA in the chromosome
    """
    def __init__(self, name=None, length=None, sequence=None):
        """
        Initialize chromosome

        :param name: Name of the chromosome
        :param length: Length of the chromosome in base-pairs
        :param sequence: Base-pair sequence of DNA in the chromosome
        """
        self.name = name.decode() if isinstance(name, bytes) else name
        self.length = length
        self.sequence = sequence.decode() if isinstance(sequence, bytes) else sequence
        if length is None and sequence is not None:
            self.length = len(sequence)
        if sequence is None and length is not None:
            self.length = length
            
    def __repr__(self):
        return "Name: %s\nLength: %d\nSequence: %s" % (self.name if self.name else '',
                                                       self.length if self.length else -1,
                                                       self.sequence[:20] + "..." if self.sequence else '')

    def __len__(self):
        """
        Get length of the chromosome.
        """
        return self.length
    
    def __getitem__(self, key):
        """
        Get object attributes by name
        """
        if key == 'name':
            return self.name
        if key == 'length':
            return self.length
        if key == 'sequence':
            return self.sequence
    
    @classmethod
    def from_fasta(cls, file_name, name=None, include_sequence=True):
        """
        Create a :class:`~Chromosome` from a FASTA file.

        This class method will load a FASTA file and convert it into
        a :class:`~Chromosome` object. If the FASTA file contains multiple
        sequences, only the first one will be read.

        :param file_name: Path to the FASTA file
        :param name: Chromosome name. If None (default), will be read
                     from the FASTA file header
        :param include_sequence: If True (default), stores the chromosome
                                 sequence in memory. Else, the sequence
                                 attribute will be set to None.
        :return: :class:`~Chromosome` if there is only a single FASTA
                 sequence in the file, list(:class:`~Chromosome`) if
                 there are multiple sequences.
        """
        with open(file_name, 'r') as fasta_file:
            fastas = SeqIO.parse(fasta_file, 'fasta')

            chromosomes = []
            for fasta in fastas:
                if include_sequence:
                    chromosome = cls(name if name else fasta.id, length=len(fasta), sequence=str(fasta.seq))
                else:
                    chromosome = cls(name if name else fasta.id, length=len(fasta))
                chromosomes.append(chromosome)

        if len(chromosomes) == 0:
            raise ValueError("File %s does not appear to be a FASTA file" % file_name)
        if len(chromosomes) == 1:
            return chromosomes[0]
        return chromosomes

    def get_restriction_sites(self, restriction_enzyme):
        """
        Find the restriction sites of a provided enzyme in this chromosome.

        Internally uses biopython to find RE sites.

        :param restriction_enzyme: The name of the restriction enzyme
                                   (e.g. HindIII)
        :return: List of RE sites in base-pairs (1-based)
        """
        logger.info("Calculating RE sites")
        try:
            re = getattr(Restriction, restriction_enzyme)
        except SyntaxError:
            raise ValueError("restriction_enzyme must be a string")
        except AttributeError:
            raise ValueError("restriction_enzyme string is not recognized: %s" % restriction_enzyme)

        return re.search(Seq.Seq(self.sequence))


class GenomicRegion(TableObject):
    """
    Class representing a genomic region.

    .. attribute:: chromosome

        Name of the chromosome this region is located on

    .. attribute:: start

        Start position of the region in base pairs

    .. attribute:: end

        End position of the region in base pairs

    .. attribute:: strand

        Strand this region is on (+1, -1)

    .. attribute:: ix

        Index of the region in the context of all genomic
        regions.

    """

    def __init__(self, start=None, end=None, chromosome=None, strand=None, ix=None, **kwargs):
        """
        Initialize this object.

        :param start: Start position of the region in base pairs
        :param end: End position of the region in base pairs
        :param chromosome: Name of the chromosome this region is located on
        :param strand: Strand this region is on (+1, -1)
        :param ix: Index of the region in the context of all genomic
                   regions.
        """
        self.start = start
        if end is None:
            end = start
        self.end = end
        if strand == "+":
            strand = 1
        elif strand == "-":
            strand = -1
        elif strand == "0" or strand == ".":
            strand = None
        self.strand = strand
        self.chromosome = chromosome.decode() if isinstance(chromosome, bytes) else chromosome
        self.ix = ix
        self.attributes = ['chromosome', 'start', 'end', 'strand', 'ix']

        for name, value in kwargs.items():
            setattr(self, name.decode() if isinstance(name, bytes) else name,
                    value.decode() if isinstance(value, bytes) else value)
            self.attributes.append(name)

    def set_attribute(self, attribute, value):
        setattr(self, attribute, value)
        if attribute not in self.attributes:
            self.attributes.append(attribute)

    @classmethod
    def from_row(cls, row):
        """
        Create a :class:`~GenomicRegion` from a PyTables row.
        """
        strand = row['strand']
        if strand == 0:
            strand = None
        return cls(start=row["start"], end=row["end"],
                   strand=strand, chromosome=row["chromosome"])

    @classmethod
    def from_string(cls, region_string):
        """
        Convert a string into a :class:`~GenomicRegion`.

        This is a very useful convenience function to quickly
        define a :class:`~GenomicRegion` object from a descriptor
        string.

        :param region_string: A string of the form
                              <chromosome>[:<start>-<end>[:<strand>]]
                              (with square brackets indicating optional
                              parts of the string). If any optional
                              part of the string is omitted, intuitive
                              defaults will be chosen.
        :return: :class:`~GenomicRegion`
        """
        chromosome = None
        start = None
        end = None
        strand = None
        
        # strip whitespace
        no_space_region_string = "".join(region_string.split())
        fields = no_space_region_string.split(':')
        
        if len(fields) > 3:
            raise ValueError("Genomic range string must be of the form <chromosome>[:<start>-<end>:[<strand>]]")
        
        # there is chromosome information
        if len(fields) > 0:
            chromosome = fields[0]
        
        # there is range information
        if len(fields) > 1 and fields[1] != '':
            start_end_bp = fields[1].split('-')
            if len(start_end_bp) > 0:
                start = str_to_int(start_end_bp[0])
            
            if len(start_end_bp) > 1:
                end = str_to_int(start_end_bp[1])

                if not end >= start:
                    raise ValueError("The end coordinate must be bigger than the start.")

        # there is strand information
        if len(fields) > 2:
            if fields[2] == '+' or fields[2] == '+1' or fields[2] == '1':
                strand = 1
            elif fields[2] == '-' or fields[2] == '-1':
                strand = -1
            else:
                raise ValueError("Strand only can be one of '+', '-', '+1', '-1', and '1'")
        return cls(start=start, end=end, chromosome=chromosome, strand=strand)
    
    def to_string(self):
        """
        Convert this :class:`~GenomicRegion` to its string representation.

        :return: str
        """
        region_string = ''
        if self.chromosome is not None:
            region_string += '%s' % self.chromosome
            
            if self.start is not None:
                region_string += ':%d' % self.start
                
                if self.end is not None:
                    region_string += '-%d' % self.end
                
                if self.strand is not None:
                    if self.strand == 1:
                        region_string += ':+'
                    else:
                        region_string += ':-'
        return region_string
    
    def __repr__(self):
        return self.to_string()

    def overlaps(self, region):
        """
        Check if this region overlaps with the specified region.

        :param region: :class:`~GenomicRegion` object or string
        """
        if isinstance(region, string_types):
            region = GenomicRegion.from_string(region)

        if region.chromosome != self.chromosome:
            return False

        if self.end is None or region.start is None or region.start <= self.end:
            if self.start is None or region.end is None or region.end >= self.start:
                return True
        return False

    def overlap(self, region):
        if region.chromosome != self.chromosome:
            return 0

        return max(0, min(self.end, region.end) - max(self.start, region.start))

    def contains(self, region):
        """
        Check if the specified region is completely contained in this region.

        :param region: :class:`~GenomicRegion` object or string
        """
        if isinstance(region, string_types):
            region = GenomicRegion.from_string(region)

        if region.chromosome != self.chromosome:
            return False

        if region.start >= self.start and region.end <= self.end:
            return True
        return False

    def _equals(self, region):
        if region.chromosome != self.chromosome:
            return False
        if region.start != self.start:
            return False
        if region.end != self.end:
            return False
        return True

    def is_reverse(self):
        if self.strand == -1 or self.strand == '-':
            return True
        return False

    def is_forward(self):
        if self.strand == 1 or self.strand == '+' or self.strand == 0:
            return True
        return False

    @property
    def strand_string(self):
        if self.is_forward():
            return '+'
        if self.is_reverse():
            return '-'
        return '.'

    @property
    def center(self):
        return self.start + (self.end - self.start)/2

    @property
    def five_prime(self):
        return self.end if self.is_reverse() else self.start

    @property
    def three_prime(self):
        return self.end if self.is_forward() else self.start

    def copy(self):
        d = {attribute: getattr(self, attribute) for attribute in self.attributes}
        return GenomicRegion(**d)

    def __eq__(self, other):
        return self._equals(other)

    def __ne__(self, other):
        return not self._equals(other)

    def __len__(self):
        return self.end - self.start

    def as_bed_line(self, score_field='score', name_field='name'):
        try:
            score = getattr(self, score_field)
        except AttributeError:
            warnings.warn("Score field {} does not exist, using '.'".format(score_field))
            score = '.'

        try:
            name = getattr(self, name_field)
        except AttributeError:
            warnings.warn("Name field {} does not exist, using '.'".format(name_field))
            name = '.'

        return "{}\t{}\t{}\t{}\t{}\t{}".format(self.chromosome, self.start, self.end,
                                               name, score, self.strand_string)

    def as_gff_line(self, source_field='source', feature_field='feature', score_field='score',
                    frame_field='frame', float_format='.2e'):
        try:
            source = getattr(self, source_field)
        except AttributeError:
            warnings.warn("Source field {} does not exist, using '.'".format(source_field))
            source = '.'

        try:
            feature = getattr(self, feature_field)
        except AttributeError:
            warnings.warn("Feature field {} does not exist, using '.'".format(feature_field))
            feature = '.'

        try:
            score = "{:{float_format}}".format(getattr(self, score_field), float_format=float_format)
        except AttributeError:
            warnings.warn("Score field {} does not exist, using '.'".format(score_field))
            score = '.'

        try:
            frame = getattr(self, frame_field)
        except AttributeError:
            warnings.warn("Frame field {} does not exist, using '.'".format(frame_field))
            frame = '.'

        no_group_items = {'start', 'end', 'chromosome', 'source', 'feature', 'score',
                          'frame', 'ix', 'strand', 'fields'}
        group = ''
        for attribute in self.attributes:
            if attribute not in no_group_items:
                a = getattr(self, attribute)
                if isinstance(a, float):
                    a = "{:{float_format}}".format(a, float_format=float_format)
                elif isinstance(a, string_types):
                    a = '"{}"'.format(a)
                group += '{} {}; '.format(attribute, a)

        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.chromosome, source,
                                                           feature, self.start + 1,
                                                           self.end, score,
                                                           self.strand_string, frame,
                                                           group)

    def expand(self,
               absolute=None, relative=None,
               absolute_left=0, absolute_right=0,
               relative_left=0.0, relative_right=0.0,
               copy=True, from_center=False):
        if absolute is not None:
            absolute_left, absolute_right = absolute, absolute
        if relative is not None:
            relative_left, relative_right = relative, relative

        extend_left_bp = absolute_left + int(relative_left * len(self))
        extend_right_bp = absolute_right + int(relative_right * len(self))

        new_region = self.copy() if copy else self
        if from_center:
            center = self.center
            new_region.start = int(center) - extend_left_bp
            new_region.end = int(center) + extend_right_bp
        else:
            new_region.start = int(self.start) - extend_left_bp
            new_region.end = int(self.end) + extend_right_bp
        return new_region

    def __add__(self, distance):
        new_region = self.copy()
        new_region.start += distance
        new_region.end += distance
        return new_region

    def __sub__(self, distance):
        return self.__add__(-distance)

    def fix_chromosome(self):
        if self.chromosome.startswith('chr'):
            self.chromosome = self.chromosome[3:]
        else:
            self.chromosome = 'chr' + self.chromosome


class LazyGenomicRegion(GenomicRegion):
    def __init__(self, row, ix=None, auto_update=True):
        self.reserved = {'_row', 'static_ix', 'strand', 'auto_update', 'attributes'}
        self._row = row
        self.static_ix = ix
        self.auto_update = auto_update

    def __getattr__(self, item):
        if item == 'reserved' or item in self.reserved:
            return object.__getattribute__(self, item)
        try:
            value = self._row[item]
            value = value.decode() if isinstance(value, bytes) else value
            return value
        except KeyError:
            raise AttributeError

    def __setattr__(self, key, value):
        if key == 'reserved' or key in self.reserved:
            super(LazyGenomicRegion, self).__setattr__(key, value)
        else:
            self._row[key] = value
            if self.auto_update:
                self.update()

    def update(self):
        self._row.update()

    @property
    def strand(self):
        try:
            return self._row["strand"]
        except KeyError:
            return None

    @property
    def ix(self):
        if self.static_ix is None:
            return self._row["ix"]
        return self.static_ix

    @property
    def attributes(self):
        return self._row.tables.colnames


class BedRegion(GenomicRegion):
    def __init__(self, bed_line, ix=None):
        try:
            self.fields = bed_line.split("\t")
        except AttributeError:
            self.fields = bed_line

        self.attributes = ('chromosome', 'start', 'end', 'name', 'score', 'strand',
                           'thick_start', 'thick_end', 'item_rgb', 'block_count',
                           'block_sizes', 'block_starts')[:len(self.fields)]
        self.ix = ix

    @property
    def chromosome(self):
        return self.fields[0]

    @property
    def start(self):
        return int(self.fields[1])

    @property
    def end(self):
        return int(self.fields[2])

    @property
    def name(self):
        try:
            return self.fields[3]
        except IndexError:
            return None

    @property
    def strand(self):
        try:
            s = self.fields[5]
            return -1 if s == '-' else 1
        except IndexError:
            return 1

    @property
    def score(self):
        try:
            return float(self.fields[4])
        except (IndexError, ValueError):
            return np.nan

    @property
    def thick_start(self):
        try:
            return int(self.fields[6])
        except IndexError:
            return None

    @property
    def thick_end(self):
        try:
            return int(self.fields[7])
        except IndexError:
            return None

    @property
    def item_rgb(self):
        try:
            return self.fields[8].split(',')
        except IndexError:
            return None

    @property
    def block_count(self):
        try:
            return int(self.fields[9])
        except IndexError:
            return None

    @property
    def block_sizes(self):
        try:
            return [int(s) for s in self.fields[10].split(',')]
        except IndexError:
            return None

    @property
    def block_sizes(self):
        try:
            return [int(s) for s in self.fields[11].split(',')]
        except IndexError:
            return None


class GffRegion(GenomicRegion):
    def __init__(self, gff_line, ix=None):
        try:
            self.fields = gff_line.split("\t")
        except AttributeError:
            self.fields = gff_line

        self.ix = ix
        self._attribute_dict = None

    @property
    def attributes(self):
        a = ['ix', 'chromosome', 'source', 'feature', 'start', 'end',
             'score', 'strand', 'frame']

        return a + list(self.attribute_dict.keys())

    @property
    def attribute_dict(self):
        if self._attribute_dict is None:
            self._attribute_dict = dict()

            attribute_fields = re.split(";\s*", self.fields[8])
            for field in attribute_fields:
                try:
                    key, value = shlex.split(field)
                except ValueError:
                    try:
                        key, value = re.split('=', field)
                    except ValueError:
                        continue
                self._attribute_dict[key] = value
        return self._attribute_dict

    def __getattr__(self, item):
        try:
            return self.attribute_dict[item]
        except (IndexError, KeyError):
            raise AttributeError("Attribute {} cannot be found".format(item))

    @property
    def seqname(self):
        return self.fields[0]

    @property
    def source(self):
        return self.fields[1]

    @property
    def feature(self):
        return self.fields[2]

    @property
    def chromosome(self):
        return self.seqname

    @property
    def start(self):
        return int(self.fields[3])

    @property
    def end(self):
        return int(self.fields[4])

    @property
    def name(self):
        return None

    @property
    def strand(self):
        try:
            s = self.fields[6]
            return -1 if s == '-' else 1
        except IndexError:
            return 1

    @property
    def score(self):
        try:
            return float(self.fields[4])
        except (IndexError, ValueError):
            return np.nan

    @property
    def frame(self):
        return self.fields[7]


class GenomicRegions(RegionBased):

    def __init__(self, regions=None):
        RegionBased.__init__(self)
        self._regions = []
        self._max_region_ix = -1

        if regions is not None:
            for region in regions:
                self.add_region(region)

    def add_region(self, region):
        """
        Add a genomic region to this object.

        This method offers some flexibility in the types of objects
        that can be loaded. See below for details.

        :param region: Can be a :class:`~GenomicRegion`, a str in the form
                       '<chromosome>:<start>-<end>[:<strand>], a dict with
                       at least the fields 'chromosome', 'start', and
                       'end', optionally 'ix', or a list of length 3
                       (chromosome, start, end) or 4 (ix, chromosome,
                       start, end).
        """
        ix = -1

        if isinstance(region, GenomicRegion):
            return self._add_region(copy.copy(region))
        elif isinstance(region, string_types):
            return self._add_region(GenomicRegion.from_string(region))
        elif type(region) is dict:
            return self._add_region(GenomicRegion(**copy.copy(region)))
        else:
            try:
                offset = 0
                if len(region) == 4:
                    ix = region[0]
                    offset += 1
                chromosome = region[offset]
                start = region[offset + 1]
                end = region[offset + 2]
                strand = 1
            except TypeError:
                raise ValueError("Node parameter has to be GenomicRegion, dict, or list")

        new_region = GenomicRegion(chromosome=chromosome, start=start, end=end, strand=strand, ix=ix)
        return self._add_region(new_region)

    def _add_region(self, region):
        region.ix = self._max_region_ix + 1

        self._regions.append(region)

        if region.ix > self._max_region_ix:
            self._max_region_ix = region.ix

        return self._region_len()

    def _region_len(self):
        return len(self._regions)

    def _region_iter(self, *args, **kwargs):
        for region in self._regions:
            yield region

    def _get_regions(self, key):
        return self._regions[key]

    @property
    def chromosome_bins(self):
        """
        Returns a dictionary of chromosomes and the start
        and end index of the bins they cover.

        Returned list is range-compatible, i.e. chromosome
        bins [0,5] cover chromosomes 1, 2, 3, and 4, not 5.
        """
        return self._chromosome_bins()

    def _chromosome_bins(self):
        chr_bins = {}
        for r in self.regions:
            if chr_bins.get(r.chromosome) is None:
                chr_bins[r.chromosome] = [r.ix, r.ix + 1]
            else:
                chr_bins[r.chromosome][1] = r.ix + 1
        return chr_bins

    @property
    def regions_dict(self):
        regions_dict = dict()
        for r in self.regions:
            regions_dict[r.ix] = r
        return regions_dict

    @property
    def bin_size(self):
        node = self.regions[0]
        return node.end - node.start + 1

    def distance_to_bins(self, distance):
        bin_size = self.bin_size
        bin_distance = int(distance/bin_size)
        if distance % bin_size > 0:
            bin_distance += 1
        return bin_distance

    def bins_to_distance(self, bins):
        return self.bin_size*bins


class RegionsTable(GenomicRegions, FileGroup):
    """
    PyTables Table wrapper for storing genomic regions.

    This class is inherited by objects working with lists of genomic
    regions, such as equi-distant bins along chromosomes in a genome
    (:class:`~Hic`) or restriction fragments of genomic DNA
    (:class:`~fanc.construct.seq.FragmentMappedReadPairs`)
    """

    _classid = 'REGIONSTABLE'

    class RegionDescription(t.IsDescription):
        """
        Description of a genomic region for PyTables Table
        """
        ix = t.Int32Col(pos=0)
        chromosome = t.StringCol(100, pos=1)
        start = t.Int64Col(pos=2)
        end = t.Int64Col(pos=3)
        strand = t.Int8Col(pos=4)
        _mask_ix = t.Int32Col(pos=5)

    def __init__(self, regions=None, file_name=None, mode='a',
                 additional_fields=None, _table_name_regions='regions',
                 tmpdir=None):
        """
        Initialize region table.

        :param data: List of regions to load in object. Can also
                     be used to load saved regions from file by
                     providing a path to an HDF5 file and setting
                     the file_name parameter to None.
        :param file_name: Path to a save file.
        :param _table_name_regions: (Internal) name of the HDF5
                                    node that stores data for this
                                    object
        """
        # parse potential unnamed argument
        if regions is not None:
            # data is file name
            if type(regions) is str or isinstance(regions, t.file.File):
                if file_name is None:
                    file_name = regions
                    regions = None

        try:
            FileGroup.__init__(self, _table_name_regions, file_name, mode=mode, tmpdir=tmpdir)
        except TypeError:
            logger.warning("RegionsTable is now a FileGroup-based object and "
                           "this object will no longer be compatible in the future")

        # check if this is an existing regions file
        try:
            group = self.file.get_node('/', _table_name_regions)

            if isinstance(group, t.table.Table):
                self._regions = group
            else:
                self._regions = self._group.regions

            if len(self._regions) > 0:
                self._max_region_ix = max(row['ix'] for row in self._regions.iterrows())
            else:
                self._max_region_ix = -1
        except t.NoSuchNodeError:
            basic_fields = dict()
            hidden_fields = dict()
            for field, description in RegionsTable.RegionDescription().columns.copy().items():
                if field.startswith('_'):
                    hidden_fields[field] = description
                else:
                    basic_fields[field] = description

            current = len(basic_fields)
            if additional_fields is not None:
                if not isinstance(additional_fields, dict) and issubclass(additional_fields, t.IsDescription):
                    # IsDescription subclass case
                    additional_fields = additional_fields.columns

                # add additional user-defined fields
                for key, value in sorted(additional_fields.items(),
                                         key=lambda x: x[1]._v_pos if x[1]._v_pos is not None else 1):
                    if key not in basic_fields:
                        value._v_pos = current
                        current += 1
                        basic_fields[key] = value
            # add hidden fields
            for key, value in sorted(hidden_fields.items(),
                                     key=lambda x: x[1]._v_pos if x[1]._v_pos is not None else 1):
                value._v_pos = current
                current += 1
                basic_fields[key] = value

            self._regions = t.Table(self._group, 'regions', basic_fields)
            self._max_region_ix = -1

        # index regions table
        create_col_index(self._regions.cols.ix)
        create_col_index(self._regions.cols.start)
        create_col_index(self._regions.cols.end)

        self._ix_to_chromosome = dict()
        self._chromosome_to_ix = dict()

        if regions is not None:
            self.add_regions(regions)
        else:
            self._update_references()

    def flush(self):
        self._regions.flush()
        self._update_references()

    def add_region(self, region, flush=True):
        # super-method, calls below '_add_region'
        ix = GenomicRegions.add_region(self, region)
        if flush:
            self.flush()
            self._update_references()
        return ix

    def _add_region(self, region):
        ix = self._max_region_ix + 1
        
        # actually append
        row = self._regions.row
        row['ix'] = ix
        row['chromosome'] = region.chromosome
        row['start'] = region.start
        row['end'] = region.end
        if hasattr(region, 'strand') and region.strand is not None:
            row['strand'] = region.strand

        for name in self._regions.colnames[5:]:
            if hasattr(region, name):
                row[name] = getattr(region, name)

        row.append()
        
        if ix > self._max_region_ix:
            self._max_region_ix = ix

        return ix

    def _update_references(self):
        chromosomes = []
        for region in self.regions(lazy=True):
            if len(chromosomes) == 0 or chromosomes[-1] != region.chromosome:
                chromosomes.append(region.chromosome)

        for i, chromosome in enumerate(chromosomes):
            self._ix_to_chromosome[i] = chromosome
            self._chromosome_to_ix[chromosome] = i

    def add_regions(self, regions):
        """
        Bulk insert multiple genomic regions.

        :param regions: List (or any iterator) with objects that
                        describe a genomic region. See
                        :class:`~RegionsTable.add_region` for options.
        """
        try:
            l = len(regions)
            _log = True
        except TypeError:
            l = None
            _log = False

        pb = RareUpdateProgressBar(max_value=l, silent=config.hide_progressbars)
        if _log:
            pb.start()

        for i, region in enumerate(regions):
            self.add_region(region, flush=False)
            if _log:
                pb.update(i)
        if _log:
            pb.finish()

        self._regions.flush()
        self._update_references()

    def data(self, key, value=None):
        """
        Retrieve or add vector-data to this object. If there is exsting data in this
        object with the same name, it will be replaced

        :param key: Name of the data column
        :param value: vector with region-based data (one entry per region)
        """
        if key not in self._regions.colnames:
            raise KeyError("%s is unknown region attribute" % key)

        if value is not None:
            for i, row in enumerate(self._regions):
                row[key] = value[i]
                row.update()
            self._regions.flush()

        return (row[key] for row in self._regions)

    def _get_region_ix(self, region):
        """
        Get index from other region properties (chromosome, start, end)
        """
        condition = "(start == %d) & (end == %d) & (chromosome == b'%s')"
        condition %= region.start, region.end, region.chromosome
        for res in self._regions.where(condition):
            return res["ix"]
        return None

    def _row_to_region(self, row, lazy=False, auto_update=True):
        if lazy:
            return LazyGenomicRegion(row, auto_update=auto_update)

        kwargs = {}
        for name in self._regions.colnames:
            if name not in RegionsTable.RegionDescription().columns.keys():
                value = row[name]
                value = value.decode() if isinstance(value, bytes) else value
                kwargs[name] = value

        try:
            mask_ix = row['_mask_ix']
        except (KeyError, ValueError):
            mask_ix = 0

        return GenomicRegion(chromosome=row["chromosome"].decode(), start=row["start"],
                             end=row["end"], ix=row["ix"], _mask_ix=mask_ix, **kwargs)

    def _region_iter(self, lazy=False, auto_update=True, *args, **kwargs):
        for row in self._regions:
            yield self._row_to_region(row, lazy=lazy, auto_update=auto_update)

    def _region_subset(self, region, lazy=False, auto_update=True, *args, **kwargs):
        for row in self._subset_rows(region):
            sub_region = self._row_to_region(row, lazy=lazy, auto_update=auto_update)
            yield sub_region

    def _get_regions(self, key):
        res = self._regions[key]

        if isinstance(res, np.ndarray):
            regions = []
            for region in res:
                regions.append(self._row_to_region(region))
            return regions
        else:
            return self._row_to_region(res)

    def _region_len(self):
        return len(self._regions)

    def _subset_rows(self, key):
        """
        Iterate over a subset of regions given the specified key.

        :param key: A :class:`~fanc.data.genomic.GenomicRegion` object,
                    or a list of the former. Also accepts slices and integers
        :return: Iterator over the specified subset of regions
        """
        if isinstance(key, slice):
            for row in self._regions.where("(ix >= {}) & (ix < {})".format(key.start, key.stop)):
                yield row
        elif isinstance(key, int):
            yield self._regions[key]
        elif isinstance(key, list) and len(key) > 0 and isinstance(key[0], int):
            for ix in key:
                yield self._regions[ix]
        else:
            if isinstance(key, string_types):
                key = GenomicRegion.from_string(key)

            if isinstance(key, GenomicRegion):
                keys = [key]
            else:
                keys = key

            for k in keys:
                if isinstance(k, string_types):
                    k = GenomicRegion.from_string(k)

                query = '('
                if k.chromosome is not None:
                    query += "(chromosome == b'%s') & " % k.chromosome
                if k.end is not None:
                    query += "(start <= %d) & " % k.end
                if k.start is not None:
                    query += "(end >= %d) & " % k.start
                if query.endswith(' & '):
                    query = query[:-3]
                query += ')'

                if len(query) == 2:
                    for row in self._regions:
                        yield row
                else:
                    for row in self._regions.where(query):
                        yield row

    def region_bins(self, region):
        start_ix = None
        end_ix = None
        for r in self.regions(region):
            if start_ix is None:
                start_ix = r.ix
            end_ix = r.ix + 1
        return slice(start_ix, end_ix, 1)


class Genome(FileGroup):
    """
    Class representing a collection of chromosomes.

    Extends the :class:`~RegionsTable` class and provides
    all the expected functionality. Provides some convenience batch
    methods that call :class:`~Chromosome` methods for every
    chromosome in this object.

    This object can be saved to file.
    """

    class ChromosomeDefinition(t.IsDescription):
        name = t.StringCol(255, pos=0)
        length = t.Int64Col(pos=1)

    def __init__(self, file_name=None, chromosomes=None, mode='a', tmpdir=None,
                 _table_name_chromosomes='chromosomes'):
        """
        Build :class:`~Genome` from a list of chromosomes or load
        previously saved object.

        :param file_name: Path of the file to load or to save to.
        :param chromosomes: List of chromosomes to load into this
                            object.
        """
        FileGroup.__init__(self, _table_name_chromosomes, file_name, mode=mode, tmpdir=tmpdir)

        # check if this is an existing regions file
        try:
            self._sequences = self._group.sequences
        except t.NoSuchNodeError:
            self._sequences = self.file.create_vlarray(self._group, 'sequences', t.VLStringAtom())

        try:
            self._chromosome_table = self._group.chromosomes
        except t.NoSuchNodeError:
            try:
                self._chromosome_table = self.file.create_table(self._group, 'chromosomes',
                                                                Genome.ChromosomeDefinition)
            except t.FileModeError:
                self._chromosome_table = None
                pass

        if chromosomes is not None:
            if isinstance(chromosomes, Chromosome):
                chromosomes = [chromosomes]

            for chromosome in chromosomes:
                self.add_chromosome(chromosome)

    def chromosomes(self):
        return self._names

    @property
    def _names(self):
        if self._chromosome_table is None:
            try:
                return self.meta['chromosome_names']
            except KeyError:
                return []
        else:
            return [row['name'].decode() if isinstance(row['name'], bytes) else row['name']
                    for row in self._chromosome_table]

    @_names.setter
    def _names(self, names):
        if self._chromosome_table is None:
            self.meta['chromosome_names'] = names
        else:
            counter = 0
            for i, row in enumerate(self._chromosome_table):
                row['name'] = names[i]
                row.update()
                counter += 1
            for i in range(counter, len(names)):
                self._chromosome_table.row['name'] = names[i]
                self._chromosome_table.row.append()
        self._chromosome_table.flush()

    @property
    def _lengths(self):
        if self._chromosome_table is None:
            try:
                return self.meta['chromosome_lengths']
            except KeyError:
                return []
        else:
            return [row['length'] for row in self._chromosome_table]

    @_lengths.setter
    def _lengths(self, lengths):
        if self._chromosome_table is None:
            self.meta['chromosome_lengths'] = lengths
        else:
            counter = 0
            for i, row in enumerate(self._chromosome_table):
                row['length'] = lengths[i]
                row.update()
                counter += 1
            for i in range(counter, len(lengths)):
                self._chromosome_table.row['length'] = lengths[i]
                self._chromosome_table.row.append()
        self._chromosome_table.flush()
        self.meta['chromosome_lengths'] = lengths

    @classmethod
    def from_folder(cls, folder_name, file_name=None, exclude=None,
                    include_sequence=True, tmpdir=None):
        """
        Load every FASTA file from a folder as a chromosome.

        :param folder_name: Path to the folder to load
        :param file_name: File to save Genome object to
        :param exclude: List or set of chromosome names that
                        should NOT be loaded
        :param include_sequence: If True, will save the
                                 chromosome sequences in the
                                 Genome object
        """
        chromosomes = []
        folder_name = os.path.expanduser(folder_name)
        for f in os.listdir(folder_name):
            try:
                chromosome = Chromosome.from_fasta(folder_name + "/" + f,
                                                   include_sequence=include_sequence)
                logger.info("Adding chromosome %s" % chromosome.name)
                if exclude is None:
                    chromosomes.append(chromosome)
                elif chromosome.name not in exclude:
                    chromosomes.append(chromosome)
            except (ValueError, IOError):
                pass

        return cls(chromosomes=chromosomes, file_name=file_name, tmpdir=tmpdir)

    @classmethod
    def from_string(cls, genome_string, file_name=None, tmpdir=None, mode='a'):
        """
        Convenience function to load a :class:`~Genome` from a string.

        :param genome_string: Path to FASTA file, path to folder with
                              FASTA files, comma-separated list of
                              paths to FASTA files, path to HDF5 file
        :param file_name: Path to save file
        :return: A :class:`~Genome` object
        """
        # case 1: FASTA file = Chromosome
        if is_fasta_file(genome_string):
            chromosomes = Chromosome.from_fasta(genome_string)
            genome = cls(chromosomes=chromosomes, file_name=file_name, tmpdir=tmpdir)
        # case 2: Folder with FASTA files
        elif os.path.isdir(genome_string):
            genome = cls.from_folder(genome_string, file_name=file_name, tmpdir=tmpdir)
        # case 3: path to HDF5 file
        elif os.path.isfile(genome_string):
            genome = cls(genome_string, tmpdir=tmpdir, mode=mode)
        # case 4: List of FASTA files
        else:
            chromosome_files = genome_string.split(',')
            chromosomes = []
            for chromosome_file in chromosome_files:
                chromosome = Chromosome.from_fasta(os.path.expanduser(chromosome_file))
                chromosomes.append(chromosome)
            genome = cls(chromosomes=chromosomes, file_name=file_name, tmpdir=tmpdir)

        return genome

    def __getitem__(self, key):
        """
        Get Genome table subsets.

        If the result is one or more rows, they will be converted to
        :class:`~Chromosome` objects, if the result is a column, it
        will be returned without conversion.
        """
        names = self._names
        lengths = self._lengths

        if isinstance(key, string_types):
            key = names.index(key)

        if isinstance(key, int):
            return Chromosome(name=names[key], length=lengths[key], sequence=self._sequences[key])
        elif isinstance(key, slice):
            l = []
            start = key.start if key.start is not None else 0
            stop = key.stop if key.stop is not None else len(names)
            step = key.step if key.step is None else 1

            for i in range(start, stop, step):
                c = Chromosome(name=names[i], length=lengths[i], sequence=self._sequences[i])
                l.append(c)
            return l
        else:
            l = []
            for i in key:
                if isinstance(i, string_types):
                    i = names.index(i)
                c = Chromosome(name=names[i], length=lengths[i], sequence=self._sequences[i])
                l.append(c)
            return l

    def __len__(self):
        return len(self._names)

    def __iter__(self):
        """
        Get iterator over :class:`~Chromosome` objects.
        """
        this = self

        class Iter(object):
            def __init__(self):
                self.current = 0

            def __iter__(self):
                self.current = 0
                return self

            def __next__(self):
                if self.current >= len(this):
                    raise StopIteration
                self.current += 1
                return this[self.current - 1]

        return Iter()

    def add_chromosome(self, chromosome):
        """
        Add a :class:`~Chromosome` to this object.

        Will choose suitable defaults for missing attributes.

        :param chromosome: :class:`~Chromosome` object or similar
                           object (e.g. dict) with the same fields
        """
        i = len(self._names)

        n = str(i)
        if chromosome.name is not None:
            n = chromosome.name

        l = 0
        if chromosome.length is not None:
            l = chromosome.length

        s = ''
        if chromosome.sequence is not None:
            s = chromosome.sequence
            if l == 0:
                l = len(s)

        self._chromosome_table.row['name'] = n
        self._chromosome_table.row['length'] = l
        self._chromosome_table.row.append()

        # self._names = self._names + [n]
        # self._lengths = self._lengths + [l]
        self._sequences.append(s)
        self._sequences.flush()
        self._chromosome_table.flush()

    def get_regions(self, split, file_name=None, chromosomes=None):
        """
        Extract genomic regions from genome.

        Provides two options:

        - Splits chromosomes at restriction sites if the split
          parameter is the name of a restriction enzyme.

        - Splits chromosomes at equi-distant points if split
          is an integer

        :param split: Name of a restriction enzyme or positive
                      integer
        :param file_name: Name of a file if the result of this
                          method should be saved to file
        :param chromosomes: List of chromosome names to include. Default: all
        :return: :class:`~GenomicRegions`
        """

        regions = RegionsTable(file_name=file_name)
        for chromosome in self:
            if chromosomes is not None and chromosome.name not in chromosomes:
                continue
            split_locations = []
            if isinstance(split, string_types):
                split_locations = chromosome.get_restriction_sites(split)
            elif isinstance(split, int):
                for i in range(split, len(chromosome) - 1, split):
                    split_locations.append(i)
            else:
                for i in split:
                    split_locations.append(i)

            for i in range(0, len(split_locations)):
                if i == 0:
                    region = GenomicRegion(start=1, end=split_locations[i], chromosome=chromosome.name)
                else:
                    region = GenomicRegion(start=split_locations[i - 1] + 1,
                                           end=split_locations[i], chromosome=chromosome.name)

                regions.add_region(region, flush=False)

            # add last node
            if len(split_locations) > 0:
                region = GenomicRegion(start=split_locations[len(split_locations) - 1] + 1,
                                       end=chromosome.length, chromosome=chromosome.name)
            else:
                region = GenomicRegion(start=1, end=chromosome.length, chromosome=chromosome.name)
            regions.add_region(region, flush=True)

        return regions

    def sub_sequence(self, chromosome, start=None, end=None):
        if start is not None:
            selection_region = GenomicRegion(chromosome=chromosome, start=start, end=end)
        elif isinstance(chromosome, GenomicRegion):
            selection_region = chromosome
        else:
            selection_region = GenomicRegion.from_string(chromosome)

        res_chromosome = self[selection_region.chromosome]
        if selection_region.start is None:
            return res_chromosome.sequence
        return res_chromosome.sequence[selection_region.start - 1:selection_region.end]


class Node(GenomicRegion, TableObject):
    """
    Class representing a node in a :class:`~Hic` object.

    Backed by a :class:`~GenomicRegion`, this class additionally
    provides methods to access the node index in the context of
    the :class:`~Hic` object.

    .. attribute:: chromosome

        Name of the chromosome this region is located on

    .. attribute:: start

        Start position of the region in base pairs

    .. attribute:: end

        End position of the region in base pairs

    .. attribute:: strand

        Strand this region is on (+1, -1)

    .. attribute:: ix

        Index of the region in the context of all genomic
        regions.
    """
    def __init__(self, chromosome=None, start=None, end=None, ix=None):
        self.ix = ix
        super(Node, self).__init__(chromosome=chromosome, start=start, end=end, ix=ix)
    
    def __repr__(self):
        if self.ix is None:
            return "%s, %d-%d" % (self.chromosome, self.start, self.end)
        else:
            return "%d: %s, %d-%d" % (self.ix, self.chromosome, self.start, self.end)


class LazyNode(LazyGenomicRegion, Node):
    def __init__(self, row, ix=None):
        LazyGenomicRegion.__init__(self, row=row, ix=ix)


class Edge(TableObject):
    """
    A contact / an Edge between two genomic regions.

    .. attribute:: source

        The index of the "source" genomic region. By convention,
        source <= sink.

    .. attribute:: sink

        The index of the "sink" genomic region.

    .. attribute:: weight

        The weight or contact strength of the edge. Can, for
        example, be the number of reads mapping to a contact.
    """
    def __init__(self, source, sink, **kwargs):
        """
        :param source: The index of the "source" genomic region
                       or :class:`~Node` object.
        :param sink: The index of the "sink" genomic region
                     or :class:`~Node` object.
        :param data: The weight or of the edge or a dictionary with
                     other fields
        """
        self._source = source
        self._sink = sink
        self.field_names = []

        for key, value in kwargs.items():
            setattr(self, key.decode() if isinstance(key, bytes) else key, value)
            self.field_names.append(key)

    @property
    def source(self):
        if isinstance(self._source, GenomicRegion):
            return self._source.ix
        return self._source

    @property
    def sink(self):
        if isinstance(self._sink, GenomicRegion):
            return self._sink.ix
        return self._sink

    @property
    def source_node(self):
        if isinstance(self._source, GenomicRegion):
            return self._source
        raise RuntimeError("Source not not provided during object initialization!")

    @property
    def sink_node(self):
        if isinstance(self._sink, GenomicRegion):
            return self._sink
        raise RuntimeError("Sink not not provided during object initialization!")

    def __repr__(self):
        base_info = "%d--%d" % (self.source, self.sink)
        for field in self.field_names:
            base_info += "\n\t%s: %s" % (field, str(getattr(self, field)))
        return base_info + "\n"


class LazyEdge(Edge):
    def __init__(self, row, nodes_table=None, auto_update=True):
        self.reserved = {'_row', '_nodes_table', 'auto_update', '_source_node', '_sink_node'}
        self._row = row
        self._nodes_table = nodes_table
        self.auto_update = auto_update
        self._source_node = None
        self._sink_node = None

    def _set_item(self, item, value):
        self._row[item] = value
        if self.auto_update:
            self.update()

    def __getattr__(self, item):
        if item == 'reserved' or item in self.reserved:
            return object.__getattribute__(self, item)
        try:
            value = self._row[item]
            value = value.decode() if isinstance(value, bytes) else value
            return value
        except KeyError:
            raise AttributeError("Attribute not supported (%s)" % str(item))

    def __setattr__(self, key, value):
        if key == 'reserved' or key in self.reserved:
            super(LazyEdge, self).__setattr__(key, value)
        else:
            self._row[key] = value
            if self.auto_update:
                self.update()

    def update(self):
        self._row.update()

    @property
    def source(self):
        return self._row['source']

    @property
    def sink(self):
        return self._row['sink']

    @property
    def source_node(self):
        if self._nodes_table is None:
            raise RuntimeError("Must set the _nodes_table attribute before calling this method!")

        if self._source_node is None:
            source_row = self._nodes_table[self.source]
            return LazyNode(source_row)
        return self._source_node

    @property
    def sink_node(self):
        if self._nodes_table is None:
            raise RuntimeError("Must set the _nodes_table attribute before calling this method!")

        if self._sink_node is None:
            sink_row = self._nodes_table[self.sink]
            return LazyNode(sink_row)
        return self._sink_node

    def __repr__(self):
        base_info = "%d--%d" % (self.source, self.sink)
        return base_info


class RegionPairs(Maskable, RegionsTable):
    """
    Class for working with data associated with pairs of regions.

    Generally, a RegionPairs object has two components:

    - Nodes or regions: (Non-overlapping) genomic regions
      obtained by splitting the genome into distinct pieces.
      See also :class:`~GenomicRegion` and :class:`~RegionsTable`

    - Edges or contacts: Pairs of genomic regions. See also
      :class:`~Edge`
    """

    _classid = 'REGIONPAIRS'

    class EntryDescription(t.IsDescription):
        source = t.Int32Col(pos=0)
        sink = t.Int32Col(pos=1)

    class EdgeIter(object):
        def __init__(self, this, _iter=None):
            self.this = this
            if _iter is None:
                self.iter = iter(this._edges)
            else:
                self.iter = iter(_iter)
            self.row_conversion_args = list()
            self.row_conversion_kwargs = dict()
            self.only_intrachromosomal = False
            self.regions_dict = None

        def __getitem__(self, item):
            res = self.this._edges[item]

            if isinstance(res, np.ndarray):
                edges = []
                for edge in res:
                    edges.append(self.this._row_to_edge(edge, *self.row_conversion_args, **self.row_conversion_kwargs))
                return edges
            else:
                edge = self.this._row_to_edge(res, *self.row_conversion_args, **self.row_conversion_kwargs)
                return edge

        def __iter__(self):
            if self.only_intrachromosomal:
                self.regions_dict = self.this.regions_dict
            return self

        def __call__(self, *args, **kwargs):
            if 'only_intrachromosomal' in kwargs:
                self.only_intrachromosomal = kwargs['only_intrachromosomal']
                del kwargs['only_intrachromosomal']
            self.row_conversion_args = args
            self.row_conversion_kwargs = kwargs
            key = kwargs.get('key', None)

            if key is not None:
                self.iter = self.this._edge_subset_rows(key=key,
                                                        only_intrachromosomal=self.only_intrachromosomal)
            return iter(self)

        def __next__(self):
            row = next(self.iter)
            if self.only_intrachromosomal:
                while self.regions_dict[row['source']].chromosome != self.regions_dict[row['sink']].chromosome:
                    row = next(self.iter)
            return self.this._row_to_edge(row, *self.row_conversion_args, **self.row_conversion_kwargs)

        def __len__(self):
            return len(self.this._edges)

    def __init__(self, file_name=None, mode='a', additional_fields=None, tmpdir=None,
                 _table_name_nodes='nodes', _table_name_edges='edges', _edge_buffer_size=1000000):

        """
        Initialize a :class:`~RegionPairs` object.

        :param file_name: Path to a save file
        :param mode: File mode to open underlying file
        :param additional_fields: Additional fields (in PyTables notation) associated with
                                  edge data, e.g. {'weight': tables.Float32Col()}
        :param _table_name_nodes: (Internal) name of the HDF5 node for regions
        :param _table_name_edges: (Internal) name of the HDF5 node for edges
        """

        # private variables
        self._max_node_ix = -1

        if file_name is not None:
            file_name = os.path.expanduser(file_name)

        # initialize inherited objects
        RegionsTable.__init__(self, file_name=file_name, _table_name_regions=_table_name_nodes,
                              mode=mode, tmpdir=tmpdir)
        Maskable.__init__(self, self.file)

        # create edge table
        if _table_name_edges in self.file.root:
            self._edges = self.file.get_node('/', _table_name_edges)
        else:
            basic_fields = self._get_field_dict(additional_fields=additional_fields)

            self._edges = MaskedTable(self.file.root, _table_name_edges, basic_fields)

        # index edge table
        create_col_index(self._edges.cols.source)
        create_col_index(self._edges.cols.sink)

        # update field names
        self._source_field_ix = 0
        self._sink_field_ix = 0
        self.field_names = []
        self._field_names_dict = dict()
        self._edge_field_defaults = dict()
        self._field_dict = self._edges.coldescrs
        for i, name in enumerate(self._edges.colnames):
            if not name.startswith("_"):
                self.field_names.append(name)
            if name == 'source':
                self._source_field_ix = i
            if name == 'sink':
                self._sink_field_ix = i
            self._field_names_dict[name] = i
            self._edge_field_defaults[name] = self._edges.coldescrs[name].dflt

        self._edge_buffer = []
        self._edge_buffer_size = _edge_buffer_size

    def disable_indexes(self):
        try:
            self._edges.cols.source.remove_index()
        except:
            pass

        try:
            self._edges.cols.sink.remove_index()
        except:
            pass
        self._edges.disable_mask_index()

    def enable_indexes(self):
        create_col_index(self._edges.cols.source)
        create_col_index(self._edges.cols.sink)
        self._edges.enable_mask_index()

    def _get_field_dict(self, additional_fields=None):
        basic_fields = RegionMatrixTable.EntryDescription().columns.copy()
        if additional_fields is not None:
            if not isinstance(additional_fields, dict) and issubclass(additional_fields, t.IsDescription):
                # IsDescription subclass case
                additional_fields = additional_fields.columns

            current = len(basic_fields)
            for key, value in sorted(additional_fields.items(), key=lambda x: x[1]._v_pos):
                if key not in basic_fields:
                    if value._v_pos is not None:
                        value._v_pos = current
                        current += 1
                    basic_fields[key] = value
        return basic_fields

    def add_node(self, node, flush=True):
        """
        Add a :class:`~Node` or :class:`~GenomicRegion`.

        :param node: :class:`~Node` or :class:`~GenomicRegion`,
                     see :func:`~RegionsTable.add_region` for details
        :param flush: Write data to file immediately after import.
        """
        return self.add_region(node, flush)

    def add_edge(self, edge, check_nodes_exist=True, flush=True, replace=False, row=None):
        """
        Add an edge to this object.

        :param edge: :class:`~Edge`, dict with at least the
                     attributes source and sink, optionally weight,
                     or a list of length 2 (source, sink) or 3
                     (source, sink, weight).
        :param check_nodes_exist: Make sure that there are nodes
                                  that match source and sink indexes
        :param flush: Write data to file immediately after import
        :param replace: If row is provided, replace values in existing edge with the ones in edge
        :param row: PyTables row object representing an edge. If provided, edge will be used to
                    modify existing row.
        """
        if not isinstance(edge, Edge):
            source = None
            sink = None

            # object
            is_object = True
            try:
                source = edge.source
                sink = edge.sink
            except AttributeError:
                is_object = False

            # dictionary
            is_dict = False
            if not is_object:
                is_dict = True
                try:
                    source = edge['source']
                    sink = edge['sink']
                except TypeError:
                    is_dict = False

            # list
            is_list = False
            if not is_object and not is_dict:
                is_list = True
                try:
                    source = edge[self._source_field_ix]
                    sink = edge[self._sink_field_ix]
                except TypeError:
                    is_list = False

            if source is None and sink is None:
                raise ValueError("Edge type not recognised (%s)" % str(type(edge)))

            if is_object:
                new_edge = self._edge_from_object(edge)
            elif is_dict:
                new_edge = self._edge_from_dict(edge)
            elif is_list:
                new_edge = self._edge_from_list(edge)
            else:
                raise ValueError("Edge type not recognised (%s)" % str(type(edge)))
        else:
            new_edge = edge

        if check_nodes_exist:
            n_regions = len(self._regions)
            if new_edge.source >= n_regions or new_edge.sink >= n_regions:
                raise ValueError("Node index exceeds number of nodes in object")

        self._add_edge(new_edge, row=row, replace=replace)

        if flush:
            self.flush()

    def _add_edge(self, edge, row, replace=False):
        source, sink = edge.source, edge.sink
        if source > sink:
            source, sink = sink, source

        if row is None:
            record = [None] * len(self._field_names_dict)
            for name, ix in self._field_names_dict.items():
                try:
                    record[ix] = getattr(edge, name)
                except AttributeError:
                    record[ix] = self._edge_field_defaults[name]
            record[self._field_names_dict['source']] = source
            record[self._field_names_dict['sink']] = sink

            self._edge_buffer.append(tuple(record))
            if len(self._edge_buffer) > self._edge_buffer_size:
                self._edges.append(self._edge_buffer)
                self._edge_buffer = []
        else:
            update = True
            if row is None:
                update = False
                row = self._edges.row
            row['source'] = source
            row['sink'] = sink
            for name in self.field_names:
                if not name == 'source' and not name == 'sink':
                    try:
                        value = getattr(edge, name)
                        if replace or not update:
                            row[name] = value
                        else:
                            row[name] += value
                    except AttributeError:
                        pass
            if update:
                row.update()
            else:
                row.append()

    def _edge_from_object(self, edge):
        return edge

    def _edge_from_dict(self, edge):
        source, sink = edge['source'], edge['sink']

        attributes = dict()
        for name, value in edge.items():
            if not name == 'source' and not name == 'sink':
                attributes[name] = value

        return Edge(source, sink, **attributes)

    def _edge_from_list(self, edge):
        source, sink = edge[self._source_field_ix], edge[self._sink_field_ix]

        attributes = dict()
        for i, name in enumerate(self.field_names):
            if not name == 'source' and not name == 'sink':
                try:
                    attributes[name] = edge[i]
                except IndexError:
                    break

        return Edge(source, sink, **attributes)

    def _default_edge_list(self):
        record = [None] * len(self._field_names_dict)
        for name, ix in self._field_names_dict.items():
            record[ix] = self._edge_field_defaults[name]
        return record

    def add_nodes(self, nodes):
        """
        Bulk-add nodes from a list.

        :param nodes: List (or iterator) of nodes. See
                      :func:`~RegionMatrixTable.add_node`
                      for details.
        """
        self.add_regions(nodes)

    def add_edges(self, edges):
        """
        Bulk-add edges from a list.

        :param edges: List (or iterator) of edges. See
                      :func:`~RegionMatrixTable.add_edge`
                      for details
        """
        for edge in edges:
            self.add_edge(edge, flush=False)
        self.flush(flush_nodes=False)

    def flush(self, flush_nodes=True, flush_edges=True,
              update_index=True, update_mappability=True):
        """
        Write data to file and flush buffers.

        :param flush_nodes: Flush nodes tables
        :param flush_edges: Flush edges table
        :param update_index: Update mask indices in edges table
        :param update_mappability: Save mappability info for fast access
        """
        if flush_nodes:
            self._regions.flush()

        if flush_edges:
            if len(self._edge_buffer) > 0:
                self._edges.append(self._edge_buffer)
                self._edge_buffer = []
            self._edges.flush(update_index=update_index)

        if update_mappability:
            self.meta['has_mappability_info'] = False
            self.mappable()

    def _edge_subset_rows(self, key, only_intrachromosomal=False):
        nodes_row, nodes_col = self._get_nodes_from_key(key, as_index=False)

        nodes_ix_row = None
        if nodes_row is not None:
            if isinstance(nodes_row, list):
                nodes_ix_row = [node.ix for node in nodes_row]
            else:
                nodes_ix_row = nodes_row.ix

        nodes_ix_col = None
        if nodes_col is not None:
            if isinstance(nodes_col, list):
                nodes_ix_col = [node.ix for node in nodes_col]
            else:
                nodes_ix_col = nodes_col.ix

        row_ranges = list(self._get_node_ix_ranges(nodes_ix_row))
        col_ranges = list(self._get_node_ix_ranges(nodes_ix_col))

        # fill matrix with weights
        for row_range in row_ranges:
            for col_range in col_ranges:
                for edge_row in self._edge_row_range(row_range[0], row_range[1],
                                                     col_range[0], col_range[1],
                                                     only_intrachromosomal=only_intrachromosomal):
                    yield edge_row

    def edge_subset(self, key=slice(0, None, None), lazy=False, auto_update=True,
                    only_intrachromosomal=False, **kwargs):
        """
        Get a subset of edges.

        :param key: Possible key types are:

                    Region types

                    - Node: Only the ix of this node will be used for
                      identification
                    - GenomicRegion: self-explanatory
                    - str: key is assumed to describe a genomic region
                      of the form: <chromosome>[:<start>-<end>:[<strand>]],
                      e.g.: 'chr1:1000-54232:+'

                    Node types

                    - int: node index
                    - slice: node range

                    List types

                    - list: This key type allows for a combination of all
                      of the above key types - the corresponding matrix
                      will be concatenated


                    If the key is a 2-tuple, each entry will be treated as the
                    row and column key, respectively,
                    e.g.: 'chr1:0-1000, chr4:2300-3000' will extract the Hi-C
                    map of the relevant regions between chromosomes 1 and 4.
        :param lazy: Enable lazy loading of edge attributes
        :param auto_update: Automatically update edge attributes on change
        :param only_intrachromosomal: Only return intra-chromosomal contacts
        :return: generator (:class:`~Edge`)
        """

        for edge_row in self._edge_subset_rows(key, only_intrachromosomal=only_intrachromosomal):
            yield self._row_to_edge(edge_row, lazy=lazy, auto_update=auto_update, **kwargs)

    def _get_nodes_from_key(self, key, as_index=False):
        if isinstance(key, tuple):
            nodes_ix_row = self._getitem_nodes(key[0], as_index=as_index)
            nodes_ix_col = self._getitem_nodes(key[1], as_index=as_index)
        else:
            nodes_ix_row = self._getitem_nodes(key, as_index=as_index)
            nodes_ix_col = []
            for region in self.regions():
                if as_index:
                    nodes_ix_col.append(region.ix)
                else:
                    nodes_ix_col.append(region)

        return nodes_ix_row, nodes_ix_col

    def _edge_row_range(self, source_start, source_end, sink_start, sink_end, only_intrachromosomal=False):
        condition = "(source > %d) & (source < %d) & (sink > %d) & (sink < %d)"
        condition1 = condition % (source_start-1, source_end+1, sink_start-1, sink_end+1)
        condition2 = condition % (sink_start-1, sink_end+1, source_start-1, source_end+1)

        if source_start > sink_start:
            condition1, condition2 = condition2, condition1

        regions_dict = None
        if only_intrachromosomal:
            regions_dict = self.regions_dict

        overlap = range_overlap(source_start, source_end, sink_start, sink_end)

        for edge_row in self._edges.where(condition1):
            if (only_intrachromosomal and
                    regions_dict[edge_row['source']].chromosome != regions_dict[edge_row['sink']].chromosome):
                continue
            yield edge_row

        for edge_row in self._edges.where(condition2):
            if overlap is not None:
                if (overlap[0] <= edge_row['source'] <= overlap[1]) and (overlap[0] <= edge_row['sink'] <= overlap[1]):
                    continue

            if (only_intrachromosomal and
                    regions_dict[edge_row['source']].chromosome != regions_dict[edge_row['sink']].chromosome):
                continue
            yield edge_row

    def _get_node_ix_ranges(self, nodes_ix=None):
        if not isinstance(nodes_ix, list):
            nodes_ix = [nodes_ix]

        # get range generator
        return ranges(nodes_ix)

    def _getitem_nodes(self, key, as_index=False):
        # None
        if key is None:
            key = slice(0, None, None)

        # 'chr1:1234:56789'
        if isinstance(key, string_types):
            key = GenomicRegion.from_string(key)

        # Node('chr1', 1234, 56789, ix=0)
        if isinstance(key, Node):
            if as_index:
                return key.ix
            else:
                return key

        # GenomicRegion('chr1', 1234, 56789)
        if isinstance(key, GenomicRegion):
            chromosome = key.chromosome
            start = key.start
            end = key.end

            # check defaults
            if chromosome is None:
                raise ValueError("Genomic region must provide chromosome name")
            if start is None:
                start = 0
            if end is None:
                end = max(row['end'] for row in self._regions.where("(chromosome == b'%s')" % chromosome))

            condition = "(chromosome == b'%s') & (end >= %d) & (start <= %d)" % (chromosome, start, end)
            if as_index:
                region_nodes = [row['ix'] for row in self._regions.where(condition)]
            else:
                region_nodes = [self._row_to_region(row) for row in self._regions.where(condition)]

            return region_nodes

        # 1:453
        if isinstance(key, slice):
            if as_index:
                return [row['ix'] for row in self._regions.iterrows(key.start, key.stop, key.step)]
            else:
                return [self._row_to_region(row) for row in self._regions.iterrows(key.start, key.stop, key.step)]

        # 432
        if isinstance(key, int):
            row = self._regions[key]
            if as_index:
                return row['ix']
            else:
                return self._row_to_region(row)

        # [item1, item2, item3]
        all_nodes_ix = []
        for item in key:
            nodes_ix = self._getitem_nodes(item, as_index=as_index)
            if isinstance(nodes_ix, list):
                all_nodes_ix += nodes_ix
            else:
                all_nodes_ix.append(nodes_ix)
        return all_nodes_ix

    def _row_to_node(self, row, lazy=False):
        if lazy:
            return LazyNode(row)
        return Node(chromosome=row["chromosome"], start=row["start"],
                    end=row["end"], ix=row["ix"])

    def get_node(self, key):
        """
        Get a single node by key.

        :param key: For possible key types see :func:`~RegionMatrixTable.__getitem__`
        :return: A :class:`~Node` matching key
        """
        found_nodes = self.get_nodes(key)
        if isinstance(found_nodes, list):
            if len(found_nodes) > 1:
                raise IndexError("More than one node found matching %s" % str(key))
            if len(found_nodes) == 1:
                return found_nodes[0]
        return found_nodes

    def get_nodes(self, key):
        """
        Get multiple nodes by key.

        :param key: For possible key types see :func:`~RegionMatrixTable.__getitem__`
        :return: A list of :class:`~Node` objects matching key
        """
        return self._getitem_nodes(key)

    def _row_to_edge(self, row, lazy=False, auto_update=True, **kwargs):
        if not lazy:
            source = row["source"]
            sink = row["sink"]
            d = dict()
            for field in self.field_names:
                if field != 'source' and field != 'sink':
                    value = row[field]
                    value = value.decode() if isinstance(value, bytes) else value
                    d[field] = value

            source_node_row = self._regions[source]
            source_node = self._row_to_node(source_node_row)
            sink_node_row = self._regions[sink]
            sink_node = self._row_to_node(sink_node_row)
            return Edge(source_node, sink_node, **d)
        else:
            return LazyEdge(row, self._regions, auto_update=auto_update)

    def get_edge(self, ix, lazy=False):
        """
        Get an edge from this object's edge list.

        :param ix: integer
        :param lazy: Use lazy loading of object attributes. Do not
                     use lazy objects outside of loop iterations!
        :return:
        """
        return self._row_to_edge(self._edges[ix], lazy=lazy)

    def nodes(self):
        """
        Iterator over this object's nodes/regions.

        See :func:`~RegionsTable.regions` for details.
        :return: Iterator over :class:`~GenomicRegions`
        """
        return self._nodes_iter()

    def _nodes_iter(self):
        return self.regions

    @property
    def edges(self):
        """
        Iterate over :class:`~Edge` objects.

        :return: Iterator over :class:`~Edge`
        """
        return self._edges_iter()

    def _edges_iter(self):
        return RegionPairs.EdgeIter(self)

    def _is_sorted(self, sortby):
        column = getattr(self._edges.cols, sortby)
        if (column.index is None or
                not column.index.is_csi):
            return False
        return True

    def create_cs_index(self, field):
        column = getattr(self._edges.cols, field)

        if not self._is_sorted(field):
            try:
                logger.info("Sorting {}...".format(field))
                if not column.is_indexed:
                    column.create_csindex()
                elif not column.index.is_csi:
                    column.reindex()
            except t.exceptions.FileModeError:
                raise RuntimeError("This object is not sorted by requested column! "
                                   "Cannot sort manually, because file is in read-only mode.")

    def edges_sorted(self, sortby, reverse=False, *args, **kwargs):
        """
        Iterate over edges sorted by a specific column.

        :param sortby: Name of column to sort over
        :param reverse: Iterate in reverse order
        :return: EdgeIter iterator
        """
        # ensure sorting on qname_ix column
        self.create_cs_index(sortby)

        if reverse:
            step = -1
        else:
            step = None
        edge_iter = RegionMatrixTable.EdgeIter(self, _iter=self._edges.itersorted(sortby, step=step))
        return edge_iter(*args, **kwargs)

    def __iter__(self):
        return self.edges

    def __len__(self):
        return len(self._edges)

    def mappable(self, sub_region=None):
        """
        Get the mappability vector of this matrix.
        """
        mappable = np.zeros(len(self.regions), dtype=bool)
        # check if mappability info already exists
        if 'has_mappability_info' in self.meta and self.meta['has_mappability_info']:
            logger.debug("Retrieving precalculated mappability...")
            for i, region in enumerate(self.regions(lazy=True)):
                mappable[i] = region._mask_ix == 0
        else:
            # prepare marginals dict
            logger.debug("Calculating mappability...")

            with RareUpdateProgressBar(max_value=len(self.edges), silent=config.hide_progressbars) as pb:
                for i, edge in enumerate(self.edges(lazy=True)):
                    mappable[edge.source] = True
                    mappable[edge.sink] = True
                    pb.update(i)

            try:
                for i, row in enumerate(self._regions):
                    if not mappable[i]:
                        row['_mask_ix'] = 1
                    else:
                        row['_mask_ix'] = 0
                    row.update()
                self._regions.flush()
                self.meta['has_mappability_info'] = True
            except (IOError, OSError, t.FileModeError, KeyError):
                logger.debug("Cannot write mappability info to read-only file")

        if sub_region is None:
            return mappable
        else:
            return mappable[self.region_bins(sub_region)]


class AccessOptimisedRegionPairs(RegionPairs):
    """
    Extends :class:`~RegionPairs` with a backend that partitions edges into chromosomes.

    This partitioning should greatly speed up edge and matrix queries for large Hi-C data sets,
    such as high-resolution (<=10kb) human Hi-C maps.

    Iterating over sorted edges performance is somewhat reduced due to the fact that we have to
    integrate tens to hundreds of tables in the sorting.
    """

    _classid = 'ACCESSOPTIMISEDREGIONPAIRS'

    class EdgeIter(RegionPairs.EdgeIter):
        """
        Class providing iterator functionality to a :class:`~RegionPairs` object.
        """
        def __init__(self, this, _iter=None):
            RegionPairs.EdgeIter.__init__(self, this, _iter=_iter)
            self.iter = _iter
            self.interchromosomal = True
            self.intrachromosomal = True

        def __getitem__(self, item):
            if isinstance(item, int):
                return self.this.get_edge(item, intrachromosomal=self.intrachromosomal,
                                          interchromosomal=self.interchromosomal)
            elif isinstance(item, slice):
                edges = []
                for row in self.get_row_range(item):
                    edge = self.this._row_to_edge(row, *self.row_conversion_args, **self.row_conversion_kwargs)
                    edges.append(edge)
                return edges

        def get_row_range(self, item):
            start = 0 if item.start is None else item.start
            stop = len(self.this.edges) if item.stop is None else item.stop
            step = 1 if item.step is None else item.step
            if step != 1:
                raise ValueError("Step sizes != 1 not currently supported in slices. %s" % str(step))

            l = 0
            for edge_table in self.this._edge_table_iter(intrachromosomal=self.intrachromosomal,
                                                         interchromosomal=self.interchromosomal):
                # not yet in range
                if start >= l + len(edge_table):
                    l += len(edge_table)
                    continue
                # over range - can stop here
                if stop < l:
                    break

                # in range, get edges
                r = (max(0, start - l), min(len(edge_table), stop - l))
                res = edge_table[r[0]:r[1]]
                for row in res:
                    yield row
                l += len(edge_table)

        def __iter__(self):
            return self

        def __call__(self, *args, **kwargs):
            if 'only_intrachromosomal' in kwargs:
                self.interchromosomal = False
                self.intrachromosomal = True
                del kwargs['only_intrachromosomal']
            if 'intrachromosomal' in kwargs:
                self.intrachromosomal = kwargs['intrachromosomal']
                del kwargs['intrachromosomal']
            if 'interchromosomal' in kwargs:
                self.interchromosomal = kwargs['interchromosomal']
                del kwargs['interchromosomal']
            self.row_conversion_args = args
            self.row_conversion_kwargs = kwargs

            key = kwargs.get('key', None)
            if key is not None:
                self.iter = self.this._edge_subset_rows(key=key,
                                                        only_intrachromosomal=self.only_intrachromosomal)
            return iter(self)

        def __next__(self):
            if self.iter is None:
                self.iter = self.this._edge_row_iter(intrachromosomal=self.intrachromosomal,
                                                     interchromosomal=self.interchromosomal)
            row = next(self.iter)
            return self.this._row_to_edge(row, *self.row_conversion_args, **self.row_conversion_kwargs)

        def __len__(self):
            return len(self.this)

    def __init__(self, file_name=None, mode='a', tmpdir=None, additional_fields=None,
                 _table_name_nodes='nodes', _table_name_edges='edges', _edge_buffer_size=1000000):
        # private variables
        self._max_node_ix = -1

        if file_name is not None:
            file_name = os.path.expanduser(file_name)

        # initialize inherited objects
        RegionsTable.__init__(self, file_name=file_name, _table_name_regions=_table_name_nodes,
                              mode=mode, tmpdir=tmpdir)
        Maskable.__init__(self, self.file)

        # create edge table
        self._field_dict = None
        self.field_names = None
        self._edge_table_dict = dict()
        self._edge_field_defaults = dict()
        self._source_field_ix = 0
        self._sink_field_ix = 1

        # existing one
        if _table_name_edges in self.file.root:
            self._edges = self.file.get_node('/', _table_name_edges)
            for edge_table in self._edges:
                if self._field_dict is None:
                    self._field_dict = edge_table.coldescrs
                    self._update_field_names(edge_table=edge_table)
                source_partition = edge_table.attrs['source_partition']
                sink_partition = edge_table.attrs['sink_partition']
                self._edge_table_dict[(source_partition, sink_partition)] = edge_table
        else:
            # create edge table definition
            self._field_dict = self._get_field_dict(additional_fields=additional_fields)

            self._edges = self.file.create_group('/', _table_name_edges)
            # will always have 0-0 partition
            edge_table = self._create_edge_table(0, 0)
            self._update_field_names(edge_table=edge_table)

        # update partitions
        self._update_partitions()

        self._edge_buffer = defaultdict(list)
        self._edge_buffer_size = _edge_buffer_size

    def _flush_table_edge_buffer(self):
        for (source_partition, sink_partition), records in self._edge_buffer.items():
            if not (source_partition, sink_partition) in self._edge_table_dict:
                self._create_edge_table(source_partition, sink_partition)
            table = self._edge_table_dict[(source_partition, sink_partition)]
            table.append(records)
        self._edge_buffer = defaultdict(list)

    def _update_field_names(self, edge_table=None):
        """
        Set internal object variables related to edge table field names.
        """
        if edge_table is None:
            for et in self._edges:
                edge_table = et
                break

        if edge_table is None:
            return

        # update field names
        self._source_field_ix = 0
        self._sink_field_ix = 0
        self.field_names = []
        self._field_names_dict = dict()
        for i, name in enumerate(edge_table.colnames):
            if not name.startswith("_"):
                self.field_names.append(name)
            if name == 'source':
                self._source_field_ix = i
            if name == 'sink':
                self._sink_field_ix = i
            self._field_names_dict[name] = i
            self._edge_field_defaults[name] = edge_table.coldescrs[name].dflt

    def _update_partitions(self):
        """
        Update the list of partition break points (split by chromosome)
        """
        self.partitions = []
        previous_chromosome = None
        for i, region in enumerate(self.regions(lazy=True)):
            if region.chromosome != previous_chromosome and previous_chromosome is not None:
                self.partitions.append(i)
            previous_chromosome = region.chromosome

    def flush(self, flush_nodes=True, flush_edges=True,
              update_index=True, update_mappability=True,
              silent=config.hide_progressbars):
        """
        Write data to file and flush buffers.

        :param flush_nodes: Flush nodes tables
        :param flush_edges: Flush edges table
        :param update_index: Update mask indices in edges table
        :param update_mappability: Save mappability info for fast access
        :param silent: do not print flush progress
        """
        if flush_nodes:
            self._regions.flush()
            # update partitions
            self._update_partitions()

        if flush_edges:
            if len(self._edge_buffer) > 0:
                self._flush_table_edge_buffer()

            if update_index:
                logger.info("Updating mask indices...")

            with RareUpdateProgressBar(max_value=sum(1 for _ in self._edges), silent=silent) as pb:
                for i, edge_table in enumerate(self._edges):
                    edge_table.flush(update_index=update_index, log_progress=False)
                    pb.update(i)

        if update_mappability:
            self.meta['has_mappability_info'] = False
            self.mappable()

    def add_regions(self, regions):
        super(AccessOptimisedRegionPairs, self).add_regions(regions)
        self._update_partitions()

    def _get_field_dict(self, additional_fields=None):
        """
        Generate a dictionary of PyTables fields to create edge table.

        Save fields dict in object variable - we need to generate a lot of tables.
        """
        if self._field_dict is not None:
            return self._field_dict
        return RegionPairs._get_field_dict(self, additional_fields=additional_fields)

    def _get_partition_ix(self, region_ix):
        """
        Bisect the partition table to get the partition index for a region index.
        """
        return bisect_right(self.partitions, region_ix)

    def _create_edge_table(self, source_partition, sink_partition):
        """
        Create and register an edge table for a partition combination.
        """
        if (source_partition, sink_partition) in self._edge_table_dict:
            return self._edge_table_dict[(source_partition, sink_partition)]

        edge_table = MaskedTable(self._edges,
                                 'chrpair_' + str(source_partition) + '_' + str(sink_partition),
                                 self._field_dict)
        edge_table.attrs['source_partition'] = source_partition
        edge_table.attrs['sink_partition'] = sink_partition

        # index
        create_col_index(edge_table.cols.source)
        create_col_index(edge_table.cols.sink)

        self._edge_table_dict[(source_partition, sink_partition)] = edge_table
        return edge_table

    def disable_indexes(self):
        for edge_table in self._edge_table_iter():
            edge_table.cols.source.remove_index()
            edge_table.cols.sink.remove_index()
            edge_table.disable_mask_index()

    def enable_indexes(self):
        for edge_table in self._edge_table_iter():
            create_col_index(edge_table.cols.source)
            create_col_index(edge_table.cols.sink)
            edge_table.enable_mask_index()

    def _get_edge_table_tuple(self, source, sink):
        if source > sink:
            source, sink = sink, source

        source_partition = self._get_partition_ix(source)
        sink_partition = self._get_partition_ix(sink)

        return source_partition, sink_partition

    def _get_edge_table(self, source, sink):
        """
        Return an edge table for this particular region index combination.
        """
        source_partition, sink_partition = self._get_edge_table_tuple(source, sink)

        if not (source_partition, sink_partition) in self._edge_table_dict:
            self._create_edge_table(source_partition, sink_partition)
        return self._edge_table_dict[(source_partition, sink_partition)]

    def _edge_table_iter(self, intrachromosomal=True, interchromosomal=True):
        """
        Iterate over internal edge tables.

        :param intrachromosomal: If true, include intra-chromosomal edge tables
        :param interchromosomal: If true, include inter-chromosomal edge tables
        :return: Edge table iterator
        """
        # intra-chromosomal
        if intrachromosomal:
            for i in range(len(self.partitions) + 1):
                if (i, i) in self._edge_table_dict:
                    yield self._edge_table_dict[(i, i)]

        # inter-chromosomal
        if interchromosomal:
            for i in range(len(self.partitions) + 1):
                for j in range(i + 1, len(self.partitions) + 1):
                    if (i, j) in self._edge_table_dict:
                        yield self._edge_table_dict[(i, j)]

    def _edge_from_list(self, edge):
        source, sink = edge[self._source_field_ix], edge[self._sink_field_ix]
        self._get_edge_table(source, sink)

        return RegionPairs._edge_from_list(self, edge)

    def _add_edge(self, edge, row, replace=False):
        """
        Add an edge to an internal edge table.
        """
        source, sink = edge.source, edge.sink
        if source > sink:
            source, sink = sink, source

        if row is None:
            record = [None] * len(self._field_names_dict)
            for name, ix in self._field_names_dict.items():
                try:
                    record[ix] = getattr(edge, name)
                except AttributeError:
                    record[ix] = self._edge_field_defaults[name]
            record[self._field_names_dict['source']] = source
            record[self._field_names_dict['sink']] = sink

            source_partition, sink_partition = self._get_edge_table_tuple(source, sink)

            self._edge_buffer[(source_partition, sink_partition)].append(tuple(record))
            if sum(len(records) for records in self._edge_buffer.values()) > self._edge_buffer_size:
                self._flush_table_edge_buffer()
        else:
            update = True
            if row is None:
                update = False
                table = self._get_edge_table(source, sink)
                row = table.row
            row['source'] = source
            row['sink'] = sink
            for name in self.field_names:
                if not name == 'source' and not name == 'sink':
                    try:
                        value = getattr(edge, name)
                        if replace or not update:
                            row[name] = value
                        else:
                            row[name] += value
                    except AttributeError:
                        pass
            if update:
                row.update()
            else:
                row.append()

    def _add_edge_from_tuple(self, edge):
        source = edge[self._source_field_ix]
        sink = edge[self._sink_field_ix]
        if source > sink:
            source, sink = sink, source

        source_partition, sink_partition = self._get_edge_table_tuple(source, sink)

        self._edge_buffer[(source_partition, sink_partition)].append(tuple(edge))
        if sum(len(records) for records in self._edge_buffer.values()) > self._edge_buffer_size:
            self._flush_table_edge_buffer()

    def get_edge(self, item, intrachromosomal=True, interchromosomal=True,
                 *row_conversion_args, **row_conversion_kwargs):
        """
        Get an edge by index.

        :param intrachromosomal: If true, include intra-chromosomal edge tables in index count
        :param interchromosomal: If true, include inter-chromosomal edge tables in index count
        :param row_conversion_args: Arguments passed to :func:`RegionPairs._row_to_edge`
        :param row_conversion_args: Keyword arguments passed to :func:`RegionPairs._row_to_edge`
        :return: :class:`~Edge`
        """
        if item < 0:
            item += len(self)

        l = 0
        for edge_table in self._edge_table_iter(intrachromosomal=intrachromosomal,
                                                interchromosomal=interchromosomal):
            if l <= item < l + len(edge_table):
                res = edge_table[item - l]
                return self._row_to_edge(res, *row_conversion_args, **row_conversion_kwargs)
            l += len(edge_table)
        raise IndexError("index out of range (%d)" % item)

    @property
    def edges(self):
        """
        Iterate over :class:`~Edge` objects.

        :return: Iterator over :class:`~Edge`
        """
        return self._edges_iter()

    def _edges_iter(self):
        return AccessOptimisedRegionPairs.EdgeIter(self)

    def _edge_row_iter(self, intrachromosomal=True, interchromosomal=True, excluded_filters=()):
        """
        Yield rows in edge tables, ordered by partition.
        """
        for edge_table in self._edge_table_iter(intrachromosomal=intrachromosomal, interchromosomal=interchromosomal):
            excluded_masks = self.get_binary_mask_from_masks(excluded_filters)
            for row in edge_table.iterrows(excluded_masks=excluded_masks):
                yield row

    def _partition_ix_range(self, start, stop):
        """
        Get a range of partitions with start and stop indices per partition from global start and stop region indices.

        :param start: Region start index
        :param stop: Region stop index
        :return: tuple, where first element is
        """
        start_partition = self._get_partition_ix(start)
        stop_partition = self._get_partition_ix(stop)

        def _is_start_of_partition(start_ix, partition_ix):
            if partition_ix == 0:
                if start_ix == 0:
                    return True
            else:
                if start_ix == self.partitions[partition_ix - 1]:
                    return True
            return False

        def _is_end_of_partition(stop_ix, partition_ix):
            if partition_ix == len(self.partitions):
                if stop_ix == len(self.regions)-1:
                    return True
            else:
                if stop_ix == self.partitions[partition_ix]-1:
                    return True
            return False

        if start_partition == stop_partition:
            complete = _is_start_of_partition(start, start_partition) and _is_end_of_partition(stop, stop_partition)
            return [(start, stop, start_partition, complete)]

        partition_ranges = []
        start_range_complete = _is_start_of_partition(start, start_partition)
        start_range = (start, self.partitions[start_partition] - 1, start_partition, start_range_complete)
        partition_ranges.append(start_range)

        for i in range(start_partition + 1, stop_partition):
            partition_ranges.append((self.partitions[i-1], self.partitions[i]-1, i, True))

        stop_range_complete = _is_end_of_partition(stop, stop_partition)
        stop_range = (self.partitions[stop_partition - 1], stop, stop_partition, stop_range_complete)
        partition_ranges.append(stop_range)
        return partition_ranges

    def _edge_row_range(self, source_start, source_end, sink_start, sink_end, only_intrachromosomal=False):
        """
        Iterate over a range of rows in this object's edge tables.

        Rows are selected based on region indices of interacting regions.
        """
        source_partition_ranges = self._partition_ix_range(source_start, source_end)
        sink_partition_ranges = self._partition_ix_range(sink_start, sink_end)

        covered = set()
        for source_partition_range in source_partition_ranges:
            source_start, source_end, source_partition, source_complete = source_partition_range

            for sink_partition_range in sink_partition_ranges:
                sink_start, sink_stop, sink_partition, sink_complete = sink_partition_range

                if only_intrachromosomal and source_partition != sink_partition:
                    continue

                if sink_partition < source_partition:
                    key = (sink_partition, source_partition)
                else:
                    key = (source_partition, sink_partition)

                if key in covered:
                    continue
                else:
                    covered.add(key)

                if key in self._edge_table_dict:
                    table = self._edge_table_dict[key]

                    # entire partition is requested, no need for where query
                    if source_complete and sink_complete:
                        for edge_row in table:
                            yield edge_row
                    else:
                        condition = "(source > %d) & (source < %d) & (sink > %d) & (sink < %d)"
                        condition1 = condition % (source_start - 1, source_end + 1, sink_start - 1, sink_end + 1)
                        condition2 = condition % (sink_start - 1, sink_end + 1, source_start - 1, source_end + 1)

                        if source_start > sink_start:
                            condition1, condition2 = condition2, condition1

                        overlap = range_overlap(source_start, source_end, sink_start, sink_end)

                        for edge_row in table.where(condition1):
                            yield edge_row

                        for edge_row in table.where(condition2):
                            if overlap is not None:
                                if (overlap[0] <= edge_row['source'] <= overlap[1]) and (overlap[0] <= edge_row['sink'] <= overlap[1]):
                                    continue

                            yield edge_row

    def _is_sorted(self, sortby):
        """
        For each edge table, check if it is sorted.
        """
        for edge_table in self._edge_table_iter():
            column = getattr(edge_table.cols, sortby)
            if (column.index is None or
                    not column.index.is_csi):
                return False
        return True

    def _edge_row_iter_sorted(self, sortby, step=None, intrachromosomal=True, interchromosomal=True):
        """
        Yield rows in edge tables, ordered by partition.
        """
        table_iterators = []
        for edge_table in self._edge_table_iter(intrachromosomal=intrachromosomal, interchromosomal=interchromosomal):
            table_iterators.append(iter(edge_table.itersorted(sortby, step=step)))

        rows = []
        for i, table_iterator in enumerate(table_iterators):
            try:
                row = next(table_iterator)
                rows.append(row)
            except StopIteration:
                del table_iterators[i]

        while len(table_iterators) > 0:
            # find current minimum or maximum
            current = None
            current_ix = None
            if step is None or step >= 0:
                for i, row in enumerate(rows):
                    if current is None or row[sortby] < current:
                        current = row[sortby]
                        current_ix = i
            else:
                for i, row in enumerate(rows):
                    if current is None or row[sortby] > current:
                        current = row[sortby]
                        current_ix = i

            yield rows[current_ix]

            try:
                rows[current_ix] = next(table_iterators[current_ix])
            except StopIteration:
                del table_iterators[current_ix]
                del rows[current_ix]

    def create_cs_index(self, field):
        for edge_table in self._edge_table_iter():
            # ensure sorting on sortby column
            column = getattr(edge_table.cols, field)

            if not self._is_sorted(field):
                try:
                    logger.info("Sorting %s..." % field)
                    if not column.is_indexed:
                        column.create_csindex()
                    elif not column.index.is_csi:
                        column.reindex()
                except t.exceptions.FileModeError:
                    raise RuntimeError("This object is not sorted by requested column! "
                                       "Cannot sort manually, because file is in read-only mode.")

    def edges_sorted(self, sortby, reverse=False, *args, **kwargs):
        """
        Iterate over edges sorted by *sortby*.
        """
        self.create_cs_index(sortby)
        if reverse:
            step = -1
        else:
            step = None
        edge_iter = AccessOptimisedRegionPairs.EdgeIter(self, _iter=self._edge_row_iter_sorted(sortby, step=step))
        return edge_iter(*args, **kwargs)

    def __len__(self):
        l = 0
        for edge_table in self._edge_table_iter():
            l += len(edge_table)
        return l

    def __iter__(self):
        return self.edges


class RegionMatrixTable(RegionPairs):
    """
    Class for working with matrix-based data.

    Generally, a RegionMatrix object has two components:

    - Nodes or regions: (Non-overlapping) genomic regions
      obtained by splitting the genome into distinct pieces.
      See also :class:`~GenomicRegion` and :class:`~RegionsTable`

    - Edges or contacts: Pairs of genomic regions with optionally
      associated weight or contact strength. See also
      :class:`~Edge`

    This is a memory-efficient implementation of a matrix data
    container. Internally, this is achieved by saving entries
    of the matrix in sparse notation, i.e. in a list of
    non-zero contacts.

    Its bracket-notation access behaves like a numpy
    array and handles data retrieval and assignment in matrix-
    fashion, e.g. m[1:3] would return rows 1 and 2 of
    the matrix m (0-based index). However, the bracket
    notation can also handle :class:`~GenomicRegion` descriptor
    strings, i.e. m['chr1','chr5'] will extract the inter-
    chromosomal matrix between chromosomes 1 and 5 only.

    Examples:

    .. code:: python

        m = RegionMatrix(file_name="/path/to/save/file")

        # load genomic regions
        genome = Genome.from_folder("/path/to/fasta/folder")
        regions = genome.get_regions("HindIII")
        m.add_regions(regions)

        # load edges
        edges = []
        edges.append(Edge(source=10, sink=23, weight=3)
        edges.append(Edge(source=8, sink=9, weight=57)
        # ...
        m.add_edges(edges)
    """

    _classid = 'REGIONMATRIXTABLE'

    def __init__(self, file_name=None, mode='a', additional_fields=None, tmpdir=None,
                 default_field=None, default_value=0.0,
                 _table_name_nodes='nodes', _table_name_edges='edges'):

        """
        Initialize a :class:`~RegionMatrixTable` object.

        :param file_name: Path to a save file
        :param mode: File mode to open underlying file
        :param _table_name_nodes: (Internal) name of the HDF5 node for regions
        :param _table_name_edges: (Internal) name of the HDF5 node for edges
        """

        # private variables
        self.default_field = default_field
        self.default_value = default_value
        self._expected_values = None
        RegionPairs.__init__(self, file_name=file_name, mode=mode, additional_fields=additional_fields, tmpdir=tmpdir,
                             _table_name_nodes=_table_name_nodes, _table_name_edges=_table_name_edges)

        if default_field is None:
            for field_name in self._edges.colnames:
                if not field_name.startswith("_") and field_name != "source" and field_name != "sink":
                    self.default_field = field_name
                    break

    def _flush_edge_buffer(self, e_buffer, replace=False,
                           update_index=True, update_mappability=True,
                           clean_zero=True, default_column=None):
        if default_column is None:
            default_column = self.default_field
        # update current rows
        for row in self._edges:
            key = (row["source"], row["sink"])

            if key in e_buffer:
                value = e_buffer[key]
                # it is a weight
                try:
                    if replace:
                        row[default_column] = float(value)
                    else:
                        row[default_column] += float(value)
                    row.update()
                except TypeError:
                    self.add_edge(value, check_nodes_exist=False, flush=False, replace=replace, row=row)
                del e_buffer[key]
        self.flush(update_index=False, update_mappability=False)

        # flush remaining buffer
        for source, sink in e_buffer:
            key = (source, sink)
            value = e_buffer[key]
            try:
                v = float(value)
                if v == 0:
                    continue
                new_edge = self._edge_from_dict({'source': source, 'sink': sink, default_column: v})
                self.add_edge(new_edge, check_nodes_exist=False, flush=False)
            except TypeError:
                self.add_edge(value, check_nodes_exist=False, flush=False)

        self.flush(update_index=update_index, update_mappability=update_mappability)
        if clean_zero:
            self._remove_zero_edges(update_index=update_index, update_mappability=update_mappability,
                                    weight_column=default_column)

    def __getitem__(self, key):
        return self.as_matrix(key)

    def as_matrix(self, key=slice(0, None, None), values_from=None, mask_missing=False,
                  default_value=None, impute_missing=False, _mappable=None):
        """
        Get a chunk of the matrix.

        :param key: Possible key types are:

                    Region types

                    - Node: Only the ix of this node will be used for
                      identification
                    - GenomicRegion: self-explanatory
                    - str: key is assumed to describe a genomic region
                      of the form: <chromosome>[:<start>-<end>:[<strand>]],
                      e.g.: 'chr1:1000-54232:+'

                    Node types

                    - int: node index
                    - slice: node range

                    List types

                    - list: This key type allows for a combination of all
                      of the above key types - the corresponding matrix
                      will be concatenated


                    If the key is a 2-tuple, each entry will be treated as the
                    row and column key, respectively,
                    e.g.: 'chr1:0-1000, chr4:2300-3000' will extract the Hi-C
                    map of the relevant regions between chromosomes 1 and 4.
        :param values_from: Determines which column will be used to populate
                            the matrix. Default is 'weight'.
        :param mask_missing: if True, will mask missing/unmappable contacts
        :param impute_missing: if True, will average missing contacts
        :return: :class:`RegionMatrix`
        """
        if values_from is None:
            values_from = self.default_field

        nodes_row, nodes_col = self._get_nodes_from_key(key, as_index=False)

        nodes_ix_row = None
        if nodes_row is not None:
            if isinstance(nodes_row, list):
                nodes_ix_row = [node.ix for node in nodes_row]
            else:
                nodes_ix_row = nodes_row.ix

        nodes_ix_col = None
        if nodes_col is not None:
            if isinstance(nodes_col, list):
                nodes_ix_col = [node.ix for node in nodes_col]
            else:
                nodes_ix_col = nodes_col.ix

        row_ranges = list(self._get_node_ix_ranges(nodes_ix_row))
        col_ranges = list(self._get_node_ix_ranges(nodes_ix_col))

        m = self._get_matrix(row_ranges, col_ranges, weight_column=values_from,
                             default_value=default_value)

        # select the correct output format
        # empty result: matrix
        rm = None
        if m.shape[0] == 0 and m.shape[1] == 0:
            rm = RegionMatrix(m, col_regions=[], row_regions=[])
        # both selectors are lists: matrix
        elif isinstance(nodes_ix_row, list) and isinstance(nodes_ix_col, list):
            rm = RegionMatrix(m, col_regions=nodes_col, row_regions=nodes_row)
        # row selector is list: vector
        elif isinstance(nodes_ix_row, list):
            rm = RegionMatrix(m[:, 0], col_regions=[nodes_col], row_regions=nodes_row)
        # column selector is list: vector
        elif isinstance(nodes_ix_col, list):
            rm = RegionMatrix(m[0, :], col_regions=nodes_col, row_regions=[nodes_row])

        if rm is not None:
            if mask_missing or impute_missing:
                if _mappable is None:
                    mappable = self.mappable()
                else:
                    mappable = _mappable
                mask = np.zeros(m.shape, dtype=bool)
                current_row = 0
                for row_range in row_ranges:
                    for i in range(row_range[0], row_range[1]+1):
                        if not mappable[i]:
                            mask[current_row] = True
                        current_row += 1

                current_col = 0
                for col_range in col_ranges:
                    for i in range(col_range[0], col_range[1]+1):
                        if not mappable[i]:
                            mask[:, current_col] = True
                        current_col += 1
                masked_rm = np.ma.MaskedArray(rm, mask=mask)

                if impute_missing:
                    return self._impute_missing_contacts(masked_rm)
                return masked_rm
            else:
                return rm

        # both must be indexes
        return m[0, 0]

    def _impute_missing_contacts(self, hic_matrix=None, _expected_contacts=None):
        """
        Impute missing contacts in a Hi-C matrix.

        :param hic_matrix: a :class:`~HicMatrix` object
        :param _expected_contacts: An ExpectedContacts object for this matrix
        :return: the input matrix with imputed values
        """
        if not hasattr(hic_matrix, "mask"):
            raise ValueError("hic_matrix must be a numpy masked array!")

        # here to avoid circular dependency
        from fanc.architecture.hic_architecture import ExpectedContacts

        if _expected_contacts is not None:
            ex = _expected_contacts
            close_ex = False
        else:
            ex = ExpectedContacts(self, smooth=True)
            close_ex = True

        intra_expected = ex.intra_expected()
        inter_expected = ex.inter_expected()

        for i in range(hic_matrix.shape[0]):
            row_region = hic_matrix.row_regions[i]
            for j in range(hic_matrix.shape[1]):
                col_region = hic_matrix.col_regions[j]

                if hic_matrix.mask[i, j]:
                    if row_region.chromosome == col_region.chromosome:
                        d = abs(row_region.ix-col_region.ix)
                        hic_matrix[i, j] = intra_expected[d]
                    else:
                        hic_matrix[i, j] = inter_expected

        if close_ex:
            ex.close()

        return hic_matrix

    def _get_matrix(self, row_ranges, col_ranges, weight_column=None, default_value=None):
        if default_value is None:
            default_value = self.default_value

        if weight_column is None:
            weight_column = self.default_field

        n_rows = 0
        for row_range in row_ranges:
            n_rows += row_range[1]-row_range[0]+1

        n_cols = 0
        for col_range in col_ranges:
            n_cols += col_range[1]-col_range[0]+1

        # create empty matrix
        m = np.full((n_rows, n_cols), default_value)

        # fill matrix with weights
        row_offset = 0
        for row_range in row_ranges:
            n_rows_sub = row_range[1] - row_range[0] + 1
            col_offset = 0
            for col_range in col_ranges:
                n_cols_sub = col_range[1] - col_range[0] + 1

                for edge_row in self._edge_row_range(row_range[0], row_range[1], col_range[0], col_range[1]):
                    source = edge_row['source']
                    sink = edge_row['sink']
                    weight = edge_row[weight_column]

                    ir = source - row_range[0] + row_offset
                    jr = sink - col_range[0] + col_offset
                    if 0 <= ir < m.shape[0] and 0 <= jr < m.shape[1]:
                        m[ir, jr] = weight

                    ir = sink - row_range[0] + row_offset
                    jr = source - col_range[0] + col_offset
                    if 0 <= ir < m.shape[0] and 0 <= jr < m.shape[1]:
                        m[ir, jr] = weight

                col_offset += n_cols_sub
            row_offset += n_rows_sub

        return m

    def as_data_frame(self, key, weight_column=None):
        """
        Get a pandas data frame by key.

        For key types see :func:`~RegionMatrixTable.__getitem__`.

        :param key: For key types see :func:`~RegionMatrixTable.__getitem__`.
        :param weight_column: Determines which column populates the DF
        :return: Pandas data frame, row and column labels are
                 corresponding node start positions
        """
        if weight_column is None:
            weight_column = self.default_field

        nodes_ix_row, nodes_ix_col = self._get_nodes_from_key(key, as_index=True)
        nodes_row, nodes_col = self._get_nodes_from_key(key, as_index=False)

        row_ranges = list(self._get_node_ix_ranges(nodes_ix_row))
        col_ranges = list(self._get_node_ix_ranges(nodes_ix_col))

        m = self._get_matrix(row_ranges, col_ranges, weight_column=weight_column)
        labels_row = []
        for node in nodes_row:
            labels_row.append(node.start)
        labels_col = []
        for node in nodes_col:
            labels_col.append(node.start)
        df = p.DataFrame(m, index=labels_row, columns=labels_col)

        return df

    def __setitem__(self, key, item):
        self.set_matrix(key, item, clean_zero=True)

    def set_matrix(self, key, item, clean_zero=False, values_to=None):
        """
        Set a chunk of the matrix.

        :param key: Possible key types are:

                    Region types

                    - Node: Only the ix of this node will be used for
                      identification
                    - GenomicRegion: self-explanatory
                    - str: key is assumed to describe a genomic region
                      of the form: <chromosome>[:<start>-<end>:[<strand>]],
                      e.g.: 'chr1:1000-54232:+'

                    Node types

                    - int: node index
                    - slice: node range

                    List types

                    - list: This key type allows for a combination of all
                      of the above key types - the corresponding matrix
                      will be concatenated

                    If the key is a 2-tuple, each entry will be treated as the
                    row and column key, respectively,
                    e.g.: 'chr1:0-1000, chr4:2300-3000' will set the entries
                    of the relevant regions between chromosomes 1 and 4.
        :param item: matrix to replace existing values
        :param clean_zero: Remove edges where 'values_to' colum is zero
        :param values_to: Determines which column is replaced by the provided
                          item. Default: weight
        """
        if values_to is None:
            values_to = self.default_field

        nodes_ix_row, nodes_ix_col = self._get_nodes_from_key(key, as_index=True)
        self._set_matrix(item, nodes_ix_row, nodes_ix_col, clean_zero=clean_zero, weight_column=values_to)

    def _set_matrix(self, item, nodes_ix_row=None, nodes_ix_col=None,
                    clean_zero=False, weight_column=None):
        if weight_column is None:
            weight_column = self.default_field

        replacement_edges = {}

        # create new edges with updated weights
        # select the correct format:
        def swap(old_source, old_sink):
            if old_source > old_sink:
                return old_sink, old_source
            return old_source, old_sink

        # both selectors are lists: matrix
        if isinstance(nodes_ix_row, list) and isinstance(nodes_ix_col, list):
            n_rows = len(nodes_ix_row)
            n_cols = len(nodes_ix_col)
            # check that we have a matrix with the correct dimensions
            if not isinstance(item, np.ndarray) or not np.array_equal(item.shape, [n_rows, n_cols]):
                raise ValueError("Item is not a numpy array with shape (%d,%d)!" % (n_rows, n_cols))

            for i in range(0, n_rows):
                for j in range(0, n_cols):
                    source = nodes_ix_row[i]
                    sink = nodes_ix_col[j]
                    source, sink = swap(source, sink)
                    weight = item[i, j]
                    key = (source, sink)
                    if key not in replacement_edges:
                        replacement_edges[key] = weight

        # row selector is list: vector
        elif isinstance(nodes_ix_row, list):
            n_rows = len(nodes_ix_row)
            if not isinstance(item, np.ndarray) or not np.array_equal(item.shape, [n_rows]):
                raise ValueError("Item is not a numpy vector of length %d!" % n_rows)

            for i, my_sink in enumerate(nodes_ix_row):
                source = nodes_ix_col
                source, sink = swap(source, my_sink)
                weight = item[i]
                key = (source, sink)
                if key not in replacement_edges:
                    replacement_edges[key] = weight

        # column selector is list: vector
        elif isinstance(nodes_ix_col, list):
            n_cols = len(nodes_ix_col)
            if not isinstance(item, np.ndarray) or not np.array_equal(item.shape, [n_cols]):
                raise ValueError("Item is not a numpy vector of length %d!" % n_cols)

            for i, my_source in enumerate(nodes_ix_col):
                sink = nodes_ix_row
                source, sink = swap(my_source, sink)
                weight = item[i]
                key = (source, sink)
                if key not in replacement_edges:
                    replacement_edges[key] = weight

        # both must be indexes
        else:
            weight = item
            source, sink = swap(nodes_ix_row, nodes_ix_col)
            key = (source, sink)
            if key not in replacement_edges:
                replacement_edges[key] = weight

        self._flush_edge_buffer(replacement_edges, replace=True,
                                clean_zero=clean_zero, default_column=weight_column)

    def _update_edge_weight(self, source, sink, weight, add=False, flush=True, weight_column=None):
        if weight_column is None:
            weight_column = self.default_field

        if source > sink:
            source, sink = sink, source

        value_set = False
        for row in self._edges.where("(source == %d) & (sink == %d)" % (source, sink)):
            original = 0
            if add:
                original = row[weight_column]
            row[weight_column] = weight + original
            row.update()
            value_set = True
            if flush:
                self.flush()
        if not value_set:
            self.add_edge(Edge(source=source, sink=sink, weight=weight), flush=flush)

    def _remove_zero_edges(self, flush=True, update_index=True,
                           update_mappability=True, weight_column=None):
        if weight_column is None:
            weight_column = self.default_field

        zero_edge_ix = []
        ix = 0
        for row in self._edges.iterrows():
            if row[weight_column] == 0:
                zero_edge_ix.append(ix)
            ix += 1

        for ix in reversed(zero_edge_ix):
            self._edges.remove_row(ix)

        if flush:
            self.flush(update_index=update_index, update_mappability=update_mappability)

    def marginals(self, weight_column=None):
        """
        Get the marginals vector of this Hic matrix.
        """
        if weight_column is None:
            weight_column = self.default_field

        # prepare marginals dict
        marginals = np.zeros(len(self.regions), float)

        logger.debug("Calculating marginals...")
        with RareUpdateProgressBar(max_value=len(self.edges), silent=config.hide_progressbars) as pb:
            for i, edge in enumerate(self.edges(lazy=True)):
                marginals[edge.source] += getattr(edge, weight_column)
                if edge.source != edge.sink:
                    marginals[edge.sink] += getattr(edge, weight_column)
                pb.update(i)

        return marginals

    def possible_contacts(self):
        logger.info("Calculating possible counts")
        regions = list(self.regions)
        chromosomes = self.chromosomes()

        cb = self.chromosome_bins
        chromosome_max_distance = defaultdict(int)
        max_distance = 0
        chromosome_subtractions = dict()
        for chromosome in chromosomes:
            start, stop = cb[chromosome]
            max_distance = max(max_distance, stop - start)
            chromosome_max_distance[chromosome] = max(chromosome_max_distance[chromosome], stop - start)
            chromosome_subtractions[chromosome] = np.zeros(stop - start,
                                                           dtype='int32')

        chromosome_mappable = defaultdict(int)
        chromosome_unmappable = defaultdict(set)
        for i, mappable in enumerate(self.mappable()):
            chromosome = regions[i].chromosome
            if not mappable:  # unmappable
                s = chromosome_subtractions[chromosome]
                o = cb[chromosome][0]
                ix = i - o
                # horizontal
                s[0: len(s) - ix] += 1
                # vertical
                for j in range(1, ix + 1):
                    if ix - j not in chromosome_unmappable[chromosome]:
                        s[j] += 1
                chromosome_unmappable[chromosome].add(ix)
            else:
                chromosome_mappable[chromosome] += 1

        inter_total = 0
        intra_total = [0] * max_distance
        chromosome_intra_total = dict()
        for chromosome, d in chromosome_max_distance.items():
            chromosome_intra_total[chromosome] = [0] * d

        for i, chromosome in enumerate(chromosomes):
            start, stop = cb[chromosome]
            count = stop - start

            # intra-chromosomal
            s = chromosome_subtractions[chromosomes[i]]
            for distance in range(0, count):
                intra_total[distance] += count - distance - s[distance]
                chromosome_intra_total[chromosome][distance] += count - distance - s[distance]

            # inter-chromosomal
            for j in range(i + 1, len(chromosomes)):
                count_mappable = chromosome_mappable[chromosomes[i]]
                count2_mappable = chromosome_mappable[chromosomes[j]]
                inter_total += count_mappable * count2_mappable

        return intra_total, chromosome_intra_total, inter_total

    def expected_values(self, selected_chromosome=None):
        if self._expected_values is not None:
            intra_expected, chromosome_intra_expected, inter_expected = self._expected_values
        else:
            # get all the bins of the different chromosomes
            chromosome_bins = self.chromosome_bins
            chromosome_dict = defaultdict(list)

            chromosome_max_distance = defaultdict(int)
            max_distance = 0
            for chromosome, (start, stop) in chromosome_bins.items():
                max_distance = max(max_distance, stop - start)
                chromosome_max_distance[chromosome] = max(chromosome_max_distance[chromosome], stop - start)

                for i in range(start, stop):
                    chromosome_dict[i] = chromosome

            chromosome_intra_sums = dict()
            chromosome_intra_expected = dict()
            for chromosome, d in chromosome_max_distance.items():
                chromosome_intra_sums[chromosome] = [0.0] * d
                chromosome_intra_expected[chromosome] = [0.0] * d

            # get the sums of edges at any given distance
            marginals = [0.0] * len(self.regions)
            inter_sums = 0.0
            intra_sums = [0.0] * max_distance
            with RareUpdateProgressBar(max_value=len(self.edges), prefix='Expected') as pb:
                for i, edge in enumerate(self.edges(lazy=True)):
                    source, sink = edge.source, edge.sink
                    weight = getattr(edge, self.default_field)

                    source_chromosome = chromosome_dict[source]
                    sink_chromosome = chromosome_dict[sink]

                    marginals[source] += weight
                    marginals[sink] += weight

                    if sink_chromosome != source_chromosome:
                        inter_sums += weight
                    else:
                        distance = sink - source
                        intra_sums[distance] += weight
                        chromosome_intra_sums[source_chromosome][distance] += weight
                    pb.update(i)

            intra_total, chromosome_intra_total, inter_total = self.possible_contacts()

            # expected values
            inter_expected = 0 if inter_total == 0 else inter_sums/inter_total

            intra_expected = [0.0] * max_distance
            bin_size = self.bin_size
            distances = []
            for d in range(max_distance):
                distances.append(bin_size * d)

                # whole genome
                count = intra_total[d]
                if count > 0:
                    intra_expected[d] = intra_sums[d] / count

            # chromosomes
            for chromosome in chromosome_intra_expected:
                for d in range(chromosome_max_distance[chromosome]):
                    chromosome_count = chromosome_intra_total[chromosome][d]
                    if chromosome_count > 0:
                        chromosome_intra_expected[chromosome][d] = chromosome_intra_sums[chromosome][d] / chromosome_count

            self._expected_values = intra_expected, chromosome_intra_expected, inter_expected

        if selected_chromosome is not None:
            return chromosome_intra_expected[selected_chromosome]

        return intra_expected, chromosome_intra_expected, inter_expected

    def iter_edge_attribute(self, attribute):
        for value in getattr(self._edges.cols, attribute):
            yield value

    def scaling_factor(self, matrix, weight_column=None):
        """
        Compute the scaling factor to another matrix.

        Calculates the ratio between the number of contacts in
        this Hic object to the number of contacts in another
        Hic object.

        :param matrix: A :class:`~Hic` object
        :param weight_column: Name of the column to calculate the scaling factor on
        :return: float
        """
        if weight_column is None:
            weight_column = self.default_field

        logger.info("Calculating scaling factor...")
        m1_sum = 0
        for v1 in self.iter_edge_attribute(weight_column):
            if np.isfinite(v1):
                m1_sum += v1

        m2_sum = 0
        for v2 in matrix.iter_edge_attribute(weight_column):
            if np.isfinite(v2):
                m2_sum += v2

        scaling_factor = m1_sum / m2_sum
        logger.debug("Scaling factor: {}".format(scaling_factor))
        return scaling_factor

    def get_combined_matrix(self, matrix, key=None, scaling_factor=None, weight_column=None):
        """
        Return a :class:`~HicMatrix` where values above the diagonal
        are from this object and values below the diagonal are from
        another :class:`~Hic` object.

        "Above the diagonal" refers to the diagonal of the complete
        Hic object, not the diagonal of the returned matrix.

        :param matrix: Another :class:`~Hic` object
        :param key: A matrix selector. Use tuple to selct row and
                    columns, also see __getitem__
        :param scaling_factor: Factor to scale the hic values. If None,
                               will be computed using
                               :func:`~Hic.scaling_factor`.
        :param weight_column: Name of the column used to combine matrices into one
        :return: :class:`~HicMatrix`
        """
        if key is None:
            key = slice(0, None, None)

        if scaling_factor is None:
            scaling_factor = self.scaling_factor(matrix, weight_column=weight_column)

        m_top = self.as_matrix(key=key, values_from=weight_column)

        # find diagonal
        row_region = m_top.row_regions[0]
        matching_index = None
        for i, col_region in enumerate(m_top.col_regions):
            if col_region == row_region:
                matching_index = i

        if matching_index is None:
            col_region = m_top.col_regions[0]
            for i, row_region in enumerate(m_top.row_regions):
                if col_region == row_region:
                    matching_index = -1 * i

        if matching_index is None:
            return m_top

        # replace diagonal
        m_bottom = matrix.as_matrix(key=key, values_from=weight_column) * scaling_factor
        top_indices = np.triu_indices(m_top.shape[0], matching_index, m_top.shape[1])
        m_bottom[top_indices] = m_top[top_indices]
        return m_bottom

    def filter(self, edge_filter, queue=False, log_progress=False):
        """
        Filter edges in this object by using a
        :class:`~HicEdgeFilter`.

        :param edge_filter: Class implementing :class:`~HicEdgeFilter`.
                            Must override valid_edge method, ideally sets mask parameter
                            during initialization.
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        :param log_progress: If true, process iterating through all edges
                             will be continuously reported.
        """
        edge_filter.set_hic_object(self)
        if not queue:
            self._edges.filter(edge_filter, _logging=log_progress)
        else:
            self._edges.queue_filter(edge_filter)

    def run_queued_filters(self, log_progress=False):
        """
        Run queued filters.

        :param log_progress: If true, process iterating through all edges
                             will be continuously reported.
        """
        self._edges.run_queued_filters(_logging=log_progress)


class AccessOptimisedRegionMatrixTable(RegionMatrixTable, AccessOptimisedRegionPairs):
    """
    Class with faster access to matrix data, based on :class:`~AccessOptimisedRegionPairs`.
    """

    _classid = 'ACCESSOPTIMISEDREGIONMATRIXTABLE'

    def __init__(self, file_name=None, mode='a', tmpdir=None, additional_fields=None,
                 default_field=None, default_value=0.0,
                 _table_name_nodes='nodes', _table_name_edges='edges'):
        """
        Initialize a :class:`~AccessOptimisedRegionMatrixTable` object.

        :param file_name: Path to a save file
        :param mode: File mode to open underlying file
        :param _table_name_nodes: (Internal) name of the HDF5 node for regions
        :param _table_name_edges: (Internal) name of the HDF5 node for edges
        """

        # private variables
        self.default_field = default_field
        self.default_value = default_value
        self._expected_values = None
        AccessOptimisedRegionPairs.__init__(self, file_name=file_name, mode=mode, additional_fields=additional_fields,
                                            tmpdir=tmpdir, _table_name_nodes=_table_name_nodes,
                                            _table_name_edges=_table_name_edges)

        if default_field is None:
            self.default_field = self.field_names[2]

    def _flush_edge_buffer(self, e_buffer, replace=False,
                           update_index=True, update_mappability=True,
                           clean_zero=True, default_column=None):
        if default_column is None:
            default_column = self.default_field

        try:
            # noinspection PyCompatibility
            e_buffer_items = e_buffer.iteritems()
        except AttributeError:
            e_buffer_items = e_buffer.items()

        # re-arrange edge buffer
        partition_e_buffer = defaultdict(dict)
        for key, edge in e_buffer_items:
            source_partition = self._get_partition_ix(key[0])
            sink_partition = self._get_partition_ix(key[1])
            partition_e_buffer[(source_partition, sink_partition)][key] = e_buffer[key]

        try:
            # noinspection PyCompatibility
            partition_e_buffer_items = partition_e_buffer.iteritems()
        except AttributeError:
            partition_e_buffer_items = partition_e_buffer.items()

        # update current rows
        for partition_key, e_buffer in partition_e_buffer_items:
            if partition_key in self._edge_table_dict:
                edge_table = self._edge_table_dict[partition_key]

                for row in edge_table:
                    key = (row["source"], row["sink"])

                    if key in e_buffer:
                        value = e_buffer[key]
                        # it is a weight
                        try:
                            if replace:
                                row[default_column] = float(value)
                            else:
                                row[default_column] += float(value)
                            row.update()
                        except TypeError:
                            self.add_edge(value, check_nodes_exist=False, flush=False, replace=replace, row=row)
                        del e_buffer[key]
                self.flush(update_index=False, update_mappability=False)

            try:
                # noinspection PyCompatibility
                e_buffer_keys = e_buffer.iterkeys()
            except AttributeError:
                e_buffer_keys = e_buffer.keys()

            # flush remaining buffer
            for source, sink in e_buffer_keys:
                key = (source, sink)
                value = e_buffer[key]
                try:
                    v = float(value)
                    if v == 0:
                        continue
                    new_edge = self._edge_from_dict({'source': source, 'sink': sink, default_column: v})
                    self.add_edge(new_edge, check_nodes_exist=False, flush=False)
                except TypeError:
                    self.add_edge(value, check_nodes_exist=False, flush=False)

        self.flush(update_index=update_index, update_mappability=update_mappability)
        if clean_zero:
            self._remove_zero_edges(update_index=update_index, update_mappability=update_mappability,
                                    weight_column=default_column)

    def _remove_zero_edges(self, flush=True, update_index=True,
                           update_mappability=True, weight_column=None):
        if weight_column is None:
            weight_column = self.default_field

        for edge_table in self._edge_table_iter():
            zero_edge_ix = []
            ix = 0
            for row in edge_table.iterrows():
                if row[weight_column] == 0:
                    zero_edge_ix.append(ix)
                ix += 1

            for ix in reversed(zero_edge_ix):
                edge_table.remove_row(ix)

        if flush:
            self.flush(update_index=update_index, update_mappability=update_mappability)

    def filter(self, edge_filter, queue=False, log_progress=not config.hide_progressbars):
        """
        Filter edges in this object by using a
        :class:`~HicEdgeFilter`.

        :param edge_filter: Class implementing :class:`~HicEdgeFilter`.
                            Must override valid_edge method, ideally sets mask parameter
                            during initialization.
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        :param log_progress: If true, process iterating through all edges
                             will be continuously reported.
        """
        total = 0
        filtered = 0
        if not queue:
            with RareUpdateProgressBar(max_value=sum(1 for _ in self._edge_table_iter()),
                                       silent=not log_progress) as pb:
                for i, edge_table in enumerate(self._edge_table_iter()):
                    stats = edge_table.filter(edge_filter, _logging=False)
                    for key, value in stats.items():
                        if key != 0:
                            filtered += stats[key]
                        total += stats[key]
                    pb.update(i)
            if log_progress:
                logger.info("Total: {}. Filtered: {}".format(total, filtered))
        else:
            for edge_table in self._edge_table_iter():
                edge_table.queue_filter(edge_filter)

    def run_queued_filters(self, log_progress=not config.hide_progressbars):
        """
        Run queued filters.

        :param log_progress: If true, process iterating through all edges
                             will be continuously reported.
        """
        total = 0
        filtered = 0
        with RareUpdateProgressBar(max_value=sum(1 for _ in self._edge_table_iter()),
                                   silent=not log_progress) as pb:
            for i, edge_table in enumerate(self._edge_table_iter()):
                stats = edge_table.run_queued_filters(_logging=False)
                for key, value in stats.items():
                    if key != 0:
                        filtered += stats[key]
                    total += stats[key]
                pb.update(i)
        if log_progress:
            logger.info("Total: {}. Filtered: {}".format(total, filtered))

    def iter_edge_attribute(self, attribute):
        for edge_table in self._edge_table_iter():
            for value in getattr(edge_table.cols, attribute):
                yield value


class Hic(RegionMatrixTable):
    """
    Class for working with Hi-C data.

    Examples:

    .. code:: python

        hic = Hic(file_name="/path/to/save/file")

        # load genomic regions
        genome = Genome.from_folder("/path/to/fasta/folder")
        regions = genome.get_regions("HindIII")
        hic.add_regions(regions)

        # load edges
        edges = []
        edges.append(HicEdge(source=10, sink=23, weight=3)
        edges.append(HicEdge(source=8, sink=9, weight=57)
        # ...
        hic.add_edges(edges)
    """

    _classid = 'HIC'

    class HicRegionAnnotationDescription(t.IsDescription):
        bias = t.Float32Col(pos=0, dflt=1)

    def __init__(self, data=None, file_name=None, mode='a', tmpdir=None,
                 _table_name_nodes='nodes', _table_name_edges='edges',
                 _table_name_node_annotations='node_annot'):

        """
        Initialize a :class:`~Hic` object.

        :param data: Can be the path to an XML file denoting a Hic object,
                     another Hic object, a :class:`~FragmentMappedReadPairs`
                     object, or a path to a save file. In the latter case,
                     this parameter may replace file_name, but only if
                     file_name is None.
        :param file_name: Path to a save file
        :param _table_name_nodes: (Internal) name of the HDF5 node for regions
        :param _table_name_edges: (Internal) name of the HDF5 node for edges
        """

        # parse potential unnamed argument
        if data is not None:
            # data is file name
            if isinstance(data, string_types) and file_name is None:
                file_name = data
                data = None

        if file_name is not None:
            file_name = os.path.expanduser(file_name)

        RegionMatrixTable.__init__(self, additional_fields={'weight': t.Float64Col(pos=0)},
                                   file_name=file_name, mode=mode, tmpdir=tmpdir,
                                   default_field='weight',
                                   _table_name_nodes=_table_name_nodes,
                                   _table_name_edges=_table_name_edges)

        if _table_name_node_annotations in self.file.root:
            self._node_annotations = self.file.get_node('/', _table_name_node_annotations)
        elif mode not in ('r', 'r+'):
            self._node_annotations = t.Table(self.file.root, _table_name_node_annotations,
                                             Hic.HicRegionAnnotationDescription)
            self._node_annotations.flush()
        else:
            # compatibility with existing objects
            self._node_annotations = None

        # add data
        self._add_data(data)

    def _add_data(self, data):
        if data is not None:
            # data is existing Hic object
            if isinstance(data, Hic):
                self.load_from_hic(data)
            else:
                try:
                    self.load_read_fragment_pairs(data)
                except AttributeError:
                    raise ValueError("Input data type not recognized")

    def load_read_fragment_pairs(self, pairs):
        """
        Load data from :class:`~fanc.construct.seq.FragmentMappedReadPairs`.

        This method automatically sums up reads mapping to the same
        fragment pairs and creates exactly one edge per fragment pair.

        :param pairs: A :class:`~fanc.construct.seq.FragmentMappedReadPairs`
                      object.
        """
        # add regions
        if len(self._regions) != 0:
            raise RuntimeError("When importing from read pairs you MUST start from an empty data set!")

        self.add_regions(pairs.regions())

        self.disable_indexes()

        l = len(pairs)

        pair_counter = 0
        with RareUpdateProgressBar(max_value=l, silent=config.hide_progressbars) as pb:
            chromosomes = self.chromosomes()
            for ix1 in range(len(chromosomes)):
                chromosome1 = chromosomes[ix1]
                for ix2 in range(ix1, len(chromosomes)):
                    chromosome2 = chromosomes[ix2]
                    logger.info("Processing pair {}-{}".format(chromosome1, chromosome2))
                    edge_buffer = defaultdict(int)
                    for pair in pairs.pairs_by_chromosomes(chromosome1, chromosome2):
                        source, sink = pair.left.fragment.ix, pair.right.fragment.ix
                        edge_buffer[(source, sink)] += 1
                        pair_counter += 1
                        pb.update(pair_counter)

                    try:
                        # noinspection PyCompatibility
                        e_buffer_items = edge_buffer.iteritems()
                    except AttributeError:
                        e_buffer_items = edge_buffer.items()

                    for (source, sink), weight in e_buffer_items:
                        self.add_edge(Edge(source=source, sink=sink, weight=weight), replace=True, flush=False,
                                      check_nodes_exist=False)
                    self.flush(update_index=False)

        self.flush(update_index=True)
        self.enable_indexes()

    def load_from_hic(self, hic, _edges_by_overlap_method=_edge_overlap_split_rao):
        """
        Load data from another :class:`~Hic` object.

        :param hic: Another :class:`~Hic` object
        :param _edges_by_overlap_method: A function that maps reads from
                                         one genomic region to others using
                                         a supplied overlap map. By default
                                         it uses the Rao et al. (2014) method.
                                         See :func:`~_edge_overlap_split_rao`
        """
        # if we do not have any nodes in this Hi-C object...
        if len(self.regions) == 0:
            logger.info("Copying Hi-C")
            # ...simply import everything
            with RareUpdateProgressBar(max_value=len(hic.regions), silent=config.hide_progressbars) as pb:
                for i, region in enumerate(hic.regions()):
                    self.add_region(region, flush=False)
                    pb.update(i)
            self.flush()

            with RareUpdateProgressBar(max_value=len(hic.edges), silent=config.hide_progressbars) as pb:
                for i, edge in enumerate(hic.edges()):
                    self.add_edge(edge, check_nodes_exist=False, flush=False)
                    pb.update(i)
            self.flush()
            self.bias_vector(hic.bias_vector())
        # if already have nodes in this HiC object...
        else:
            logger.info("Binning Hi-C contacts")
            # create region "overlap map"
            overlap_map = _get_overlap_map(hic.regions(), self.regions())

            self.disable_indexes()
            edge_counter = 0
            with RareUpdateProgressBar(max_value=len(hic.edges), silent=config.hide_progressbars) as pb:
                chromosomes = hic.chromosomes()
                for i in range(len(chromosomes)):
                    for j in range(i, len(chromosomes)):
                        logger.info("Chromosomes: {}-{}".format(chromosomes[i], chromosomes[j]))
                        edges = defaultdict(int)
                        for edge in hic.edge_subset(key=(chromosomes[i], chromosomes[j])):
                            old_source, old_sink = edge.source, edge.sink
                            old_weight = getattr(edge, hic.default_field)

                            for new_edge in _edges_by_overlap_method([old_source, old_sink, old_weight], overlap_map):
                                if new_edge[2] != 0:
                                    edges[(new_edge[0], new_edge[1])] += new_edge[2]

                            edge_counter += 1
                            pb.update(edge_counter)

                        for (source, sink), weight in viewitems(edges):
                            self.add_edge(Edge(source=source, sink=sink, weight=weight), flush=False)

            self.flush(update_index=True)
            self.enable_indexes()

    def copy(self, file_name, tmpdir=None):
        return Hic(data=self, file_name=file_name, tmpdir=tmpdir, mode='w')

    def bin(self, bin_size, file_name=None):
        """
        Map edges in this object to equi-distant bins.

        :param bin_size: Bin size in base pairs
        :param file_name: File name of the new, binned Hic object
        :return: :class:`~Hic` object
        """
        # find chromosome lengths
        logger.info("Constructing binned genome...")
        chromosomes = self.chromosomes()
        chromosome_sizes = {chromosome: 0 for chromosome in chromosomes}
        for region in self.regions():
            if chromosome_sizes[region.chromosome] < region.end:
                chromosome_sizes[region.chromosome] = region.end

        chromosome_list = []
        for chromosome in chromosomes:
            chromosome_list.append(Chromosome(name=chromosome, length=self.chromosome_lens[chromosome]))

        genome = Genome(chromosomes=chromosome_list)
        regions = genome.get_regions(bin_size)
        genome.close()

        logger.info("Binning edges...")
        hic = self.__class__(file_name=file_name, mode='w')
        hic.add_regions(regions)
        regions.close()
        hic.load_from_hic(self)

        return hic

    @classmethod
    def from_hic(cls, hics, file_name=None, tmpdir=None, only_intrachromosomal=False):
        if isinstance(hics, Hic):
            hics = [hics]

        logger.info("Checking if regions are identical")
        identical = True
        for i in range(1, len(hics)):
            if len(hics[i].regions) != len(hics[0].regions):
                identical = False
                break

            for self_region, hic_region in zip(hics[0].regions, hics[i].regions):
                if self_region.chromosome != hic_region.chromosome:
                    identical = False
                    break
                if self_region.start != hic_region.start:
                    identical = False
                    break
                if self_region.end != hic_region.end:
                    identical = False
                    break

        merged_hic = cls(file_name=file_name, tmpdir=tmpdir, mode='w')
        merged_hic.disable_indexes()
        if not identical:
            logger.warning("Regions in your Hi-C objects are not identical. Attempting a merge, "
                           "but it will probably be painfully slow. Ensure identical regions before a"
                           "merge by using the same FASTA/genome object for building both Hi-C objects.")
            merged_hic.merge(hics)
        else:
            merged_hic.add_regions(hics[0].regions)

            chromosomes = hics[0].chromosomes()
            for i in range(len(chromosomes)):
                r2 = range(i, i + 1) if only_intrachromosomal else range(i, len(chromosomes))
                for j in r2:
                    logger.info("Chromosomes: {}-{}".format(chromosomes[i], chromosomes[j]))
                    edges = defaultdict(int)
                    for hic in hics:
                        for edge in hic.edge_subset(key=(chromosomes[i], chromosomes[j])):
                            key = (edge.source, edge.sink)
                            edges[key] += edge.weight

                    for (source, sink), weight in viewitems(edges):
                        merged_hic.add_edge(Edge(source=source, sink=sink, weight=weight), flush=False)
            merged_hic.flush()
        merged_hic.enable_indexes()
        return merged_hic

    def _merge(self, hic, _edge_buffer_size=5000000):
        """
        Merge this object with another :class:`~Hic` object.

        First merges genomic regions, then merges edges.
        It is strongly advised that the genomic regions in
        both objects are the same, although this method will attempt to
        "translate" regions from one object to the other if
        this is not the case.

        :param hic: :class:`~Hic` object to be merged into this one
        """
        ix_conversion = {}

        # check if regions are identical (saves a lot of time)
        logger.info("Checking if regions are identical")
        identical = True
        region_counter = 0
        for self_region, hic_region in zip(self.regions(), hic.regions()):
            if self_region.chromosome != hic_region.chromosome:
                identical = False
                break
            if self_region.start != hic_region.start:
                identical = False
                break
            if self_region.end != hic_region.end:
                identical = False
                break
            ix_conversion[region_counter] = region_counter
            region_counter += 1

        if region_counter < len(hic.regions):
            identical = False

        if not identical:
            ix_conversion = {}
            # merge genomic regions
            logger.info("Merging genomic regions...")

            l = len(hic.regions)

            with RareUpdateProgressBar(max_value=l, silent=config.hide_progressbars) as pb:
                for i, region in enumerate(hic.regions):
                    ix = self._get_region_ix(region)
                    if ix is None:
                        ix = self.add_region([region.chromosome, region.start, region.end], flush=False)
                    ix_conversion[region.ix] = ix
                    pb.update(i)
                self._regions.flush()

        # merge edges
        logger.info("Merging contacts...")
        edge_buffer = {}
        l = len(hic.edges)
        with RareUpdateProgressBar(max_value=l, silent=config.hide_progressbars) as pb:
            for i, merge_edge in enumerate(hic.edges):
                merge_source = ix_conversion[merge_edge.source]
                merge_sink = ix_conversion[merge_edge.sink]
                merge_weight = merge_edge.weight

                if merge_source > merge_sink:
                    merge_source, merge_sink = merge_sink, merge_source

                edge_buffer[(merge_source, merge_sink)] = merge_weight

                pb.update(i)

                if len(edge_buffer) > _edge_buffer_size:
                    logger.info("Flushing buffer...")
                    self._flush_edge_buffer(edge_buffer, replace=False,
                                            update_index=False, update_mappability=False)
                    edge_buffer = {}
        logger.info("Final flush...")
        self._flush_edge_buffer(edge_buffer, replace=False,
                                update_index=False, update_mappability=False)

    def merge(self, hic_or_hics, _edge_buffer_size=5000000):
        """
        Merge this object with other :class:`~Hic` objects.

        First merges genomic regions, then merges edges.
        It is strongly advised that the genomic regions in
        both objects are the same, although this method will attempt to
        "translate" regions from one object to the other if
        this is not the case.

        :param hic_or_hics: :class:`~Hic` object or a list
                            of :class:`~Hic` objects to be
                            merged into this one
        """
        if isinstance(hic_or_hics, Hic):
            hic_or_hics = (hic_or_hics,)

        try:
            for hic in hic_or_hics:
                logger.info("Merging {}".format(hic.file_name))
                self._merge(hic, _edge_buffer_size=_edge_buffer_size)
        except TypeError:
            logger.error('{} is not a Hic object or an iterable'.format(hic_or_hics))
            raise

        logger.info("Removing zero edges")
        self._remove_zero_edges(update_index=True, update_mappability=True)

    def to_cooler(self, path):
        """
        Export Hi-C data as cooler file. Only contacts that have not been
        filtered are exported.
        https://github.com/mirnylab/cooler/
        If input Hi-C matrix is uncorrected, the uncorrected matrix is stored.
        If it is corrected, the uncorrected matrix is stored and the bias vector.
        Cooler always calculates corrected matrix on-the-fly from the uncorrected
        matrix and the bias vector.
        :param path: Output path for cooler file
        """
        try:
            import cooler
        except ImportError:
            logger.error("Cannot import cooler. Install cooler 'pip install cooler'.")
            raise
        from cooler.io import parse_cooler_uri
        import h5py
        import itertools as it
        n_contacts = len(self)
        contact_dtype = [("source", np.int_), ("sink", np.int_), ("weight", np.float_)]
        bias = self.bias_vector()
        is_corrected = ~np.all(np.isclose(bias, 1.))
        def contact_generator():
            '''
            Generator yielding all unfiltered contacts in Hi-C object as tuples
            (source node id, sink node id, edge weight).
            '''
            for node in self.edges(lazy=True):
                yield (node.source, node.sink, node.weight)
        logging.info("Loading contacts")
        contact_array = np.fromiter(contact_generator(), dtype=contact_dtype, count=n_contacts)
        logging.info("Sorting contacts")
        order = np.argsort(contact_array, order=("source", "sink"))
        # Cooler stores uncorrected counts and bias vector, so if matrix is corrected
        # we have to "uncorrect" it here
        if is_corrected:
            bias_flat = np.take(bias, contact_array["source"])
            bias_flat *= np.take(bias, contact_array["sink"])
            counts = np.rint(contact_array["weight"]/bias_flat).astype(np.int_)
        else:
            counts = np.rint(contact_array["weight"]).astype(np.int_)
        contact_dict = {
            "bin1_id": contact_array["source"][order],
            "bin2_id": contact_array["sink"][order],
            "count": counts[order],
        }
        region_dicts = [{"chrom": r.chromosome, "start": r.start - 1, "end": r.end} for r in self.regions()]
        region_df = p.DataFrame(region_dicts)
        logging.info("Writing cooler")
        cooler.io.create(cool_uri=path, bins=region_df, pixels=contact_dict)
        if is_corrected:
            cool_path, group_path = parse_cooler_uri(path)
            logging.info("Writing bias vector")
            # Copied this section from
            # https://github.com/mirnylab/cooler/blob/356a89f6a62e2565f42ff13ec103352f20d251be/cooler/cli/balance.py#L195
            with h5py.File(cool_path, 'r+') as h5:
                grp = h5[group_path]
                # add the bias column to the file
                h5opts = dict(compression='gzip', compression_opts=6)
                grp['bins'].create_dataset("weight", data=bias, **h5opts)

    def flush(self, flush_nodes=True, flush_edges=True,
              update_index=True, update_mappability=True):
        """
        Write data to file and flush buffers.

        :param flush_nodes: Flush nodes tables
        :param flush_edges: Flush edges table
        :param update_index: Update mask indices in edges table
        """
        RegionMatrixTable.flush(self, flush_nodes=flush_nodes, flush_edges=flush_edges,
                                update_index=update_index, update_mappability=update_mappability)
        self._node_annotations.flush()

    def filter_diagonal(self, distance=0, queue=False):
        """
        Convenience function that applies a :class:`~DiagonalFilter`.

        :param distance: Distance from the diagonal up to which matrix entries
                         will be filtered. The default, 0, filters only the
                         diagonal itself.
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        """
        mask = self.add_mask_description('diagonal',
                                         'Mask the diagonal of the Hic matrix (up to distance %d)' % distance)
        diagonal_filter = DiagonalFilter(distance=distance, mask=mask)
        self.filter(diagonal_filter, queue)

    def filter_low_coverage_regions(self, rel_cutoff=None, cutoff=None, queue=False):
        """
        Convenience function that applies a :class:`~LowCoverageFilter`.

        The cutoff can be provided in two ways:
        1. As an absolute threshold. Regions with contact count below this
        absolute threshold are filtered
        2. As a fraction relative to the median contact count of all regions.

        If both is supplied, whichever threshold is lower will be selected.

        If no parameter is supplied, rel_cutoff will be chosen as 0.1.

        :param rel_cutoff: A cutoff as a fraction (0-1) of the median contact count of all
                           regions.
        :param cutoff: A cutoff in absolute contact counts (can be float) below
                       which regions are considered "low coverage"
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        """
        if cutoff is None and rel_cutoff is None:
            rel_cutoff = 0.1

        mask = self.add_mask_description('low_coverage',
                                         'Mask low coverage regions in the Hic matrix '
                                         '(absolute cutoff {:.4}, relative '
                                         'cutoff {:.1%}'.format(float(cutoff) if cutoff else 0., float(rel_cutoff) if rel_cutoff else 0.))

        low_coverage_filter = LowCoverageFilter(self, rel_cutoff=rel_cutoff, cutoff=cutoff, mask=mask)
        self.filter(low_coverage_filter, queue)
    
    def bias_vector(self, vector=None):
        """
        Get or set the bias vector of this Hic matrix.

        :param vector: Numpy float vector. If provided, sets the
                       the bias vector of the object.
        """

        if vector is not None:
            if len(vector) != len(self.regions):
                raise ValueError("Bias vector must be the same length as number of regions "
                                 "(is: {}, should: {})".format(len(vector), len(self.regions)))

            # overwrite biases
            if len(self._node_annotations) == len(vector):
                for i, row in enumerate(self._node_annotations):
                    row['bias'] = vector[i]
                    row.update()
            # create new biases
            else:
                row = self._node_annotations.row
                for value in vector:
                    row['bias'] = value
                    row.append()
            self._node_annotations.flush()
            return vector

        vector = np.ones(len(self.regions))
        if len(self._node_annotations) > 0:
            for i, row in enumerate(self._node_annotations):
                vector[i] = row['bias']

        return vector

    def mappable_regions(self):
        mappability = self.mappable()
        mappable = defaultdict(int)
        for i, region in enumerate(self.regions()):
            if mappability[i]:
                mappable[region.chromosome] += 1
        return mappable

    @property
    def architecture(self):
        import fanc.legacy.architecture.hic_architecture as ha
        return ha.HicArchitecture(self)

    @classmethod
    def from_juicer(cls, juicer_file, juicer_tools_jar_path, genome_file, resolution,
                    norm='NONE', output_file=None, inter_chromosomal=True,
                    chromosomes=None):
        hic = cls(file_name=output_file, mode='w')

        # regions
        logger.info("Building genome {}".format(genome_file))
        genome = Genome.from_string(genome_file)
        regions = genome.get_regions(resolution, chromosomes=chromosomes)
        hic.add_regions(regions)
        genome.close()
        regions.close()

        if chromosomes is None:
            chromosomes = hic.chromosomes()
        chromosome_bins = hic.chromosome_bins

        logger.info("Extracting edges from hic file")
        nan_counter = 0
        for i in range(len(chromosomes)):
            chromosome1 = chromosomes[i]
            offset1 = chromosome_bins[chromosome1][0]

            for j in range(i, len(chromosomes)):
                if i != j and not inter_chromosomal:
                    continue
                chromosome2 = chromosomes[j]
                offset2 = chromosome_bins[chromosome2][0]

                logger.info("{} -- {}".format(chromosome1, chromosome2))

                juicer_command = ['java', '-jar', juicer_tools_jar_path,
                                  'dump', 'observed', norm, juicer_file,
                                  chromosome1, chromosome2, 'BP', str(resolution)]

                juicer_process = subprocess.Popen(juicer_command, stdout=subprocess.PIPE)

                for line in juicer_process.stdout:
                    fields = line.rstrip().split()

                    try:
                        start, end, weight = int(fields[0]), int(fields[1]), float(fields[2])
                    except ValueError:
                        continue

                    start_ix = int(start / resolution) + offset1
                    end_ix = int(end / resolution) + offset2
                    if start_ix > end_ix:
                        start_ix, end_ix = end_ix, start_ix

                    if not np.isfinite(weight):
                        nan_counter += 1
                        continue

                    hic.add_edge([start_ix, end_ix, weight], check_nodes_exist=False, flush=False,
                                 replace=True)

        hic.flush()

        if nan_counter > 0:
            logger.warning("{} contacts could not be imported, "
                           "because they had non-finite values.".format(nan_counter))

        return hic

    def sample(self, n, exact=False, file_name=None):
        if isinstance(n, Hic):
            n = len(n.edges)

        region_pairs = []
        if not exact:
            weights = []
            logger.info("Using sampling with replacement")
            for edge in self.edges(lazy=True):
                region_pairs.append((edge.source, edge.sink))
                weights.append(edge.weight)
            s = sum(weights)
            p = [w / s for w in weights]
        else:
            p = None
            logger.info("Using sampling without replacement")
            for edge in self.edges(lazy=True):
                for i in range(int(edge.weight)):
                    region_pairs.append((edge.source, edge.sink))

        new_hic = self.__class__(file_name=file_name, mode='w')
        new_hic.add_regions(self.regions)
        new_edges = defaultdict(int)
        for new_pair_ix in np.random.choice(len(region_pairs), size=n, replace=not exact, p=p):
            new_edges[region_pairs[new_pair_ix]] += 1
        new_edges = [[source, sink, weight] for (source, sink), weight in new_edges.items()]
        new_hic.add_edges(new_edges)

        return new_hic

    def subset_hic(self, *regions, **kwargs):
        """
        Subset a Hic object by specifying one or more subset regions.

        :param regions: string or GenomicRegion object(s)
        :param kwargs: Supports
                       file_name: destination file name of subset Hic object;
                       tmpdir: if True works in tmp until object is closed
        :return: Hic
        """
        file_name = kwargs.get("file_name", None)
        tmpdir = kwargs.get('tmpdir', None)

        new_hic = AccessOptimisedHic(file_name=file_name, mode='w', tmpdir=tmpdir)

        bias_vector = self.bias_vector()

        new_bias_vector = []
        ix_converter = {}
        ix = 0
        for region_string in regions:
            for region in self.regions(region_string):
                ix_converter[region.ix] = ix
                ix += 1

                new_hic.add_region(region, flush=False)
                new_bias_vector.append(bias_vector[region.ix])
        new_hic.flush()

        for i, region_string1 in enumerate(regions):
            for j in range(i, len(regions)):
                region_string2 = regions[j]
                for edge in self.edge_subset(key=(region_string1, region_string2), lazy=True):
                    source = ix_converter[edge.source]
                    sink = ix_converter[edge.sink]
                    new_hic.add_edge([source, sink, edge.weight], flush=False)
        new_hic.flush()

        new_hic.bias_vector(new_bias_vector)
        return new_hic


class AccessOptimisedHic(Hic, AccessOptimisedRegionMatrixTable):
    """
    Class with faster access to matrix data, based on :class:`~AccessOptimisedRegionPairs`.
    """

    _classid = 'ACCESSOPTIMISEDHIC'

    def __init__(self, data=None, file_name=None, mode='a', tmpdir=None,
                 _table_name_nodes='nodes', _table_name_edges='edges',
                 _table_name_node_annotations='node_annot'):
        """
        Initialize a :class:`~AccessOptimisedHic` object.

        :param data: Can be the path to an XML file denoting a Hic object,
                     another Hic object, a :class:`~FragmentMappedReadPairs`
                     object, or a path to a save file. In the latter case,
                     this parameter may replace file_name, but only if
                     file_name is None.
        :param file_name: Path to a save file
        :param _table_name_nodes: (Internal) name of the HDF5 node for regions
        :param _table_name_edges: (Internal) name of the HDF5 node for edges
        """

        # parse potential unnamed argument
        if data is not None:
            # data is file name
            if isinstance(data, string_types) and file_name is None:
                file_name = data
                data = None

        if file_name is not None:
            file_name = os.path.expanduser(file_name)

        AccessOptimisedRegionMatrixTable.__init__(self, additional_fields={'weight': t.Float64Col(pos=0)},
                                                  file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                  default_field='weight',
                                                  _table_name_nodes=_table_name_nodes,
                                                  _table_name_edges=_table_name_edges)

        if _table_name_node_annotations in self.file.root:
            self._node_annotations = self.file.get_node('/', _table_name_node_annotations)
        elif mode not in ('r', 'r+'):
            self._node_annotations = t.Table(self.file.root, _table_name_node_annotations,
                                             Hic.HicRegionAnnotationDescription)
            self._node_annotations.flush()
        else:
            # compatibility with existing objects
            self._node_annotations = None

        # add data
        self._add_data(data)

    def flush(self, flush_nodes=True, flush_edges=True,
              update_index=True, update_mappability=True,
              silent=False):
        """
        Write data to file and flush buffers.

        :param flush_nodes: Flush nodes tables
        :param flush_edges: Flush edges table
        :param update_index: Update mask indices in edges table
        """
        AccessOptimisedRegionMatrixTable.flush(self, flush_nodes=flush_nodes, flush_edges=flush_edges,
                                               update_index=update_index, update_mappability=update_mappability,
                                               silent=silent)
        self._node_annotations.flush()

    @classmethod
    def from_hic(cls, hics, file_name=None, tmpdir=None, only_intrachromosomal=False):
        if isinstance(hics, Hic):
            hics = [hics]

        logger.info("Checking if regions are identical")
        identical = True
        for i in range(1, len(hics)):
            if len(hics[i].regions) != len(hics[0].regions):
                identical = False
                break

            for self_region, hic_region in zip(hics[0].regions, hics[i].regions):
                if self_region.chromosome != hic_region.chromosome:
                    identical = False
                    break
                if self_region.start != hic_region.start:
                    identical = False
                    break
                if self_region.end != hic_region.end:
                    identical = False
                    break

        if not identical:
            logger.warn("Regions in your Hi-C objects are not identical. Attempting a merge, "
                        "but it will probably be painfully slow. Ensure identical regions before a"
                        "merge by using the same FASTA/genome object for building both Hi-C objects.")
            merged_hic = cls(file_name=file_name, tmpdir=tmpdir, mode='w')
            merged_hic.merge(hics)
            return merged_hic

        merged_hic = cls(file_name=file_name, tmpdir=tmpdir, mode='w')
        merged_hic.add_regions(hics[0].regions)

        # check if edge table partitions are also identical
        try:
            merged_partitions = merged_hic.partitions
            partitions_identical = True
            for hic in hics:
                if hic.partitions != merged_partitions:
                    partitions_identical = False
                    break
        except AttributeError:
            # this is an old-style Hi-C object
            partitions_identical = False

        partition_keys = set()
        for hic in hics:
            for partition_key in hic._edge_table_dict.keys():
                partition_keys.add(partition_key)

        merged_hic.disable_indexes()
        if partitions_identical:
            logger.info("Partitions identical, performing fast merge.")
            for partition_key in partition_keys:

                edge_buffer = defaultdict(int)
                for hic in hics:
                    if partition_key in hic._edge_table_dict:
                        for row in hic._edge_table_dict[partition_key]:
                            source, sink, weight = row['source'], row['sink'], row[hic.default_field]
                            edge_buffer[(source, sink)] += weight

                merged_edge_table = merged_hic._create_edge_table(partition_key[0], partition_key[1])
                merged_row = merged_edge_table.row
                for (source, sink), weight in viewitems(edge_buffer):
                    merged_row['source'] = source
                    merged_row['sink'] = sink
                    merged_row[merged_hic.default_field] = weight
                    merged_row.append()
                merged_edge_table.flush(update_index=False)
        else:
            logger.info("Partition tables not identical, this will be a slower merge.")
            chromosomes = merged_hic.chromosomes()
            for i in range(len(chromosomes)):
                r2 = [i] if only_intrachromosomal else range(i, len(chromosomes))
                for j in r2:
                    logger.info("Chromosomes: {}-{}".format(chromosomes[i], chromosomes[j]))
                    edge_buffer = defaultdict(int)
                    for hic in hics:
                        for edge in hic.edge_subset(key=(chromosomes[i], chromosomes[j]), lazy=True):
                            edge_buffer[(edge.source, edge.sink)] += edge.weight

                    # re-arrange edge buffer
                    partition_edge_buffer = defaultdict(dict)
                    for key, edge in viewitems(edge_buffer):
                        source_partition = merged_hic._get_partition_ix(key[0])
                        sink_partition = merged_hic._get_partition_ix(key[1])
                        partition_edge_buffer[(source_partition, sink_partition)][key] = edge_buffer[key]

                    # update current rows
                    for partition_key, e_buffer in viewitems(partition_edge_buffer):
                        merged_edge_table = merged_hic._create_edge_table(partition_key[0], partition_key[1])

                        merged_row = merged_edge_table.row
                        for (source, sink), weight in viewitems(e_buffer):
                            merged_row['source'] = source
                            merged_row['sink'] = sink
                            merged_row[merged_hic.default_field] = weight
                            merged_row.append()
                        merged_edge_table.flush(update_index=False)
        merged_hic.flush()
        merged_hic.enable_indexes()

        return merged_hic

    def load_from_hic(self, hic, _edges_by_overlap_method=_edge_overlap_split_rao):
        """
        Load data from another :class:`~Hic` object.

        :param hic: Another :class:`~Hic` object
        :param _edges_by_overlap_method: A function that maps reads from
                                         one genomic region to others using
                                         a supplied overlap map. By default
                                         it uses the Rao et al. (2014) method.
                                         See :func:`~_edge_overlap_split_rao`
        """

        try:
            hic.partitions
        except AttributeError:
            return Hic.load_from_hic(self, hic, _edges_by_overlap_method=_edge_overlap_split_rao)

        # if we do not have any nodes in this Hi-C object...
        if len(self.regions) == 0:
            logger.info("Copying Hi-C")
            self.add_regions(hic.regions)
            self.add_edges(hic.edges)
            self.bias_vector(hic.bias_vector())
        # if already have nodes in this HiC object...
        else:
            logger.info("Binning Hi-C contacts")
            # create region "overlap map"
            overlap_map = _get_overlap_map(hic.regions(), self.regions())

            self.disable_indexes()
            edge_counter = 0
            with RareUpdateProgressBar(max_value=len(hic.edges), silent=config.hide_progressbars) as pb:

                for edge_table in hic._edge_table_dict.values():
                    edge_buffer = defaultdict(int)
                    for row in edge_table:
                        for new_edge in _edges_by_overlap_method((row['source'], row['sink'],
                                                                 row[hic.default_field]),
                                                                 overlap_map):
                            if new_edge[2] != 0:
                                edge_buffer[(new_edge[0], new_edge[1])] += new_edge[2]
                        edge_counter += 1
                        pb.update(edge_counter)

                    # re-arrange edge buffer
                    partition_edge_buffer = defaultdict(dict)
                    for key, edge in viewitems(edge_buffer):
                        source_partition = self._get_partition_ix(key[0])
                        sink_partition = self._get_partition_ix(key[1])
                        partition_edge_buffer[(source_partition, sink_partition)][key] = edge_buffer[key]

                    # update current rows
                    for partition_key, e_buffer in viewitems(partition_edge_buffer):
                        this_edge_table = self._create_edge_table(partition_key[0], partition_key[1])

                        merged_row = this_edge_table.row
                        for (source, sink), weight in viewitems(e_buffer):
                            merged_row['source'] = source
                            merged_row['sink'] = sink
                            merged_row[self.default_field] = weight
                            merged_row.append()
                        this_edge_table.flush(update_index=False)

            self.flush(update_index=True, update_mappability=True)
            self.enable_indexes()

    def filter(self, edge_filter, queue=False, log_progress=not config.hide_progressbars):
        """
        Filter edges in this object by using a
        :class:`~HicEdgeFilter`.

        :param edge_filter: Class implementing :class:`~HicEdgeFilter`.
                            Must override valid_edge method, ideally sets mask parameter
                            during initialization.
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        :param log_progress: If true, process iterating through all edges
                             will be continuously reported.
        """
        edge_filter.set_hic_object(self)

        AccessOptimisedRegionMatrixTable.filter(self, edge_filter, queue=queue, log_progress=log_progress)


def load_hic(file_name, mode='r', tmpdir=None, _edge_table_name='edges'):
    f = t.open_file(file_name, mode='r')
    n = f.get_node('/' + _edge_table_name)
    if isinstance(n, MaskedTable):
        hic_class = Hic
    elif isinstance(n, t.group.Group):
        hic_class = AccessOptimisedHic
    else:
        raise ValueError("%s is not a valid Hi-C object file" % file_name)

    f.close()
    return hic_class(file_name=file_name, mode=mode, tmpdir=tmpdir)


class HicEdgeFilter(with_metaclass(ABCMeta, MaskFilter)):
    """
    Abstract class that provides filtering functionality for the
    edges/contacts in a :class:`~Hic` object.

    Extends MaskFilter and overrides valid(self, row) to make
    :class:`~HicEdge` filtering more "natural".

    To create custom filters for the :class:`~Hic` object, extend this
    class and override the valid_edge(self, edge) method.
    valid_edge should return False for a specific :class:`~HicEdge` object
    if the object is supposed to be filtered/masked and True
    otherwise. See :class:`~DiagonalFilter` for an example.

    Pass a custom filter to the :func:`~Hic.filter` method in :class:`~Hic`
    to apply it.
    """

    def __init__(self, mask=None):
        """
        Initialize HicEdgeFilter.

        :param mask: The Mask object that should be used to mask
                     filtered :class:`~HicEdge` objects. If None the default
                     Mask will be used.
        """
        super(HicEdgeFilter, self).__init__(mask)
        self._hic = None

    @abstractmethod
    def valid_edge(self, edge):
        """
        Determine if a :class:`~HicEdge` object is valid or should
        be filtered.

        When implementing custom HicEdgeFilter this method must be
        overridden. It should return False for :class:`~HicEdge` objects that
        are to be fitered and True otherwise.

        Internally, the :class:`~Hic` object will iterate over all HicEdge
        instances to determine their validity on an individual
        basis.

        :param edge: A :class:`~HicEdge` object
        :return: True if :class:`~HicEdge` is valid, False otherwise
        """
        pass

    def set_hic_object(self, hic_object):
        """
        Set the :class:`~Hic` instance to be filtered by this
        HicEdgeFilter.

        Used internally by :class:`~Hic` instance.

        :param hic_object: :class:`~Hic` object
        """
        self._hic = hic_object

    def valid(self, row):
        """
        Map valid_edge to MaskFilter.valid(self, row).

        :param row: A pytables Table row.
        :return: The boolean value returned by valid_edge.
        """
        edge = self._hic._row_to_edge(row, lazy=True)
        return self.valid_edge(edge)


class DiagonalFilter(HicEdgeFilter):
    """
    Filter contacts in the diagonal of a :class:`~Hic` matrix.
    """
    def __init__(self, distance=0, mask=None):
        """
        Initialize filter with chosen parameters.

        :param distance: Distance from the diagonal up to which
                         contacts will be filtered
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered edges.
        """
        HicEdgeFilter.__init__(self, mask=mask)
        self.distance = distance

    def valid_edge(self, edge):
        """
        Check if an edge is on (or near) the diagonal of the :class:`~Hic` matrix.
        """
        if abs(edge.source-edge.sink) <= self.distance:
            return False
        return True


class LowCoverageFilter(HicEdgeFilter):
    """
    Filter a :class:`~HicEdge` if it connects a region that
    does not have a contact count larger than a specified
    cutoff.

    If the cutoff is not provided, it is automatically
    chosen at 10% of the mean contact count of all regions.
    """
    def __init__(self, hic_object, cutoff=None, rel_cutoff=None, mask=None):
        """
        Initialize filter with these settings.

        The cutoff can be provided in two ways:
        1. As an absolute threshold. Regions with contact count below this
        absolute threshold are filtered
        2. As a fraction relative to the median contact count of all regions.

        If both is supplied, whichever threshold is lower will be selected.

        :param hic_object: The :class:`~Hic` object that this
                           filter will be called on. Needed for
                           contact count calculation.
        :param rel_cutoff: A cutoff as a fraction (0-1) of the median contact count of all
                           regions. If cutoff and rel_cutoff are None, will be set to 10%
        :param cutoff: A cutoff in absolute contact counts (can be float) below
                       which regions are considered "low coverage"
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered edges.
        """
        HicEdgeFilter.__init__(self, mask=mask)

        self._marginals = hic_object.marginals()
        if cutoff is None and rel_cutoff is None:
            rel_cutoff = 0.1
            logger.info("Using default 10 percent relative coverage as cutoff")

        if cutoff is not None and rel_cutoff is not None:
            cutoff = min(cutoff if cutoff else float("inf"),
                         self.calculate_cutoffs(rel_cutoff)[0] if rel_cutoff else float("inf"))
        elif rel_cutoff is not None:
            cutoff = self.calculate_cutoffs(rel_cutoff)[0]
        logger.info("Final absolute cutoff threshold is {:.4}".format(float(cutoff)))

        self._regions_to_mask = set()
        for i, contacts in enumerate(self._marginals):
            if contacts < cutoff:
                self._regions_to_mask.add(i)
        logger.info("Selected a total of {} ({:.1%}) regions to be masked".format(
            len(self._regions_to_mask), len(self._regions_to_mask)/len(hic_object.regions)))

    def calculate_cutoffs(self, fraction_threshold=0.05):
        lower = np.median(self._marginals[self._marginals > 0])*fraction_threshold
        upper = np.median(self._marginals[self._marginals > 0])+lower
        return lower, upper

    def valid_edge(self, edge):
        """
        Check if an edge falls into a low-coverage region.
        """
        if edge.source in self._regions_to_mask:
            return False
        if edge.sink in self._regions_to_mask:
            return False
        return True


class RegionMatrix(np.ndarray):
    def __new__(cls, input_matrix, col_regions=None, row_regions=None):
        obj = np.asarray(input_matrix).view(cls)
        obj._row_region_trees = None
        obj._col_region_trees = None
        obj.set_col_regions(col_regions)
        obj.set_row_regions(row_regions)
        return obj

    def _interval_tree_regions(self, regions):
        intervals = defaultdict(list)
        for i, region in enumerate(regions):
            interval = intervaltree.Interval(region.start - 1, region.end,
                                             data=i)
            intervals[region.chromosome].append(interval)

        interval_trees = {chromosome: intervaltree.IntervalTree(intervals)
                          for chromosome, intervals in intervals.items()}
        return interval_trees

    def set_row_regions(self, regions):
        self.row_regions = regions
        if regions is not None:
            self._row_region_trees = self._interval_tree_regions(regions)
        else:
            self._row_region_trees = None

    def set_col_regions(self, regions):
        self.col_regions = regions
        if regions is not None:
            self._col_region_trees = self._interval_tree_regions(regions)
        else:
            self._col_region_trees = None

    def __array_finalize__(self, obj):
        if obj is None:
            return

        self.set_row_regions(getattr(obj, 'row_regions', None))
        self.set_col_regions(getattr(obj, 'col_regions', None))

    def __setitem__(self, key, item):
        self._setitem = True
        try:
            if isinstance(self, np.ma.core.MaskedArray):
                out = np.ma.MaskedArray.__setitem__(self, key, item)
            else:
                out = np.ndarray.__setitem__(self, key, item)
        finally:
            self._setitem = False

    def __getitem__(self, index):
        self._getitem = True

        # convert string types into region indexes
        if isinstance(index, tuple):
            if len(index) == 2:
                row_key = self._convert_key(
                    index[0],
                    self._row_region_trees if hasattr(self, '_row_region_trees') else None
                )
                col_key = self._convert_key(
                    index[1],
                    self._col_region_trees if hasattr(self, '_col_region_trees') else None
                )
                index = (row_key, col_key)
            elif len(index) == 1:
                row_key = self._convert_key(index[0], self._row_region_trees)
                col_key = slice(0, len(self.col_regions), 1)
                index = (row_key, )
            else:
                col_key = slice(0, len(self.col_regions), 1)
                row_key = index
                index = row_key
        else:
            row_key = self._convert_key(index, self._row_region_trees)
            try:
                col_key = slice(0, len(self.col_regions), 1)
            except TypeError:
                col_key = None
            index = row_key

        try:
            if isinstance(self, np.ma.core.MaskedArray):
                out = np.ma.MaskedArray.__getitem__(self, index)
            else:
                out = np.ndarray.__getitem__(self, index)
        finally:
            self._getitem = False

        if not isinstance(out, np.ndarray):
            return out

        # get regions
        try:
            row_regions = self.row_regions[row_key]
        except TypeError:
            row_regions = None

        try:
            col_regions = self.col_regions[col_key]
        except TypeError:
            col_regions = None

        if isinstance(row_regions, GenomicRegion):
            out.row_regions = [row_regions]
        else:
            out.row_regions = row_regions

        if isinstance(col_regions, GenomicRegion):
            out.col_regions = [col_regions]
        else:
            out.col_regions = col_regions

        return out

    def __getslice__(self, start, stop):
        return self.__getitem__(slice(start, stop))

    def _convert_key(self, key, region_trees):
        if isinstance(key, string_types):
            key = GenomicRegion.from_string(key)

        if isinstance(key, GenomicRegion):
            start = None
            stop = None
            try:
                key_start = 0 if key.start is None else max(0, key.start - 1)
                key_end = key.end
                for interval in region_trees[key.chromosome][key_start:key_end]:
                    i = interval.data
                    start = min(i, start) if start is not None else i
                    stop = max(i + 1, stop) if stop is not None else i + 1
            except KeyError:
                raise ValueError("Requested chromosome {} was not "
                                 "found in this matrix.".format(key.chromosome))

            if start is None or stop is None:
                raise ValueError("Requested region {} was not found in this matrix.".format(key))

            return slice(start, stop, 1)
        return key

    # def _convert_key(self, key, regions):
    #     if isinstance(key, string_types):
    #         key = GenomicRegion.from_string(key)
    #     if isinstance(key, GenomicRegion):
    #         key_start = 0 if key.start is None else max(0, key.start)
    #         key_end = key.end
    #         start = None
    #         stop = None
    #         for i, region in enumerate(regions):
    #             if region.chromosome == key.chromosome:
    #                 if (key_end is None or region.start <= key_end) and region.end >= key_start:
    #                     if start is None:
    #                         start = i
    #                     stop = i
    #         if start is None or stop is None:
    #             raise ValueError("Requested region {} was not found in this matrix.".format(key))
    #         return slice(start, stop+1, 1)
    #     return key

    def __reduce__(self):
        # Get the parent's __reduce__ tuple
        pickled_state = super(RegionMatrix, self).__reduce__()
        # Create our own tuple to pass to __setstate__
        new_state = pickled_state[2] + (pickle.dumps(self.row_regions), pickle.dumps(self.col_regions))
        # Return a tuple that replaces the parent's __setstate__ tuple with our own
        return pickled_state[0], pickled_state[1], new_state

    def __setstate__(self, state):
        self.row_regions = pickle.loads(state[-2])
        self.col_regions = pickle.loads(state[-1])
        # Call the parent's __setstate__ with the other tuple elements.
        super(RegionMatrix, self).__setstate__(state[0:-2])


def _get_overlap_map(old_regions, new_regions):
    # 1. organize regions in self by chromosome
    new_region_map = {}
    for i, new_region in enumerate(new_regions):
        if not new_region.chromosome in new_region_map:
            new_region_map[new_region.chromosome] = []
        new_region_map[new_region.chromosome].append([new_region.start,new_region.end,i])
        
    # 2. iterate over regions in hic to find overlap
    def _get_overlap(new_region, old_region):
        new_region_length = new_region[1] - new_region[0] + 1
        overlap = min(old_region[1], new_region[1]) - max(old_region[0],new_region[0]) + 1
        return max(0,overlap/new_region_length)
        
    old_to_new = {}
    current_chromosome = ''
    current_ix = 0
    for i, old_region in enumerate(old_regions):
        old_to_new[i] = []
        if current_chromosome != old_region.chromosome:
            current_ix = 0
            current_chromosome = old_region.chromosome
        
        found_overlap = True
        while found_overlap:
            found_overlap = False
            if current_ix < len(new_region_map[current_chromosome]):
                new_region = new_region_map[current_chromosome][current_ix]
                overlap = _get_overlap(new_region, [old_region.start,old_region.end,i])
                if overlap > 0:
                    old_to_new[i].append([new_region[2], overlap])
                    current_ix += 1
                    found_overlap = True
                elif old_region.start > new_region[1]:
                    current_ix += 1
                    found_overlap = True
            
        current_ix -= 1
    
    return old_to_new


def merge_regions(regions):
    sorted_regions = sorted(regions, key=lambda r: (r.chromosome, r.start))

    merged_regions = []
    current_regions = []
    last_end = None
    for region in sorted_regions:
        if len(current_regions) == 0:
            current_regions.append(region)
            last_end = region.end
        elif region.chromosome == current_regions[0].chromosome and region.start < last_end:
            current_regions.append(region)
            last_end = max(last_end, region.end)
        else:
            merged_region = GenomicRegion(chromosome=current_regions[0].chromosome,
                                          start=current_regions[0].start, end=last_end,
                                          strand=current_regions[0].strand)
            merged_regions.append(merged_region)
            current_regions = [region]
            last_end = region.end

    merged_region = GenomicRegion(chromosome=current_regions[0].chromosome,
                                  start=current_regions[0].start, end=last_end,
                                          strand=current_regions[0].strand)
    merged_regions.append(merged_region)

    return merged_regions
