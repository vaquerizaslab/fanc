"""
This module provides a number of classes for the calculation of Hi-C (and other matrices) architectural
features, such as AB domains, methods for boundary/TAD identification, and lots more.

The classes are built in a way, that only the first method call requesting the architectural feature
will trigger the calculation of all necessary values. Every subsequent call will retrieve the pre-calculated
values from memory/file, unless explicitly providing the 'force=True' parameter, which will always
trigger a recalculation. This behavior leads to a very natural syntax, without the need for explicit calls
to a "calculate" method, but with the possibility of buffering data.

Example:

.. code:: python

   hic = kaic.sample_hic()
   ex = ExpectedContacts(hic)  # will not trigger calculation of results yet
   intra_expected_foo = ex.intra_expected()  # triggers calculation
   intra_expected_bar = ex.intra_expected()  # simply retrieves data from memory


A recalculation is also avoided when restoring data from file.

"""


from __future__ import division
from kaic.config import config
from kaic.architecture.architecture import TableArchitecturalFeature, calculateondemand, ArchitecturalFeature
from kaic.architecture.genome_architecture import MatrixArchitecturalRegionFeature, VectorArchitecturalRegionFeature, \
    MatrixArchitecturalRegionFeatureFilter
from kaic.architecture.maxima_callers import MaximaCallerDelta
from kaic.data.genomic import GenomicRegion, HicEdgeFilter, Edge, Hic, Genome, Bedpe
from collections import defaultdict
from kaic.tools.general import ranges, to_slice
from kaic.tools.matrix import apply_sliding_func, kth_diag_indices, trim_stats
from kaic.data.general import FileGroup
from Bio.SeqUtils import GC as calculate_gc_content
import numpy as np
import tables as t
import itertools
from scipy.misc import imresize
from scipy.stats import trim_mean
from scipy.stats.mstats import gmean
from bisect import bisect_left
from kaic.tools.general import RareUpdateProgressBar
from future.utils import string_types
import warnings
import logging
logger = logging.getLogger(__name__)


def cis_trans_ratio(hic, normalise=False):
    """
    Calculate the cis/trans ratio for a Hic object.
    
    :param hic: :class:`~kaic,data.genomic.Hic` object
    :param normalise: If True, will normalise ratio to the possible number of cis/trans contacts
                      in this genome. Makes ratio comparable across different genomes
    :return: tuple (ratio, cis, trans, factor)
    """
    cis = 0
    trans = 0
    regions_dict = hic.regions_dict
    for edge in hic.edges(lazy=True):
        if regions_dict[edge.source].chromosome == regions_dict[edge.sink].chromosome:
            cis += edge.weight
        else:
            trans += edge.weight
    if not normalise:
        return cis / (cis + trans), cis, trans, 1.0
    with PossibleContacts(hic) as possible:
        possible_intra = possible.intra_possible()
        possible_inter = possible.inter_possible()
        f = possible_intra / possible_inter

        return cis / (cis + trans * f), cis, trans, f


class HicArchitecture(object):
    """
    Convenience class to access Hi-C architectural features.
    """
    def __init__(self, hic):
        self.hic = hic

    @property
    def expected_contacts(self, regions=None):
        return ExpectedContacts(self.hic, smooth=False, regions=regions)

    @property
    def possible_contacts(self, regions=None):
        return PossibleContacts(self.hic, regions=regions)

    @property
    def directionality_index(self, window_sizes=2000000):
        if isinstance(window_sizes, int):
            window_sizes = (window_sizes, )
        return DirectionalityIndex(self.hic, window_sizes=window_sizes)


class HicEdgeCollection(MatrixArchitecturalRegionFeature):
    """
    Collection of edges from multiple matrices, accessible via the same regions.
    """
    _classid = 'HICEDGECOLLECTION'

    def __init__(self, hics=None, additional_fields=None, file_name=None, mode='a', tmpdir=None,
                 only_intra_chromosomal=False):
        """
        Initialize :class:`~HicEdgeCollection`.

        :param hics: Iterable of :class:`~kaic,data.genomic.Hic` objects
        :param additional_fields: Any additional meta fields to include in the Hi-C collection
                                  (edge field). Must be PyTables Column description(s)
        :param file_name: Path to save file
        :param mode: File mode ('r' = read-only, 'w' = (over)write, 'a' = append)
        :param tmpdir: Temporary directory
        :param only_intra_chromosomal: If True, ignore inter-chromosomal edges.
        """
        if not isinstance(hics, string_types) and hics is not None:

            if additional_fields is not None:
                if not isinstance(additional_fields, dict) and issubclass(additional_fields, t.IsDescription):
                    # IsDescription subclass case
                    additional_fields = additional_fields.columns
            else:
                additional_fields = dict()

            original_fields = additional_fields.copy()

            self.shared_base_field_names = []
            for i, hic in enumerate(hics):
                field_descriptions = hic._field_dict
                for name, description in field_descriptions.items():

                    if name.startswith("_") or name in {'source', 'sink'}:
                        continue

                    for j in range(len(hics)):
                        new_name = "%s_%d" % (name, j)
                        if new_name in original_fields:
                            raise ValueError("%s already in hic object, please choose different name." % new_name)
                        if new_name in additional_fields:
                            continue

                        if name not in self.shared_base_field_names:
                            self.shared_base_field_names.append(name)

                        description_args = {'dflt': description.dflt,
                                            'pos': len(additional_fields)}
                        if description.__class__ == t.description.StringCol:
                            description_args['itemsize'] = description.itemsize

                        additional_fields[new_name] = description.__class__(**description_args)
            MatrixArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode,
                                                      data_fields=additional_fields,
                                                      regions=hics[0].regions, tmpdir=tmpdir)
        elif (hics is None and file_name is not None) or file_name is None:
            if hics is not None and file_name is None:
                file_name = hics
            MatrixArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir)
            self._calculated = True

        self.hics = hics
        self.only_intra_chromosomal = only_intra_chromosomal

    def _calculate(self, *args, **kwargs):
        chromosomes = self.hics[0].chromosomes()

        # step 1: combine all hic objects into one
        #         and calculate variance
        for chr_i, chromosome1 in enumerate(chromosomes):
            for chr_j in range(chr_i, len(chromosomes)):
                chromosome2 = chromosomes[chr_j]

                if self.only_intra_chromosomal and chromosome1 != chromosome2:
                    continue

                logger.info("Processing chromosomes %s-%s" % (chromosome1, chromosome2))

                edges = dict()
                for i, hic in enumerate(self.hics):
                    logger.info("Processing Hic %d/%d" % (i, len(self.hics)))

                    for edge in hic.edge_subset(key=(chromosome1, chromosome2), lazy=True):
                        key = (edge.source, edge.sink)
                        if key not in edges:
                            edges[key] = {}
                        for field in self.shared_base_field_names:
                            if field not in edges[key]:
                                edges[key][field] = [None] * len(self.hics)
                            edges[key][field][i] = getattr(edge, field, None)

                try:
                    # noinspection PyCompatibility
                    edge_items = edges.iteritems()
                except AttributeError:
                    edge_items = edges.items()

                for key, d in edge_items:
                    try:
                        # noinspection PyCompatibility
                        d_items = d.iteritems()
                    except AttributeError:
                        d_items = d.items()
                    for field, values in d_items:
                        source = key[0]
                        sink = key[1]

                        d = dict()
                        for i, value in enumerate(values):
                            d["%s_%d" % (field, i)] = value

                        self.add_edge(Edge(source=source, sink=sink, **d), flush=False)

        self.flush()


class ExpectedContacts(TableArchitecturalFeature):
    """
    Calculate the expected number of contacts, intra- and inter-chromosomal.

    Intra-chromosomal contacts take into account the distance between two regions,
    inter-chromosomal expected contacts are a genome-wide average.

    :param hic: A :class:`~kaic.data.genomic.RegionMatrixTable` object
                    such as :class:`~kaic.data.genomic.Hic`
    :param file_name: Path to save file location
    :param mode: File mode ('r' = read-only, 'w' = (over)write, 'a' = append)
    :param tmpdir: Temporary directory
    :param smooth: If True, applies curve smoothing (sliding window average)
                   to expected intra-chromosomal contacts.
    :param min_reads: Minimum number of reads in sliding window to apply smooting
                      function
    :param regions: A region selector string, :class:`~kaic.data.genomic.GenomicRegion`, or lists thereof.
    :param weight_column: Name of the column containing the weights/values for the
                          expected value calculation. If None, this will be the default field
                          in the provided :class:`~kaic.data.genomic.RegionMatrixTable`

    """
    _classid = 'EXPECTEDCONTACTS'

    def __init__(self, hic=None, file_name=None, mode='a', tmpdir=None, smooth=False, min_reads=400,
                 regions=None, weight_column=None, _table_name='distance_decay'):
        """
        Initialize an :class:`~ExpectedContacts` object.
        """
        if isinstance(hic, string_types) and file_name is None:
            file_name = hic
            hic = None

        TableArchitecturalFeature.__init__(self, _table_name,
                                           {'distance': t.Int64Col(), 'intra': t.Float64Col(),
                                            'contacts': t.Float64Col(), 'pixels': t.Float64Col()},
                                           file_name=file_name, mode=mode, tmpdir=tmpdir)

        self.hic = hic
        self.smooth = smooth
        self.min_reads = min_reads
        self.regions = regions
        if weight_column is None and hic is not None:
            self.weight_column = self.hic.default_field
        else:
            self.weight_column = weight_column

        if self._calculated:
            self._marginals_array = self._group.marginals

    def _calculate(self):
        """
        Get intra- and inter-chromosomal expected contact counts.
        """

        if self.regions is None:
            regions = [region for region in self.hic.regions]
            edges_iter = self.hic.edges(lazy=True)
        else:
            regions = [region for region in self.hic.subset(self.regions)]
            edges_iter = self.hic.edge_subset(key=(self.regions, self.regions))

        chromosome_offsets = dict()
        chromosome_regions = defaultdict(list)
        regions_dict = dict()
        for i, r in enumerate(regions):
            regions_dict[r.ix] = (i, r.chromosome)
            chromosome_regions[r.chromosome].append(r)
            if r.chromosome not in chromosome_offsets:
                chromosome_offsets[r.chromosome] = i
        max_distance = max([len(chromosome_regions[chromosome]) for chromosome in chromosome_regions])

        try:
            bias_vector = self.hic.bias_vector()
        except AttributeError:
            bias_vector = [1.0] * len(self.hic.regions)
            if self.smooth:
                warnings.warn("This object does not support smoothing, returning unsmoothed values.")
                self.smooth = False

        # get the sums of edges at any given distance
        marginals = [0.0] * len(regions)
        inter_sums = 0.0
        intra_sums = [0.0] * max_distance
        intra_uncorrected = [0] * max_distance
        for edge in edges_iter:
            source, sink = edge.source, edge.sink
            weight = getattr(edge, self.weight_column)

            source_i, source_chromosome = regions_dict[source]
            sink_i, sink_chromosome = regions_dict[sink]

            marginals[source_i] += weight
            marginals[sink_i] += weight

            if sink_chromosome != source_chromosome:
                inter_sums += weight
            else:
                distance = sink - source
                intra_sums[distance] += weight
                intra_uncorrected[distance] += int(weight / (bias_vector[source] * bias_vector[sink]) + 0.5)

        logger.info("Calculating possible counts")
        # getting per-chromosome subtractions
        chromosome_subtractions = dict()
        chromosomes = list(chromosome_regions.keys())
        for chromosome in chromosomes:
            chromosome_subtractions[chromosome] = np.zeros(len(chromosome_regions[chromosome]), dtype='int32')

        chromosome_mappable = defaultdict(int)
        chromosome_unmappable = defaultdict(set)
        for i, marginal in enumerate(marginals):
            chromosome = regions[i].chromosome
            if marginal < 10e-10:  # unmappable
                s = chromosome_subtractions[chromosome]
                o = chromosome_offsets[chromosome]
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
        for i in range(len(chromosomes)):
            count = len(chromosome_regions[chromosomes[i]])

            # intra-chromosomal
            s = chromosome_subtractions[chromosomes[i]]
            for distance in range(0, count):
                intra_total[distance] += count - distance - s[distance]

            # inter-chromosomal
            for j in range(i + 1, len(chromosomes)):
                count_mappable = chromosome_mappable[chromosomes[i]]
                count2_mappable = chromosome_mappable[chromosomes[j]]
                inter_total += count_mappable * count2_mappable

        try:
            inter_expected = inter_sums / inter_total
        except ZeroDivisionError:
            inter_expected = 0.0

        intra_expected = [0.0] * max_distance
        bin_size = self.hic.bin_size
        distances = []
        for distance in range(max_distance):
            distances.append(bin_size * distance)
            count = intra_total[distance]
            if count > 0:
                intra_expected[distance] = intra_sums[distance] / count

        # save marginals in object
        marginals = np.array(marginals)
        self._marginals_array = self.file.create_carray(self._group, 'marginals',
                                                        t.Atom.from_dtype(marginals.dtype),
                                                        marginals.shape)
        self._marginals_array[:] = marginals

        self.data('distance', distances)
        self.meta['inter'] = inter_expected

        # return here if smoothing not requested
        if not self.smooth:
            self.data('intra', intra_expected)
            self.data('contacts', intra_sums)
            self.data('pixels', intra_total)
            return

        # smoothing
        smoothed_intra_sums = np.zeros(len(intra_sums))
        smoothed_intra_total = np.zeros(len(intra_sums))
        for i in range(len(intra_sums)):
            uncorrected_reads = intra_uncorrected[i]
            smoothed_reads = intra_sums[i]
            smoothed_pixels = intra_total[i]
            window_size = 0
            can_extend = True
            # smooth to a minimum number of reads per distance
            while uncorrected_reads < self.min_reads and can_extend:
                window_size += 1
                can_extend = False
                # check if we can increase the window to the left
                if i - window_size >= 0:
                    uncorrected_reads += intra_uncorrected[i-window_size]
                    smoothed_reads += intra_sums[i-window_size]
                    smoothed_pixels += intra_total[i-window_size]
                    can_extend = True
                # check if we can increase the window to the right
                if i + window_size < len(intra_sums):
                    uncorrected_reads += intra_uncorrected[i + window_size]
                    smoothed_reads += intra_sums[i+window_size]
                    smoothed_pixels += intra_total[i+window_size]
                    can_extend = True
            smoothed_intra_sums[i] = smoothed_reads
            smoothed_intra_total[i] = smoothed_pixels

        intra_expected = smoothed_intra_sums/smoothed_intra_total

        self.data('intra', intra_expected)
        self.data('contacts', smoothed_intra_sums)
        self.data('pixels', smoothed_intra_total)

    @calculateondemand
    def marginals(self):
        return self._marginals_array[:]

    @calculateondemand
    def intra_expected(self):
        """
        Get the list of expected intra-chromosomal contacts.
        :return: list of floats
        """
        return self[:, 'intra']

    @calculateondemand
    def inter_expected(self):
        """
        Get the inter-chromosomal expected (normalized) contact count
        :return: float
        """
        return self.meta['inter']

    @calculateondemand
    def distance(self):
        """
        Get a list of distances between regions matching the :func:`~ExpectedContacts.intra_expected`.
        :return: list of ints
        """
        return self[:, 'distance']

    @calculateondemand
    def intra_contacts(self):
        """
        Get a list of observed intra-chromosomal number of contacts grouped by region distance.
        :return: list of floats
        """
        return self[:, 'contacts']

    @calculateondemand
    def intra_pixels(self):
        """
        Get the number of pixels/pairs of bins that correspond to intra-chromosomal contacts at a certain distance.
        :return: list of ints
        """
        return self[:, 'pixels']


class ObservedExpectedRatio(MatrixArchitecturalRegionFeature):
    """
    Calculate the ratio of observed over expected contacts in a :class:`~kaic.data.genomic.RegionMatrixTable`.

    :param hic: A :class:`~kaic.data.genomic.RegionMatrixTable` object
                such as :class:`~kaic.data.genomic.Hic`
    :param file_name: Path to save file location
    :param mode: File mode ('r' = read-only, 'w' = (over)write, 'a' = append)
    :param tmpdir: Temporary directory
    :param regions: A region selector string, :class:`~kaic.data.genomic.GenomicRegion`, or lists thereof.
    :param weight_column: Name of the column containing the weights/values for the
                          expected value calculation. If None, this will be the default field
                          in the provided :class:`~kaic.data.genomic.RegionMatrixTable`
    :return: :class:`~kaic.architecture.genome_architecture.MatrixArchitecturalRegionFeature`
    """
    _classid = 'OBSERVEDEXPECTEDRATIO'

    def __init__(self, hic=None, file_name=None, mode='a', tmpdir=None, regions=None,
                 weight_column='weight', per_chromosome=True, _table_name='observed_expected'):
        self.region_selection = regions

        # are we retrieving an existing object?
        if isinstance(hic, string_types) and file_name is None:
            file_name = hic
            hic = None

        if hic is None and file_name is not None:
            MatrixArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                      _table_name_edges=_table_name)
        else:
            if regions is None:
                regions = hic.regions
                self.region_conversion = {region.ix: region.ix for region in hic.regions}
            else:
                self.region_conversion = {region.ix: i for i, region in enumerate(hic.subset(regions))}
                regions = hic.subset(regions)
            MatrixArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                      data_fields={'ratio': t.Float32Col()}, regions=regions,
                                                      default_field='ratio',
                                                      _table_name_edges=_table_name)
        self.default_value = 1.0
        self.hic = hic
        self.weight_column = weight_column
        self.per_chromosome = per_chromosome
        self._expected_by_chromosome = dict()
        self._expected_inter = None
        self._expected_intra = None

    def _expected(self, chromosome1, chromosome2):
        try:
            return self._expected_by_chromosome[(chromosome1, chromosome2)]
        except KeyError:
            if self.per_chromosome:
                if chromosome1 == chromosome2:
                    with ExpectedContacts(self.hic, regions=chromosome1, weight_column=self.weight_column) as ex:
                        e = ex.intra_expected()
                else:
                    if self._expected_inter is None:
                        with ExpectedContacts(self.hic, weight_column=self.weight_column) as ex_inter:
                            self._expected_inter = ex_inter.inter_expected()
                    e = self._expected_inter
            else:
                if self._expected_inter is None or self._expected_intra is None:
                    with ExpectedContacts(self.hic, regions=self.region_selection,
                                          weight_column=self.weight_column) as ex_inter:
                        self._expected_inter = ex_inter.inter_expected()
                        self._expected_intra = ex_inter.intra_expected()
                if chromosome1 == chromosome2:
                    e = self._expected_intra
                else:
                    e = self._expected_inter
            self._expected_by_chromosome[(chromosome1, chromosome2)] = e

            return e

    def _calculate(self):
        region_selection = self.region_selection if self.region_selection is not None else slice(0, None, None)
        regions_dict = self.hic.regions_dict

        for edge in self.hic.edge_subset(key=(region_selection, region_selection), lazy=True):
            source = edge.source
            new_source = self.region_conversion[source]
            sink = edge.sink
            new_sink = self.region_conversion[sink]
            weight = getattr(edge, self.weight_column)

            e = self._expected(regions_dict[source].chromosome, regions_dict[sink].chromosome)
            try:
                d = new_sink - new_source
                e = e[d]
            except (TypeError, IndexError):
                pass

            if e != 0:
                self.add_edge(Edge(new_source, new_sink, ratio=weight/e), flush=False)
        self.flush()


class ComparisonMatrix(MatrixArchitecturalRegionFeature):
    """
    Compare two :class:`~kaic.data.genomic.RegionMatrixTable` objects.

    Define edge comparison function manually.

    :param matrix1: :class:`~kaic.data.genomic.RegionMatrixTable`, such as :class:`~kaic.data.genomic.Hic`
    :param matrix2: :class:`~kaic.data.genomic.RegionMatrixTable`, such as :class:`~kaic.data.genomic.Hic`
    :param file_name: Path to save file location
    :param mode: File mode ('r' = read-only, 'w' = (over)write, 'a' = append)
    :param tmpdir: Temporary directory
    :param regions: A region selector string, :class:`~kaic.data.genomic.GenomicRegion`, or lists thereof.
                    Will subset both matrices using these region(s) before the calculation
    :param scale_matrices: If True, will scale the matrices naively by artificially increasing the number of
                           reads in one matrix uniformly, so that it has the same number of total contacts
                           as the other matrix.
    :param log2: If True, will log2-transform the output values.
    :param weight_column: Name of the column containing the weights/values for the
                          expected value calculation. If None, this will be the default field
                          in the provided :class:`~kaic.data.genomic.RegionMatrixTable`
    """
    _classid = 'COMPARISONMATRIX'

    def __init__(self, matrix1=None, matrix2=None, comparison_function=lambda x: x[0] / x[1],
                 file_name=None, mode='a', tmpdir=None,
                 regions=None, scale_matrices=False, log2=True,
                 weight_column='weight', ignore_zero=False,
                 _table_name='expected_contacts'):
        self.region_selection = regions

        # are we retrieving an existing object?
        if isinstance(matrix1, string_types) and matrix2 is None and file_name is None:
            file_name = matrix1
            matrix1 = None

        if matrix1 is None and matrix2 is None and file_name is not None:
            MatrixArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                      default_field='fc', _table_name_edges=_table_name)
        else:
            if regions is None:
                regions = matrix1.regions
                self.region_conversion = {region.ix: region.ix for region in matrix1.regions}
            else:
                self.region_conversion = {region.ix: i for i, region in enumerate(matrix1.subset(regions))}
                regions = matrix1.subset(regions)
            MatrixArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                      data_fields={'fc': t.Float32Col()}, regions=regions,
                                                      default_field='fc', _table_name_edges=_table_name)
        self.matrix1 = matrix1
        self.matrix2 = matrix2
        self.weight_column = weight_column
        self.scale_matrices = scale_matrices
        self.log2 = log2
        self.compare = comparison_function
        self.ignore_zero = ignore_zero

    def _calculate(self):
        if self.scale_matrices:
            scaling_factor = self.matrix1.scaling_factor(self.matrix2)
        else:
            scaling_factor = 1.

        chromosomes = self.chromosomes()

        for i in range(len(chromosomes)):
            chromosome1 = chromosomes[i]
            for j in range(i, len(chromosomes)):
                chromosome2 = chromosomes[j]

                edges1 = dict()
                for edge in self.matrix1.edge_subset(key=(chromosome1, chromosome2), lazy=True):
                    try:
                        source = self.region_conversion[edge.source]
                        sink = self.region_conversion[edge.sink]
                    except KeyError:
                        continue

                    edges1[(source, sink)] = getattr(edge, self.weight_column)

                for edge in self.matrix2.edge_subset(key=(chromosome1, chromosome2), lazy=True):
                    try:
                        source = self.region_conversion[edge.source]
                        sink = self.region_conversion[edge.sink]
                    except KeyError:
                        continue

                    edge1_weight = edges1[(source, sink)] if (source, sink) in edges1 else self.matrix1.default_value
                    if self.ignore_zero and edge1_weight == 0 or edge.weight == 0:
                        continue

                    weight = self.compare([edge1_weight, scaling_factor*edge.weight])
                    if self.log2:
                        weight = np.log2(weight)
                    self.add_edge([source, sink, weight], flush=False)

        self.flush()


class FoldChangeMatrix(ComparisonMatrix):
    """
    Calculate the fold-change matrix of two :class:`~kaic.data.genomic.RegionMatrixTable` objects.

    fc_ij = matrix1_ij/matrix2_ij

    :param matrix1: :class:`~kaic.data.genomic.RegionMatrixTable`, such as :class:`~kaic.data.genomic.Hic`
    :param matrix2: :class:`~kaic.data.genomic.RegionMatrixTable`, such as :class:`~kaic.data.genomic.Hic`
    :param file_name: Path to save file location
    :param mode: File mode ('r' = read-only, 'w' = (over)write, 'a' = append)
    :param tmpdir: Temporary directory
    :param regions: A region selector string, :class:`~kaic.data.genomic.GenomicRegion`, or lists thereof.
                    Will subset both matrices using these region(s) before the calculation
    :param scale_matrices: If True, will scale the matrices naively by artificially increasing the number of
                           reads in one matrix uniformly, so that it has the same number of total contacts
                           as the other matrix.
    :param log2: If True, will log2-transform the output values.
    :param weight_column: Name of the column containing the weights/values for the
                          expected value calculation. If None, this will be the default field
                          in the provided :class:`~kaic.data.genomic.RegionMatrixTable`
    """
    _classid = 'FOLDCHANGEMATRIX'

    def __init__(self, matrix1=None, matrix2=None, file_name=None, mode='a', tmpdir=None,
                 regions=None, scale_matrices=False, log2=True,
                 weight_column='weight', ignore_zero=False,
                 _table_name='expected_contacts'):

        ComparisonMatrix.__init__(self, matrix1=matrix1, matrix2=matrix2,
                                  comparison_function=lambda x: x[0] / x[1],
                                  file_name=file_name, mode=mode, tmpdir=tmpdir,
                                  regions=regions, scale_matrices=scale_matrices,
                                  log2=log2, weight_column=weight_column,
                                  ignore_zero=ignore_zero,
                                  _table_name=_table_name)


class DifferenceMatrix(ComparisonMatrix):
    """
    Calculate the fold-change matrix of two :class:`~kaic.data.genomic.RegionMatrixTable` objects.

    fc_ij = matrix1_ij/matrix2_ij

    :param matrix1: :class:`~kaic.data.genomic.RegionMatrixTable`, such as :class:`~kaic.data.genomic.Hic`
    :param matrix2: :class:`~kaic.data.genomic.RegionMatrixTable`, such as :class:`~kaic.data.genomic.Hic`
    :param file_name: Path to save file location
    :param mode: File mode ('r' = read-only, 'w' = (over)write, 'a' = append)
    :param tmpdir: Temporary directory
    :param regions: A region selector string, :class:`~kaic.data.genomic.GenomicRegion`, or lists thereof.
                    Will subset both matrices using these region(s) before the calculation
    :param scale_matrices: If True, will scale the matrices naively by artificially increasing the number of
                           reads in one matrix uniformly, so that it has the same number of total contacts
                           as the other matrix.
    :param log2: If True, will log2-transform the output values.
    :param weight_column: Name of the column containing the weights/values for the
                          expected value calculation. If None, this will be the default field
                          in the provided :class:`~kaic.data.genomic.RegionMatrixTable`
    """
    _classid = 'DIFFERENCEMATRIX'

    def __init__(self, matrix1=None, matrix2=None, file_name=None, mode='a', tmpdir=None,
                 regions=None, scale_matrices=False, log2=True,
                 weight_column='weight', ignore_zero=False,
                 _table_name='expected_contacts'):

        ComparisonMatrix.__init__(self, matrix1=matrix1, matrix2=matrix2,
                                  comparison_function=lambda x: x[0] - x[1],
                                  file_name=file_name, mode=mode, tmpdir=tmpdir,
                                  regions=regions, scale_matrices=scale_matrices,
                                  log2=log2, weight_column=weight_column,
                                  ignore_zero=ignore_zero,
                                  _table_name=_table_name)


class ABDomainMatrix(MatrixArchitecturalRegionFeature):
    """
    Calculate the AB domain/compartment matrix by first calculating the observed/expected matrix
    and then returning a matrix in which every entry m_ij is the correlation between row i and row j
    of the observed/expected matrix.
    You can also directly calculate the Hi-C correlation matrix by setting 'ratio' to False.

    :param hic: A :class:`~kaic.data.genomic.Hic` matrix
    :param file_name: Path to save file location
    :param mode: File mode ('r' = read-only, 'w' = (over)write, 'a' = append)
    :param tmpdir: Temporary directory
    :param regions: A region selector string, :class:`~kaic.data.genomic.GenomicRegion`, or lists thereof.
                    Will subset both matrices using these region(s) before the calculation
    :param ratio: If False, will omit the step of calculating the observed/expected matrix and return
                  a Hi-C correlation matrix directly.
    :param weight_column: Name of the column containing the weights/values for the
                          expected value calculation. If None, this will be the default field
                          in the provided :class:`~kaic.data.genomic.RegionMatrixTable`
    :param per_chromosome: If True, will only calculate the intra-chromosomal ABDomainMatrix,
                           which will save computation time and memory.
    """
    _classid = 'ABDOMAINMATRIX'

    def __init__(self, hic=None, file_name=None, mode='a', tmpdir=None, regions=None,
                 ratio=True, weight_column='weight', per_chromosome=True, _table_name='ab_domain_matrix'):
        self.region_selection = regions

        # are we retrieving an existing object?
        if isinstance(hic, string_types) and file_name is None:
            file_name = hic
            hic = None

        if hic is None and file_name is not None:
            MatrixArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                      _table_name_edges=_table_name)
        else:
            self._region_selection = regions
            if regions is None:
                regions = hic.regions
                self.region_conversion = {region.ix: region.ix for region in hic.regions}
            else:
                self.region_conversion = {region.ix: i for i, region in enumerate(hic.subset(regions))}
                regions = hic.subset(regions)
            MatrixArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                      data_fields={'correlation': t.Float32Col()}, regions=regions,
                                                      default_field='correlation', _table_name_edges=_table_name)
        self.hic = hic
        self.weight_column = weight_column
        self.ratio = ratio
        self.per_chromosome = per_chromosome

    def _calculate(self):
        oer = self.hic
        if self.ratio:
            oer = ObservedExpectedRatio(self.hic, weight_column=self.weight_column,
                                        regions=self._region_selection, per_chromosome=self.per_chromosome)

        if self.per_chromosome:
            chromosomes = self.chromosomes()
            for chromosome in chromosomes:
                m = oer[chromosome, chromosome]
                corr_m = np.corrcoef(m)

                logger.info("Chromosome {}".format(chromosome))
                with RareUpdateProgressBar(max_value=m.shape[0], silent=config.hide_progressbars) as pb:
                    for i, row_region in enumerate(m.row_regions):
                        for j, col_region in enumerate(m.col_regions):
                            if j < i:
                                continue
                            source = self.region_conversion[row_region.ix]
                            sink = self.region_conversion[col_region.ix]
                            try:
                                if np.isnan(corr_m[i, j]):
                                    continue
                                self.add_edge([source, sink, corr_m[i, j]], flush=False)
                            except IndexError:
                                pass
                        pb.update(i)
        else:
            m = oer[:]
            corr_m = np.corrcoef(m)
            with RareUpdateProgressBar(max_value=m.shape[0], silent=config.hide_progressbars) as pb:
                for i, row_region in enumerate(m.row_regions):
                    for j in range(i, len(m.row_regions)):
                        col_region = m.row_regions[j]
                        source = self.region_conversion[row_region.ix]
                        sink = self.region_conversion[col_region.ix]
                        self.add_edge([source, sink, corr_m[i, j]], flush=False)
                    pb.update(i)
        if self.ratio:
            oer.close()

        self.flush()


class ABDomains(VectorArchitecturalRegionFeature):
    """
    Calculate AB domain membership per region in a Hi-C matrix. This analysis is based on the
    sign of the first eigenvector of the :class:`~ABDomainMatrix` - a region will be assigned
    domain 'A' if its corresponding value in the eigenvector is >= 0, and 'B' otherwise.

    :param data: A :class:`~kaic.data.genomic.Hic` object or :class:`~ABDomainMatrix`
    :param file_name: Path to save file location
    :param mode: File mode ('r' = read-only, 'w' = (over)write, 'a' = append)
    :param tmpdir: Temporary directory
    :param genome: A Genome for GC content calculation, used to decide if negative
                   eigenvector means A or B
    :param regions: A region selector string, :class:`~kaic.data.genomic.GenomicRegion`, or lists thereof.
                    Will subset both matrices using these region(s) before the calculation
    :param per_chromosome: If True, will only calculate the intra-chromosomal ABDomainMatrix,
                           which will save computation time and memory.
    """
    _classid = 'ABDOMAINS'

    def __init__(self, data=None, file_name=None, mode='a', tmpdir=None, genome=None,
                 per_chromosome=True, regions=None, eigenvector=0, _table_name='ab_domains'):
        self.region_selection = regions

        # are we retrieving an existing object?
        if isinstance(data, string_types) and file_name is None:
            file_name = data
            data = None

        if data is None and file_name is not None:
            VectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                      _table_name_data=_table_name)
        else:
            if regions is None:
                regions = data.regions
            else:
                regions = data.subset(regions)

            fields = {'ev': t.Float32Col()}

            VectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                      data_fields=fields, regions=regions,
                                                      _table_name_data=_table_name)

        self.per_chromosome = per_chromosome
        self.data = data
        try:
            self._eigenvector = self.meta.eigenvector
        except AttributeError:
            self._eigenvector = eigenvector

        self.genome = genome

    def _calculate(self):
        if isinstance(self.data, Hic):
            ab_data = ABDomainMatrix(self.data, regions=self.region_selection, per_chromosome=self.per_chromosome)
            close_data = True
        else:
            close_data = False
            ab_data = self.data

        ev = np.zeros(len(self.regions))
        if self.per_chromosome:
            for chromosome in self.chromosomes():
                m = ab_data[chromosome, chromosome]
                m[np.isnan(m)] = 0
                w, v = np.linalg.eig(m)
                ab_vector = v[:, self._eigenvector]
                for i, region in enumerate(m.row_regions):
                    ev[region.ix] = ab_vector[i]
        else:
            m = ab_data[:]
            m[np.isnan(m)] = 0
            w, v = np.linalg.eig(m)
            ab_vector = v[:, self._eigenvector]
            for i, region in enumerate(m.row_regions):
                ev[region.ix] = ab_vector[i]

        if self.genome is not None:
            logger.info("Using GC content to orient eigenvector...")
            if isinstance(self.genome, string_types):
                genome = Genome.from_string(self.genome, mode='r')
            else:
                genome = self.genome

            gc_content = [np.nan] * len(self.regions)
            for chromosome in self.chromosomes():
                logger.info("{}".format(chromosome))
                chromosome_sequence = genome[chromosome].sequence
                for region in self.regions(chromosome):
                    s = chromosome_sequence[region.start - 1:region.end]
                    gc_content[region.ix] = calculate_gc_content(s)
            gc_content = np.array(gc_content)

            # use gc content to orient AB domain vector per chromosome
            cb = self.chromosome_bins
            for chromosome, (start, end) in cb.items():
                ev_sub = ev[start:end]
                gc_sub = gc_content[start:end]
                a_ixs = np.where(ev_sub >= 0.)
                b_ixs = np.where(ev_sub < 0.)
                gc_a = np.nanmean(gc_sub[a_ixs])
                gc_b = np.nanmean(gc_sub[b_ixs])

                if gc_a < gc_b:  # AB compartments are reversed!
                    ev[start:end] = -1 * ev_sub

        for region in self.regions(lazy=True):
            region.ev = ev[region.ix]
        self.flush()

        if close_data:
            ab_data.close()

    @calculateondemand
    def ab_domain_eigenvector(self, region=None):
        """
        Get the eigenvector of the :class:`~ABDomainMatrix`
        :return: list of floats
        """
        regions = self.regions(region, lazy=False)
        return [r.ev for r in regions]

    @calculateondemand
    def ab_regions(self):
        """
        Get a list of regions, each with a 'type' attribute that is either 'A' or 'B',
        depending on the sign of the matrix eigenvector.
        :return: list of regions
        """
        ev = self.ab_domain_eigenvector()

        domains = []
        current_domain = None
        last_region = None
        current_scores = []
        for i, region in enumerate(self.regions(lazy=False)):
            domain_type = 'A' if ev[i] >= 0 else 'B'

            if last_region is not None and region.chromosome != last_region.chromosome:
                current_domain = None
                current_scores = []

            if current_domain is None:
                current_domain = GenomicRegion(chromosome=region.chromosome, start=region.start, end=region.end,
                                               score=ev[i], type=domain_type, name=domain_type)
                current_scores = [ev[i]]
            else:
                if (region.ev < 0 and last_region.ev < 0) or (region.ev >= 0 and last_region.ev >= 0):
                    current_domain.end = region.end
                    current_scores.append(ev[i])
                    current_domain.score = np.nanmean(current_scores)
                else:
                    domains.append(current_domain)
                    current_domain = GenomicRegion(chromosome=region.chromosome, start=region.start, end=region.end,
                                                   score=ev[i], type=domain_type, name=domain_type)
                    current_scores = [ev[i]]
            last_region = region
        return domains


class PossibleContacts(TableArchitecturalFeature):
    """
    Calculate the possible number of intra- and inter-chromosomal contacts in a
    :class:`~kaic.data.genomic.RegionMatrixTable`. This is a combinatorial approach
    that also takes into account unmappable/masked regions in the genome.

    :param hic: :class:`~kaic.data.genomic.RegionMatrixTable`
    :param file_name: Path to save file location
    :param mode: File mode ('r' = read-only, 'w' = (over)write, 'a' = append)
    :param tmpdir: Temporary directory
    :param regions: A region selector string, :class:`~kaic.data.genomic.GenomicRegion`, or lists thereof.
                    Will subset both matrices using these region(s) before the calculation
    :param weight_column: Name of the column containing the weights/values for the
                          expected value calculation. If None, this will be the default field
                          in the provided :class:`~kaic.data.genomic.RegionMatrixTable`
    """
    _classid = 'POSSIBLECONTACTS'

    def __init__(self, hic=None, file_name=None, mode='a', tmpdir=None, regions=None,
                 weight_column='weight', _table_name='possible_contacts'):
        if isinstance(hic, string_types) and file_name is None:
            file_name = hic
            hic = None

        TableArchitecturalFeature.__init__(self, _table_name,
                                           {'intra': t.Int32Col(), 'inter': t.Int32Col()},
                                           file_name=file_name, mode=mode, tmpdir=tmpdir)

        self.hic = hic
        self.regions = regions
        self.weight_column = weight_column

    def _calculate(self):
        marginals = self.hic.marginals(weight_column=self.weight_column)

        mappable = defaultdict(int)
        if self.regions is None:
            for r in self.hic.regions(lazy=True):
                if marginals[r.ix] > 0:
                    mappable[r.chromosome] += 1
        else:
            if isinstance(self.regions, string_types) or isinstance(self.regions, GenomicRegion):
                self.regions = [self.regions]

            for region in self.regions:
                if isinstance(region, string_types):
                    region = GenomicRegion.from_string(region)

                for r in self.hic.subset(region, lazy=True):
                    if marginals[r.ix] > 0:
                        mappable[r.chromosome] += 1

        # calculate possible combinations
        intra_possible = 0
        inter_possible = 0
        chromosomes = list(mappable.keys())
        for i in range(len(chromosomes)):
            chromosome1 = chromosomes[i]
            n1 = mappable[chromosome1]
            intra_possible += n1**2/2 + n1/2
            for j in range(i+1, len(chromosomes)):
                chromosome2 = chromosomes[j]
                n2 = mappable[chromosome2]
                inter_possible += n1*n2

        self.data('intra', [intra_possible])
        self.data('inter', [inter_possible])

    @calculateondemand
    def intra_possible(self):
        """
        Get the number of theoretically possible intra-chromosomal region pars for this matrix.
        :return: int
        """
        return self[0, 'intra']

    @calculateondemand
    def inter_possible(self):
        """
        Get the number of theoretically possible inter-chromosomal region pars for this matrix.
        :return: int
        """
        return self[0, 'inter']


class RowRegionMatrix(np.ndarray):
    """
    Covenience class to add row information to a matrix.
    """
    def __new__(cls, input_matrix, regions=None, fields=None, y_values=None):
        obj = np.asarray(input_matrix).view(cls)
        obj.regions = regions
        obj.fields = fields
        obj.y_values = y_values
        obj._chromosome_index = None
        obj._region_index = None
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return

        self.regions = getattr(obj, 'regions', None)
        self.fields = getattr(obj, 'fields', None)
        self.y_values = getattr(obj, 'y_values', None)
        self._chromosome_index = None
        self._region_index = None

    def _build_region_index(self):
        self._chromosome_index = defaultdict(list)
        self._region_index = defaultdict(list)
        for i, region in enumerate(self.regions):
            self._region_index[region.chromosome].append(i)
            self._chromosome_index[region.chromosome].append(region.end)

    def region_bins(self, region):
        """
        Get a slice of indices for regions in this matrix spanned by 'region'.
        :param region: A :class:`~kaic.data.genomic.GenomicRegion` object or region selector string
        :return: slice
        """
        if self._region_index is None or self._chromosome_index is None:
            self._build_region_index()
        if isinstance(region, string_types):
            region = GenomicRegion.from_string(region)
        start_ix = bisect_left(self._chromosome_index[region.chromosome], region.start)
        end_ix = bisect_left(self._chromosome_index[region.chromosome], region.end)
        start_region_ix = self._region_index[region.chromosome][start_ix]
        end_region_ix = self._region_index[region.chromosome][end_ix]
        return slice(start_region_ix, end_region_ix + 1)

    def __getitem__(self, index):
        self._getitem = True

        # convert string types into region indexes
        if isinstance(index, tuple):
            row_key = self._convert_region_key(index[0])
            col_key = self._convert_field_key(index[1])
            index = (row_key, col_key)
        else:
            row_key = self._convert_region_key(index)
            try:
                col_key = slice(0, len(self.fields), 1)
            except TypeError:
                col_key = None
            index = row_key

        try:
            out = np.ndarray.__getitem__(self, index)
        finally:
            self._getitem = False

        if not isinstance(out, np.ndarray):
            return out

        # get regions
        try:
            row_regions = self.regions[row_key]
        except TypeError:
            row_regions = None

        if not isinstance(col_key, list):
            try:
                col_fields = self.fields[col_key]
            except TypeError:
                col_fields = None
        else:
            try:
                col_fields = [self.fields[key] for key in col_key]
            except TypeError:
                col_fields = None

        out.fields = col_fields
        out.regions = row_regions

        return out

    def __getslice__(self, start, stop):
        return self.__getitem__(slice(start, stop))

    def _convert_region_key(self, key):
        if isinstance(key, string_types):
            key = GenomicRegion.from_string(key)
        if isinstance(key, GenomicRegion):
            return self.region_bins(key)
        return key

    def _convert_field_key(self, key):
        if isinstance(key, string_types):
            return self.fields.index(key)

        if isinstance(key, list):
            if len(key) == 0:
                raise ValueError("Length of supplied list is 0.")
            l = []
            for k in key:
                if isinstance(k, string_types):
                    k = self.fields.index(k)
                l.append(k)
            return l

        return key


class MultiVectorArchitecturalRegionFeature(VectorArchitecturalRegionFeature):
    """
    Base class for numeric vector-based objects. Vectors must all be of the same type for this
    class to function properly. Provides a convenience method to convert
    multiple vectors into a matrix.

    :param file_name: Path to save file location
    :param mode: File mode ('r' = read-only, 'w' = (over)write, 'a' = append)
    :param tmpdir: Path to temporary directory
    :param data_fields: dict or class with PyTables data types
    :param regions: A region selector string, :class:`~kaic.data.genomic.GenomicRegion`, or lists thereof.
                    Will subset both matrices using these region(s) before the calculation
    :param data: dict with data to load while initializing the matrix (only here for convenience)

    """
    _classid = 'MULTIVECTORARCHITECTURALREGIONFEATURE'

    def __init__(self, file_name=None, mode='a', data_fields=None,
                 regions=None, data=None, _table_name_data='array_region_data',
                 tmpdir=None):
        VectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, data_fields=data_fields,
                                                  regions=regions, data=data, _table_name_data=_table_name_data,
                                                  tmpdir=tmpdir)
        self._y_values = None

    def _fields(self, keys=None):
        if keys is None:
            keys = self.data_field_names
        elif isinstance(keys, slice) or isinstance(keys, int):
            keys = self.data_field_names[keys]
        elif isinstance(keys, list) and len(keys) > 0 and isinstance(keys[0], int):
            keys = [self.data_field_names[i] for i in keys]
        return keys

    def as_matrix(self, regions=None, keys=None):
        """
        Get parts or all of this object's vectors concatenated into a matrix.

        :param regions: If None, will default to all regions in the matrix. Else provide an iterator
                        over :class:`~kaic.data.genomic.GenomicRegion`s
        :param keys: Keys of vectors to include in matrix. Will by default include all vectors.
        :return: :class:`~RowRegionMatrix`
        """
        if regions is None:
            region_iter = self.regions(lazy=False)
        else:
            region_iter = self.subset(regions, lazy=False)

        keys = self._fields(keys)

        array = []
        array_regions = []
        array_keys = [key for key in keys]
        for region in region_iter:
            array_regions.append(region)
            row = []
            for key in array_keys:
                row.append(getattr(region, key))
            array.append(row)
        array = np.array(array)

        if self.y_values is not None:
            y_values = []
            keyset = set(keys)
            for i, field in enumerate(self.data_field_names):
                if field in keyset:
                    y_values.append(self.y_values[i])
        else:
            y_values = np.arange(array.shape[1] + 1)

        return RowRegionMatrix(array, regions=array_regions, fields=array_keys, y_values=y_values)

    @property
    def y_values(self):
        """
        Get the y values corresponding to vector keys.
        :return: list of floats or ints
        """
        return self._y_values

    @y_values.setter
    def y_values(self, values):
        """
        Set the y values of vectors in this matrix.

        :param values: list or numpy array of floats or ints.
        """
        if len(values) != len(self.data_field_names):
            raise ValueError("Length of y-values "
                             "({}) must be the same as length of data fields ({})".format(len(values),
                                                                                          len(self.data_field_names)))
        self._y_values = values

    def _calculate(self, *args, **kwargs):
        pass


class DirectionalityIndex(MultiVectorArchitecturalRegionFeature):
    """
    Calculate the directionality index for a given set of window sizes on a Hi-C matrix.

    The directionality index (Dixon 2012 et al.) is a measure for up-/downstream biases of contact counts any
    given region displays.

    :param hic: :class:`~kaic.data.genomic.RegionMatrixTable`, typically
                a :class:`~kaic.data.genomic.Hic`object
    :param file_name: Path to save file location
    :param mode: File mode ('r' = read-only, 'w' = (over)write, 'a' = append)
    :param tmpdir: Path to temporary directory
    :param weight_column: Name of the column containing the weights/values for the
                          expected value calculation. If None, this will be the default field
                          in the provided :class:`~kaic.data.genomic.RegionMatrixTable`
    :param regions: A region selector string, :class:`~kaic.data.genomic.GenomicRegion`, or lists thereof.
                    Will subset both matrices using these region(s) before the calculation
    :param window_sizes: A list of intergers with the window sizes (in base pairs) to use for
                         directionality index calculations.
    """
    _classid = 'DIRECTIONALITYINDEX'

    def __init__(self, hic=None, file_name=None, mode='a', tmpdir=None,
                 weight_column=None, regions=None, window_sizes=(2000000,),
                 _table_name='directionality_index'):

        self.region_selection = regions

        # are we retrieving an existing object?
        if isinstance(hic, string_types) and file_name is None:
            file_name = hic
            hic = None

        if hic is None and file_name is not None:
            MultiVectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                           _table_name_data=_table_name)
        else:
            di_fields = {}
            self.window_sizes = []
            for i, window_size in enumerate(window_sizes):
                di_fields['di_%d' % window_size] = t.Float32Col(pos=i)
                self.window_sizes.append(window_size)

            if regions is None:
                regions = hic.regions
            else:
                regions = hic.subset(regions)

            MultiVectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                           data_fields=di_fields, regions=regions,
                                                           _table_name_data=_table_name)

            self.hic = hic
            if weight_column is None:
                self.weight_column = self.hic.default_field
            else:
                self.weight_column = weight_column

        self.window_sizes = []
        for colname in self._regions.colnames:
            if colname.startswith("di_"):
                window_size = int(colname[3:])
                self.window_sizes.append(window_size)
        self.y_values = self.window_sizes

    def _get_boundary_distances(self):
        n_bins = len(self.hic.regions)
        # find distances to chromosome boundaries in bins
        boundary_dist = np.zeros(n_bins, dtype=int)
        last_chromosome = None
        last_chromosome_index = 0
        for i, region in enumerate(self.hic.regions(lazy=True)):
            chromosome = region.chromosome
            if last_chromosome is not None and chromosome != last_chromosome:
                chromosome_length = i-last_chromosome_index
                for j in range(chromosome_length):
                    boundary_dist[last_chromosome_index+j] = min(j, i-last_chromosome_index-1-j)
                last_chromosome_index = i
            last_chromosome = chromosome
        chromosome_length = n_bins-last_chromosome_index
        for j in range(chromosome_length):
            boundary_dist[last_chromosome_index+j] = min(j, n_bins-last_chromosome_index-1-j)

        return boundary_dist

    def _directionality_index(self, window_size=2000000):
        bin_window_size = self.hic.distance_to_bins(window_size)

        n_bins = len(self.hic.regions)
        boundary_dist = self._get_boundary_distances()

        if self.region_selection is not None:
            edge_iter = self.hic.edge_subset((self.region_selection, self.region_selection),
                                             only_intrachromosomal=True)
        else:
            edge_iter = self.hic.edges(lazy=True, only_intrachromosomal=True)

        left_sums = np.zeros(n_bins)
        right_sums = np.zeros(n_bins)
        directionality_index = np.zeros(n_bins)
        for edge in edge_iter:
            source = edge.source
            sink = edge.sink
            weight = getattr(edge, self.weight_column)
            if source == sink:
                continue
            if sink - source <= bin_window_size:
                if boundary_dist[sink] >= sink-source:
                    left_sums[sink] += weight
                if boundary_dist[source] >= sink-source:
                    right_sums[source] += weight

        for i in range(n_bins):
            A = left_sums[i]
            B = right_sums[i]
            E = (A+B)/2
            if E != 0 and B-A != 0:
                directionality_index[i] = ((B-A)/abs(B-A)) * ((((A-E)**2)/E) + (((B-E)**2)/E))

        if self.region_selection is not None:
            nodes_ix = self.hic._getitem_nodes(key=self.region_selection, as_index=True)

            if not isinstance(nodes_ix, list):
                nodes_ix = [nodes_ix]

            di_sub = []
            for node_range in ranges(nodes_ix):
                di_sub += list(directionality_index[node_range[0]:node_range[1]+1])

            return np.array(di_sub)
        return directionality_index

    def _calculate(self):
        for window_size in self.window_sizes:
            logger.info("Calculating directionality index for window size {}".format(window_size))
            directionality_index = self._directionality_index(window_size)
            self.data("di_%d" % window_size, directionality_index)

    @calculateondemand
    def directionality_index(self, window_size=None):
        """
        Get the directionality index for the given window size.

        :param window_size: Window size in base pairs. If None, will select the first window size
                            used in the initialization of this object
        :return: list of floats
        """
        if window_size is None:
            window_size = self.window_sizes[0]
        return self[:, 'di_%d' % window_size]


class InsulationIndex(MultiVectorArchitecturalRegionFeature):
    """
    Calculate the insulation index for a list of window sizes.

    The insulation index (Crane et al. 2015) is a measure of the number of interactions that cross a region
    in the genome. If the value is very low, a region effectively acts as a "contact boundary", while
    highl values may, for example, reflect that a region is part of a highly self-interacting region.

    This class offers a number of different calculation and normalisation options. While the default
    should work well in many cases, it is highly recommended that you read through the options listed below
    to get the most out of your analysis.

    :param hic: :class:`~kaic.data.genomic.RegionMatrixTable`, typically
                a :class:`~kaic.data.genomic.Hic`object
    :param file_name: Path to save file location
    :param mode: File mode ('r' = read-only, 'w' = (over)write, 'a' = append)
    :param tmpdir: Path to temporary directory
    :param regions: A region selector string, :class:`~kaic.data.genomic.GenomicRegion`, or lists thereof.
                    Will subset both matrices using these region(s) before the calculation
    :param offset: Offset of the insulation square from the diagonal. Can be useful to avoid biases
                   stemming from very bright diagonals
    :param normalise: If True, will normalise the insulation values by the chromosome average. You can change
                      this behavior to division by a smaller sliding window average by specifiying a number of
                      bins using the _normalisation_window argument
    :param impute_missing: If True, will do a very simplistic missing value imputation by replacing missing
                           values with the expected intra-chromosomal value given average contact values at
                           this region separation. This may be useful if too many regins in the genome are
                           missing or you observe artefacts.
    :param window_sizes: List of window sizes in base pairs to calculate the insulation index
    :param log: If True, log2-transform insulation values. Particularly useful in combination with the
                'normalise' option.
    :param subtract_mean: Instead of dividing by the mean when 'normalise' is True, subtract the mean.
                          This is useful when calculating the insulation index on already-log-transformed values,
                          such as in a :class:`~FoldChangeMatrix`
    """
    _classid = 'INSULATIONINDEX'

    def __init__(self, hic=None, file_name=None, mode='a', tmpdir=None,
                 regions=None, offset=0, normalise=False, impute_missing=False,
                 window_sizes=(200000,), log=False, _normalisation_window=300,
                 subtract_mean=False, trim_mean_proportion=0.0, geometric_mean=False,
                 _table_name='insulation_index'):
        self.region_selection = regions

        # are we retrieving an existing object?
        if isinstance(hic, string_types) and file_name is None:
            file_name = hic
            hic = None

        if hic is None and file_name is not None:
            MultiVectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                           _table_name_data=_table_name)
        else:
            if regions is None:
                regions = hic.regions
            else:
                regions = hic.subset(regions)

            ii_fields = {}
            self.window_sizes = []
            for i, window_size in enumerate(window_sizes):
                ii_fields['ii_%d' % window_size] = t.Float32Col(pos=i)
                self.window_sizes.append(window_size)

            MultiVectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                           data_fields=ii_fields, regions=regions,
                                                           _table_name_data=_table_name)

        self.window_sizes = []
        for colname in self._regions.colnames:
            if colname.startswith("ii_"):
                window_size = int(colname[3:])
                self.window_sizes.append(window_size)
        self.y_values = self.window_sizes

        self.offset = offset
        self.hic = hic
        self.impute_missing = impute_missing
        self.normalise = normalise
        self.normalisation_window = _normalisation_window
        self.subtract_mean = subtract_mean
        self.trim_mean_proportion = trim_mean_proportion
        self.log = log
        if geometric_mean:
            self._stat = gmean
        else:
            self._stat = np.nanmean

    def _insulation_index_lowmem(self, window_size, window_offset=0, aggr_func=None, _mappable=None,
                                 _na_threshold=0.5, _expected=None):
        lowmem = False
        if aggr_func is None:
            lowmem = True
            aggr_func = np.ma.mean

        chromosome_bins = self.hic.chromosome_bins
        window_offset += 1

        if _mappable is None:
            _mappable = self.hic.mappable()

        if self.impute_missing and _expected is None:
            _expected = ExpectedContacts(self.hic, smooth=True)

        if _expected is not None:
            intra_expected = _expected.intra_expected()
        else:
            intra_expected = None

        def _pair_to_bins(source, sink):
            # correct index by chromosome and window offset
            i = source - chromosome_start + window_offset
            j = sink - chromosome_start - window_offset

            start = max(i, j - window_size + 1)
            stop = min(j + 1, i + window_size)
            for ii_bin in range(start, stop):
                yield ii_bin

        ii_list = []
        chromosomes = self.hic.chromosomes()
        for chromosome in chromosomes:
            chromosome_start, chromosome_stop = chromosome_bins[chromosome]
            if lowmem:
                values_by_chromosome = [0 for _ in range(chromosome_start, chromosome_stop)]
            else:
                values_by_chromosome = [list() for _ in range(chromosome_start, chromosome_stop)]

            # add each edge weight to every insulation window that contains it
            for edge in self.hic.edge_subset((chromosome, chromosome), lazy=True):
                for ii_bin in _pair_to_bins(edge.source, edge.sink):
                    weight = getattr(edge, self.hic.default_field)
                    if lowmem:
                        values_by_chromosome[ii_bin] += weight
                    else:
                        values_by_chromosome[ii_bin].append(weight)

            # add imputed values, if requested
            if intra_expected is not None:
                covered = set()
                for ix in range(chromosome_start, chromosome_stop):
                    if not _mappable[ix]:
                        sink = ix
                        for source in range(max(chromosome_start, ix - window_size - window_offset - 2), ix):
                            if (source, sink) in covered:
                                continue
                            covered.add((source, sink))
                            weight = intra_expected[sink - source]
                            for ii_bin in _pair_to_bins(source, sink):
                                if lowmem:
                                    values_by_chromosome[ii_bin] += weight
                                else:
                                    values_by_chromosome[ii_bin].append(weight)

                        source = ix
                        for sink in range(ix, min(chromosome_stop, ix + window_offset + window_size + 2)):
                            if (source, sink) in covered:
                                continue
                            covered.add((source, sink))
                            weight = intra_expected[sink - source]
                            for ii_bin in _pair_to_bins(source, sink):
                                if lowmem:
                                    values_by_chromosome[ii_bin] += weight
                                else:
                                    values_by_chromosome[ii_bin].append(weight)

                for k in range(len(values_by_chromosome)):
                    if (k - window_offset < window_size - 1
                            or k + window_offset > len(values_by_chromosome) - window_size):
                        if lowmem:
                            values_by_chromosome[k] = np.nan
                        else:
                            values_by_chromosome[k] = []
                    else:
                        if not lowmem:
                            for _ in range(len(values_by_chromosome[k]), window_size**2):
                                values_by_chromosome[k].append(0)
                        else:
                            values_by_chromosome[k] /= window_size**2
            # count unmappable bins in every window
            else:
                unmappable_horizontal = [0 for _ in range(chromosome_start, chromosome_stop)]
                unmappable_vertical = [0 for _ in range(chromosome_start, chromosome_stop)]

                for ix in range(chromosome_start, chromosome_stop):
                    if not _mappable[ix]:
                        # horizontal
                        start_bin = ix - chromosome_start + window_offset
                        for ii_bin in range(start_bin, start_bin + window_size):
                            if 0 <= ii_bin < len(unmappable_horizontal):
                                unmappable_horizontal[ii_bin] += 1
                        # vertical
                        start_bin = ix - chromosome_start - window_offset
                        for ii_bin in range(start_bin - window_size + 1, start_bin + 1):
                            if 0 <= ii_bin < len(unmappable_vertical):
                                unmappable_vertical[ii_bin] += 1

                for k in range(len(values_by_chromosome)):
                    values = values_by_chromosome[k]
                    na_vertical = unmappable_vertical[k] * window_size
                    na_horizontal = unmappable_horizontal[k] * window_size
                    na_overlap = unmappable_horizontal[k] * unmappable_vertical[k]
                    na_total = na_vertical + na_horizontal - na_overlap

                    # take into account nan values when adding zeros
                    if ((na_total > (window_size**2 * _na_threshold))
                            or k - window_offset < window_size - 1
                            or k + window_offset > len(values_by_chromosome) - window_size):
                        if lowmem:
                            values_by_chromosome[k] = np.nan
                        else:
                            values_by_chromosome[k] = []
                    else:
                        if not lowmem:
                            for _ in range(len(values), window_size**2 - na_total):
                                values_by_chromosome[k].append(0)
                        else:
                            values_by_chromosome[k] /= (window_size**2 - na_total)

            if lowmem:
                ii_by_chromosome = values_by_chromosome
            else:
                ii_by_chromosome = []
                for values in values_by_chromosome:
                    ii_by_chromosome.append(np.nan) if len(values) == 0 else ii_by_chromosome.append(aggr_func(values))

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)

                ii_by_chromosome = np.array(ii_by_chromosome)
                if self.normalise:
                    logger.debug("Normalising insulation index")
                    if self.normalisation_window is not None:
                        logger.debug("Sliding window average")
                        mean_ins = apply_sliding_func(ii_by_chromosome, self.normalisation_window,
                                                      func=lambda x: trim_stats(x[np.isfinite(x)],
                                                                                self.trim_mean_proportion,
                                                                                stat=self._stat))
                    else:
                        logger.debug("Whole chromosome mean")
                        mean_ins = trim_stats(ii_by_chromosome[np.isfinite(ii_by_chromosome)],
                                              self.trim_mean_proportion, stat=self._stat)

                    if not self.subtract_mean:
                        logger.debug("Dividing by mean")
                        ii_by_chromosome = ii_by_chromosome / mean_ins
                    else:
                        logger.debug("Subtracting mean")
                        ii_by_chromosome = ii_by_chromosome - mean_ins
            ii_list.append(ii_by_chromosome)

        ins_matrix = np.array(list(itertools.chain.from_iterable(ii_list)))

        if self.region_selection is not None:
            ins_matrix = ins_matrix[self.hic.region_bins(self.region_selection)]

        if self.log:
            logger.debug("Log-transforming insulation index")
            return np.log2(ins_matrix)
        else:
            return ins_matrix

    def _calculate(self):
        mappable = self.hic.mappable()
        if self.impute_missing:
            ex = ExpectedContacts(self.hic, smooth=True)
        else:
            ex = None
        offset_bins = self.hic.distance_to_bins(self.offset)

        for window_size in self.window_sizes:
            logger.info("Calculating insulation index for window size {}".format(window_size))
            bins = self.hic.distance_to_bins(window_size)
            insulation_index = self._insulation_index_lowmem(bins, offset_bins,
                                                             _mappable=mappable, _expected=ex)

            self.data("ii_%d" % window_size, insulation_index)
        if ex is not None:
            ex.close()

    @calculateondemand
    def insulation_index(self, window_size, nan_window=None, nan_threshold=0.33):
        """
        Get the insulation vector for the given window size.

        It is possible to filter insulation values on the fly that are in areas
        of high NaN-density (generally highly repetitive regions) using the
        'nan_' parameters. If nan_window is not None, and more than a fraction
        of nan_treshold values in nan_window upstream OR downstream of a region
        is NaN, the region's insulation value will also be set to NaN.

        :param window_size: Window size in base pairs
        :param nan_window: Window size in bins to filter insulation values based
                           on the NaN density in surrounding regions (see above)
        :param nan_threshold: Fraction of allowed NaNs in nan_window up- or downstream
                              of a region
        :return: list of insulation values for each genomic region
        """
        if window_size is None:
            window_size = self.window_sizes[0]

        if isinstance(window_size, int):
            window_size = 'ii_%d' % window_size

        ii = self[:, window_size]
        # filter vector on neighboring nan regions
        if nan_window is not None:
            ii_nan = np.array(ii).copy()
            for i, value in enumerate(ii):
                # if already nan, no need to calculate
                if np.isnan(value):
                    continue

                left = max(0, i - nan_window)
                right = min(len(ii), i + nan_window + 1)
                ii_left = ii[left:i]
                ii_right = ii[i+1:right]

                nans_left = np.sum(np.isnan(ii_left)) + (nan_window - len(ii_left))
                nans_right = np.sum(np.isnan(ii_right)) + (nan_window - len(ii_right))

                if nans_left/nan_window > nan_threshold or nans_right/nan_window > nan_threshold:
                    ii_nan[i] = np.nan
            ii = list(ii_nan)

        return ii

    def boundaries(self, window_size, min_score=None, delta_window=3, log=False, sub_bin_precision=False,
                   call_maxima=False):
        """
        Call insulation boundaries based on minima in an insulation vector of this object.

        :param window_size: Window size in base pairs to select insulation vector
        :param min_score: Minimum difference between minimum and the closest maximum
                          in the insulation vector for a region to be considered a
                          boundary
        :param delta_window: Window size in bins to control smoothing of the delta function
                             used to calculate the derivative of the insulation index.
                             Calculation takes into account d bins upstream and d
                             bins downstream for a total window size of 2*d + 1 bins.
        :param log: Log2-transform insulation index before boundary calls
        :param sub_bin_precision: Call boundaries with sub bin precision, by taking
                                  into account the precise zero transition of the delta vector.
        :param call_maxima: Call maxima instead of minima as boundaries
        :return: list of :class:`~kaic.data.genomic.GenomicRegion`
        """
        index = self.insulation_index(window_size)
        if log:
            index = np.log2(index)
        peaks = MaximaCallerDelta(index, window_size=delta_window, sub_bin_precision=sub_bin_precision)
        if call_maxima:
            minima, scores = peaks.get_maxima()
        else:
            minima, scores = peaks.get_minima()
        regions = list(self.regions)

        boundaries = []
        for i, ix in enumerate(minima):
            if min_score is not None and scores[i] < min_score:
                continue
            if sub_bin_precision:
                region = regions[int(ix)]
                frac = ix % 1
                shift = int((frac - .5)*(region.end - region.start))
                b = GenomicRegion(chromosome=region.chromosome, start=region.start + shift, end=region.end + shift, score=scores[i])
            else:
                region = regions[ix]
                b = GenomicRegion(chromosome=region.chromosome, start=region.start, end=region.end, score=scores[i])
            boundaries.append(b)

        logger.info("Found {} boundaries for window size {}".format(len(boundaries), window_size))
        return boundaries


class RegionContactAverage(MultiVectorArchitecturalRegionFeature):
    """
    Calculates the average number of contacts for a given region within a certain window.

    The size and shape of the window can be controlled with the below parameters. In principle,
    the window is very similar to the insulation index window, but has added padding (directionality
    index window has a width of 1, region contact average can be as wide as required) and can be offset
    from the diagonal.

    :param matrix: :class:`~kaic.data.genomic.RegionMatrixTable`, typically
                   a :class:`~kaic.data.genomic.Hic`object
    :param file_name: Path to save file location
    :param mode: File mode ('r' = read-only, 'w' = (over)write, 'a' = append)
    :param tmpdir: Path to temporary directory
    :param window_sizes: List of window sizes in base pairs to calculate the insulation index
    :param regions: A region selector string, :class:`~kaic.data.genomic.GenomicRegion`, or lists thereof.
                    Will subset both matrices using these region(s) before the calculation
    :param offset: Offset of window from the diagonal of the matrix
    :param padding: Padding of window (width = 1 + 2*padding)
    :param impute_missing: If True, will do a very simplistic missing value imputation by replacing missing
                           values with the expected intra-chromosomal value given average contact values at
                           this region separation. This may be useful if too many regins in the genome are
                           missing or you observe artefacts.
    """
    _classid = 'REGIONCONTACTAVERAGE'

    def __init__(self, matrix=None, file_name=None, mode='a', tmpdir=None,
                 window_sizes=(200000,), regions=None, offset=0, padding=1, impute_missing=True,
                 _table_name='contact_average'):
        self.region_selection = regions

        # are we retrieving an existing object?
        if isinstance(matrix, string_types) and file_name is None:
            file_name = matrix
            matrix = None

        if matrix is None and file_name is not None:
            MultiVectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                           _table_name_data=_table_name)
        else:
            if regions is None:
                regions = matrix.regions
                self.region_conversion = {region.ix: region.ix for region in matrix.regions}
            else:
                self.region_conversion = {region.ix: i for i, region in enumerate(matrix.subset(regions))}
                regions = matrix.subset(regions)

            av_fields = {}
            self.window_sizes = []
            n = 0
            for i, window_size in enumerate(window_sizes):
                av_fields['av_%d' % window_size] = t.Float32Col(pos=n)
                av_fields['avl_%d' % window_size] = t.Float32Col(pos=n+1)
                av_fields['avr_%d' % window_size] = t.Float32Col(pos=n+2)
                n += 3
                self.window_sizes.append(window_size)

            MultiVectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                           data_fields=av_fields, regions=regions,
                                                           _table_name_data=_table_name)

        self.window_sizes = []
        for colname in self._regions.colnames:
            if colname.startswith("av_"):
                window_size = int(colname[3:])
                self.window_sizes.append(window_size)

        # self.y_values = self.window_sizes
        self.offset = offset
        self.padding = padding
        self.matrix = matrix
        self.impute_missing = impute_missing

    def _contact_average(self, window_size, offset=0, padding=1, _aggr_func=np.ma.mean):
        av_values = dict()
        for chromosome in self.chromosomes():
            matrix = self.matrix.as_matrix(key=(chromosome, chromosome), mask_missing=True,
                                           impute_missing=self.impute_missing)

            # region index
            for i, region in enumerate(matrix.row_regions):
                ix = self.region_conversion[region.ix]
                slice_left = slice(max(0, i-offset-window_size), max(0, i-offset+1))
                slice_right = slice(min(i+offset, matrix.shape[0]), min(i+offset+window_size+1, matrix.shape[0]))
                slice_vertical = slice(max(0, i-padding), min(i+padding, matrix.shape[0]))

                value_left = _aggr_func(matrix[slice_vertical, slice_left])
                value_right = _aggr_func(matrix[slice_vertical, slice_right])

                av_values[ix] = (value_left, value_right)

        av_left = np.zeros(len(self.regions))
        av_right = np.zeros(len(self.regions))
        for region in self.regions(lazy=True):
            if region.ix in av_values:
                av_left[region.ix], av_right[region.ix] = av_values[region.ix]
        return av_left, av_right

    def _calculate(self):
        offset_bins = self.matrix.distance_to_bins(self.offset)
        for window_size in self.window_sizes:
            bins = self.matrix.distance_to_bins(window_size)
            av_values_left, av_values_right = self._contact_average(bins, offset_bins)
            av_values = (av_values_left + av_values_right)/2
            self.data("av_%d" % window_size, av_values)
            self.data("avl_%d" % window_size, av_values_left)
            self.data("avr_%d" % window_size, av_values_right)

    @calculateondemand
    def average_contacts(self, window_size):
        """
        Get the list of region contact averages for the given window size.

        :param window_size: Window size in base pairs.
        :return: list of floats
        """
        if window_size is None:
            window_size = self.window_sizes[0]
        return self[:, 'av_%d' % window_size]


class VectorDifference(MultiVectorArchitecturalRegionFeature):
    """
    Calculate the difference between two :class:`~MultiVectorArchitecturalRegionFeature` objects.

    Will substract any two vectors in both objects that share the same name.

    :param vector1: :class:`~MultiVectorArchitecturalRegionFeature`
    :param vector2: :class:`~MultiVectorArchitecturalRegionFeature`
    :param absolute: Transform differences into absolute values
    :param file_name: Path to save file location
    :param mode: File mode ('r' = read-only, 'w' = (over)write, 'a' = append)
    :param tmpdir: Path to temporary directory
    """
    _classid = 'VECTORDIFF'

    def __init__(self, vector1=None, vector2=None, absolute=False, file_name=None, mode='a', tmpdir=None,
                 _table_name='vector_diff'):
        if isinstance(vector1, string_types) and file_name is None and vector2 is None:
            file_name = vector1
            vector1 = None

        if vector1 is None and file_name is not None:
            MultiVectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                           _table_name_data=_table_name)
        else:
            if vector1 is None or vector2 is None:
                raise ValueError("Must provide both vector1 and vector2!")

            self.vector1 = vector1
            self.vector2 = vector2

            vector2_values = set(vector2.y_values)
            diff_fields = {}
            self.external_fields = []
            n = 0
            for i, field in enumerate(vector1.y_values):
                if field in vector2_values:
                    diff_fields['diff_{}'.format(field)] = t.Float32Col(pos=n)
                    n += 1
                    self.external_fields.append(vector1.data_field_names[i])

            MultiVectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                           data_fields=diff_fields, regions=vector1.regions,
                                                           _table_name_data=_table_name)

        self.absolute = absolute

        self.window_sizes = []
        for colname in self._regions.colnames:
            if colname.startswith("diff_"):
                window_size = int(colname[5:])
                self.window_sizes.append(window_size)

        self.y_values = self.window_sizes

    def _calculate(self, *args, **kwargs):
        for i, field in enumerate(self.external_fields):
            logger.info("Calculating difference for {}".format(field))
            v1 = np.array(self.vector1[:, field])
            v2 = np.array(self.vector2[:, field])
            d = v1-v2
            if self.absolute:
                d = abs(d)
            self.data("diff_{}".format(self.y_values[i]), d)

    @calculateondemand
    def difference(self, window_size=None):
        """
        Get the vector difference for the given window size
        :param window_size: Window size in base pairs
        :return: list of floats
        """
        if window_size is None:
            window_size = self.window_sizes[0]
        return self[:, 'diff_{}'.format(window_size)]


class MetaMatrixBase(ArchitecturalFeature, FileGroup):
    """
    Meta class for the extraction of submatrices from a matrix using a list of regions.

    :param array: :class:`~MultiVectorArchitecturalRegionFeature`
    :param regions: A region selector string, :class:`~kaic.data.genomic.GenomicRegion`, or lists thereof.
                    Will subset both matrices using these region(s) before the calculation
    :param window_width: Width of the extracted sub-matrix in bins
    :param data_selection: Names or indexes of the vectors to extract submatrix from.
    :param file_name: Path to save file location
    :param mode: File mode ('r' = read-only, 'w' = (over)write, 'a' = append)
    :param tmpdir: Path to temporary directory
    :param orient_strand: flip submatrix vertically, if region is on reverse strand
    """
    _classid = 'METAMATRIXBASE'

    def __init__(self, array=None, regions=None, window_width=50, data_selection=None,
                 file_name=None, mode='a', tmpdir=None, orient_strand=False,
                 _group_name='meta_base'):

        ArchitecturalFeature.__init__(self)

        if isinstance(array, string_types) and file_name is None:
            file_name = array
            array = None

        FileGroup.__init__(self, _group_name, file_name, mode=mode, tmpdir=tmpdir)

        try:
            self.meta_matrix = self.file.get_node(self._group, 'meta_matrix')
            self._calculated = True
        except t.NoSuchNodeError:
            self.meta_matrix = None

        self.array = array
        self.window_width = self.array.distance_to_bins(window_width)
        self._data_selection = None
        self._matrix_shape = None
        self.data_selection = data_selection
        self.regions = regions
        self.orient_strand = orient_strand

    @property
    def data_selection(self):
        """
        Get the name of the selected data vector.
        """
        return self._data_selection

    @data_selection.setter
    def data_selection(self, selection):
        matrix_shape = [self.window_width * 2 + 1, 0]
        if selection is None:
            selection = range(0, len(self.array.data_field_names))

        if isinstance(selection, range):
            selection = list(selection)

        if isinstance(selection, int) or isinstance(selection, string_types):
            selection = [selection]

        if isinstance(selection, list) or isinstance(selection, tuple):
            new_list = []
            for i in selection:
                ix = i
                if isinstance(i, string_types):
                    ix = self.array.data_field_names.index(i)
                if ix <= len(self.array.data_field_names):
                    new_list.append(ix)

            matrix_shape[1] = len(new_list)
            self._data_selection = new_list
        else:
            raise ValueError("Unsupported data_selection type({})".format(type(selection)))

        self._matrix_shape = tuple(matrix_shape)

    def _sub_matrices(self):
        chromosome_regions = defaultdict(list)
        chromosomes = self.array.chromosomes()
        chromosome_set = set(chromosomes)

        region_counter = 0
        for i, region in enumerate(self.regions):
            if region.chromosome in chromosome_set:
                midpoint = (region.start + region.end) / 2
                chromosome_regions[region.chromosome].append((midpoint, region, i))
                region_counter += 1

        ds = self.data_selection
        try:
            ds = to_slice(ds)
        except ValueError:
            pass

        with RareUpdateProgressBar(max_value=region_counter, silent=config.hide_progressbars) as pb:
            counter = 0
            for chromosome in chromosomes:
                matrix = self.array.as_matrix(chromosome)
                for pos, region, i in chromosome_regions[chromosome]:
                    counter += 1
                    pb.update(counter)
                    try:
                        bin_range = matrix.region_bins(GenomicRegion(start=pos, end=pos, chromosome=chromosome))
                    except IndexError:
                        logger.error("Cannot find bin range for {}:{}".format(chromosome, pos))
                        sub_matrix = np.zeros((2*self.window_width+1, ds))
                        sub_matrix[sub_matrix == 0] = np.nan
                        yield i, region, sub_matrix
                        continue

                    for region_ix in range(bin_range.start, bin_range.stop):
                        sub_matrix = matrix[region_ix - self.window_width:region_ix + self.window_width + 1, ds]
                        if self.orient_strand and hasattr(region, 'strand') and region.is_reverse():
                            sub_matrix = np.fliplr(sub_matrix)
                        yield i, region, sub_matrix

    @calculateondemand
    def _calculate(self, *args, **kwargs):
        pass


class MetaArray(MetaMatrixBase):
    """
    Calculate an average array matrix from a list of regions. This can be used, for example, to
    create a pile-up flame plot of insulation values at regions of interest, such as all boundaries in the
    genome.

    :param array: :class:`~MultiVectorArchitecturalRegionFeature`
    :param regions: A region selector string, :class:`~kaic.data.genomic.GenomicRegion`, or lists thereof.
                    Will subset both matrices using these region(s) before the calculation
    :param window_width: Width of the extracted sub-matrix in bins
    :param data_selection: Names or indexes of the vectors to extract submatrix from.
    :param file_name: Path to save file location
    :param mode: File mode ('r' = read-only, 'w' = (over)write, 'a' = append)
    :param tmpdir: Path to temporary directory
    :param orient_strand: flip submatrix vertically, if region is on reverse strand
    """
    _classid = 'METAARRAY'

    def __init__(self, array=None, regions=None, window_width=50000, data_selection=None,
                 file_name=None, mode='a', tmpdir=None, orient_strand=False,
                 _group_name='meta_matrix'):
        MetaMatrixBase.__init__(self, array=array, regions=regions, window_width=window_width,
                                data_selection=data_selection, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                orient_strand=orient_strand, _group_name=_group_name)

    def _calculate(self):
        avg_matrix = np.zeros(self._matrix_shape)
        count_matrix = np.zeros(self._matrix_shape)
        for _, _, m_sub in self._sub_matrices():
            if m_sub.shape == self._matrix_shape:
                count_matrix += np.isnan(m_sub) is False
                m_sub[np.isnan(m_sub)] = 0
                avg_matrix += m_sub

        avg_matrix /= count_matrix

        self.meta_matrix = self.file.create_carray(self._group, 'meta_matrix', t.Float32Atom(),
                                                   tuple(reversed(self._matrix_shape)))
        self.meta_matrix[:] = avg_matrix.T

    @calculateondemand
    def matrix(self):
        """
        Get the meta matrix
        :return: numpy array
        """
        return self.meta_matrix[:]

    def x(self):
        """
        A list of (relative) x values for the meta matrix
        :return: list of ints
        """
        x = []
        for i in range(-1*self.window_width, self.window_width + 1):
            d = self.array.bins_to_distance(i)
            x.append(d)
        return x

    def y(self):
        """
        Get a list of y-values (such as window sizes) for this meta matrix.
        :return:
        """
        y = []
        for i in self.data_selection:
            y.append(self.array.y_values[i])
        return y


class MetaHeatmap(MetaMatrixBase):
    """
    Extract sub-rows from array by a list of regions and concatenate into a heatmap array.

    :param array: :class:`~MultiVectorArchitecturalRegionFeature`
    :param regions: A region selector string, :class:`~kaic.data.genomic.GenomicRegion`, or lists thereof.
                    Will subset both matrices using these region(s) before the calculation
    :param window_width: Width of the extracted sub-matrix in bins
    :param data_selection: Names or indexes of the vectors to extract submatrix from.
    :param file_name: Path to save file location
    :param mode: File mode ('r' = read-only, 'w' = (over)write, 'a' = append)
    :param tmpdir: Path to temporary directory
    :param orient_strand: flip submatrix vertically, if region is on reverse strand
    """
    _classid = 'METAHEATMAP'

    def __init__(self, array=None, regions=None, window_width=50000, data_selection=None,
                 file_name=None, mode='a', tmpdir=None, orient_strand=False,
                 _group_name='meta_heatmap'):
        if data_selection is None and array is not None:
            data_selection = array.data_field_names[0]

        if not isinstance(data_selection, string_types) and not isinstance(data_selection, int):
            raise ValueError("data_selection parameter must be "
                             "int, str, or None, but is {}".format(type(data_selection)))

        if isinstance(array, string_types) and file_name is None:
            file_name = array
            array = None

        if file_name is not None and array is None:
            MetaMatrixBase.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                    _group_name=_group_name, orient_strand=orient_strand)
        else:
            MetaMatrixBase.__init__(self, array=array, regions=regions, window_width=window_width,
                                    data_selection=data_selection, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                    _group_name=_group_name, orient_strand=orient_strand)
            self.regions_subset = []

    def _calculate(self):
        order = []
        heatmap = []
        regions = []
        for i, region, m_sub in self._sub_matrices():
            order.append(i)
            regions.append(region)
            if m_sub.shape == self._matrix_shape:
                heatmap.append(m_sub[:, 0])
            else:
                m = np.zeros((self._matrix_shape[0],))
                m[m == 0] = np.nan
                heatmap.append(m)

        order = np.array(order).argsort()
        heatmap = np.array(heatmap)[order]
        self.regions_subset = [regions[i] for i in order]

        self.meta_matrix = self.file.create_carray(self._group, 'meta_matrix', t.Float32Atom(),
                                                   heatmap.shape)
        self.meta_matrix[:] = heatmap

    @calculateondemand
    def heatmap(self):
        """
        Get array heatmap.
        :return: numpy array
        """
        return self.meta_matrix[:]

    def x(self):
        """
        Get a list of (relative) x values for this meta matrix.
        :return: list of ints
        """
        x = []
        for i in range(-1*self.window_width, self.window_width + 1):
            d = self.array.bins_to_distance(i)
            x.append(d)
        return x


class MetaRegionAverage(MetaMatrixBase):
    """
    Calculate an average profile of array values for a list of regions.

    :param array: :class:`~MultiVectorArchitecturalRegionFeature`
    :param regions: A region selector string, :class:`~kaic.data.genomic.GenomicRegion`, or lists thereof.
                    Will subset both matrices using these region(s) before the calculation
    :param window_width: Width of the extracted sub-matrix in bins
    :param data_selection: Names or indexes of the vectors to extract submatrix from.
    :param file_name: Path to save file location
    :param mode: File mode ('r' = read-only, 'w' = (over)write, 'a' = append)
    :param tmpdir: Path to temporary directory
    :param orient_strand: flip submatrix vertically, if region is on reverse strand
    """
    _classid = 'METAREGIONAVG'

    def __init__(self, array=None, regions=None, window_width=50000, data_selection=None,
                 file_name=None, mode='a', tmpdir=None, orient_strand=False,
                 _group_name='meta_region_avg'):

        if data_selection is None and array is not None:
            data_selection = array.data_field_names[0]

        if not isinstance(data_selection, string_types) and not isinstance(data_selection, int):
            raise ValueError("data_selection parameter must be "
                             "int, str, or None, but is {}".format(type(data_selection)))

        if isinstance(array, string_types) and file_name is None:
            file_name = array
            array = None

        if file_name is not None and array is None:
            MetaMatrixBase.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                    _group_name=_group_name, orient_strand=orient_strand)
        else:
            MetaMatrixBase.__init__(self, array=array, regions=regions, window_width=window_width,
                                    data_selection=data_selection, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                    orient_strand=orient_strand, _group_name=_group_name)

    def _calculate(self):
        averages = []
        for i, region, m_sub in self._sub_matrices():
            if m_sub.shape == self._matrix_shape:
                v = m_sub[:, 0]
                averages.append(np.nansum(v)/np.sum(~np.isnan(v)))
        averages = np.array(averages)

        self.meta_matrix = self.file.create_carray(self._group, 'meta_matrix', t.Float32Atom(),
                                                   averages.shape)
        self.meta_matrix[:] = averages

    @calculateondemand
    def averages(self):
        """
        Get average profiles.
        :return: numpy array
        """
        return self.meta_matrix[:]


def cumulative_matrix(hic, regions, window, cache_matrix=False, norm=False, silent=False,
                      _mappable=None):
    """
    Construct a matrix by superimposing subsets of the Hi-C matrix from different regions.

    For each region, a matrix of a certain window size (region at the center) will be
    extracted. All matrices will be added and divided by the number of regions.

    :param hic: Hic object (or any RegionMatrix compatible)
    :param regions: Iterable with :class:`GenomicRegion` objects
    :param window: Window size in base pairs around each region
    :param cache_matrix: Load chromosome matrices into memory to speed up matrix generation
    :param norm: If True, will normalise the averaged maps to the expected map (per chromosome)
    :param silent: Suppress progress bar output
    :return: numpy array
    """
    bins = window/hic.bin_size
    bins_half = max(1, int(bins/2))
    shape = bins_half * 2 + 1

    chromosome_bins = hic.chromosome_bins
    regions_by_chromosome = defaultdict(list)
    region_counter = 0
    for region in regions:
        regions_by_chromosome[region.chromosome].append(region)
        region_counter += 1

    cmatrix = np.zeros((shape, shape))
    counter = 0
    counter_matrix = np.zeros((shape, shape))
    with RareUpdateProgressBar(max_value=len(regions), silent=silent) as pb:
        i = 0
        for chromosome, chromosome_regions in regions_by_chromosome.items():
            if chromosome not in chromosome_bins:
                continue

            matrix_expected = np.ones((shape, shape))
            if norm:
                with ExpectedContacts(hic, regions=chromosome) as ex:
                    intra_expected = ex.intra_expected()

                    for j in range(shape):
                        matrix_expected[kth_diag_indices(shape, j)] = intra_expected[j]
                        matrix_expected[kth_diag_indices(shape, -1 * j)] = intra_expected[j]

            if cache_matrix:
                chromosome_matrix = hic.as_matrix((chromosome, chromosome), mask_missing=True,
                                                  _mappable=_mappable)
                offset = chromosome_bins[chromosome][0]
            else:
                chromosome_matrix = hic
                offset = 0

            for region in chromosome_regions:
                i += 1

                center_region = GenomicRegion(chromosome=chromosome, start=region.center, end=region.center)
                center_bin = list(hic.subset(center_region))[0].ix
                if center_bin - bins_half < chromosome_bins[chromosome][0] \
                        or chromosome_bins[chromosome][1] <= center_bin + bins_half + 1:
                    continue
                center_bin -= offset

                s = slice(center_bin-bins_half, center_bin+bins_half+1)
                if cache_matrix:
                    matrix = chromosome_matrix[s, s]
                else:
                    matrix = chromosome_matrix.as_matrix((s, s), mask_missing=True,
                                                         _mappable=_mappable)

                if matrix.shape[0] != matrix.shape[1] or matrix.shape[0] != shape:
                    continue

                cmatrix += matrix / matrix_expected
                inverted_mask = ~matrix.mask
                counter_matrix += inverted_mask.astype('int')
                counter += 1
                pb.update(i)

    if counter is None:
        raise ValueError("No valid regions found!")

    cmatrix /= counter_matrix

    if norm:
        cmatrix = np.log2(cmatrix)

    return cmatrix


def _aggregate_region_bins(hic, region, offset=0):
    region_slice = hic.region_bins(region)
    return region_slice.start-offset, region_slice.stop-offset


def extract_submatrices(hic, region_pairs, norm=False,
                        log=True, cache=True, mask=True, mask_inf=True,
                        keep_invalid=False, orient_strand=False):
    cl = hic.chromosome_lengths
    cb = hic.chromosome_bins

    valid_region_pairs = defaultdict(list)
    invalid_region_pairs = list()
    logger.info("Checking region pair validity...")
    valid, invalid = 0, 0
    for ix, (region1, region2) in enumerate(region_pairs):
        is_invalid = False
        if region1 is None or region2 is None:
            is_invalid = True
        elif region1.start < 1 or region1.chromosome not in cl or region1.end > cl[region1.chromosome]:
            is_invalid = True
        elif region2.start < 1 or region2.chromosome not in cl or region2.end > cl[region2.chromosome]:
            is_invalid = True
        if is_invalid:
            invalid += 1
            invalid_region_pairs.append(ix)
        else:
            valid += 1
            valid_region_pairs[(region1.chromosome, region2.chromosome)].append((ix, region1, region2))
    logger.info("{}/{} region pairs are invalid".format(invalid, valid))

    logger.info("Calculating mappability...")
    mappable = hic.mappable()

    intra_expected, inter_expected = dict(), None
    if norm:
        logger.info("Calculating expected values...")
        _, intra_expected, inter_expected = hic.expected_values()

    order = []
    matrices = []
    with RareUpdateProgressBar(max_value=valid, prefix='Matrices') as pb:
        current_matrix = 0
        for (chromosome1, chromosome2), regions_pairs_by_chromosome in valid_region_pairs.items():
            if cache:
                matrix = hic.as_matrix((chromosome1, chromosome2), mask_missing=mask,
                                       _mappable=mappable)
                offset1 = cb[chromosome1][0]
                offset2 = cb[chromosome2][0]
            else:
                matrix = hic
                offset1 = 0
                offset2 = 0

            for (region_ix, region1, region2) in regions_pairs_by_chromosome:
                current_matrix += 1
                region1_bins = _aggregate_region_bins(hic, region1, offset1)
                region2_bins = _aggregate_region_bins(hic, region2, offset2)

                if cache:
                    ms = matrix[region1_bins[0]:region1_bins[1], region2_bins[0]: region2_bins[1]]
                    m = ms.copy()
                    del ms
                else:
                    s1 = slice(region1_bins[0], region1_bins[1])
                    s2 = slice(region2_bins[0], region2_bins[1])
                    m = hic.as_matrix((s1, s2), mask_missing=mask, _mappable=mappable)

                if norm:
                    e = np.ones(m.shape)
                    if chromosome1 != chromosome2:
                        e.fill(inter_expected)
                    else:
                        for i, row in enumerate(range(region1_bins[0], region1_bins[1])):
                            for j, col in enumerate(range(region2_bins[0], region2_bins[1])):
                                ix = abs(col - row)
                                e[i, j] = intra_expected[chromosome1][ix]

                    if log:
                        m = np.log2(m/e)
                        m[np.isnan(m)] = 0.
                    else:
                        m = m/e
                        m[np.isnan(m)] = 1

                if mask_inf:
                    m_mask = np.isinf(m)
                    if not hasattr(m, 'mask'):
                        m = np.ma.masked_where(m_mask, m)
                    m.mask += m_mask

                if orient_strand and region1.is_reverse() and region2.is_reverse():
                    m = np.flip(np.flip(m, 0), 1)

                pb.update(current_matrix)
                matrices.append(m)
                order.append(region_ix)

            if cache:
                del matrix

    if keep_invalid:
        for region_ix in invalid_region_pairs:
            matrices.append(None)
            order.append(region_ix)

    return [matrices[ix] for ix in np.argsort(order)]


def _rescale_oe_matrix(matrix, bin_size, scaling_exponent=-0.25):
    rm = np.zeros(matrix.shape)
    b = bin_size
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            v = (abs(i - j) * b + b) ** scaling_exponent
            try:
                if matrix.mask[i, j]:
                    continue
            except AttributeError:
                pass
            rm[i, j] = v * matrix[i, j]
            rm[j, i] = v * matrix[j, i]
    return rm


def aggregate_boundaries(hic, boundary_regions, window=200000,
                         rescale=False, scaling_exponent=-0.25,
                         **kwargs):
    kwargs.setdefault('norm', True)
    kwargs.setdefault('keep_invalid', False)
    kwargs.setdefault('log', True)

    region_pairs = []
    for region in boundary_regions:
        new_start = int(region.center - int(window / 2))
        new_end = int(region.center + int(window / 2))
        new_region = GenomicRegion(chromosome=region.chromosome, start=new_start, end=new_end,
                                   strand=region.strand)
        region_pairs.append((new_region, new_region))

    counter_matrix = None
    matrix_sum = None
    for m in extract_submatrices(hic, region_pairs, **kwargs):
        if counter_matrix is None:
            shape = m.shape
            counter_matrix = np.zeros(shape)
            matrix_sum = np.zeros(shape)

        if hasattr(m, 'mask'):
            inverted_mask = ~m.mask
            counter_matrix += inverted_mask.astype('int')
        else:
            counter_matrix += np.ones(counter_matrix.shape)

        matrix_sum += m

    am = matrix_sum / counter_matrix

    if rescale:
        am = _rescale_oe_matrix(am, hic.bin_size, scaling_exponent=scaling_exponent)

    return am


def _tad_matrix_iterator(hic, tad_regions, absolute_extension=0, relative_extension=1., **kwargs):
    region_pairs = []
    for region in tad_regions:
        new_region = region.expand(absolute=absolute_extension, relative=relative_extension)
        region_pairs.append((new_region, new_region))

    for m in extract_submatrices(hic, region_pairs, **kwargs):
        yield m


def aggregate_tads(hic, tad_regions, pixels=90, rescale=False, scaling_exponent=-0.25,
                   interpolation='nearest', keep_mask=True,
                   absolute_extension=0, relative_extension=1.0,
                   **kwargs):
    kwargs.setdefault('norm', True)
    kwargs.setdefault('keep_invalid', False)
    kwargs.setdefault('log', True)

    shape = (pixels, pixels)
    counter_matrix = np.zeros(shape)
    matrix_sum = np.zeros(shape)
    for m in _tad_matrix_iterator(hic, tad_regions,
                                  absolute_extension=absolute_extension,
                                  relative_extension=relative_extension, **kwargs):
        ms = imresize(m, shape, interp=interpolation, mode='F')

        if keep_mask and hasattr(ms, 'mask'):
            mask = imresize(m.mask, shape, interp='nearest').astype('bool')
            ms = np.ma.masked_where(mask, ms)
            inverted_mask = ~mask
            counter_matrix += inverted_mask.astype('int')
        else:
            counter_matrix += np.ones(shape)

        matrix_sum += ms

    am = matrix_sum/counter_matrix

    if rescale:
        am = _rescale_oe_matrix(am, hic.bin_size, scaling_exponent=scaling_exponent)

    return am


def tad_strength(hic, tad_regions, **kwargs):
    kwargs.setdefault('norm', True)
    kwargs.setdefault('keep_invalid', False)
    kwargs.setdefault('log', False)
    kwargs['relative_extension'] = 1.
    kwargs['absolute_extension'] = 0

    tad_strengths = []
    for m in _tad_matrix_iterator(hic, tad_regions, **kwargs):
        tl = int(m.shape[0]/3)
        upper_third = slice(0, tl)
        middle_third = slice(tl, 2*tl)
        lower_third = slice(2*tl, m.shape[0])
        tad_sum = m[middle_third, middle_third].sum() / np.logical_not(m.mask)[middle_third, middle_third].sum()
        upper_sum = m[upper_third, middle_third].sum() / np.logical_not(m.mask)[upper_third, middle_third].sum()
        lower_sum = m[lower_third, upper_third].sum() / np.logical_not(m.mask)[lower_third, upper_third].sum()
        ts = float(tad_sum / ((upper_sum + lower_sum) / 2))
        tad_strengths.append(np.log2(ts))
    return tad_strengths


def _loop_regions_from_bedpe(bedpe):
    anchors = []
    for region in bedpe.regions:
        a1 = GenomicRegion(chromosome=region.chromosome1, start=region.start1, end=region.end1)
        a2 = GenomicRegion(chromosome=region.chromosome2, start=region.start2, end=region.end2)
        anchors.append((a1, a2))
    return anchors


def _loop_matrix_iterator(hic, loop_regions, pixels=16, **kwargs):
    left = int(pixels / 2)
    right = left if pixels % 2 == 1 else left - 1

    if isinstance(loop_regions, Bedpe):
        loop_regions = _loop_regions_from_bedpe(loop_regions)

    bin_size = hic.bin_size
    region_pairs = []
    invalid = 0
    for (anchor1, anchor2) in loop_regions:
        a1 = GenomicRegion(chromosome=anchor1.chromosome, start=anchor1.center, end=anchor1.center)
        a2 = GenomicRegion(chromosome=anchor2.chromosome, start=anchor2.center, end=anchor2.center)

        try:
            r1 = list(hic.regions(a1))[0].copy()
            r2 = list(hic.regions(a2))[0].copy()
            r1.start -= left * bin_size
            r1.end += right * bin_size
            r2.start -= left * bin_size
            r2.end += right * bin_size
            region_pairs.append((r1, r2))
        except IndexError:
            invalid += 1
            region_pairs.append((None, None))

    if invalid > 0:
        logger.warning("{} region pairs invalid, most likely due to missing chromosome data".format(invalid))

    for m in extract_submatrices(hic, region_pairs, **kwargs):
        yield m


def aggregate_loops(hic, loop_regions, pixels=16, **kwargs):
    kwargs.setdefault('norm', True)
    kwargs.setdefault('keep_invalid', False)
    kwargs.setdefault('log', True)

    shape = (pixels, pixels)
    counter_matrix = np.zeros(shape)
    matrix_sum = np.zeros(shape)
    for m in _loop_matrix_iterator(hic, loop_regions, pixels=pixels, **kwargs):
        if hasattr(m, 'mask'):
            inverted_mask = ~m.mask
            counter_matrix += inverted_mask.astype('int')
        else:
            counter_matrix += np.ones(shape)
        matrix_sum += m

    return matrix_sum/counter_matrix


def loop_strength(hic, loop_regions, pixels=16, **kwargs):
    kwargs.setdefault('log', False)
    kwargs.setdefault('norm', True)
    kwargs['keep_invalid'] = True

    if isinstance(loop_regions, Bedpe):
        loop_regions = _loop_regions_from_bedpe(loop_regions)

    # generating new regions
    new_region_pairs = []  # 0: original, 1: control left, 2: control right
    for (region1, region2) in loop_regions:
        d = int(abs(region1.center - region2.center))
        new_left = GenomicRegion(chromosome=region1.chromosome, start=region1.start - d, end=region1.end - d)
        new_right = GenomicRegion(chromosome=region1.chromosome, start=region2.start + d, end=region2.end + d)
        new_region_pairs.append((region1, region2))
        new_region_pairs.append((new_left, region1))
        new_region_pairs.append((region2, new_right))

    original, left, right = [], [], []
    for i, m in enumerate(_loop_matrix_iterator(hic, new_region_pairs, pixels=pixels, **kwargs)):
        if m is not None:
            value = float(m.sum()/np.logical_not(m.mask).sum())
        else:
            value = None

        if i % 3 == 0:
            original.append(value)
        elif i % 3 == 1:
            left.append(value)
        else:
            right.append(value)

    ratios = []
    for i in range(len(original)):
        if original[i] is None:
            continue
        if left[i] is None and right[i] is None:
            continue

        if left[i] is None:
            r = original[i]/right[i]
        elif right[i] is None:
            r = original[i]/left[i]
        else:
            r = original[i]/((left[i]+right[i])/2)
        ratios.append(np.log2(r))
    return ratios


def vector_enrichment_profile(oe, vector, mappable=None, per_chromosome=True,
                              percentiles=(20.0, 40.0, 60.0, 80.0, 100.0),
                              symmetric_at=None, exclude_chromosomes=()):
    if len(exclude_chromosomes) > 0:
        chromosome_bins = oe.chromosome_bins
        exclude_vector = []
        for chromosome in oe.chromosomes():
            if chromosome not in exclude_chromosomes:
                b = chromosome_bins[chromosome]
                for v in vector[b[0]:b[1]]:
                    exclude_vector.append(v)
                # exclude_vector += vector[b[0]:b[1]]
    else:
        exclude_vector = vector
    exclude_vector = np.array(exclude_vector)

    if symmetric_at is not None:
        lv = exclude_vector[exclude_vector <= symmetric_at]
        gv = exclude_vector[exclude_vector > symmetric_at]
        lv_cutoffs = np.nanpercentile(lv, percentiles)
        gv_cutoffs = np.nanpercentile(gv, percentiles)
        bin_cutoffs = np.concatenate((lv_cutoffs, gv_cutoffs))
    else:
        bin_cutoffs = np.nanpercentile(exclude_vector, percentiles)

    bins = []
    for value in vector:
        bins.append(bisect_left(bin_cutoffs, value))

    s = len(bin_cutoffs)
    m = np.zeros((s, s))
    c = np.zeros((s, s))

    if mappable is None:
        mappable = oe.mappable()

    if per_chromosome:
        for chromosome in oe.chromosomes():
            if chromosome in exclude_chromosomes:
                continue
            oem = oe[chromosome, chromosome]
            for i, row_region in enumerate(oem.row_regions):
                if not mappable[row_region.ix]:
                    continue
                i_bin = s - bins[row_region.ix] - 1
                for j, col_region in enumerate(oem.col_regions):
                    if not mappable[col_region.ix]:
                        continue
                    j_bin = s - bins[col_region.ix] - 1
                    value = oem[i, j]

                    m[i_bin, j_bin] += value
                    c[i_bin, j_bin] += 1
                    m[j_bin, i_bin] += value
                    c[j_bin, i_bin] += 1
    else:
        oem = oe[:]
        for i in range(oem.shape):
            if not mappable[i]:
                continue
            i_bin = s - bins[i] - 1
            for j in range(i, oem.shape):
                if not mappable[j]:
                    continue
                j_bin = s - bins[j] - 1
                value = oem[i, j]

                m[i_bin, j_bin] += value
                c[i_bin, j_bin] += 1
                m[j_bin, i_bin] += value
                c[j_bin, i_bin] += 1

    m /= c
    # m[c == 0] = 0
    rev_cutoffs = bin_cutoffs[::-1]
    for i in range(len(rev_cutoffs) - 1, 0, -1):
        if np.isclose(rev_cutoffs[i - 1], rev_cutoffs[i]):
            m[:, i - 1] = m[:, i]
            m[i - 1, :] = m[i, :]

    return np.log2(m), rev_cutoffs


def region_score_enrichment_profile(hic, regions, per_chromosome=True,
                                    percentiles=(20.0, 40.0, 60.0, 80.0, 100.0)):

    regions_chromosomes = set(hic.chromosomes())
    vector = []
    for chromosome in hic.chromosomes():
        if chromosome not in regions_chromosomes:
            continue

        v = [r.score for r in regions.subset(chromosome)]

        lh = len(list(hic.regions(chromosome)))
        if lh != len(v):
            raise ValueError("The number of values in chromosome {} is not equal to the "
                             "number of regions in the Hi-C object ({} vs {})!".format(chromosome,
                                                                                       len(v), lh))

        vector += v

    with ObservedExpectedRatio(hic, per_chromosome=per_chromosome) as oe:
        mappable = hic.mappable()
        return vector_enrichment_profile(oe, vector, mappable=mappable, per_chromosome=per_chromosome,
                                         percentiles=percentiles)


def ab_enrichment_profile(hic, genome, percentiles=(20.0, 40.0, 60.0, 80.0, 100.0),
                          per_chromosome=True, only_gc=False, symmetric_at=None,
                          exclude_chromosomes=()):

    logger.info("Generating profile...")
    with ObservedExpectedRatio(hic, per_chromosome=per_chromosome) as oe:
        if only_gc:
            # calculate GC content
            if isinstance(genome, string_types):
                logger.info("Loading genome...")
                genome = Genome.from_string(genome)

            logger.info("Calculating GC content...")
            gc_content = [np.nan] * len(hic.regions)
            for chromosome in hic.chromosomes():
                logger.info("{}".format(chromosome))
                chromosome_sequence = genome[chromosome].sequence
                for region in hic.regions(chromosome):
                    s = chromosome_sequence[region.start - 1:region.end]
                    gc_content[region.ix] = calculate_gc_content(s)
            gc_content = np.array(gc_content)
            ev = gc_content
        else:
            with ABDomainMatrix(oe, ratio=False, per_chromosome=per_chromosome) as ab:
                with ABDomains(ab, genome=genome) as abd:
                    ev = np.array(abd.ab_domain_eigenvector())

        mappable = hic.mappable()
        return vector_enrichment_profile(oe, ev, mappable=mappable, per_chromosome=per_chromosome,
                                         percentiles=percentiles, symmetric_at=symmetric_at,
                                         exclude_chromosomes=exclude_chromosomes)


def contact_directionality_bias(hic, regions, distance=1000000, region_anchor='center', **kwargs):
    forward_region_pairs = []
    reverse_region_pairs = []
    for region in regions:
        pos = int(getattr(region, region_anchor))
        new_region = GenomicRegion(chromosome=region.chromosome, start=pos, end=pos, strand=region.strand)
        if region.is_forward():
            forward_region_pairs.append((new_region, new_region.expand(absolute=distance)))
        else:
            reverse_region_pairs.append((new_region, new_region.expand(absolute=distance)))

    cumulative_forward = np.zeros(int(distance / hic.bin_size) * 2 + 1)
    count_forward = np.zeros(int(distance / hic.bin_size) * 2 + 1)
    for matrix in extract_submatrices(hic, forward_region_pairs, **kwargs):
        cumulative_forward += matrix[0, :]
        if hasattr(matrix, 'mask'):
            inverted_mask = ~matrix.mask
            count_forward += inverted_mask.astype('int')[0, :]
        else:
            count_forward += np.ones(count_forward.shape)

    cumulative_reverse = np.zeros(int(distance / hic.bin_size) * 2 + 1)
    count_reverse = np.zeros(int(distance / hic.bin_size) * 2 + 1)
    for matrix in extract_submatrices(hic, reverse_region_pairs, **kwargs):
        cumulative_reverse += matrix[0, :]
        if hasattr(matrix, 'mask'):
            inverted_mask = ~matrix.mask
            count_reverse += inverted_mask.astype('int')[0, :]
        else:
            count_reverse += np.ones(count_reverse.shape)

    avg_forward = cumulative_forward / count_forward
    avg_reverse = cumulative_reverse / count_reverse

    bin_size = hic.bin_size
    d = []
    ratio_forward = []
    ratio_reverse = []
    center = int(len(avg_forward)/2)
    for ix in range(center + 1):
        d.append(ix * bin_size)
        ratio_forward.append(avg_forward[center + ix] / avg_forward[center - ix])
        ratio_reverse.append(avg_reverse[center + ix] / avg_reverse[center - ix])

    return d, ratio_forward, ratio_reverse


"""
Filters for architecture objects
"""


class ZeroWeightFilter(MatrixArchitecturalRegionFeatureFilter):
    """
    Filter edges where every associated weight is 0.
    """
    def __init__(self, mask=None):
        """
        Initialize filter with chosen parameters.

        :param mask: Optional Mask object describing the mask
                     that is applied to filtered edges.
        """
        MatrixArchitecturalRegionFeatureFilter.__init__(self, mask=mask)

    def valid_edge(self, edge):
        """
        Check if an edge is on (or near) the diagonal of the :class:`~Hic` matrix.
        """
        i = 0
        while True:
            try:
                weight = getattr(edge, 'weight_' + str(i))
                if weight == 0:
                    return False
                i += 1
            except AttributeError:
                return True


class ExpectedObservedCollectionFilter(MatrixArchitecturalRegionFeatureFilter):
    """
    Filter edges where every associated weight is 0.
    """
    def __init__(self, collection, fold_change=1, filter_when_single_invalid=False, mask=None):
        """
        Initialize filter with chosen parameters.

        :param mask: Optional Mask object describing the mask
                     that is applied to filtered edges.
        """
        MatrixArchitecturalRegionFeatureFilter.__init__(self, mask=mask)

        # count number of matrices
        n_hic = 0
        while 'weight_' + str(n_hic) in collection.field_names:
            n_hic += 1

        self.intra_expected = dict()
        self.inter_expected = dict()
        for i in range(n_hic):
            with ExpectedContacts(collection, weight_column='weight_' + str(i)) as ex:
                self.intra_expected[i] = ex.intra_expected()
                self.inter_expected[i] = ex.inter_expected()

        self.n_hic = n_hic
        self.fold_change = fold_change
        self.filter_single = filter_when_single_invalid
        self.regions_dict = collection.regions_dict

    def valid_edge(self, edge):
        """
        Check if an edge weight is at least fold_change above
        the expected weight for this contact.
        """
        source = edge.source
        sink = edge.sink
        intra = False
        if self.regions_dict[source].chromosome == self.regions_dict[sink].chromosome:
            intra = True
        n_failed = 0
        for i in range(self.n_hic):
            if intra:
                expected = self.intra_expected[i][abs(sink-source)]
            else:
                expected = self.inter_expected[i]

            if getattr(edge, 'weight_' + str(i)) < self.fold_change*expected:
                if self.filter_single:
                    return False
                else:
                    n_failed += 1
        if n_failed == self.n_hic:
            return False
        return True


class BackgroundLigationCollectionFilter(MatrixArchitecturalRegionFeatureFilter):
    """
    Filter edges where every associated weight is 0.
    """
    def __init__(self, collection, fold_change=1, filter_when_single_invalid=False,
                 all_contacts=True, mask=None):
        """
        Initialize filter with chosen parameters.

        :param mask: Optional Mask object describing the mask
                     that is applied to filtered edges.
        """
        MatrixArchitecturalRegionFeatureFilter.__init__(self, mask=mask)

        regions_dict = collection.regions_dict
        # count number of matrices
        n_hic = 0
        while 'weight_' + str(n_hic) in collection.field_names:
            n_hic += 1

        inter_count = defaultdict(int)
        inter_sum = defaultdict(int)
        for edge in collection.edges(lazy=True):
            intra = regions_dict[edge.source].chromosome == regions_dict[edge.sink].chromosome
            for i in range(n_hic):
                if intra:
                    inter_count[i] += 1
                    inter_sum[i] += getattr(edge, 'weight_' + str(i))

        if all_contacts:
            with PossibleContacts(collection, weight_column='weight_0') as pc:
                for i in range(n_hic):
                    inter_count[i] = pc.inter_possible()

        self.cutoff = dict()
        for i in range(n_hic):
            if inter_count[i] == 0:
                self.cutoff[i] = 0
            else:
                self.cutoff[i] = fold_change*(inter_sum[i]/inter_count[i])

        self.n_hic = n_hic
        self.filter_single = filter_when_single_invalid

    def valid_edge(self, edge):
        """
        Check if an edge weight is at least fold_change above
        the expected weight for this contact.
        """
        n_failed = 0
        for i in range(self.n_hic):
            if getattr(edge, 'weight_' + str(i)) < self.cutoff[i]:
                if self.filter_single:
                    return False
                else:
                    n_failed += 1
        if n_failed == self.n_hic:
            return False
        return True


class MinMaxDistanceCollectionFilter(MatrixArchitecturalRegionFeatureFilter):
    """
    Filter edges where every associated weight is 0.
    """
    def __init__(self, collection, min_distance, max_distance, mask=None):
        """
        Initialize filter with chosen parameters.

        :param mask: Optional Mask object describing the mask
                     that is applied to filtered edges.
        """
        MatrixArchitecturalRegionFeatureFilter.__init__(self, mask=mask)

        self.min_distance_bins = collection.distance_to_bins(min_distance)
        self.max_distance_bins = collection.distance_to_bins(max_distance)
        self.regions_dict = collection.regions_dict

    def valid_edge(self, edge):
        """
        Check if an edge weight is at least fold_change above
        the expected weight for this contact.
        """
        source = edge.source
        sink = edge.sink
        distance_bins = sink-source

        # inter-chromosomal are valid by default
        if self.regions_dict[source].chromosome != self.regions_dict[sink].chromosome:
            return True

        if self.min_distance_bins <= distance_bins <= self.max_distance_bins:
            return True

        return False


class DiagonalCollectionFilter(MatrixArchitecturalRegionFeatureFilter):
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
        MatrixArchitecturalRegionFeatureFilter.__init__(self, mask=mask)
        self.distance = distance

    def valid_edge(self, edge):
        """
        Check if an edge is on (or near) the diagonal of the :class:`~Hic` matrix.
        """
        if abs(edge.source-edge.sink) <= self.distance:
            return False
        return True


class BackgroundLigationFilter(HicEdgeFilter):
    """
    Filter a :class:`~HicEdge` if it does not have a weight
    larger than  [fold_change*background ligation frequency].

    Background ligation frequency is estimated as the average
    of all non-zero inter-chromosomal contacts of this Hic object.
    """
    def __init__(self, hic, fold_change=5, all_contacts=False, mask=None):
        """
        Initialize filter with these settings.

        :param hic: The :class:`~Hic` object that this
                    filter will be called on. Needed for
                    contact count calculation.
        :param fold_change: Lowest acceptable edge weight is calculated
                            as fold_change*(inter_sum/inter_count)
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered edges.
        """
        HicEdgeFilter.__init__(self, mask=mask)

        regions_dict = hic.regions_dict

        inter_count = 0
        inter_sum = 0
        for edge in hic.edges(lazy=True):
            if regions_dict[edge.source].chromosome != regions_dict[edge.sink].chromosome:
                inter_count += 1
                inter_sum += edge.weight

        if all_contacts:
            with PossibleContacts(hic) as pc:
                inter_count = pc.inter_possible()

        if inter_count == 0:
            self.cutoff = 0
        else:
            self.cutoff = fold_change*(inter_sum/inter_count)

    def valid_edge(self, edge):
        """
        Check if an edge weight is below background ligation frequency.
        """
        if edge.weight < self.cutoff:
            return False
        return True


class ExpectedObservedEnrichmentFilter(HicEdgeFilter):
    """
    Filter a :class:`~HicEdge` if it does not have a weight
    larger than fold_change times its expected value.
    """
    def __init__(self, hic, fold_change=2, mask=None):
        """
        Initialize filter with these settings.

        :param hic: The :class:`~Hic` object that this
                    filter will be called on. Needed for
                    expected contact count calculation.
        :param fold_change: Lowest acceptable edge weight is calculated
                            as fold_change*expected
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered edges.
        """
        HicEdgeFilter.__init__(self, mask=mask)

        with ExpectedContacts(hic) as ex:
            self.intra_expected = ex.intra_expected()
            self.inter_expected = ex.inter_expected()
        self.regions_dict = hic.regions_dict
        self.fold_change = fold_change

    def valid_edge(self, edge):
        """
        Check if an edge weight is at least fold_change above
        the expected weight for this contact.
        """
        source = edge.source
        sink = edge.sink
        if self.regions_dict[source].chromosome == self.regions_dict[sink].chromosome:
            expected = self.intra_expected[abs(sink-source)]
        else:
            expected = self.inter_expected

        if edge.weight < self.fold_change*expected:
            return False
        return True
