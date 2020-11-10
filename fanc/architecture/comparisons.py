from __future__ import division

import tables
import logging
import operator
from collections import defaultdict
from genomic_regions import GenomicRegion, as_region
from future.utils import string_types

from ..tools.general import RareUpdateProgressBar
from ..regions import RegionsTable
from ..matrix import RegionMatrixTable, RegionMatrix
from .domains import RegionScoreParameterTable

import numpy as np
from sklearn.decomposition import PCA

logger = logging.getLogger(__name__)


def _edge_collection(*hics, region=None, scale=True,
                     filters=None, **kwargs):
    """
    Get all weights from the same edge across different Hi-C matrices.

    :param hics: Hic/matrix objects or list of Hic/Matrix objects
    :param region: Optional region subset (e.g. a single chromosome)
    :param scale: True if Hic matrices should be scaled to the
                  number of valid pairs
    :return: dict of lists
    """
    inter_chromosomal = kwargs.pop("inter_chromosomal", False)

    if len(hics) == 1:
        hics = hics[0]

    if filters is None:
        filters = []

    scaling_factors = [1.0] * len(hics)
    if scale:
        reference_hic = None
        for i, hic in enumerate(hics):
            if i == 0:
                reference_hic = hic
            else:
                s = hic.scaling_factor(reference_hic)
                scaling_factors[i] = s

    total_edges = sum(len(hic.edges) for hic in hics)

    if isinstance(region, GenomicRegion) or isinstance(region, string_types) or region is None:
        regions = [region]
    else:
        regions = region

    region_pairs = []
    for i in range(len(regions)):
        region1 = as_region(regions[i]) if regions[i] is not None else None
        if region1 is None:
            region_pairs.append(None)
            continue

        for j in range(i, len(regions)):
            region2 = as_region(regions[j]) if regions[j] is not None else None
            if not inter_chromosomal and region1.chromosome != region2.chromosome:
                continue
            region_pairs.append((region1, region2))

    edges = defaultdict(list)
    j = 0
    with RareUpdateProgressBar(max_value=total_edges, prefix='Edge collection') as pb:
        for i, hic in enumerate(hics):
            logger.debug("Adding Hic {} ({}) to edge collection".format(i, region))

            for region_pair in region_pairs:
                for edge in hic.edges(region_pair, lazy=True,
                                      inter_chromosomal=inter_chromosomal, **kwargs):
                    weight = edge.weight * scaling_factors[i]
                    source, sink = edge.source, edge.sink

                    weight_list = edges[(source, sink)]
                    while len(weight_list) < i:
                        weight_list.append(0)

                    if len(weight_list) <= i:
                        weight_list.append(weight)

                    pb.update(j)
                    j += 1

    total = len(hics)
    for _, weights in edges.items():
        while len(weights) < total:
            weights.append(0)

    before_filtering = len(edges)
    for key in list(edges.keys()):
        for f in filters:
            if not f.valid(key[0], key[1], edges[key]):
                del edges[key]
    after_filtering = len(edges)

    logger.info("Valid edges: {}/{}".format(after_filtering, before_filtering))

    return edges


class EdgeCollectionFilter(object):
    def __init__(self):
        pass

    def valid(self, source, sink, weights):
        raise NotImplementedError("Subclasses of EdgeCollectionFilter must implement "
                                  "'valid'!")


class ObservedExpectedFilter(EdgeCollectionFilter):
    def __init__(self, hics, fold_change_cutoff=1.0, min_valid_edges=1,
                 oe_per_chromosome=True):
        EdgeCollectionFilter.__init__(self)

        self._cutoff = fold_change_cutoff
        self._min_valid = min_valid_edges
        self._oe_per_chromosome = oe_per_chromosome

        self._intra_expected = {}
        self._inter_expected = {}
        self._intra_expected_chromosome = {}
        for i, hic in enumerate(hics):
            intra_expected, intra_expected_chromosome, inter_expected = hic.expected_values()
            self._intra_expected[i] = intra_expected
            self._inter_expected[i] = inter_expected
            self._intra_expected_chromosome[i] = intra_expected_chromosome

        self._chromosome_dict = {}
        for region in hics[0].regions:
            self._chromosome_dict[region.ix] = region.chromosome

    def valid(self, source, sink, weights):
        source_chromosome = self._chromosome_dict[source]
        sink_chromosome = self._chromosome_dict[sink]

        count_valid = 0
        for i, weight in enumerate(weights):
            if source_chromosome == sink_chromosome:
                d = abs(sink - source)
                if self._oe_per_chromosome:
                    e = self._intra_expected_chromosome[i][source_chromosome][d]
                else:
                    e = self._intra_expected[i][d]
            else:
                e = self._inter_expected[i]

            fc = weight / e
            if np.isfinite(fc) and fc > self._cutoff:
                count_valid += 1

        if count_valid > self._min_valid:
            return True
        return False


class NonzeroFilter(EdgeCollectionFilter):
    def __init__(self):
        EdgeCollectionFilter.__init__(self)

    def valid(self, source, sink, weights):
        if 0 in weights:
            return False
        return True


class AbsoluteWeightFilter(EdgeCollectionFilter):
    def __init__(self, cutoff, min_above_cutoff=1):
        EdgeCollectionFilter.__init__(self)

        self.cutoff = cutoff
        self._min_valid = min_above_cutoff

    def valid(self, source, sink, weights):
        counter = 0
        for weight in weights:
            if weight >= self.cutoff:
                counter += 1

        if counter < self._min_valid:
            return False
        return True


class MinMaxDistanceFilter(EdgeCollectionFilter):
    def __init__(self, hics, min_distance=None, max_distance=None):
        self.bin_size = None
        for hic in hics:
            if self.bin_size is None:
                self.bin_size = hic.bin_size
            else:
                if self.bin_size != hic.bin_size:
                    raise ValueError("Hic objects must have same bin size! "
                                     "({}/{})".format(self.bin_size, hic.bin_size))
        EdgeCollectionFilter.__init__(self)
        self.min_distance = np.round(min_distance/self.bin_size) if min_distance is not None else None
        self.max_distance = np.round(max_distance/self.bin_size) if max_distance is not None else None

    def valid(self, source, sink, weights):
        d = abs(sink - source)
        if self.min_distance is not None and d < self.min_distance:
            return False
        if self.max_distance is not None and d > self.max_distance:
            return False
        return True


class ComparisonMatrix(RegionMatrixTable):

    _classid = 'COMPARISONMATRIX'

    def __init__(self, *args, **kwargs):
        RegionMatrixTable.__init__(self, *args, **kwargs)

    def compare(self, weight1, weight2):
        """
        Compare two edge weights.

        :param weight1: float
        :param weight2: float
        :return: float
        """
        raise NotImplementedError("Subclasses of ComparisonMatrix must implement "
                                  "'compare'")

    @classmethod
    def from_matrices(cls, matrix1, matrix2, file_name=None, tmpdir=None, mode='w',
                      log_cmp=False, ignore_infinite=True, ignore_zeros=False,
                      scale=True, **kwargs):
        """
        Create a comparison matrix from two compatible matrix objects.

        The resulting object can be treated like any other matrix in
        FAN-C, offering the same convenience functions for regions and
        edges.

        :param matrix1: First matrix object, such as a Hi-C matrix
        :param matrix2: Second matrix object, such as a Hi-C matrix
        :param file_name: Path to the comparison output file
        :param tmpdir: Optional. If ``True``, will work in temporary
                       directory until file is closed
        :param mode: Write mode of the output file. Only change this if you
                     know what you are doing - setting this to 'a' could lead
                     to unexpected consequences!
        :param log_cmp: If ``True``, log2-transform the comparison matrix
                        value after the comparison has been performed. Useful,
                        for example, for fold-change matrices
        :param ignore_infinite: If ``True``, will remove infinite values from
                                the final comparison matrix
        :param ignore_zeros: If ``True``, will only compare edge weights when
                             both are non-zero.
        :param scale: Scale matrices to the same sequencing depth (sum of all
                      edge weights) before the comparison. You can set this
                      to ``False`` if you know the type of normalisation you
                      performed already takes care of this.
        :param kwargs: Keyword arguments passed to
                       :func:`~fanc.matrix.RegionPairsContainer.edges`
        :return: :class:`~fanc.architecture.comparisons.ComparisonMatrix`
        """
        kwargs['lazy'] = True

        comparison_matrix = cls(file_name=file_name, mode=mode, tmpdir=tmpdir)
        comparison_matrix.add_regions(matrix1.regions, preserve_attributes=False)

        sf = 1.0
        if scale:
            if kwargs.get('oe', False):
                logger.warning('Not computing scaling factor due to O/E conversion.')
            else:
                sf = matrix2.scaling_factor(matrix1)

        compare = comparison_matrix.compare
        chromosomes = matrix1.chromosomes()
        n_chromosome_pairs = int(np.round(len(chromosomes)**2/2 + len(chromosomes)/2))
        current_chromosome_pair = 0
        with RareUpdateProgressBar(max_value=n_chromosome_pairs, prefix='Compare') as pb:
            for chr_i in range(len(chromosomes)):
                chromosome1 = chromosomes[chr_i]
                for chr_j in range(chr_i, len(chromosomes)):
                    chromosome2 = chromosomes[chr_j]

                    edges1 = dict()
                    for edge in matrix1.edges((chromosome1, chromosome2), **kwargs):
                        edges1[(edge.source, edge.sink)] = edge.weight

                    edges2 = dict()
                    for edge in matrix2.edges((chromosome1, chromosome2), **kwargs):
                        edges2[(edge.source, edge.sink)] = edge.weight * sf

                    for key, w1 in edges1.items():
                        try:
                            w2 = edges2[key]
                        except KeyError:
                            w2 = matrix2._default_value

                        if ignore_zeros and (w1 == 0 or w2 == 0):
                            continue
                        weight = compare(w1, w2)
                        if log_cmp:
                            weight = np.log2(weight)
                        if ignore_infinite and not np.isfinite(weight):
                            continue

                        comparison_matrix.add_edge_simple(key[0], key[1], weight=weight)

                    for key, w2 in edges2.items():
                        if key not in edges1:
                            w1 = matrix1._default_value

                            if ignore_zeros and (w1 == 0 or w2 == 0):
                                continue
                            weight = compare(w1, w2)
                            if log_cmp:
                                weight = np.log2(weight)
                            if ignore_infinite and not np.isfinite(weight):
                                continue

                            comparison_matrix.add_edge_simple(key[0], key[1], weight=weight)

                    current_chromosome_pair += 1
                    pb.update(current_chromosome_pair)

        comparison_matrix.flush()
        return comparison_matrix


class FoldChangeMatrix(ComparisonMatrix):

    _classid = 'FOLDCHANGEMATRIX'

    def __init__(self, *args, **kwargs):
        ComparisonMatrix.__init__(self, *args, **kwargs)

    def compare(self, weight1, weight2):
        try:
            return weight1 / weight2
        except ZeroDivisionError:
            return np.nan


class DifferenceMatrix(ComparisonMatrix):

    _classid = 'DIFFERENCEMATRIX'

    def __init__(self, *args, **kwargs):
        ComparisonMatrix.__init__(self, *args, **kwargs)

    def compare(self, weight1, weight2):
        return weight1 - weight2


class ComparisonScores(RegionScoreParameterTable):

    _classid = 'COMPARISONSCORES'

    def __init__(self, *args, **kwargs):
        RegionScoreParameterTable.__init__(self, *args, **kwargs)

    def compare(self, score1, score2):
        """
        Compare two edge weights.

        :param score1: float
        :param score2: float
        :return: float
        """
        raise NotImplementedError("Subclasses of ComparisonScores must implement "
                                  "'compare'")

    @classmethod
    def from_scores(cls, scores1, scores2, attributes=None, file_name=None, tmpdir=None,
                    log=False, field_prefix='cmp_', **kwargs):
        """
        Compare parameter-based scores in a :class:`~fanc.architecture.domains.RegionScoreParameterTable`.

        :param scores1: First :class:`~fanc.architecture.domains.RegionScoreParameterTable`
        :param scores2: Second :class:`~fanc.architecture.domains.RegionScoreParameterTable`
        :param attributes: If ``None``, will do all possible comparisons. Provide a list of
                           region attributes (e.g. ["insulation_1000000", "insulation_2000000"])
                           for specific comparisons.
        :param file_name: Optional path to an output file
        :param tmpdir: Optional. If ``True``, will work in temporary
                       directory until file is closed
        :param log: log2-transform values after comparison
        :param field_prefix: Prefix of the output field
        :param kwargs: Keyword arguments passed on to :func:`~genomic_regions.RegionBased.regions`
        :return: :class:`~fanc.architecture.comparisons.ComparisonScores`
        """
        # all matching parameters
        if attributes is None:
            attributes = []
            attributes1 = scores1._score_fields
            for a, p in zip(scores2._score_fields, scores2._parameters):
                if a in attributes1:
                    attributes.append(p)

        comparison_scores = cls(parameter_values=attributes, parameter_prefix=field_prefix,
                                file_name=file_name, mode='w', tmpdir=tmpdir)
        comparison_scores.add_regions(scores1.regions, preserve_attributes=False)

        region_pairs = list()
        region_ixs = dict()
        for region in scores1.regions:
            region_ixs[(region.chromosome, region.start, region.end)] = len(region_pairs)
            region_pairs.append([region, None])

        for region in scores2.regions(**kwargs):
            ix = region_ixs[(region.chromosome, region.start, region.end)]
            region_pairs[ix][1] = region

        for attribute in attributes:
            field = scores1._score_field_converter([attribute])[0]
            cmp_scores = []
            for r1, r2 in region_pairs:
                v1 = getattr(r1, field)
                v2 = getattr(r2, field)
                v_cmp = comparison_scores.compare(v1, v2)
                if log:
                    v_cmp = np.log2(v_cmp)
                cmp_scores.append(v_cmp)
            comparison_scores.scores(attribute, cmp_scores)

        return comparison_scores


class FoldChangeScores(ComparisonScores):

    _classid = 'FOLDCHANGESCORES'

    def __init__(self, *args, **kwargs):
        ComparisonScores.__init__(self, *args, **kwargs)

    def compare(self, score1, score2):
        return score1 / score2


class DifferenceScores(ComparisonScores):

    _classid = 'DIFFERENCESCORES'

    def __init__(self, *args, **kwargs):
        ComparisonScores.__init__(self, *args, **kwargs)

    def compare(self, score1, score2):
        return score1 - score2


class ComparisonRegions(RegionsTable):

    _classid = 'COMPARISONREGIONS'

    def __init__(self, *args, **kwargs):
        RegionsTable.__init__(self, *args, **kwargs)

    def compare(self, score1, score2):
        """
        Compare two edge weights.

        :param score1: float
        :param score2: float
        :return: float
        """
        raise NotImplementedError("Subclasses of ComparisonRegions must implement "
                                  "'compare'")

    @classmethod
    def from_regions(cls, region_based1, region_based2, attribute='score',
                     file_name=None, tmpdir=None, log=False, score_field=None,
                     **kwargs):
        """
        Compare genomic tracks with region-associated scores.

        All scores are assumed to be floats.

        :param region_based1: First :class:`~genomic_regions.RegionBased` object
        :param region_based2: Second :class:`~genomic_regions.RegionBased` object
        :param attribute: Name of the attribute to be compared. Typically "score"
        :param file_name: Optional path to an output file
        :param tmpdir: Optional. If ``True``, will work in temporary
                       directory until file is closed
        :param log: If ``True``, will log2-transform values after comparison
        :param score_field: Name of the attribute comparison scores will be saved to.
                            Will use ``attribute`` if not provided.
        :param kwargs: Keyword arguments passed on to :func:`~genomic_regions.RegionBased.regions`
        :return: :class:`~fanc.architecture.comparisons.ComparisonRegions`
        """
        if score_field is None:
            score_field = attribute

        logger.debug("Using scores from '{}' field, writing to '{}' field".format(attribute, score_field))

        comparison_regions = cls(file_name=file_name, mode='w', tmpdir=tmpdir,
                                 additional_fields={attribute: tables.Float32Col()})
        comparison_regions.add_regions(region_based1.regions, preserve_attributes=False)

        regions = dict()
        for region in region_based1.regions(**kwargs):
            regions[(region.chromosome, region.start, region.end)] = region

        scores = []
        for region2 in region_based2.regions:
            region1 = regions[(region2.chromosome, region2.start, region2.end)]
            v1 = getattr(region1, attribute)
            v2 = getattr(region2, attribute)

            v = comparison_regions.compare(v1, v2)
            if log:
                v = np.log2(v)
            scores.append(v)
        comparison_regions.region_data(score_field, scores)
        return comparison_regions


class FoldChangeRegions(ComparisonRegions):

    _classid = 'FOLDCHANGEREGIONS'

    def __init__(self, *args, **kwargs):
        ComparisonRegions.__init__(self, *args, **kwargs)

    def compare(self, score1, score2):
        return score1 / score2


class DifferenceRegions(ComparisonRegions):

    _classid = 'DIFFERENCEREGIONS'

    def __init__(self, *args, **kwargs):
        ComparisonRegions.__init__(self, *args, **kwargs)

    def compare(self, score1, score2):
        return score1 - score2


class SplitMatrix(object):
    def __init__(self, matrix_top, matrix_bottom, scaling_factor=1.):
        self.matrix_top = matrix_top
        self.matrix_bottom = matrix_bottom
        self.scaling_factor = scaling_factor

    def matrix(self, *args, **kwargs):
        """
        Return a :class:`~HicMatrix` where values above the diagonal
        are from this object and values below the diagonal are from
        another :class:`~Hic` object.

        "Above the diagonal" refers to the diagonal of the complete
        Hic object, not the diagonal of the returned matrix.

        :param key: A matrix selector. Use tuple to select row and
                    columns
        :param scaling_factor: Factor to scale the hic values. If None,
                               will be computed using
                               :func:`~Hic.scaling_factor`.
        :param kwargs: Keyword arguments passed to both matrices
                       :func:`~fanc.matrix.RegionMatrixContainer.matrix`
                       functions.
        :return: :class:`~RegionMatrix`
        """
        m_top = self.matrix_top.matrix(*args, **kwargs)

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
        m_bottom = self.matrix_bottom.matrix(*args, **kwargs) * self.scaling_factor
        top_indices = np.triu_indices(m_top.shape[0], matching_index, m_top.shape[1])
        m_bottom[top_indices] = m_top[top_indices]

        return RegionMatrix(m_bottom, row_regions=m_top.row_regions, col_regions=m_top.col_regions)


class EdgeCollectionSelector(object):
    def __init__(self):
        pass

    def filter_edge_collection(self, edge_collection, sample_size, *args, **kwargs):
        raise NotImplementedError("Subclasses of EdgeCollectionSelector must implement "
                                  "filter_edge_collection!")


class LargestVarianceSelector(EdgeCollectionSelector):
    def __init__(self):
        EdgeCollectionSelector.__init__(self)

    def filter_edge_collection(self, edge_collection, sample_size, *args, **kwargs):
        if sample_size is None:
            sample_size = len(edge_collection)

        variances = dict()
        for key, weights in edge_collection.items():
            v = np.nanvar(weights)
            if np.isfinite(v):
                variances[key] = v

        variances_sorted = sorted(variances.items(), key=operator.itemgetter(1), reverse=True)

        for i in range(min(len(variances_sorted), sample_size)):
            source, sink = variances_sorted[i][0]
            weights = edge_collection[(source, sink)]
            yield source, sink, weights


class PassthroughSelector(EdgeCollectionSelector):
    def __init__(self):
        EdgeCollectionSelector.__init__(self)

    def filter_edge_collection(self, edge_collection, sample_size, *args, **kwargs):
        if sample_size is None:
            sample_size = len(edge_collection)

        for i, (key, weights) in enumerate(edge_collection.items()):
            if i > sample_size:
                break
            yield key[0], key[1], weights


class LargestFoldChangeSelector(EdgeCollectionSelector):
    def __init__(self):
        EdgeCollectionSelector.__init__(self)

    def filter_edge_collection(self, edge_collection, sample_size=None, *args, **kwargs):
        if sample_size is None:
            sample_size = len(edge_collection)

        max_fold_changes = dict()
        for key, weights in edge_collection.items():
            min_weight = np.nanmin(weights)
            max_weight = np.nanmax(weights)

            fc = max_weight / min_weight
            if np.isfinite(fc):
                max_fold_changes[key] = fc

        fc_sorted = sorted(max_fold_changes.items(), key=operator.itemgetter(1), reverse=True)

        for i in range(min(len(fc_sorted), sample_size)):
            source, sink = fc_sorted[i][0]
            weights = edge_collection[(source, sink)]
            yield source, sink, weights


def hic_pca(*hics, sample_size=None, region=None, strategy='variance', scale=True, log=False,
            ignore_zeros=False, oe_enrichment=None, min_distance=None, max_distance=None,
            background_ligation=False, min_libraries_above_background=1,
            **kwargs):
    """
    Run a PCA analysis on a set of Hi-C matrices.

    Note: this is not a compartment analysis. Use :class:`~ABCompartmentMatrix` for
    that purpose.

    :param hics: Two or more Hi-C objects
    :param sample_size: Optional. Set an upper limit on the number of edges sampled
                        for this analysis. If not specified, will use all applicable
                        edges. Used in conjunction with ``strategy`` to determine
                        which edges to prioritise
    :param region: Optionally specify a region string to limit the PCA to that region.
    :param strategy: Sort order of edges. Used in conjunction with ``sample_size``.
                     One of "variance" (default), "fold-change", or "passthrough".
                     variance: edges sorted by variance of contact strength
                     across samples; fold-change: edges sorted by size of fold-change
                     of contact strength across samples; passthrough: unsorted, edges
                     appear in order they are stored in the object
    :param scale: If ``True`` (default), the matrix values are scaled to their sequencing
                  depth before running PCA. If you are using the default normalisation,
                  matrix entries correspond to contact probabilities and the margins are
                  equal to 1 and there is not Need for scaling, so you can set this to
                  ``False`` in order to save computational time.
    :param log: Log-transform contact strength prior to PCA
    :param ignore_zeros: Only use contacts that are non-zero in all samples
    :param oe_enrichment: Used for "fold-change" ``strategy``, at least on edge must
                          have an O/E of this value or larger.
    :param min_distance: regions must be at least this far apart (in base pairs) to be
                         used for PCA
    :param max_distance: regions must be at least this close together (in base pairs) to be
                         used for PCA
    :param background_ligation: Use the average inter-chromosomal contact strength as
                                background ligation signal and only use pixels where at
                                least ``min_libraries_above_background`` samples have
                                an O/E signal above background
    :param min_libraries_above_background: Minimum number of libraries/samples that
                                           must have an O/E above background ligation
                                           signal for each pixel.
    :param kwargs: Keyword arguments passed to :func:`~fanc.matrix.RegionMatrixTable.edges`
    :return: sklearn PCA object, PCA result
    """

    strategies = {
        'variance': LargestVarianceSelector(),
        'fold-change': LargestFoldChangeSelector(),
        'passthrough': PassthroughSelector(),
        None: PassthroughSelector(),
    }

    selector = strategies[strategy]

    filters = []
    if ignore_zeros:
        filters.append(NonzeroFilter())
    if oe_enrichment is not None:
        filters.append(ObservedExpectedFilter(hics,
                                              fold_change_cutoff=oe_enrichment,
                                              min_valid_edges=min_libraries_above_background))
    if min_distance is not None or max_distance is not None:
        filters.append(MinMaxDistanceFilter(hics, min_distance=min_distance, max_distance=max_distance))
    if background_ligation:
        hic = hics[0]
        _, _, inter_expected = hic.expected_values()
        filters.append(AbsoluteWeightFilter(cutoff=inter_expected,
                                            min_above_cutoff=min_libraries_above_background))

    edge_collection = _edge_collection(*hics, region=region, scale=scale,
                                       filters=filters, **kwargs)

    pca_edges = selector.filter_edge_collection(edge_collection, sample_size)

    values = []
    for source, sink, weights in pca_edges:
        values.append(weights)

    if not log:
        y = np.array(values)
    else:
        y = np.log(np.array(values))

    pca = PCA()
    pca_res = pca.fit_transform(y.T)
    logger.info("Variance explained: %s" % str(pca.explained_variance_ratio_))

    return pca, pca_res
