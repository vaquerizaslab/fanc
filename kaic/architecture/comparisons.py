from __future__ import division

import tables
import logging
import operator
from collections import defaultdict

from ..tools.general import RareUpdateProgressBar
from ..regions import RegionsTable
from ..matrix import RegionMatrixTable, Edge
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
    if len(hics) == 1:
        hics = hics[0]

    if filters is None:
        filters = []

    scaling_factors = defaultdict(lambda: 1.0)
    if scale:
        reference_hic = None
        for i, hic in enumerate(hics):
            if i == 0:
                reference_hic = hic
            else:
                s = hic.scaling_factor(reference_hic)
                scaling_factors[i] = s

    edges = defaultdict(list)
    for i, hic in enumerate(hics):
        logger.debug("Adding Hic {} ({}) to edge collection".format(i, region))

        with RareUpdateProgressBar(max_value=len(hic.edges)) as pb:
            for j, edge in enumerate(hic.edges(region, lazy=True, **kwargs)):
                weight = edge.weight * scaling_factors[i]
                source, sink = edge.source, edge.sink
                edges[(source, sink)].append(weight)

                pb.update(j)

        for k, v in edges.items():
            if len(v) < i + 1:
                v.append(0)

    for key in list(edges.keys()):
        for f in filters:
            if not f.valid(key[0], key[1], edges[key]):
                del edges[key]

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
        if self.min_distance is None and d < self.min_distance:
            return False
        if self.max_distance is None and d > self.max_distance:
            return False


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
    def from_matrices(cls, matrix1, matrix2, file_name=None, tmpdir=None,
                      log=False, ignore_infinite=True, *args, **kwargs):
        comparison_matrix = cls(file_name=file_name, mode='w', tmpdir=tmpdir)
        comparison_matrix.add_regions(matrix1.regions, preserve_attributes=False)

        chromosomes = matrix1.chromosomes()
        for chr_i in range(len(chromosomes)):
            chromosome1 = chromosomes[chr_i]
            for chr_j in range(chr_i, len(chromosomes)):
                chromosome2 = chromosomes[chr_j]

                edges = _edge_collection(matrix1, matrix2, region=(chromosome1, chromosome2),
                                         *args, **kwargs)
                for (source, sink), weights in edges.items():
                    weight = comparison_matrix.compare(*weights)
                    if log:
                        weight = np.log2(weight)
                    if ignore_infinite and not np.isfinite(weight):
                        continue
                    comparison_matrix.add_edge([source, sink, weight])
        comparison_matrix.flush()
        return comparison_matrix


class FoldChangeMatrix(ComparisonMatrix):

    _classid = 'FOLDCHANGEMATRIX'

    def __init__(self, *args, **kwargs):
        ComparisonMatrix.__init__(self, *args, **kwargs)

    def compare(self, weight1, weight2):
        return weight1 / weight2


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
                    log=False, field_prefix='cmp_', *args, **kwargs):
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

        for region in scores2.regions:
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
                     file_name=None, tmpdir=None, log=False, score_field='score',
                     *args, **kwargs):
        comparison_regions = cls(file_name=file_name, mode='w', tmpdir=tmpdir,
                                 additional_fields={attribute: tables.Float32Col()})
        comparison_regions.add_regions(region_based1.regions, preserve_attributes=False)

        regions = dict()
        for region in region_based1.regions:
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

    if log:
        y = np.array(values)
    else:
        y = np.log(np.array(values))

    pca = PCA()
    pca_res = pca.fit_transform(y.T)
    logger.info("Variance explained: %s" % str(pca.explained_variance_ratio_))

    return pca, pca_res
