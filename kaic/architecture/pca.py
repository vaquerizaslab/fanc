from __future__ import division

import logging
import operator
from collections import defaultdict

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
        for edge in hic.edges(region, lazy=True, **kwargs):
            weight = edge.weight * scaling_factors[i]
            edges[(edge.source, edge.sink)].append(weight)

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


class BackgroundLigationFilter(EdgeCollectionFilter):
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


class EdgeCollectionSelector(object):
    def __init__(self):
        pass

    def filter_edge_collection(self, edge_collection, sample_size, *args, **kwargs):
        raise NotImplementedError("Subclasses of EdgeCollectionSelector must implement "
                                  "filter_edge_collection!")


class LargestVarianceSelector(EdgeCollectionFilter):
    def __init__(self):
        EdgeCollectionFilter.__init__(self)

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


class PassthroughSelector(EdgeCollectionFilter):
    def __init__(self):
        EdgeCollectionFilter.__init__(self)

    def filter_edge_collection(self, edge_collection, sample_size, *args, **kwargs):
        if sample_size is None:
            sample_size = len(edge_collection)

        for i, (key, weights) in enumerate(edge_collection.items()):
            if i > sample_size:
                break
            yield key[0], key[1], weights


class LargestFoldChangeSelector(EdgeCollectionFilter):
    def __init__(self):
        EdgeCollectionFilter.__init__(self)

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
            ignore_zeros=False, background_enrichment=None,
            min_libraries_above_background=1,
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
    if background_enrichment is not None:
        filters.append(BackgroundLigationFilter(hics,
                                                fold_change_cutoff=background_enrichment,
                                                min_valid_edges=min_libraries_above_background))

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
