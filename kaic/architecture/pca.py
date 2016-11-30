from __future__ import division
from kaic.architecture.hic_architecture import ExpectedObservedCollectionFilter, \
    BackgroundLigationCollectionFilter, HicEdgeCollection
from kaic.data.genomic import GenomicRegion
from sklearn.decomposition import PCA
from abc import ABCMeta, abstractmethod
import tables as t
import numpy as np
import tempfile
import os.path
import logging
logger = logging.getLogger(__name__)


class HicCollectionWeightMeanVariance(HicEdgeCollection):
    _classid = 'HICCOLLECTIONWEIGHTMEANVARIANCE'

    def __init__(self, hics=None, file_name=None, mode='a', tmpdir=None,
                 only_intra_chromosomal=False, scale_libraries=False):
        additional_fields = {'var': t.Float32Col(pos=0), 'mean': t.Float32Col(pos=1)}
        HicEdgeCollection.__init__(self, hics, additional_fields=additional_fields, file_name=file_name,
                                   mode=mode, tmpdir=tmpdir, only_intra_chromosomal=only_intra_chromosomal)
        self.scale_libraries = scale_libraries

    def _calculate(self, *args, **kwargs):
        HicEdgeCollection._calculate(self, *args, **kwargs)

        weight_sums = np.zeros(len(self.hics))
        for edge in self.edges(lazy=True):
            weights = np.zeros(len(self.hics))
            for i in range(len(self.hics)):
                weight = getattr(edge, 'weight_' + str(i), 0.0)
                if np.isnan(weight):
                    weight = 0.0
                    setattr(edge, 'weight_' + str(i), 0.0)
                weights[i] = weight
                weight_sums[i] += weight
            edge.var = np.var(weights)
            edge.mean = np.mean(weights)
        self.flush()

        if self.scale_libraries:
            weight_ratios = weight_sums/weight_sums[0]
            for edge in self.edges(lazy=True):
                for i in range(len(self.hics)):
                    weight = getattr(edge, 'weight_' + str(i))
                    setattr(edge, 'weight_' + str(i), weight/weight_ratios[i])
        self.flush()

        self._edges.cols.var.create_csindex()


class PairSelection(object):
    __metaclass__ = ABCMeta

    def __init__(self):
        self.collection = None
        pass

    def set_collection(self, collection):
        self.collection = collection

    @abstractmethod
    def pair_selection(self, **kwargs):
        pass


class LargestVariancePairSelection(PairSelection):
    def __init__(self, sample_size=100000, regions=None, lazy=False):
        PairSelection.__init__(self)
        self.sample_size = sample_size
        self.regions = regions
        self.lazy = lazy

    def pair_selection(self, sample_size=None, lazy=None, regions=None):
        if lazy is None:
            lazy = self.lazy
        if sample_size is None:
            sample_size = self.sample_size
        if regions is None:
            regions = self.regions

        if isinstance(regions, str):
            regions = GenomicRegion.from_string(regions)

        regions_dict = None
        if regions is not None:
            regions_dict = self.collection.regions_dict

        for j, edge in enumerate(self.collection.edges_sorted('var', reverse=True, lazy=lazy)):
            if regions is None:
                yield edge
            else:
                source_region = regions_dict[edge.source]
                sink_region = regions_dict[edge.sink]
                if source_region.overlaps(regions) and sink_region.overlaps(regions):
                    yield edge

            if j >= sample_size:
                # raise StopIteration
                break


class PassthroughPairSelection(PairSelection):
    def __init__(self, sample_size=None, regions=None, lazy=False):
        PairSelection.__init__(self)
        self.sample_size = sample_size
        self.regions = regions
        self.lazy = lazy

    def pair_selection(self, sample_size=None, lazy=None, regions=None):
        if lazy is None:
            lazy = self.lazy
        if sample_size is None:
            sample_size = self.sample_size
        if regions is None:
            regions = self.regions

        if regions is None:
            edge_iterator = self.collection.edges(lazy=lazy)
        else:
            edge_iterator = self.collection.edge_subset(key=regions, lazy=lazy)

        for j, edge in enumerate(edge_iterator):
            yield edge

            if sample_size is not None and j >= sample_size:
                # raise StopIteration
                break


class LargestFoldChangePairSelection(PairSelection):
    def __init__(self, require_enriched=2, require_nonenriched=2, fold_change=1.5, sample_size=100000, regions=None,
                 lazy=False):
        PairSelection.__init__(self)
        self.sample_size = sample_size
        self.lazy = lazy
        self.fold_change = fold_change
        self.require_enriched = require_enriched
        self.require_nonenriched = require_nonenriched
        self.regions = regions

    def pair_selection(self, require_enriched=2, require_nonenriched=None,
                       fold_change=None, sample_size=None, lazy=None, regions=None):
        # count weights
        n_hic = 0
        while 'weight_' + str(n_hic) in self.collection.field_names:
            n_hic += 1

        # update parameters
        if fold_change is None:
            fold_change = self.fold_change
        if require_enriched is None:
            require_enriched = self.require_enriched
        if require_nonenriched is None:
            require_nonenriched = min(n_hic-require_enriched, self.require_nonenriched)
        if lazy is None:
            lazy = self.lazy
        if sample_size is None:
            sample_size = self.sample_size
        if regions is None:
            regions = self.regions

        if isinstance(regions, str):
            regions = GenomicRegion.from_string(regions)

        regions_dict = None
        if regions is not None:
            regions_dict = self.collection.regions_dict

        # select edges
        for edge_counter, edge in enumerate(self.collection.edges_sorted('var', reverse=True, lazy=lazy)):
            if edge_counter >= sample_size:
                # raise StopIteration
                break

            weights = []
            for i in range(n_hic):
                weights.append(getattr(edge, 'weight_' + str(i)))

            # require at least require_enriched weights to be fold_change
            # higher than each of the remaining require_nonenriched samples
            n_enriched = 0
            for i in range(len(weights)):
                local_enriched = 0
                for j in range(len(weights)):
                    if i == j:
                        continue
                    if weights[i] >= fold_change*weights[j]:
                        local_enriched += 1
                if local_enriched >= require_nonenriched:
                    n_enriched += 1

            if n_enriched >= require_enriched:
                if regions is None:
                    yield edge
                else:
                    source_region = regions_dict[edge.source]
                    sink_region = regions_dict[edge.sink]
                    if source_region.overlaps(regions) and sink_region.overlaps(regions):
                        yield edge


def do_pca(hics, pair_selection=None, tmpdir=None, eo_cutoff=0.0, bg_cutoff=1.0,
           log=True, only_intra_chromosomal=False, scale_libraries=True,
           regions=None, **kwargs):
    if pair_selection is None:
        pair_selection = LargestVariancePairSelection()

    pair_selection.regions = regions

    if tmpdir is not None:
        tmpdir = tempfile.mkdtemp(dir=os.path.expanduser(tmpdir))
    else:
        tmpdir = tempfile.mkdtemp()
    if not tmpdir.endswith('/'):
        tmpdir += '/'

    logger.info("Joining objects")
    existing_coll = False
    if isinstance(hics, HicCollectionWeightMeanVariance):
        logger.info("Found collection.")
        coll = hics
        n_hics = 0
        while 'weight_' + str(n_hics) in coll.field_names:
            n_hics += 1
        existing_coll = True
    else:
        n_hics = len(hics)
        coll = HicCollectionWeightMeanVariance(hics, file_name=tmpdir + 'coll.m',
                                               only_intra_chromosomal=only_intra_chromosomal,
                                               scale_libraries=scale_libraries)
        coll.calculate()

        if eo_cutoff != 0.0:
            eof = ExpectedObservedCollectionFilter(coll)
            coll.filter(eof, queue=True)

        if bg_cutoff != 1.0:
            bgf = BackgroundLigationCollectionFilter(coll)
            coll.filter(bgf, queue=True)
        coll.run_queued_filters()

    pair_selection.set_collection(coll)

    values = []
    for edge in pair_selection.pair_selection(**kwargs):
        weights = []
        for i in range(n_hics):
            weights.append(getattr(edge, 'weight_' + str(i)))
        values.append(weights)

    if log:
        y = np.array(values)
    else:
        y = np.log(np.array(values))

    pca = PCA()
    pca_res = pca.fit_transform(y.T)
    logger.info("Variance explained: %s" % str(pca.explained_variance_ratio_))

    if not existing_coll:
        coll.close()
    return pca, pca_res
