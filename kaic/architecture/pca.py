from __future__ import division
from kaic.data.genomic import Edge
from kaic.architecture.hic_architecture import MatrixArchitecturalRegionFeature, \
    BackgroundLigationFilter, ExpectedObservedEnrichmentFilter, ZeroWeightFilter, \
    ExpectedObservedCollectionFilter, BackgroundLigationCollectionFilter, \
    HicEdgeCollection
import tables as t
import numpy as np
import tempfile
import os.path
import logging
logging.basicConfig(level=logging.INFO)


class HicCollectionWeightMeanVariance(HicEdgeCollection):
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
            for i in xrange(len(self.hics)):
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
                for i in xrange(len(self.hics)):
                    weight = getattr(edge, 'weight_' + str(i))
                    setattr(edge, 'weight_' + str(i), weight/weight_ratios[i])
        self.flush()

        self._edges.cols.var.create_csindex()


def do_pca(hics, tmpdir=None, sample_size=10000):
    if tmpdir is not None:
        tmpdir = tempfile.mkdtemp(dir=os.path.expanduser(tmpdir))
    else:
        tmpdir = tempfile.mkdtemp()
    if not tmpdir.endswith('/'):
        tmpdir += '/'

    logging.info("Joining objects")
    coll = HicCollectionWeightMeanVariance(hics, file_name=tmpdir + 'coll.m')
    coll.calculate()

    eof = ExpectedObservedCollectionFilter(coll)
    bgf = BackgroundLigationCollectionFilter(coll)
    coll.filter(eof, queue=True)
    coll.filter(bgf, queue=True)
    coll.run_queued_filters()

    # filter "uninteresting" contacts
    logging.info("Filtering uninteresting contacts")
    for hic in hics_copy:
        print len(hic.edges)
        blf = BackgroundLigationFilter(hic, fold_change=5)
        eof = ExpectedObservedEnrichmentFilter(hic, fold_change=1.5)
        hic.filter(blf, queue=True)
        hic.filter(eof, queue=True)
        hic.run_queued_filters()
        print len(hic.edges)

    # join into one object
    logging.info("Joining objects")
    coll = HicWeightMeanVariance(hics, file_name=tmpdir + 'coll.m')
    zwf = ZeroWeightFilter()
    logging.info("Filtering zero-weights")
    coll.filter(zwf)

    values = []
    current = 0
    for edge in coll.edges_sorted('var', reverse=True, lazy=True):
        if np.isnan(edge.var):
            continue

        current += 1
        if current >= sample_size:
            break

        weights = []
        for i in xrange(0, len(hics)):
            weights.append(edge['weight_' + str(i)])
        values.append(weights)
