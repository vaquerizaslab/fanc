from fanc.legacy.data.genomic import AccessOptimisedHic, Hic
import numpy as np
from collections import defaultdict
import logging
logger = logging.getLogger(__name__)


def correct(hic, tolerance=1e-2, max_iterations=500, whole_matrix=True,
            inter_chromosomal=True, intra_chromosomal=True, copy=False, file_name=None,
            optimise=True):
    if copy:
        if optimise:
            hic_new = AccessOptimisedHic(file_name=file_name, mode='w')
        else:
            hic_new = Hic(file_name=file_name, mode='w')

        hic_new.add_regions(hic.regions())
        hic_new.flush()

        hic_new.disable_indexes()
        for edge in hic.edges(intra_chromosomal=intra_chromosomal,
                              inter_chromosomal=inter_chromosomal):
            hic_new.add_edge(edge, flush=False)
        hic_new.enable_indexes()
        hic_new.flush()
        hic = hic_new

    if not whole_matrix:
        bias_vectors = []
        for chromosome in hic.chromosomes():
            region_converter = dict()
            bias_vector = []
            for i, region in enumerate(hic.subset(chromosome)):
                region_converter[region.ix] = i
                bias_vector.append(1)
            bias_vector = np.array(bias_vector, dtype='float64')

            marginal_error = tolerance + 1
            current_iteration = 0
            logger.info("Starting iterations for chromosome {} ".format(chromosome))
            while (marginal_error > tolerance and
                   current_iteration <= max_iterations):
                m = np.zeros(len(bias_vector), dtype='float64')
                for edge in hic.edge_subset(key=(chromosome, chromosome), lazy=True):
                    source = region_converter[edge.source]
                    sink = region_converter[edge.sink]
                    m[source] += edge.weight
                    if source != sink:
                        m[sink] += edge.weight

                bias_vector *= m
                marginal_error = _marginal_error(m)
                for edge in hic.edge_subset(key=(chromosome, chromosome), lazy=True):
                    source = region_converter[edge.source]
                    sink = region_converter[edge.sink]
                    weight = edge.weight
                    edge.weight = 0 if m[sink] == 0 else weight / np.sqrt(m[source]) / np.sqrt(m[sink])
                hic.flush(silent=True)
                current_iteration += 1
                logger.info("Iteration: %d, error: %lf" % (current_iteration, marginal_error))
            bias_vectors.append(bias_vector)
        logger.info("Done.")
        logger.info("Adding bias vector...")
        hic.bias_vector(np.concatenate(bias_vectors))

        logger.info("Done.")
    else:
        bias_vector = np.ones(len(hic.regions), float)
        marginal_error = tolerance + 1
        current_iteration = 0
        logger.info("Starting iterations")
        while (marginal_error > tolerance and
               current_iteration <= max_iterations):
            m = hic.marginals()
            bias_vector *= m
            marginal_error = _marginal_error(m)
            for edge in hic.edges(lazy=True, intra_chromosomal=intra_chromosomal,
                                  inter_chromosomal=inter_chromosomal):
                source = edge.source
                sink = edge.sink
                weight = edge.weight
                edge.weight = 0 if m[sink] == 0 else weight/np.sqrt(m[source])/np.sqrt(m[sink])
            hic.flush()
            current_iteration += 1
            logger.info("Iteration: %d, error: %lf" % (current_iteration, marginal_error))
        hic.bias_vector(bias_vector)


def _marginal_error(marginals, percentile=99.9):
    marginals = marginals[marginals != 0]
    error = np.percentile(np.abs(marginals - marginals.mean()), percentile)
    return error / marginals.mean()
