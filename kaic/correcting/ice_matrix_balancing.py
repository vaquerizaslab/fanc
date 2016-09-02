import numpy as np
import logging
logger = logging.getLogger(__name__)


def correct(hic, tolerance=1e-2, max_iterations=500):
    bias_vector = np.ones(len(hic.regions()), float)
    marginal_error = tolerance + 1
    current_iteration = 0
    logger.info("Starting iterations")
    while (marginal_error > tolerance and
           current_iteration <= max_iterations):
        m = hic.marginals()
        bias_vector *= m
        marginal_error = _marginal_error(m)
        for edge in hic.edges(lazy=True):
            source = edge.source
            sink = edge.sink
            weight = edge.weight
            edge.weight = weight/np.sqrt(m[source])/np.sqrt(m[sink])
        hic.flush()
        current_iteration += 1
        logger.info("Iteration: %d, error: %lf" % (current_iteration, marginal_error))
    hic.bias_vector(bias_vector)


def _marginal_error(marginals, percentile=99.9):
    marginals = marginals[marginals != 0]
    error = np.percentile(np.abs(marginals - marginals.mean()), percentile)
    return error / marginals.mean()
