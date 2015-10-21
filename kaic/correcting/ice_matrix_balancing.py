import numpy as np
import logging


def correct(hic, tolerance=1e-2, max_iterations=500):
    bias_vector = np.ones(len(hic.regions()), float)
    marginal_error = tolerance + 1
    current_iteration = 0
    logging.info("Starting iterations")
    while (marginal_error > tolerance and
           current_iteration <= max_iterations):
        m = get_marginals(hic)
        bias_vector *= m
        marginal_error = _marginal_error(m)
        for row in hic._edges._iter_visible_and_masked():
            source = row['source']
            sink = row['sink']
            weight = row['weight']
            row['weight'] = weight/np.sqrt(m[source])/np.sqrt(m[sink])
            row.update()
        hic.flush()
        current_iteration += 1
        logging.info("Iteration: %d, error: %lf" % (current_iteration, marginal_error))
    hic.bias_vector(bias_vector)


def get_marginals(hic):
    # prepare marginals dict
    marginals = np.zeros(len(hic.regions()), float)
    
    for edge in hic.edges():
        marginals[edge.source] += edge.weight
        if edge.source != edge.sink:
            marginals[edge.sink] += edge.weight

    return marginals


def _marginal_error(marginals, percentile=99.9):
    marginals = marginals[marginals != 0]
    error = np.percentile(np.abs(marginals - marginals.mean()), percentile)
    return error / marginals.mean()
