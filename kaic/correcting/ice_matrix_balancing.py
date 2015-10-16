import numpy as np
import logging


def correct(hic, tolerance=1e-2, max_iterations=500):
    bias_vector = np.ones(len(hic.regions()), float)
    marginal_error = tolerance + 1
    current_iteration = 0
    logging.info("Starting iterations")
    while (marginal_error > tolerance
           and current_iteration <= max_iterations):
        m = get_marginals(hic, normalize_to_mean_one=True)
        bias_vector *= m
        marginal_error = _marginal_error(m)
        for row in hic._edges._iter_visible_and_masked():
            source = row['source']
            sink = row['sink']
            weight = row['weight']
            row['weight'] = weight/m[source]/m[sink]
            row.update()
        hic.flush()
        current_iteration += 1
        logging.info("Iteration: %d, error: %lf" % (current_iteration, marginal_error))
    hic.bias_vector(bias_vector)

def get_marginals(hic, normalize_to_mean_one=False):
    # prepare marginals dict
    marginals = np.zeros(len(hic.regions()), float)
    
    for edge in hic.edges():
        marginals[edge.source] += edge.weight
        if edge.source != edge.sink:
            marginals[edge.sink] += edge.weight
        
    if normalize_to_mean_one:
        divisor = marginals[marginals > 0].mean()
        marginals /= divisor
        marginals[marginals == 0] = 1
        marginals -= 1
        marginals *= .6
        marginals += 1
    
    return marginals


def _marginal_error(marginals):
    marginals = marginals[marginals != 0]
    error = np.max(np.abs(marginals - marginals.mean()))
    return error / marginals.mean()
