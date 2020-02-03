import warnings
import numpy as np
from fanc.tools.matrix import remove_sparse_rows, restore_sparse_rows
from fanc.legacy.data.genomic import Hic, AccessOptimisedHic, Edge
import logging
logger = logging.getLogger(__name__)


def correct(hic, whole_matrix=True, intra_chromosomal=True, inter_chromosomal=True,
            copy=False, file_name=None, optimise=True,
            restore_coverage=False):
    hic_new = None
    chromosome_starts = dict()
    last_chromosome = None
    if copy:
        if optimise:
            hic_new = AccessOptimisedHic(file_name=file_name, mode='w')
        else:
            hic_new = Hic(file_name=file_name, mode='w')

        for i, region in enumerate(hic.regions()):
            hic_new.add_region(region, flush=False)
            if region.chromosome != last_chromosome:
                chromosome_starts[region.chromosome] = i
            last_chromosome = region.chromosome
        hic_new.flush()
        hic_new.disable_indexes()

    if not whole_matrix:
        bias_vectors = []
        for chromosome in hic.chromosomes():
            m = hic[chromosome, chromosome]
            m_corrected, bias_vector_chromosome = correct_matrix(m, restore_coverage=restore_coverage)
            if hic_new is None:
                logger.debug("Replacing corrected edges in existing Hic object...")
                hic[chromosome, chromosome] = m_corrected
            else:
                logger.debug("Adding corrected edges in new Hic object ...")
                chromosome_offset = chromosome_starts[chromosome]
                for i in range(m_corrected.shape[0]):
                    i_region = i + chromosome_offset
                    nonzero_idx = np.nonzero(m_corrected[i])[0]
                    for j in nonzero_idx[nonzero_idx >= i]:
                        j_region = j + chromosome_offset
                        weight = m_corrected[i, j]
                        hic_new.add_edge(Edge(source=i_region, sink=j_region, weight=weight), flush=False)
            bias_vectors.append(bias_vector_chromosome)
        logger.debug("Done.")
        logger.debug("Adding bias vector...")
        if hic_new is None:
            hic.bias_vector(np.concatenate(bias_vectors))
        else:
            hic_new.bias_vector(np.concatenate(bias_vectors))
        logger.debug("Done.")
    else:
        logger.debug("Fetching whole genome matrix")
        m = hic[:, :]
        cb = hic.chromosome_bins
        if not intra_chromosomal:
            for chromosome, bins in cb.items():
                m[bins[0]:bins[1], bins[0]:bins[1]] = 0
        if not inter_chromosomal:
            for chromosome, bins in cb:
                m[0:bins[0], bins[0]:bins[1]] = 0
                m[bins[1]:, bins[0]:bins[1]] = 0

        m_corrected, bias_vector = correct_matrix(m, restore_coverage=restore_coverage)
        if hic_new is None:
            logger.debug("Replacing corrected edges in existing Hic object...")
            hic[:, :] = m_corrected
        else:
            logger.debug("Adding corrected edges in new Hic object ...")
            for i in range(m_corrected.shape[0]):
                nonzero_idx = np.nonzero(m_corrected[i])[0]
                for j in nonzero_idx[nonzero_idx >= i]:
                    weight = m_corrected[i, j]
                    hic_new.add_edge([i, j, weight], flush=False)
        logger.debug("Done.")
        logger.debug("Adding bias vector...")
        if hic_new is None:
            hic.bias_vector(bias_vector)
        else:
            hic_new.bias_vector(bias_vector)
        logger.debug("Done.")

    if hic_new is None:
        return hic
    hic_new.flush()
    hic_new.enable_indexes()
    return hic_new


def correct_matrix(m, max_attempts=50, restore_coverage=False):
    # remove zero-sum rows
    removed_rows = []
    m_nonzero, ixs = remove_sparse_rows(m, cutoff=0)
    removed_rows.append(ixs)

    has_errors = True
    iterations = 0
    x = None
    while has_errors:
        has_errors = False

        try:
            x = get_bias_vector(m_nonzero)
        except ValueError as e:
            logger.debug("Matrix balancing failed (this can happen!), \
                          removing sparsest rows to try again. Error: \
                          %s" % str(e))
            m_nonzero, ixs = remove_sparse_rows(m_nonzero)
            removed_rows.append(ixs)
            has_errors = True

        iterations += 1
        if iterations > max_attempts:
            raise RuntimeError("Exceeded maximum attempts (%d)" % max_attempts)

    if restore_coverage:
        x = x*np.sqrt(np.sum(m_nonzero)/m_nonzero.shape[0])

    logger.debug("Applying bias vector")
    m_nonzero = x*m_nonzero*x[:, np.newaxis]

    logger.debug(removed_rows)
    logger.debug("Restoring {} sets ({} total) sparse rows.".format(
        len(removed_rows), sum(len(x) for x in removed_rows)))
    # restore zero rows
    m_nonzero = restore_sparse_rows(m_nonzero, removed_rows)
    x = restore_sparse_rows(x, removed_rows)

    return m_nonzero, x


def get_bias_vector(A, x0=None, tol=1e-06, delta=0.1, Delta=3, fl=0, high_precision=False, outer_limit=300):
    logger.debug("Starting matrix balancing")

    with warnings.catch_warnings():
        warnings.filterwarnings('error')

        try:
            # basic variables
            # n=size_(A,1)
            if not isinstance(A, np.ndarray):
                try:
                    if high_precision:
                        A = np.array(A, dtype=np.float128)
                    else:
                        A = np.array(A)
                except AttributeError:
                    A = np.array(A)
            n = A.shape[0]
            # e=ones_(n,1)
            try:
                if high_precision:
                    e = np.ones(n, dtype=np.float128)
                else:
                    e = np.ones(n)
            except AttributeError:
                e = np.ones(n)

            if not x0:
                try:
                    if high_precision:
                        x0 = np.ones(n, np.float128)
                    else:
                        x0 = np.ones(n)
                except AttributeError:
                    x0 = np.ones(n)
            else:
                try:
                    if high_precision:
                        x0 = np.array(x0, np.float128)
                    else:
                        x0 = np.array(x0)
                except AttributeError:
                    x0 = np.array(x0)
            res = np.array([])
            g = 0.9
            etamax = 0.1
            eta = etamax
            stop_tol = tol * 0.5
            x = x0.copy()
            rt = tol ** 2
            # v=x.dot((A * x))
            v = x*A.dot(x)
            rk = 1 - v
            # rho_km1=rk.T * rk
            rho_km1 = rk.T.dot(rk)
            rout = rho_km1
            rho_km2 = rho_km1
            rold = rout
            MVP = 0
            i = 0

            n_iterations_outer = 0
            while rout > rt:
                n_iterations_outer += 1

                if n_iterations_outer > outer_limit:
                    raise ValueError("Number of iterations has exceeded the limit (%d)." % outer_limit)

                i += 1
                k = 0
                y = e.copy()
                innertol = max(eta ** 2 * rout, rt)
                n_iterations_inner = 0
                while rho_km1 > innertol:
                    n_iterations_inner += 1

                    k += 1
                    if k == 1:
                        try:
                            Z = rk / v
                        except Warning:
                            raise ValueError("v=0; Remove zero or sparse rows")
                        p = Z.copy()
                        rho_km1 = rk.T.dot(Z)
                    else:
                        beta = rho_km1 / rho_km2
                        p = Z + beta * p
                    # w = x.*(A*(x.*p)) + v.*p;
                    w = x*A.dot(x*p) + v*p
                    alpha = rho_km1 / p.T.dot(w)
                    ap = alpha * p
                    ynew = y + ap
                    if min(ynew) <= delta:
                        if delta == 0:
                            break
                        ind = np.where(ap < 0)[0]
                        # gamma = min((delta  - y(ind))./ap(ind));
                        gamma = min((delta-y[ind])/ap[ind])
                        y = y + gamma * ap
                        break
                    if max(ynew) >= Delta:
                        ind = np.where(ynew > Delta)[0]
                        gamma = min((Delta-y[ind])/ap[ind])
                        y = y + gamma * ap
                        break
                    y = ynew.copy()
                    rk = rk - alpha * w
                    rho_km2 = rho_km1.copy()
                    Z = rk / v
                    rho_km1 = rk.T.dot(Z)

                try:
                    x = x*y
                except Warning:
                    raise ValueError("Value in x or y too small to represent numerically. Try removing sparse rows")

                v = x*A.dot(x)
                rk = 1 - v
                rho_km1 = rk.T.dot(rk)
                rout = rho_km1.copy()
                MVP = MVP + k + 1
                rat = rout / rold
                rold = rout.copy()
                res_norm = np.sqrt(rout)
                eta_o = eta
                eta = g * rat
                if g * eta_o ** 2 > 0.1:
                    eta = max(eta, g * eta_o ** 2)
                eta = max(min(eta, etamax), stop_tol / res_norm)
                if fl == 1:
                    res = np.array([[res], [res_norm]])

            logger.debug("Matrix-vector products = %d\n" % MVP)
            logger.debug("Outer iterations: %d" % n_iterations_outer)
        except Warning as e:
            logger.error(str(e))
            raise ValueError("Generic catch all warnings")

    return x

