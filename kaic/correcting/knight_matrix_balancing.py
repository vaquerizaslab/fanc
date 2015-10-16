import logging
import warnings
import numpy as np
from kaic.tools.matrix import remove_sparse_rows, restore_sparse_rows


def correct(hic, only_intra_chromosomal=False):
    if only_intra_chromosomal:
        for chromosome in hic.chromosomes():
            m = hic[chromosome, chromosome]
            m_corrected = correct_matrix(m)
            hic[chromosome, chromosome] = m_corrected
    else:
        m = hic[:, :]
        m_corrected = correct_matrix(m)
        hic[:, :] = m_corrected


def correct_matrix(m, max_attempts=50):
    # remove zero-sum rows
    removed_rows = []
    m_nonzero, ixs = remove_sparse_rows(m, cutoff=0)
    removed_rows.append(ixs)

    has_errors = True
    iterations = 0
    while has_errors:
        has_errors = False

        try:
            x = get_bias_vector(m_nonzero)
        except ValueError, e:
            logging.info("Matrix balancing failed (this can happen!), \
                          removing sparsest rows to try again. Error: \
                          %s" % str(e))
            m_nonzero, ixs = remove_sparse_rows(m_nonzero)
            removed_rows.append(ixs)
            has_errors = True

        iterations += 1
        if iterations > max_attempts:
            raise RuntimeError("Exceeded maximum attempts (%d)" % max_attempts)

    for i in range(0, m_nonzero.shape[0]):
        for j in range(0, m_nonzero.shape[1]):
            m_nonzero[i, j] = x[i]*m_nonzero[i, j]*x[j]
    
    # restore zero rows
    for idx in reversed(removed_rows):
        m_nonzero = restore_sparse_rows(m_nonzero, idx)
        x = restore_sparse_rows(x, idx)

    return m_nonzero


def get_bias_vector(A, x0=None, tol=1e-06, delta=0.1, Delta=3, fl=0, high_precision=False, outer_limit=300):
    logging.info("Starting matrix balancing")

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
                    raise StopIteration("Number of iterations has exceeded the limit (%d)." % outer_limit)

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
                    print "%d %d   %.3e %.3e %.3e \n" % (i, k, res_norm, min(y), min(x))
                    res = np.array([[res], [res_norm]])

            logging.info("Matrix-vector products = %d\n" % MVP)
            logging.info("Outer iterations: %d" % n_iterations_outer)
        except Warning, e:
            logging.error(str(e))
            raise ValueError("Generic catch all warnings")

    return x

