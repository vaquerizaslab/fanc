import numpy as np
from scipy.stats.mstats import gmean


def remove_sparse_rows(m, cutoff=None):
    s = np.sum(m, 0)
    
    if cutoff is None:
        cutoff = min(s)
    
    idxs = np.where(s <= cutoff)[0]
    m_removed = np.delete(m, idxs, 0)
    m_removed = np.delete(m_removed, idxs, 1)
    
    return m_removed, idxs
    
    
def restore_sparse_rows(m, idx_sets, rows=None):
    abs_idx = []
    for idxs in reversed(idx_sets):
        for i in sorted(idxs):
            shift = 0
            for j in sorted(abs_idx):
                if j + shift < i:
                    shift += 1
            abs_idx.append(i - shift)
    abs_idx.sort()
    a = np.insert(m, abs_idx, 0, axis=0)
    if len(m.shape) > 1:
        a = np.insert(a, abs_idx, 0, axis=1)
    return a


def is_symmetric(m, tol=1e-10):
    for i in range(0, m.shape[0]):
        for j in range(i, m.shape[1]):
            if abs(m[i, j]-m[j, i]) > tol:
                return False
    return True


def apply_sliding_func(a, window, func=np.ma.mean):
    """
    Apply function on a sliding window over an array, ignoring Numpy NaN values.

    :param a: Numpy array on which function is applied
    :param window: The sliding window is i - window:i + window + 1
                   so total window is twice this parameter.
    :param func: Function to apply
    """
    out = np.empty(a.shape)
    for i in range(len(a)):
        window_start = max(0, i - window)
        window_end = min(len(a), i + window + 1)
        cur_window = a[window_start:window_end]
        out[i] = func(cur_window[~np.isnan(cur_window)])
    return out


def delta_window(x, window, ignore_mask=False, mask_thresh=.5):
    try:
        x.mask = np.logical_or(x.mask, ~np.isfinite(x))
    except AttributeError:
        x = np.ma.masked_invalid(x)
    n = len(x)
    delta = np.empty(n)
    for i in range(n):
        if i < window or n - i <= window - 1:
            delta[i] = np.nan
            continue
        down_slice = slice(i + 1, i + window + 1)
        up_slice = slice(i - window, i)
        if not ignore_mask and (np.sum(x.mask[down_slice]) > window*mask_thresh or
                np.sum(x.mask[up_slice]) > window*mask_thresh):
            delta[i] = np.nan
        delta[i] = np.ma.mean(x[down_slice] - x[i]) - np.ma.mean(x[up_slice] - x[i])
    return delta


def kth_diag_indices(n, k):
    """
    Return indices of bins k steps away from the diagonal.
    (from http://stackoverflow.com/questions/10925671/numpy-k-th-diagonal-indices)
    """

    rows, cols = np.diag_indices(n)
    if k < 0:
        return rows[-k:], cols[:k]
    elif k > 0:
        return rows[:-k], cols[k:]
    else:
        return rows, cols


def trim_stats(a, proportiontocut=0.0, axis=0, stat=np.nanmean):
    if proportiontocut >= 0.5:
        raise ValueError("Cannot cut more than 50% of values off distribution tails!")

    s = np.sort(a, axis=axis)
    ix = int(s.shape[axis] * proportiontocut)

    sl = [slice(None)] * s.ndim
    sl[axis] = slice(ix, len(s)-ix, 1)
    s_sub = s[tuple(sl)]

    return stat(s_sub)


def nangmean(a, axis=0, dtype=None):
    return gmean(a[np.isfinite(a)], axis=axis, dtype=dtype)
