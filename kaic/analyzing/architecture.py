from kaic.data.genomic import Hic
import numpy as np
import logging
import ipdb
from scipy.signal import savgol_filter

log = logging.getLogger(__name__)
log.setLevel(10)

def kth_diag_indices(n, k):
    #http://stackoverflow.com/questions/10925671/numpy-k-th-diagonal-indices
    rows, cols = np.diag_indices(n)
    if k < 0:
        return rows[:k], cols[-k:]
    elif k > 0:
        return rows[k:], cols[:-k]
    else:
        return rows, cols

def insulation_index(hic, d, hic_matrix=None):
    chr_bins = hic.chromosome_bins
    n = len(hic.regions())
    if hic_matrix is None:
        log.debug("Fetching matrix")
        hic_matrix = hic[:,:]
    ins_matrix = np.empty(n)
    log.debug("Starting processing")
    for i, r in enumerate(hic.regions()):
        if i - chr_bins[r.chromosome][0] < d:
            ins_matrix[i] = np.nan
            continue
        if chr_bins[r.chromosome][1] - i < d:
            ins_matrix[i] = np.nan
            continue
        ins_matrix[i] = np.ma.sum(hic_matrix[i: i + d, i - d:i])
    return ins_matrix

def rel_insulation_index(hic, d, hic_matrix=None):
    chr_bins = hic.chromosome_bins
    n = len(hic.regions())
    if hic_matrix is None:
        log.debug("Fetching matrix")
        hic_matrix = hic[:,:]
    rel_ins_matrix = np.empty(n)
    log.debug("Starting processing")
    for i, r in enumerate(hic.regions()):
        if i - chr_bins[r.chromosome][0] < d:
            rel_ins_matrix[i] = np.nan
            continue
        if chr_bins[r.chromosome][1] - i < d:
            rel_ins_matrix[i] = np.nan
            continue
        rel_ins_matrix[i] = ((np.ma.sum(hic_matrix[i - d:i, i - d:i]) +
            np.ma.sum(hic_matrix[i:i + d, i:i + d])) /
            (2*np.ma.sum(hic_matrix[i: i + d, i - d:i])))
    return rel_ins_matrix

def contact_band(hic, d1, d2, hic_matrix=None, use_oe_ratio=False):
    chr_bins = hic.chromosome_bins
    n = len(hic.regions())
    if hic_matrix is None:
        log.debug("Fetching matrix")
        hic_matrix = hic[:,:]
    band = np.empty(n)
    if use_oe_ratio:
        raise NotImplementedError("oe not implemented yet :(")
    log.debug("Starting processing")
    for i, r in enumerate(hic.regions()):
        if i - chr_bins[r.chromosome][0] < d2:
            band[i] = np.nan
            continue
        if chr_bins[r.chromosome][1] - i < d2:
            band[i] = np.nan
            continue
        band[i] = np.ma.sum(hic_matrix[i - d2:i - d1, i + d1:i + d2])
    return band

def observed_expected_ratio(hic, hic_matrix=None, per_chromosome=True):
    chr_bins = hic.chromosome_bins
    n = len(hic.regions())
    if hic_matrix is None:
        log.debug("Fetching matrix")
        hic_matrix = hic[:,:]
    oe = hic_matrix.copy()
    log.debug("starting processing")
    if per_chromosome:
        for so, eo in chr_bins.itervalues():
            # Correcting intrachromosomal contacts by mean contact count at each diagonal
            for i in range(n):
                oe[so:eo, so:eo][kth_diag_indices(eo - so, -i)] /= np.ma.mean(hic_matrix[so:eo, so:eo][kth_diag_indices(eo - so, -i)])
            # Correcting interchromosomal contacts by mean of all contact counts between
            # each set of chromosomes
            for si, ei in chr_bins.itervalues():
                # Only correct upper triangle
                if si <= so:
                    continue
                oe[so:eo, si:ei] /= np.ma.mean(hic_matrix[so:eo, si:ei])
    else:
        for i in range(n):
            oe[kth_diag_indices(n, i)] /= np.ma.mean(hic_matrix.diagonal(i))
    # Copying upper triangle to lower triangle
    oe[np.tril_indices(n)] = oe.T[np.tril_indices(n)]
    return oe

def call_peaks_deriv(x, window_length=9, polyorder=2, as_range=True, thresh=None):
    # Calculate Savitzky-Golay smoothed 1st and 2nd order derivative of input
    d1 = savgol_filter(x=x, window_length=window_length, polyorder=polyorder, deriv=1)
    d2 = savgol_filter(x=x, window_length=window_length, polyorder=polyorder, deriv=2)
    # Find indices of zero crossings of 1st derivative
    zc1 = np.nonzero(np.diff(np.signbit(d1)))[0]
    if thresh:
        # If threshold is set, only retain peaks where 2nd derivative
        # is above threshold (sharpness of peak)
        is_above = np.abs(d2[zc1]) > thresh
        zc1 = zc1[is_above]
    # Minima have negative 2nd order derivative at minimum
    is_minimum = d2[zc1] > 0
    if as_range:
        # Boundaries of peaks can be defined as the inflection point up-
        # and downstream of the peak. (zero crossings of 2nd derivative)
        zc2 = np.nonzero(np.diff(np.signbit(d2)))[0]
        boundary_ix = np.searchsorted(zc2, zc1, side="right")
        extrema = np.empty((len(zc1), 2), dtype=np.int_)
        # Have to catch boundary condition if inflection point before first peak
        # or after last peak is not known
        try:
            extrema[:, 0] = zc2[boundary_ix - 1]
        except IndexError:
            boundary_ix[0] = 1
            extrema[:, 0] = zc2[boundary_ix - 1]
            extrema[0][0] = zc1[0]
        try:
            extrema[:, 1] = zc2[boundary_ix]
        except IndexError:
            boundary_ix[-1] = 1
            extrema[:, 1] = zc2[boundary_ix]
            extrema[-1][1] = zc1[-1]
        return (extrema[is_minimum], extrema[~is_minimum])
    return (zc1[is_minimum], zc1[~is_minimum])
