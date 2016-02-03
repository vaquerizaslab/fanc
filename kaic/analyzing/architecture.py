from __future__ import division
from kaic.data.genomic import Hic
from kaic.plotting.plotter import GenomicTrack
import numpy as np
import logging
import itertools as it
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

# class HicArchitecture(object):
#     def __init__(self, hic):
#         self.hic = hic

def insulation_index(hic, d, mask_thresh=.5, hic_matrix=None, aggr_func=np.ma.median):
    chr_bins = hic.chromosome_bins
    n = len(hic.regions())
    if hic_matrix is None:
        log.debug("Fetching matrix")
        hic_matrix = hic[:,:]
    hic_matrix = hic_matrix.masked_matrix
    ins_matrix = np.empty(n)
    log.debug("Starting processing")
    skipped = 0
    for i, r in enumerate(hic.regions()):
        if (i - chr_bins[r.chromosome][0] < d or
            chr_bins[r.chromosome][1] - i < d):
            ins_matrix[i] = np.nan
            continue
        if np.sum(hic_matrix.mask[i: i + d, i - d:i]) > d*d*mask_thresh:
            # If too close to the edge of chromosome or
            # if more than half of the entries in this quadrant are masked (unmappable)
            # exclude it from the analysis
            skipped += 1
            ins_matrix[i] = np.nan
            continue
        ins_matrix[i] = aggr_func(hic_matrix[i: i + d, i - d:i])
    log.info("Skipped {} regions because >{:.1%} of matrix positions were masked".format(skipped, mask_thresh))
    return ins_matrix

def rel_insulation_index(hic, d, mask_thresh=.5, hic_matrix=None, aggr_func=np.ma.median):
    chr_bins = hic.chromosome_bins
    n = len(hic.regions())
    if hic_matrix is None:
        log.debug("Fetching matrix")
        hic_matrix = hic[:,:]
    hic_matrix = hic_matrix.masked_matrix
    rel_ins_matrix = np.empty(n)
    log.debug("Starting processing")
    skipped = 0
    for i, r in enumerate(hic.regions()):
        if (i - chr_bins[r.chromosome][0] < d or
            chr_bins[r.chromosome][1] - i < d):
            rel_ins_matrix[i] = np.nan
            continue
        if np.sum(hic_matrix.mask[i - d: i + d, i - d:i + d]) > 4*d*d*mask_thresh:
            # If too close to the edge of chromosome or
            # if more than half of the entries in this quadrant are masked (unmappable)
            # exclude it from the analysis
            skipped += 1
            rel_ins_matrix[i] = np.nan
            continue
        rel_ins_matrix[i] = (aggr_func(hic_matrix[i: i + d, i - d:i]) /
            aggr_func([hic_matrix[i - d:i, i - d:i], hic_matrix[i:i + d, i:i + d]]))
    log.info("Skipped {} regions because >{:.1%} of matrix positions were masked".format(skipped, mask_thresh))
    return rel_ins_matrix

def contact_band(hic, d1, d2, mask_thresh=.5, hic_matrix=None, use_oe_ratio=False, aggr_func=np.ma.median):
    chr_bins = hic.chromosome_bins
    n = len(hic.regions())
    if hic_matrix is None:
        log.debug("Fetching matrix")
        hic_matrix = hic[:,:]
    hic_matrix = hic_matrix.masked_matrix
    band = np.empty(n)
    if use_oe_ratio:
        raise NotImplementedError("oe not implemented yet :(")
    log.debug("Starting processing")
    skipped = 0
    for i, r in enumerate(hic.regions()):
        if (i - chr_bins[r.chromosome][0] < d2 or
            chr_bins[r.chromosome][1] - i < d2):
            band[i] = np.nan
            continue
        if np.sum(hic_matrix.mask[i - d2:i - d1, i + d1:i + d2]) > ((d2 - d1)**2)*mask_thresh:
            # If too close to the edge of chromosome or
            # if more than half of the entries in this quadrant are masked (unmappable)
            # exclude it from the analysis
            skipped += 1
            band[i] = np.nan
            continue
        band[i] = aggr_func(hic_matrix[i - d2:i - d1, i + d1:i + d2])
    log.info("Skipped {} regions because >{:.1%} of matrix positions were masked".format(skipped, mask_thresh))
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

def create_gtf_from_region_ix(regions, region_ix):
    for s, e in region_ix:
        yield "{}\tkaic\tboundary\t{}\t{}\t.\t.\t.\n".format(
            regions[s].chromosome, regions[s].start, regions[e].end)

class PeakCallerDeriv(object):
    def __init__(self, x, window_size=15, poly_order=2):
        self.x = x
        self.window_size = window_size
        self.poly_order = poly_order
        self._call_peaks()

    def _call_peaks(self):
        # Calculate Savitzky-Golay smoothed x and 1st and 2nd order derivative of input
        self.x_smooth = savgol_filter(self.x, window_length=self.window_size, polyorder=self.poly_order, deriv=0)
        self.d1 = savgol_filter(self.x, window_length=self.window_size, polyorder=self.poly_order, deriv=1)
        self.d2 = savgol_filter(self.x, window_length=self.window_size, polyorder=self.poly_order, deriv=2)
        # Find indices of zero crossings of 1st derivative
        self._peaks = np.nonzero(np.diff(np.signbit(self.d1)))[0]
        log.info("Found {} raw peaks".format(len(self._peaks)))
        # Minima have negative 2nd order derivative at minimum
        self._min_mask = self.d2[self._peaks] > 0
        # Boundaries of peaks can be defined as the inflection point up-
        # and downstream of the peak. (zero crossings of 2nd derivative)
        self.infl_points = np.nonzero(np.diff(np.signbit(self.d2)))[0]
        infl_ix = np.searchsorted(self.infl_points, self._peaks, side="right")
        self.peak_bounds = np.empty((len(self._peaks), 2), dtype=np.int_)
        # Have to catch boundary condition if inflection point before first peak
        # or after last peak is not known
        try:
            self.peak_bounds[:, 0] = self.infl_points[infl_ix - 1]
        except IndexError:
            infl_ix[0] = 1
            self.peak_bounds[:, 0] = self.infl_points[infl_ix - 1]
            self.peak_bounds[0][0] = self._peaks[0]
        try:
            self.peak_bounds[:, 1] = self.infl_points[infl_ix]
        except IndexError:
            infl_ix[-1] = 1
            self.peak_bounds[:, 1] = self.infl_points[infl_ix]
            self.peak_bounds[-1][1] = self._peaks[-1]
        peak_vals = self.x_smooth[self._peaks]
        # For each peak get average value of inflection points left and right
        # Maybe what we want instead is the inflection point that has largest distance to peak summit?
        ip_means = np.ma.mean([self.x_smooth[self.peak_bounds[:, 0]], self.x_smooth[self.peak_bounds[:, 1]]], axis=0)
        self._heights = peak_vals - ip_means
        self._peak_mask = np.full(self._peaks.shape, True, dtype=np.bool_)

    def get_peaks(self, as_range=False):
        if as_range:
            return self.peak_bounds[self._peak_mask]
        return self._peaks[self._peak_mask]

    def get_minima(self, as_range=False):
        if as_range:
            return self.peak_bounds[np.logical_and(self._peak_mask, self._min_mask)]
        return self._peaks[np.logical_and(self._peak_mask, self._min_mask)]

    def get_maxima(self, as_range=False):
        if as_range:
            return self.peak_bounds[np.logical_and(self._peak_mask, ~self._min_mask)]
        return self._peaks[np.logical_and(self._peak_mask, ~self._min_mask)]

    @property
    def heights(self):
        return self._heights[self._peak_mask]

    def filter(self, z_score_thresh=None, z_score_window=50, steepness_thresh=None, height_thresh=None, abs_thresh=None):
        def rolling_func_nan(a, window, func=np.mean):
            out = np.empty(a.shape)
            for i in range(len(a)):
                if i < window or i >= len(a) - window:
                    out[i] = np.nan
                    continue
                cur_window = a[i - window:i + window + 1]
                out[i] = func(cur_window[~np.isnan(cur_window)])
            return out
        if abs_thresh:
            abs_pass = np.where(self._min_mask, self.x_smooth < abs_thresh[0], self.x_smooth > abs_thresh[1])
            self._peak_mask = np.logical_and(self._peak_mask, abs_pass)
            log.info("Discarding {}({:.1%}) of total peaks due to absolute value threshold ({})".format(np.sum(~abs_pass), np.sum(~abs_pass)/len(self._peaks), abs_thresh))
        if z_score_thresh:
            rolling_mean = rolling_func_nan(self.x, z_score_window, func=np.mean)
            z_trans = (self.x - rolling_mean)/np.std(self.x[~np.isnan(self.x)])
            z_pass = np.logical_or(z_trans[self._peaks] < -z_score_thresh, z_trans[self._peaks] > z_score_thresh)
            self._peak_mask = np.logical_and(self._peak_mask, z_pass)
            log.info("Discarding {}({:.1%}) of total peaks due to z-score threshold ({})".format(np.sum(~z_pass), np.sum(~z_pass)/len(self._peaks), z_score_thresh))
        if steepness_thresh:
            # If threshold is set, only retain peaks where 2nd derivative
            # is above the required percentile for steepness
            #d2_thresh = np.nanpercentile(np.abs(self.d2[self._peaks]), 100 - steepness_thresh)
            steep_pass = np.abs(self.d2[zc1]) > steepness_thresh
            self._peak_mask = np.logical_and([self._peak_mask, steep_pass], axis=0)
            log.info("Discarding {}({:.1%}) of total peaks due to steepness threshold ({})".format(np.sum(~steep_pass), np.sum(~steep_pass)/len(self._peaks), steepness_thresh))
        if height_thresh:
            #height_thresh_minima = np.nanpercentile(self._heights[self._min_mask], height_thresh)
            #height_thresh_maxima = np.nanpercentile(self._heights[~self._min_mask], 100 - height_thresh)
            #height_pass = np.where(self._min_mask, self._heights <= height_thresh_minima, self._heights >= height_thresh_maxima)
            height_pass = np.abs(self._heights[self._peaks]) > height_thresh
            self._peak_mask = np.logical_and(self._peak_mask, height_pass)
            log.info("Discarding {}({:.1%}) of total peaks due to peak height threshold ({})".format(np.sum(~height_pass), np.sum(~height_pass)/len(self._peaks), height_thresh))
