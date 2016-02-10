from __future__ import division
from kaic.data.genomic import Hic
from kaic.plotting.plotter import GenomicTrack
import numpy as np
import logging
import itertools as it
from scipy.signal import savgol_filter

log = logging.getLogger(__name__)
log.setLevel(10)

def rolling_func_nan(a, window, func=np.ma.mean):
    out = np.empty(a.shape)
    for i in range(len(a)):
        window_start = max(0, i - window)
        window_end = min(len(a), i + window + 1)
        cur_window = a[window_start:window_end]
        out[i] = func(cur_window[~np.isnan(cur_window)])
    return out

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

def insulation_index(hic, d, mask_thresh=.5, hic_matrix=None, aggr_func=np.ma.mean):
    chr_bins = hic.chromosome_bins
    n = len(hic.regions())
    if hic_matrix is None:
        log.debug("Fetching matrix")
        hic_matrix = hic[:,:]
    if not hasattr(hic_matrix, "mask"):
        hic_matrix = hic_matrix.masked_matrix
    ins_matrix = np.empty(n)
    log.debug("Starting processing")
    skipped = 0
    for i, r in enumerate(hic.regions()):
        if (i - chr_bins[r.chromosome][0] < d or
            chr_bins[r.chromosome][1] - i <= d + 1):
            ins_matrix[i] = np.nan
            continue
        if hic_matrix.mask[i, i]:
            ins_matrix[i] = np.nan
            continue
        ins_slice = (slice(i + 1, i + d + 1), slice(i - d, i))
        if np.sum(hic_matrix.mask[ins_slice]) > d*d*mask_thresh:
            # If too close to the edge of chromosome or
            # if more than half of the entries in this quadrant are masked (unmappable)
            # exclude it from the analysis
            skipped += 1
            ins_matrix[i] = np.nan
            continue
        ins_matrix[i] = aggr_func(hic_matrix[ins_slice])
    log.info("Skipped {} regions because >{:.1%} of matrix positions were masked".format(skipped, mask_thresh))
    return ins_matrix

def rel_insulation_index(hic, d, mask_thresh=.5, hic_matrix=None, aggr_func=np.ma.mean):
    chr_bins = hic.chromosome_bins
    n = len(hic.regions())
    if hic_matrix is None:
        log.debug("Fetching matrix")
        hic_matrix = hic[:,:]
    if not hasattr(hic_matrix, "mask"):
        hic_matrix = hic_matrix.masked_matrix
    rel_ins_matrix = np.empty(n)
    log.debug("Starting processing")
    skipped = 0
    for i, r in enumerate(hic.regions()):
        if (i - chr_bins[r.chromosome][0] < d or
            chr_bins[r.chromosome][1] - i <= d + 1):
            rel_ins_matrix[i] = np.nan
            continue
        if hic_matrix.mask[i, i]:
            rel_ins_matrix[i] = np.nan
            continue
        up_rel_slice = (slice(i - d, i), slice(i - d, i))
        down_rel_slice = (slice(i + 1, i + d + 1), slice(i + 1, i + d + 1))
        ins_slice = (slice(i + 1, i + d + 1), slice(i - d, i))
        if (np.sum(hic_matrix.mask[up_rel_slice]) > d*d*mask_thresh or
            np.sum(hic_matrix.mask[down_rel_slice]) > d*d*mask_thresh or
            np.sum(hic_matrix.mask[ins_slice]) > d*d*mask_thresh):
            # If too close to the edge of chromosome or
            # if more than half of the entries in this quadrant are masked (unmappable)
            # exclude it from the analysis
            skipped += 1
            rel_ins_matrix[i] = np.nan
            continue
        rel_ins_matrix[i] = (aggr_func(hic_matrix[ins_slice]) /
            aggr_func(np.ma.dstack((hic_matrix[up_rel_slice], hic_matrix[down_rel_slice]))))
    log.info("Skipped {} regions because >{:.1%} of matrix positions were masked".format(skipped, mask_thresh))
    return rel_ins_matrix

def contact_band(hic, d1, d2, mask_thresh=.5, hic_matrix=None, use_oe_ratio=False, aggr_func=np.ma.mean):
    chr_bins = hic.chromosome_bins
    n = len(hic.regions())
    if hic_matrix is None:
        log.debug("Fetching matrix")
        hic_matrix = hic[:,:]
    if not hasattr(hic_matrix, "mask"):
        hic_matrix = hic_matrix.masked_matrix
    band = np.empty(n)
    if use_oe_ratio:
        raise NotImplementedError("oe not implemented yet :(")
    log.debug("Starting processing")
    skipped = 0
    for i, r in enumerate(hic.regions()):
        if (i - chr_bins[r.chromosome][0] < d2 or
            chr_bins[r.chromosome][1] - i <= d2 + 1):
            band[i] = np.nan
            continue
        if hic_matrix.mask[i, i]:
            band[i] = np.nan
            continue
        band_slice = (slice(i + d1 + 1, i + d2 + 1), slice(i - d2, i - d1))
        if np.sum(hic_matrix.mask[band_slice]) > ((d2 - d1)**2)*mask_thresh:
            # If too close to the edge of chromosome or
            # if more than half of the entries in this quadrant are masked (unmappable)
            # exclude it from the analysis
            skipped += 1
            band[i] = np.nan
            continue
        band[i] = aggr_func(hic_matrix[band_slice])
    log.info("Skipped {} regions because >{:.1%} of matrix positions were masked".format(skipped, mask_thresh))
    return band

def directionality_index(hic, d, mask_thresh=.5, hic_matrix=None, correct_sum=False, stat=np.ma.sum):
    chr_bins = hic.chromosome_bins
    n = len(hic.regions())
    if hic_matrix is None:
        log.debug("Fetching matrix")
        hic_matrix = hic[:,:]
    if not hasattr(hic_matrix, "mask"):
        hic_matrix = hic_matrix.masked_matrix
    di = np.empty(n)
    log.debug("Starting processing")
    skipped = 0
    for i, r in enumerate(hic.regions()):
        if (i - chr_bins[r.chromosome][0] < d or
            chr_bins[r.chromosome][1] - i <= d + 1):
            di[i] = np.nan
            continue
        if hic_matrix.mask[i, i]:
            di[i] = np.nan
            continue
        up_slice = (i, slice(i - d, i))
        up_masked = np.sum(hic_matrix.mask[up_slice])
        down_slice = (i, slice(i + 1, i + d + 1))
        down_masked = np.sum(hic_matrix.mask[down_slice])
        if (up_masked > d*mask_thresh or
            down_masked > d*mask_thresh):
            # If too close to the edge of chromosome or
            # if more than half of the entries in this quadrant are masked (unmappable)
            # exclude it from the analysis
            skipped += 1
            di[i] = np.nan
            continue
        up_sum = stat(hic_matrix[up_slice])
        down_sum = stat(hic_matrix[down_slice])
        expected_sum = (up_sum + down_sum)/2
        if correct_sum:
            up_sum /= 1 - (up_masked/d)
            down_sum /= 1 - (down_masked/d)
        di[i] = ((down_sum - up_sum)/abs(down_sum - up_sum))*(((up_sum - expected_sum)**2)/expected_sum + ((down_sum - expected_sum)**2)/expected_sum)
    log.info("Skipped {} regions because >{:.1%} of matrix positions were masked".format(skipped, mask_thresh))
    return di

def observed_expected_ratio(hic, hic_matrix=None, per_chromosome=True, stat=np.ma.mean):
    chr_bins = hic.chromosome_bins
    n = len(hic.regions())
    if hic_matrix is None:
        log.debug("Fetching matrix")
        hic_matrix = hic[:,:]
    hic_matrix = hic_matrix.masked_matrix
    oe = hic_matrix.copy()
    log.debug("starting processing")
    if per_chromosome:
        for c_start, c_end in chr_bins.itervalues():
            # Correcting intrachromoc_startmal contacts by mean contact count at each diagonal
            for i in range(c_end - c_start):
                ind = kth_diag_indices(c_end - c_start, -i)
                oe[c_start:c_end, c_start:c_end][ind] /= stat(hic_matrix[c_start:c_end, c_start:c_end][ind])
            # Correcting interchromoc_startmal contacts by mean of all contact counts between
            # each set of chromoc_startmes
            for other_start, other_end in chr_bins.itervalues():
                # Only correct upper triangle
                if other_start <= c_start:
                    continue
                oe[c_start:c_end, other_start:other_end] /= stat(hic_matrix[c_start:c_end, other_start:other_end])
    else:
        for i in range(n):
            oe[kth_diag_indices(n, -i)] /= stat(hic_matrix.diagonal(i))
    # Copying upper triangle to lower triangle
    oe[np.tril_indices(n)] = oe.T[np.tril_indices(n)]
    return oe

def delta_window(x, window):
    n = len(x)
    delta = np.empty(n)
    for i in range(n):
        if (i < window or n - i <= window - 1):
            delta[i] = np.nan
            continue
        delta[i] = np.mean(x[i + 1:i + window + 1] - x[i]) - np.mean(x[i - window:i] - x[i])
    return delta

def expected(hic, hic_matrix=None, per_chromosome=True, stat=np.ma.mean):
    chr_bins = hic.chromosome_bins
    n = len(hic.regions())
    if hic_matrix is None:
        log.debug("Fetching matrix")
        hic_matrix = hic[:,:]
    hic_matrix = hic_matrix.masked_matrix
    log.debug("starting processing")
    expected = np.zeros(hic_matrix.shape)
    if per_chromosome:
        for so, c_end in chr_bins.itervalues():
            for i in range(eo - so):
                expected[so:eo, so:eo][kth_diag_indices(eo - so, -i)] = stat(hic_matrix[so:eo, so:eo][kth_diag_indices(eo - so, -i)])
            # Correcting interchromosomal contacts by mean of all contact counts between
            # each set of chromosomes
            for si, ei in chr_bins.itervalues():
                # Only correct upper triangle
                if si <= so:
                    continue
                expected[so:eo, si:ei] = stat(hic_matrix[so:eo, si:ei])
    else:
        for i in range(n):
            expected[kth_diag_indices(n, -i)] = stat(hic_matrix.diagonal(-i))
    # Copying upper triangle to lower triangle
    expected[np.tril_indices(n)] = expected.T[np.tril_indices(n)]
    return expected

def create_gtf_from_region_ix(regions, region_ix):
    for s, e in region_ix:
        yield "{}\tkaic\tboundary\t{}\t{}\t.\t.\t.\n".format(
            regions[s].chromosome, regions[s].start, regions[e].end)

class PeakCallerDelta(object):
    def __init__(self, x, window_size=7):
        self.x = x
        self.window_size = window_size
        self._call_peaks()

    def _call_peaks(self):
        self.delta = delta_window(self.x, self.window_size)
        self._peaks = np.nonzero(np.diff(np.signbit(self.delta)))[0]
        log.info("Found {} raw peaks".format(len(self._peaks)))
        self.delta_d1 = savgol_filter(self.delta, window_length=self.window_size, polyorder=2, deriv=1)
        self._delta_peaks = np.nonzero(np.diff(np.signbit(self.delta_d1)))[0]
        self.delta_d2 = savgol_filter(self.delta, window_length=self.window_size, polyorder=2, deriv=2)
        # Figure out which delta zero crossings are minima in x
        self._min_mask = self.delta_d1[self._peaks] > 0
        # Figure out which delta peaks are minima, have d2 > 0
        self._delta_min_mask = self.delta_d2[self._delta_peaks] > 0
        # Find local extrema in delta on left and right side of x extrema
        self._left_value = np.full(self._peaks.shape, np.nan)
        self._right_value = np.full(self._peaks.shape, np.nan)
        self._right_ix = np.searchsorted(self._delta_peaks, self._peaks, side="right")
        for i, p in enumerate(self._peaks):
            try:
                self._left_value[i] = self.delta[self._delta_peaks[self._right_ix[i] - 1]]
            except IndexError:
                pass
            try:
                self._right_value[i] = self.delta[self._delta_peaks[self._right_ix[i]]]
            except IndexError:
                pass
        #left_ix = np.searchsorted(self._delta_peaks, self._peaks, side="right")
        #left_peak = self.delta[self._delta_peaks[left_ix]]
        #right_peak = self.delta[self._delta_peaks[left_ix + 1]]
        # Score
        self._scores = np.abs(self._left_value - self._right_value)
        self._peak_mask = np.full(self._peaks.shape, True, dtype=np.bool_)

    def get_peaks(self):
        return self._peaks[self._peak_mask]

    def get_minima(self):
        return self._peaks[np.logical_and(self._peak_mask, self._min_mask)]

    def get_maxima(self):
        return self._peaks[np.logical_and(self._peak_mask, ~self._min_mask)]

    @property
    def scores(self):
        return self._scores[self._peak_mask]

    def filter(self, delta_score_thresh=None):
        if delta_score_thresh:
            delta_score_pass = self._scores > delta_score_thresh
            self._peak_mask = np.logical_and(self._peak_mask, delta_score_pass)
            log.info("Discarding {}({:.1%}) of total peaks due to delta score threshold ({})".format(np.sum(~delta_score_pass), np.sum(~delta_score_pass)/len(self._peaks), delta_score_thresh))

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
        # Minima have positive 2nd order derivative at minimum
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
        if abs_thresh:
            abs_pass = np.where(self._min_mask, self.x_smooth[self._peaks] < abs_thresh[0], self.x_smooth[self._peaks] > abs_thresh[1])
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
