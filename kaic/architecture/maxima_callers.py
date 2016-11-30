import numpy as np
from scipy.signal import savgol_filter
from kaic.tools.matrix import delta_window
import logging
logger = logging.getLogger(__name__)


class BaseMaximaCaller(object):
    def __init__(self):
        self._peaks = None
        self._min_mask = None
        self._peak_mask = None
        self._scores = None

    def get_peaks(self):
        return self._peaks[self._peak_mask], self._scores[self._peak_mask]

    def get_minima(self):
        return (self._peaks[np.logical_and(self._min_mask, self._peak_mask)],
                self._scores[np.logical_and(self._min_mask, self._peak_mask)])

    def get_maxima(self):
        return (self._peaks[np.logical_and(~self._min_mask, self._peak_mask)],
                self._scores[np.logical_and(~self._min_mask, self._peak_mask)])


class MaximaCallerMatrix(BaseMaximaCaller):
    def __init__(self, x, seed_delta_ix=10, seed_delta_window=7):
        BaseMaximaCaller.__init__(self)
        self.x = x
        self.seed_delta_ix = seed_delta_ix
        self.seed_delta_window = seed_delta_window
        self._call_peaks()
        self.count = 0
        self.below = 0

    def _call_peaks(self):
        self.delta = delta_window(self.x[:, self.seed_delta_ix], self.seed_delta_window)
        self.delta_d1 = savgol_filter(self.delta, window_length=self.seed_delta_window, polyorder=2, deriv=1)
        # Figure out which delta zero crossings are minima in x
        self._peaks = np.nonzero(np.diff(np.signbit(self.delta)))[0]
        self._min_mask = self.delta_d1[self._peaks] > 0
        self._peak_mask = np.full(self._peaks.shape, True, dtype=np.bool_)

    def filter(self, min_threshold, max_threshold, start=0, end=None, lenience=.1):
        self._scores = np.ma.abs(np.ma.sum(self.x[self._peaks, start:end], axis=1))
        self.below = np.where(self._min_mask[:, np.newaxis], self.x[self._peaks] < min_threshold,
                              self.x[self._peaks] > max_threshold)
        self.count = np.sum(self.below[:, start:end], axis=1)
        pass_mask = self.count > (1 - lenience)*(end - start)
        logger.info("Discarding {}({:.1%}) of total peaks due to thresholds ({}, {})".format(np.sum(~pass_mask),
                                                                                             np.sum(~pass_mask) /
                                                                                             len(self._peaks),
                                                                                             min_threshold,
                                                                                             max_threshold))
        self._peak_mask = np.logical_and(self._peak_mask, pass_mask)


class MaximaCallerDelta(BaseMaximaCaller):
    def __init__(self, x, window_size=7):
        BaseMaximaCaller.__init__(self)
        self.x = x
        self.window_size = window_size
        self._call_peaks()

    def _call_peaks(self):
        self.delta = delta_window(self.x, self.window_size)
        self._peaks = np.nonzero(np.diff(np.signbit(self.delta)))[0]
        logger.info("Found {} raw peaks".format(len(self._peaks)))
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

        # left_ix = np.searchsorted(self._delta_peaks, self._peaks, side="right")
        # left_peak = self.delta[self._delta_peaks[left_ix]]
        # right_peak = self.delta[self._delta_peaks[left_ix + 1]]
        # Score
        self._scores = np.abs(self._left_value - self._right_value)
        self._peak_mask = np.full(self._peaks.shape, True, dtype=np.bool_)

    def filter(self, delta_score_thresh=None):
        if delta_score_thresh:
            delta_score_pass = self._scores > delta_score_thresh
            self._peak_mask = np.logical_and(self._peak_mask, delta_score_pass)
            logger.info("Discarding {}({:.1%}) of total peaks due to delta score threshold ({})".format(
                np.sum(~delta_score_pass), np.sum(~delta_score_pass)/len(self._peaks), delta_score_thresh)
            )
