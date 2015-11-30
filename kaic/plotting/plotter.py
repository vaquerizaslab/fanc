from matplotlib.ticker import Formatter, MaxNLocator
from kaic.data.genomic import GenomicRegion, HicMatrix
import numpy as np
import math
import matplotlib as mpl
import seaborn as sns
plt = sns.plt
import logging
log = logging.getLogger(__name__)
log.setLevel(10)

def millify(n, precision=1):
    """Take input float and return human readable string.
    E.g.:
    millify(1000f0) -> "10k"
    millify(2300000) -> "2M"

    Parameters
    ----------
    n : int, float
        Number to be converted
    precision : int
        Number of decimals displayed in output string

    Returns
    -------
    str : Human readable string representation of n
    """
    millnames = ["","k","M","B","T"]
    if n == 0:
        return 0
    n = float(n)
    millidx = max(0, min(len(millnames) - 1,
                      int(math.floor(math.log10(abs(n))/3))))
    return "{:.{prec}f}{}".format(n/10**(3*millidx), millnames[millidx], prec=precision)

class GenomeCoordFormatter(Formatter):
    def __init__(self, chromosome, start, end):
        self.chromosome = chromosome
        self.start = start
        self.end = end

    def __call__(self, x, pos=None):
        if pos == 0 or x == 0:
            return "{}:{}".format(self.chromosome, self.start)
        return millify(x)

class GenomeCoordLocator(MaxNLocator):
    def __init__(self, chromosome, start, end, **kwargs):
        self.chromsome = chromosome
        self.start = start
        self.end = end
        super(GenomeCoordLocator, self).__init__(**kwargs)

    def __call__(self):
        vmin, vmax = self.axis.get_view_interval()
        return self.tick_values(max(self.start, vmin), min(self.end, vmax))

class BasePlotter(object):
    def __init__(self):
        self._fig = None
        self._ax = None

    def add_fig_ax(self, fig, ax):
        self._fig = fig
        self._ax = ax

    @property
    def fig(self):
        if self._fig is None:
            self._make_fig_ax()
        return self._fig
    
    @property
    def ax(self):
        if self._ax is None:
             self._make_fig_ax()
        return self._ax

    def _make_fig_ax(self):
        self._fig, self._ax = plt.subplots()
        return self._fig, self._ax

    def _plot(self, region):
        raise NotImplementedError("Subclasses need to override _plot function")
    
    def plot(self, region):
        if isinstance(region, basestring):
            region = GenomicRegion.from_string(region)
        self._plot(region)
        return self.fig, self.ax

class HicPlot(BasePlotter):
    def __init__(self, hic_data, colormap='viridis', max_height=None, norm="log",
                 vmin=None, vmax=None):
        super(HicPlot, self).__init__()
        self.hic_data = hic_data
        self.colormap = colormap
        self.max_height = max_height
        self.vmin = vmin
        self.vmax = vmax
        if norm == "log":
            self.norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
        elif norm == "lin":
            self.norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        else:
            raise ValueError("'{}'' not a valid normalization method.".format(norm))

    def _plot(self, region):
        log.debug("Generating matrix from hic object")
        hm = self.hic_data[region, region]
        hm[np.tril_indices(hm.shape[0])] = np.nan
        # Remove part of matrix further away than max_height
        if self.max_height:
            for i, r in enumerate(hm.row_regions):
                if r.start - region.start > self.max_height:
                    hm[np.triu_indices(hm.shape[0], k=i)] = np.nan
                    break
        hm_masked = np.ma.MaskedArray(hm, mask=np.isnan(hm))
        log.debug("Rotating matrix")
        # prepare an array of the corner coordinates of the Hic-matrix
        # Distances have to be scaled by sqrt(2), because the diagonals of the bins
        # are sqrt(2)*len(bin_size)
        sqrt2 = math.sqrt(2)
        bin_coords = np.r_[[(x.start - 1) for x in hm.row_regions], (hm.row_regions[-1].end)]/sqrt2
        X, Y = np.meshgrid(bin_coords, bin_coords)
        # rotatate coordinate matrix 45 degrees
        sin45 = math.sin(math.radians(45))
        X_, Y_ = X*sin45 + Y*sin45, X*sin45 - Y*sin45
        # shift x coords to correct start coordinate and center the diagonal directly on the 
        # x-axis
        X_ -= np.min(X_) - (hm.row_regions[0].start - 1)
        Y_ -= .5*np.min(Y_) + .5*np.max(Y_)
        with sns.axes_style("ticks"):
            # normalize colors
            cmap = mpl.cm.get_cmap(self.colormap)
            log.debug("Plotting matrix")
            # create plot
            sns.plt.pcolormesh(X_, Y_, hm_masked, axes=self.ax, cmap=cmap, norm=self.norm)
            # set limits and aspect ratio
            self.ax.set_aspect(aspect="equal")
            self.ax.set_xlim(hm.row_regions[0].start - 1, hm.row_regions[-1].end)
            self.ax.set_ylim(0, self.max_height if self.max_height else 0.5*(region.end-region.start))
            log.debug("Setting custom x tick formatter")
            # set genome tick formatter
            self.ax.xaxis.set_major_formatter(GenomeCoordFormatter(region.chromosome, region.start, region.end))
            # remove y ticks
            self.ax.set_yticks([])
            # Hide the left, right and top spines
            sns.despine(left=True)
            # hide background patch
            self.ax.patch.set_visible(False)
            # Only show ticks on the left and bottom spines
            self.ax.xaxis.set_ticks_position('bottom')
            log.debug("Setting tight layout")
            # make figure margins accommodate labels
            sns.plt.tight_layout()
        cmap_data = mpl.cm.ScalarMappable(norm=self.norm, cmap=cmap)
        cmap_data.set_array([self.vmin if self.vmin else np.ma.min(hm_masked), self.vmax if self.vmax else np.ma.max(hm_masked)])
        cax, kw = mpl.colorbar.make_axes(self.ax, location="top", shrink=0.4)
        plt.colorbar(cmap_data, cax=cax, **kw)
        import ipdb
        ipdb.set_trace()
