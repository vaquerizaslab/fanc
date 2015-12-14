from matplotlib.ticker import Formatter, MaxNLocator
from kaic.data.genomic import GenomicRegion
from abc import abstractmethod, ABCMeta
import numpy as np
import math
import matplotlib as mpl
import logging
import seaborn as sns
plt = sns.plt
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
    millnames = ["", "k", "M", "B", "T"]
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


class BasePlotter1D(object):

    __metaclass__ = ABCMeta

    def __init__(self):
        self.fig = None
        self.ax = None

    @abstractmethod
    def _plot(self, region=None):
        raise NotImplementedError("Subclasses need to override _plot function")

    @abstractmethod
    def _refresh(self, region=None):
        raise NotImplementedError("Subclasses need to override _plot function")
    
    def plot(self, region=None, ax=None):
        if isinstance(region, basestring):
            region = GenomicRegion.from_string(region)

        if ax is not None:
            self.ax = ax
            self.fig = ax.figure
        else:
            self.fig, self.ax = plt.subplots()

        self._plot(region)
        return self.fig, self.ax


class BasePlotter2D(object):

    __metaclass__ = ABCMeta

    def __init__(self):
        self.fig = None
        self.ax = None
        self.cid = None
        self.current_chromosome_x = None
        self.current_chromosome_y = None

    @abstractmethod
    def _plot(self, x_region=None, y_region=None):
        raise NotImplementedError("Subclasses need to override _plot function")

    @abstractmethod
    def _refresh(self, x_region=None, y_region=None):
        raise NotImplementedError("Subclasses need to override _refresh function")

    def mouse_release_refresh(self, _):
        xlim = self.ax.get_xlim()
        x_start, x_end = (xlim[0], xlim[1]) if xlim[0] < xlim[1] else (xlim[1], xlim[0])
        x_region = GenomicRegion(x_start, x_end, self.current_chromosome_x)
        ylim = self.ax.get_ylim()
        y_start, y_end = (ylim[0], ylim[1]) if ylim[0] < ylim[1] else (ylim[1], ylim[0])
        y_region = GenomicRegion(y_start, y_end, self.current_chromosome_y)

        self._refresh(x_region, y_region)

    def plot(self, x_region=None, y_region=None, ax=None, interactive=True):
        if isinstance(x_region, basestring):
            x_region = GenomicRegion.from_string(x_region)
            self.current_chromosome_x = x_region.chromosome

        if isinstance(y_region, basestring):
            y_region = GenomicRegion.from_string(y_region)
            self.current_chromosome_y = y_region.chromosome

        if ax is not None:
            self.ax = ax
            self.fig = ax.figure
        else:
            self.fig, self.ax = plt.subplots()

        if interactive:
            plt.ion()
            self.cid = self.fig.canvas.mpl_connect('button_release_event', self.mouse_release_refresh)

        self._plot(x_region, y_region)
        return self.fig, self.ax


class HicPlot2D(BasePlotter2D):
    def __init__(self, hic):
        super(HicPlot2D, self).__init__()
        self.hic = hic

    def _plot(self, x_region=None, y_region=None):
        m = self.hic[y_region, x_region]
        self.im = self.ax.imshow(m, interpolation='none',
                                 extent=[m.col_regions[0].start, m.col_regions[-1].end,
                                         m.row_regions[-1].end, m.row_regions[0].start])

    def _refresh(self, x_region=None, y_region=None):
        m = self.hic[x_region, y_region]

        self.im.set_data(m)
        self.im.set_extent([m.col_regions[0].start, m.col_regions[-1].end,
                            m.row_regions[-1].end, m.row_regions[0].start])


class HicPlot(BasePlotter1D):
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

    def _plot(self, region=None):
        log.debug("Generating matrix from hic object")
        if region is None:
            hm = self.hic_data[:]
        else:
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
        #import ipdb
        #ipdb.set_trace()
