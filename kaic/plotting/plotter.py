from matplotlib.ticker import Formatter, MaxNLocator, LinearLocator
from matplotlib.widgets import Slider
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

def prepare_normalization(norm="lin", vmin=None, vmax=None):
    if norm == "log":
        return mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
    elif norm == "lin":
        return mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    else:
        raise ValueError("'{}'' not a valid normalization method.".format(norm))

def plot(plots, region):
    n = len(plots)
    fig, axes = plt.subplots(n, sharex=True)
    for p, a in zip(plots, axes):
        p.plot(region, a)
    return fig, axes

class GenomeCoordFormatter(Formatter):
    def __init__(self, chromosome=None, start=None):
        if isinstance(chromosome, GenomicRegion):
            self.chromosome = chromosome.chromosome
            self.start = chromosome.start
        else:
            self.chromosome = chromosome
            self.start = start

    def __call__(self, x, pos=None):
        if self.start is not None and self.chromosome is not None and (pos == 0 or x == 0):
            return "{}:{}".format(self.chromosome, self.start)
        return millify(x)

class BufferedMatrix(object):
    def __init__(self, hic_data):
        self.hic_data = hic_data
        self.buffered_x_region = None
        self.buffered_y_region = None
        self.buffered_matrix = None

    def is_buffered_region(self, x_region, y_region):
        if (self.buffered_y_region is None or not self.buffered_y_region.contains(y_region) or
                self.buffered_x_region is None or not self.buffered_x_region.contains(x_region) or
                self.buffered_matrix is None):
            return False
        return True

    def get_matrix(self, x_region=None, y_region=None):
        if not self.is_buffered_region(x_region, y_region):
            log.info("Buffering matrix")
            if x_region.start is not None and x_region.end is not None:
                x_region_size = x_region.end-x_region.start
                new_x_start = max(1, x_region.start-x_region_size)
                new_x_end = x_region.end + x_region_size
                self.buffered_x_region = GenomicRegion(new_x_start, new_x_end, x_region.chromosome)
            else:
                self.buffered_x_region = GenomicRegion(None, None, x_region.chromosome)

            if y_region.start is not None and y_region.end is not None:
                y_region_size = y_region.end-y_region.start
                new_y_start = max(1, y_region.start-y_region_size)
                new_y_end = y_region.end + y_region_size
                self.buffered_y_region = GenomicRegion(new_y_start, new_y_end, y_region.chromosome)
            else:
                self.buffered_y_region = GenomicRegion(None, None, y_region.chromosome)
            self.buffered_matrix = self.hic_data[self.buffered_x_region, self.buffered_y_region]
        return self.buffered_matrix[y_region, x_region]

    @property
    def buffered_min(self):
        return np.ma.min(self.buffered_matrix) if self.buffered_matrix is not None else None

    @property
    def buffered_max(self):
        return np.ma.max(self.buffered_matrix) if self.buffered_matrix is not None else None

class BasePlotter(object):

    __metaclass__ = ABCMeta

    def __init__(self):
        self._ax = None

    @abstractmethod
    def _plot(self, region=None):
        raise NotImplementedError("Subclasses need to override _plot function")

    @abstractmethod
    def _refresh(self, region=None):
        raise NotImplementedError("Subclasses need to override _refresh function")

    @abstractmethod
    def plot(self, region=None):
        raise NotImplementedError("Subclasses need to override plot function")
    
    @property
    def fig(self):
        return self._ax.figure

    @property
    def ax(self):
        if not self._ax:
            log.debug("Creating new figure object.")
            _, self._ax = plt.subplots()
        return self._ax

    @ax.setter
    def ax(self, value):
        self._ax = value

class BasePlotter1D(BasePlotter):

    __metaclass__ = ABCMeta

    def plot(self, region=None, ax=None):
        if isinstance(region, basestring):
            region = GenomicRegion.from_string(region)
        if ax:
            self.ax = ax
        # set genome tick formatter
        self.ax.xaxis.set_major_formatter(GenomeCoordFormatter(region))

        self._plot(region)
        return self.fig, self.ax

class BasePlotterHic(object):

    __metaclass__ = ABCMeta

    def __init__(self, hic_data, colormap='viridis', norm="log",
                 vmin=None, vmax=None, show_colorbar=True, adjust_range=True):
        self.hic_data = hic_data
        self.hic_buffer = BufferedMatrix(hic_data)
        self.colormap = mpl.cm.get_cmap(colormap)
        self._vmin = vmin
        self._vmax = vmax
        self.norm = prepare_normalization(norm=norm, vmin=vmin, vmax=vmax)
        self.cax = None
        self.colorbar = None
        self.slider = None
        self.show_colorbar = show_colorbar
        self.adjust_range = adjust_range

    def add_colorbar(self):
        cmap_data = mpl.cm.ScalarMappable(norm=self.norm, cmap=self.colormap)
        cmap_data.set_array([self.vmin, self.vmax])
        self.cax, kw = mpl.colorbar.make_axes(self.ax, location="top", shrink=0.4)
        self.colorbar = plt.colorbar(cmap_data, cax=self.cax, **kw)

    def add_adj_slider(self):
        plot_position = self.cax.get_position()
        vmin_axs = plt.axes([plot_position.x0, 0.05, plot_position.width, 0.03], axisbg='#f3f3f3')
        self.vmin_slider = Slider(vmin_axs, 'vmin', self.vmin, self.vmax, valinit=self.vmin,
                                  facecolor='#dddddd', edgecolor='none')
        vmax_axs = plt.axes([plot_position.x0, 0.02, plot_position.width, 0.03], axisbg='#f3f3f3')
        self.vmax_slider = Slider(vmax_axs, 'vmax', self.vmin, self.vmax, valinit=self.vmax,
                                  facecolor='#dddddd', edgecolor='none')
        self.fig.subplots_adjust(top=0.90, bottom=0.15)
        self.vmin_slider.on_changed(self._slider_refresh)
        self.vmax_slider.on_changed(self._slider_refresh)

    def _slider_refresh(self, val):
        new_vmin = self.vmin_slider.val
        new_vmax = self.vmax_slider.val
        self.im.set_clim(vmin=new_vmin, vmax=new_vmax)

    @property
    def vmin(self):
        return self._vmin if self._vmin else self.hic_buffer.buffered_min

    @property
    def vmax(self):
        return self._vmax if self._vmax else self.hic_buffer.buffered_max
 
class BasePlotter2D(BasePlotter):

    __metaclass__ = ABCMeta

    def __init__(self):
        BasePlotter.__init__(self)
        self.cid = None
        self.current_chromosome_x = None
        self.current_chromosome_y = None
        self.last_ylim = None
        self.last_xlim = None

    def mouse_release_refresh(self, _):
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()

        if xlim != self.last_xlim or ylim != self.last_ylim:
            self.last_xlim = xlim
            self.last_ylim = ylim
            x_start, x_end = (xlim[0], xlim[1]) if xlim[0] < xlim[1] else (xlim[1], xlim[0])
            x_region = GenomicRegion(x_start, x_end, self.current_chromosome_x)

            y_start, y_end = (ylim[0], ylim[1]) if ylim[0] < ylim[1] else (ylim[1], ylim[0])
            y_region = GenomicRegion(y_start, y_end, self.current_chromosome_y)

            self._refresh(x_region, y_region)

            # this should take care of any unwanted ylim changes
            # from custom _refresh methods
            self.ax.set_ylim(self.last_ylim)
            self.ax.set_xlim(self.last_xlim)

    def plot(self, x_region=None, y_region=None, ax=None):
        if isinstance(x_region, basestring):
            x_region = GenomicRegion.from_string(x_region)
        if ax:
            self.ax = ax

        self.current_chromosome_x = x_region.chromosome

        if y_region is None:
            y_region = x_region

        if isinstance(y_region, basestring):
            y_region = GenomicRegion.from_string(y_region)

        self.current_chromosome_y = y_region.chromosome

        # set base-pair formatters
        self.ax.xaxis.set_major_formatter(GenomeCoordFormatter(x_region))
        self.ax.yaxis.set_major_formatter(GenomeCoordFormatter(y_region))
        # set release event callback
        self.cid = self.fig.canvas.mpl_connect('button_release_event', self.mouse_release_refresh)

        self._plot(x_region, y_region)
        return self.fig, self.ax

class HicPlot2D(BasePlotter2D, BasePlotterHic):
    def __init__(self, hic_data, colormap='viridis', norm="log",
                 vmin=None, vmax=None, show_colorbar=True,
                 adjust_range=True):
        BasePlotter2D.__init__(self)
        BasePlotterHic.__init__(self, hic_data=hic_data, colormap=colormap,
                                norm=norm, vmin=vmin, vmax=vmax, show_colorbar=show_colorbar,
                                adjust_range=adjust_range)

    def _plot(self, x_region=None, y_region=None):
        m = self.hic_buffer.get_matrix(x_region=x_region, y_region=y_region)
        self.im = self.ax.imshow(m, interpolation='nearest', cmap=self.colormap, norm=self.norm,
                                 extent=[m.col_regions[0].start, m.col_regions[-1].end,
                                         m.row_regions[-1].end, m.row_regions[0].start])
        self.last_ylim = self.ax.get_ylim()
        self.last_xlim = self.ax.get_xlim()

        if self.show_colorbar:
            self.add_colorbar()
            if self.adjust_range:
                self.add_adj_slider()

    def _refresh(self, x_region=None, y_region=None):
        print "refreshing"
        m = self.hic_buffer.get_matrix(x_region=x_region, y_region=y_region)

        self.im.set_data(m)
        self.im.set_extent([m.col_regions[0].start, m.col_regions[-1].end,
                            m.row_regions[-1].end, m.row_regions[0].start])


class HicSideBySidePlot2D(object):
    def __init__(self, hic1, hic2, colormap='viridis', norm="log",
                 vmin=None, vmax=None):
        self.hic_plotter1 = HicPlot2D(hic1, colormap=colormap, norm=norm, vmin=vmin, vmax=vmax)
        self.hic_plotter2 = HicPlot2D(hic2, colormap=colormap, norm=norm, vmin=vmin, vmax=vmax)

    def plot(self, region):
        fig = plt.figure()
        ax1 = plt.subplot(121)
        ax2 = plt.subplot(122, sharex=ax1, sharey=ax1)

        self.hic_plotter1.plot(x_region=region, y_region=region, ax=ax1)
        self.hic_plotter2.plot(x_region=region, y_region=region, ax=ax2)

        return fig, ax1, ax2


class HicComparisonPlot2D(HicPlot2D):
    def __init__(self, hic_top, hic_bottom, colormap='viridis', norm='log',
                 vmin=None, vmax=None, scale_matrices=True):
        super(HicComparisonPlot2D, self).__init__(hic_top, colormap=colormap, norm=norm, vmin=vmin, vmax=vmax)
        self.hic_top = hic_top
        self.hic_bottom = hic_bottom
        self.scaling_factor = 1
        if scale_matrices:
            self.scaling_factor = hic_top.scaling_factor(hic_bottom)

    def _get_matrix(self, x_region, y_region):
        print x_region, y_region
        return self.hic_top.get_combined_matrix(self.hic_bottom, key=(y_region, x_region),
                                                scaling_factor=self.scaling_factor)


class HicPlot(BasePlotter1D, BasePlotterHic):
    def __init__(self, hic_data, colormap='viridis', max_dist=None, norm="log",
                 vmin=None, vmax=None, show_colorbar=True, adjust_range=False):
        BasePlotter1D.__init__(self)
        BasePlotterHic.__init__(self, hic_data, colormap=colormap, vmin=vmin, vmax=vmax,
                                show_colorbar=show_colorbar, adjust_range=adjust_range)
        self.max_dist = max_dist

    def _plot(self, region=None):
        log.debug("Generating matrix from hic object")
        if region is None:
            raise ValueError("Cannot plot triangle plot for whole genome.")
        hm = self.hic_buffer.get_matrix(region, region)
        hm[np.tril_indices(hm.shape[0])] = np.nan
        # Remove part of matrix further away than max_dist
        if self.max_dist:
            for i, r in enumerate(hm.row_regions):
                if r.start - region.start > self.max_dist:
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
        log.debug("Plotting matrix")
        # create plot
        self.ax.pcolormesh(X_, Y_, hm_masked, cmap=self.colormap, norm=self.norm)
        # set limits and aspect ratio
        self.ax.set_aspect(aspect="equal")
        self.ax.set_xlim(hm.row_regions[0].start - 1, hm.row_regions[-1].end)
        self.ax.set_ylim(0, self.max_dist if self.max_dist else 0.5*(region.end-region.start))
        # remove y ticks
        self.ax.set_yticks([])
        # Hide the left, right and top spines
        sns.despine(left=True)
        # hide background patch
        self.ax.patch.set_visible(False)
        # Only show ticks on the left and bottom spines
        self.ax.xaxis.set_ticks_position('bottom')

    def _refresh(self, region=None):
        pass
