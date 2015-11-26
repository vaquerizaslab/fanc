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
    """Take input float and return string which is turned human readable.
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
        if x == 0:
            return self.chromosome
        return millify(x + self.start)

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
        self.plotter = None
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
        hm = self.hic_data[region.to_string(), region.to_string()]
        n = hm.shape[0]
        log.debug("Taking upper triangle.")
        # mask areas we don't want to plot
        # lower triangle
        mask_lower = np.tril_indices(n, k=-1)
        hm[mask_lower] = np.nan
        if self.max_height is not None:
            # upper right corner
            mask_upper = np.triu_indices(n, k=self.max_height)
            hm[mask_upper] = np.nan
        triangle = np.ma.masked_array(hm, np.isnan(hm))
        log.debug("Rotating matrix")
        # prepare an array of the corner coordinates of the Hic-matrix
        sqrt2 = math.sqrt(2)
        sin45 = math.sin(math.radians(45))
        bin_coords = np.r_[[sqrt2*(x.start - 1) for x in hm.row_regions], sqrt2*(hm.row_regions[-1].end -1)]
        X, Y = np.meshgrid(bin_coords, bin_coords)
        # rotation coordinate matrix 45 degrees
        t = np.array([[sin45, sin45], [-sin45, sin45]])
        # "rotate" A
        corner_coords = np.dot(np.stack((X, Y), axis=-1), t)
        # and shift x coords to correct start
        corner_coords[:,:,0] += (np.abs(np.min(corner_coords[:,:,0])) + hm.row_regions[0].start)
        import ipdb
        ipdb.set_trace()

        # flip triangle (because that pcolormesh works opposite than imshow)
        flip_triangle = np.flipud(triangle)

        with sns.axes_style("ticks"):
            # normalize colors
            cmap = mpl.cm.get_cmap(self.colormap)

            log.debug("Plotting matrix")
            # create plot
            caxes = sns.plt.pcolormesh(corner_coords[:,:,0], corner_coords[:,:,1], flip_triangle, axes=self.ax, cmap=cmap, norm=self.norm)

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

