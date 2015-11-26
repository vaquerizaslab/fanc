from matplotlib.ticker import Formatter, MaxNLocator
from kaic.data.genomic import GenomicRegion, HicMatrix
import numpy as np
import math
import matpotlib as mpl
import seaborn as sns
plt = sns.plt

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
        self.chromsome = chromosome
        self.start = start
        self.end = end

    def __call__(self, x, pos=None):
        if x == 0:
            return chromosome
        return millify(x + start)

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
        if self.fig is None:
            self.fig, self.ax = self._make_fig_ax()
        return self._fig
    
    @property
    def ax(self):
        if self.ax is None:
            self.fig, self.ax = self._make_fig_ax()
        return self._ax

    def _make_fig_ax(self):
        return plt.subplots()

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
        hm = hic_data[region.to_string(), region.to_string()]
        n = hm.shape[0]

        # mask areas we don't want to plot
        # lower triangle
        mask_lower = np.tril_indices(n, k=-1)
        hm[mask_lower] = np.nan
        if max_height is not None:
            # upper right corner
            mask_upper = np.triu_indices(n, k=max_height)
            hm[mask_upper] = np.nan
        triangle = np.ma.masked_array(hm, np.isnan(hm))

        # prepare an array of tuples that will be used to rotate triangle
        A = np.array([(y, x) for x in range(n, -1, -1) for y in range(n + 1)])
        # rotation matrix 45 degrees
        t = np.array([[0.707, 0.707], [-0.707, 0.707]])
        # "rotate" A
        A = np.dot(A, t)
        # transform A into x and y values
        X = A[:, 1].reshape((n + 1, n + 1))
        Y = A[:, 0].reshape((n + 1, n + 1))

        # flip triangle (because that pcolormesh works opposite than imshow)
        flip_triangle = np.flipud(triangle)

        with sns.axes_style("ticks"):
            # normalize colors
            cmap = mpl.cm.get_cmap(self.colormap)

            # create plot
            caxes = sns.plt.pcolormesh(X, Y, flip_triangle, axes=self.ax, cmap=cmap, norm=self.norm)

            # re-calculate and reset axis limits
            max_x = max(A[:, 1])
            max_y = 0
            for i in xrange(flip_triangle.shape[0]):
                for j in xrange(flip_triangle.shape[1]):
                    if flip_triangle[i, j] is not np.ma.masked:
                        max_y = max(max_y, Y[i, j]+2)
            self.ax.set_ylim((-1, max_y))
            self.ax.set_xlim((0, max_x))

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

            # make figure margins accommodate labels
            sns.plt.tight_layout()

