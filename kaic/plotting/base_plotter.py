from __future__ import division, print_function
from matplotlib.ticker import MaxNLocator, Formatter, Locator
from kaic.data.genomic import GenomicRegion
from matplotlib.widgets import Slider
from abc import abstractmethod, ABCMeta
import numpy as np
import matplotlib as mpl
import logging
import seaborn as sns
import math

plt = sns.plt
log = logging.getLogger(__name__)
log.setLevel(10)


def append_axes(parent, side, thickness, padding, length=None, shrink=1., **kwargs):
    """
    Add an axes on any side of parent ax without resizing the parent ax.

    :param parent: Parent axes
    :param side: Side on which axes is appended ("top", "bottom", "left" or "right")
    :param thickness: Thickness of newly created axes, in inches. Measured in
                      direction from the parent axes to the side where axes is created
    :param padding: Padding between parent and new axes
    :param length: Length of new axes perpendicular thickness. By default same length
                   as parent axes
    :param shrink: Set length to a certain fraction of parent axes length. No effect
                   if length is set explicitely.
    :param kwargs: Additional keyword args passed to figure.add_axes method
    :return: Axes instance
    """
    figsize = parent.figure.get_size_inches()
    bbox = parent.get_position()
    if side in ("top", "bottom"):
        thickness = thickness/figsize[1]
        padding = padding/figsize[1]
        length = length/figsize[0] if length is not None else shrink*bbox.width
        if side == "top":
            hor = bbox.x0 + (bbox.width - length)/2
            vert = bbox.y0 + bbox.height + padding
            width = length
            height = thickness
        else:
            length = length if length is not None else shrink*bbox.width
            hor = bbox.x0 + (bbox.width - length)/2
            vert = bbox.y0 - padding - thickness
            width = length
            height = thickness
    elif side in ("right", "left"):
        thickness = thickness/figsize[0]
        padding = padding/figsize[0]
        length = length/figsize[1] if length is not None else shrink*bbox.height
        if side == "right":
            hor = bbox.x0 + bbox.width + padding
            vert = bbox.y0 + (bbox.height - length)/2
            width = thickness
            height = length
        else:
            hor = bbox.x0 - padding - thickness
            vert = bbox.y0 + (bbox.height - length)/2
            width = thickness
            height = length
    else:
        raise ValueError("Illegal parameter side '{}'".format(side))
    return parent.figure.add_axes([hor, vert, width, height], **kwargs)


def force_aspect(ax, aspect):
    """
    Force aspect ratio of a matplotlib axes

    :param ax: axes whose aspect ratio will be forced
    :param aspect: Aspect ratio (width/height)
    """
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    data_aspect = (xlim[1] - xlim[0])/(ylim[1] - ylim[0])
    aspect = abs(data_aspect/aspect)
    ax.set_aspect(aspect)


class SymmetricNorm(mpl.colors.Normalize):
    """
    Normalizes data for plotting on a divergent colormap.
    Automatically chooses vmin and vmax so that the colormap is
    centered at zero.
    """
    def __init__(self, vmin=None, vmax=None, percentile=None):
        """
        :param vmin: Choose vmin manually
        :param vmax: Choose vmax manually
        :param percentile: Instead of taking the minmum or maximum to
                           automatically determine vmin/vmax, take the
                           percentile.
                           eg. with 5, vmin is 5%-ile, vmax 95%-ile
        """
        mpl.colors.Normalize.__init__(self, vmin=vmin, vmax=vmax, clip=False)
        self.percentile = percentile
        self.vmin = vmin
        self.vmax = vmax

    def _get_min(self, A):
        if self.percentile:
            return np.nanpercentile(A, 100 - self.percentile)
        else:
            return np.ma.min(A[~np.isnan(A)])

    def _get_max(self, A):
        if self.percentile:
            return np.nanpercentile(A, self.percentile)
        else:
            return np.ma.max(A[~np.isnan(A)])

    def autoscale(self, A):
        vmin = self._get_min(A)
        vmax = self._get_max(A)
        abs_max = max(abs(vmin), abs(vmax))
        self.vmin = -1.*abs_max
        self.vmax = abs_max

    def autoscale_None(self, A):
        vmin = self.vmin if self.vmin else self._get_min(A)
        vmax = self.vmax if self.vmax else self._get_max(A)
        abs_max = max(abs(vmin), abs(vmax))
        self.vmin = -1.*abs_max
        self.vmax = abs_max


def millify(n, precision=1):
    """
    Take input float and return human readable string.
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
    millidx = max(0, min(len(millnames) - 1, int(math.floor(math.log10(abs(n))/3))))

    return "{:.{prec}f}{}".format(n/10**(3*millidx), millnames[millidx], prec=precision)


class GenomeCoordFormatter(Formatter):
    """
    Process axis tick labels to give nice representations
    of genomic coordinates
    """
    def __init__(self, chromosome, display_scale=True):
        """
        :param chromosome: :class:`~kaic.data.genomic.GenomicRegion` or string
        :param display_scale: Boolean
                              Display distance scale at bottom right
        """
        if isinstance(chromosome, GenomicRegion):
            self.chromosome = chromosome.chromosome
        else:
            self.chromosome = chromosome
        self.display_scale = display_scale

    def _format_val(self, x, prec_offset=0):
        if x == 0:
            oom_loc = 0
        else:
            oom_loc = int(math.floor(math.log10(abs(x))))
        view_range = self.axis.axes.get_xlim()
        oom_range = int(math.floor(math.log10(abs(view_range[1] - view_range[0]))))
        if oom_loc >= 3:
            return "{:.{prec}f}kb".format(x/1000, prec=max(0, 3 + prec_offset - oom_range))
        return "{:.0f}b".format(x)

    def __call__(self, x, pos=None):
        """
        Return label for tick at coordinate x. Relative position of
        ticks can be specified with pos. First tick gets chromosome name.
        """
        s = self._format_val(x, prec_offset=1)
        if pos == 0 or x == 0:
            return "{}:{}".format(self.chromosome, s)
        return s

    def get_offset(self):
        """
        Return information about the distances between
        tick bars and the size of the view window.
        Is called by matplotlib and displayed in lower right corner
        of plots.
        """
        if not self.display_scale or self.locs is None or len(self.locs) == 0:
            return ""
        view_range = self.axis.axes.get_xlim()
        view_dist = abs(view_range[1] - view_range[0])
        tick_dist = self.locs[2] - self.locs[1]
        minor_tick_dist = tick_dist/5
        minor_tick_dist_str = self._format_val(minor_tick_dist, prec_offset=2)
        tick_dist_str = self._format_val(tick_dist, prec_offset=1)
        view_dist_str = self._format_val(view_dist)
        return "{}|{}|{}".format(minor_tick_dist_str, tick_dist_str, view_dist_str)


class GenomeCoordLocator(MaxNLocator):
    """
    Choose locations of genomic coordinate ticks on the plot axis.
    Behaves like default Matplotlib locator, except that it always
    places a tick at the start and the end of the window.
    """
    def __call__(self):
        vmin, vmax = self.axis.get_view_interval()
        ticks = self.tick_values(vmin, vmax)
        # Make sure that first and last tick are the start
        # and the end of the genomic range plotted. If next
        # ticks are too close, remove them.
        ticks[0] = vmin
        ticks[-1] = vmax
        if ticks[1] - vmin < (vmax - vmin)/(self._nbins*3):
            ticks = np.delete(ticks, 1)
        if vmax - ticks[-2] < (vmax - vmin)/(self._nbins*3):
            ticks = np.delete(ticks, -2)
        return self.raise_if_exceeds(np.array(ticks))


class MinorGenomeCoordLocator(Locator):
    """
    Choose locations of minor tick marks between major
    tick labels. Modification of the Matplotlib AutoMinorLocator,
    except that it uses the distance between 2nd and 3rd major
    mark as reference, instead of 2nd and 3rd.
    """
    def __init__(self, n):
        self.ndivs = n

    def __call__(self):
        majorlocs = self.axis.get_majorticklocs()
        try:
            majorstep = majorlocs[2] - majorlocs[1]
        except IndexError:
            # Need at least two major ticks to find minor tick locations
            # TODO: Figure out a way to still be able to display minor
            # ticks without two major ticks visible. For now, just display
            # no ticks at all.
            majorstep = 0
        if self.ndivs is None:
            if majorstep == 0:
                # TODO: Need a better way to figure out ndivs
                ndivs = 1
            else:
                x = int(np.round(10 ** (np.log10(majorstep) % 1)))
                if x in [1, 5, 10]:
                    ndivs = 5
                else:
                    ndivs = 4
        else:
            ndivs = self.ndivs
        minorstep = majorstep / ndivs
        vmin, vmax = self.axis.get_view_interval()
        if vmin > vmax:
            vmin, vmax = vmax, vmin
        if len(majorlocs) > 0:
            t0 = majorlocs[1]
            tmin = ((vmin - t0) // minorstep + 1) * minorstep
            tmax = ((vmax - t0) // minorstep + 1) * minorstep
            locs = np.arange(tmin, tmax, minorstep) + t0
            cond = np.abs((locs - t0) % majorstep) > minorstep / 10.0
            locs = locs.compress(cond)
        else:
            locs = []
        return self.raise_if_exceeds(np.array(locs))


def _prepare_normalization(norm="lin", vmin=None, vmax=None):
    if isinstance(norm, mpl.colors.Normalize):
        norm.vmin = vmin
        norm.vmax = vmax
        return norm
    if norm == "log":
        return mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
    elif norm == "lin":
        return mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    else:
        raise ValueError("'{}'' not a valid normalization method.".format(norm))


class BasePlotter(object):

    __metaclass__ = ABCMeta

    def __init__(self, title='', aspect=1., axes_style="ticks"):
        self.ax = None
        self.cax = None
        self.title = title
        self.has_legend = False
        self._aspect = aspect
        self.axes_style = axes_style

    def _before_plot(self, region=None, *args, **kwargs):
        self.ax.set_title(self.title)

    def _after_plot(self, region=None, *args, **kwargs):
        pass

    @abstractmethod
    def _plot(self, region=None, *args, **kwargs):
        raise NotImplementedError("Subclasses need to override _plot function")

    def plot(self, region=None, ax=None, *args, **kwargs):
        if ax is None:
            self.ax = plt.gca()
        else:
            self.ax = ax

        if isinstance(region, basestring):
            region = GenomicRegion.from_string(region)

        self._before_plot(region=region, *args, **kwargs)
        plot_output = self._plot(region=region, *args, **kwargs)
        self._after_plot(region=region, *args, **kwargs)

        if plot_output is None:
            return self.fig, self.ax
        return plot_output

    @property
    def fig(self):
        return self.ax.figure

    def add_legend(self, *args, **kwargs):
        if not self.has_legend:
            self.ax.legend(*args, **kwargs)

    def get_default_aspect(self):
        return self._aspect

    def remove_genome_ticks(self):
        plt.setp(self.ax.get_xticklabels(), visible=False)
        self.ax.xaxis.offsetText.set_visible(False)


class BasePlotter1D(BasePlotter):

    __metaclass__ = ABCMeta

    def __init__(self, title='', aspect=1., axes_style="ticks"):
        BasePlotter.__init__(self, title=title, aspect=aspect,
                             axes_style=axes_style)
        self._mouse_release_handler = None
        self._last_xlim = None
        self.current_chromosome = None

    def _before_plot(self, region=None, *args, **kwargs):
        BasePlotter._before_plot(self, region=region, *args, **kwargs)
        self.ax.xaxis.set_major_formatter(GenomeCoordFormatter(region))
        self.ax.xaxis.set_major_locator(GenomeCoordLocator(nbins=5))
        self.ax.xaxis.set_minor_locator(MinorGenomeCoordLocator(n=5))
        self.current_chromosome = region.chromosome

    def _after_plot(self, region=None, *args, **kwargs):
        BasePlotter._after_plot(self, region=region, *args, **kwargs)
        self.ax.set_xlim(region.start, region.end)
        self._mouse_release_handler = self.fig.canvas.mpl_connect('button_release_event', self._mouse_release_event)

    def refresh(self, region=None, *args, **kwargs):
        self._refresh(region, *args, **kwargs)

        # this should take care of any unwanted ylim changes
        # from custom _refresh methods
        self.ax.set_xlim(self._last_xlim)

    @abstractmethod
    def _refresh(self, region=None, *args, **kwargs):
        raise NotImplementedError("Subclasses need to override _refresh function")

    def _mouse_release_event(self, event):
        xlim = self.ax.get_xlim()

        if xlim != self._last_xlim:
            self._last_xlim = xlim
            x_start, x_end = (xlim[0], xlim[1]) if xlim[0] < xlim[1] else (xlim[1], xlim[0])
            x_region = GenomicRegion(x_start, x_end, self.current_chromosome)
            self.refresh(region=x_region)


class BasePlotterMatrix(BasePlotter):
    """
    Mix-in class to provide methods for mapping colorvalues
    in special areas in the plots etc.
    """

    __metaclass__ = ABCMeta

    def __init__(self, colormap='viridis', norm="log", vmin=None, vmax=None,
                 show_colorbar=True, blend_zero=True, title='',
                 unmappable_color=".9", illegal_color=None):
        BasePlotter.__init__(self, title=title)

        if isinstance(colormap, basestring):
            colormap = mpl.cm.get_cmap(colormap)

        self.colormap = colormap
        self._vmin = vmin
        self._vmax = vmax
        self.norm = _prepare_normalization(norm=norm, vmin=vmin, vmax=vmax)
        self.unmappable_color = unmappable_color
        self.blend_zero = blend_zero
        self.illegal_color = illegal_color
        self.show_colorbar = show_colorbar
        self.colorbar = None

    def get_color_matrix(self, matrix):
        """
        Convert matrix of scalar values to final
        matrix of RGB values suitable for plotting.
        """
        color_matrix = self.colormap(self.norm(matrix))
        if self.blend_zero or self.unmappable_color:
            zero_mask = np.isclose(matrix, 0.)
        if self.blend_zero:
            color_matrix[zero_mask] = self.colormap(0)
        if self.illegal_color:
            color_matrix[~np.isfinite(matrix)] = mpl.colors.colorConverter.to_rgba(self.illegal_color)
        if self.unmappable_color:
            color_matrix[np.all(zero_mask, axis=0), :] = mpl.colors.colorConverter.to_rgba(self.unmappable_color)
            color_matrix[:, np.all(zero_mask, axis=1)] = mpl.colors.colorConverter.to_rgba(self.unmappable_color)
        return color_matrix

    def add_colorbar(self, ax=None):
        """
        Add colorbar to the plot.

        :param ax: Optional axis on which to draw the colorbar
        """
        cmap_data = mpl.cm.ScalarMappable(norm=self.norm, cmap=self.colormap)
        cmap_data.set_array(np.array([self.vmin, self.vmax]))
        self.colorbar = plt.colorbar(cmap_data, cax=None, orientation="vertical")

    def remove_colorbar_ax(self):
        try:
            self.fig.delaxes(self.colorbar.ax)
        except KeyError:
            pass

    @property
    def vmin(self):
        return self._vmin if self._vmin else self.norm.vmin

    @property
    def vmax(self):
        return self._vmax if self._vmax else self.norm.vmax

    def _update_norm(self, norm=None, vmin=None, vmax=None):
        if vmin is None:
            vmin = self.norm.vmin
        if vmax is None:
            vmax = self.norm.vmax
        if norm is None:
            norm = self.norm
        self.norm = _prepare_normalization(norm, vmin, vmax)

#
# class IntensitySlider(object):
#     def __init__(self, vmin, vmax, init_value=0):
#         self.vmin = vmin
#         self.vmax = vmax
#         self.init_value = init_value
#
#     def show(self, ax):

class BasePlotter2D(BasePlotter):

    __metaclass__ = ABCMeta

    def __init__(self, title, aspect=1., axes_style="ticks"):
        BasePlotter.__init__(self, title=title, aspect=aspect,
                             axes_style=axes_style)
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