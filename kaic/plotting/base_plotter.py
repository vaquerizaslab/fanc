from __future__ import division, print_function
from kaic.config import config
from kaic.plotting.helpers import style_ticks_whitegrid, check_kwargs
from matplotlib.ticker import MaxNLocator, Formatter, Locator
from kaic.data.genomic import GenomicRegion
from abc import abstractmethod, abstractproperty, ABCMeta
import numpy as np
import matplotlib as mpl
import seaborn as sns
import math
from future.utils import with_metaclass, string_types
import logging
logger = logging.getLogger(__name__)

plt = sns.plt


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

        if oom_loc >= 6:
            return "{:.{prec}f}Mb".format(x / 1000000, prec=max(0, 6 + prec_offset - oom_range))
        elif oom_loc >= 3:
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


class PlotMeta(ABCMeta):
    """
    Metaclass for all plotting classes. Automatically adds
    all parent classes' __init__ method docstrings to the subclasses'
    __init__ method docstring.
    """
    def __new__(cls, clsname, bases, dct):
        new_init_doc = dct["__init__"].__doc__
        if new_init_doc is None:
            new_init_doc = ""
        for b in bases:
            if b.__name__ == "object":
                continue
            if b.__init__.__doc__ is not None and len(b.__init__.__doc__) > 1:
                new_init_doc += "\nArguments from {}:".format(b.__name__)
                new_init_doc += b.__init__.__doc__
        if len(new_init_doc) > 0:
            dct["__init__"].__doc__ = new_init_doc
        return super(PlotMeta, cls).__new__(cls, clsname, bases, dct)

class BasePlotter(with_metaclass(PlotMeta, object)):

    def __init__(self, title='', aspect=1., axes_style="ticks",
                 draw_ticks=True, draw_tick_labels=True, **kwargs):
        """
        :param title: Title drawn on top of the figure panel
        :param aspect: Aspect ratio of the plot. Can be overriden by setting
                       the height_ratios in class:`~GenomicFigure`
        :param axes_style: Set styling of the axes, can be anything
                           that seaborn supports. See
                           http://seaborn.pydata.org/tutorial/aesthetics.html#styling-figures-with-axes-style-and-set-style
        :param draw_ticks: Draw tickmarks. Default: True
        :param draw_tick_labels: Draw tick labels. Default: True
        """
        check_kwargs(self, kwargs)
        self.ax = None
        self.cax = None
        self.title = title
        self.has_legend = False
        self._has_ticks = draw_ticks
        self._has_tick_labels = draw_tick_labels
        self._aspect = aspect
        self.axes_style = axes_style
        self.overlays = []

    def _before_plot(self, region=None, *args, **kwargs):
        self.ax.set_title(self.title)

    def _after_plot(self, region=None, *args, **kwargs):
        for o in self.overlays:
            o.plot(self, region)
        if not self._has_ticks:
            self.remove_genome_ticks()
        if not self._has_tick_labels:
            self.remove_genome_labels()

    @abstractmethod
    def _plot(self, region=None, *args, **kwargs):
        raise NotImplementedError("Subclasses need to override _plot function")

    def plot(self, region=None, ax=None, *args, **kwargs):
        if ax is None:
            self.ax = plt.gca()
        else:
            self.ax = ax

        if isinstance(region, string_types):
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
        """
        Remove all genome coordinate tickmarks.
        """
        if self.ax:
            plt.setp(self.ax.get_xticklines(), visible=False)
            self.ax.xaxis.set_visible(False)
            self.ax.spines["bottom"].set_visible(False)
        else:
            self._has_ticks = False

    def remove_genome_labels(self):
        """
        Remove all genome coordinate labels.
        """
        if self.ax:
            plt.setp(self.ax.get_xticklabels(), visible=False)
            self.ax.xaxis.offsetText.set_visible(False)
        else:
            self._has_ticklabels = False

    def remove_colorbar_ax(self):
        if self.cax is None:
            return
        try:
            self.fig.delaxes(self.cax)
            self.cax = None
        except KeyError:
            pass

    def add_overlay(self, overlay):
        plot_class = self.__class__.__name__
        if plot_class not in overlay.compatibility:
            raise ValueError("Overlay type {} is not compatible with plotter type {}.".format(
                overlay.__class__.__name__, plot_class))
        self.overlays.append(overlay)


class BasePlotter1D(with_metaclass(PlotMeta, BasePlotter)):

    def __init__(self, n_ticks=3, n_minor_ticks=5, **kwargs):
        """
        :param n_ticks: Number of major ticks
        :param n_minor_ticks: Number of minor ticks per major tick
        """
        BasePlotter.__init__(self, **kwargs)
        self.n_tick_bins = n_ticks + 2
        self.n_minor_ticks = n_minor_ticks
        self._mouse_release_handler = None
        self._last_xlim = None
        self._current_chromosome = None

    def _before_plot(self, region=None, *args, **kwargs):
        BasePlotter._before_plot(self, region=region, *args, **kwargs)
        self.ax.xaxis.set_major_formatter(GenomeCoordFormatter(region))
        self.ax.xaxis.set_major_locator(GenomeCoordLocator(nbins=self.n_tick_bins))
        self.ax.xaxis.set_minor_locator(MinorGenomeCoordLocator(n=self.n_minor_ticks))
        self._current_chromosome = region.chromosome

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
            x_region = GenomicRegion(x_start, x_end, self._current_chromosome)
            self.refresh(region=x_region)


class ScalarDataPlot(BasePlotter1D):
    """
    Base class for plotting scalar values like ChIP-seq signal.
    Provides methods for converting lists of values and regions
    to plotting coordinates.
    """
    _STYLE_STEP = "step"
    _STYLE_MID = "mid"

    def __init__(self, style="step", ylim=None, yscale="linear",
                 condensed=False, n_yticks=2, **kwargs):
        """
        :param style: 'step' Draw values in a step-wise manner for each bin
                      'mid' Draw values connecting mid-points of bins
        :param ylim: Set y-axis limits as tuple. Can leave upper or lower
                     limit undetermined by setting None, e.g. (2.5, None).
                     Default: Automatically determined by data limits
        :param yscale: Scale of y-axis. Is passed to Matplotlib set_yscale,
                       so any valid argument ("linear", "log", etc.) works
                       Default: "linear"
        :param condensed: Only show maximum y-axis tick.
                          Default: False
        :param n_yticks: Number of y-axis ticks. If only the maximum
                         tick should be displayed set condensed to True.
                         Default: 2
        """
        kwargs.setdefault("aspect", .2)
        kwargs.setdefault("axes_style", style_ticks_whitegrid)
        BasePlotter1D.__init__(self, **kwargs)
        self.style = style
        self.ylim = ylim
        self.yscale = yscale
        self.condensed = condensed
        if n_yticks < 2:
            raise ValueError("At least 2 ticks needed. Use condensed argument for only one.")
        self.n_ybins = n_yticks - 1
        if style not in self._STYLES:
            raise ValueError("Only the styles {} are supported.".format(list(self._STYLES.keys())))

    def _before_plot(self, region=None, *args, **kwargs):
        BasePlotter1D._before_plot(self, region=region, *args, **kwargs)
        self.ax.yaxis.set_major_locator(MaxNLocator(self.n_ybins))

    def _after_plot(self, region=None, *args, **kwargs):
        BasePlotter1D._after_plot(self, region=region, *args, **kwargs)
        if self.ylim:
            self.ax.set_ylim(self.ylim)
        self.ax.set_yscale(self.yscale)
        if self.condensed:
            low, high = self.ax.get_ylim()
            self.ax.set_yticks([high])
            self.ax.set_yticklabels([high], va='top', size='large')

    def _get_values_per_step(self, values, region_list):
        x = np.empty(len(region_list)*2)
        y = np.empty(len(region_list)*2)
        for i, r in enumerate(region_list):
            j = i*2
            x[j], x[j + 1] = r.start, r.end
            y[j:j + 2] = values[i]
        return x, y

    def _get_values_per_mid(self, values, region_list):
        x = np.empty(len(values), dtype=np.int_)
        for i, r in enumerate(region_list):
            x[i] = int((r.end + r.start)/2 + 0.5)
        return x, values

    def get_plot_values(self, values, region_list):
        """
        Convert values and regions to final x- and y-
        coordinates for plotting, based on the selected style.

        :param values: List or array of scalar values
        :param region_list: List of class:`~kaic.data.genomic.GenomicRegion`,
                            one for each value
        """
        return self._STYLES[self.style](self, values, region_list)

    _STYLES = {_STYLE_STEP: _get_values_per_step,
               _STYLE_MID: _get_values_per_mid}


class BasePlotterMatrix(with_metaclass(PlotMeta, object)):
    """
    Mix-in class to provide methods for mapping colorvalues
    in special areas in the plots etc.
    Does not inherit from BasePlotter, since it's meant to be inherited
    together with another BasePlotter class.
    """

    def __init__(self, colormap=config.colormap_hic, norm="log", vmin=None, vmax=None,
                 show_colorbar=True, blend_zero=True, replacement_color=None,
                 unmappable_color=".9", illegal_color=None, colorbar_symmetry=None, **kwargs):
        """
        :param colormap: Can be the name of a colormap or a Matplotlib colormap instance
        :param norm: Can be "lin", "log" or any Matplotlib Normalization instance
        :param vmin: Clip interactions below this value
        :param vmax: Clip interactions above this value
        :param show_colorbar: Draw a colorbar
        :param blend_zero: If True then zero count bins will be drawn using replacement_color
        :param replacement_color: If None use the lowest color in the colormap, otherwise
                                  use the specified color. Can be any valid matplotlib
                                  color name or specification.
        :param unmappable_color: Draw unmappable bins using this color. Defaults to
                                 light gray (".9")
        :param illegal_color: Draw non-finite (NaN, +inf, -inf) bins using this color. Defaults to
                                 None (no special color).
        :param colorbar_symmetry: Set to enforce that the colorbar is symemtric around
                                  this value. Default: None
        """
        check_kwargs(self, kwargs)
        if isinstance(colormap, string_types):
            colormap = mpl.cm.get_cmap(colormap)
        self.colormap = colormap
        self._vmin = vmin
        self._vmax = vmax
        self.norm = _prepare_normalization(norm=norm, vmin=vmin, vmax=vmax)
        self.unmappable_color = unmappable_color
        self.blend_zero = blend_zero
        self.illegal_color = illegal_color
        self.show_colorbar = show_colorbar
        self.colorbar_symmetry = colorbar_symmetry
        self.colorbar = None
        self.replacement_color = replacement_color
        self.cax = None
        if isinstance(self.show_colorbar, mpl.axes.Axes):
            self.cax = self.show_colorbar

    def get_color_matrix(self, matrix):
        """
        Convert matrix of scalar values to final
        matrix of RGB values suitable for plotting.

        :param matrix: Data matrix with values to be plotted
        """
        if self.replacement_color is not None:
            cv = mpl.colors.ColorConverter()
            replacement_color = cv.to_rgba(self.replacement_color)
        else:
            replacement_color = self.colormap(0)

        color_matrix = self.colormap(self.norm(matrix))
        if self.blend_zero or self.unmappable_color:
            if not hasattr(matrix, 'mask'):
                zero_mask = np.isclose(matrix, 0.) | np.isnan(matrix)
            else:
                zero_mask = np.ma.getmaskarray(matrix)

            if self.blend_zero:
                color_matrix[zero_mask] = replacement_color
        if self.illegal_color:
            color_matrix[~np.isfinite(matrix)] = mpl.colors.colorConverter.to_rgba(self.illegal_color)
        if self.unmappable_color:
            color_matrix[np.all(zero_mask, axis=1), :] = mpl.colors.colorConverter.to_rgba(self.unmappable_color)
            if matrix.shape[0] == matrix.shape[1]:
                color_matrix[:, np.all(zero_mask, axis=0)] = mpl.colors.colorConverter.to_rgba(self.unmappable_color)
        return color_matrix

    def add_colorbar(self, ax=None, baseline=None):
        """
        Add colorbar to the plot.

        :param ax: Optional axis on which to draw the colorbar. Default: colorbar ax
        :param baseline: symmetric axis around this value. Asymmetric if None.
        """
        if baseline is None and self.colorbar_symmetry is not None:
            baseline = self.colorbar_symmetry

        if ax is None and self.cax is not None:
            ax = self.cax
        cmap_data = mpl.cm.ScalarMappable(norm=self.norm, cmap=self.colormap)

        cmap_data.set_array(np.array([self.vmin, self.vmax]))
        self.colorbar = plt.colorbar(cmap_data, cax=ax, orientation="vertical")
        self.update_colorbar(vmin=self.vmin, vmax=self.vmax, baseline=baseline)

    def update_colorbar(self, vmin=None, vmax=None, baseline=None):
        if baseline is None and self.colorbar_symmetry is not None:
            baseline = self.colorbar_symmetry

        if baseline is None:
            self.colorbar.set_clim(vmax=vmax, vmin=vmin)
        else:
            old_vmin, old_vmax = self.colorbar.get_clim()
            if vmin is None:
                vmin = old_vmin
            if vmax is None:
                vmax = old_vmax

            max_diff = max(abs(baseline - vmax), abs(baseline - vmin))
            self.colorbar.set_clim(vmin=baseline-max_diff, vmax=baseline+max_diff)
        self.colorbar.draw_all()

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


class BasePlotter2D(with_metaclass(PlotMeta, BasePlotter)):

    def __init__(self, n_ticks=3, n_minor_ticks=5, **kwargs):
        """
        :param n_ticks: Number of major ticks
        :param n_minor_ticks: Number of minor ticks per major tick
        """
        kwargs.setdefault("aspect", 1.)
        BasePlotter.__init__(self, **kwargs)
        self.n_tick_bins = n_ticks + 2
        self.n_minor_ticks = n_minor_ticks
        self._mouse_release_handler = None
        self._current_chromosome_x = None
        self._current_chromosome_y = None
        self._last_ylim = None
        self._last_xlim = None

    def _before_plot(self, region=None, *args, **kwargs):
        x_region, y_region = region
        BasePlotter._before_plot(self, region=x_region, *args, **kwargs)
        self.ax.xaxis.set_major_formatter(GenomeCoordFormatter(x_region))
        self.ax.xaxis.set_major_locator(GenomeCoordLocator(nbins=self.n_tick_bins))
        self.ax.xaxis.set_minor_locator(MinorGenomeCoordLocator(n=self.n_minor_ticks))
        self.ax.yaxis.set_major_formatter(GenomeCoordFormatter(y_region))
        self.ax.yaxis.set_major_locator(GenomeCoordLocator(nbins=self.n_tick_bins))
        self.ax.yaxis.set_minor_locator(MinorGenomeCoordLocator(n=self.n_minor_ticks))
        self._current_chromosome_x = x_region.chromosome
        self._current_chromosome_y = y_region.chromosome

    def _after_plot(self, region=None, *args, **kwargs):
        x_region, y_region = region
        BasePlotter._after_plot(self, region=x_region, *args, **kwargs)
        self.ax.set_xlim(x_region.start, x_region.end)
        self.ax.set_ylim(y_region.start, y_region.end)
        self._mouse_release_handler = self.fig.canvas.mpl_connect('button_release_event', self._mouse_release_event)

    def refresh(self, region=None, *args, **kwargs):
        self._refresh(region, *args, **kwargs)

        # this should take care of any unwanted ylim changes
        # from custom _refresh methods
        self.ax.set_xlim(self._last_xlim)
        self.ax.set_ylim(self._last_ylim)

    @abstractmethod
    def _refresh(self, region=None, *args, **kwargs):
        raise NotImplementedError("Subclasses need to override _refresh function")

    def _mouse_release_event(self, event):
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()

        if xlim != self._last_xlim or ylim != self._last_ylim:
            self._last_xlim = xlim
            self._last_ylim = ylim
            x_start, x_end = (xlim[0], xlim[1]) if xlim[0] < xlim[1] else (xlim[1], xlim[0])
            x_region = GenomicRegion(x_start, x_end, self._current_chromosome_x)
            y_start, y_end = (ylim[0], ylim[1]) if ylim[0] < ylim[1] else (ylim[1], ylim[0])
            y_region = GenomicRegion(y_start, y_end, self._current_chromosome_y)
            self.refresh(region=(x_region, y_region))

    def plot(self, regions=None, ax=None, *args, **kwargs):
        if ax is None:
            self.ax = plt.gca()
        else:
            self.ax = ax

        if isinstance(regions, tuple):
            x_region, y_region = regions
        else:
            x_region = regions
            y_region = x_region

        if isinstance(x_region, string_types):
            x_region = GenomicRegion.from_string(x_region)

        if isinstance(y_region, string_types):
            y_region = GenomicRegion.from_string(y_region)

        self._current_chromosome_x = x_region.chromosome
        self._current_chromosome_y = y_region.chromosome

        self._before_plot(region=(x_region, y_region), *args, **kwargs)
        plot_output = self._plot(region=(x_region, y_region), *args, **kwargs)
        self._after_plot(region=(x_region, y_region), *args, **kwargs)

        if plot_output is None:
            return self.fig, self.ax
        return plot_output


class BaseOverlayPlotter(with_metaclass(PlotMeta, object)):

    def __init__(self, **kwargs):
        check_kwargs(self, kwargs)

    @abstractproperty
    def compatibility(self):
        return None

    @abstractmethod
    def _plot(self, base_plot, region):
        raise NotImplementedError("Subclasses need to override _plot function")

    def plot(self, base_plot, region):
        self._plot(base_plot, region)
