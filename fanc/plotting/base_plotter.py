from __future__ import division, print_function
from fanc.config import config
from fanc.plotting.helpers import style_ticks_whitegrid, LimitGroup
from matplotlib.ticker import MaxNLocator, Formatter, Locator
from genomic_regions import GenomicRegion
from abc import abstractmethod, abstractproperty, ABCMeta
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
from future.utils import with_metaclass, string_types
import logging
logger = logging.getLogger(__name__)


class GenomeCoordFormatter(Formatter):
    """
    Process axis tick labels to give nice representations
    of genomic coordinates
    """
    def __init__(self, chromosome, minor_div=None, display_scale=True, display_chromosome=True):
        """
        :param chromosome: :class:`~fanc.data.genomic.GenomicRegion` or string
        :param minor_div: Divide each major tick by this many minor ticks.
        :param display_scale: Boolean
                              Display distance scale at bottom right
        """
        if isinstance(chromosome, GenomicRegion):
            self.chromosome = chromosome.chromosome
        else:
            self.chromosome = chromosome
        self.minor_div = minor_div
        self.display_scale = display_scale
        self.display_chromosome = display_chromosome

    def _format_val(self, x, prec_offset=0):
        if x == 0:
            oom_loc = 0
        else:
            oom_loc = int(math.floor(math.log10(abs(x))))
        view_range = self.axis.axes.get_xlim()
        oom_range = int(math.floor(math.log10(abs(view_range[1] - view_range[0]))))

        if oom_loc >= 6:
            return "{:.{prec}f}Mb".format(x / 1000000, prec=max(1, 6 + prec_offset - oom_range))
        elif oom_loc >= 3:
            return "{:.{prec}f}kb".format(x/1000, prec=max(1, 3 + prec_offset - oom_range))
        return "{:.0f}b".format(x)

    def __call__(self, x, pos=None):
        """
        Return label for tick at coordinate x. Relative position of
        ticks can be specified with pos. First tick gets chromosome name.
        """
        s = self._format_val(x, prec_offset=1)
        if self.display_chromosome and (pos == 0 or x == 0):
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
        tick_dist = abs(self.locs[1] - self.locs[0])
        minor_tick_dist = abs(tick_dist/self.minor_div)
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


class MinorGenomeCoordLocator(mpl.ticker.AutoMinorLocator):
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
            majorstep = majorlocs[1] - majorlocs[0]
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
            t0 = majorlocs[0]
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
                # new_init_doc += "\nArguments from {}:".format(b.__name__)
                new_init_doc += b.__init__.__doc__
        if len(new_init_doc) > 0:
            dct["__init__"].__doc__ = new_init_doc
        return super(PlotMeta, cls).__new__(cls, clsname, bases, dct)


class BasePlotter(with_metaclass(PlotMeta, object)):

    def __init__(self, title='', aspect=1., axes_style="ticks", ylabel=None,
                 draw_ticks=True, draw_minor_ticks=False, draw_major_ticks=True,
                 draw_tick_labels=True, draw_tick_legend=False, draw_x_axis=True,
                 draw_chromosome_label=True, padding=None, extra_padding=0,
                 fix_chromosome=False, invert_x=False, ax=None, **kwargs):
        """
        :param str title: Title drawn on top of the figure panel.
        :param float aspect: Aspect ratio of the plot, height/width.
                             So 0.5 means half as high as wide.
                             Default: Plot type specific sensible value
        :param axes_style: Set styling of the axes, can be anything
                           that seaborn supports. See
                           http://seaborn.pydata.org/tutorial/aesthetics.html#styling-figures-with-axes-style-and-set-style
        :param str ylabel: Label for y-axis. Default: None
        :param bool draw_ticks: Draw tickmarks. Default: True
        :param bool draw_major_ticks: Draw major tickmarks. Default: True
        :param bool draw_minor_ticks: Draw minor tickmarks. Default: True
        :param bool draw_tick_labels: Draw tick labels. Default: True
        :param bool draw_x_axis: If False, remove genome x-axis completely.
                                 Default: True
        :param bool draw_tick_legend: Draw legend for the tick distances. Default: True
        :param float padding: Padding in inches to the next plot. Default: None,
                              automatically calculated.
        :param float extra_padding: Add or subtract the specified inches from
                                    the automatically calculated padding from
                                    this panel to the next.
        :param bool fix_chromosome: If True modify chromosome identifiers for this plot,
                                    removing or adding 'chr' as necessary. Default: False
        :param bool invert_x: Invert x-axis. Default=False
                              Caution: This only works reliably when ``independent_x=True``
                              in the GenomicFigure instance, since x-axis of all plots
                              are linked otherwise.
        :param Axes ax: Matplotlib axes instance that the plot will be drawn on.
                        Only necessary if you don't intend to use a 
                        :class:`~fanc.plotting.plotter.GenomicFigure` to draw the
                        plot. Default: None
        """
        if len(kwargs) > 0:
            raise TypeError("Unexpected keyword argument used: {}".format(kwargs))
        self.ax = ax
        self.cax = None
        self.title = title
        self._has_legend = False

        self._draw_ticks = draw_ticks
        self._draw_major_ticks = draw_major_ticks
        self._draw_minor_ticks = draw_minor_ticks
        if not self._draw_ticks:
            self._draw_major_ticks = False
            self._draw_minor_ticks = False
        elif not self._draw_major_ticks or not self._draw_minor_ticks:
            self._draw_ticks = False
        self._draw_chromosome_label = draw_chromosome_label

        self._draw_tick_labels = draw_tick_labels
        self._draw_x_axis = draw_x_axis
        self._draw_tick_legend = draw_tick_legend
        self.aspect = aspect
        self.padding = padding
        # Total padding excl title and adj slider padding, set from GenomicFigure
        self._total_padding = None
        self.extra_padding = extra_padding
        self.axes_style = axes_style
        self.overlays = []
        self.ylabel = ylabel
        self.fix_chromosome = fix_chromosome
        self._dimensions_stale = False
        self.invert_x = invert_x
        self._drawn = False

    def _before_plot(self, region):
        if self._drawn:
            self._clear()
        self._drawn = True
        self.ax.set_title(self.title)
        if self.ylabel and len(self.ylabel) > 0:
            self.ax.set_ylabel(self.ylabel)

    def _after_plot(self, region):
        for o in self.overlays:
            o.plot(self, region)
        if not self._draw_ticks:
            self.remove_genome_ticks(major=not self._draw_major_ticks,
                                     minor=not self._draw_minor_ticks)
        if not self._draw_tick_labels:
            self.remove_genome_labels()
        if not self._draw_x_axis:
            self.remove_genome_axis()
        if not self._draw_tick_legend:
            self.remove_tick_legend()
        if self.invert_x and not self.ax.xaxis_inverted():
            self.ax.invert_xaxis()

    @abstractmethod
    def _plot(self, region):
        raise NotImplementedError("Subclasses need to override _plot function")

    def plot(self, region):
        if isinstance(region, string_types):
            region = GenomicRegion.from_string(region)
        if self.fix_chromosome:
            region = region.fix_chromosome(copy=True)
        if self.ax is None:
            self.ax = plt.gca()
        self._before_plot(region)
        plot_output = self._plot(region)
        self._after_plot(region)

        if plot_output is None:
            return self.fig, self.ax
        return plot_output

    @property
    def fig(self):
        return self.ax.figure

    def add_legend(self, *args, **kwargs):
        if not self._has_legend:
            self.ax.legend(*args, **kwargs)
            self._has_legend = True

    def remove_genome_ticks(self, major=True, minor=True):
        """
        Remove all genome coordinate tickmarks.
        """
        if major:
            plt.setp(self.ax.xaxis.get_majorticklines(), visible=False)
            self.ax.set_xticks([])
        if minor:
            plt.setp(self.ax.xaxis.get_minorticklines(), visible=False)
        if major or minor:
            self._draw_ticks = False

    def remove_genome_labels(self):
        """
        Remove all genome coordinate labels.
        """
        if self.ax:
            plt.setp(self.ax.get_xticklabels(), visible=False)
        self._draw_tick_labels = False

    def remove_genome_axis(self):
        """
        Remove the genome x-axis completely.
        """
        if self.ax:
            self.ax.xaxis.set_visible(False)
            self.ax.spines["bottom"].set_visible(False)
        self._draw_x_axis = False
        self.remove_genome_labels()
        self.remove_genome_ticks()
        self.remove_tick_legend()

    def remove_tick_legend(self):
        """
        Remove the tick mark legend.
        """
        if self.ax:
            self.ax.xaxis.offsetText.set_visible(False)
            self.ax.yaxis.offsetText.set_visible(False)
        self._draw_tick_legend = False

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

    def _clear(self):
        """
        Clear all axes of this plot.

        Called when the plot() is called a second time in order to clear
        the axes for the second plotting.
        """
        if self.ax:
            self.ax.cla()
        if self.cax:
            self.cax.cla()

    def show(self, *args, **kwargs):
        plt.show(*args, **kwargs)

    def save(self, *args, **kwargs):
        self.ax.figure.savefig(*args, **kwargs)


class BasePlotter1D(BasePlotter):

    def __init__(self, n_ticks=5, n_minor_ticks=5, **kwargs):
        """
        :param int n_ticks: Number of major x-axis genome coordinate ticks
        :param int n_minor_ticks: Number of minor ticks per major tick
        """
        super(BasePlotter1D, self).__init__(**kwargs)
        if n_ticks < 2:
            raise ValueError("Need at least two ticks. Set draw_ticks to False to hide all ticks.")
        self.n_tick_bins = n_ticks - 1
        self.n_minor_ticks = n_minor_ticks
        self._mouse_release_handler = None
        self._last_xlim = None
        self._current_chromosome = None

    def _before_plot(self, region):
        super(BasePlotter1D, self)._before_plot(region)
        self.ax.xaxis.set_major_formatter(GenomeCoordFormatter(region,
                                                               minor_div=self.n_minor_ticks,
                                                               display_chromosome=self._draw_chromosome_label))
        self.ax.xaxis.set_major_locator(GenomeCoordLocator(nbins=self.n_tick_bins))
        self.ax.xaxis.set_minor_locator(MinorGenomeCoordLocator(n=self.n_minor_ticks))
        self.ax.set_xlim(region.start, region.end)
        self._current_chromosome = region.chromosome

    def _after_plot(self, region):
        super(BasePlotter1D, self)._after_plot(region)
        self._mouse_release_handler = self.fig.canvas.mpl_connect('button_release_event', self._mouse_release_event)

    def refresh(self, region):
        self._refresh(region)
        # this should take care of any unwanted ylim changes
        # from custom _refresh methods
        self.ax.set_xlim(self._last_xlim)

    @abstractmethod
    def _refresh(self, region):
        raise NotImplementedError("Subclasses need to override _refresh function")

    def _mouse_release_event(self, event):
        xlim = self.ax.get_xlim()

        if xlim != self._last_xlim:
            self._last_xlim = xlim
            x_start, x_end = (max(xlim[0], 0), xlim[1]) if xlim[0] < xlim[1] else (xlim[1], xlim[0])
            x_region = GenomicRegion(chromosome=self._current_chromosome, start=int(x_start), end=int(x_end))
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
                 condensed=False, n_yticks=3, **kwargs):
        """
        :param str style: 'step' Draw values in a step-wise manner for each bin
                          'mid' Draw values connecting mid-points of bins
        :param tupleylim: Set y-axis limits as tuple. Can leave upper or lower
                          limit undetermined by setting None, e.g. (2.5, None).
                          Alternatively, a :class:`~fanc.plotting.helpers.LimitGroup` instance can
                          be passed to synchronize limits across multiple plots.
                          Default: Automatically determined by data limits
        :type lim: tuple(float, float)
        :param str yscale: Scale of y-axis. Is passed to Matplotlib set_yscale,
                           so any valid argument ("linear", "log", etc.) works
                           Default: "linear"
        :param bool condensed: Only show maximum y-axis tick.
                          Default: False
        :param int n_yticks: Number of y-axis ticks. If only the maximum
                             tick should be displayed set condensed to True.
                             Default: 3
        """
        kwargs.setdefault("aspect", .2)
        kwargs.setdefault("axes_style", style_ticks_whitegrid)
        super(ScalarDataPlot, self).__init__(**kwargs)
        self.style = style
        if isinstance(ylim, LimitGroup):
            self.ylim = None
            self.ylim_group = ylim
        else:
            self.ylim = ylim
            self.ylim_group = None
        self.yscale = yscale
        self.condensed = condensed
        if n_yticks < 2:
            raise ValueError("At least 2 ticks needed. Use condensed argument for only one.")
        self.n_yticks = n_yticks
        if style not in self._STYLES:
            raise ValueError("Only the styles {} are supported.".format(list(self._STYLES.keys())))

    def _before_plot(self, region):
        super(ScalarDataPlot, self)._before_plot(region)
        self.ax.set_yscale(self.yscale)
        self.ax.yaxis.set_major_locator(MaxNLocator(nbins=self.n_yticks - 1,
                                                    min_n_ticks=self.n_yticks))

    def _after_plot(self, region):
        super(ScalarDataPlot, self)._after_plot(region)
        if self.ylim:
            self.ax.set_ylim(self.ylim)
        if self.condensed:
            low, high = self.ax.get_ylim()
            self.ax.set_yticks([high])
            self.ax.set_yticklabels([high], va='top', size='large')

    def _get_values_per_step(self, regions, attribute='score'):
        x = []
        y = []
        for region in regions:
            x.append(region.start)
            x.append(region.end)
            v = getattr(region, attribute)
            y.append(v)
            y.append(v)
        return np.array(x), np.array(y)

    def _get_values_per_mid(self, regions, attribute='score'):
        x = []
        y = []
        for region in regions:
            x.append(int(region.center))
            y.append(getattr(region, attribute))
        return np.array(x), np.array(y)

    def values_from_region_iter(self, regions, attribute='score'):
        """
        Get x, y coordinates from a region iterator.

        :param regions: iterable over :class:`~genomic_regions.GenomicRegion` objects
        :param attribute: attribute to extract from each region as score (default: 'score')
        :return: x, y
        """
        return self._STYLES[self.style](self, regions, attribute)

    _STYLES = {_STYLE_STEP: _get_values_per_step,
               _STYLE_MID: _get_values_per_mid}


class BasePlotterMatrix(with_metaclass(PlotMeta, object)):
    """
    Mix-in class to provide methods for mapping colorvalues
    in special areas in the plots etc.
    Does not inherit from BasePlotter, since it's meant to be inherited
    together with another BasePlotter class.
    """

    def __init__(self, colormap=config.colormap_hic, norm="lin", vmin=None, vmax=None,
                 show_colorbar=True, blend_zero=False, replacement_color=None,
                 unmappable_color=".9", illegal_color=None, colorbar_symmetry=None,
                 cax=None, **kwargs):
        """
        :param colormap: Can be the name of a colormap or a Matplotlib colormap instance
        :param norm: Can be "lin", "log" or any Matplotlib Normalization instance
        :param float vmin: Clip interactions below this value
        :param float vmax: Clip interactions above this value
        :param bool show_colorbar: Draw a colorbar. Default: True
        :param bool blend_zero: If True then zero count bins will be drawn using replacement_color
        :param str replacement_color: If None use the lowest color in the colormap, otherwise
                                      use the specified color. Can be any valid matplotlib
                                      color name or specification.
        :param unmappable_color: Draw unmappable bins using this color. Defaults to
                                 light gray (".9")
        :param illegal_color: Draw non-finite (NaN, +inf, -inf) bins using this color. Defaults to
                                 None (no special color).
        :param float colorbar_symmetry: Set to enforce that the colorbar is symmetric around
                                        this value. Default: None
        """
        super(BasePlotterMatrix, self).__init__(**kwargs)
        if isinstance(colormap, string_types):
            colormap = mpl.cm.get_cmap(colormap)
        self.colormap = colormap
        self._vmin = vmin
        self._vmax = vmax
        self.unmappable_color = unmappable_color
        self.blend_zero = blend_zero
        self.illegal_color = illegal_color
        self.show_colorbar = show_colorbar
        self.colorbar_symmetry = colorbar_symmetry
        self.colorbar = None
        self.replacement_color = replacement_color
        self.cax = cax
        self.norm = norm
        self._map_norm = None

    def _before_plot(self, region):
        super(BasePlotterMatrix, self)._before_plot(region)

    def _after_plot(self, region):
        super(BasePlotterMatrix, self)._after_plot(region)
        if self.show_colorbar:
            self.add_colorbar()
        else:
            self.remove_colorbar_ax()

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

        if not hasattr(matrix, 'mask'):
            zero_mask = np.isclose(matrix, 0.) | np.isnan(matrix)
        else:
            zero_mask = np.ma.getmaskarray(matrix)

        if self._vmin is not None:
            vmin = self._vmin
        else:
            if self.norm != 'log':
                vmin = np.nanmin(matrix)
            else:
                vmin = np.nanmin(matrix[matrix > 0])
        if self._vmax is not None:
            vmax = self._vmax
        else:
            vmax = np.nanmax(matrix)

        if self.colorbar_symmetry is not None:
            max_diff = max(abs(self.colorbar_symmetry - vmax), abs(self.colorbar_symmetry - vmin))
            vmin = self.colorbar_symmetry - max_diff
            vmax = self.colorbar_symmetry + max_diff

        if self._map_norm is None:
            self._map_norm = _prepare_normalization(norm=self.norm, vmin=vmin, vmax=vmax)
        color_matrix = self.colormap(self._map_norm(matrix))

        if self.blend_zero:
            color_matrix[zero_mask] = replacement_color
        if self.illegal_color is not None:
            color_matrix[~np.isfinite(matrix)] = mpl.colors.colorConverter.to_rgba(self.illegal_color)
        if self.unmappable_color is not None:
            color_matrix[np.all(zero_mask, axis=1), :] = mpl.colors.colorConverter.to_rgba(self.unmappable_color)
            if matrix.shape[0] == matrix.shape[1]:
                color_matrix[:, np.all(zero_mask, axis=0)] = mpl.colors.colorConverter.to_rgba(self.unmappable_color)
        return color_matrix

    def add_colorbar(self, ax=None, baseline=None, **kwargs):
        """
        Add colorbar to the plot.

        :param ax: Optional axis on which to draw the colorbar. Default: colorbar ax
        :param baseline: symmetric axis around this value. Asymmetric if None.
        """
        kwargs.setdefault('orientation', 'vertical')
        if baseline is None and self.colorbar_symmetry is not None:
            baseline = self.colorbar_symmetry

        if ax is None and self.cax is not None:
            ax = self.cax
        cmap_data = mpl.cm.ScalarMappable(norm=self._map_norm, cmap=self.colormap)

        cmap_data.set_array(np.array([self.vmin, self.vmax]))
        self.colorbar = plt.colorbar(cmap_data, cax=ax, **kwargs)
        self.update_colorbar(vmin=self.vmin, vmax=self.vmax, baseline=baseline)
        return self.colorbar

    def update_colorbar(self, vmin=None, vmax=None, baseline=None):
        if baseline is None and self.colorbar_symmetry is not None:
            baseline = self.colorbar_symmetry

        if baseline is None:
            self.colorbar.mappable.set_clim(vmax=vmax, vmin=vmin)
        else:
            old_vmin, old_vmax = self.colorbar.mappable.get_clim()
            if vmin is None:
                vmin = old_vmin
            if vmax is None:
                vmax = old_vmax

            max_diff = max(abs(baseline - vmax), abs(baseline - vmin))
            self.colorbar.mappable.set_clim(vmin=baseline-max_diff, vmax=baseline+max_diff)
        self.colorbar.draw_all()

    @property
    def vmin(self):
        return self._vmin if self._vmin else self._map_norm.vmin

    @property
    def vmax(self):
        return self._vmax if self._vmax else self._map_norm.vmax

    def _update_norm(self, norm=None, vmin=None, vmax=None):
        if vmin is None:
            vmin = self._map_norm.vmin
        if vmax is None:
            vmax = self._map_norm.vmax
        if norm is None:
            norm = self._map_norm
        self._map_norm = _prepare_normalization(norm, vmin, vmax)


class BasePlotter2D(BasePlotter):

    def __init__(self, n_ticks=3, n_minor_ticks=5, **kwargs):
        """
        :param int n_ticks: Number of major ticks
        :param int n_minor_ticks: Number of minor ticks per major tick
        """
        kwargs.setdefault("aspect", 1.)
        super(BasePlotter2D, self).__init__(**kwargs)
        self.n_tick_bins = n_ticks + 2
        self.n_minor_ticks = n_minor_ticks
        self._mouse_release_handler = None
        self._current_chromosome_x = None
        self._current_chromosome_y = None
        self._last_ylim = None
        self._last_xlim = None

    def _before_plot(self, region):
        x_region, y_region = region
        super(BasePlotter2D, self)._before_plot(x_region)
        self.ax.xaxis.set_major_formatter(GenomeCoordFormatter(x_region, minor_div=self.n_minor_ticks,
                                                               display_chromosome=self._draw_chromosome_label))
        self.ax.xaxis.set_major_locator(GenomeCoordLocator(nbins=self.n_tick_bins))
        self.ax.xaxis.set_minor_locator(MinorGenomeCoordLocator(n=self.n_minor_ticks))
        self.ax.yaxis.set_major_formatter(GenomeCoordFormatter(y_region, minor_div=self.n_minor_ticks,
                                                               display_chromosome=self._draw_chromosome_label))
        self.ax.yaxis.set_major_locator(GenomeCoordLocator(nbins=self.n_tick_bins))
        self.ax.yaxis.set_minor_locator(MinorGenomeCoordLocator(n=self.n_minor_ticks))

        self._current_chromosome_x = x_region.chromosome
        self._current_chromosome_y = y_region.chromosome

    def _after_plot(self, region):
        x_region, y_region = region
        super(BasePlotter2D, self)._after_plot(region)
        self._mouse_release_handler = self.fig.canvas.mpl_connect('button_release_event', self._mouse_release_event)

    def refresh(self, region):
        self._refresh(region)

        # this should take care of any unwanted ylim changes
        # from custom _refresh method
        self.ax.set_xlim(self._last_xlim)
        self.ax.set_ylim(self._last_ylim)

    def remove_genome_ticks(self, major=True, minor=True):
        """
        Remove all genome coordinate tickmarks.
        """
        if major:
            plt.setp(self.ax.xaxis.get_majorticklines(), visible=False)
            plt.setp(self.ax.yaxis.get_majorticklines(), visible=False)
            self.ax.set_yticks([], minor=False)
        if minor:
            plt.setp(self.ax.xaxis.get_minorticklines(), visible=False)
            plt.setp(self.ax.yaxis.get_minorticklines(), visible=False)
            self.ax.set_yticks([], minor=True)
        if major or minor:
            self._draw_ticks = False

    @abstractmethod
    def _refresh(self, region):
        raise NotImplementedError("Subclasses need to override _refresh function")

    def _mouse_release_event(self, event):
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()

        if xlim != self._last_xlim or ylim != self._last_ylim:
            self._last_xlim = xlim
            self._last_ylim = ylim
            x_start, x_end = (xlim[0], xlim[1]) if xlim[0] < xlim[1] else (xlim[1], xlim[0])
            x_region = GenomicRegion(chromosome=self._current_chromosome_x, start=x_start, end=x_end)
            y_start, y_end = (ylim[0], ylim[1]) if ylim[0] < ylim[1] else (ylim[1], ylim[0])
            y_region = GenomicRegion(chromosome=self._current_chromosome_y, start=y_start, end=y_end)
            self.refresh(region=(x_region, y_region))

    def plot(self, regions):
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

        if self.ax is None:
            self.ax = plt.gca()
        self._before_plot((x_region, y_region))
        plot_output = self._plot((x_region, y_region))
        self._after_plot((x_region, y_region))

        if plot_output is None:
            return self.fig, self.ax
        return plot_output


class BaseOverlayPlotter(with_metaclass(PlotMeta, object)):

    def __init__(self, **kwargs):
        pass

    @abstractproperty
    def compatibility(self):
        return None

    @abstractmethod
    def _plot(self, base_plot, region):
        raise NotImplementedError("Subclasses need to override _plot function")

    def plot(self, base_plot, region):
        self._plot(base_plot, region)


class BaseAnnotation(with_metaclass(PlotMeta, object)):

    def __init__(self, fix_chromosome=False):
        """
        :param bool fix_chromosome: If True modify chromosome identifiers for this plot,
                                    removing or adding 'chr' as necessary. Default: False
        """
        self.fix_chromosome = fix_chromosome

    def plot(self, region):
        if isinstance(region, string_types):
            region = GenomicRegion.from_string(region)
        if self.fix_chromosome:
            chromosome = region.chromosome
            if chromosome.startswith('chr'):
                chromosome = chromosome[3:]
            else:
                chromosome = 'chr' + chromosome
            region = GenomicRegion(chromosome=chromosome, start=region.start, end=region.end)
        self._plot(region)

    @abstractmethod
    def _plot(self, region):
        raise NotImplementedError("Subclasses need to override _plot function")
