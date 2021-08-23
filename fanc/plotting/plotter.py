from __future__ import division, print_function
from fanc.config import config
import matplotlib as mpl
from matplotlib.ticker import NullLocator, MaxNLocator
from fanc import load
from genomic_regions import GenomicRegion, GenomicDataFrame, merge_overlapping_regions
from fanc.plotting.base_plotter import BasePlotter1D, ScalarDataPlot, BaseOverlayPlotter, \
                                       BasePlotter, BaseAnnotation
from fanc.plotting.hic_plotter import BasePlotterMatrix
from fanc.plotting.helpers import append_axes, style_ticks_whitegrid, get_region_field, \
                                  region_to_pbt_interval, absolute_wspace_hspace, \
                                  box_coords_abs_to_rel, figure_line, figure_rectangle, \
                                  parse_bedtool_input, get_region_based_object, \
                                  load_score_data
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import types
import numpy as np
import seaborn as sns
import pybedtools as pbt
try:
    from itertools import izip as zip
except ImportError:
    pass
import itertools
import re
import warnings
from collections import defaultdict
from future.utils import string_types
import fanc
from fanc.tools.general import human_format
import logging
logger = logging.getLogger(__name__)


def hide_axis(ax):
    """
    Hide the axis bar, ticks and tick labels on all sides.

    :param ax: matplotlib Axes instance
    """
    sns.despine(ax=ax, top=True, left=True, bottom=True, right=True)
    try:
        plt.setp(ax.get_xticklabels(), visible=False)
    except IndexError:
        pass
    try:
        plt.setp(ax.get_yticklabels(), visible=False)
    except IndexError:
        pass
    try:
        plt.setp(ax.get_xticklines(), visible=False)
    except IndexError:
        pass
    try:
        plt.setp(ax.get_yticklines(), visible=False)
    except IndexError:
        pass
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.xaxis.offsetText.set_visible(False)
    ax.yaxis.offsetText.set_visible(False)


class GenomicFigure(object):
    """
    GenomicFigure composed of one or more plots.
    All plots are arranged in a single column, their genomic coordinates aligned.
    """

    def __init__(self, plots, width=4., ticks_last=False,
                 invert_x=False, cax_padding=.3, cax_width=.3, fig_padding=(.5, .5, 1., 1.),
                 independent_x=False):
        """
        :param list plots: List of plot instances each will form a separate panel in the figure.
                      Should inherit from :class:`~fanc.plotting.baseplotter.BasePlotter` or
                      :class:`~fanc.plotting.baseplotter.BaseAnnotation`
        :param float width: Width of the plots in inches. Height is automatically determined
                      from the specified aspect ratios of the Plots.
                      Default: 5.
        :param bool ticks_last: Only draw genomic coordinate tick labels on last (bottom) plot
        :param bool invert_x: Invert x-axis for on all plots. Default: False
        :param float cax_padding: Distance between plots and the colorbar in inches. Default: 0.3
        :param float cax_width: Width of colorbar in inches. Default: 0.5
        :param fig_padding: Distance between the edges of the plots and the figure borders
                            in inches (bottom, top, left, right). Default: (.5, .5, .5, .5)
        :type fig_padding: tuple(float, float, float, float)
        :param bool independent_x: When plotting, supply separate coordinates for each plot
                              in the figure. Default: False
        """
        self.plots = []
        self.annotations = []
        for p in plots:
            if isinstance(p, BasePlotter):
                self.plots.append(p)
            elif isinstance(p, BaseAnnotation):
                self.annotations.append(p)
            else:
                raise ValueError("Incompatible plot {}.".format(p))
        for p in self.annotations:
            p._verify(self)
        self.n = len(self.plots)
        self.ticks_last = ticks_last
        self._width = width
        self._cax_padding = cax_padding
        self._cax_width = cax_width
        self._fig_padding = fig_padding
        self.invert_x = invert_x
        self._independent_x = independent_x
        self._figure_setup()

    def _calc_figure_setup(self):
        aspects = [p.aspect if p.aspect is not None else 1. for p in self.plots]
        pad_b, pad_t, pad_l, pad_r = self._fig_padding
        total_width = pad_l + pad_r + self._width + self._cax_width + self._cax_padding
        plot_heights = [a*self._width for a in aspects]
        plot_pads = []
        for i, p in enumerate(self.plots):
            if self.ticks_last and i < len(self.axes) - 1:
                p.remove_genome_labels()
                p.remove_tick_legend()
            if p.padding is not None:
                pad = p.padding
            else:
                pad = config.pad_empty_axis + p.extra_padding
                if p._draw_tick_labels:
                    pad += config.pad_with_label
                if p._draw_ticks:
                    pad += config.pad_with_ticks
                if p._draw_tick_legend:
                    pad += config.pad_with_tick_legend
                p._total_padding = pad
                # Add a bit of space if adjustment slider is present
                pad_adj_slider = 0.
                if getattr(p, "adjust_range", False):
                    pad_adj_slider = config.adjustment_slider_height + config.pad_empty_axis
                # Add a bit of space if next plot has title
                pad_title = 0.
                if i < self.n - 1 and len(self.plots[i + 1].title) > 0:
                    pad_title = config.pad_next_title
                pad += pad_adj_slider + pad_title
            plot_pads.append(pad)
        total_height = pad_t + pad_b + sum(plot_heights) + sum(plot_pads)
        figsize = (total_width, total_height)
        cax_l = pad_l + self._width + self._cax_padding
        cur_top = pad_t
        ax_specs = []
        for i in range(self.n):
            specs = {}
            specs["ax"] = list(box_coords_abs_to_rel(cur_top, pad_l, self._width, plot_heights[i], figsize))
            specs["cax"] = list(box_coords_abs_to_rel(cur_top, cax_l, self._cax_width, plot_heights[i], figsize))
            ax_specs.append(specs)
            cur_top += plot_heights[i] + plot_pads[i]
        return ax_specs, figsize

    def _figure_setup(self):
        ax_specs, figsize = self._calc_figure_setup()
        fig = plt.figure(figsize=figsize)
        sharey_ax = None
        for i in range(self.n):
            plot = self.plots[i]
            with sns.axes_style("ticks" if plot.axes_style is None else
                                plot.axes_style):
                ax = fig.add_axes(ax_specs[i]["ax"], sharex=self.axes[0] if i > 0 and not self._independent_x else None,
                                  sharey=sharey_ax if getattr(plot, "sharey", False) else None)
            if not sharey_ax and getattr(plot, "sharey", False):
                sharey_ax = ax
            cax = fig.add_axes(ax_specs[i]["cax"])
            plot.ax = ax
            plot.cax = cax
            if self.invert_x:
                plot.invert_x = True

    def _update_figure_setup(self):
        ax_specs, figsize = self._calc_figure_setup()
        self.fig.set_size_inches(figsize)
        for i in range(self.n):
            self.plots[i].ax.set_position(ax_specs[i]["ax"])
            if self.plots[i].cax is not None:
                self.plots[i].cax.set_position(ax_specs[i]["cax"])

    @property
    def fig(self):
        return self.axes[0].figure

    def plot(self, region):
        """
        Make a plot of the specified region.

        :param region: A string describing a region "2L:10000000-12000000" or
                       a :class:`~fanc.data.genomic.GenomicRegion`.
                       If ``independent_x`` was set, a list of regions
                       equal to the length of plots + the length of
                       annotations (if present) must be supplied.
                       The order of the regions here must be the same
                       as the order in which the plots were supplied to
                       the GenomicFigure constructor.
        :type region: string or ~fanc.data.genomic.GenomicRegion
        :return: A matplotlib Figure instance and a list of figure axes
        """
        if self._independent_x:
            plot_regions = region
            if len(plot_regions) != len(self.plots) + len(self.annotations):
                raise ValueError("{} regions supplied. Figure has {} plots and "
                    "{} annotations.".format(len(plot_regions), len(self.plots), len(self.annotations)))
        else:
            plot_regions = [region]*(len(self.plots) + len(self.annotations))
        # First iteration plot each panel and record ylims for ylim_groups
        for r, p in zip(plot_regions, self.plots):
            p.plot(r)
            if getattr(p, "ylim_group", None) is not None:
                p.ylim_group.add_limit(p.ax.get_ylim())
        # Second iteration apply ylim_group limits
        for p in self.plots:
            if getattr(p, "ylim_group", None) is not None:
                p.ax.set_ylim(p.ylim_group.get_limit())
                p.ax.yaxis.reset_ticks()
        # Third iteration check if any axes dimensions need updating
        # and reset ylim_groups
        dimensions_stale = False
        for p in self.plots:
            if p._dimensions_stale:
                dimensions_stale = True
                p._dimensions_stale = False
            if getattr(p, "ylim_group", None) is not None:
                p.ylim_group.reset_limit()
        if dimensions_stale:
            self._update_figure_setup()
        # Plot annotations
        for p in self.annotations:
            p._verify(self)
        for r, p in zip(plot_regions[len(self.plots):], self.annotations):
            p.plot(r)
        return self.fig, self.axes

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        plt.close(self.fig)

    @property
    def axes(self):
        return [p.ax for p in self.plots]

    @property
    def caxes(self):
        return [p.cax for p in self.plots]


class HighlightAnnotation(BaseAnnotation):
    """
    Vertical lines or rectangles (shaded regions) which can be
    positioned at specific genomic coordinates from a BED file
    and can span multiple panels of the GenomicFigure.

    Useful for highlighting specific regions across multiple
    panels of figures.
    """

    def __init__(self, bed, plot1=None, plot2=None, plot_kwargs=None,
                 **kwargs):
        """
        :param bed: Anything pybedtools can parse. Path to BED-file
                    GTF-file, or a list of tuples [(chr, start, end), ...]
                    If features are 1bp long, lines are drawn. If they
                    are > 1bp rectangles are drawn. Their appearance
                    can be controlled using the plot_kwargs.
        :type bed: string or pybedtool.BedTool
        :param plot1: First plot where line should start. Default: First
        :param plot2: Second plot where line should end. Default: Last
        :param plot_kwargs: Dictionary of properties which are passed
                            to matplotlib.lines.Line2D or
                            matplotlib.patches.Rectangle constructor
        """
        super(HighlightAnnotation, self).__init__(**kwargs)
        self.plot_kwargs = {
            "linewidth": 1.,
            "color": "black",
            "linestyle": "solid",
            "alpha": .5,
        }
        if plot_kwargs is not None:
            self.plot_kwargs.update(plot_kwargs)
        self.bedtool = parse_bedtool_input(bed)
        self.plot1 = self._plot1 = plot1
        self.plot2 = self._plot2 = plot2
        self.patches = []
        self.lines = []

    def _plot(self, region):
        x_trans = self._plot1.ax.transData
        y_trans1 = self._plot1.ax.transAxes + self._plot1.ax.figure.transFigure.inverted()
        y_trans2 = self._plot2.ax.transAxes + self._plot2.ax.figure.transFigure.inverted()
        blended_trans = mpl.transforms.blended_transform_factory(x_trans, self._plot1.ax.figure.transFigure)
        interval = region_to_pbt_interval(region)
        hits = self.bedtool.all_hits(interval)
        for r in hits:
            if len(r) > 1:
                self._draw_rectangle(r, x_trans, y_trans1, y_trans2, blended_trans)
            else:
                self._draw_line(r, x_trans, y_trans1, y_trans2, blended_trans)

    def _draw_rectangle(self, r, x_trans, y_trans1, y_trans2, plot_trans):
        s, e = r.start, r.end
        y1_t = y_trans1.transform((0, 1))[1]
        y2_t = y_trans2.transform((0, 0))[1]
        y1_t, y2_t = sorted([y1_t, y2_t])
        patch = figure_rectangle(self._plot1.ax.figure, xy=(s, y1_t),
                                 width=e - s, height=y2_t - y1_t,
                                 transform=plot_trans, **self.plot_kwargs)
        patch.set_transform(plot_trans)
        self.patches.append(patch)

    def _draw_line(self, r, x_trans, trans1, trans2, plot_trans):
        s = r.start
        y1_t = trans1.transform((0, 1))[1]
        y2_t = trans2.transform((0, 0))[1]
        l = figure_line(self._plot1.ax.figure, xdata=[s, s],
                        ydata=[y1_t, y2_t], transform=plot_trans,
                        **self.plot_kwargs)
        l.set_transform(plot_trans)
        self.lines.append(l)

    def _verify(self, gfig):
        if self.plot1 is None:
            self._plot1 = gfig.plots[0]
        if self.plot2 is None:
            self._plot2 = gfig.plots[-1]
        if not all([self._plot1 in gfig.plots, self._plot2 in gfig.plots]):
            raise ValueError("At least one plot in the HighlightAnnotation is"
                             "not part of the GenomicFigure")
        # Make sure plot1 comes first in plot list
        if gfig.plots.index(self._plot2) < gfig.plots.index(self._plot1):
            self._plot1, self._plot2 = self._plot2, self._plot1
        return True

    def _refresh(self, region):
        for a in itertools.chain(self.patches, self.lines):
            a.remove()
        self.lines = []
        self.patches = []
        self._plot(region)


class GenomicRegionsPlot(ScalarDataPlot):
    """
    Plot scalar values from one or more :class:`~GenomicRegions` objects
    """

    def __init__(self, regions, attributes=None, names=None,
                 colors=None, linestyles=None, legend=True, **kwargs):
        """
        :param regions: :class:`~GenomicRegions`
        :param attributes: Only draw attributes from the track objects
                           which match this description.
                           Should be a list of names. Supports wildcard matching
                           and regex.
        :param names: Supply list of names for each track.
        :param legend: Draw legend. Default: True
        """
        kwargs.setdefault("aspect", .2)
        super(GenomicRegionsPlot, self).__init__(**kwargs)

        self.regions = regions
        self.attributes = attributes
        self.lines = []
        self.names = names
        if colors is None:
            self.colors = itertools.cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
        else:
            self.colors = itertools.cycle(colors)

        if linestyles is None:
            self.linestyles = itertools.cycle(['-'])
        else:
            self.linestyles = itertools.cycle(linestyles)
        self.legend = legend

    def _plot(self, region):
        line_counter = 0
        for i, name in enumerate(self.regions.data_field_names):
            if not self.attributes or any(re.match("^" + a.replace("*", ".*") + "$", name) for a in self.attributes):
                regions = []
                values = []
                for r in self.regions.subset(region):
                    regions.append(r)
                    values.append(getattr(r, name))
                x, y = self.get_plot_values(values, regions)
                if self.names is not None:
                    if not self.attributes:
                        label = self.names[line_counter]
                    else:
                        label = self.names[self.attributes.index(name)]
                elif self.regions.y_values is not None:
                    label = "{}".format(self.regions.y_values[i])
                else:
                    label = "{}".format(name)
                l = self.ax.plot(x, y, label=label,
                                 color=next(self.colors),
                                 linestyle=next(self.linestyles))
                self.lines.append(l[0])
                line_counter += 1

        if self.legend:
            self.add_legend()
        self.remove_colorbar_ax()

    def _refresh(self, region):
        line_counter = 0
        for i, name in enumerate(self.regions.data_field_names):
            if not self.attributes or any(re.match("^" + a.replace("*", ".*") + "$", name) for a in self.attributes):
                regions = []
                values = []
                for r in self.regions.subset(region):
                    regions.append(r)
                    values.append(getattr(r, name))
                x, y = self.get_plot_values(values, regions)
                self.lines[line_counter].set_xdata(x)
                self.lines[line_counter].set_ydata(y)
                line_counter += 1


class RegionsValuesPlot(ScalarDataPlot):
    """
    Plot scalar values from one or more :class:`~GenomicTrack` objects
    """
    def __init__(self, regions, values, symmetry=None, **kwargs):
        """
        :param regions: :class:`~GenomicRegions`
        :param style: 'step' Draw values in a step-wise manner for each bin
              'mid' Draw values connecting mid-points of bins
        :param attributes: Only draw attributes from the track objects
                           which match this description.
                           Should be a list of names. Supports wildcard matching
                           and regex.
        :param title: Used as title for plot
        :param aspect: Default aspect ratio of the plot. Can be overriden by setting
               the height_ratios in :class:`~fanc.plotting.GenomicFigure`
        """
        kwargs.setdefault("axes_style", style_ticks_whitegrid)
        kwargs.setdefault("aspect", .2)
        super(RegionsValuesPlot, self).__init__(**kwargs)

        self.regions = regions
        self.legend = True
        if not isinstance(values, dict):
            self.values = {'data': values}
            self.legend = False
        else:
            self.values = values
        self.lines = []
        self.symmetry = symmetry

    def _plot_values(self, region):
        for label, region_values in self.values.items():
            regions = []
            values = []
            for i, r in enumerate(self.regions):
                if r.chromosome != region.chromosome:
                    continue
                if r.end < region.start or r.start > region.end:
                    continue

                regions.append(r)
                values.append(region_values[i])

            x, y = self.get_plot_values(values, regions)
            yield label, x, y

    def _plot(self, region):
        for label, x, y in self._plot_values(region):
            l = self.ax.plot(x, y, label=label)
            self.lines.append(l[0])

        if self.symmetry is not None:
            ylim = self.ax.get_ylim()
            d = max(abs(np.array(ylim) - self.symmetry))
            new_ylim = (self.symmetry - d, self.symmetry + d)
            self.ax.set_ylim(new_ylim)

        if self.legend:
            self.add_legend()

        self.remove_colorbar_ax()

    def _refresh(self, region):
        for i, (label, x, y) in enumerate(self._plot_values(region)):
            self.lines[i].set_xdata(x)
            self.lines[i].set_ydata(y)


class GenomicVectorArrayPlot(BasePlotterMatrix, BasePlotter1D):
    """
    Plot matrix from a :class:`~fanc.architecture.hic_architecture.MultiVectorArchitecturalRegionFeature` object.
    """

    def __init__(self, array, parameters=None, y_coords=None, y_scale='linear',
                 plot_kwargs=None, genomic_format=False, **kwargs):
        """
        :param array: :class:`~fanc.architecture.hic_architecture.MultiVectorArchitecturalRegionFeature`
        :param keys: keys for which vectors to use for array. None indicates all vectors will be used.
        :param y_coords: Matrices in the :class:`~GenomicTrack` object are
                         unitless. Can provide the coordinates for the
                         y-direction here. Matrix has shape (X, Y) must
                         have shape Y or Y + 1
        :param y_scale: Set scale of the y-axis, is passed to Matplotlib set_yscale, so any
                        valid argument ("linear", "log", etc.) works
        :param plot_kwargs: Keyword-arguments passed on to pcolormesh
        """
        kwargs.setdefault("aspect", .3)
        kwargs.setdefault("norm", 'lin')

        super(GenomicVectorArrayPlot, self).__init__(**kwargs)
        self.array = array
        if parameters is None:
            self.parameters = array._parameters
        else:
            self.parameters = parameters
        if plot_kwargs is None:
            plot_kwargs = {}
        self.plot_kwargs = plot_kwargs
        if y_coords is None:
            if self.array._parameters is not None:
                self.y_coords = self.array._parameters
            else:
                self.y_coords = None
        self.hm = None
        self.y_scale = y_scale
        self._genomic_format = genomic_format

    def _plot(self, region):
        x, y, self.hm = self._mesh_data(region=region)
        self.collection = self.ax.pcolormesh(x, y, np.ma.masked_invalid(self.hm.T),
                                             rasterized=True, cmap=self.colormap,
                                             norm=self._map_norm, **self.plot_kwargs)

        self.collection._A = None
        self._update_mesh_colors()
        self.ax.set_yscale(self.y_scale)
        self.ax.set_ylim(self.parameters[0], self.parameters[-1])

        if self._genomic_format:
            self.ax.set_yticklabels([human_format(tick) for tick in self.ax.get_yticks()])

        def drag_pan(self, button, key, x, y):
            mpl.axes.Axes.drag_pan(self, button, 'x', x, y)  # pretend key=='x'

        self.ax.drag_pan = types.MethodType(drag_pan, self.ax)

    def _mesh_data(self, region):
        hm = self.array.score_matrix(region=region, parameters=self.parameters)
        bin_coords = np.r_[[(x.start - 1) for x in hm.row_regions], hm.row_regions[-1].end]
        x, y = np.meshgrid(bin_coords, self.parameters)
        return x, y, hm

    def _update_mesh_colors(self):
        # pcolormesh doesn't support plotting RGB arrays directly like imshow, have to workaround
        # See https://github.com/matplotlib/matplotlib/issues/4277
        # http://stackoverflow.com/questions/29232439/plotting-an-irregularly-spaced-rgb-image-in-python/29232668?noredirect=1#comment46710586_29232668
        color_matrix = self.get_color_matrix(np.ma.masked_invalid(self.hm))
        color_tuple = color_matrix.transpose((1, 0, 2)).reshape(
            (color_matrix.shape[0] * color_matrix.shape[1], color_matrix.shape[2]))
        self.collection.set_facecolor(color_tuple)

    def _refresh(self, region):
        x, y, self.hm = self._mesh_data(region)

        self.collection._coordinates[:, :, 0] = x
        # update matrix data
        self.collection.set_array(self.hm.T.ravel())
        self._update_mesh_colors()


class MirrorMatrixPlot(BasePlotter1D):
    """
    Stack two plots on top of each other, bottom plot inverted.
    Especially suited to stacking two Hic plots (triangles) on top
    of each other.
    """
    def __init__(self, top_plot, bottom_plot, gap=0, cax_gap=.05, **kwargs):
        """
        :param top_plot: Plot instance on top
        :param bottom_plot: Plot instance on bottom
        :param gap: Gap between plots in inches
        :param cax_gap: Gap between colorbars in inches
        """
        kwargs.setdefault("aspect", 1.)
        super(MirrorMatrixPlot, self).__init__(**kwargs)
        self.top_plot = top_plot
        self.bottom_plot = bottom_plot
        self.parent_ax = None
        self.parent_cax = None
        self.top_ax = None
        self.bottom_ax = None
        self.gap = gap
        self.cax_gap = cax_gap
        self.cax = None

    def _add_split_ax(self, ax, gap, sharex=False):
        bbox = ax.get_position()
        figsize = ax.figure.get_size_inches()
        gap = gap/figsize[1]
        top_ax = ax.figure.add_axes([bbox.x0, bbox.y0 + gap/2 + bbox.height/2,
                                     bbox.width, bbox.height/2 - gap/2], sharex=ax if sharex else None)
        bottom_ax = ax.figure.add_axes([bbox.x0, bbox.y0,
                                        bbox.width, bbox.height/2 - gap/2], sharex=ax if sharex else None)
        return top_ax, bottom_ax

    def _plot(self, region):
        if self.cax is None:
            self.cax = append_axes(self.ax, 'right', 0.3, 0.05)
        # Check if ax has already been split
        if self.parent_ax is not self.ax:
            self.parent_ax = self.ax
            self.top_ax, self.bottom_ax = self._add_split_ax(self.ax, self.gap, sharex=True)
            sns.despine(ax=self.ax, top=True, left=True, bottom=False, right=True)
            self.ax.yaxis.set_major_locator(NullLocator())
            self.ax.yaxis.set_minor_locator(NullLocator())
        if self.parent_cax is not self.cax:
            self.parent_cax = self.cax
            self.top_plot.cax, self.bottom_plot.cax = self._add_split_ax(self.cax, self.cax_gap)
            self.cax.set_visible(False)
        self.top_plot.ax = self.top_ax
        self.top_plot.plot(region)
        self.bottom_plot.ax = self.bottom_ax
        self.bottom_plot.plot(region)
        self.bottom_ax.invert_yaxis()
        hide_axis(self.top_ax)
        hide_axis(self.bottom_ax)

        # if not hasattr(self.top_plot, 'colorbar') or self.top_plot.colorbar is None:
        #     sns.despine(ax=self.top_plot.cax, top=True, left=True, bottom=True, right=True)
        #     self.top_plot.cax.xaxis.set_visible(False)
        #     self.top_plot.cax.yaxis.set_visible(False)
        #
        # if not hasattr(self.bottom_plot, 'colorbar') or self.bottom_plot.colorbar is None:
        #     sns.despine(ax=self.bottom_plot.cax, top=True, left=True, bottom=True, right=True)
        #     self.bottom_plot.cax.xaxis.set_visible(False)
        #     self.bottom_plot.cax.yaxis.set_visible(False)

    def _refresh(self, region):
        self.top_plot.refresh(region)
        self.bottom_plot.refresh(region)

    def _clear(self):
        self.top_plot._clear()
        self.bottom_plot._clear()


VerticalSplitPlot = MirrorMatrixPlot


class GenomicFeaturePlot(BasePlotter1D):
    """
    Plot discrete genomic regions from BED, GTF files or similar.
    Just draws a black box where the feature is located.
    """
    def __init__(self, regions, feature_types=False, label_field="gene_symbol",
                 label_func=None, color='black', **kwargs):
        """
        :param regions: Any input that pybedtools can parse. Can be a path to a
                        GTF/BED file or a list of tuples [(2L, 500, 1000), (3R, 400, 600), ...]
        :param feature_types: If the input file is a GTF, only draw certain feature types (3rd column)
                              If False, draw all features on a common track
                              If None, automatically draw different feature types on separate tracks
                              If a list, draw only the feature types in the list on separate tracks,
                              don't draw the rest.
        :param label_field: Use this field as a label for each feature drawn.
                         Can be an integer to select a specific column or the name of an attribute
                         in the GTF file. If None or False no label is drawn. E.g. 2 for the second
                         column or "score" for the score attribute.
        :param label_func: Alternatively, label can be generated by calling this function which
                           takes pybedtools.Interval als argument and returns label string
        """
        kwargs.setdefault("aspect", .2)
        super(GenomicFeaturePlot, self).__init__(**kwargs)
        self.bedtool = parse_bedtool_input(regions)
        if feature_types is None and self.bedtool.file_type == "gff":
            feature_types = set(f[2] for f in self.bedtool)
        elif isinstance(feature_types, string_types):
            feature_types = [feature_types]
        self.feature_types = feature_types
        self._n_tracks = 1 if not self.feature_types else len(self.feature_types)
        self.label_field = label_field
        self.label_func = label_func
        self.color = color

    def _plot(self, region):
        interval = region_to_pbt_interval(region)
        genes = self.bedtool.all_hits(interval)
        trans = self.ax.get_xaxis_transform()
        pos = {k: i/(1 + self._n_tracks)
               for i, k in (enumerate(self.feature_types) if self.feature_types else [(0, "")])}
        stroke_length = max(0.03, 1/self._n_tracks - .05)
        for k, p in pos.items():
            self.ax.text(0, p, k, transform=self.ax.transAxes, ha="left", size=5)
        for g in genes:
            feature_type = g[2] if self.feature_types else ""
            if self.feature_types and feature_type not in self.feature_types:
                continue
            label = get_region_field(g, self.label_field, return_default="") if self.label_field else ""
            gene_patch = patches.Rectangle(
                (g.start, pos[feature_type]),
                width=abs(g.end - g.start), height=stroke_length,
                transform=trans, color=self.color
            )
            self.ax.add_patch(gene_patch)
            label_x = .5*(max(region.start, g.start) + min(region.end, g.end))
            self.ax.text(label_x, pos[feature_type] + stroke_length + .05,
                         label if not self.label_func else self.label_func(g),
                         transform=trans, ha="center", size="small")

        sns.despine(ax=self.ax, top=True, left=True, right=True)
        self.ax.yaxis.set_major_locator(NullLocator())
        self.remove_colorbar_ax()

    def _refresh(self, region):
        pass


class GenomicFeatureScorePlot(BasePlotter1D):
    """
    Plot discrete genomic regions from BED, GTF files or similar.

    Regions will be plotted as bars with the height equal to the score provided in the file.
    """
    def __init__(self, regions, attribute='score', feature_types=None,
                 color_neutral='grey', color_forward='red',
                 color_reverse='blue', annotation_field=None, ylim=None,
                 **kwargs):
        """
        :param regions: Any input that pybedtools can parse. Can be a path to a
                        GTF/BED file or a list of tuples [(2L, 500, 1000), (3R, 400, 600), ...]
        :param attribute: Field in the fiel which is plotted as score. Can be integer or attribute name
                          (if it is a GTF/GFF file)
        :param feature_types: If the input file is a GTF, only draw certain feature types (3rd column)
                              If False, draw all features on a common track
                              If None, automatically draw different feature types on separate tracks
                              If a list, draw only the feature types in the list on separate tracks,
                              don't draw the rest.
        """
        kwargs.setdefault("aspect", .2)
        kwargs.setdefault("axes_style", "ticks")
        super(GenomicFeatureScorePlot, self).__init__(**kwargs)
        self.regions = get_region_based_object(regions)
        if isinstance(feature_types, string_types):
            feature_types = [feature_types]
        self.feature_types = feature_types
        self.color_forward = color_forward
        self.color_reverse = color_reverse
        self.color_neutral = color_neutral
        self.attribute = attribute
        self.annotation_field = annotation_field
        self.ylim = ylim

        self._n_tracks = 1 if not self.feature_types else len(self.feature_types)

    def _plot(self, region):
        x = []
        y = []
        width = []
        colors = []
        annotations = []
        for f in self.regions.regions(region):
            if self.feature_types is not None:
                try:
                    if not f[2] in self.feature_types:
                        continue
                except ValueError:
                    pass

            x.append(f.start)
            width.append(f.end-f.start)
            if isinstance(self.attribute, int):
                score = f.fields[self.attribute]
            else:
                try:
                    score = getattr(f, self.attribute)
                except AttributeError:
                    score = 1
            try:
                score = float(score)
            except ValueError:
                score = 1
            y.append(score)

            if self.annotation_field is not None:
                if self.annotation_field != 'strand':
                    try:
                        annotation = getattr(f, self.annotation_field)
                        if annotation is None or annotation == '.':
                            annotation = ''
                        else:
                            try:
                                annotation = "{:.3f}".format(float(annotation))
                            except ValueError:
                                annotation = "{}".format(annotation)
                    except AttributeError:
                        annotation = ''

                else:
                    annotation = ''
                    if f.strand == '+' or f.strand == 1:
                        colors.append(self.color_forward)
                        annotation += ' (+)' if annotation != '' else '+'
                    elif f.strand == '-' or f.strand == -1:
                        colors.append(self.color_reverse)
                        annotation += ' (-)' if annotation != '' else '-'
                    else:
                        colors.append(self.color_neutral)
            else:
                annotation = ''
            annotations.append(annotation)

        for i in range(len(x)):
            x[i] += width[i]/2
        x = np.array(x)
        y = np.array(y)
        width = np.array(width)
        colors = np.array(colors)
        annotations = np.array(annotations)

        ixs = np.where(np.isfinite(y))[0]
        if len(colors) == len(x):
            colors = colors[ixs]
        else:
            colors = None
        y = y[ixs]
        x = x[ixs]
        width = width[ixs]
        annotations = annotations[ixs]

        rects = self.ax.bar(x, y, width=width, color=colors, edgecolor=colors, alpha=0.5)

        sns.despine(ax=self.ax, top=True, right=True)
        self.remove_colorbar_ax()
        if self.ylim is not None:
            self.ax.set_ylim(self.ylim)

        for i, rect in enumerate(rects):
            if i % 2 == 0:
                offset = 0.85
            else:
                offset = 1.0
            self.ax.text(rect.get_x() + rect.get_width() / 2., self.ax.get_ylim()[1] * offset,
                         annotations[i], ha='center', va='bottom', fontsize='x-small', alpha=0.5)
        # self.ax.yaxis.set_major_locator(NullLocator())
        # self.remove_colorbar_ax()

    def _refresh(self, region):
        pass


class RegionPlotBase(ScalarDataPlot):
    """
    Base class for region panels.
    """

    def __init__(self, data, plot_kwargs=None, labels=None, exclude=None, include=None,
                 bin_size=None, bins=None, **kwargs):
        """
        :param data: Data or list of data. Or dictionary, where keys represent
                     data labels. Data can be paths to files on the disk or anthing
                     that pybedtools can parse [(chr, start, end, score), ...].
                     If a list of data or dict is provided multiple lines are drawn.
                     Examples:
                     data=["x.bigwig", [("chr11", 40, 50, 1.8), ("chr11", 50, 70, 4.3)]]
                     data=["y.bedgraph", "z.bigwig"]
                     data={"x_chip": "x.bedgraph", "y_chip": "y.bigwig"}
        :param labels: List of labels for each data file. Used as label in the legend.
                       Ignored if labels are specified in data dictionary.
        :param bin_size: Bin values using fixed size bins of the given size.
                         If None, will plot values as they are in the data.
        :param fill: Fill space between x-axis and data line. Default: True
        :param plot_kwargs: Dictionary of additional keyword arguments passed to the plot function
        """
        kwargs.setdefault("aspect", .2)
        ScalarDataPlot.__init__(self, **kwargs)
        self.data = []
        self.labels = labels

        if isinstance(data, string_types):
            data = [data]
        # If data has attribute keys, assume it's dictionary
        elif hasattr(data, "keys"):
            if labels is not None:
                raise ValueError("Labels cannot not be assigned accurately if input is dict of data!")
            self.labels = list(data.keys())
            data = list(data.values())
        # First assume that input is an iterable with multiple datasets
        try:
            for d in data:
                self.data.append(load_score_data(d))
        except (ValueError, TypeError):
            # Assume input is a single data item
            self.data = [load_score_data(data)]
        self.plot_kwargs = {} if plot_kwargs is None else plot_kwargs
        self.exclude = exclude
        self.include = include
        self.bin_size = bin_size
        self.bins = bins

    def _region_valid(self, region):
        if self.exclude is not None:
            for key, value in self.exclude.items():
                a = getattr(region, key)
                if isinstance(value, string_types):
                    if a == string_types:
                        return False
                elif isinstance(value, set) or isinstance(value, list):
                    if a in value:
                        return False
        if self.include is not None:
            for key, value in self.include.items():
                a = getattr(region, key)
                if isinstance(value, string_types):
                    if a != string_types:
                        return False
                elif isinstance(value, set) or isinstance(value, list):
                    if a not in value:
                        return False
        return True

    def _region_iter(self, data, region=None, **kwargs):
        kwargs.setdefault('lazy', False)
        region_valid = self._region_valid
        if self.bin_size is not None or self.bins is not None:
            for region in data.binned_regions(region, bins=self.bins, bin_size=self.bin_size, **kwargs):
                if region_valid(region):
                    yield region
        else:
            for region in data.regions(region, **kwargs):
                if region_valid(region):
                    yield region

    def _data_iter(self, region, **kwargs):
        kwargs.setdefault('lazy', False)
        for data in self.data:
            yield self._region_iter(data, region, **kwargs)


class LinePlot(RegionPlotBase):
    """
    Plot data as line. Data must be :class:`~genomic_regions.RegionBased`
    """

    def __init__(self, data, fill=True, attribute='score',
                 colors=None, show_legend=None, legend_location='best', **kwargs):
        """
        :param data: Data or list of data. Or dictionary, where keys represent
                     data labels. Data can be paths to files on the disk or anything
                     that pybedtools can parse [(chr, start, end, score), ...].
                     If a list of data or dict is provided multiple lines are drawn.
                     Examples:
                     data=["x.bigwig", [("chr11", 40, 50, 1.8), ("chr11", 50, 70, 4.3)]]
                     data=["y.bedgraph", "z.bigwig"]
                     data={"x_chip": "x.bedgraph", "y_chip": "y.bigwig"}
        :param labels: List of labels for each data file. Used as label in the legend.
                       Ignored if labels are specified in data dictionary.
        :param bin_size: Bin values using fixed size bins of the given size.
                         If None, will plot values as they are in the data.
        :param fill: Fill space between x-axis and data line. Default: True
        :param plot_kwargs: Dictionary of additional keyword arguments passed to the plot function
        """
        super(LinePlot, self).__init__(data, **kwargs)
        self.lines = []
        self.fills = []
        self.fill = fill
        if isinstance(colors, string_types):
            colors = [colors]
        elif isinstance(colors, dict):
            if self.labels is None:
                raise ValueError("Colors can only be assigned as dict of labels are present")
            colors = [colors[l] for l in self.labels]
        elif colors is None:
            colors = ('red', 'blue', 'green', 'purple', 'yellow', 'black', 'orange', 'pink', 'cyan', 'lawngreen')

        self.colors = itertools.cycle(colors)
        self.attribute = attribute
        self.show_legend = show_legend
        self.legend_location = legend_location

    def _line_values(self, region):
        for i, region_iter in enumerate(self._data_iter(region)):
            x, y = self.values_from_region_iter(region_iter, self.attribute)
            yield i, x, y

    def _plot(self, region):
        limits = [999999999999999, 0]

        for i, x, y in self._line_values(region):
            limits[0] = min(limits[0], x[0])
            limits[1] = max(limits[1], x[-1])
            kwargs = self.plot_kwargs.copy()
            kwargs.setdefault('color', next(self.colors))
            kwargs.setdefault('label', self.labels[i] if self.labels else "")
            l = self.ax.plot(x, y, **kwargs)[0]

            self.lines.append(l)
            if self.fill:
                f = self.ax.fill_between(x, [0] * len(y), y, color=l.get_color(), alpha=l.get_alpha())
                self.fills.append(f)
        if self.labels is not None:
            if self.show_legend or (self.show_legend is None and len(self.labels) > 1):
                self.ax.legend(self.lines, self.labels, loc=self.legend_location)

        self.remove_colorbar_ax()
        sns.despine(ax=self.ax, top=True, right=True)
        
        if region.start is not None and region.end is not None:
            self.ax.set_xlim([region.start, region.end])

    def _refresh(self, region):
        for f in self.fills:
            f.remove()

        self.fills = []
        for i, x, y in self._line_values(region):
            self.lines[i].set_xdata(x)
            self.lines[i].set_ydata(y)
            if self.fill:
                f = self.ax.fill_between(x, [0] * len(y), y, color=self.lines[i].get_color(),
                                         alpha=self.lines[i].get_alpha())
                self.fills.append(f)


class BarPlot(RegionPlotBase):
    """
    Plot data as line. Data can be from BigWig or bedgraph files or anything pyBedTools can parse.
    """

    def __init__(self, data, labels=None, colors=None, alpha=0.5, min_score=None,
                 min_bar_width=0.0, attribute='score', show_legend=None,
                 legend_location='best', **kwargs):
        """
        :param data: Data or list of data. Or dictionary, where keys represent
                     data labels. Data can be paths to files on the disk or anthing
                     that pybedtools can parse [(chr, start, end, score), ...].
                     If a list of data or dict is provided multiple lines are drawn.
                     Examples:
                     data=["x.bigwig", [("chr11", 40, 50, 1.8), ("chr11", 50, 70, 4.3)]]
                     data=["y.bedgraph", "z.bigwig"]
                     data={"x_chip": "x.bedgraph", "y_chip": "y.bigwig"}
        :param labels: List of labels for each data file. Used as label in the legend.
                       Ignored if labels are specified in data dictionary.
        :param colors: List of matplotlib-compatible colors. If this list is shorter than
                       data, it will recycle from the beginning of the list
        :param alpha: float, alpha (transparency) value of each bar (useful for overlapping bars)
        :param min_score: Minimum score of a region to be plotted as a bar
        :param min_bar_width: Minimum plotting width of a bar in fraction of the plotting window.
                              Useful if regions in data are smaller than a pixel, making them
                              invisible.
        :param plot_kwargs: Dictionary of additional keyword arguments passed to the bar function
        """
        super(BarPlot, self).__init__(data, **kwargs)
        self.labels = labels
        self.lines = []

        if isinstance(colors, dict):
            if self.labels is None:
                raise ValueError("Colors can only be assigned as dict of labels are present")
            colors = [colors[l] for l in self.labels]
        if colors is None:
            colors = ('red', 'blue', 'green', 'purple', 'yellow', 'black', 'orange', 'pink', 'cyan', 'lawngreen')

        self.colors = itertools.cycle(colors)
        self.alpha = alpha
        self.min_score = min_score
        self.min_bar_width = min_bar_width
        self.attribute = attribute
        self.show_legend = show_legend
        self.legend_location = legend_location

    def _bar_values(self, region):
        min_width = self.min_bar_width * len(region)
        for i, region_iter in enumerate(self._data_iter(region)):
            x, w, h = [], [], []
            for region in region_iter:
                score = getattr(region, self.attribute)
                if self.min_score is not None and score < self.min_score:
                    continue
                x.append(region.start)
                width = len(region)
                if width < min_width:
                    width = min_width
                w.append(width)
                h.append(region.score)
            color = next(self.colors)
            yield x, w, h, color

    def _plot(self, region):
        bars = []
        for i, (x, w, h, c) in enumerate(self._bar_values(region)):
            kwargs = self.plot_kwargs.copy()
            kwargs.setdefault('color', c)
            kwargs.setdefault('align', 'edge')
            kwargs.setdefault('label', self.labels[i] if self.labels else "")
            kwargs.setdefault('alpha', self.alpha)

            b = self.ax.bar(x, h, w, **kwargs)
            if b is not None and len(b) > 0:
                bars.append(b[0])

        if self.labels is not None:
            if self.show_legend or (self.show_legend is None and len(self.labels) > 1):
                self.ax.legend(bars, self.labels, loc=self.legend_location)
        self.remove_colorbar_ax()
        sns.despine(ax=self.ax, top=True, right=True)
        
        if region.start is not None and region.end is not None:
            self.ax.set_xlim([region.start, region.end])

    def _refresh(self, region):
        self.ax.clear()
        self._plot(region)


class BigWigPlot(ScalarDataPlot):
    """
    Plot data from on or more BigWig or Bedgraph files. *Deprecated, use LinePlot instead*.
    """

    def __init__(self, bigwigs, names=None, bin_size=None, fill=True,
                 plot_kwargs=None, **kwargs):
        """
        :param bigwigs: Path or list of paths to bigwig or bedgraph files
        :param names: List of names for each bigwig. Used as label in the legend.
        :param bin_size: Bin BigWig values using fixed size bins of the given size.
                         If None, will plot values as they are in the BigWig file
        :param fill: Fill space between x-axis and data line. Default: True
        :param plot_kwargs: Dictionary of additional keyword arguments passed to the plot function
        """
        warnings.warn("BigWigPlot is deprecated, use LinePlot instead.")
        kwargs.setdefault("aspect", .2)
        super(BigWigPlot, self).__init__(**kwargs)
        if not isinstance(bigwigs, (list, tuple, types.GeneratorType)):
            bigwigs = [bigwigs]
        self.plot_kwargs = {} if plot_kwargs is None else plot_kwargs
        self.bigwigs = []
        for bw in bigwigs:
            if isinstance(bw, string_types):
                bw = fanc.load(bw)
            self.bigwigs.append(bw)
        self.names = names
        self.bin_size = bin_size
        self.lines = []
        self.fill = fill

    def _bin_intervals(self, region, intervals):
        """
        Bin values fearch interval from Bigwig in evenly spaced bins
        suitable for plotting.
        """
        bin_coords = np.r_[slice(region.start, region.end, self.bin_size), region.end]
        interval_records = np.core.records.fromrecords(intervals, names=["start", "end", "value"])
        out_values = np.full(len(bin_coords) - 1, np.nan, dtype=np.float_)
        start_overlap = np.searchsorted(interval_records["end"], bin_coords[:-1], side="right")
        end_overlap = np.searchsorted(interval_records["start"], bin_coords[1:], side="left")
        for i, (s, e) in enumerate(zip(start_overlap, end_overlap)):
            assert e >= s
            if s == e:
                # In this case no suitable values found in bigwig, leave nan
                continue
            weights = interval_records["end"][s:e] - interval_records["start"][s:e] + 1
            out_values[i] = np.average(interval_records["value"][s:e], weights=weights)
        bin_regions = [GenomicRegion(chromosome=region.chromosome, start=s, end=e)
                       for s, e in zip(bin_coords[:-1], bin_coords[1:])]
        return bin_regions, out_values

    def _line_values(self, region):
        for i, b in enumerate(self.bigwigs):
            if isinstance(b, fanc.data.genomic.RegionBased):
                intervals = b.region_intervals(region)
            else:
                intervals = b.intervals(region.chromosome, region.start - 1, region.end)

            if self.bin_size:
                regions, bw_values = self._bin_intervals(region, intervals)
            else:
                regions = [GenomicRegion(chromosome=region.chromosome, start=s, end=e) for s, e, v in intervals]
                bw_values = [v for s, e, v in intervals]
            x, y = self.get_plot_values(bw_values, regions)
            yield i, x, y

    def _plot(self, region):
        for i, x, y in self._line_values(region):
            l = self.ax.plot(x, y, label=self.names[i] if self.names else "",
                             **self.plot_kwargs)[0]
            self.lines.append(l)
            if self.fill:
                self.ax.fill_between(x, [0] * len(y), y, color=l.get_color())
        if self.names:
            self.add_legend()
        self.remove_colorbar_ax()
        sns.despine(ax=self.ax, top=True, right=True)

    def _refresh(self, region):
        for i, x, y in self._line_values(region):
            self.lines[i].set_xdata(x)
            self.lines[i].set_ydata(y)


class OldGenePlot(BasePlotter1D):
    """
    Plot genes including exon/intron structure from BED, GTF files or similar.
    """

    def __init__(self, genes, feature_types=('exon',), color_neutral='gray',
                 color_forward='orangered', color_reverse='darkturquoise',
                 color_score=False, vdist=0.2, box_height=0.1, show_labels=True,
                 label_field='name', font_size=9, arrow_size=8,
                 line_width=1, group_by='transcript_id', text_position='alternate',
                 collapse=False, squash=False, min_gene_size=None,
                 lookahead=2000000, **kwargs):
        """
        :param genes: Any input that pybedtools can parse. Can be a path to a
                      GTF/BED file
        :param feature_types: If the input file is a GTF, only draw certain feature types (3rd column)
                              If False, draw all features on a common track
                              If None, automatically draw different feature types on separate tracks
                              If a list, draw only the feature types in the list on separate tracks,
                              don't draw the rest.
        :param label_field: Field of input file for labelling transcripts/genes. Default: name
        :param collapse: Draw all transcripts on a single row. Everything will overlap. Default: False
        :param squash: Squash all exons belonging to a single grouping unit (merging overlapping exons).
                       Useful especially when setting group_by="gene_id" or "gene_symbol".
                       Genes will still draw on separate rows, if necessary. Default: False
        """
        kwargs.setdefault("aspect", .5)
        kwargs.setdefault("axes_style", "ticks")
        super(GenePlot, self).__init__(**kwargs)
        self.genes = get_region_based_object(genes)
        # ignore feature types if input is not GFF or GTF
        if self.genes.file_type != "gff" and self.genes.file_type != "gtf":
            feature_types = None

        self.color_score = color_score
        if color_score:
            scores = []
            for region in self.genes.regions:
                try:
                    scores.append(region.score)
                except ValueError:
                    scores.append(0.0)
            self.abs_max_score = max([min(scores), max(scores)])

        if isinstance(feature_types, string_types):
            feature_types = [feature_types]
        self.feature_types = feature_types
        self.color_forward = color_forward
        self.color_reverse = color_reverse
        self.color_neutral = color_neutral
        self.vdist = vdist
        self.box_height = box_height
        self.font_size = font_size
        self.group_by = group_by
        self.arrow_size = arrow_size
        self.text_position = text_position
        self.line_width = line_width
        self.show_labels = show_labels
        self.label_field = label_field
        self.collapse = collapse
        self.squash = squash
        self.min_gene_size = min_gene_size
        self.lookahead = lookahead

        self.lines = []
        self.patches = []
        self.texts = []

        self._n_tracks = 1 if not self.feature_types else len(self.feature_types)

    def _plot_genes(self, region):
        plot_range = region.end - region.start
        exon_hits = self.genes.regions(region.expand(absolute=self.lookahead))
        # trans = self.ax.get_xaxis_transform()

        gene_number = 0
        genes = defaultdict(list)
        for exon in exon_hits:
            if self.feature_types is not None:
                try:
                    if exon.feature not in self.feature_types:
                        continue
                except ValueError:
                    pass

            # get gene name
            try:
                name = getattr(exon, self.label_field)
            except AttributeError:
                name = None

            # get transcript id for grouping
            try:
                transcript_id = getattr(exon, self.group_by)
            except AttributeError:
                transcript_id = name

            if name is None and transcript_id is None:
                warnings.warn("Could not find either gene name or {}".format(self.group_by))
                name = str(gene_number)
                transcript_id = str(gene_number)
            elif name is None:
                name = transcript_id
            elif transcript_id is None:
                transcript_id = name

            exon_region = GenomicRegion(chromosome=region.chromosome, start=exon.start, end=exon.end,
                                        name=name, id=transcript_id, strand=exon.strand)

            genes[transcript_id].append(exon_region)
            gene_number = len(genes)

        # Squash transcripts
        if self.squash:
            for group_id, exons in genes.items():
                merged_exons = merge_overlapping_regions(exons)
                for exon in merged_exons:
                    exon.name = group_id
                genes[group_id] = merged_exons

        # sort exons
        for transcript_id, exons in genes.items():
            exons.sort(key=lambda x: (x.chromosome, x.start))

        # sort transcripts
        genes = [(name, exons) for name, exons in genes.items()]
        genes.sort(key=lambda x: x[1][0].start)

        genes_by_row = []
        for gene, exons in genes:
            # gene region - only used for calculating avoidance
            start = exons[0].start - 0.02 * plot_range
            end = exons[-1].end + 0.02 * plot_range
            gene_region = GenomicRegion(chromosome=region.chromosome, start=start, end=end)

            if self.collapse:
                if len(genes_by_row) == 0:
                    genes_by_row.append([(gene, gene_region, exons)])
                else:
                    genes_by_row[0].append((gene, gene_region, exons))
            else:
                # find empty spot in row
                spot_found = False
                for i, row in enumerate(genes_by_row):
                    overlaps = False
                    for row_gene in row:
                        if gene_region.overlaps(row_gene[1]):
                            overlaps = True
                            break
                    if not overlaps:
                        row.append((gene, gene_region, exons))
                        spot_found = True
                        break

                if not spot_found:
                    genes_by_row.append([(gene, gene_region, exons)])

        def _plot_gene(name, gene_region, exons, offset, text_position='top'):
            if exons[0].strand == 1:
                bar_marker = '$>$' if not self.invert_x else '$<$'
            elif exons[0].strand == -1:
                bar_marker = '$<$' if not self.invert_x else '$>$'
            else:
                bar_marker = 0

            gene_color = _set_gene_color(exons, self.color_score)

            bar_start, bar_end = exons[0].start, exons[-1].end
            bar_step_size = int(0.02 * plot_range)
            marker_correction = -1

            if bar_start == bar_end:
                bar_end += 1

            bar_x = list(range(bar_start, bar_end, bar_step_size))
            if bar_x[-1] != bar_end:
                bar_x += [bar_end]
                marker_correction -= 1
            bar_y = [offset] * len(bar_x)

            bar, = self.ax.plot(bar_x, bar_y, c=gene_color, linewidth=self.line_width)
            # transparent markers
            marker_bar, = self.ax.plot(bar_x[1:marker_correction], bar_y[1:marker_correction], marker=bar_marker,
                                       markersize=self.arrow_size, c=gene_color)
            self.lines.append(bar)
            self.lines.append(marker_bar)

            # plot exons
            for exon in exons:
                actual_width = exon.end - exon.start + 1
                if self.min_gene_size is None:
                    width = actual_width
                    exon_start = exon.start
                else:
                    width = max(self.min_gene_size, actual_width)
                    start_offset = (width - actual_width) / 2
                    exon_start = exon.start - start_offset
                width = width if not self.min_gene_size else max(self.min_gene_size, width)
                patch = self.ax.add_patch(
                    patches.Rectangle(
                        (exon_start, offset - self.box_height/2),  # (x,y)
                        width,  # width
                        self.box_height,  # height
                        facecolor=gene_color,
                        alpha=0.5,
                        edgecolor="none"
                    )
                )
                self.patches.append(patch)

            if self.show_labels and exons[0].start < region.end and exons[-1].end > region.start:
                if text_position == 'top':
                    text_y = offset - (self.box_height/2)*1.05
                    text_valign = 'bottom'
                elif text_position == 'bottom':
                    text_y = offset + (self.box_height / 2) * 1.05
                    text_valign = 'top'
                else:
                    raise ValueError("Text position '{}' not supported".format(text_position))

                text_x = max(exons[0].start, region.start)
                text_x = min(text_x, region.end)

                text = self.ax.text(text_x, text_y, exons[0].name,
                                    verticalalignment=text_valign, horizontalalignment='left',
                                    fontsize=self.font_size, family='monospace', color='gray')
                self.texts.append(text)

        def _set_gene_color(exons, color_score=False):
            cmap = plt.get_cmap('RdBu')
            if color_score:
                norm_score = (exons[0].score - (-self.abs_max_score)) / (self.abs_max_score - (-self.abs_max_score))
                gene_color = cmap(norm_score)
            else:
                if exons[0].strand == 1:
                    gene_color = self.color_forward
                elif exons[0].strand == -1:
                    gene_color = self.color_reverse
                else:
                    gene_color = self.color_neutral
            return gene_color

        for offset, row in enumerate(genes_by_row):
            for i, (name, gene_region, exons) in enumerate(row):
                if self.text_position == 'alternate':
                    if i % 2 == 0:
                        text_position = 'top'
                    else:
                        text_position = 'bottom'
                else:
                    text_position = self.text_position

                _plot_gene(name, gene_region, exons, offset*self.vdist, text_position)

        self.ax.set_ylim((len(genes_by_row)-1)*self.vdist+self.box_height/2*1.5, -1*self.box_height/2*1.5)

    def _plot(self, region):
        self._plot_genes(region=region)

        def drag_pan(self, button, key, x, y):
            mpl.axes.Axes.drag_pan(self, button, 'x', x, y)  # pretend key=='x'

        self.ax.drag_pan = types.MethodType(drag_pan, self.ax)
        self.ax.get_yaxis().set_visible(False)
        sns.despine(ax=self.ax, top=True, right=True, left=True)
        self.remove_colorbar_ax()

    def _refresh(self, region):
        while len(self.lines) > 0:
            el = self.lines.pop()
            el.remove()
            del el
        while len(self.patches) > 0:
            el = self.patches.pop()
            el.remove()
            del el
        while len(self.texts) > 0:
            el = self.texts.pop()
            el.remove()
            del el

        self._plot_genes(region)


class FeatureLayerPlot(BasePlotter1D):
    """
    Plot genomic features in layers grouped by name/type.::

        B1         -   -      -
        B2      -   -     - -   -
        L1       ---   ---     --
            _____________________
            0    10    20    30
    """

    def __init__(self, features, gff_grouping_attribute=None,
                 element_height=0.8, include=None, exclude=None,
                 color_by='strand', colors=((1, 'orangered'), (-1, 'darkturquoise')),
                 shadow=True, shadow_width=0.005,
                 collapse=False, **kwargs):
        """
        :param features: Any input that pybedtools can parse. Can be a path to a
                         GTF/BED file. If BED, elements will be grouped by name,
                         if GFF will be grouped by feature type
        :param gff_grouping_attribute: By default, GFF entries are grouped by feature type,
                                       change this to any attribute using this parameter
        :param element_height: Height of an individual element in the plot. A row's height
                               is 1.0, so you should choose a value smaller than that.
        :param color_by: element attribute to color the element by. Currently, only categorical
                         values are supported.
        :param colors: List of (attribute, color) pairs to color elements according to some attribute
        :param shadow: Draws a translucent box under each element the is min_element_width wide.
                       Useful if the size of elements is very small compared to plotting region
        :param shadow_width: Width of the shadow of an element in fraction of plotting region.
                             Some very small features won't be visible in this plot unless
                             you increase this parameter
        :param collapse: Collapse all rows onto a single one (ignore grouping)
        """
        kwargs.setdefault("aspect", 0.1)
        super(FeatureLayerPlot, self).__init__(**kwargs)
        self.features = get_region_based_object(features)
        if gff_grouping_attribute is None:
            self.grouping_attribute = 'feature' if hasattr(self.features, 'file_type') and \
                                                   self.features.file_type == 'gff' else 'name'
        else:
            self.grouping_attribute = gff_grouping_attribute
        self.element_height = element_height
        self.min_element_width = shadow_width
        self.draw_shadow = shadow
        self.top_offset = (1. - self.element_height) / 2 if self.element_height < 1. else 0
        self._patches = []
        self._color_by = color_by
        if colors is not None:
            self.colors = dict(colors)
        else:
            self.colors = dict()
        self._collapse = collapse
        self.include = set(include) if include is not None else None
        self.exclude = set(exclude) if exclude is not None else None

    def _plot_elements(self, region):
        groups = defaultdict(list)
        for element in self.features.regions(region):
            group = ' '
            try:
                if not self._collapse:
                    group = getattr(element, self.grouping_attribute)
            except AttributeError:
                pass

            if self.include is not None:
                if group not in self.include:
                    continue

            if self.exclude is not None:
                if group in self.exclude:
                    continue

            groups[group].append(element)

        region_width = region.end - region.start + 1
        tick_positions = []
        tick_labels = []
        for i, name in enumerate(sorted(groups)):
            y_offset = len(groups) - i
            tick_positions.append(y_offset - .5)
            tick_labels.append(name) if name != '.' else ''
            for element in groups[name]:
                try:
                    ca = getattr(element, self._color_by)
                    color = self.colors[ca]
                except KeyError:
                    color = 'grey'

                element_width = element.end - element.start + 1
                shadow_width = self.min_element_width * region_width
                if self.draw_shadow and shadow_width > element_width:
                    shadow_start = element.start - (shadow_width - element_width) / 2
                    shadow_patch = self.ax.add_patch(
                        patches.Rectangle(
                            (shadow_start, y_offset - self.top_offset - self.element_height),  # (x,y)
                            shadow_width,  # width
                            self.element_height,  # height
                            facecolor=color,
                            alpha=0.3,
                            edgecolor="none"
                        )
                    )
                    self._patches.append(shadow_patch)

                patch = self.ax.add_patch(
                    patches.Rectangle(
                        (element.start, y_offset - self.top_offset - self.element_height),  # (x,y)
                        element_width,  # width
                        self.element_height,  # height
                        facecolor=color,
                        alpha=1.,
                        edgecolor="none"
                    )
                )
                self._patches.append(patch)

        self.ax.set_ylim(0, len(groups))
        self.ax.set_yticks(tick_positions)
        self.ax.set_yticklabels(tick_labels)
        self.ax.tick_params(axis=u'y', which=u'both', length=0)

    def _plot(self, region):
        self._plot_elements(region)
        self.ax.spines['right'].set_visible(False)
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['left'].set_visible(False)
        self.ax.yaxis.set_ticks_position('left')
        self.ax.xaxis.set_ticks_position('bottom')
        self.remove_colorbar_ax()

    def _refresh(self, region):
        while len(self._patches) > 0:
            patch = self._patches.pop()
            patch.remove()
            del patch
        self._plot_elements(region)


class GenomicDataFramePlot(ScalarDataPlot):
    """
    Plot data from a table.
    """

    def __init__(self, genomic_data_frame, names=None,
                 plot_kwargs=None, **kwargs):
        """
        :param genomic_data_frame: :class:`~GenomicDataFrame`
        :param names: List of column names to plot (on the same axis)
        :param plot_kwargs: Dictionary of additional keyword arguments passed to the plot function
        """
        kwargs.setdefault("axes_style", style_ticks_whitegrid)
        kwargs.setdefault("aspect", .2)
        super(GenomicDataFramePlot, self).__init__(**kwargs)
        if isinstance(genomic_data_frame, string_types):
            genomic_data_frame = GenomicDataFrame.read_table(genomic_data_frame)
        self.plot_kwargs = {} if plot_kwargs is None else plot_kwargs
        self.genomic_data_frame = genomic_data_frame
        if names is None:
            names = genomic_data_frame.columns[3]
        if isinstance(names, string_types):
            names = [names]
        self.names = names
        self.lines = []

    def _draw_lines(self, region):
        x = []
        ys = [[] for _ in self.names]
        for r in self.genomic_data_frame.subset(region):
            x.append(r.center)
            for i, name in enumerate(self.names):
                value = getattr(r, name)
                ys[i].append(value)

        for i, name in enumerate(self.names):
            for line in self.ax.plot(x, ys[i], label=name, **self.plot_kwargs):
                self.lines.append(line)

    def _plot(self, region):
        self._draw_lines(region)

        self.remove_colorbar_ax()
        sns.despine(ax=self.ax, top=True, right=True)

    def _refresh(self, region):
        while len(self.lines) > 0:
            self.lines.pop(0).remove()
        plt.gca().set_prop_cycle(plt.matplotlib.rcParams['axes.prop_cycle'])
        self._draw_lines(region)


class Virtual4CPlot(BasePlotter1D):
    """
    Plot a 'virtual 4C' plot of the interactions of a specific genomic viewpoint as a 
    line plot. Extracts interactions from regions which overlap the viewpoint and 
    plots the mean of these interactions.

    This plot is an alternative to the :class:`~fanc.plotting.HicSlicePlot` with
    slightly different usage.
    """
    def __init__(self, hic, viewpoint, color='blue', alpha=1.0,
                 norm=True, oe=False, check_valid=True, excluded_filters=0, mask=True,
                 *args, **kwargs):
        """
        :param hic: :class:`~fanc.Hic` or :class:`~fanc.RegionMatrix`.
        :param viewpoint: Viewpoint to use for virtual 4C. String, e.g. 
                "2L:1000000-1500000" or :class:`~fanc.GenomicRegion`.
        :param color: Line colour to use for plotting. Default: blue.
        :param alpha: Transparency to use for plotting. Default: 1.0.
        :param norm: Use normalised values from Hi-C matrix. Default: True.
        :param oe: Use observed/expected values. Default: False.
        :param mask: Use values from masked Hi-C matrix. Default: True.
        :param check_valid: Parameter for Hi-C matrix extraction. Default: True.
        :param excluded_filters: Parameter for Hi-C matrix extraction. Default: 0.
        """
        BasePlotter1D.__init__(self, *args, **kwargs)
        self.hic = hic

        if isinstance(viewpoint, string_types):
            viewpoint = fanc.GenomicRegion.from_string(viewpoint)

        self.viewpoint = viewpoint
        self.lines = []
        self.color = color
        self.alpha = alpha
        self.norm = norm
        self.oe = oe
        self.check_valid = check_valid
        self.excluded_filters = excluded_filters
        self.mask = mask

    def _plot(self, region):
        submatrix = self.hic.matrix((self.viewpoint, region), norm=self.norm, oe=self.oe,
                                    check_valid=self.check_valid, excluded_filters=self.excluded_filters,
                                    mask=self.mask)
        v4c_signal = np.nanmean(submatrix, axis=0)
        x = []
        for r in self.hic.regions(region):
            x.append(r.center)

        line = self.ax.plot(x, v4c_signal, color=self.color, alpha=self.alpha)
        self.lines.append(line)
        self.remove_colorbar_ax()

    def _refresh(self, region):
        while len(self.lines) > 0:
            self.lines.pop(0).remove()

        self._plot(region)


class GenePlot(BasePlotter1D):
    """
    Plot genes including exon/intron structure from BED, GTF files or similar.
    """

    def __init__(self, genes, feature_types=('exon',), color_neutral='gray',
                 color_forward='orangered', color_reverse='darkturquoise',
                 color_score=False, box_height=0.1, show_labels=True,
                 label_field='name', font_size=9, arrow_size=8, show_arrows=True,
                 line_width=1, group_by='transcript_id',
                 collapse=False, squash=False, min_gene_size=None,
                 lookahead=2000000, score_colormap='RdBu', relative_text_offset=0.01,
                 relative_marker_step=0.03, no_labels_outside_plot=False,
                 include=None, exclude=None, **kwargs):
        """
        :param genes: Any input that pybedtools can parse. Can be a path to a
                      GTF/BED file
        :param feature_types: If the input file is a GTF, only draw certain feature types (3rd column)
                              If False, draw all features on a common track
                              If None, automatically draw different feature types on separate tracks
                              If a list, draw only the feature types in the list on separate tracks,
                              don't draw the rest.
        :param label_field: Field of input file for labelling transcripts/genes. Default: name
        :param collapse: Draw all transcripts on a single row. Everything will overlap. Default: False
        :param squash: Squash all exons belonging to a single grouping unit (merging overlapping exons).
                       Useful especially when setting group_by="gene_id" or "gene_symbol".
                       Genes will still draw on separate rows, if necessary. Default: False
        """
        kwargs.setdefault("aspect", .5)
        kwargs.setdefault("axes_style", "ticks")
        super(GenePlot, self).__init__(**kwargs)
        self.genes = get_region_based_object(genes)
        # ignore feature types if input is not GFF or GTF
        if self.genes.file_type != "gff" and self.genes.file_type != "gtf":
            feature_types = None

        self.color_score = color_score
        if color_score:
            scores = []
            for region in self.genes.regions:
                try:
                    scores.append(region.score)
                except ValueError:
                    scores.append(0.0)
            self.abs_max_score = max([min(scores), max(scores)])

        if isinstance(feature_types, string_types):
            feature_types = [feature_types]
        self.feature_types = feature_types
        self.color_forward = color_forward
        self.color_reverse = color_reverse
        self.color_neutral = color_neutral
        self.box_height = box_height
        self.font_size = font_size if font_size is not None else 9
        self.group_by = group_by
        self.arrow_size = arrow_size
        self.show_arrows = show_arrows
        self.line_width = line_width
        self.show_labels = self._normalise_include_exclude(show_labels)
        self.label_field = label_field
        self.collapse = collapse
        self.squash = squash
        self.min_gene_size = min_gene_size
        self.lookahead = lookahead
        self.score_colormap = score_colormap
        self.relative_text_offset = relative_text_offset
        self.relative_marker_step = relative_marker_step
        self.include = self._normalise_include_exclude(include)
        self.exclude = self._normalise_include_exclude(exclude)
        if isinstance(no_labels_outside_plot, bool):
            if no_labels_outside_plot:
                self.no_labels_outside_plot = 0.05
            else:
                self.no_labels_outside_plot = 2.0
        else:
            self.no_labels_outside_plot = float(no_labels_outside_plot)

        self.lines = []
        self.patches = []
        self.texts = []

        self._n_tracks = 1 if not self.feature_types else len(self.feature_types)

    def _normalise_include_exclude(self, ie_dict):
        if ie_dict is None or isinstance(ie_dict, bool):
            return ie_dict

        new_ie = dict()
        for key in ie_dict.keys():
            if isinstance(ie_dict[key], string_types):
                new_ie[key] = {ie_dict[key]}
            else:
                new_ie[key] = set(ie_dict[key])
        return new_ie

    def _plot_genes(self, region):
        fig = self.ax.figure

        plot_range = region.end - region.start
        exon_hits = self.genes.regions(region.expand(absolute=self.lookahead))
        # trans = self.ax.get_xaxis_transform()

        gene_number = 0
        genes = defaultdict(list)
        for i, exon in enumerate(exon_hits):
            if self.exclude is not None:
                exclude_exon = False
                for attribute_name, disallowed in self.exclude.items():
                    a = getattr(exon, attribute_name, None)
                    if a in disallowed:
                        exclude_exon = True
                        break

                if exclude_exon:
                    continue

            if self.include is not None:
                exclude_exon = False
                for attribute_name, allowed in self.include.items():
                    a = getattr(exon, attribute_name, None)
                    if a not in allowed:
                        exclude_exon = True
                        break

                if exclude_exon:
                    continue

            if self.feature_types is not None:
                try:
                    if exon.feature not in self.feature_types:
                        continue
                except ValueError:
                    pass

            # get gene name
            try:
                name = getattr(exon, self.label_field)
            except AttributeError:
                name = None

            # get transcript id for grouping
            try:
                transcript_id = getattr(exon, self.group_by)
            except AttributeError:
                transcript_id = i

            if name is None and transcript_id is None:
                warnings.warn("Could not find either gene name or {}".format(self.group_by))
                name = str(gene_number)
                transcript_id = str(gene_number)
            elif name is None:
                name = transcript_id
            elif transcript_id is None:
                transcript_id = name

            # exon_region = GenomicRegion(chromosome=region.chromosome, start=exon.start, end=exon.end,
            #                             name=name, id=transcript_id, strand=exon.strand)
            exon_region = exon.copy()
            exon_region.set_attribute('name', name)
            exon_region.set_attribute('id', transcript_id)

            genes[transcript_id].append(exon_region)
            gene_number = len(genes)

        plot_labels = dict()
        for gene, exons in genes.items():
            plot_labels[gene] = True
            if isinstance(self.show_labels, bool) and not self.show_labels:
                plot_labels[gene] = False
            elif isinstance(self.show_labels, dict):
                for attribute_name, allowed in self.show_labels.items():
                    a = getattr(exons[0], attribute_name, None)
                    if a is None or a not in allowed:
                        plot_labels[gene] = False
                        break

        # Squash transcripts
        if self.squash:
            for group_id, exons in genes.items():
                names = [exon.name for exon in exons]
                name = max(set(names), key=names.count)
                merged_exons = merge_overlapping_regions(exons)
                for exon in merged_exons:
                    exon.name = name
                genes[group_id] = merged_exons

        # sort exons
        for transcript_id, exons in genes.items():
            exons.sort(key=lambda x: (x.chromosome, x.start))

        # sort transcripts
        genes = [(name, exons) for name, exons in genes.items()]
        genes.sort(key=lambda x: x[1][0].start)

        # first plot all of the genes and labels in one place,
        # then figure out their non-overlapping positions
        def _plot_label(label):
            text = self.ax.text(0, 0, label,
                                verticalalignment='center', horizontalalignment='right',
                                fontsize=self.font_size,
                                family='monospace', color='gray')
            self.texts.append(text)
            return text

        def _set_gene_color(exons, color_score=False):
            cmap = plt.get_cmap(self.score_colormap)
            if color_score:
                norm_score = (exons[0].score - (-self.abs_max_score)) / (self.abs_max_score - (-self.abs_max_score))
                gene_color = cmap(norm_score)
            else:
                if exons[0].strand == 1:
                    gene_color = self.color_forward
                elif exons[0].strand == -1:
                    gene_color = self.color_reverse
                else:
                    gene_color = self.color_neutral
            return gene_color

        def _add_patch(x, y, width, height, color, alpha=1.0):
            return self.ax.add_patch(
                patches.Rectangle(
                    (x, y),  # (x,y)
                    width,  # width
                    height,  # height
                    facecolor=color,
                    alpha=alpha,
                    edgecolor="none"
                )
            )

        def _plot_gene(exons, offset):
            if self.show_arrows:
                if exons[0].is_forward():
                    bar_marker = '$>$' if not self.invert_x else '$<$'
                elif exons[0].is_reverse():
                    bar_marker = '$<$' if not self.invert_x else '$>$'
                else:
                    bar_marker = ''
            else:
                bar_marker = ''

            gene_color = _set_gene_color(exons, self.color_score)

            bar_start, bar_end = exons[0].start, exons[-1].end
            bar_step_size = int(self.relative_marker_step * plot_range)
            marker_correction = -1

            if bar_start == bar_end:
                bar_end += 1

            bar_x = list(range(bar_start, bar_end, bar_step_size))
            if bar_x[-1] != bar_end:
                bar_x += [bar_end]
                marker_correction -= 1
            bar_y = [offset] * len(bar_x)

            bar, = self.ax.plot(bar_x, bar_y, c=gene_color, linewidth=self.line_width)
            #bar = _add_patch()
            # transparent markers
            marker_bar, = self.ax.plot(bar_x[1:marker_correction], bar_y[1:marker_correction], marker=bar_marker,
                                       markersize=self.arrow_size, c=gene_color)
            self.lines.append(bar)
            self.lines.append(marker_bar)

            # plot exons
            exon_patches = []
            for exon in exons:
                actual_width = exon.end - exon.start + 1
                if self.min_gene_size is None:
                    width = actual_width
                    exon_start = exon.start
                else:
                    width = max(self.min_gene_size, actual_width)
                    start_offset = (width - actual_width) / 2
                    exon_start = exon.start - start_offset
                width = width if not self.min_gene_size else max(self.min_gene_size, width)
                patch = _add_patch(exon_start, offset - self.box_height / 2,
                                   width, self.box_height, gene_color, alpha=0.5)
                self.patches.append(patch)
                exon_patches.append(patch)

            # if self.show_labels and exons[0].start < region.end and exons[-1].end > region.start:
            #     _plot_label(exons[0].name)
            return [bar, marker_bar] + exon_patches

        gene_elements = []
        for gene, exons in genes:
            if exons[0].start < region.end and exons[-1].end > region.start:
                if plot_labels[gene]:
                    text = _plot_label(exons[0].name)
                else:
                    text = _plot_label('')

                bars = _plot_gene(exons, 0)
                gene_elements.append((text, bars))

        renderer = fig.canvas.get_renderer()

        genes_by_row = []
        inv = self.ax.transData.inverted()
        for text, bars in gene_elements:
            bb = text.get_window_extent(renderer=renderer)
            text_start, text_end = inv.transform(bb.intervalx)
            text_width = text_end - text_start

            x = bars[0].get_xdata(orig=True)
            gene_start, gene_end = x[0], x[-1]

            new_text_start = gene_start - self.relative_text_offset * plot_range
            start = new_text_start - text_width
            end = gene_end + self.relative_text_offset * plot_range

            gene_region = GenomicRegion(chromosome=region.chromosome, start=start, end=end)
            gene_fraction = (start - region.start) / len(region)
            if gene_fraction < -1*self.no_labels_outside_plot:
                text.set_text('')

            if self.collapse:
                offset = 0
            else:
                # find empty spot in row
                offset = -1
                for i, row in enumerate(genes_by_row):
                    overlaps = False
                    for row_gene in row:
                        if gene_region.overlaps(row_gene):
                            overlaps = True
                            break
                    if not overlaps:
                        row.append(gene_region)
                        offset = i
                        break

                if offset < 0:
                    offset = len(genes_by_row)
                    genes_by_row.append([gene_region])

            # move text and bars to new position(s)
            bars[0].set_ydata([offset] * len(bars[0].get_ydata()))
            bars[1].set_ydata([offset] * len(bars[1].get_ydata()))
            for bar in bars[2:]:
                bar.set_y(offset - self.box_height/2)
            text.set_x(new_text_start)
            text.set_y(offset)

        # if self.font_size is None:
        #     bbox = self.ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        #     font_size = (bbox.height*10 / (len(genes_by_row)-1))
        #     print(font_size, bbox.height, (bbox.y1-bbox.y0), self.ax.get_window_extent())
        #     for text, _ in gene_elements:
        #         text.set_size(font_size)

        self.ax.set_ylim(len(genes_by_row), -1)

    def _plot(self, region):
        self._plot_genes(region=region)

        def drag_pan(self, button, key, x, y):
            mpl.axes.Axes.drag_pan(self, button, 'x', x, y)  # pretend key=='x'

        self.ax.drag_pan = types.MethodType(drag_pan, self.ax)
        self.ax.get_yaxis().set_visible(False)
        sns.despine(ax=self.ax, top=True, right=True, left=True)
        self.remove_colorbar_ax()

    def _refresh(self, region):
        while len(self.lines) > 0:
            el = self.lines.pop()
            el.remove()
            del el
        while len(self.patches) > 0:
            el = self.patches.pop()
            el.remove()
            del el
        while len(self.texts) > 0:
            el = self.texts.pop()
            el.remove()
            del el

        self._plot_genes(region)
