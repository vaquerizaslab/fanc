from __future__ import division, print_function
from kaic.config import config
import matplotlib as mpl
from matplotlib.ticker import NullLocator
from kaic.data.genomic import GenomicRegion
from kaic.plotting.base_plotter import BasePlotter1D, ScalarDataPlot, BaseOverlayPlotter
from kaic.plotting.hic_plotter import BasePlotterMatrix
from kaic.plotting.helpers import append_axes, style_ticks_whitegrid, get_region_field, \
                                  region_to_pbt_interval, absolute_wspace_hspace
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import types
import numpy as np
import seaborn as sns
import pybedtools as pbt
try:
    from itertools import izip as zip
except ImportError:
    pass
import re
import warnings
from collections import defaultdict
from future.utils import string_types
import logging
logger = logging.getLogger(__name__)

plt = sns.plt


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
    def __init__(self, plots, height_ratios=None, figsize=None, hspace=.5,
                 gridspec_args=None, ticks_last=False, fix_chromosome=None,
                 invert_x=False, hide_x=None):
        """
        Creates a GenomicFigure composed of one or more plots.
        All plots are arranged in a single column, their genomic coordinates aligned.

        :param plots: List of plot instances, which should inherit
                      from :class:`~BasePlotter`
        :param height_ratios: Usually the aspect ratio of all plots is already specified
                              using the "aspect" argument when constructing the plot instances.
                              Alternatively, the aspect ratios can be overriden here by supplying
                              a list of aspect ratios, one for each plot.
                              The sum of this list is also used to calculate the total height of
                              the plot in inches height = width*sum(height_ratios).
                              If any or all entries are None, the aspect ratio of these
                              plots default to the value specified in the plot instances
                              using the aspect argument.
        :param figsize: Specify figure size directly (width, height) of figure in inches.
                        Defaults is (6, 6*sum(height_ratios))
                        None can be used to as a placeholder for the default value, eg.
                        (8, None) is converted to (8, 8*sum(height_rations))
        :param hspace: Distance between plot panels in inches
        :param gridspec_args: Optional keyword-arguments passed directly to GridSpec constructor
        :param ticks_last: Only draw genomic coordinate tick labels on last (bottom) plot
        :param fix_chromosome: boolean list, same length as plots. If an element is True, the corresponding plot
                               will receive a modified chromosome identifier (omitting or adding 'chr' as necessary)
        """
        self.plots = plots
        self.n = len(plots)
        self.ticks_last = ticks_last
        if height_ratios is None:
            height_ratios = [None]*self.n
        for i in range(self.n):
            if height_ratios[i] is None:
                height_ratios[i] = self.plots[i].get_default_aspect()
        self.height_ratios = height_ratios
        if figsize is None:
            figsize = (None, None)
        width, height = figsize
        if width is None:
            width = 6
        if height is None:
            height = width*sum(height_ratios) + hspace*self.n
        self.figsize = width, height
        if not gridspec_args:
            gridspec_args = {}
        gs = gridspec.GridSpec(self.n, 2, height_ratios=self.height_ratios, width_ratios=[1, .05], **gridspec_args)
        fig = plt.figure(figsize=self.figsize)
        absolute_wspace_hspace(fig, gs, .3, hspace)
        for i in range(self.n):
            with sns.axes_style("ticks" if plots[i].axes_style is None else
                                plots[i].axes_style):
                if i > 0:
                    ax = plt.subplot(gs[i, 0], sharex=self.axes[0])
                else:
                    ax = plt.subplot(gs[i, 0])
            plots[i].ax = ax
            plots[i].cax = plt.subplot(gs[i, 1])

        # fix chromosome identifiers
        if fix_chromosome is None:
            self.fix_chromosome = [False] * self.n
        else:
            self.fix_chromosome = fix_chromosome
        if len(self.fix_chromosome) != self.n:
            raise ValueError("fix_chromosome ({}) must be the same length "
                             "as plots ({})".format(len(self.fix_chromosome), self.n))

        self.invert_x = invert_x

        # hide x axes
        if hide_x is None:
            self.hide_x = [False] * self.n
        else:
            self.hide_x = hide_x
        if len(self.hide_x) != self.n:
            raise ValueError("hide_x ({}) must be the same length "
                             "as plots ({})".format(len(self.hide_x), self.n))

    @property
    def fig(self):
        return self.axes[0].figure
    
    def plot(self, region):
        """
        Make a plot of the specified region.
        :param region: A string describing a region "2L:10000000-12000000" or
                       a :class:`~GenomicRegion`
        :return: A matplotlib Figure instance and a list of figure axes
        """
        for i, (p, a) in enumerate(zip(self.plots, self.axes)):

            plot_region = region
            if self.fix_chromosome[i]:
                chromosome = region.chromosome
                if chromosome.startswith('chr'):
                    chromosome = chromosome[3:]
                else:
                    chromosome = 'chr' + chromosome
                plot_region = GenomicRegion(chromosome=chromosome, start=region.start, end=region.end)

            p.plot(plot_region, ax=a)
            if self.ticks_last and i < len(self.axes) - 1:
                p.remove_genome_ticks()

            if self.invert_x:
                a.invert_xaxis()

            if self.hide_x[i]:
                a.xaxis.set_visible(False)

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


class GenomicTrackPlot(ScalarDataPlot):
    def __init__(self, tracks, style="step", attributes=None, title='', aspect=.2,
                 axes_style=style_ticks_whitegrid):
        """
        Plot scalar values from one or more class:`~GenomicTrack` objects

        :param tracks: class:`~GenomicTrack`
        :param style: 'step' Draw values in a step-wise manner for each bin
              'mid' Draw values connecting mid-points of bins
        :param attributes: Only draw attributes from the track objects
                           which match this description.
                           Should be a list of names. Supports wildcard matching
                           and regex.
        :param title: Used as title for plot
        :param aspect: Default aspect ratio of the plot. Can be overriden by setting
               the height_ratios in class:`~GenomicFigure`
        """
        ScalarDataPlot.__init__(self, style=style, title=title, aspect=aspect,
                                axes_style=axes_style)
        if not isinstance(tracks, list):
            tracks = [tracks]
        self.tracks = tracks
        self.attributes = attributes
        self.lines = []

    def _plot(self, region=None, ax=None, *args, **kwargs):
        for track in self.tracks:
            bins = track.region_bins(region)
            values = track[bins]
            regions = track.regions[bins]
            for k, v in values.items():
                if not self.attributes or any(re.match(a.replace("*", ".*"), k) for a in self.attributes):
                    x, y = self.get_plot_values(v, regions)
                    l = self.ax.plot(x, y,
                                     label="{}{}".format(track.title + "_"
                                                         if track.title and len(self.tracks) > 1
                                                         else "", k))
                    self.lines.append(l[0])
        self.add_legend()
        self.remove_colorbar_ax()

    def _refresh(self, region=None, ax=None, *args, **kwargs):
        for track in self.tracks:
            bins = track.region_bins(region)
            values = track[bins]
            regions = track.regions[bins]
            current_line = 0
            for k, v in values.items():
                if not self.attributes or any(re.match(a.replace("*", ".*"), k) for a in self.attributes):
                    x, y = self.get_plot_values(v, regions)
                    self.lines[current_line].set_xdata(x)
                    self.lines[current_line].set_ydata(y)
                    current_line += 1


class GenomicRegionsPlot(ScalarDataPlot):
    def __init__(self, regions, style="step", attributes=None, title='', aspect=.2,
                 axes_style=style_ticks_whitegrid, ylim=None):
        """
        Plot scalar values from one or more class:`~GenomicTrack` objects

        :param regions: class:`~GenomicRegions`
        :param style: 'step' Draw values in a step-wise manner for each bin
              'mid' Draw values connecting mid-points of bins
        :param attributes: Only draw attributes from the track objects
                           which match this description.
                           Should be a list of names. Supports wildcard matching
                           and regex.
        :param title: Used as title for plot
        :param aspect: Default aspect ratio of the plot. Can be overriden by setting
               the height_ratios in class:`~GenomicFigure`
        """
        ScalarDataPlot.__init__(self, style=style, title=title, aspect=aspect,
                                axes_style=axes_style)

        self.regions = regions
        self.attributes = attributes
        self.lines = []
        self.ylim = ylim

    def _plot(self, region=None, ax=None, *args, **kwargs):
        for i, name in enumerate(self.regions.data_field_names):
            if not self.attributes or any(re.match("^" + a.replace("*", ".*") + "$", name) for a in self.attributes):
                regions = []
                values = []
                for r in self.regions.subset(region):
                    regions.append(r)
                    values.append(getattr(r, name))
                x, y = self.get_plot_values(values, regions)
                if self.regions.y_values is not None:
                    l = self.ax.plot(x, y, label="{}".format(self.regions.y_values[i]))
                else:
                    l = self.ax.plot(x, y, label="{}".format(name))
                self.lines.append(l[0])

        if self.ylim is not None:
            self.ax.set_ylim(self.ylim)

        self.add_legend()
        self.remove_colorbar_ax()

    def _refresh(self, region=None, ax=None, *args, **kwargs):
        current_line = 0
        for name in self.regions.data_field_names:
            if not self.attributes or any(re.match(a.replace("*", ".*"), name) for a in self.attributes):
                regions = []
                values = []
                for r in self.regions.subset(region):
                    regions.append(r)
                    values.append(getattr(r, name))
                x, y = self.get_plot_values(values, regions)
                self.lines[current_line].set_xdata(x)
                self.lines[current_line].set_ydata(y)
                current_line += 1


class GenomicMatrixPlot(BasePlotter1D, BasePlotterMatrix):
    def __init__(self, track, attribute, y_coords=None, y_scale='linear', plot_kwargs=None, title='',
                 colormap=config.colormap_hic, norm="lin", vmin=None, vmax=None,
                 show_colorbar=True, blend_zero=False,
                 unmappable_color=".9", illegal_color=None, aspect=.3,
                 axes_style="ticks"):
        """
        Plot matrix from a class:`~GenomicTrack` objects

        :param track: class:`~GenomicTrack` containing the matrix
        :param attribute: Which matrix from the track object to draw
        :param y_coords: Matrices in the class:`~GenomicTrack` object are
                         unitless. Can provide the coordinates for the
                         y-direction here. Matrix has shape (X, Y) must
                         have shape Y or Y + 1
        :param y_scale: Set scale of the y-axis, is passed to Matplotlib set_yscale, so any
                        valid argument ("linear", "log", etc.) works
        :param plot_kwargs: Keyword-arguments passed on to pcolormesh
        :param title: Used as title for plot
        :param aspect: Default aspect ratio of the plot. Can be overriden by setting
                       the height_ratios in class:`~GenomicFigure`
        :param blend_zero: If True then zero count bins will be drawn using the minimum
                   value in the colormap, otherwise transparent
        :param unmappable_color: Draw unmappable bins using this color. Defaults to
                                 light gray (".9")
        :param illegal_color: Draw non-finite (NaN, +inf, -inf) bins using this color. Defaults to
                         None (no special color).
        """
        BasePlotter1D.__init__(self, title=title, aspect=aspect, axes_style=axes_style)
        BasePlotterMatrix.__init__(self, colormap=colormap, norm=norm,
                                   vmin=vmin, vmax=vmax, show_colorbar=show_colorbar,
                                   blend_zero=blend_zero, unmappable_color=unmappable_color,
                                   illegal_color=illegal_color)
        self.track = track
        self.attribute = attribute
        if plot_kwargs is None:
            plot_kwargs = {}
        self.plot_kwargs = plot_kwargs
        self.y_coords = y_coords
        self.hm = None
        self.y_scale = y_scale

    def _plot(self, region=None, ax=None, *args, **kwargs):

        x, y, self.hm = self._mesh_data(region=region)
        self.collection = self.ax.pcolormesh(x, y, self.hm.T, rasterized=True, cmap=self.colormap,
                                             norm=self.norm, **self.plot_kwargs)
        self.collection._A = None
        self._update_mesh_colors()
        self.ax.set_yscale(self.y_scale)
        if self.y_coords is not None:
            self.ax.set_ylim(self.y_coords[0], self.y_coords[-1])

        if self.show_colorbar:
            self.add_colorbar()

    def _mesh_data(self, region):
        bins = self.track.region_bins(region)
        hm = self.track[bins][self.attribute]
        regions = self.track.regions()[bins]
        bin_coords = np.r_[[(x.start - 1) for x in regions], regions[-1].end]
        x, y = np.meshgrid(bin_coords, (self.y_coords if self.y_coords is not None
                                        else np.arange(hm.shape[1] + 1)))
        return x, y, hm

    def _update_mesh_colors(self):
        # pcolormesh doesn't support plotting RGB arrays directly like imshow, have to workaround
        # See https://github.com/matplotlib/matplotlib/issues/4277
        # http://stackoverflow.com/questions/29232439/plotting-an-irregularly-spaced-rgb-image-in-python/29232668?noredirect=1#comment46710586_29232668
        color_matrix = self.get_color_matrix(self.hm)
        color_tuple = color_matrix.transpose((1, 0, 2)).reshape(
            (color_matrix.shape[0] * color_matrix.shape[1], color_matrix.shape[2]))
        self.collection.set_color(color_tuple)

    def _refresh(self, region=None, *args, **kwargs):
        x, y, self.hm = self._mesh_data(region)

        self.collection._coordinates[:, :, 0] = x
        # update matrix data
        self.collection.set_array(self.hm.T.ravel())
        self._update_mesh_colors()


class GenomicVectorArrayPlot(BasePlotter1D, BasePlotterMatrix):
    def __init__(self, array, keys=None, y_coords=None, y_scale='linear', plot_kwargs=None, title='',
                 colormap=config.colormap_hic, colorbar_symmetry=None, norm="lin", vmin=None, vmax=None,
                 show_colorbar=True, blend_zero=True, replacement_color=None,
                 unmappable_color=".9", illegal_color=None, aspect=.3,
                 axes_style="ticks"):
        """
        Plot matrix from a class:`~GenomicTrack` objects

        :param array: class:`~MultiVectorArchitecturalRegionFeature`
        :param keys: keys for which vectors to use for array. None indicates all vectors will be used.
        :param y_coords: Matrices in the class:`~GenomicTrack` object are
                         unitless. Can provide the coordinates for the
                         y-direction here. Matrix has shape (X, Y) must
                         have shape Y or Y + 1
        :param y_scale: Set scale of the y-axis, is passed to Matplotlib set_yscale, so any
                        valid argument ("linear", "log", etc.) works
        :param plot_kwargs: Keyword-arguments passed on to pcolormesh
        :param title: Used as title for plot
        :param aspect: Default aspect ratio of the plot. Can be overriden by setting
                       the height_ratios in class:`~GenomicFigure`
        :param blend_zero: If True then zero count bins will be drawn using the minimum
                   value in the colormap, otherwise transparent
        :param unmappable_color: Draw unmappable bins using this color. Defaults to
                                 light gray (".9")
        :param illegal_color: Draw non-finite (NaN, +inf, -inf) bins using this color. Defaults to
                         None (no special color).
        """
        BasePlotter1D.__init__(self, title=title, aspect=aspect, axes_style=axes_style)
        BasePlotterMatrix.__init__(self, colormap=colormap, norm=norm, colorbar_symmetry=colorbar_symmetry,
                                   vmin=vmin, vmax=vmax, show_colorbar=show_colorbar,
                                   blend_zero=blend_zero, unmappable_color=unmappable_color,
                                   illegal_color=illegal_color, replacement_color=replacement_color)
        self.array = array
        self.keys = keys
        if plot_kwargs is None:
            plot_kwargs = {}
        self.plot_kwargs = plot_kwargs
        if y_coords is None:
            if self.array.y_values is not None:
                self.y_coords = self.array.y_values
            else:
                self.y_coords = None
        self.hm = None
        self.y_scale = y_scale

    def _plot(self, region=None, ax=None, *args, **kwargs):
        x, y, self.hm = self._mesh_data(region=region)
        self.collection = self.ax.pcolormesh(x, y, np.ma.masked_invalid(self.hm.T), rasterized=True, cmap=self.colormap,
                                             norm=self.norm, **self.plot_kwargs)

        self.collection._A = None
        self._update_mesh_colors()
        self.ax.set_yscale(self.y_scale)
        if self.y_coords is not None:
            self.ax.set_ylim(self.y_coords[0], self.y_coords[-1])

        if self.show_colorbar:
            self.add_colorbar()

        def drag_pan(self, button, key, x, y):
            mpl.axes.Axes.drag_pan(self, button, 'x', x, y)  # pretend key=='x'

        self.ax.drag_pan = types.MethodType(drag_pan, self.ax)

    def _mesh_data(self, region):
        hm = self.array.as_matrix(regions=region, keys=self.keys)
        regions = hm.regions
        bin_coords = np.r_[[(x.start - 1) for x in regions], regions[-1].end]
        x, y = np.meshgrid(bin_coords, (self.y_coords if self.y_coords is not None
                                        else np.arange(hm.shape[1] + 1)))
        return x, y, hm

    def _update_mesh_colors(self):
        # pcolormesh doesn't support plotting RGB arrays directly like imshow, have to workaround
        # See https://github.com/matplotlib/matplotlib/issues/4277
        # http://stackoverflow.com/questions/29232439/plotting-an-irregularly-spaced-rgb-image-in-python/29232668?noredirect=1#comment46710586_29232668
        color_matrix = self.get_color_matrix(np.ma.masked_invalid(self.hm))
        color_tuple = color_matrix.transpose((1, 0, 2)).reshape(
            (color_matrix.shape[0] * color_matrix.shape[1], color_matrix.shape[2]))
        self.collection.set_color(color_tuple)

    def _refresh(self, region=None, *args, **kwargs):
        x, y, self.hm = self._mesh_data(region)

        self.collection._coordinates[:, :, 0] = x
        # update matrix data
        self.collection.set_array(self.hm.T.ravel())
        self._update_mesh_colors()


class HicPeakPlot(BaseOverlayPlotter):
    """
    Overlay peaks onto Hicplot or HicPlot2D
    """
    def __init__(self, peaks, radius):
        """
        :param peaks: Kaic peaks instance
        :param radius: Radius in bp for plotted circles
        """
        BaseOverlayPlotter.__init__(self)
        self.peaks = peaks
        self.radius = radius
        self.circle_props = {"edgecolor": "black", "fill": False}

    @property
    def compatibility(self):
        return ["HicPlot", "HicPlot2D"]

    def _plot(self, base_plot, region):
        def plot_hicplot(start, end, radius):
            x = .5*(start + end)
            y = .5*(end - start)
            circle = patches.Circle((x, y), radius, **self.circle_props)
            base_plot.ax.add_patch(circle)

        def plot_hicplot2d(start, end, radius):
            circle = patches.Circle((start, end), radius, **self.circle_props)
            base_plot.ax.add_patch(circle)
            circle = patches.Circle((end, start), radius, **self.circle_props)
            base_plot.ax.add_patch(circle)

        plot_dispatch = {
            "HicPlot": plot_hicplot,
            "HicPlot2D": plot_hicplot2d
        }

        base_plot_class = base_plot.__class__.__name__
        peaks_gen = self.peaks.edge_subset((region, region))
        for p in peaks_gen:
            plot_dispatch[base_plot_class](p.source_node.start, p.sink_node.end, self.radius)


class VerticalSplitPlot(BasePlotter1D):
    """
    Stack two plots on top of each other, bottom plot inverted.
    Especially suited to stacking two Hic plots (triangles) on top
    of each other.
    """
    def __init__(self, top_plot, bottom_plot, gap=0, cax_gap=.05, title="", aspect=1.,
                 axes_style="ticks"):
        """
        Create split plot.

        :param top_plot: Plot instance on top
        :param bottom_plot: Plot instace on bottom
        :param gap: Gap between plots in inches
        :param cax_gap: Gap between colorbars in inches
        :param title: Used as title for plot
        :param aspect: Default aspect ratio of the plot. Can be overriden by setting
                       the height_ratios in class:`~GenomicFigure`
        """
        BasePlotter1D.__init__(self, title=title, aspect=aspect, axes_style=axes_style)
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

    def _plot(self, region=None, ax=None, *args, **kwargs):
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
        self.top_plot.plot(region, ax=self.top_ax)
        self.bottom_plot.plot(region, ax=self.bottom_ax)
        self.bottom_ax.invert_yaxis()
        hide_axis(self.top_ax)
        hide_axis(self.bottom_ax)

        if not hasattr(self.top_plot, 'colorbar') or self.top_plot.colorbar is None:
            sns.despine(ax=self.top_plot.cax, top=True, left=True, bottom=True, right=True)
            self.top_plot.cax.xaxis.set_visible(False)
            self.top_plot.cax.yaxis.set_visible(False)

        if not hasattr(self.bottom_plot, 'colorbar') or self.bottom_plot.colorbar is None:
            sns.despine(ax=self.bottom_plot.cax, top=True, left=True, bottom=True, right=True)
            self.bottom_plot.cax.xaxis.set_visible(False)
            self.bottom_plot.cax.yaxis.set_visible(False)

    def _refresh(self, region=None, ax=None, *args, **kwargs):
        pass

    def remove_genome_ticks(self):
        plt.setp(self.ax.get_xticklabels(), visible=False)
        self.ax.xaxis.offsetText.set_visible(False)


class GenomicFeaturePlot(BasePlotter1D):
    """
    Plot discrete genomic regions from BED, GTF files or similar.
    Just draws a black box where the feature is located.
    """
    def __init__(self, regions, title="", feature_types=False, label_field="gene_symbol",
                 label_func=None, aspect=.2, axes_style="ticks"):
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
        :param title: Used as title for plot
        :param aspect: Default aspect ratio of the plot. Can be overriden by setting
                       the height_ratios in class:`~GenomicFigure`
        """
        BasePlotter1D.__init__(self, title=title, aspect=aspect, axes_style=axes_style)
        if isinstance(regions, pbt.BedTool):
            self.bedtool = regions
        else:
            self.bedtool = pbt.BedTool(regions)
        if feature_types is None and self.bedtool.file_type == "gff":
            feature_types = set(f[2] for f in self.bedtool)
        elif isinstance(feature_types, string_types):
            feature_types = [feature_types]
        self.feature_types = feature_types
        self._n_tracks = 1 if not self.feature_types else len(self.feature_types)
        self.label_field = label_field
        self.label_func = label_func

    def _plot(self, region=None, ax=None, *args, **kwargs):
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
                transform=trans, color="black"
            )
            self.ax.add_patch(gene_patch)
            self.ax.text((g.start + g.end)/2, pos[feature_type] + stroke_length + .05,
                         label if not self.label_func else self.label_func(g),
                         transform=trans, ha="center", size="small")

        sns.despine(ax=self.ax, top=True, left=True, right=True)
        self.ax.yaxis.set_major_locator(NullLocator())
        self.remove_colorbar_ax()

    def _refresh(self, region=None, ax=None, *args, **kwargs):
        pass


class GenomicFeatureScorePlot(BasePlotter1D):
    """
    Plot discrete genomic regions from BED, GTF files or similar.

    Regions will be plotted as bars with the height equal to the score provided in the file.
    """
    def __init__(self, regions, title="", feature_types=None, aspect=.2, axes_style="ticks",
                 color_neutral='grey', color_forward='red', color_reverse='blue', show_labels=True):
        """
        :param regions: Any input that pybedtools can parse. Can be a path to a
                        GTF/BED file or a list of tuples [(2L, 500, 1000), (3R, 400, 600), ...]
        :param feature_types: If the input file is a GTF, only draw certain feature types (3rd column)
                              If False, draw all features on a common track
                              If None, automatically draw different feature types on separate tracks
                              If a list, draw only the feature types in the list on separate tracks,
                              don't draw the rest.
        :param title: Used as title for plot
        :param aspect: Default aspect ratio of the plot. Can be overriden by setting
                       the height_ratios in class:`~GenomicFigure`
        """
        BasePlotter1D.__init__(self, title=title, aspect=aspect, axes_style=axes_style)
        if not isinstance(regions, pbt.BedTool):
            self.bedtool = pbt.BedTool(regions)
        else:
            self.bedtool = regions
        if feature_types is None and self.bedtool.file_type == "gff":
            feature_types = set(f[2] for f in self.bedtool)
        elif isinstance(feature_types, string_types):
            feature_types = [feature_types]
        self.feature_types = feature_types
        self.color_forward = color_forward
        self.color_reverse = color_reverse
        self.color_neutral = color_neutral
        self.show_labels = show_labels

        self._n_tracks = 1 if not self.feature_types else len(self.feature_types)

    def _plot(self, region=None, ax=None, *args, **kwargs):
        interval = region_to_pbt_interval(region)
        features = self.bedtool.all_hits(interval)
        # trans = self.ax.get_xaxis_transform()

        x = []
        y = []
        width = []
        colors = []
        annotations = []
        for f in features:
            if self.feature_types is not None:
                try:
                    if not f[2] in self.feature_types:
                        continue
                except ValueError:
                    pass

            x.append(f.start)
            width.append(f.end-f.start)
            try:
                score = float(f.score)
            except ValueError:
                score = 1
            y.append(score)

            try:
                if f.name != '.':
                    annotation = f.name
                else:
                    annotation = ''
            except ValueError:
                annotation = ''

            if f.strand == '+' or f.strand == 1:
                colors.append(self.color_forward)
                annotation += ' (+)' if annotation != '' else '+'
            elif f.strand == '-' or f.strand == -1:
                colors.append(self.color_reverse)
                annotation += ' (-)' if annotation != '' else '-'
            else:
                colors.append(self.color_neutral)

            if not self.show_labels:
                annotation = ''

            annotations.append(annotation)

        rects = self.ax.bar(x, y, width=width, color=colors, edgecolor=colors)

        for i, rect in enumerate(rects):
            if i % 2 == 0:
                offset = 1.05
            else:
                offset = 1.25
            self.ax.text(rect.get_x() + rect.get_width() / 2., offset * rect.get_height(),
                         annotations[i], ha='center', va='bottom')

        sns.despine(ax=self.ax, top=True, right=True)
        self.remove_colorbar_ax()
        # self.ax.yaxis.set_major_locator(NullLocator())
        # self.remove_colorbar_ax()

    def _refresh(self, region=None, ax=None, *args, **kwargs):
        pass


class BigWigPlot(ScalarDataPlot):
    def __init__(self, bigwigs, names=None, style="step", title='', bin_size=None,
                 plot_kwargs=None, ylim=None, aspect=.2, axes_style=style_ticks_whitegrid):
        """
        Plot data from on or more BigWig files.

        :param bigwigs: Path or list of paths to bigwig files
        :param names: List of names for each bigwig. Used as label in the legend.
        :param style: 'step' Draw values in a step-wise manner for each bin
                      'mid' Draw values connecting mid-points of bins
        :param title: Title of the plot
        :param bin_size: Bin BigWig values using fixed size bins of the given size.
                         If None, will plot values as they are in the BigWig file
        :param plot_kwargs: Dictionary of additional keyword arguments passed to the plot function
        :param ylim: Tuple to set y-axis limits
        :param aspect: Default aspect ratio of the plot. Can be overriden by setting
               the height_ratios in class:`~GenomicFigure`
        """
        ScalarDataPlot.__init__(self, style=style, title=title, aspect=aspect,
                                axes_style=axes_style)
        if isinstance(bigwigs, string_types):
            bigwigs = [bigwigs]
        self.plot_kwargs = {} if plot_kwargs is None else plot_kwargs
        self.bigwigs = bigwigs
        self.names = names
        self.bin_size = bin_size
        self.ylim = ylim
        self.x = None
        self.y = None

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
            weights = interval_records["end"][s:e] - interval_records["start"][s:e]
            out_values[i] = np.average(interval_records["value"][s:e], weights=weights)
        bin_regions = [GenomicRegion(chromosome=region.chromosome, start=s, end=e)
                       for s, e in zip(bin_coords[:-1], bin_coords[1:])]
        return bin_regions, out_values

    def _plot(self, region=None, ax=None, *args, **kwargs):
        import pyBigWig

        for i, b in enumerate(self.bigwigs):
            if isinstance(b, string_types):
                try:
                    bw = pyBigWig.open(b)
                    intervals = bw.intervals(region.chromosome, region.start - 1, region.end)
                finally:
                    bw.close()
            else:
                intervals = b.intervals(region.chromosome, region.start - 1, region.end)

            if self.bin_size:
                regions, bw_values = self._bin_intervals(region, intervals)
            else:
                regions = [GenomicRegion(chromosome=region.chromosome, start=s, end=e) for s, e, v in intervals]
                bw_values = [v for s, e, v in intervals]
            self.x, self.y = self.get_plot_values(bw_values, regions)
            self.ax.plot(self.x, self.y, label=self.names[i] if self.names else "", **self.plot_kwargs)
        if self.names:
            self.add_legend()
        self.remove_colorbar_ax()
        sns.despine(ax=self.ax, top=True, right=True)
        if self.ylim:
            self.ax.set_ylim(self.ylim)

    def _refresh(self, region=None, ax=None, *args, **kwargs):
        pass


class GenePlot(BasePlotter1D):
    """
    Plot genes including exon/intron structure from BED, GTF files or similar.
    """
    def __init__(self, genes, title="", feature_types=('exon',), aspect=.5, axes_style="ticks",
                 color_neutral='gray', color_forward='orangered', color_reverse='darkturquoise',
                 color_score=False, vdist=0.2, box_height=0.1, show_labels=True, font_size=9, arrow_size=8,
                 line_width=1, group_by='transcript_id', text_position='alternate', collapse=False,
                 min_gene_size=None):
        """
        :param genes: Any input that pybedtools can parse. Can be a path to a
                      GTF/BED file
        :param feature_types: If the input file is a GTF, only draw certain feature types (3rd column)
                              If False, draw all features on a common track
                              If None, automatically draw different feature types on separate tracks
                              If a list, draw only the feature types in the list on separate tracks,
                              don't draw the rest.
        :param title: Used as title for plot
        :param aspect: Default aspect ratio of the plot. Can be overriden by setting
                       the height_ratios in class:`~GenomicFigure`
        """
        BasePlotter1D.__init__(self, title=title, aspect=aspect, axes_style=axes_style)
        if not isinstance(genes, pbt.BedTool):
            self.bedtool = pbt.BedTool(genes)
        else:
            self.bedtool = genes

        # ignore feature types if inout is not GFF or GTF
        if self.bedtool.file_type != "gff" and self.bedtool.file_type != "gtf":
            feature_types = None

        self.color_score = color_score
        if color_score:
            scores = []
            for interval in self.bedtool:
                try:
                    scores.append(float(interval.score))
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
        self.collapse = collapse
        self.min_gene_size = min_gene_size

        self.lines = []
        self.patches = []
        self.texts = []

        self._n_tracks = 1 if not self.feature_types else len(self.feature_types)

    def _plot_genes(self, region=None):
        plot_range = region.end - region.start
        interval = region_to_pbt_interval(region)
        exon_hits = self.bedtool.all_hits(interval)
        # trans = self.ax.get_xaxis_transform()

        gene_number = 0
        genes = defaultdict(list)
        for exon in exon_hits:
            if self.feature_types is not None:
                try:
                    if not exon[2] in self.feature_types:
                        continue
                except ValueError:
                    pass

            # get gene name
            try:
                name = exon.name
            except ValueError:
                name = None

            # get transcript id for grouping
            try:
                transcript_id = exon.attrs[self.group_by]
            except KeyError:
                transcript_id = name

            if name is None and transcript_id is None:
                warnings.warn("Could not find either gene name or {}".format(self.group_by))
                name = str(gene_number)
                transcript_id = str(gene_number)
            elif name is None:
                name = transcript_id
            elif transcript_id is None:
                transcript_id = name

            exon_region = GenomicRegion(chromosome=region.chromosome, start=exon.start + 1, end=exon.end,
                                        name=name, id=transcript_id, strand=exon.strand)

            # get gene score
            if self.color_score:
                try:
                    exon_region.score = float(exon.score)
                except ValueError:
                    exon_region.score = None

            genes[transcript_id].append(exon_region)
            gene_number = len(genes)

        # sort exons
        for transcript_id, exons in genes.items():
            exons.sort(key=lambda x: x.start)

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
                bar_marker = '$>$'
            elif exons[0].strand == -1:
                bar_marker = '$<$'
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

            if self.show_labels:
                if text_position == 'top':
                    text_y = offset - (self.box_height/2)*1.05
                    text_valign = 'bottom'
                elif text_position == 'bottom':
                    text_y = offset + (self.box_height / 2) * 1.05
                    text_valign = 'top'
                else:
                    raise ValueError("Text position '{}' not supported".format(text_position))

                text = self.ax.text(exons[0].start, text_y, exons[0].name, verticalalignment=text_valign,
                                    horizontalalignment='left', fontsize=self.font_size, family='monospace',
                                    color='gray')
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

    def _plot(self, region=None, ax=None, *args, **kwargs):
        self._plot_genes(region=region)

        def drag_pan(self, button, key, x, y):
            mpl.axes.Axes.drag_pan(self, button, 'x', x, y)  # pretend key=='x'

        self.ax.drag_pan = types.MethodType(drag_pan, self.ax)
        self.ax.get_yaxis().set_visible(False)
        sns.despine(ax=self.ax, top=True, right=True, left=True)
        self.remove_colorbar_ax()

    def _refresh(self, region=None, ax=None, *args, **kwargs):
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
