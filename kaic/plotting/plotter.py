from __future__ import division, print_function
from kaic.config import config
import matplotlib as mpl
from matplotlib.ticker import NullLocator, MaxNLocator
from kaic import load
from kaic.data.genomic import GenomicRegion, GenomicDataFrame
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
import pyBigWig
import kaic
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
    """
    GenomicFigure composed of one or more plots.
    All plots are arranged in a single column, their genomic coordinates aligned.
    """

    def __init__(self, plots, height_ratios=None, figsize=None, hspace=.5,
                 gridspec_args=None, ticks_last=False, fix_chromosome=None,
                 invert_x=False, hide_x=None):
        """
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
        self.gs = gs
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
    """
    Plot scalar values from one or more class:`~GenomicTrack` objects
    """

    def __init__(self, tracks, attributes=None, **kwargs):
        """
        :param tracks: class:`~GenomicTrack`
        :param attributes: Only draw attributes from the track objects
                           which match this description.
                           Should be a list of names. Supports wildcard matching
                           and regex.
        """
        kwargs.setdefault("aspect", .2)
        kwargs.setdefault("axes_style", style_ticks_whitegrid)
        ScalarDataPlot.__init__(self, **kwargs)
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
    """
    Plot scalar values from one or more class:`~GenomicTrack` objects
    """

    def __init__(self, regions, attributes=None, names=None, **kwargs):
        """
        :param regions: class:`~GenomicRegions`
        :param attributes: Only draw attributes from the track objects
                           which match this description.
                           Should be a list of names. Supports wildcard matching
                           and regex.
        :param names: Supply list of names for each track.
        :param legend: Draw legend. Default: True
        """
        kwargs.setdefault("aspect", .2)
        ScalarDataPlot.__init__(self, **kwargs)

        self.regions = regions
        self.attributes = attributes
        self.lines = []
        self.names = names
        self.legend = legend

    def _plot(self, region=None, ax=None, *args, **kwargs):
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
                l = self.ax.plot(x, y, label=label)
                self.lines.append(l[0])
                line_counter += 1

        if self.legend:
            self.add_legend()
        self.remove_colorbar_ax()

    def _refresh(self, region=None, ax=None, *args, **kwargs):
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
    Plot scalar values from one or more class:`~GenomicTrack` objects
    """
    def __init__(self, regions, values, symmetry=None, **kwargs):
        """
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
        kwargs.setdefault("axes_style", style_ticks_whitegrid)
        kwargs.setdefault("aspect", .2)
        ScalarDataPlot.__init__(self, **kwargs)

        self.regions = regions
        self.legend = True
        if not isinstance(values, dict):
            self.values = {'data': values}
            self.legend = False
        else:
            self.values = values
        self.lines = []
        self.symmetry = symmetry

    def _plot_values(self, region=None):
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

    def _plot(self, region=None, ax=None, *args, **kwargs):
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

    def _refresh(self, region=None, ax=None, *args, **kwargs):
        for i, (label, x, y) in enumerate(self._plot_values(region)):
            self.lines[i].set_xdata(x)
            self.lines[i].set_ydata(y)


class GenomicMatrixPlot(BasePlotter1D, BasePlotterMatrix):
    """
    Plot matrix from a class:`~GenomicTrack` objects.
    """

    def __init__(self, track, attribute, y_coords=None, y_scale='linear', plot_kwargs=None,
                 **kwargs):
        """
        :param track: class:`~GenomicTrack` containing the matrix
        :param attribute: Which matrix from the track object to draw
        :param y_coords: Matrices in the class:`~GenomicTrack` object are
                         unitless. Can provide the coordinates for the
                         y-direction here. Matrix has shape (X, Y) must
                         have shape Y or Y + 1
        :param y_scale: Set scale of the y-axis, is passed to Matplotlib set_yscale, so any
                        valid argument ("linear", "log", etc.) works
        :param plot_kwargs: Keyword-arguments passed on to pcolormesh
        """
        kwargs.setdefault("aspect", .3)
        BasePlotter1D.__init__(self, **kwargs)
        BasePlotterMatrix.__init__(self, **kwargs)
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
    """
    Plot matrix from a class:`~MultiVectorArchitecturalRegionFeature` object.
    """

    def __init__(self, array, keys=None, y_coords=None, y_scale='linear', plot_kwargs=None,
                 **kwargs):
        """
        :param array: class:`~MultiVectorArchitecturalRegionFeature`
        :param keys: keys for which vectors to use for array. None indicates all vectors will be used.
        :param y_coords: Matrices in the class:`~GenomicTrack` object are
                         unitless. Can provide the coordinates for the
                         y-direction here. Matrix has shape (X, Y) must
                         have shape Y or Y + 1
        :param y_scale: Set scale of the y-axis, is passed to Matplotlib set_yscale, so any
                        valid argument ("linear", "log", etc.) works
        :param plot_kwargs: Keyword-arguments passed on to pcolormesh
        """
        kwargs.setdefault("aspect", .3)
        BasePlotter1D.__init__(self, **kwargs)
        BasePlotterMatrix.__init__(self, **kwargs)
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
        self.ax.set_ylim(self.hm.y_values[0], self.hm.y_values[-1])

        if self.show_colorbar:
            self.add_colorbar()

        def drag_pan(self, button, key, x, y):
            mpl.axes.Axes.drag_pan(self, button, 'x', x, y)  # pretend key=='x'

        self.ax.drag_pan = types.MethodType(drag_pan, self.ax)

    def _mesh_data(self, region):
        hm = self.array.as_matrix(regions=region, keys=self.keys)
        bin_coords = np.r_[[(x.start - 1) for x in hm.regions], hm.regions[-1].end]
        x, y = np.meshgrid(bin_coords, hm.y_values)
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
    def __init__(self, peaks, radius=None, **kwargs):
        """
        :param peaks: Kaic peaks instance
        :param radius: Radius in bp for plotted circles
        """
        BaseOverlayPlotter.__init__(self, **kwargs)
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
        peaks_gen = self.peaks.edge_subset((region, region), distances_in_bp=True)
        for p in peaks_gen:
            plot_dispatch[base_plot_class](p.source_node.start, p.sink_node.end, p.radius)


class VerticalSplitPlot(BasePlotter1D):
    """
    Stack two plots on top of each other, bottom plot inverted.
    Especially suited to stacking two Hic plots (triangles) on top
    of each other.
    """
    def __init__(self, top_plot, bottom_plot, gap=0, cax_gap=.05, **kwargs):
        """
        :param top_plot: Plot instance on top
        :param bottom_plot: Plot instace on bottom
        :param gap: Gap between plots in inches
        :param cax_gap: Gap between colorbars in inches
        """
        kwargs.setdefault("aspect", 1.)
        BasePlotter1D.__init__(self, **kwargs)
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


class GenomicFeaturePlot(BasePlotter1D):
    """
    Plot discrete genomic regions from BED, GTF files or similar.
    Just draws a black box where the feature is located.
    """
    def __init__(self, regions, feature_types=False, label_field="gene_symbol",
                 label_func=None, **kwargs):
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
        BasePlotter1D.__init__(self, **kwargs)
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
    def __init__(self, regions, attribute='score', feature_types=None, color_neutral='grey', color_forward='red',
                 color_reverse='blue', annotation_field=None, **kwargs):
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
        BasePlotter1D.__init__(self, **kwargs)
        if isinstance(regions, string_types):
            self.regions = kaic.load(regions)
        else:
            self.regions = regions

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

    def _plot(self, region=None, ax=None, *args, **kwargs):
        x = []
        y = []
        width = []
        colors = []
        annotations = []
        for f in self.regions[region]:
            if self.feature_types is not None:
                try:
                    if not f[2] in self.feature_types:
                        continue
                except ValueError:
                    pass

            x.append(f.start)
            width.append(f.end-f.start)
            try:
                try:
                    score = getattr(f, self.attribute)
                except AttributeError:
                    score = 1
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

    def _refresh(self, region=None, ax=None, *args, **kwargs):
        pass


class BigWigPlot(ScalarDataPlot):
    """
    Plot data from on or more BigWig files.
    """

    def __init__(self, bigwigs, names=None, bin_size=None, fill=True,
                 plot_kwargs=None, **kwargs):
        """
        :param bigwigs: Path or list of paths to bigwig files
        :param names: List of names for each bigwig. Used as label in the legend.
        :param bin_size: Bin BigWig values using fixed size bins of the given size.
                         If None, will plot values as they are in the BigWig file
        :param fill: Fill space between x-axis and data line. Default: True
        :param plot_kwargs: Dictionary of additional keyword arguments passed to the plot function
        """
        kwargs.setdefault("aspect", .2)
        ScalarDataPlot.__init__(self, **kwargs)
        if isinstance(bigwigs, string_types):
            bigwigs = [bigwigs]
        self.plot_kwargs = {} if plot_kwargs is None else plot_kwargs
        self.bigwigs = []
        for bw in bigwigs:
            if isinstance(bw, string_types):
                bw = kaic.load(bw)
            if isinstance(bw, kaic.BigWig):
                self.bigwigs.append(bw.bw)
            else:
                self.bigwigs.append(bw)
        self.names = names
        self.bin_size = bin_size
        self.lines = []
        self.condensed = condensed
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
            weights = interval_records["end"][s:e] - interval_records["start"][s:e]
            out_values[i] = np.average(interval_records["value"][s:e], weights=weights)
        bin_regions = [GenomicRegion(chromosome=region.chromosome, start=s, end=e)
                       for s, e in zip(bin_coords[:-1], bin_coords[1:])]
        return bin_regions, out_values

    def _line_values(self, region):
        for i, b in enumerate(self.bigwigs):
            intervals = b.intervals(region.chromosome, region.start - 1, region.end)

            if self.bin_size:
                regions, bw_values = self._bin_intervals(region, intervals)
            else:
                regions = [GenomicRegion(chromosome=region.chromosome, start=s, end=e) for s, e, v in intervals]
                bw_values = [v for s, e, v in intervals]
            x, y = self.get_plot_values(bw_values, regions)
            yield i, x, y

    def _plot(self, region=None, ax=None, *args, **kwargs):
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
        if self.condensed:
            low, high = self.ax.get_ylim()
            # self.ax.set_yticks([low, high])
            # self.ax.set_yticklabels([self.title, high], va='top', size='large')
            self.ax.set_yticks([high])
            self.ax.set_yticklabels([high], va='top', size='large')

    def _refresh(self, region=None, ax=None, *args, **kwargs):
        for i, x, y in self._line_values(region):
            self.lines[i].set_xdata(x)
            self.lines[i].set_ydata(y)


class GenePlot(BasePlotter1D):
    """
    Plot genes including exon/intron structure from BED, GTF files or similar.
    """

    def __init__(self, genes, title="", feature_types=('exon',), color_neutral='gray',
                 color_forward='orangered', color_reverse='darkturquoise',
                 color_score=False, vdist=0.2, box_height=0.1, show_labels=True, font_size=9, arrow_size=8,
                 line_width=1, group_by='transcript_id', text_position='alternate', collapse=False,
                 min_gene_size=None, **kwargs):
        """
        :param genes: Any input that pybedtools can parse. Can be a path to a
                      GTF/BED file
        :param feature_types: If the input file is a GTF, only draw certain feature types (3rd column)
                              If False, draw all features on a common track
                              If None, automatically draw different feature types on separate tracks
                              If a list, draw only the feature types in the list on separate tracks,
                              don't draw the rest.
        """
        kwargs.setdefault("aspect", .5)
        kwargs.setdefault("axes_style", "ticks")
        BasePlotter1D.__init__(self, **kwargs)
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


class FeatureLayerPlot(BasePlotter1D):
    """
    Plot genomic features in layers grouped by name/type.

        B1         -   -      -
        B2      -   -     - -   -
        L1       ---   ---     --
            _____________________
            0    10    20    30
    """

    def __init__(self, features, gff_grouping_attribute=None,
                 element_height=0.8, include=None, exclude=None,
                 color_by='strand', colors=((1, 'r'), (-1, 'b')),
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
                         values are supported
        :param colors: List of (attribute, color) pairs to color elements according to some attribute
        :param shadow: Draws a translucent box under each element the is min_element_width wide.
                       Useful if the size of elements is very small compared to plotting region
        :param shadow_width: Width of the shadow of an element in fraction of plotting region.
                             Some very small features won't be visible in this plot unless
                             you increase this parameter
        :param collapse: Collapse all rows onto a single one (ignore grouping)
        """
        kwargs.setdefault("aspect", 1.)
        BasePlotter1D.__init__(self, **kwargs)

        if isinstance(features, string_types):
            self.features = load(features)
        else:
            self.features = features
        if gff_grouping_attribute is None:
            self.grouping_attribute = 'feature' if self.features.file_type == 'gff' else 'name'
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
        elements = self.features[region]
        groups = defaultdict(list)
        for element in elements:
            if not self._collapse:
                group = getattr(element, self.grouping_attribute)
            else:
                group = ' '

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
            tick_labels.append(name)
            for element in groups[name]:
                try:
                    color = self.colors[element.strand]
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

    def _plot(self, region=None, ax=None, *args, **kwargs):
        self._plot_elements(region)
        self.ax.spines['right'].set_visible(False)
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['left'].set_visible(False)
        self.ax.yaxis.set_ticks_position('left')
        self.ax.xaxis.set_ticks_position('bottom')
        self.remove_colorbar_ax()

    def _refresh(self, region=None, ax=None, *args, **kwargs):
        while len(self._patches) > 0:
            patch = self._patches.pop()
            patch.remove()
            del patch
        self._plot_elements(region)


class GenomicDataFramePlot(ScalarDataPlot):
    def __init__(self, genomic_data_frame, names=None,
                 plot_kwargs=None, **kwargs):
        """
        Plot data from a table.

        :param genomic_data_frame: :class:`~GenomicDataFrame`
        :param names: List of column names to plot (on the same axis)
        :param plot_kwargs: Dictionary of additional keyword arguments passed to the plot function
        """
        kwargs.setdefault("axes_style", style_ticks_whitegrid)
        kwargs.setdefault("aspect", .2)
        ScalarDataPlot.__init__(self, **kwargs)
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

    def _plot(self, region=None, ax=None, *args, **kwargs):
        self._draw_lines(region)

        self.remove_colorbar_ax()
        sns.despine(ax=self.ax, top=True, right=True)

    def _refresh(self, region=None, ax=None, *args, **kwargs):
        while len(self.lines) > 0:
            self.lines.pop(0).remove()
        plt.gca().set_prop_cycle(plt.matplotlib.rcParams['axes.prop_cycle'])
        self._draw_lines(region)
