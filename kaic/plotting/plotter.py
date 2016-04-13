from __future__ import division, print_function
from matplotlib.ticker import NullLocator
from kaic.data.genomic import GenomicRegion
from kaic.plotting import style_ticks_whitegrid
from kaic.plotting.base_plotter import BasePlotter1D
from kaic.plotting.hic_plotter import BasePlotterMatrix
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import numpy as np
import logging
import seaborn as sns
import pybedtools as pbt
import itertools as it
import re

plt = sns.plt
log = logging.getLogger(__name__)
log.setLevel(10)


def region_to_pbt_interval(region):
    return pbt.cbedtools.Interval(chrom=region.chromosome, start=region.start - 1, end=region.end)


def get_region_field(interval, field, return_default=False):
    """
    Take BedTool region and return value stored in the specified field.
    Will try to fetch field from specific integer index, from Interval attribute and 
    lastly from the BedTool attributes dictionary present for GTF files

    :param return_default: If False, raise ValueError if field cannot be found. If anything
                           else return this value instead.
    """
    try:
        return interval[int(field)]
    except ValueError:
        pass
    try:
        return getattr(interval, field)
    except AttributeError:
        pass
    if interval.file_type == "gff":
        try:
            return interval.attrs[field]
        except KeyError:
            pass
    if return_default != False:
        return return_default
    else:
        raise ValueError("Field {} can't be found in inteval {}".format(field, interval))


def absolute_wspace_hspace(fig, gs, wspace, hspace):
    """
    Set distance between subplots of a GridSpec instance in inches. Updates the
    GridSpec instance and returns the calculated relative (as required by GridSpec) as tuple.

    :param fig: Figure instance
    :param gs: GridSpec instance
    :param wspace: Distance in inches horizontal
    :param hspace: Distance in inches vertical
    :return: (wspace, hspace) as a fraction of axes size
    """
    figsize = fig.get_size_inches()
    sp_params = gs.get_subplot_params(fig)
    wspace = wspace/figsize[0]
    hspace = hspace/figsize[1]
    tot_width = sp_params.right - sp_params.left
    tot_height = sp_params.top - sp_params.bottom
    nrows, ncols = gs.get_geometry()
    wspace = wspace*ncols/(tot_width - wspace*ncols + wspace)
    hspace = hspace*nrows/(tot_height - hspace*nrows + hspace)
    if not wspace > 0 or not hspace > 0:
        raise ValueError("Invalid relative spacing ({}, {}) calculated, "
                         "Probably distance set too large.".format(wspace, hspace))
    gs.update(wspace=wspace, hspace=hspace)
    return wspace, hspace


class GenomicFigure(object):
    def __init__(self, plots, height_ratios=None, figsize=None, gridspec_args=None,
                 ticks_last=False):
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
        :param gridspec_args: Optional keyword-arguments passed directly to GridSpec constructor
        :param ticks_last: Only draw genomic coordinate tick labels on last (bottom) plot
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
            height = width*sum(height_ratios)
        self.figsize = width, height
        if not gridspec_args:
            gridspec_args = {}
        gridspec_args["wspace"] = gridspec_args.get("wspace", .1)
        gridspec_args["hspace"] = gridspec_args.get("hspace", .2)
        gs = gridspec.GridSpec(self.n, 2, height_ratios=self.height_ratios, width_ratios=[1, .05], **gridspec_args)
        self.axes = []
        plt.figure(figsize=self.figsize)
        for i in xrange(self.n):
            with sns.axes_style("ticks" if plots[i].axes_style is None else
                                plots[i].axes_style):
                if i > 0:
                    ax = plt.subplot(gs[i, 0], sharex=self.axes[0])
                else:
                    ax = plt.subplot(gs[i, 0])

            if hasattr(plots[i], 'cax'):
                plots[i].cax = plt.subplot(gs[i, 1])
            else:
                cax = plt.subplot(gs[i, 1])
                sns.despine(ax=cax, top=True, left=True, bottom=True, right=True)
                cax.xaxis.set_visible(False)
                cax.yaxis.set_visible(False)
            plots[i].ax = ax
            self.axes.append(ax)

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
            p.plot(region, ax=a)
            if self.ticks_last and i < len(self.axes) - 1:
                p.remove_genome_ticks()
        return self.fig, self.axes

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        plt.close(self.fig)


class ScalarDataPlot(BasePlotter1D):
    """
    Base class for plotting scalar values like ChIP-seq signal.
    Provides methods for converting lists of values and regions
    to plotting coordinates.
    """
    _STYLE_STEP = "step"
    _STYLE_MID = "mid"

    def __init__(self, style="step", title='', aspect=.2, axes_style=style_ticks_whitegrid):
        BasePlotter1D.__init__(self, title=title, aspect=aspect, axes_style=axes_style)
        self.style = style
        if style not in self._STYLES:
            raise ValueError("Only the styles {} are supported.".format(self._STYLES.keys()))

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
            x[i] = int(round((r.end + r.start)/2))
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

    def remove_colorbar_ax(self):
        if self.cax is None:
            return
        try:
            self.fig.delaxes(self.cax)
        except KeyError:
            pass

    _STYLES = {_STYLE_STEP: _get_values_per_step,
               _STYLE_MID: _get_values_per_mid}


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
            for k, v in values.iteritems():
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
            for k, v in values.iteritems():
                if not self.attributes or any(re.match(a.replace("*", ".*"), k) for a in self.attributes):
                    x, y = self.get_plot_values(v, regions)
                    self.lines[current_line].set_xdata(x)
                    self.lines[current_line].set_ydata(y)
                    current_line += 1


class GenomicMatrixPlot(BasePlotter1D, BasePlotterMatrix):
    def __init__(self, track, attribute, y_coords=None, y_scale='linear', plot_kwargs=None, title='',
                 colormap='viridis', norm="lin", vmin=None, vmax=None,
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

    def _add_split_ax(self, ax, gap):
        bbox = ax.get_position()
        figsize = ax.figure.get_size_inches()
        gap = gap/figsize[1]
        top_ax = ax.figure.add_axes([bbox.x0, bbox.y0 + gap/2 + bbox.height/2,
                                     bbox.width, bbox.height/2 - gap/2])
        bottom_ax = ax.figure.add_axes([bbox.x0, bbox.y0,
                                        bbox.width, bbox.height/2 - gap/2], sharex=top_ax)
        return top_ax, bottom_ax

    def _plot(self, region=None, ax=None, *args, **kwargs):
        # Check if ax has already been split
        if self.parent_ax is not self.ax:
            self.parent_ax = self.ax
            self.top_ax, self.bottom_ax = self._add_split_ax(self.ax, self.gap)
            sns.despine(ax=self.ax, top=True, left=True, bottom=True, right=True)
            self.ax.xaxis.set_major_locator(NullLocator())
            self.ax.yaxis.set_major_locator(NullLocator())
        if self.parent_cax is not self.cax:
            self.parent_cax = self.cax
            self.top_plot.cax, self.bottom_plot.cax = self._add_split_ax(self.cax, self.cax_gap)
            self.cax.set_visible(False)
        self.top_plot.plot(region, ax=self.top_ax)
        self.bottom_plot.plot(region, ax=self.bottom_ax)
        self.bottom_plot.ax.invert_yaxis()
        sns.despine(ax=self.top_ax, top=True, left=True, bottom=True, right=True)
        self.top_ax.xaxis.set_major_locator(NullLocator())
        self.top_ax.xaxis.set_minor_locator(NullLocator())

    def _refresh(self, region=None, ax=None, *args, **kwargs):
        pass

    def remove_genome_ticks(self):
        plt.setp(self.bottom_ax.get_xticklabels(), visible=False)
        self.bottom_ax.xaxis.offsetText.set_visible(False)


class GenomicFeaturePlot(BasePlotter1D):
    """
    Plot discrete genomic regions from BED, GTF files or similar.
    Just draws a black box where the feature is located.
    """
    def __init__(self, regions, title="", feature_types=False, label_field="gene_symbol",
                 label_format=None, label_cast=None, aspect=.2, axes_style="ticks"):
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
        :param label_format: If a label field is specified, can also supply a python format string
                             to customize how the label is formatted. E.g. "{:.2}"
        :param label_cast: Supply a function to cast the label to a specific type, e.g. int or float
        :param title: Used as title for plot
        :param aspect: Default aspect ratio of the plot. Can be overriden by setting
                       the height_ratios in class:`~GenomicFigure`
        """
        BasePlotter1D.__init__(self, title=title, aspect=aspect, axes_style=axes_style)
        self.bedtool = pbt.BedTool(regions)
        if feature_types is None and self.bedtool.file_type == "gff":
            feature_types = set(f[2] for f in self.bedtool)
        elif isinstance(feature_types, (str, unicode)):
            feature_types = [feature_types]
        self.feature_types = feature_types
        self._n_tracks = 1 if not self.feature_types else len(self.feature_types)
        self.label_field = label_field
        self.label_format = label_format
        self.label_cast = label_cast

    def _plot(self, region=None, ax=None, *args, **kwargs):
        interval = region_to_pbt_interval(region)
        genes = self.bedtool.all_hits(interval)
        trans = self.ax.get_xaxis_transform()
        pos = {k: i/(1 + self._n_tracks)
               for i, k in (enumerate(self.feature_types) if self.feature_types else [(0, "")])}
        stroke_length = max(0.03, 1/self._n_tracks - .05)
        for k, p in pos.iteritems():
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
                         self.label_format.format(self.label_cast(label)
                                                  if self.label_cast
                                                  else label)
                         if self.label_format else label,
                         transform=trans, ha="center", size="small")
        sns.despine(ax=self.ax, top=True, left=True, right=True)
        self.ax.yaxis.set_major_locator(NullLocator())
        # self.remove_colorbar_ax()

    def _refresh(self, region=None, ax=None, *args, **kwargs):
        pass


class BigWigPlot(ScalarDataPlot):
    # TODO make this work - wWigIO just won't install
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
        if isinstance(bigwigs, basestring):
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
        bin_regions = [GenomicRegion(chromosome=region.chromosome, start=s, end=e)
                       for s, e in it.izip(bin_coords[:-1], bin_coords[1:])]
        interval_records = np.core.records.fromrecords(intervals, names=["start", "end", "value"])
        out_values = np.full(len(bin_coords) - 1, np.nan, dtype=np.float_)
        start_overlap = np.searchsorted(interval_records["start"], bin_coords[:-1], side="right") - 1
        end_overlap = np.searchsorted(interval_records["end"], bin_coords[1:], side="left")
        for i, (s, e) in enumerate(it.izip(start_overlap, end_overlap)):
            assert e >= s
            if s == e:
                out_values[i] = interval_records["value"][s]
                continue
            total_range = bin_coords[i + 1] - bin_coords[i]
            weighted_value = 0
            # Have to control for edge cases where first and/or last bin only partially overlap with
            # interval
            weighted_value += (min(interval_records["end"][s], bin_coords[i + 1]) -
                               max(interval_records["start"][s], bin_coords[i]))*interval_records["value"][s]
            weighted_value += (min(interval_records["end"][e], bin_coords[i + 1]) -
                               max(interval_records["start"][e], bin_coords[i]))*interval_records["value"][e]
            # Once edge case is taken care of the rest of the intervals can be binned evenly
            if e - s > 1:
                weighted_value += np.sum((interval_records["end"][s + 1:e] - interval_records["start"][s + 1:e]) *
                                         interval_records["value"][s + 1:e])
            out_values[i] = weighted_value/total_range
        return bin_regions, out_values

    def _plot(self, region=None, ax=None, *args, **kwargs):
        import wWigIO
        for i, b in enumerate(self.bigwigs):
            try:
                bw = wWigIO.open(b)
                intervals = wWigIO.getIntervals(b, region.chromosome, region.start - 1, region.end)
            finally:
                wWigIO.close(b)
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


class GenomicFeaturePlot(BasePlotter1D):
    def __init__(self, regions, labels, title="", color='black'):
        BasePlotter1D.__init__(self, title=title)
        sorted_regions = sorted(zip(regions, labels), key=lambda x: (x[0].chromosome, x[0].start))
        regions, self.labels = zip(*sorted_regions)

        self.color = color
        self.regions = GenomicRegions(regions=regions)

        for region in regions:
            print(region)

    def _plot(self, region=None, ax=None):
        trans = self.ax.get_xaxis_transform()
        overlap_regions = self.regions.range(region)
        for r in overlap_regions:
            print(r.ix)
            region_patch = patches.Rectangle(
                (r.start, 0.05),
                width=abs(r.end - r.start), height=0.6,
                transform=trans,
                color=self.color
            )
            self.ax.add_patch(region_patch)
            self.ax.text((r.start + r.end)/2, 0.8, self.labels[r.ix], transform=trans,
                         ha="center", size="small")

        # self.ax.spines['right'].set_visible(False)
        # self.ax.spines['top'].set_visible(False)
        # self.ax.spines['left'].set_visible(False)
        # self.ax.spines['bottom'].set_visible(False)
        # self.ax.xaxis.set_ticks_position('bottom')
        # self.ax.yaxis.set_visible(False)
        # self.ax.xaxis.set_visible(False)

    def _refresh(self, **kwargs):
        pass