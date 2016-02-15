from matplotlib.ticker import MaxNLocator, ScalarFormatter
from matplotlib.widgets import Slider
import kaic
from kaic.data.genomic import GenomicRegion, RegionsTable, GenomicRegions
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from abc import abstractmethod, ABCMeta
import numpy as np
import math
import matplotlib as mpl
import logging
import seaborn as sns
import pybedtools as pbt
import itertools as it
import tables
import re
import track as track_module
import warnings
import copy
plt = sns.plt
log = logging.getLogger(__name__)
log.setLevel(10)

sns.set_style("ticks")

class SymmetricNorm(mpl.colors.Normalize):
    def autoscale(self, A):
        vmin = np.ma.min(A)
        vmax = np.ma.max(A)
        abs_max = max(abs(vmin), abs(vmax))
        self.vmin = -1.*abs_max
        self.vmax = abs_max

    autoscale_None = autoscale

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

def region_to_pbt_interval(region):
    return pbt.cbedtools.Interval(chrom=region.chromosome, start=region.start - 1, end=region.end)

def get_typed_array(input_iterable, nan_strings, count=-1):
    try:
        return np.fromiter((0 if x in nan_strings else x for x in input_iterable), int, count)
    except ValueError:
        pass
    try:
        return np.fromiter((np.nan if x in nan_strings else x for x in input_iterable), float, count)
    except ValueError:
        pass
    return np.fromiter(input_iterable, str, count)

class NewGenomicTrack(object):
    def __init__(self, file_name, title=None):
        """
        Initialize a genomic track.

        :param file_name: Genomic track file. Can be BED, BEDGRAPH, WIG, BIGWIG etc.
        :param title: The overall title of the track.
        """
        self.file_name = file_name
        self.title = title
        self.track = track.load(file_name)

    def __getitem__(self, key):
        if isinstance(key, int) or isinstance(key, slice):
            return {t.name: t[key] for t in self._tracks}
        if isinstance(key, basestring):
            region = GenomicRegion.from_string(key)
        elif isinstance(key, GenomicRegion):
            region = key

        return self[self.region_bins(region)]

class GenomicTrack(RegionsTable):
    def __init__(self, file_name, title=None, data_dict=None, regions=None, _table_name_tracks='tracks'):
        """
        Initialize a genomic track.

        :param file_name: Storage location of the genomic track HDF5 file
        :param data_dict: Dictionary containing data tracks as numpy arrays.
                          The arrays must have as many elements in the first
                          dimension as there are regions.
        :param title: The overall title of the track.
        :param regions: An iterable of (:class: `~kaic.data.genomic.GenomicRegion~)
                        or String elemnts that describe regions.
        """
        RegionsTable.__init__(self, file_name=file_name)
        # Check if file already exist, then don't allow to add more
        if regions:
            if not len(self._regions) > 0:
                self.add_regions(regions)
            else:
                warnings.warn("Existing GenomicTrack object already contains regions, ignoring regions to be added")
        if _table_name_tracks in self.file.root:
            self._tracks = self.file.get_node('/', _table_name_tracks)
        else:
            self._tracks = self.file.create_group('/', _table_name_tracks, "Genomic tracks")
        if title:
            self.title = title
        if data_dict:
            for k, v in data_dict.iteritems():
                self.add_data(k, v)

    @classmethod
    def from_gtf(cls, file_name, gtf_file, store_attrs=None, nan_strings=[".", ""]):
        """
        Import a GTF file as GenomicTrack.

        :param file_name: Storage location of the genomic track HDF5 file
        :param gtf_file: Location of GTF file_name
        :param store_attrs: List or listlike
                            Only store attributes in the list
        :param nan_strings: These characters will be considered NaN for parsing.
                            Will become 0 for int arrays, np.nan for float arrays
                            and left as is for string arrays.
        """
        gtf = pbt.BedTool(gtf_file)
        n = len(gtf)
        regions = []
        values = {}
        for i, f in enumerate(gtf.sort()):
            regions.append(GenomicRegion(chromosome=f.chrom, start=f.start, end=f.end, strand=f.strand))
            # If input is a GTF file, also store the type and source fields
            if f.file_type in ("gff", "gtf"):
                f.attrs["source"] = a.fields[1]
                f.attrs["feature"] = a.fields[2]
            # Check if there is a new attribute that hasn't occured before
            for k in f.attrs.keys():
                if not k in values and (not store_attrs or k in store_attrs):
                    if i > 0:
                        # Fill up values for this attribute with nan
                        values[k] = [nan_strings[0]]*i
                    else:
                        values[k] = []
            for k in values.keys():
                values[k].append(f.attrs.get(k, nan_strings[0]))
        for k, v in values.iteritems():
            values[k] = get_typed_array(v, nan_strings=nan_strings, count=n)
        return cls(file_name=file_name, data_dict=values, regions=regions)

    def to_bedgraph(self, prefix, tracks=None, skip_nan=True):
        if tracks is None:
            tracks = [t.name for t in self._tracks]
        elif not isinstance(tracks, list) and not isinstance(tracks, tuple):
            tracks = [tracks]
        for t in tracks:
            if self[t].ndim > 1:
                continue
            log.info("Writing track {}".format(t))
            with open("{}{}.bedgraph".format(prefix, t), "w") as f:
                for r, v in it.izip(self.regions, self[t]):
                    if skip_nan and np.isnan(v):
                        continue
                    f.write("{}\t{}\t{}\t{}\n".format(r.chromosome, r.start, r.end, v))

    def add_data(self, name, values, description=None):
        """
        Add a single genomic track to the object

        :param name: A string representing the name or name of the track
        :param values: A numpy array of values for each region in the object
        :param description: Longer description of track contents.
        """
        if values.shape[0] != len(self._regions):
            raise ValueError("First dimension of values must have as many elements "
                             "({}) as there are regions ({})".format(values.shape, len(self._regions)))
        self.file.create_array(self._tracks, name, values, description if description else "")

    def __getitem__(self, key):
        if isinstance(key, int) or isinstance(key, slice):
            return {t.name: t[key] for t in self._tracks}
        if isinstance(key, basestring):
            if key in self._tracks:
                return getattr(self._tracks, key)[:]
            region = GenomicRegion.from_string(key)
            return self[self.region_bins(region)]

    @property
    def tracks(self):
        # this = self
        # class TrackList(object):
        #     def __init__(self, track_object):
        #         self.track_object = track_object

        #     def __str__(self):
        #         return str(self.track_object._tracks)

        #     def __
        return {t.name: t[:] for t in self._tracks}

    @property
    def title(self):
        try:
            return self.file.get_node("/title").read()
        except tables.NoSuchNodeError:
            return None

    @title.setter
    def title(self, value):
        if value is None:
            self.file.remove_node("/title")
            return None
        try:
            self.file.create_array("/", "title", value)
        except tables.NodeError:
            self.file.remove_node("/title")
            self.file.create_array("/", "title", value)

class GenomicFigure(object):
    def __init__(self, plots, height_ratios=None, figsize=None, gridspec_args=None):
        self.plots = plots
        self.n = len(plots)
        if not gridspec_args:
            gridspec_args = {}
        gs = gridspec.GridSpec(self.n, 2, wspace=.1, hspace=.1, height_ratios=height_ratios, width_ratios=[1, .05], **gridspec_args)
        if figsize is None:
            figsize = (6, 6*self.n)
        self.axes = []
        fig = plt.figure(figsize=figsize)
        for i in xrange(self.n):
            if i > 0:
                ax = plt.subplot(gs[i, 0], sharex=self.axes[0])
            else:
                ax = plt.subplot(gs[i, 0])
            plots[i].cax = plt.subplot(gs[i, 1])
            self.axes.append(ax)
        #if figsize is None:
        #    figsize = (8, 4*self.n)
        #_, self.axes = plt.subplots(self.n, sharex=True, figsize=figsize)

    @property
    def fig(self):
        return self.axes[0].figure
    
    def plot(self, region):
        for p, a in zip(self.plots, self.axes):
            p.plot(region, ax=a)
        #self.fig.tight_layout()
        return self.fig, self.axes

    # def add_colorbar(self):
    #     vmin, vmax = float("inf"), float("-inf")
    #     for p in self.plots:
    #         if p.vmin < vmin:
    #             vmin = p.vmin
    #         if p.vmax > vmax:
    #             vmax = p.vmax
    #     cmap_data = mpl.cm.ScalarMappable(norm=self.norm, cmap=self.colormap)
    #     cmap_data.set_array([self.vmin, self.vmax])
    #     self.cax, kw = mpl.colorbar.make_axes(self.ax, location="top", shrink=0.4)
    #     self.colorbar = plt.colorbar(cmap_data, cax=self.cax, **kw)

    # @property
    # def norm(self):
    #     return self.p
    

class GenomeCoordFormatter(ScalarFormatter):
    def __init__(self, chromosome=None, start=None):
        ScalarFormatter.__init__(self, useOffset=False)
        if isinstance(chromosome, GenomicRegion):
            self.chromosome = chromosome.chromosome
            self.start = chromosome.start
        else:
            self.chromosome = chromosome
            self.start = start

    def __call__(self, x, pos=None):
        s = ScalarFormatter.__call__(self, x=x, pos=pos)
        if pos == 0 or x == 0:
            return "{}:{}".format(self.chromosome, s)
        return s

    def get_offset(self):
        """
        Returns little offset string that is written in bottom right corner
        of plot by default.
        """
        if len(self.locs) == 0:
            return ""
        s = ""
        if self.orderOfMagnitude:
            s = millify(10**self.orderOfMagnitude, precision=0)
        return self.fix_minus(s) + "b"

class GenomeCoordLocator(MaxNLocator):
    def __call__(self):
        vmin, vmax = self.axis.get_view_interval()
        ticks = self.tick_values(vmin, vmax)
        # Make sure that first and last tick are the start
        # and the end of the genomic range plotted. If next
        # ticks are too close, remove them.
        if ticks[0] - vmin < (vmax - vmin)/(self._nbins*3):
            ticks = ticks[1:]
        if vmax - ticks[-1] < (vmax - vmin)/(self._nbins*3):
            ticks = ticks[:-1]
        ticks = np.r_[vmin, ticks, vmax]
        return ticks

class BufferedMatrix(object):
    _STRATEGY_ALL = "all"
    _STRATEGY_FIXED = "fixed"
    _STRATEGY_RELATIVE ="relative"

    def __init__(self, data, buffering_strategy="relative", buffering_arg=1):
        self.data = data
        if buffering_strategy not in self._BUFFERING_STRATEGIES:
            raise ValueError("Only support the buffering strategies {}".format(self._BUFFERING_STRATEGIES.keys()))
        self.buffering_strategy = buffering_strategy
        self.buffering_arg = buffering_arg
        self.buffered_region = None
        self.buffered_matrix = None

    @classmethod
    def from_hic_matrix(cls, hic_matrix):
        bm = cls(data=None, buffering_strategy="all")
        bm.buffered_region = bm._STRATEGY_ALL
        bm.buffered_matrix = hic_matrix
        return bm

    def is_buffered_region(self, *regions):
        if (self.buffered_region is None or self.buffered_matrix is None or
                (not self.buffered_region == self._STRATEGY_ALL and not all(rb.contains(rq) for rb, rq in it.izip(self.buffered_region, regions)))):
            return False
        return True

    def get_matrix(self, *regions):
        regions = tuple(reversed([r for r in regions]))
        if not self.is_buffered_region(*regions):
            log.info("Buffering matrix")
            self._BUFFERING_STRATEGIES[self.buffering_strategy](self, *regions)
        return self.buffered_matrix[tuple(regions)]

    def _buffer_all(self, *regions):
        self.buffered_region = self._STRATEGY_ALL
        self.buffered_matrix = self.data[tuple([slice(None, None)]*len(regions))]

    def _buffer_relative(self, *regions):
        self.buffered_region = []
        for rq in regions:
            if rq.start is not None and rq.end is not None:
                rq_size = rq.end - rq.start
                new_start = max(1, rq.start - rq_size*self.buffering_arg)
                new_end = rq.end + rq_size*self.buffering_arg
                self.buffered_region.append(GenomicRegion(start=new_start, end=new_end, chromosome=rq.chromosome))
            else:
                self.buffered_region.append(GenomicRegion(start=None, end=None, chromosome=rq.chromosome))
        self.buffered_matrix = self.data[tuple(self.buffered_region)]

    def _buffer_fixed(self, *regions):
        self.buffered_region = []
        for rq in regions:
            if rq.start is not None and rq.end is not None:
                new_start = max(1, rq.start - self.buffering_arg)
                new_end = rq.end + self.buffering_arg
                self.buffered_region.append(GenomicRegion(start=new_start, end=new_end, chromosome=rq.chromosome))
            else:
                self.buffered_region.append(GenomicRegion(start=None, end=None, chromosome=rq.chromosome))
        self.buffered_matrix = self.data[tuple(self.buffered_region)]

    @property
    def buffered_min(self):
        return np.ma.min(self.buffered_matrix) if self.buffered_matrix is not None else None

    @property
    def buffered_max(self):
        return np.ma.max(self.buffered_matrix) if self.buffered_matrix is not None else None

    _BUFFERING_STRATEGIES = {_STRATEGY_ALL: _buffer_all, 
                             _STRATEGY_RELATIVE: _buffer_relative,
                             _STRATEGY_FIXED: _buffer_fixed}


class BufferedCombinedMatrix(BufferedMatrix):
    def __init__(self, hic_top, hic_bottom, scale_matrices=True, buffering_strategy="relative", buffering_arg=1):
        super(BufferedCombinedMatrix, self).__init__(None, buffering_strategy, buffering_arg)

        scaling_factor = 1
        if scale_matrices:
            scaling_factor = hic_top.scaling_factor(hic_bottom)

        class CombinedData(object):
            def __init__(self, hic_top, hic_bottom, scaling_factor=1):
                self.hic_top = hic_top
                self.hic_bottom = hic_bottom
                self.scaling_factor = scaling_factor

            def __getitem__(self, item):
                return hic_top.get_combined_matrix(self.hic_bottom, key=item, scaling_factor=self.scaling_factor)

        self.data = CombinedData(hic_top, hic_bottom, scaling_factor)

class BasePlotter(object):

    __metaclass__ = ABCMeta

    def __init__(self, title):
        self._ax = None
        self.cax = None
        self.title = title

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

    def __init__(self, title):
        BasePlotter.__init__(self, title=title)

    def plot(self, region=None, ax=None):
        if isinstance(region, basestring):
            region = GenomicRegion.from_string(region)
        if ax:
            self.ax = ax
        # set genome tick formatter
        self.ax.xaxis.set_major_formatter(GenomeCoordFormatter(region))
        self.ax.xaxis.set_major_locator(GenomeCoordLocator(nbins=10))
        self.ax.set_title(self.title)
        self._plot(region)
        self.ax.set_xlim(region.start, region.end)
        return self.fig, self.ax

class BasePlotterHic(object):

    __metaclass__ = ABCMeta

    def __init__(self, hic_data, colormap='viridis', norm="log",
                 vmin=None, vmax=None, show_colorbar=True, adjust_range=True,
                 buffering_strategy="relative", buffering_arg=1, blend_masked=False):
        self.hic_data = hic_data
        if isinstance(hic_data, kaic.Hic):
            self.hic_buffer = BufferedMatrix(hic_data, buffering_strategy=buffering_strategy, buffering_arg=buffering_arg)
        elif isinstance(hic_data, kaic.data.genomic.HicMatrix):
            self.hic_buffer = BufferedMatrix.from_hic_matrix(hic_data)
        else:
            raise ValueError("Unknown type for hic_data")
        self.colormap = copy.copy(mpl.cm.get_cmap(colormap))
        if blend_masked:
            self.colormap.set_bad(self.colormap(0))
        self._vmin = vmin
        self._vmax = vmax
        self.norm = prepare_normalization(norm=norm, vmin=vmin, vmax=vmax)
        self.colorbar = None
        self.slider = None
        self.show_colorbar = show_colorbar
        self.adjust_range = adjust_range

    def add_colorbar(self):
        cmap_data = mpl.cm.ScalarMappable(norm=self.norm, cmap=self.colormap)
        cmap_data.set_array([self.vmin, self.vmax])
        #self.cax, kw = mpl.colorbar.make_axes(self.ax, location="top", shrink=0.4)
        #self.cax, kw = mpl.colorbar.make_axes_gridspec(self.ax, orientation="horizontal",
        #                                      aspect=30, shrink=0.6)
        #divider = make_axes_locatable(self.ax)
        #self.cax = divider.append_axes("top", "5%", pad="5%")
        self.colorbar = plt.colorbar(cmap_data, cax=self.cax, orientation="vertical")

    def add_adj_slider(self):
        plot_position = self.cax.get_position()
        vmin_axs = plt.axes([plot_position.x0, 0.05, plot_position.width, 0.03], axisbg='#f3f3f3')
        self.vmin_slider = Slider(vmin_axs, 'vmin', self.hic_buffer.buffered_min,
                                  self.hic_buffer.buffered_max, valinit=self.vmin,
                                  facecolor='#dddddd', edgecolor='none')
        vmax_axs = plt.axes([plot_position.x0, 0.02, plot_position.width, 0.03], axisbg='#f3f3f3')
        self.vmax_slider = Slider(vmax_axs, 'vmax', self.hic_buffer.buffered_min,
                                  self.hic_buffer.buffered_max, valinit=self.vmax,
                                  facecolor='#dddddd', edgecolor='none')
        self.fig.subplots_adjust(top=0.90, bottom=0.15)
        self.vmin_slider.on_changed(self._slider_refresh)
        self.vmax_slider.on_changed(self._slider_refresh)

    def _slider_refresh(self, val):
        new_vmin = self.vmin_slider.val
        new_vmax = self.vmax_slider.val
        self.im.set_clim(vmin=new_vmin, vmax=new_vmax)
        self.colorbar.set_clim(vmin=new_vmin, vmax=new_vmax)
        self.colorbar.draw_all()

    @property
    def vmin(self):
        return self._vmin if self._vmin else self.hic_buffer.buffered_min

    @property
    def vmax(self):
        return self._vmax if self._vmax else self.hic_buffer.buffered_max
 
class BasePlotter2D(BasePlotter):

    __metaclass__ = ABCMeta

    def __init__(self, title):
        BasePlotter.__init__(self, title=title)
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
    def __init__(self, hic_data, title='', colormap='viridis', norm="log",
                 vmin=None, vmax=None, show_colorbar=True,
                 adjust_range=True, buffering_strategy="relative", buffering_arg=1,
                 blend_masked=False):
        BasePlotter2D.__init__(self, title=title)
        BasePlotterHic.__init__(self, hic_data=hic_data, colormap=colormap,
                                norm=norm, vmin=vmin, vmax=vmax, show_colorbar=show_colorbar,
                                adjust_range=adjust_range, buffering_strategy=buffering_strategy,
                                buffering_arg=buffering_arg, blend_masked=blend_masked)

    def _plot(self, x_region=None, y_region=None):
        m = self.hic_buffer.get_matrix(x_region, y_region)
        self.im = self.ax.imshow(m, interpolation='none', cmap=self.colormap, norm=self.norm, origin="upper",
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
        m = self.hic_buffer.get_matrix(x_region, y_region)

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
                 vmin=None, vmax=None, scale_matrices=True, show_colorbar=True,
                 buffering_strategy="relative", buffering_arg=1):
        super(HicComparisonPlot2D, self).__init__(hic_top, colormap=colormap, norm=norm, vmin=vmin, vmax=vmax,
                                                  show_colorbar=show_colorbar)
        self.hic_top = hic_top
        self.hic_bottom = hic_bottom
        self.hic_buffer = BufferedCombinedMatrix(hic_top, hic_bottom, scale_matrices, buffering_strategy, buffering_arg)


class HicPlot(BasePlotter1D, BasePlotterHic):
    def __init__(self, hic_data, title='', colormap='viridis', max_dist=None, norm="log",
                 vmin=None, vmax=None, show_colorbar=True, adjust_range=False,
                 buffering_strategy="relative", buffering_arg=1, blend_masked=False):
        BasePlotter1D.__init__(self, title=title)
        BasePlotterHic.__init__(self, hic_data, colormap=colormap, vmin=vmin, vmax=vmax,
                                show_colorbar=show_colorbar, adjust_range=adjust_range,
                                buffering_strategy=buffering_strategy, buffering_arg=buffering_arg,
                                norm=norm, blend_masked=blend_masked)
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
        X_ -= X_[1, 0] - (hm.row_regions[0].start - 1)
        Y_ -= .5*np.min(Y_) + .5*np.max(Y_)
        # create plot
        artist = self.ax.pcolormesh(X_, Y_, hm_masked, cmap=self.colormap, norm=self.norm)
        #import mpldatacursor
        #mpldatacursor.datacursor(artist)
        # set limits and aspect ratio
        self.ax.set_aspect(aspect="equal")
        self.ax.set_ylim(0, self.max_dist if self.max_dist else 0.5*(region.end-region.start))
        # remove y ticks
        self.ax.set_yticks([])
        # hide background patch
        self.ax.patch.set_visible(False)

        if self.show_colorbar:
            self.add_colorbar()
            if self.adjust_range:
                self.add_adj_slider()

    def _refresh(self, region=None):
        pass

class ScalarPlot(BasePlotter1D):
    def __init__(self, values, regions, title=''):
        BasePlotter1D.__init__(self, title=title)
        self.values = values
        self.regions = regions

    def _get_values_per_bp(self, region_list):
        v = np.empty(region_list[-1].end - region_list[0].start + 1)
        n = 0
        for r in region_list:
            v[n:n + r.end - r.start + 1] = self.values[r.ix]
            n += r.end - r.start + 1
        return v

    def _plot(self, region):
        region_list = list(self.regions.intersect(region))
        v = get_values_per_bp(region_list, self.values)
        self.ax.plot(np.arange(region_list[0].start, region_list[-1].end + 1), v)

    def _refresh(self):
        pass

class GenomicTrackPlot(BasePlotter1D):
    _STYLE_STEP = "step"
    _STYLE_MID = "mid"

    def __init__(self, tracks, style="step", attributes=None, title=''):
        BasePlotter1D.__init__(self, title=title)
        if not isinstance(tracks, list):
            tracks = [tracks]
        self.tracks = tracks
        self.attributes = attributes
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

    def _plot(self, region):
        for track in self.tracks:
            bins = track.region_bins(region)
            values = track[bins]
            regions = track.regions()[bins]
            for k, v in values.iteritems():
                if not self.attributes or any(re.match(a.replace("*", ".*"), k) for a in self.attributes):
                    x, y = self._STYLES[self.style](self, v, regions)
                    self.ax.plot(x, y, label="{}{}".format(track.title + "_" if track.title and len(self.tracks) > 1 else "", k))
        self.ax.legend()

    def _refresh(self):
        pass

    _STYLES = {_STYLE_STEP: _get_values_per_step,
               _STYLE_MID: _get_values_per_mid}

class GenomicArrayPlot(BasePlotter1D):
    def __init__(self, track, attribute, plot_kwargs=None, title=''):
        BasePlotter1D.__init__(self, title=title)
        self.track = track
        self.attribute = attribute
        if plot_kwargs is None:
            plot_kwargs = {}
        self.plot_kwargs = plot_kwargs

    def _plot(self, region):
        bins = self.track.region_bins(region)
        values = self.track[bins][self.attribute]
        regions = self.track.regions()[bins]
        bin_coords = np.r_[[(x.start - 1) for x in regions], (regions[-1].end)]
        X, Y = np.meshgrid(bin_coords, np.arange(values.shape[1] + 1))
        mesh = self.ax.pcolormesh(X, Y, values.T, **self.plot_kwargs)
        self.colorbar = plt.colorbar(mesh, cax=self.cax, orientation="vertical")
        #import mpldatacursor
        #mpldatacursor.datacursor(mesh)

    def _refresh(self):
        pass

class GeneModelPlot(BasePlotter1D):
    def __init__(self, gtf, title="", feature_types=None, id_field="gene_symbol"):
        import pybedtools as pbt
        BasePlotter1D.__init__(self, title=title)
        self.gtf = pbt.BedTool(gtf)
        if feature_types is None:
            feature_types = list(set(f[2] for f in self.gtf))
        elif isinstance(feature_types, (str, unicode)):
            feature_types = [feature_types]
        self.feature_types = feature_types
        self.id_field = id_field

    def _plot(self, region):
        interval = region_to_pbt_interval(region)
        genes = self.gtf.all_hits(interval)
        trans = self.ax.get_xaxis_transform()
        pos = {k: (i + 1)/len(self.feature_types) for i, k in enumerate(self.feature_types.iterkeys())}
        for g in genes:
            if not g[2]in self.feature_types:
                continue
            try:
                label = g.attrs[self.id_field]
            except KeyError:
                label = ""
            gene_patch = patches.Rectangle(
                (g.start, pos[g[2]]),
                width=abs(g.end - g.start), height=0.03,
                transform=trans, color="black"
            )
            self.ax.add_patch(gene_patch)
            self.ax.text((g.start + g.end)/2, pos[g[2]] + .035, label, transform=trans,
                         ha="center", size="small")
            self.ax.spines['right'].set_visible(False)
            self.ax.spines['top'].set_visible(False)
            self.ax.spines['left'].set_visible(False)
            self.ax.spines['bottom'].set_visible(False)
            self.ax.xaxis.set_ticks_position('bottom')
            self.ax.yaxis.set_visible(False)
            #self.ax.xaxis.set_visible(False)

    def _refresh(self):
        pass


class GenomicFeaturePlot(BasePlotter1D):
    def __init__(self, regions, labels, title="", color='black'):
        BasePlotter1D.__init__(self, title=title)
        sorted_regions = sorted(zip(regions, labels), key=lambda x: (x[0].chromosome, x[0].start))
        regions, self.labels = zip(*sorted_regions)

        self.color = color
        self.regions = GenomicRegions(regions=regions)

        for region in regions:
            print region

    def _plot(self, region):
        trans = self.ax.get_xaxis_transform()
        overlap_regions = self.regions.range(region)
        for r in overlap_regions:
            region_patch = patches.Rectangle(
                (r.start, 0.05),
                width=abs(r.end - r.start), height=0.6,
                transform=trans, color=self.color
            )
            self.ax.add_patch(region_patch)
            self.ax.text((r.start + r.end)/2, 0.8, self.labels[r.ix], transform=trans,
                         ha="center", size="small")

        self.ax.spines['right'].set_visible(False)
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['left'].set_visible(False)
        self.ax.spines['bottom'].set_visible(False)
        self.ax.xaxis.set_ticks_position('bottom')
        self.ax.yaxis.set_visible(False)
        self.ax.xaxis.set_visible(False)

    def _refresh(self):
        pass
