from __future__ import division, print_function
from matplotlib.ticker import MaxNLocator, Formatter, Locator
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
import warnings
import pyBigWig
plt = sns.plt
log = logging.getLogger(__name__)
log.setLevel(10)

sns.set_style("ticks")


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


def get_region_field(interval, field, return_none=False):
    """
    Take BedTool region and return value stored in the specified field.
    Will try to fetch field from specific integer index, from Interval attribute and 
    lastly from the BedTool attributes dictionary present for GTF files

    :param return_none: Instead of raising exception, return none if value can't be found
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
    if return_none:
        return None
    else:
        raise ValueError("Field {} can't be found in inteval {}".format(field, interval))


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
    def from_gtf(cls, file_name, gtf_file, store_attrs=None, nan_strings=(".", "")):
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
                f.attrs["source"] = f.fields[1]
                f.attrs["feature"] = f.fields[2]
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
                    f.write("{}\t{}\t{}\t{}\n".format(r.chromosome, r.start - 1, r.end, v))

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
            return
        try:
            self.file.create_array("/", "title", value)
        except tables.NodeError:
            self.file.remove_node("/title")
            self.file.create_array("/", "title", value)


class GenomicFigure(object):
    def __init__(self, plots, height_ratios=None, figsize=None, gridspec_args=None, ticks_last=False):
        """
        Creates a GenomicFigure composed of one or more plots.
        All plots are arranged in a single column, their genomic coordinates aligned.

        :param plots: List of plot instances, which should inherit
                      from :class:`~BasePlotter`
        :param height_ratios: A list of heights, one for each plot. The ratio of these
                              numbers represents the height ratios of the plots.
                              Also, the sum of this list is used as the height of
                              the plot in inches. The width is fixed at 6 inches at
                              the moment.
                              If any or all entries are None the aspect ratio of these
                              plots defaut to a plot-type specific value.
        :param figsize: Specify figure size directly (width, height) of figure in inches.
                        Defaults is (6, sum(height_rations))
                        None can be used to as a placeholder for the default value, eg.
                        (8, None) is converted to (8, sum(height_rations))
        :param gridspec_args: Optional keyword-arguments passed directly to GridSpec constructor
        :param ticks_last: Only draw genomic coordinate tick labels on last bottom plot
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
            if i > 0:
                ax = plt.subplot(gs[i, 0], sharex=self.axes[0])
            else:
                ax = plt.subplot(gs[i, 0])
            plots[i].cax = plt.subplot(gs[i, 1])
            plots[i].ax = ax
            self.axes.append(ax)

    @property
    def fig(self):
        return self.axes[0].figure
    
    def plot(self, region):
        for i, (p, a) in enumerate(zip(self.plots, self.axes)):
            p.plot(region, ax=a)
            # if self.height_ratios is not None:
            #     force_aspect(a, self.height_ratios[i])
            if self.ticks_last and i < len(self.axes) - 1:
                plt.setp(a.get_xticklabels(), visible=False)
                a.xaxis.offsetText.set_visible(False)
        return self.fig, self.axes

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        plt.close(self.fig)


class GenomeCoordFormatter(Formatter):
    """
    Process axis tick labels to give nice reprensations
    of genomic coordinates
    """
    def __init__(self, chromosome, display_scale=True):
        """
        :param chromosome: :class:`~kaic.GenomicRegion` or string
        :param display_scale: Boolean
                              Display distance scale at bottom right
        """
        if isinstance(chromosome, GenomicRegion):
            self.chromosome = chromosome.chromosome
        else:
            self.chromosome = chromosome
        self.display_scale = display_scale

    def _format_val(self, x, prec_offset=0):
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
        if not self.display_scale:
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


class BufferedMatrix(object):
    """
    Buffer contents of any :class:`~kaic.Hic` like objects. Matrix is
    prefetched and stored in memory. Buffer contents can quickly be fetched
    from memory. Different buffering strategies allow buffering of nearby
    regions so that adjacent parts of the matrix can quickly be fetched.
    """
    _STRATEGY_ALL = "all"
    _STRATEGY_FIXED = "fixed"
    _STRATEGY_RELATIVE = "relative"

    def __init__(self, data, buffering_strategy="relative", buffering_arg=1):
        """
        Initialize a buffer for Matrix-like objects that support
        indexing using class:`~GenomicRegion` objects, such as class:`~kaic.Hic`
        or class:`~kaic.RegionMatrix` objects.

        :param data: Data to be buffered
        :param buffering_strategy: "all", "fixed" or "relative"
                                   "all" buffers the whole matrix
                                   "fixed" buffers a fixed area, specified by buffering_arg
                                           around the query area
                                   "relative" buffers a multiple of the query area.
                                              With buffering_arg=1 the query area plus
                                              the same amount upstream and downstream
                                              are buffered
        :param buffering_arg: Number specifying how much around the query area is buffered
        """
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
                (not self.buffered_region == self._STRATEGY_ALL and not
                 all(rb.contains(rq) for rb, rq in it.izip(self.buffered_region, regions)))):
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
        self.buffered_matrix = self.data[tuple([slice(0, None, None)]*len(regions))]

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

    def __init__(self, title, aspect=1.):
        self._ax = None
        self.cax = None
        self.title = title
        self.has_legend = False
        self._aspect = aspect

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

    def remove_colorbar_ax(self):
        if self.cax is None:
            return
        try:
            self.fig.delaxes(self.cax)
        except KeyError:
            pass

    def add_legend(self, *args, **kwargs):
        if not self.has_legend:
            self.ax.legend(*args, **kwargs)

    def get_default_aspect(self):
        return self._aspect


class BasePlotter1D(BasePlotter):

    __metaclass__ = ABCMeta

    def __init__(self, title, aspect=1.):
        BasePlotter.__init__(self, title=title, aspect=aspect)

    def plot(self, region=None, ax=None):
        if isinstance(region, basestring):
            region = GenomicRegion.from_string(region)
        if ax:
            self.ax = ax
        # set genome tick formatter
        self.ax.xaxis.set_major_formatter(GenomeCoordFormatter(region))
        self.ax.xaxis.set_major_locator(GenomeCoordLocator(nbins=5))
        self.ax.xaxis.set_minor_locator(MinorGenomeCoordLocator(n=5))
        self.ax.set_title(self.title)
        self._plot(region)
        self.ax.set_xlim(region.start, region.end)
        return self.fig, self.ax


class BasePlotterMatrix(object):
    """
    Mix-in class to provide methods for mapping colorvalues
    in special areas in the plots etc.
    """

    __metaclass__ = ABCMeta

    def __init__(self, colormap='viridis', norm="log", vmin=None, vmax=None,
                 show_colorbar=True, blend_zero=True,
                 unmappable_color=".9", illegal_color=None):
        if isinstance(colormap, basestring):
            colormap = mpl.cm.get_cmap(colormap)
        self.colormap = colormap
        self._vmin = vmin
        self._vmax = vmax
        self.norm = prepare_normalization(norm=norm, vmin=vmin, vmax=vmax)
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
            unmappable = np.all(zero_mask, axis=0)
            color_matrix[unmappable, :] = mpl.colors.colorConverter.to_rgba(self.unmappable_color)
            color_matrix[:, unmappable] = mpl.colors.colorConverter.to_rgba(self.unmappable_color)
        return color_matrix

    def add_colorbar(self):
        """
        Add colorbar to the plot. Draws on the colorbar axes cax.
        """
        cmap_data = mpl.cm.ScalarMappable(norm=self.norm, cmap=self.colormap)
        cmap_data.set_array(np.array([self.vmin, self.vmax]))
        self.colorbar = plt.colorbar(cmap_data, cax=self.cax, orientation="vertical")

    @property
    def vmin(self):
        return self._vmin if self._vmin else self.norm.vmin

    @property
    def vmax(self):
        return self._vmax if self._vmax else self.norm.vmax


class BasePlotterHic(BasePlotterMatrix):

    __metaclass__ = ABCMeta

    def __init__(self, hic_data, colormap='viridis', norm="log",
                 vmin=None, vmax=None, show_colorbar=True, adjust_range=True,
                 buffering_strategy="relative", buffering_arg=1, blend_zero=True,
                 unmappable_color=".9", illegal_color=None):
        BasePlotterMatrix.__init__(self, colormap=colormap, norm=norm,
                                   vmin=vmin, vmax=vmax, show_colorbar=show_colorbar,
                                   blend_zero=blend_zero, unmappable_color=unmappable_color,
                                   illegal_color=illegal_color)
        self.hic_data = hic_data
        if isinstance(hic_data, kaic.Hic):
            self.hic_buffer = BufferedMatrix(hic_data, buffering_strategy=buffering_strategy,
                                             buffering_arg=buffering_arg)
        elif isinstance(hic_data, kaic.data.genomic.RegionMatrix):
            self.hic_buffer = BufferedMatrix.from_hic_matrix(hic_data)
        else:
            raise ValueError("Unknown type for hic_data")
        self.slider = None
        self.adjust_range = adjust_range

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


class BasePlotter2D(BasePlotter):

    __metaclass__ = ABCMeta

    def __init__(self, title, aspect=1.):
        BasePlotter.__init__(self, title=title, aspect=aspect)
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
                 blend_zero=True, unmappable_color=".9",
                 aspect=1.):
        """
        Initialize a 2D Hi-C heatmap plot.

        :param hic_data: class:`~kaic.Hic` or class:`~kaic.RegionMatrix`
        :param title: Title drawn on top of the figure panel
        :param colormap: Can be the name of a colormap or a Matplotlib colormap instance
        :param norm: Can be "lin", "log" or any Matplotlib Normalization instance
        :param vmin: Clip interactions below this value
        :param vmax: Clip interactions above this value
        :param show_colorbar: Draw a colorbar
        :param adjust_range: Draw a slider to adjust vmin/vmax interactively
        :param buffering_strategy: A valid buffering strategy for class:`~BufferedMatrix`
        :param buffering_arg: Adjust range of buffering for class:`~BufferedMatrix`
        :param blend_zero: If True then zero count bins will be drawn using the minimum
                           value in the colormap, otherwise transparent
        :param unmappable_color: Draw unmappable bins using this color. Defaults to
                                 light gray (".9")
        :param aspect: Default aspect ratio of the plot. Can be overriden by setting
                       the height_ratios in class:`~GenomicFigure`
        """
        BasePlotter2D.__init__(self, title=title, aspect=aspect)
        BasePlotterHic.__init__(self, hic_data=hic_data, colormap=colormap,
                                norm=norm, vmin=vmin, vmax=vmax, show_colorbar=show_colorbar,
                                adjust_range=adjust_range, buffering_strategy=buffering_strategy,
                                buffering_arg=buffering_arg, blend_zero=blend_zero,
                                unmappable_color=unmappable_color)

    def _plot(self, x_region=None, y_region=None):
        m = self.hic_buffer.get_matrix(x_region, y_region)
        self.im = self.ax.imshow(self.get_color_matrix(m), interpolation='none',
                                 cmap=self.colormap, norm=self.norm, origin="upper",
                                 extent=[m.col_regions[0].start, m.col_regions[-1].end,
                                         m.row_regions[-1].end, m.row_regions[0].start])
        self.last_ylim = self.ax.get_ylim()
        self.last_xlim = self.ax.get_xlim()

        if self.show_colorbar:
            self.add_colorbar()
            if self.adjust_range:
                self.add_adj_slider()

    def _refresh(self, x_region=None, y_region=None):
        print("refreshing")
        m = self.hic_buffer.get_matrix(x_region, y_region)
        self.im.set_data(self.get_color_matrix(m))
        self.im.set_extent([m.col_regions[0].start, m.col_regions[-1].end,
                            m.row_regions[-1].end, m.row_regions[0].start])


class HicSideBySidePlot2D(object):
    def __init__(self, hic1, hic2, colormap='viridis', norm="log",
                 vmin=None, vmax=None, aspect=1.):
        self.hic_plotter1 = HicPlot2D(hic1, colormap=colormap, norm=norm, vmin=vmin, vmax=vmax, aspect=aspect)
        self.hic_plotter2 = HicPlot2D(hic2, colormap=colormap, norm=norm, vmin=vmin, vmax=vmax, aspect=aspect)

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
                 buffering_strategy="relative", buffering_arg=1, aspect=1.):
        super(HicComparisonPlot2D, self).__init__(hic_top, colormap=colormap, norm=norm, vmin=vmin, vmax=vmax,
                                                  show_colorbar=show_colorbar, aspect=aspect)
        self.hic_top = hic_top
        self.hic_bottom = hic_bottom
        self.hic_buffer = BufferedCombinedMatrix(hic_top, hic_bottom, scale_matrices, buffering_strategy, buffering_arg)


class HicPlot(BasePlotter1D, BasePlotterHic):
    def __init__(self, hic_data, title='', colormap='viridis', max_dist=None, norm="log",
                 vmin=None, vmax=None, show_colorbar=True, adjust_range=False,
                 buffering_strategy="relative", buffering_arg=1, blend_zero=True,
                 unmappable_color=".9", illegal_color=None, aspect=.5):
        """
        Initialize a triangle Hi-C heatmap plot.

        :param hic_data: class:`~kaic.Hic` or class:`~kaic.RegionMatrix`
        :param title: Title drawn on top of the figure panel
        :param colormap: Can be the name of a colormap or a Matplotlib colormap instance
        :param norm: Can be "lin", "log" or any Matplotlib Normalization instance
        :param max_dist: Only draw interactions up to this distance
        :param vmin: Clip interactions below this value
        :param vmax: Clip interactions above this value
        :param show_colorbar: Draw a colorbar
        :param adjust_range: Draw a slider to adjust vmin/vmax interactively
        :param buffering_strategy: A valid buffering strategy for class:`~BufferedMatrix`
        :param buffering_arg: Adjust range of buffering for class:`~BufferedMatrix`
        :param blend_zero: If True then zero count bins will be drawn using the minimum
                           value in the colormap, otherwise transparent
        :param unmappable_color: Draw unmappable bins using this color. Defaults to
                                 light gray (".9")
        :param aspect: Default aspect ratio of the plot. Can be overriden by setting
                       the height_ratios in class:`~GenomicFigure`
        """
        BasePlotter1D.__init__(self, title=title, aspect=aspect)
        BasePlotterHic.__init__(self, hic_data, colormap=colormap, vmin=vmin, vmax=vmax,
                                show_colorbar=show_colorbar, adjust_range=adjust_range,
                                buffering_strategy=buffering_strategy, buffering_arg=buffering_arg,
                                norm=norm, blend_zero=blend_zero, unmappable_color=unmappable_color,
                                illegal_color=illegal_color)
        self.max_dist = max_dist

    def _plot(self, region=None):
        log.debug("Generating matrix from hic object")
        if region is None:
            raise ValueError("Cannot plot triangle plot for whole genome.")
        # Have to copy unfortunately, otherwise modifying matrix in buffer
        hm = self.hic_buffer.get_matrix(region, region)
        hm = kaic.data.genomic.RegionMatrix(np.copy(hm), col_regions=hm.col_regions, row_regions=hm.row_regions)
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
        # pcolormesh doesn't support plotting RGB arrays directly like imshow, have to workaround
        # See https://github.com/matplotlib/matplotlib/issues/4277
        # http://stackoverflow.com/questions/29232439/plotting-an-irregularly-spaced-rgb-image-in-python/29232668?noredirect=1#comment46710586_29232668
        color_matrix = self.get_color_matrix(hm)
        color_tuple = color_matrix.transpose((1,0,2)).reshape(
            (color_matrix.shape[0]*color_matrix.shape[1],color_matrix.shape[2]))
        self.collection = self.ax.pcolormesh(X_, Y_, hm, cmap=self.colormap, norm=self.norm, rasterized=True)
        self.collection._A = None
        self.collection.set_color(color_tuple)
        # set limits and aspect ratio
        #self.ax.set_aspect(aspect="equal")
        self.ax.set_ylim(0, self.max_dist/2 if self.max_dist else (region.end-region.start)/2)
        # remove outline everywhere except at bottom
        sns.despine(ax=self.ax, top=True, right=True, left=True)
        self.ax.set_yticks([])
        # hide background patch
        self.ax.patch.set_visible(False)
        if self.show_colorbar:
            self.add_colorbar()
            if self.adjust_range:
                self.add_adj_slider()

    def _refresh(self, region=None):
        pass


class ScalarDataPlot(BasePlotter1D):
    """
    Base class for plotting scalar values like ChIP-seq signal.
    Provides methods for converting lists of values and regions
    to plotting coordinates.
    """
    _STYLE_STEP = "step"
    _STYLE_MID = "mid"

    def __init__(self, style="step", title='', aspect=.2):
        BasePlotter1D.__init__(self, title=title, aspect=aspect)
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

    _STYLES = {_STYLE_STEP: _get_values_per_step,
               _STYLE_MID: _get_values_per_mid}


class BigWigPlot(ScalarDataPlot):
    def __init__(self, bigwigs, names=None, style="step", title='', bin_size=None,
                 plot_kwargs=None, ylim=None, aspect=.2):
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
        ScalarDataPlot.__init__(self, style=style, title=title, aspect=aspect)
        if isinstance(bigwigs, basestring):
            bigwigs = [bigwigs]
        self.plot_kwargs = {} if plot_kwargs is None else plot_kwargs
        self.bigwigs = bigwigs
        self.names = names
        self.bin_size = bin_size
        self.ylim = ylim

    def _bin_intervals(self, region, intervals):
        """
        Bin values fearch interval from Bigwig in evenly spaced bins
        suitable for plotting.
        """
        bin_coords = np.r_[slice(region.start, region.end, self.bin_size), region.end]
        bin_regions = [GenomicRegion(chromosome=region.chromosome, start=s, end=e) for s, e in it.izip(bin_coords[:-1], bin_coords[1:])]
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
                weighted_value += np.sum((interval_records["end"][s + 1:e] - interval_records["start"][s + 1:e])*
                                          interval_records["value"][s + 1:e])
            out_values[i] = weighted_value/total_range
        return bin_regions, out_values

    def _plot(self, region):
        for i, b in enumerate(self.bigwigs):
            try:
                bw = pyBigWig.open(b)
                intervals = bw.intervals(region.chromosome, region.start, region.end)
            finally:
                bw.close()
            if self.bin_size:
                regions, bw_values = self._bin_intervals(region, intervals)
            else:
                regions = [GenomicRegion(chromosome=region.chromosome, start=s, end=e) for s, e, v in intervals]
                bw_values = [v for s, e, v in intervals]
            x, y = self.get_plot_values(bw_values, regions)
            self.ax.plot(x, y, label=self.names[i] if self.names else "", **self.plot_kwargs)
        if self.names:
            self.add_legend()
        self.remove_colorbar_ax()
        sns.despine(ax=self.ax, top=True, right=True)
        if self.ylim:
            self.ax.set_ylim(self.ylim)

    def _refresh(self, region):
        pass


class GenomicTrackPlot(ScalarDataPlot):
    def __init__(self, tracks, style="step", attributes=None, title='', aspect=.2):
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
        ScalarDataPlot.__init__(self, style=style, title=title, aspect=aspect)
        if not isinstance(tracks, list):
            tracks = [tracks]
        self.tracks = tracks
        self.attributes = attributes

    def _plot(self, region=None, ax=None):
        for track in self.tracks:
            bins = track.region_bins(region)
            values = track[bins]
            regions = track.regions[bins]
            for k, v in values.iteritems():
                if not self.attributes or any(re.match(a.replace("*", ".*"), k) for a in self.attributes):
                    x, y = self.get_plot_values(v, regions)
                    self.ax.plot(x, y,
                                 label="{}{}".format(track.title + "_"
                                                     if track.title and len(self.tracks) > 1
                                                     else "", k))
        self.add_legend()
        self.remove_colorbar_ax()

    def _refresh(self, **kwargs):
        pass


class GenomicMatrixPlot(BasePlotter1D, BasePlotterMatrix):
    def __init__(self, track, attribute, y_coords=None, plot_kwargs=None, title='',
                 colormap='viridis', norm="lin", vmin=None, vmax=None,
                 show_colorbar=True, blend_zero=False,
                 unmappable_color=".9", illegal_color=None, aspect=.3):
        """
        Plot matrix from a class:`~GenomicTrack` objects

        :param track: class:`~GenomicTrack` containing the matrix
        :param attribute: Which matrix from the track object to draw
        :param y_coords: Matrices in the class:`~GenomicTrack` object are
                         unitless. Can provide the coordinates for the
                         y-direction here. Matrix has shape (X, Y) must
                         have shape Y or Y + 1
        :param plot_kwargs: Keyword-arguments passed on to pcolormesh
        :param title: Used as title for plot
        :param aspect: Default aspect ratio of the plot. Can be overriden by setting
                       the height_ratios in class:`~GenomicFigure`
        """
        BasePlotter1D.__init__(self, title=title, aspect=aspect)
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

    def _plot(self, region=None, ax=None):
        bins = self.track.region_bins(region)
        values = self.track[bins][self.attribute]
        regions = self.track.regions()[bins]
        bin_coords = np.r_[[(x.start - 1) for x in regions], (regions[-1].end)]
        X, Y = np.meshgrid(bin_coords, (self.y_coords if self.y_coords is not None else np.arange(values.shape[1] + 1)))
        color_matrix = self.get_color_matrix(values.T)
        color_tuple = color_matrix.transpose((1,0,2)).reshape(
            (color_matrix.shape[0]*color_matrix.shape[1],color_matrix.shape[2]))
        self.collection = self.ax.pcolormesh(X, Y, values.T, rasterized=True, cmap=self.colormap,
                                             norm=self.norm, **self.plot_kwargs)
        self.collection._A = None
        self.collection.set_color(color_tuple)
        if self.show_colorbar:
            self.add_colorbar()

    def _refresh(self, **kwargs):
        pass


class GeneModelPlot(BasePlotter1D):
    def __init__(self, gtf, title="", feature_types=None, id_field="gene_symbol",
                 label_format=None, label_cast=None, aspect=.2):
        import pybedtools as pbt
        BasePlotter1D.__init__(self, title=title, aspect=aspect)
        self.gtf = pbt.BedTool(gtf)
        if feature_types is None:
            feature_types = set(f[2] for f in self.gtf)
        elif isinstance(feature_types, (str, unicode)):
            feature_types = [feature_types]
        self.feature_types = feature_types
        self.id_field = id_field
        self.label_format = label_format
        self.label_cast = label_cast

    def _plot(self, region=None, ax=None):
        interval = region_to_pbt_interval(region)
        genes = self.gtf.all_hits(interval)
        trans = self.ax.get_xaxis_transform()
        pos = {k: i/(1 + len(self.feature_types)) for i, k in enumerate(self.feature_types)}
        stroke_length = max(0.03, 1/len(self.feature_types) - .05)
        for k, p in pos.iteritems():
            self.ax.text(0, p, k, transform=self.ax.transAxes, ha="left", size=5)
        for g in genes:
            if not g[2]in self.feature_types:
                continue
            if isinstance(self.id_field, int):
                label = g[self.id_field]
            else:
                try:
                    label = g.attrs[self.id_field]
                except KeyError:
                    label = ""
            gene_patch = patches.Rectangle(
                (g.start, pos[g[2]]),
                width=abs(g.end - g.start), height=stroke_length,
                transform=trans, color="black"
            )
            self.ax.add_patch(gene_patch)
            self.ax.text((g.start + g.end)/2, pos[g[2]] + stroke_length + .05,
                         self.label_format.format(self.label_cast(label)
                                                  if self.label_cast
                                                  else label)
                         if self.label_format else label,
                         transform=trans, ha="center", size="small")
            self.ax.spines['right'].set_visible(False)
            self.ax.spines['top'].set_visible(False)
            self.ax.spines['left'].set_visible(False)
            self.ax.spines['bottom'].set_visible(False)
            self.ax.xaxis.set_ticks_position('bottom')
            self.ax.yaxis.set_visible(False)
        self.remove_colorbar_ax()

    def _refresh(self, **kwargs):
        pass


class GenomicFeaturePlot(BasePlotter1D):
    def __init__(self, regions, labels, title="", color='black', aspect=.2):
        BasePlotter1D.__init__(self, title=title, aspect=aspect)
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

    def _refresh(self, **kwargs):
        pass
