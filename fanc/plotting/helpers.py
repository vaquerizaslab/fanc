from __future__ import division
import matplotlib as mpl
import numpy as np
from math import log10, floor
import pybedtools as pbt
import fanc
from genomic_regions import RegionBased


style_ticks_whitegrid = {
    'axes.axisbelow': True,
    'axes.edgecolor': '.15',
    'axes.facecolor': 'white',
    'axes.grid': True,
    'figure.facecolor': 'white',
    'grid.color': '.8',
    'grid.linestyle': '-',
    'lines.solid_capstyle': 'round',
}


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
        width = length
        height = thickness
        if side == "top":
            hor = bbox.x0 + (bbox.width - length)/2
            vert = bbox.y0 + bbox.height + padding
        else:
            hor = bbox.x0 + (bbox.width - length)/2
            vert = bbox.y0 - padding - thickness
    elif side in ("right", "left"):
        thickness = thickness/figsize[0]
        padding = padding/figsize[0]
        length = length/figsize[1] if length is not None else shrink*bbox.height
        width = thickness
        height = length
        if side == "right":
            hor = bbox.x0 + bbox.width + padding
            vert = bbox.y0 + (bbox.height - length)/2
        else:
            hor = bbox.x0 - padding - thickness
            vert = bbox.y0 + (bbox.height - length)/2
    else:
        raise ValueError("Illegal parameter side '{}'".format(side))
    return parent.figure.add_axes([hor, vert, width, height], **kwargs)


def absolute_wspace_hspace(fig, gs, wspace=None, hspace=None):
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
    tot_width = sp_params.right - sp_params.left
    tot_height = sp_params.top - sp_params.bottom
    nrows, ncols = gs.get_geometry()
    if wspace is not None:
        wspace = wspace/figsize[0]
        wspace = wspace*ncols/(tot_width - wspace*ncols + wspace)
    else:
        wspace = gs.wspace
    if hspace is not None:
        hspace = hspace/figsize[1]
        hspace = hspace*nrows/(tot_height - hspace*nrows + hspace)
    else:
        hspace = gs.hspace
    if not wspace > 0 or not hspace > 0:
        raise ValueError("Invalid relative spacing ({}, {}) calculated, "
                         "Probably distance set too large.".format(wspace, hspace))
    gs.update(wspace=wspace, hspace=hspace)
    return wspace, hspace


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


def box_coords_abs_to_rel(top, left, width, height, figsize):
    f_width, f_height = figsize
    rel_bottom = (f_height - top - height)/f_height
    rel_left = left/f_width
    rel_width = width/f_width
    rel_height = height/f_height
    return rel_left, rel_bottom, rel_width, rel_height


# Borrowed from figure.text method
# https://github.com/matplotlib/matplotlib/blob/a4999acbbf6ebd6fa211f70becd49887dce663ab/lib/matplotlib/figure.py#L1495
def figure_line(fig, xdata, ydata, **kwargs):
    """
    Add a line to the figure, independent of axes.
    Coordinates in (0, 1) relative to bottom left of the figure.
    All kwargs are passed to Line2D constructor.
    """
    l = mpl.lines.Line2D(xdata, ydata, **kwargs)
    fig._set_artist_props(l)
    fig.lines.append(l)
    l._remove_method = lambda h: fig.lines.remove(h)
    fig.stale = True
    return l


# Borrowed from figure.text method
# https://github.com/matplotlib/matplotlib/blob/a4999acbbf6ebd6fa211f70becd49887dce663ab/lib/matplotlib/figure.py#L1495
def figure_rectangle(fig, xy, width, height, **kwargs):
    """
    Add a rectangle to the given figure independent of axes.
    Coordinates in (0, 1) relative to bottom left of the figure.
    All kwargs are passed to Rectangle constructor.
    """
    p = mpl.patches.Rectangle(xy, width, height, **kwargs)
    fig._set_artist_props(p)
    fig.patches.append(p)
    p._remove_method = lambda h: fig.patches.remove(h)
    fig.stale = True
    return p


# From https://stackoverflow.com/a/3413529/4603385
def round_sig(x, sig=2):
    if x == 0:
        return 0.
    return round(x, sig - int(floor(log10(abs(x)))) - 1)


class LimitGroup(object):
    """
    Can be used for synchronizing axis limits across multiple
    plots. Pass the same instance of this class to all plots
    that should have synchronized axis limits.
    """

    def __init__(self, limit=None, sig=2):
        """
        :param limit: tuple (vmin, vmax) to set absolute limits
                      for axis. Final limits will be chosen
                      within this absolute limit. Pass None
                      for vmin or vmax to set no limit.
                      Default: (None, None)
        :param sig: Round limits to sig significant digits.
                    If None, don't round.
                    Default: 2
        """
        self.limit = limit
        if self.limit is None:
            self.limit = (None, None)
        self.sig = sig
        self.limit_list = []

    def add_limit(self, limit):
        self.limit_list.append(limit)

    def get_limit(self):
        if len(self.limit_list) < 1:
            return self.limit
        vmin = min(x[0] if x[0] is not None else float("+Inf") for x in self.limit_list if x is not None)
        vmax = max(x[1] if x[1] is not None else float("-Inf") for x in self.limit_list if x is not None)
        if vmin is None or (self.limit[0] is not None and vmin < self.limit[0]):
            vmin = self.limit[0]
        if vmax is None or (self.limit[1] is not None and vmax > self.limit[1]):
            vmax = self.limit[1]
        if self.sig is not None:
            vmin = round_sig(vmin, self.sig)
            vmax = round_sig(vmax, self.sig)
        return (vmin, vmax)

    def reset_limit(self):
        self.limit_list = []


def parse_bedtool_input(x):
    """
    Pass x to BedTool constructor and return.
    If x is already BedTool return withouth modification.
    """
    if isinstance(x, pbt.BedTool):
        return x
    return pbt.BedTool(x)


def get_region_based_object(input_object):
    if isinstance(input_object, RegionBased):
        return input_object
    return fanc.load(input_object)


def load_score_data(data):
    # If it's already an instance of fanc data, just return it
    if isinstance(data, RegionBased):
        return data

    # If it's anything else let's hope fanc.load can deal with it
    return fanc.load(data)

