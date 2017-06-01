from __future__ import division
import matplotlib as mpl
import numpy as np
import pybedtools as pbt

style_ticks_whitegrid = {
    'axes.axisbelow': True,
    'axes.edgecolor': '.15',
    'axes.facecolor': 'white',
    'axes.grid': True,
    'axes.labelcolor': '.15',
    'axes.linewidth': 1.25,
    'figure.facecolor': 'white',
    'font.family': ['sans-serif'],
    'grid.color': '.8',
    'grid.linestyle': '-',
    'image.cmap': 'Greys',
    'legend.frameon': False,
    'legend.numpoints': 1,
    'legend.scatterpoints': 1,
    'lines.solid_capstyle': 'round',
    'text.color': '.15',
    'xtick.color': '.15',
    'xtick.direction': 'out',
    'xtick.major.size': 6,
    'xtick.minor.size': 3,
    'ytick.color': '.15',
    'ytick.direction': 'out',
    'ytick.major.size': 6,
    'ytick.minor.size': 3}


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
