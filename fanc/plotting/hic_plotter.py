import fanc
from fanc.config import config
from fanc.plotting.base_plotter import BasePlotterMatrix, BasePlotter1D, BasePlotter2D, ScalarDataPlot, \
                                       PlotMeta, BaseOverlayPlotter
from fanc.plotting.helpers import append_axes, style_ticks_whitegrid
from fanc.tools.general import str_to_int
from genomic_regions import GenomicRegion, as_region
from ..matrix import RegionMatrixTable, RegionMatrixContainer
from ..architecture.comparisons import SplitMatrix
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons
import matplotlib.gridspec as grd
import matplotlib.patches as patches
from abc import ABCMeta
import numpy as np
from scipy.ndimage.filters import gaussian_filter
import itertools as it
import types
import seaborn as sns
from future.utils import with_metaclass, string_types
from collections import defaultdict
from itertools import cycle
from ..matrix import RegionMatrix
from ..peaks import ObservedPeakFilter, FdrPeakFilter, EnrichmentPeakFilter, MappabilityPeakFilter
import logging
logger = logging.getLogger(__name__)


def prepare_hic_buffer(hic_data, buffering_strategy="relative", buffering_arg=1,
                       weight_field=None, default_value=None, smooth_sigma=None,
                       norm=True, oe=False, log=False):
    """
    Prepare :class:`~BufferedMatrix` from hic data.

    :param hic_data: :class:`~fanc.data.genomic.RegionMatrixTable` or
                     :class:`~fanc.data.genomic.RegionMatrix`
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
    if isinstance(hic_data, RegionMatrixContainer):
        return BufferedMatrix(hic_data, buffering_strategy=buffering_strategy,
                              buffering_arg=buffering_arg, weight_field=weight_field,
                              default_value=default_value, smooth_sigma=smooth_sigma,
                              norm=norm, oe=oe, log=log)
    else:
        raise ValueError("Unknown type for hic_data")


class BufferedMatrix(object):
    """
    Buffer contents of any :class:`~fanc.Hic` like objects. Matrix is
    prefetched and stored in memory. Buffer contents can quickly be fetched
    from memory. Different buffering strategies allow buffering of nearby
    regions so that adjacent parts of the matrix can quickly be fetched.
    """
    _STRATEGY_ALL = "all"
    _STRATEGY_FIXED = "fixed"
    _STRATEGY_RELATIVE = "relative"

    def __init__(self, data, buffering_strategy="relative", buffering_arg=1,
                 weight_field=None, default_value=None, smooth_sigma=None,
                 norm=True, oe=False, log=False):
        """
        Initialize a buffer for Matrix-like objects that support
        indexing using class:`~GenomicRegion` objects, such as class:`~fanc.Hic`
        or class:`~fanc.RegionMatrix` objects.

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
            raise ValueError("Only support the buffering strategies {}".format(list(self._BUFFERING_STRATEGIES.keys())))
        self.buffering_strategy = buffering_strategy
        self.buffering_arg = buffering_arg
        self.buffered_region = None
        self.buffered_matrix = None
        self.weight_field = weight_field
        self.default_value = default_value
        self.smooth_sigma = smooth_sigma
        self.norm = norm
        self.oe = oe
        self.log = log

    @classmethod
    def from_hic_matrix(cls, hic_matrix, weight_field=None, default_value=None,
                        smooth_sigma=None, norm=True, oe=False, log=False):
        """
        Wrap a :class:`~HicMatrix` in a :class:`~BufferedMatrix` container.
        Default buffering strategy is set to "all" by default.

        :param hic_matrix: :class:`~HicMatrix`
        :return: :class:`~BufferedMatrix`
        """
        bm = cls(data=None, buffering_strategy="all", weight_field=weight_field,
                 default_value=default_value, smooth_sigma=smooth_sigma,
                 norm=norm, oe=oe, log=log)
        bm.buffered_region = bm._STRATEGY_ALL
        bm.buffered_matrix = hic_matrix
        return bm

    def is_buffered_region(self, *regions):
        """
        Check if set of :class:`~GenomicRegion`s is already buffered in this matrix.

        :param regions: :class:`~GenomicRegion` object(s)
        :return:
        """
        if (self.buffered_region is None or self.buffered_matrix is None or
                (not self.buffered_region == self._STRATEGY_ALL and not
                 all(rb.contains(rq) for rb, rq in zip(self.buffered_region, regions)))):
            return False
        return True

    def get_matrix(self, *regions):
        """
        Retrieve a sub-matrix by the given :class:`~GenomicRegion` object(s).

        Will automatically load data if a non-buffered region is requested.

        :param regions: :class:`~GenomicRegion` object(s)
        :return: :class:`~HicMatrix`
        """
        # regions = tuple(reversed([r for r in regions]))
        if not self.is_buffered_region(*regions):
            logger.debug("Buffering matrix")
            self._BUFFERING_STRATEGIES[self.buffering_strategy](self, *regions)
        m = self.buffered_matrix[tuple(regions)]
        if self.smooth_sigma is not None:
            mf = gaussian_filter(m, self.smooth_sigma)
            m = fanc.matrix.RegionMatrix(mf, row_regions=m.row_regions, col_regions=m.col_regions)
        return m

    def _buffer_all(self, *regions):
        """
        No buffering, just loads everything in the object into memory.

        Obviously very memory intensive.

        :param regions: :class:`~GenomicRegion` objects
        :return: :class:`~HicMatrix`
        """
        self.buffered_region = self._STRATEGY_ALL
        self.buffered_matrix = self.data.matrix(key=tuple([slice(0, None, None)]*len(regions)),
                                                score_field=self.weight_field,
                                                default_value=self.default_value,
                                                norm=self.norm, oe=self.oe, log=self.log)

    def _buffer_relative(self, *regions):
        """
        Load the requested :class:`~GenomicRegion` and buffer an additional fration
        of the matrix given by buffering_arg*len(region)

        :param regions: :class:`~GenomicRegion` objects
        :return: :class:`~HicMatrix`
        """
        self.buffered_region = []
        for rq in regions:
            if rq.start is not None and rq.end is not None:
                rq_size = rq.end - rq.start
                new_start = max(1, rq.start - rq_size*self.buffering_arg)
                new_end = rq.end + rq_size*self.buffering_arg
                self.buffered_region.append(GenomicRegion(start=new_start, end=new_end,
                                                          chromosome=rq.chromosome))
            else:
                self.buffered_region.append(GenomicRegion(start=None, end=None, chromosome=rq.chromosome))
        self.buffered_matrix = self.data.matrix(tuple(self.buffered_region),
                                                score_field=self.weight_field,
                                                default_value=self.default_value,
                                                norm=self.norm, oe=self.oe, log=self.log)

    def _buffer_fixed(self, *regions):
        """
        Load the requested :class:`~GenomicRegion` and buffer an additional
        fixed part of the matrix given by buffering_arg

        :param regions: :class:`~GenomicRegion` objects
        :return: :class:`~HicMatrix`
        """
        self.buffered_region = []
        for rq in regions:
            if rq.start is not None and rq.end is not None:
                new_start = max(1, rq.start - self.buffering_arg)
                new_end = rq.end + self.buffering_arg
                self.buffered_region.append(GenomicRegion(start=new_start, end=new_end,
                                                          chromosome=rq.chromosome))
            else:
                self.buffered_region.append(GenomicRegion(start=None, end=None, chromosome=rq.chromosome))
        self.buffered_matrix = self.data.matrix(tuple(self.buffered_region),
                                                score_field=self.weight_field,
                                                default_value=self.default_value,
                                                norm=self.norm, oe=self.oe, log=self.log)

    @property
    def buffered_min(self):
        """
        Find the smallest non-zero buffered matrix value.

        :return: float or None if nothing is buffered
        """
        return float(np.nanmin(self.buffered_matrix[np.ma.nonzero(self.buffered_matrix)]))\
            if self.buffered_matrix is not None else None

    @property
    def buffered_max(self):
        """
        Find the largest buffered matrix value
        :return: float or None if nothing is buffered
        """
        return float(np.nanmax(self.buffered_matrix)) if self.buffered_matrix is not None else None

    _BUFFERING_STRATEGIES = {_STRATEGY_ALL: _buffer_all,
                             _STRATEGY_RELATIVE: _buffer_relative,
                             _STRATEGY_FIXED: _buffer_fixed}


class BufferedCombinedMatrix(BufferedMatrix):
    """
    A buffered square matrix where values above and below the diagonal
    come from different matrices.
    """
    def __init__(self, top_matrix, bottom_matrix, scale_matrices=True,
                 **kwargs):
        super(BufferedCombinedMatrix, self).__init__(None, **kwargs)

        if top_matrix.bin_size != bottom_matrix.bin_size:
            raise ValueError("Can only combine matrices with equal binning!")

        scaling_factor = top_matrix.scaling_factor(bottom_matrix) if scale_matrices else 1.
        self.data = SplitMatrix(top_matrix, bottom_matrix, scaling_factor=scaling_factor)


class BasePlotterHic(BasePlotterMatrix):
    """
    Base class for plotting Hi-C data.

    **warning: should always be inherited from first, before the other
    base classed**

    Makes use of matrix buffering by :class:`~BufferedMatrix` internally.
    """

    def __init__(self, hic_data, adjust_range=False, buffering_strategy="relative",
                 buffering_arg=1, weight_field=None, default_value=None, smooth_sigma=None,
                 matrix_norm=True, oe=False, log=False, **kwargs):
        """
        :param hic_data: Path to Hi-C data on disk or
                        :class:`~fanc.data.genomic.Hic` or :class:`~fanc.data.genomic.RegionMatrix`
        :param adjust_range: Draw a slider to adjust vmin/vmax interactively. Default: False
        :param buffering_strategy: A valid buffering strategy for :class:`~BufferedMatrix`
        :param buffering_arg: Adjust range of buffering for :class:`~BufferedMatrix`
        """
        super(BasePlotterHic, self).__init__(**kwargs)
        if isinstance(hic_data, string_types):
            hic_data = fanc.load(hic_data, mode="r")
        self.hic_data = hic_data
        self.hic_buffer = prepare_hic_buffer(hic_data, buffering_strategy=buffering_strategy,
                                             buffering_arg=buffering_arg, weight_field=weight_field,
                                             default_value=default_value, smooth_sigma=smooth_sigma,
                                             norm=matrix_norm, oe=oe, log=log)
        self.slider = None
        self.adjust_range = adjust_range
        self.vmax_slider = None


class SquareMatrixPlot(BasePlotterHic, BasePlotter2D):
    """
    Plot Hi-C map as a square matrix.
    """
    sharey = True

    def __init__(self, hic_data, flip=False, **kwargs):
        """
        :param flip: Transpose matrix before plotting
        """
        super(SquareMatrixPlot, self).__init__(hic_data=hic_data, **kwargs)
        self.vmax_slider = None
        self.current_matrix = None
        self.flip = flip
        self.im = None

    def _plot(self, region):
        self._refresh(region)

        self.ax.spines['right'].set_visible(False)
        self.ax.spines['top'].set_visible(False)
        self.ax.xaxis.set_ticks_position('bottom')
        self.ax.yaxis.set_ticks_position('left')

        if self.adjust_range:
            self.add_adj_slider()

        return self.ax

    def add_adj_slider(self, ax=None):
        if ax is None:
            pad = self._total_padding if self._total_padding is not None else 0.
            ax = append_axes(self.ax, 'bottom', config.adjustment_slider_height, pad)

        vmin = self.hic_buffer.buffered_min
        vmax = self.hic_buffer.buffered_max
        self.vmax_slider = Slider(ax, 'vmax', vmin,
                                  vmax, valinit=self.vmax,
                                  facecolor='#dddddd', edgecolor='none')

        self.vmax_slider.on_changed(self._slider_refresh)

        self.ax.patch.set_visible(False)

    def _slider_refresh(self, val):
        # new_vmin = self.vmin_slider.val
        new_vmax = self.vmax_slider.val
        if self.colorbar_symmetry is not None:
            diff = abs(self.colorbar_symmetry - new_vmax)
            self.im.set_clim(vmin=self.colorbar_symmetry - diff, vmax=self.colorbar_symmetry + diff)
        else:
            self.im.set_clim(vmin=self.hic_buffer.buffered_min, vmax=new_vmax)
        # Hack to force redraw of image data
        self._refresh(None)

        if self.colorbar is not None:
            self.update_colorbar(vmax=new_vmax)

    def _refresh(self, region):
        if region is not None:
            self.current_matrix = self.hic_buffer.get_matrix(*region)
        old_image = self.im

        m = np.transpose(self.current_matrix)
        if self.flip:
            m = np.flipud(m)
            extent = [self.current_matrix.row_regions[0].start, self.current_matrix.row_regions[-1].end,
                      self.current_matrix.col_regions[0].start, self.current_matrix.col_regions[-1].end]
        else:
            extent = [self.current_matrix.row_regions[0].start, self.current_matrix.row_regions[-1].end,
                      self.current_matrix.col_regions[-1].end, self.current_matrix.col_regions[0].start]

        color_matrix = self.get_color_matrix(m)
        self.im = self.ax.imshow(color_matrix, interpolation='none', aspect='auto',
                                 cmap=self.colormap, norm=self._map_norm, extent=extent)

        self.ax.set_xlim(extent[0], extent[1])
        self.ax.set_ylim(extent[2], extent[3])

        if old_image is not None:
            old_image.remove()


HicPlot2D = SquareMatrixPlot


class SplitMatrixPlot(HicPlot2D):
    def __init__(self, hic_top, hic_bottom, scale_matrices=True,
                 weight_field=None, default_value=None, smooth_sigma=None,
                 matrix_norm=True, oe=False, log=False, **kwargs):
        super(SplitMatrixPlot, self).__init__(hic_top, **kwargs)
        self.hic_top = hic_top
        self.hic_bottom = hic_bottom
        self.hic_buffer = BufferedCombinedMatrix(hic_bottom, hic_top, scale_matrices=scale_matrices,
                                                 weight_field=weight_field, default_value=default_value,
                                                 smooth_sigma=smooth_sigma, norm=matrix_norm, oe=oe,
                                                 log=log)


HicComparisonPlot2D = SplitMatrixPlot


class HicSlicePlot(ScalarDataPlot):
    """
    Draw Hi-C data as virtual 4C-plot. All interactions that
    involve the slice region are shown.
    """

    def __init__(self, hic_data, slice_region, names=None,
                 colors=None, fill=None,
                 buffering_strategy="relative", buffering_arg=1,
                 weight_field=None, default_value=None, **kwargs):
        """
        :param hic_data: :class:`~fanc.Hic` or :class:`~fanc.RegionMatrix`. Can be list of
                         multiple Hi-C datasets.
        :param slice_region: String ("2L:1000000-1500000") or :class:`~GenomicRegion`.
                             All interactions involving this region are shown.
        :param names: If multiple Hi-C datasets are provided, can pass a list of names.
                      Are used as names in the legend of the plot.
        :param buffering_strategy: A valid buffering strategy for class:`~BufferedMatrix`
        :param buffering_arg: Adjust range of buffering for class:`~BufferedMatrix`
        """
        kwargs.setdefault("aspect", .3)
        super(HicSlicePlot, self).__init__(**kwargs)
        if not isinstance(hic_data, (list, tuple)):
            hic_data = [hic_data]
        self.hic_buffers = []
        for h in hic_data:
            hb = prepare_hic_buffer(h,
                                    buffering_strategy=buffering_strategy,
                                    buffering_arg=buffering_arg,
                                    weight_field=weight_field,
                                    default_value=default_value)
            self.hic_buffers.append(hb)
        self.names = names
        if isinstance(slice_region, string_types):
            slice_region = GenomicRegion.from_string(slice_region)
        self.slice_region = slice_region
        self.x = None
        self.y = None
        self.lines = []
        self.fill = fill

        if colors is None:
            prop_cycle = plt.rcParams['axes.prop_cycle']
            colors = prop_cycle.by_key()['color']
        elif isinstance(colors, string_types):
            colors = [colors]

        self.colors = colors

    def _refresh(self, region):
        for line in self.lines:
            line.remove()
        self.lines = []

        color_cycle = cycle(self.colors)
        for i, b in enumerate(self.hic_buffers):
            hm = b.get_matrix(self.slice_region, region.expand(relative=0.5))
            m = np.mean(hm, axis=0)
            bin_coords = np.r_[[x.start for x in hm.col_regions], hm.col_regions[-1].end]
            bin_coords = (bin_coords[1:] + bin_coords[:-1])/2
            line = self.ax.plot(bin_coords, m,
                                color=next(color_cycle),
                                label=self.names[i] if self.names else "")[0]
            self.lines.append(line)

            if self.fill:
                self.ax.fill_between(bin_coords, [0] * len(m), m, color=line.get_color())

    def _plot(self, region):
        self._refresh(region)

        # plot patch at slice region
        bottom, top = self.ax.get_ylim()
        patch_region = self.slice_region.copy()

        rect = patches.Rectangle((patch_region.start, bottom),
                                 len(patch_region), (top - bottom)/20,
                                 facecolor='grey', edgecolor='grey')
        self.ax.add_patch(rect)

        if self.names:
            self.add_legend()
        self.remove_colorbar_ax()
        sns.despine(ax=self.ax, top=True, right=True)


class TriangularMatrixPlot(BasePlotterHic, BasePlotter1D):
    """
    A triangle Hi-C heatmap plot.
    """

    def __init__(self, hic_data, max_dist=None, proportional=True,
                 rasterized=True, **kwargs):
        """
        :param max_dist: Only draw interactions up to this distance
        :param proportional: Automatically determine aspect ratio of plot
                             so that x- and y-axis are proportional. Default: True
        :param rasterized: Draw map as image (True) or vector graphic (False).
                           Default: True
        """
        kwargs.setdefault("aspect", None)
        super(TriangularMatrixPlot, self).__init__(hic_data=hic_data, **kwargs)
        self.proportional = proportional
        self.max_dist = str_to_int(max_dist)
        self.hm = None
        self.rasterized = rasterized
        self.collection = None

    def _plot(self, region):
        logger.debug("Generating matrix from hic object")
        if region is None:
            raise ValueError("Cannot plot triangle plot for whole genome.")
        if region.start is None:
            region.start = 1
        if region.end is None:
            region.end = self.hic_data.chromosome_lengths[region.chromosome]
        if self.aspect is None and self.proportional:
            if self.max_dist is None:
                self.aspect = .5
            else:
                rl = region.end - region.start
                self.aspect = .5*min(self.max_dist, rl)/rl
            self._dimensions_stale = True

        # Have to copy unfortunately, otherwise modifying matrix in buffer
        self._refresh(region)

        # set limits and aspect ratio
        # self.ax.set_aspect(aspect="equal")
        self.ax.set_ylim(0, self.max_dist/2 if self.max_dist else (region.end-region.start)/2)
        # remove outline everywhere except at bottom
        sns.despine(ax=self.ax, top=True, right=True, left=True)
        self.ax.set_yticks([])
        # hide background patch
        self.ax.patch.set_visible(False)
        if self.adjust_range:
            self.add_adj_slider()

        def drag_pan(self, button, key, x, y):
            mpl.axes.Axes.drag_pan(self, button, 'x', x, y)  # pretend key=='x'

        self.ax.drag_pan = types.MethodType(drag_pan, self.ax)

    def _mesh_data(self, region):
        hm = self.hic_buffer.get_matrix(region, region)
        hm_copy = fanc.matrix.RegionMatrix(np.copy(hm), col_regions=hm.col_regions,
                                           row_regions=hm.row_regions)
        # update coordinates
        bin_coords = np.r_[[x.start for x in hm_copy.row_regions], hm_copy.row_regions[-1].end]
        # Make sure the matrix is not protruding over the end of the requested plotting region
        if bin_coords[0] < region.start <= bin_coords[1]:
            bin_coords[0] = region.start
        if bin_coords[-1] > region.end >= bin_coords[-2]:
            bin_coords[-1] = region.end
        bin_coords = np.true_divide(bin_coords, np.sqrt(2))
        x, y = np.meshgrid(bin_coords, bin_coords)
        # rotatate coordinate matrix 45 degrees
        sin45 = np.sin(np.radians(45))
        x_, y_ = x * sin45 + y * sin45, x * sin45 - y * sin45

        return x_, y_, hm_copy

    def _update_mesh_colors(self):
        # pcolormesh doesn't support plotting RGB arrays directly like imshow, have to workaround
        # See https://github.com/matplotlib/matplotlib/issues/4277
        # http://stackoverflow.com/questions/29232439/plotting-an-irregularly-spaced-rgb-image-in-python/29232668?noredirect=1#comment46710586_29232668
        color_matrix = self.get_color_matrix(self.hm)
        color_tuple = color_matrix.transpose((1, 0, 2)).reshape(
            (color_matrix.shape[0] * color_matrix.shape[1], color_matrix.shape[2]))
        self.collection.set_facecolor(color_tuple)

    def _refresh(self, region):
        x_, y_, hm = self._mesh_data(region)
        self.hm = hm

        old_collection = self.collection
        self.collection = self.ax.pcolormesh(x_, y_, hm, cmap=self.colormap, norm=self._map_norm,
                                             rasterized=self.rasterized)
        self.collection._A = None
        self._update_mesh_colors()

        if old_collection is not None:
            old_collection.remove()

    def add_adj_slider(self, ax=None):
        if ax is None:
            pad = self._total_padding if self._total_padding is not None else 0.
            ax = append_axes(self.ax, 'bottom', config.adjustment_slider_height, pad)

        self.vmax_slider = Slider(ax, 'vmax', self.hic_buffer.buffered_min,
                                  self.hic_buffer.buffered_max, valinit=self.vmax,
                                  facecolor='#dddddd', edgecolor='none')

        self.vmax_slider.on_changed(self._slider_refresh)

    def _slider_refresh(self, val):
        # new_vmin = self.vmin_slider.val
        new_vmax = self.vmax_slider.val
        if self.colorbar_symmetry is not None:
            diff = abs(self.colorbar_symmetry-new_vmax)
            self._update_norm(vmin=self.colorbar_symmetry-diff, vmax=self.colorbar_symmetry+diff)
        else:
            self._update_norm(vmax=new_vmax)
        self._update_mesh_colors()
        if self.colorbar is not None:
            self.update_colorbar(vmax=new_vmax)


HicPlot = TriangularMatrixPlot


class HicPeakPlot(BaseOverlayPlotter):
    """
    Overlay peaks onto Hicplot or HicPlot2D. Accepts
    :class:`~fanc.data.network.PeakInfo`.
    Add to HicPlot or HicPlot2D using add_overlay method.
    """
    def __init__(self, peaks, radius=None, circle_props={}, **kwargs):
        """
        :param peaks: fanc peaks instance
        :param radius: Radius in bp for plotted circles.
                       If not specified (default), use the radius of the
                       peak itself. This is often too small to see,
                       providing a value like 50000 helps. Default: None
        :param circe_props: Dictionary with properties for the plotted circles
                            for the matplotlib.patches.Circle constructor.
                            Default: Black edges, no fill, linewidth 3 pt
        """
        super(HicPeakPlot, self).__init__(**kwargs)
        self.peaks = peaks
        self.radius = radius
        self.circle_props = {
            "edgecolor": "black",
            "fill": False,
            "linewidth": 3,
        }
        self.circle_props.update(circle_props)

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
        plot_func = plot_dispatch[base_plot_class]
        peaks_gen = self.peaks.edge_subset((region, region), distances_in_bp=True)
        for p in peaks_gen:
            plot_func((p.source_node.start + p.source_node.end)/2,
                      (p.sink_node.start + p.sink_node.end)/2,
                      self.radius if self.radius is not None else getattr(p, "radius", self.peaks.bin_size*3)*self.peaks.bin_size)


class EdgeFilterBuffer(object):
    def __init__(self, hic, plot_field='weight', default_value=0, log2=False):
        self.hic = hic
        self._region_buffer = dict()
        self._edge_buffer = dict()
        self._filters = []
        self._plot_field = plot_field
        self._default_value = default_value
        self._log2 = log2
        self._last_matrix = None

    def add_filter(self, filter):
        self._filters.append(filter)

    @property
    def buffered_min(self):
        """
        Find the smallest non-zero buffered matrix value.

        :return: float or None if nothing is buffered
        """
        return float(np.nanmin(self._last_matrix[np.ma.nonzero(self._last_matrix)])) \
            if self._last_matrix is not None else None

    @property
    def buffered_max(self):
        """
        Find the largest buffered matrix value
        :return: float or None if nothing is buffered
        """
        return float(np.nanmax(self._last_matrix)) if self._last_matrix is not None else None

    def get_matrix(self, *regions):
        row_region = regions[0]
        col_region = regions[1]
        row_key = (row_region.chromosome, row_region.start, row_region.end)
        col_key = (col_region.chromosome, col_region.start, col_region.end)

        if row_key in self._region_buffer:
            sub_row_regions = self._region_buffer[row_key]
        else:
            sub_row_regions = list(self.hic.regions(row_region))
            self._region_buffer[row_key] = sub_row_regions

        if col_key in self._region_buffer:
            sub_col_regions = self._region_buffer[col_key]
        else:
            sub_col_regions = list(self.hic.regions(col_region))
            self._region_buffer[col_key] = sub_col_regions

        row_offset = sub_row_regions[0].ix
        col_offset = sub_col_regions[0].ix

        if (row_key, col_key) in self._edge_buffer:
            edges = self._edge_buffer[(row_key, col_key)]
        else:
            edges = list(self.hic.edge_subset((row_region, col_region)))
            self._edge_buffer[(row_key, col_key)] = edges

        m = np.full((len(sub_row_regions), len(sub_col_regions)), self._default_value,
                    dtype=np.float64)

        n_filtered = 0
        n_filtered_by_filter = defaultdict(int)
        for edge in edges:
            valid = True
            for filter_ix, filter in enumerate(self._filters):
                if not filter.valid_peak(edge):
                    n_filtered_by_filter[filter_ix] += 1
                    valid = False
                    break

            if not valid:
                n_filtered += 1
                continue

            source = edge.source
            sink = edge.sink
            weight = getattr(edge, self._plot_field)
            if self._log2:
                weight = np.log2(weight)

            ir = source - row_offset
            jr = sink - col_offset
            if 0 <= ir < m.shape[0] and 0 <= jr < m.shape[1]:
                m[ir, jr] = weight

            ir = sink - row_offset
            jr = source - col_offset
            if 0 <= ir < m.shape[0] and 0 <= jr < m.shape[1]:
                m[ir, jr] = weight

        rm = RegionMatrix(m, row_regions=sub_row_regions, col_regions=sub_col_regions)
        self._last_matrix = rm
        return rm


class EdgeHicPlot(HicPlot2D):
    def __init__(self, hic, plot_field='weight', default_value=0,
                 log2=False, highlight_edges=False,
                 highlight_limit=100, *args, **kwargs):
        HicPlot2D.__init__(self, hic, *args, **kwargs)

        self.hic_buffer = EdgeFilterBuffer(hic, plot_field=plot_field,
                                           default_value=default_value,
                                           log2=log2)
        self._highlight_edges = highlight_edges
        self._highlight_circles = []
        self._highlight_limit = highlight_limit

    def update_highlights(self, m=None):
        for circle in self._highlight_circles:
            circle.remove()

        self._highlight_circles = []

        if self._highlight_edges and m is not None:
            all_x, all_y = np.where(m > 0)
            if not len(all_x) > self._highlight_limit:
                for x, y in zip(all_x, all_y):
                    region_x = m.row_regions[x]
                    region_y = m.col_regions[y]

                    circle = patches.Circle((region_x.center, region_y.center), len(region_x)*3,
                                            facecolor='red')
                    self.ax.add_patch(circle)
                    self._highlight_circles.append(circle)

    def _plot(self, region):
        HicPlot2D._plot(self, region)

    def _refresh(self, region):
        HicPlot2D._refresh(self, region)


class PeakParameterPlot(object):
    def __init__(self, peaks, font_size=5, observed_init=1.,
                 oe_d_init=1., oe_h_init=1., oe_v_init=1., oe_l_init=1.,
                 fdr_d_init=.1, fdr_h_init=.1, fdr_v_init=.1, fdr_l_init=.1,
                 mappability_d_init=0., mappability_h_init=0.,
                 mappability_v_init=0., mappability_l_init=0.,
                 oe_slider_range=(0, 5), oe_slider_step=0.1,
                 fdr_slider_range=(0, .1), fdr_slider_step=0.001,
                 mappability_slider_range=(0, 1.), mappability_slider_step=0.05,
                 observed_range=(0, 50), observed_step=1,
                 **kwargs):
        self.peaks = peaks
        self.font_size = font_size
        self.hic_args = kwargs
        
        self.observed_init = observed_init

        self.oe_init = {
            'd': oe_d_init,
            'v': oe_v_init,
            'h': oe_h_init,
            'l': oe_l_init,
        }

        self.fdr_init = {
            'd': fdr_d_init,
            'v': fdr_v_init,
            'h': fdr_h_init,
            'l': fdr_l_init,
        }

        self.mappability_init = {
            'd': mappability_d_init,
            'v': mappability_v_init,
            'h': mappability_h_init,
            'l': mappability_l_init,
        }

        self.fdr_range = fdr_slider_range
        self.fdr_step = fdr_slider_step
        self.oe_range = oe_slider_range
        self.oe_step = oe_slider_step
        self.mappability_range = mappability_slider_range
        self.mappability_step = mappability_slider_step
        self.observed_range = observed_range
        self.observed_step = observed_step

        self.observed_cutoff = self.observed_init
        self.oe_cutoffs = self.oe_init.copy()
        self.fdr_cutoffs = self.fdr_init.copy()
        self.mappability_cutoffs = self.mappability_init.copy()

        self.ax_fdr_sliders = {}
        self.fdr_sliders = {}
        self.oe_sliders = {}
        self.ax_oe_sliders = {}
        self.ax_mappability_sliders = {}
        self.mappability_sliders = {}
        self.observed_slider = None

        self.observed_filter = None
        self.fdr_filter = None
        self.mappability_filter = None
        self.oe_filter = None

        self.filtered_plots = []
        self.hic_plots = []

        self.button = None
        self.region_pairs = []

        self.fig = None

    def plot(self, *regions):
        plt.rcParams.update({'font.size': self.font_size})

        # process region pairs
        self.region_pairs = []
        for pair in regions:
            if isinstance(pair, string_types):
                r = as_region(pair)
                pair = (r, r)
            elif isinstance(pair, GenomicRegion):
                pair = (pair, pair)

            try:
                r1, r2, vmax = pair
            except ValueError:
                r1, r2 = pair
                if isinstance(r2, float) or isinstance(r2, int):
                    vmax = r2
                    r2 = r1
                else:
                    vmax = None
            r1 = as_region(r1)
            r2 = as_region(r2)
            self.region_pairs.append((r1, r2, vmax))

        #
        # necessary plots: hic, oe_d, fdr_d, uncorrected,
        #
        gs = grd.GridSpec(len(self.region_pairs) + 4, 4,
                          height_ratios=[10] * len(self.region_pairs) + [1, 1, 1, 3],
                          wspace=0.3, hspace=0.5)

        self.fig = plt.figure(figsize=(10, len(self.region_pairs) * 2 + 2), dpi=150)

        # sliders
        inner_observed_gs = grd.GridSpecFromSubplotSpec(3, 1,
                                                        subplot_spec=gs[len(self.region_pairs) + 3, 0],
                                                        wspace=0.0, hspace=0.0)

        ax_observed_slider = plt.subplot(inner_observed_gs[0, 0])
        self.observed_slider = Slider(ax_observed_slider, 'uncorrected',
                                      self.observed_range[0], self.observed_range[1],
                                      valinit=self.observed_init, valstep=self.observed_step)

        self.ax_oe_sliders = dict()
        self.oe_sliders = dict()
        self.ax_fdr_sliders = dict()
        self.fdr_sliders = dict()
        self.ax_mappability_sliders = dict()
        self.mappability_sliders = dict()
        for i, neighborhood in enumerate(['d', 'h', 'v', 'l']):
            # O/E
            self.ax_oe_sliders[neighborhood] = plt.subplot(gs[len(self.region_pairs) + 0, i])
            self.oe_sliders[neighborhood] = Slider(self.ax_oe_sliders[neighborhood],
                                                   'O/E {}'.format(neighborhood.upper()),
                                                   self.oe_range[0], self.oe_range[1],
                                                   valinit=self.oe_init[neighborhood],
                                                   valstep=self.oe_step)
            # FDR
            self.ax_fdr_sliders[neighborhood] = plt.subplot(gs[len(self.region_pairs) + 1, i])
            self.fdr_sliders[neighborhood] = Slider(self.ax_fdr_sliders[neighborhood],
                                                    'FDR {}'.format(neighborhood.upper()),
                                                    self.fdr_range[0], self.fdr_range[1],
                                                    valinit=self.fdr_init[neighborhood],
                                                    valstep=self.fdr_step)
            # Mappability
            self.ax_mappability_sliders[neighborhood] = plt.subplot(gs[len(self.region_pairs) + 2, i])
            self.mappability_sliders[neighborhood] = Slider(self.ax_mappability_sliders[neighborhood],
                                                            'Map {}'.format(neighborhood.upper()),
                                                            self.mappability_range[0],
                                                            self.mappability_range[1],
                                                            valinit=self.mappability_init[neighborhood],
                                                            valstep=self.mappability_step)

        # check button
        inner_button_gs = grd.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[len(self.region_pairs) + 3, 1],
                                                      wspace=0.0, hspace=0.0)
        ax_button = plt.subplot(inner_button_gs[0, 0])
        self.button = CheckButtons(ax_button, ['Show loops'], [False])

        # filters
        self.observed_filter = ObservedPeakFilter(cutoff=self.observed_init)
        self.oe_filter = EnrichmentPeakFilter(enrichment_d_cutoff=self.oe_init['d'],
                                              enrichment_h_cutoff=self.oe_init['h'],
                                              enrichment_v_cutoff=self.oe_init['v'],
                                              enrichment_ll_cutoff=self.oe_init['l'])
        self.fdr_filter = FdrPeakFilter(fdr_d_cutoff=self.fdr_init['d'],
                                        fdr_ll_cutoff=self.fdr_init['h'],
                                        fdr_h_cutoff=self.fdr_init['v'],
                                        fdr_v_cutoff=self.fdr_init['l'])
        self.mappability_filter = MappabilityPeakFilter(mappability_d_cutoff=self.mappability_init['d'],
                                                        mappability_h_cutoff=self.mappability_init['h'],
                                                        mappability_v_cutoff=self.mappability_init['v'],
                                                        mappability_ll_cutoff=self.mappability_init['l'])

        self.hic_plots = []
        self.filtered_plots = []
        for i, (r1, r2, vmax) in enumerate(self.region_pairs):
            ax_hic = plt.subplot(gs[i, 0])
            ax_filtered = plt.subplot(gs[i, 1])
            ax_oe = plt.subplot(gs[i, 2])
            ax_fdr = plt.subplot(gs[i, 3])

            hic_plot = EdgeHicPlot(self.peaks,
                                   ax=ax_hic, show_colorbar=False, adjust_range=False,
                                   unmappable_color='white', vmax=vmax,
                                   highlight_edges=True, highlight_limit=200,
                                   draw_tick_legend=False,
                                   **self.hic_args)
            filtered_plot = EdgeHicPlot(self.peaks,
                                        ax=ax_filtered, show_colorbar=False, adjust_range=False,
                                        unmappable_color='white', vmax=vmax,
                                        draw_tick_legend=False,
                                        **self.hic_args)
            oe_plot = EdgeHicPlot(self.peaks, plot_field='oe_d', default_value=0, norm='lin',
                                  ax=ax_oe, show_colorbar=False, colormap='white_red', adjust_range=False,
                                  vmin=1, vmax=4, log2=False, unmappable_color='white',
                                  draw_tick_legend=False,)
            fdr_plot = EdgeHicPlot(self.peaks, plot_field='fdr_d', default_value=1, norm='lin',
                                   ax=ax_fdr, show_colorbar=False, adjust_range=False,
                                   vmin=0, vmax=0.05, colormap='Greys_r', unmappable_color='white',
                                   draw_tick_legend=False)

            filtered_plot.hic_buffer.add_filter(self.observed_filter)
            filtered_plot.hic_buffer.add_filter(self.oe_filter)
            filtered_plot.hic_buffer.add_filter(self.fdr_filter)
            filtered_plot.hic_buffer.add_filter(self.mappability_filter)

            self.hic_plots.append(hic_plot)
            self.filtered_plots.append(filtered_plot)

            hic_plot.plot((r1, r2))
            filtered_plot.plot((r1, r2))
            oe_plot.plot((r1, r2))
            fdr_plot.plot((r1, r2))

            ax_filtered.set_yticklabels([])
            ax_oe.set_yticklabels([])
            ax_fdr.set_yticklabels([])

        for neighborhood in ['d', 'h', 'v', 'l']:
            self.oe_sliders[neighborhood].on_changed(self.update_oe_filter)
            self.fdr_sliders[neighborhood].on_changed(self.update_fdr_filter)
            self.mappability_sliders[neighborhood].on_changed(self.update_mappability_filter)

        self.observed_slider.on_changed(self.update_observed_filter)

        self.button.on_clicked(self.refresh_plots)
        self.refresh_plots()
        return self.fig

    def refresh_plots(self, event=None):
        for i, (r1, r2, vmax) in enumerate(self.region_pairs):
            self.filtered_plots[i].refresh((r1, r2))

            if self.button.get_status()[0]:
                self.hic_plots[i].update_highlights(self.filtered_plots[i].hic_buffer._last_matrix)
            else:
                self.hic_plots[i].update_highlights()

        self.fig.canvas.draw()

    def update_observed_filter(self, event):
        self.observed_cutoff = self.observed_slider.val
        self.observed_filter.cutoff = self.observed_cutoff
        logger.info("Observed cutoff set to {}".format(self.observed_cutoff))
        self.refresh_plots()

    def update_oe_filter(self, event):
        self.oe_cutoffs = {
            'd': self.oe_sliders['d'].val,
            'h': self.oe_sliders['h'].val,
            'v': self.oe_sliders['v'].val,
            'l': self.oe_sliders['l'].val
        }

        self.oe_filter.enrichment_d_cutoff = self.oe_cutoffs['d']
        self.oe_filter.enrichment_h_cutoff = self.oe_cutoffs['h']
        self.oe_filter.enrichment_v_cutoff = self.oe_cutoffs['v']
        self.oe_filter.enrichment_ll_cutoff = self.oe_cutoffs['l']
        logger.info("O/E cutoffs set to: {}".format(self.oe_cutoffs))
        self.refresh_plots()

    def update_fdr_filter(self, event):
        self.fdr_cutoffs = {
            'd': self.fdr_sliders['d'].val,
            'h': self.fdr_sliders['h'].val,
            'v': self.fdr_sliders['v'].val,
            'l': self.fdr_sliders['l'].val
        }
        self.fdr_filter.fdr_d_cutoff = self.fdr_cutoffs['d']
        self.fdr_filter.fdr_h_cutoff = self.fdr_cutoffs['h']
        self.fdr_filter.fdr_v_cutoff = self.fdr_cutoffs['v']
        self.fdr_filter.fdr_ll_cutoff = self.fdr_cutoffs['l']

        logger.info("FDR cutoffs set to: {}".format(self.fdr_cutoffs))
        self.refresh_plots()

    def update_mappability_filter(self, event):
        self.mappability_cutoffs = {
            'd': self.mappability_sliders['d'].val,
            'h': self.mappability_sliders['h'].val,
            'v': self.mappability_sliders['v'].val,
            'l': self.mappability_sliders['l'].val
        }
        self.mappability_filter.mappability_d_cutoff = self.mappability_cutoffs['d']
        self.mappability_filter.mappability_h_cutoff = self.mappability_cutoffs['h']
        self.mappability_filter.mappability_v_cutoff = self.mappability_cutoffs['v']
        self.mappability_filter.mappability_ll_cutoff = self.mappability_cutoffs['l']

        logger.info("Mappability cutoffs set to: {}".format(self.mappability_cutoffs))
        self.refresh_plots()
