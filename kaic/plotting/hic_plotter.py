import kaic
from kaic.config import config
from kaic.plotting.base_plotter import BasePlotterMatrix, BasePlotter1D, BasePlotter2D, ScalarDataPlot, \
                                       PlotMeta, BaseOverlayPlotter
from kaic.plotting.helpers import append_axes, style_ticks_whitegrid
from genomic_regions import GenomicRegion
from ..matrix import RegionMatrixTable, RegionMatrixContainer
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import matplotlib.patches as patches
from abc import ABCMeta
import numpy as np
from scipy.ndimage.filters import gaussian_filter
import itertools as it
import types
import seaborn as sns
from future.utils import with_metaclass, string_types
import logging
logger = logging.getLogger(__name__)


def prepare_hic_buffer(hic_data, buffering_strategy="relative", buffering_arg=1,
                       weight_field=None, default_value=None, smooth_sigma=None,
                       norm=True, oe=False, log=False):
    """
    Prepare :class:`~BufferedMatrix` from hic data.

    :param hic_data: :class:`~kaic.data.genomic.RegionMatrixTable` or
                     :class:`~kaic.data.genomic.RegionMatrix`
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
    if isinstance(hic_data, RegionMatrixTable):
        return BufferedMatrix(hic_data, buffering_strategy=buffering_strategy,
                              buffering_arg=buffering_arg, weight_field=weight_field,
                              default_value=default_value, smooth_sigma=smooth_sigma,
                              norm=norm, oe=oe, log=log)
    elif isinstance(hic_data, RegionMatrixContainer):
        return BufferedMatrix.from_hic_matrix(hic_data, weight_field=weight_field,
                                              default_value=default_value,
                                              smooth_sigma=smooth_sigma,
                                              norm=norm, oe=oe, log=log)
    else:
        raise ValueError("Unknown type for hic_data")


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

    def __init__(self, data, buffering_strategy="relative", buffering_arg=1,
                 weight_field=None, default_value=None, smooth_sigma=None,
                 norm=True, oe=False, log=False):
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
        regions = tuple(reversed([r for r in regions]))
        if not self.is_buffered_region(*regions):
            logger.debug("Buffering matrix")
            self._BUFFERING_STRATEGIES[self.buffering_strategy](self, *regions)
        m = self.buffered_matrix[tuple(regions)]
        if self.smooth_sigma is not None:
            mf = gaussian_filter(m, self.smooth_sigma)
            m = kaic.matrix.RegionMatrix(mf, row_regions=m.row_regions, col_regions=m.col_regions)
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
    def __init__(self, top_matrix, bottom_matrix, scale_matrices=True, buffering_strategy="relative", buffering_arg=1):
        super(BufferedCombinedMatrix, self).__init__(None, buffering_strategy, buffering_arg)

        scaling_factor = 1
        if scale_matrices:
            scaling_factor = top_matrix.scaling_factor(bottom_matrix)

        class CombinedData(object):
            def __init__(self, hic_top, hic_bottom, scaling_factor=1):
                self.hic_top = hic_top
                self.hic_bottom = hic_bottom
                self.scaling_factor = scaling_factor

            def __getitem__(self, item):
                return top_matrix.get_combined_matrix(self.hic_bottom, key=item, scaling_factor=self.scaling_factor)

        self.data = CombinedData(top_matrix, bottom_matrix, scaling_factor)


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
                        :class:`~kaic.data.genomic.Hic` or :class:`~kaic.data.genomic.RegionMatrix`
        :param adjust_range: Draw a slider to adjust vmin/vmax interactively. Default: False
        :param buffering_strategy: A valid buffering strategy for :class:`~BufferedMatrix`
        :param buffering_arg: Adjust range of buffering for :class:`~BufferedMatrix`
        """
        super(BasePlotterHic, self).__init__(**kwargs)
        if isinstance(hic_data, string_types):
            hic_data = kaic.load(hic_data, mode="r")
        self.hic_data = hic_data
        self.hic_buffer = prepare_hic_buffer(hic_data, buffering_strategy=buffering_strategy,
                                             buffering_arg=buffering_arg, weight_field=weight_field,
                                             default_value=default_value, smooth_sigma=smooth_sigma,
                                             norm=matrix_norm, oe=oe, log=log)
        self.slider = None
        self.adjust_range = adjust_range
        self.vmax_slider = None


class HicPlot2D(BasePlotterHic, BasePlotter2D):
    """
    Plot Hi-C map as a square matrix.
    """
    sharey = True

    def __init__(self, hic_data, flip=False, **kwargs):
        """
        :param flip: Transpose matrix before plotting
        """
        super(HicPlot2D, self).__init__(hic_data=hic_data, **kwargs)
        self.vmax_slider = None
        self.current_matrix = None
        self.flip = flip

    def _plot(self, region):
        self.current_matrix = self.hic_buffer.get_matrix(*region)
        m = np.transpose(self.current_matrix)

        if self.flip:
            m = np.flipud(m)
            extent = [self.current_matrix.row_regions[0].end, self.current_matrix.row_regions[-1].end,
                      self.current_matrix.col_regions[0].start, self.current_matrix.col_regions[-1].end]
        else:
            extent = [self.current_matrix.row_regions[0].end, self.current_matrix.row_regions[-1].end,
                      self.current_matrix.col_regions[-1].end, self.current_matrix.col_regions[0].start]

        color_matrix = self.get_color_matrix(m)
        self.im = self.ax.imshow(color_matrix, interpolation='none', aspect='auto',
                                 cmap=self.colormap, norm=self.norm, extent=extent)

        self.ax.spines['right'].set_visible(False)
        self.ax.spines['top'].set_visible(False)
        self.ax.xaxis.set_ticks_position('bottom')
        self.ax.yaxis.set_ticks_position('left')

        self.ax.set_xlim(extent[0], extent[1])
        self.ax.set_ylim(extent[2], extent[3])

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
        if self.adjust_range:
            self.add_adj_slider()

    def _slider_refresh(self, val):
        # new_vmin = self.vmin_slider.val
        new_vmax = self.vmax_slider.val
        if self.colorbar_symmetry is not None:
            diff = abs(self.colorbar_symmetry - new_vmax)
            self.im.set_clim(vmin=self.colorbar_symmetry - diff, vmax=self.colorbar_symmetry + diff)
        else:
            self.im.set_clim(vmin=self.hic_buffer.buffered_min, vmax=new_vmax)
        # Hack to force redraw of image data
        self.im.set_data(self.current_matrix)

        if self.colorbar is not None:
            self.update_colorbar(vmax=new_vmax)

    def _refresh(self, region):
        self.current_matrix = self.hic_buffer.get_matrix(*region)
        m = np.transpose(self.current_matrix)

        if self.flip:
            m = np.flipud(m)
            extent = [self.current_matrix.row_regions[0].end, self.current_matrix.row_regions[-1].end,
                      self.current_matrix.col_regions[0].start, self.current_matrix.col_regions[-1].end]
        else:
            extent = [self.current_matrix.row_regions[0].end, self.current_matrix.row_regions[-1].end,
                      self.current_matrix.col_regions[-1].end, self.current_matrix.col_regions[0].start]

        color_matrix = self.get_color_matrix(m)
        self.im.set_data(color_matrix)
        self.im.set_extent(extent)


class HicSideBySidePlot2D(object):
    def __init__(self, hic1, hic2, colormap=config.colormap_hic, norm="log",
                 vmin=None, vmax=None, aspect=1., axes_style="ticks"):
        self.hic_plotter1 = HicPlot2D(hic1, colormap=colormap, norm=norm,
                                      vmin=vmin, vmax=vmax, aspect=aspect, axes_style=axes_style)
        self.hic_plotter2 = HicPlot2D(hic2, colormap=colormap, norm=norm,
                                      vmin=vmin, vmax=vmax, aspect=aspect, axes_style=axes_style)

    def plot(self, region):
        fig = plt.figure()
        ax1 = plt.subplot(121)
        ax2 = plt.subplot(122, sharex=ax1, sharey=ax1)

        self.hic_plotter1.plot(x_region=region, y_region=region, ax=ax1)
        self.hic_plotter2.plot(x_region=region, y_region=region, ax=ax2)

        return fig, ax1, ax2


class HicComparisonPlot2D(HicPlot2D):
    def __init__(self, hic_top, hic_bottom, buffering_strategy="relative", 
                 buffering_arg=1., scale_matrices=True, **kwargs):
        super(HicComparisonPlot2D, self).__init__(hic_top, **kwargs)
        self.hic_top = hic_top
        self.hic_bottom = hic_bottom
        self.hic_buffer = BufferedCombinedMatrix(hic_bottom, hic_top, scale_matrices,
                                                 buffering_strategy, buffering_arg)


class HicSlicePlot(ScalarDataPlot):
    """
    Draw a Hi-C data as virtual 4C-plot. All interactions that 
    involve the slice region are shown.
    """

    def __init__(self, hic_data, slice_region, names=None,
                 buffering_strategy="relative", buffering_arg=1,
                 weight_field=None, default_value=None, **kwargs):
        """
        :param hic_data: :class:`~kaic.Hic` or :class:`~kaic.RegionMatrix`. Can be list of
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

    def _plot(self, region):
        for i, b in enumerate(self.hic_buffers):
            hm = np.mean(b.get_matrix(self.slice_region, region).T, axis=0)
            bin_coords = np.r_[[x.start for x in hm.row_regions], hm.row_regions[-1].end]
            bin_coords = (bin_coords[1:] + bin_coords[:-1])/2
            self.ax.plot(bin_coords, hm, label=self.names[i] if self.names else "")
        if self.names:
            self.add_legend()
        self.remove_colorbar_ax()
        sns.despine(ax=self.ax, top=True, right=True)

    def _refresh(self, region):
        pass


class HicPlot(BasePlotterHic, BasePlotter1D):
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
        kwargs.setdefault("aspect", .5)
        super(HicPlot, self).__init__(hic_data=hic_data, **kwargs)
        self.proportional = proportional
        self.max_dist = max_dist
        self.hm = None
        self.rasterized = rasterized

    def _plot(self, region):
        logger.debug("Generating matrix from hic object")
        if region is None:
            raise ValueError("Cannot plot triangle plot for whole genome.")
        if region.start is None:
            region.start = 1
        if region.end is None:
            region.end = self.hic_data.chromosome_lengths[region.chromosome]
        if self.proportional:
            if self.max_dist is None:
                self.aspect = .5
            else:
                rl = region.end - region.start
                self.aspect = .5*min(self.max_dist, rl)/rl
            self._dimensions_stale = True

        # Have to copy unfortunately, otherwise modifying matrix in buffer
        x_, y_, hm = self._mesh_data(region)
        self.hm = hm

        self.collection = self.ax.pcolormesh(x_, y_, hm, cmap=self.colormap, norm=self.norm,
                                             rasterized=self.rasterized)
        self.collection._A = None
        self._update_mesh_colors()

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
        hm_copy = kaic.matrix.RegionMatrix(np.copy(hm), col_regions=hm.col_regions,
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

        self.collection._coordinates[:, :, 0] = x_
        # update matrix data
        self.collection.set_array(self.hm.ravel())
        self._update_mesh_colors()

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


class HicPeakPlot(BaseOverlayPlotter):
    """
    Overlay peaks onto Hicplot or HicPlot2D. Accepts
    :class:`~kaic.data.network.PeakInfo`.
    Add to HicPlot or HicPlot2D using add_overlay method.
    """
    def __init__(self, peaks, radius=None, circle_props={}, **kwargs):
        """
        :param peaks: Kaic peaks instance
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
