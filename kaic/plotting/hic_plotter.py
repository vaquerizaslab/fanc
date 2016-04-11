import kaic
from kaic.plotting.base_plotter import BasePlotterMatrix, BasePlotter1D, BasePlotter2D, append_axes
from kaic.data.genomic import GenomicRegion
import matplotlib as mpl
from matplotlib.widgets import Slider
from abc import ABCMeta
import numpy as np
import itertools as it
import types
import logging
import seaborn as sns
plt = sns.plt


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
        """
        Wrap a :class:`~HicMatrix` in a :class:`~BufferedMatrix` container.
        Default buffering strategy is set to "all" by default.

        :param hic_matrix: :class:`~HicMatrix`
        :return: :class:`~BufferedMatrix`
        """
        bm = cls(data=None, buffering_strategy="all")
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
                 all(rb.contains(rq) for rb, rq in it.izip(self.buffered_region, regions)))):
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
            logging.info("Buffering matrix")
            self._BUFFERING_STRATEGIES[self.buffering_strategy](self, *regions)
        return self.buffered_matrix[tuple(regions)]

    def _buffer_all(self, *regions):
        """
        No buffering, just loads everything in the object into memory.

        Obviously very memory intensive.

        :param regions: :class:`~GenomicRegion` objects
        :return: :class:`~HicMatrix`
        """
        self.buffered_region = self._STRATEGY_ALL
        self.buffered_matrix = self.data[tuple([slice(0, None, None)]*len(regions))]

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
                self.buffered_region.append(GenomicRegion(start=new_start, end=new_end, chromosome=rq.chromosome))
            else:
                self.buffered_region.append(GenomicRegion(start=None, end=None, chromosome=rq.chromosome))
        self.buffered_matrix = self.data[tuple(self.buffered_region)]

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
                self.buffered_region.append(GenomicRegion(start=new_start, end=new_end, chromosome=rq.chromosome))
            else:
                self.buffered_region.append(GenomicRegion(start=None, end=None, chromosome=rq.chromosome))
        self.buffered_matrix = self.data[tuple(self.buffered_region)]

    @property
    def buffered_min(self):
        """
        Find the smallest non-zero buffered matrix value.

        :return: float or None if nothing is buffered
        """
        return np.ma.min(self.buffered_matrix[np.ma.nonzero(self.buffered_matrix)])\
            if self.buffered_matrix is not None else None

    @property
    def buffered_max(self):
        """
        Find the largest buffered matrix value
        :return: float or None if nothing is buffered
        """
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
        self.vmax_slider = None


class HicPlot2D(BasePlotter2D, BasePlotterHic):
    def __init__(self, hic_data, title='', colormap='viridis', norm="log",
                 vmin=None, vmax=None, show_colorbar=True,
                 adjust_range=True, buffering_strategy="relative", buffering_arg=1,
                 blend_zero=True, unmappable_color=".9",
                 aspect=1., axes_style="ticks"):
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
        BasePlotter2D.__init__(self, title=title, aspect=aspect, axes_style=axes_style)
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
        self.ax.spines['right'].set_visible(False)
        self.ax.spines['top'].set_visible(False)
        self.ax.xaxis.set_ticks_position('bottom')
        self.ax.yaxis.set_ticks_position('left')

        self.last_ylim = self.ax.get_ylim()
        self.last_xlim = self.ax.get_xlim()

        if self.show_colorbar:
            self.add_colorbar()
            if self.adjust_range:
                self.add_adj_slider()

    def _refresh(self, x_region=None, y_region=None):
        m = self.hic_buffer.get_matrix(x_region, y_region)
        self.im.set_data(self.get_color_matrix(m))
        self.im.set_extent([m.col_regions[0].start, m.col_regions[-1].end,
                            m.row_regions[-1].end, m.row_regions[0].start])


class HicSideBySidePlot2D(object):
    def __init__(self, hic1, hic2, colormap='viridis', norm="log",
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
    def __init__(self, hic_top, hic_bottom, colormap='viridis', norm='log',
                 vmin=None, vmax=None, scale_matrices=True, show_colorbar=True,
                 buffering_strategy="relative", buffering_arg=1, aspect=1.,
                 axes_style="ticks"):
        super(HicComparisonPlot2D, self).__init__(hic_top, colormap=colormap, norm=norm,
                                                  vmin=vmin, vmax=vmax,
                                                  show_colorbar=show_colorbar, aspect=aspect,
                                                  axes_style=axes_style)
        self.hic_top = hic_top
        self.hic_bottom = hic_bottom
        self.hic_buffer = BufferedCombinedMatrix(hic_top, hic_bottom, scale_matrices, buffering_strategy, buffering_arg)


class HicPlot(BasePlotter1D, BasePlotterHic):
    def __init__(self, hic_data, title='', colormap='viridis', max_dist=None, norm="log",
                 vmin=None, vmax=None, show_colorbar=True, adjust_range=False,
                 buffering_strategy="relative", buffering_arg=1, blend_zero=True,
                 unmappable_color=".9", illegal_color=None, aspect=.5,
                 axes_style="ticks"):
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
        :param illegal_color: Draw non-finite (NaN, +inf, -inf) bins using this color. Defaults to
                                 None (no special color).
        :param aspect: Default aspect ratio of the plot. Can be overriden by setting
                       the height_ratios in class:`~GenomicFigure`
        """
        BasePlotter1D.__init__(self, title=title, aspect=aspect, axes_style=axes_style)
        BasePlotterHic.__init__(self, hic_data, colormap=colormap, vmin=vmin, vmax=vmax,
                                show_colorbar=show_colorbar, adjust_range=adjust_range,
                                buffering_strategy=buffering_strategy, buffering_arg=buffering_arg,
                                norm=norm, blend_zero=blend_zero, unmappable_color=unmappable_color,
                                illegal_color=illegal_color)
        self.max_dist = max_dist
        self.hm = None

    def _plot(self, region=None, *args, **kwargs):
        logging.debug("Generating matrix from hic object")
        if region is None:
            raise ValueError("Cannot plot triangle plot for whole genome.")
        # Have to copy unfortunately, otherwise modifying matrix in buffer
        x_, y_, hm = self._mesh_data(region)
        self.hm = hm

        self.collection = self.ax.pcolormesh(x_, y_, hm, cmap=self.colormap, norm=self.norm, rasterized=True)
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
        if self.show_colorbar:
            cax = None
            if isinstance(self.show_colorbar, mpl.axes.Axes):
                cax = self.show_colorbar
            self.add_colorbar(ax=cax)
        if self.adjust_range:
            self.add_adj_slider()

        def drag_pan(self, button, key, x, y):
            mpl.axes.Axes.drag_pan(self, button, 'x', x, y)  # pretend key=='x'

        self.ax.drag_pan = types.MethodType(drag_pan, self.ax)

    def _mesh_data(self, region):
        hm = self.hic_buffer.get_matrix(region, region)
        hm_copy = kaic.data.genomic.RegionMatrix(np.copy(hm), col_regions=hm.col_regions,
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
        self.collection.set_color(color_tuple)

    def _refresh(self, region=None, *args, **kwargs):
        x_, y_, hm = self._mesh_data(region)
        self.hm = hm

        self.collection._coordinates[:, :, 0] = x_
        # update matrix data
        self.collection.set_array(self.hm.ravel())
        self._update_mesh_colors()

    def add_adj_slider(self, ax=None):
        if ax is None:
            ax = append_axes(self.ax, 'top', 1, 0.05)

        self.vmax_slider = Slider(ax, 'vmax', self.hic_buffer.buffered_min,
                                  self.hic_buffer.buffered_max, valinit=self.vmax,
                                  facecolor='#dddddd', edgecolor='none')

        self.vmax_slider.on_changed(self._slider_refresh)

    def _slider_refresh(self, val):
        # new_vmin = self.vmin_slider.val
        new_vmax = self.vmax_slider.val
        self._update_norm(vmax=new_vmax)
        self._update_mesh_colors()
        self.colorbar.set_clim(vmax=new_vmax)
        self.colorbar.draw_all()
