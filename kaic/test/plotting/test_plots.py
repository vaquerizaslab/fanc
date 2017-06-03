from __future__ import division
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import kaic
import kaic.plotting as kplot
import os.path
import pytest
import numpy as np


def get_example_hic():
    directory = os.path.dirname(os.path.realpath(__file__))
    h = kaic.Hic(directory + "/../data/test_network/rao2014.chr11_77400000_78600000.hic", mode='r')
    return h

def get_example_bigwig():
    directory = os.path.dirname(os.path.realpath(__file__))
    h = directory + "/../data/test_plotting/sample_bigwig.bigwig"
    return h

def get_example_bedgraph():
    directory = os.path.dirname(os.path.realpath(__file__))
    h = directory + "/../data/test_plotting/sample_bigwig.bedgraph"
    return h

def abs_ax_aspect(ax):
    bbox = ax.get_position()
    fig = ax.figure
    aspect = (bbox.height*fig.get_figheight())/(bbox.width*fig.get_figwidth())
    return aspect

@pytest.mark.plotting
class TestHicPlot:
    def setup_method(self, method):
        self.hic = get_example_hic()
        self.hic_matrix = self.hic[:]
        self.hic_matrix[10, :] = 0
        self.hic_matrix[:, 10] = 0

    def teardown_method(self, method):
        self.hic.close()

    @pytest.mark.parametrize("norm", ["lin", "log"])
    @pytest.mark.parametrize("max_dist", [None, 400000])
    @pytest.mark.parametrize("colormap", ["viridis", mpl.cm.get_cmap("Reds")])
    @pytest.mark.parametrize("colorbar", [True, False])
    @pytest.mark.parametrize("blend_zero", [True, False])
    @pytest.mark.parametrize("aspect", [.345])
    @pytest.mark.parametrize("vrange", [(None, .3), (.01, .4)])
    @pytest.mark.parametrize("crange", [(77390001, 78600000), (77800000, 78000000)])
    @pytest.mark.parametrize("unmappable_color", [".345"])
    def test_hicplot_inputs(self, norm, max_dist, colormap, colorbar,
                            blend_zero, aspect, vrange, crange, unmappable_color):
        start, end = crange
        vmin, vmax = vrange
        hplot = kplot.HicPlot(hic_data=self.hic_matrix, title="quark", norm=norm, vmin=vmin, vmax=vmax,
                              max_dist=max_dist, colormap=colormap, show_colorbar=colorbar,
                              blend_zero=blend_zero, aspect=aspect, unmappable_color=unmappable_color)
        gfig = kplot.GenomicFigure([hplot])
        selector = "chr11:{}-{}".format(start, end)
        fig, axes = gfig.plot(selector)
        assert axes[0].get_title() == "quark"
        norm_values = {"lin": mpl.colors.Normalize,
                       "log": mpl.colors.LogNorm}
        assert isinstance(hplot.collection.norm, norm_values[norm])
        assert axes[0].get_xlim() == (start, end)
        if vmin is not None:
            assert hplot.collection.norm.vmin == vmin
        if vmax is not None:
            assert hplot.collection.norm.vmax == vmax
        assert hplot.aspect == aspect
        assert abs_ax_aspect(axes[0]) == pytest.approx(aspect)
        colorbar_values = {True: mpl.colorbar.Colorbar,
                           False: type(None)}
        assert isinstance(hplot.colorbar, colorbar_values[colorbar])
        hic_matrix = self.hic_matrix[selector, selector]
        zero_mask = np.isclose(hic_matrix, 0.)
        unmap_mask = np.all(zero_mask, axis=0)
        unmap_mask = np.logical_or(unmap_mask, unmap_mask[:, np.newaxis])
        if unmappable_color is not None:
            unmap_color = mpl.colors.colorConverter.to_rgba(unmappable_color)
            assert np.all(np.isclose(hplot.collection.get_facecolors()[unmap_mask.T.ravel(), :], unmap_color))
        # blend zero only really makes sense with log normalization, with 
        # linear normalization 0 values map to the first colormap value anyway
        if norm == "log":
            exp_zero = np.logical_and(zero_mask, np.logical_not(unmap_mask))
            if np.sum(exp_zero) > 0:
                zero_value = hplot.colormap(0)
                zero_blended = np.all(np.isclose(hplot.collection.get_facecolors()[exp_zero.T.ravel(), :], zero_value))
                assert blend_zero == zero_blended
        plt.close(fig)

    @pytest.mark.parametrize("norm", ["lin", "log"])
    @pytest.mark.parametrize("colormap", ["viridis", mpl.cm.get_cmap("Reds")])
    @pytest.mark.parametrize("colorbar", [True, False])
    @pytest.mark.parametrize("blend_zero", [True, False])
    @pytest.mark.parametrize("aspect", [.345])
    @pytest.mark.parametrize("vrange", [(None, .3), (.01, .4)])
    @pytest.mark.parametrize("crange", [(77390001, 78600000), (77800000, 78000000)])
    @pytest.mark.parametrize("unmappable_color", [".345"])
    def test_hicplot2d_inputs(self, norm, colormap, colorbar,
                            blend_zero, aspect, vrange, crange, unmappable_color):
        start, end = crange
        vmin, vmax = vrange
        hplot = kplot.HicPlot2D(hic_data=self.hic_matrix, title="quark", norm=norm, vmin=vmin, vmax=vmax,
                                colormap=colormap, show_colorbar=colorbar,
                                blend_zero=blend_zero, aspect=aspect, unmappable_color=unmappable_color)
        gfig = kplot.GenomicFigure([hplot])
        selector = "chr11:{}-{}".format(start, end)
        fig, axes = gfig.plot(selector)
        assert axes[0].get_title() == "quark"
        norm_values = {"lin": mpl.colors.Normalize,
                       "log": mpl.colors.LogNorm}
        assert isinstance(hplot.norm, norm_values[norm])
        assert axes[0].get_xlim() == (start, end)
        assert axes[0].get_ylim() == (start, end)
        if vmin is not None:
            assert hplot.im.norm.vmin == vmin
        if vmax is not None:
            assert hplot.im.norm.vmax == vmax
        assert hplot.aspect == aspect
        assert abs_ax_aspect(axes[0]) == pytest.approx(aspect)
        colorbar_values = {True: mpl.colorbar.Colorbar,
                           False: type(None)}
        assert isinstance(hplot.colorbar, colorbar_values[colorbar])
        hic_matrix = self.hic_matrix[selector, selector]
        zero_mask = np.isclose(hic_matrix, 0.)
        unmap_mask = np.all(zero_mask, axis=0)
        unmap_mask = np.logical_or(unmap_mask, unmap_mask[:, np.newaxis])
        if unmappable_color is not None:
            unmap_color = mpl.colors.colorConverter.to_rgba(unmappable_color)
            assert np.all(np.isclose(hplot.im.get_array()[unmap_mask, :], unmap_color))
        # blend zero only really makes sense with log normalization, with
        # linear normalization 0 values map to the first colormap value anyway
        if norm == "log":
            exp_zero = np.logical_and(zero_mask, np.logical_not(unmap_mask))
            if np.sum(exp_zero) > 0:
                zero_value = hplot.colormap(0)
                zero_blended = np.all(np.isclose(hplot.im.get_array()[exp_zero, :], zero_value))
                assert blend_zero == zero_blended
        plt.close(fig)

    @pytest.mark.parametrize("n_hic", [1, 3])
    @pytest.mark.parametrize("yscale", ["linear", "log"])
    @pytest.mark.parametrize("ylim", [(0, 10), (None, 5)])
    @pytest.mark.parametrize("slice_region", ["chr11:77800000-77850000"])
    @pytest.mark.parametrize("names", [None, True])
    def test_hicsliceplot(self, n_hic, yscale, ylim, slice_region, names):
        hic_data = self.hic_matrix if n_hic == 1 else [self.hic_matrix]*n_hic
        names_passed = None if names is None else ["hic{}".format(i) for i in range(n_hic)]
        splot = kplot.HicSlicePlot(hic_data, slice_region, names=names_passed,
                                   ylim=ylim, yscale=yscale)
        gfig = kplot.GenomicFigure([splot])
        selector = "chr11:{}-{}".format(77400000, 78600000)
        fig, axes = gfig.plot(selector)
        assert axes[0].get_yscale() == yscale
        for a, b in zip(ylim, axes[0].get_ylim()):
            assert a is None or a == b
        if names is not None:
            assert all(l.get_label() == n_p for l, n_p in zip(axes[0].get_lines(), names_passed))


@pytest.mark.plotting
class TestPlots:
    def setup_method(self, method):
        self.bigwig_path = get_example_bigwig()
        self.pyBigWig = pytest.importorskip("pyBigWig")
        self.bigwig = self.pyBigWig.open(self.bigwig_path)
        self.bedgraph_path = get_example_bedgraph()

    def teardown_method(self, method):
        self.bigwig.close()

    @pytest.mark.parametrize("file", ["bigwig", "bedgraph"])
    @pytest.mark.parametrize("crange", [(10, 700)])
    @pytest.mark.parametrize("n_bw", [1, 3])
    @pytest.mark.parametrize("ylim", [(0, 10), (None, 5)])
    @pytest.mark.parametrize("yscale", ["linear", "log"])
    @pytest.mark.parametrize("names", [None, True])
    @pytest.mark.parametrize("bin_size", [50, 2])
    def test_bigwig_plot(self, file, crange, n_bw, ylim, yscale, names, bin_size):
        path = getattr(self, "{}_path".format(file))
        bw_data = path if n_bw == 1 else [path]*n_bw
        names_passed = None if names is None else ["bw{}".format(i) for i in range(n_bw)]
        bplot = kplot.BigWigPlot(bw_data, names=names_passed, bin_size=bin_size,
                                  ylim=ylim, yscale=yscale)
        gfig = kplot.GenomicFigure([bplot])
        fig, axes = gfig.plot("chr2:{}-{}".format(*crange))
        assert axes[0].get_yscale() == yscale
        for a, b in zip(ylim, axes[0].get_ylim()):
            assert a is None or a == b
        if names is not None:
            assert all(l.get_label() == n_p for l, n_p in zip(axes[0].get_lines(), names_passed))
        # Should figure out exactly how many data points I expect instead...
        assert all(len(l.get_ydata()) > 5 for l in axes[0].get_lines())
