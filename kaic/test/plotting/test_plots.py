from __future__ import division
import kaic
import kaic.plotting as kplot
import os.path
import pytest
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# @pytest.fixture
# def example_hic(request):
#     directory = os.path.dirname(os.path.realpath(__file__))
#     h = kaic.Hic(directory + "/../data/test_network/rao2014.chr11_77400000_78600000.hic", mode='r')
#     def fin():
#         h.close()
#     request.addfinalizer(fin)
#     return h

def get_example_hic():
    directory = os.path.dirname(os.path.realpath(__file__))
    h = kaic.Hic(directory + "/../data/test_network/rao2014.chr11_77400000_78600000.hic", mode='r')
    return h

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
        assert hplot.get_default_aspect() == aspect
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
        assert hplot.get_default_aspect() == aspect
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
