from __future__ import division
import kaic
import kaic.plotting as kplot
import os.path
import pytest
import matplotlib as mpl
import numpy as np

@pytest.fixture
def example_hic(request):
    directory = os.path.dirname(os.path.realpath(__file__))
    h = kaic.Hic(directory + "../test_network/rao2014.chr11_77400000_78600000.hic", mode='r')
    def fin():
        h.close()
    request.addfinalizer(fin)
    return h

def get_example_hic():
    directory = os.path.dirname(os.path.realpath(__file__))
    h = kaic.Hic(directory + "../test_network/rao2014.chr11_77400000_78600000.hic", mode='r')
    return h

class TestHicPlot:
    def setup_method(self, method):
        self.hic = get_example_hic()

    @pytest.mark.parametrize("norm", ["lin", "log"])
    @pytest.mark.parametrize("max_dist", [None, 400000, 800000])
    @pytest.mark.parametrize("colormap", ["viridis", mpl.cm.get_cmap("Reds")])
    @pytest.mark.parametrize("colorbar", [True, False])
    @pytest.mark.parametrize("blend_zero", [True, False])
    @pytest.mark.parametrize("aspect", [.2, 1.])
    @pytest.mark.parametrize("start", [77390001, 77800000])
    @pytest.mark.parametrize("end", [78600000, 78000000])
    def test_hicplot_inputs(self, **kwargs):
        start, end = kwargs.pop("start"), kwargs.pop("end")
        hplot = kplot.HicPlot(hic_data=self.hic, title="quark", **kwargs)
        gfig = kplot.GenomicFigure([hplot])
        selector = "chr11:{}-{}".format(start, end)
        fig, axes = gfig.plot(selector)
        assert axes[0].get_title() == "quark"
        norm_values = {"lin": mpl.colors.Normalize,
                       "log": mpl.colors.LogNorm}
        assert isinstance(hplot.norm, norms[kwargs[norm]])
        assert axes[0].get_xlim() == (kwargs["start"], kwargs["end"])
        assert hplot.get_default_aspect() == kwargs["aspect"]
        colorbar_values = {True: mpl.colorbar.Colorbar,
                           False: NoneType}
        assert isinstance(hplot.colorbar, colorbar_values[kwargs["colorbar"]])
        hic_matrix = self.hic[selector, selector]
        zero_mask = np.isclose(hic_matrix, 0.).T.ravel()
        zero_value = hplot.colormap(0)
        assert all(c == zero_value for c in hplot.collection.get_facecolors()[zero_mask, :])
