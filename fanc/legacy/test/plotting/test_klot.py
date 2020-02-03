from __future__ import division, print_function
import fanc
import pytest
import subprocess as sp
import uuid
import os

def run_klot(*args):
    args = ["klot"] + [str(x) for x in args]
    return sp.check_call(args)

@pytest.fixture
def output_filename():
    return "{}.pdf".format(uuid.uuid4().hex)

def get_filesize(path):
    return os.stat(str(path)).st_size

@pytest.mark.plotting
@pytest.mark.klot
class TestHicPlot:
    def setup_method(self, method):
        self.hic_path = fanc.example_data["hic"]

    @pytest.mark.parametrize("crange", [(77390001, 78600000)])
    def test_hicplot(self, crange, tmpdir, output_filename):
        start, end = crange
        out_path = tmpdir.join(output_filename)
        # vmin, vmax = vrange
        args = [
            "-o", out_path,
            "chr11:{}-{}".format(start, end),
            "-p", "-t", "hic", self.hic_path,
        ]
        result = run_klot(*args)
        assert get_filesize(out_path) > 10000

    @pytest.mark.parametrize("crange", [(77390001, 78600000)])
    def test_hicplot_2d(self, crange, tmpdir, output_filename):
        start, end = crange
        out_path = tmpdir.join(output_filename)
        # vmin, vmax = vrange
        args = [
            "-o", out_path,
            "chr11:{}-{}".format(start, end),
            "-p", "-t", "hic2d", self.hic_path,
        ]
        result = run_klot(*args)
        assert get_filesize(out_path) > 10000

@pytest.mark.plotting
@pytest.mark.klot
class TestScorePlots:
    def setup_method(self, method):
        self.bigwig_path = self.bigwig_path = fanc.example_data["chip_bigwig"]
        self.bedgraph_path = fanc.example_data["chip_bedgraph"]
        self.peak_path = fanc.example_data["chip_peak_bed"]

    @pytest.mark.parametrize("crange", [(77390001, 78600000)])
    def test_region(self, crange, tmpdir, output_filename):
        start, end = crange
        out_path = tmpdir.join(output_filename)
        # vmin, vmax = vrange
        args = [
            "-o", out_path,
            "chr11:{}-{}".format(start, end),
            "-p", "-t", "region", self.peak_path,
        ]
        result = run_klot(*args)
        assert get_filesize(out_path) > 10000

    @pytest.mark.parametrize("crange", [(77390001, 78600000)])
    def test_bigwig_bedgraph(self, crange, tmpdir, output_filename):
        start, end = crange
        out_path = tmpdir.join(output_filename)
        # vmin, vmax = vrange
        args = [
            "-o", out_path,
            "chr11:{}-{}".format(start, end),
            "-p", "-t", "bigwig", self.bigwig_path,
            "-p", "-t", "bigwig", self.bedgraph_path,
        ]
        result = run_klot(*args)
        assert get_filesize(out_path) > 10000
