from __future__ import division, print_function
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import kaic
import kaic.plotting as kplot
import os.path
import pytest
import numpy as np
import pybedtools as pbt
import subprocess as sp
import uuid
import os

@pytest.fixture(scope="session")
def run_klot():
    def do_run(*args):
        args = ["klot"] + list(args)
        return sp.run(args, check=True)
    return do_run

@pytest.fixture
def output_filename():
    return "{}.pdf".format(uuid.uuid4().hex)

@pytest.mark.plotting
@pytest.mark.klot
class TestHicPlot:
    def setup_method(self, method):
        self.hic_path = kaic.example_data["hic"]

    @pytest.mark.parametrize("crange", [(77390001, 78600000)])
    def test_hicplot(self, crange, tmpdir, run_klot, output_filename):
        start, end = crange
        out_path = tmpdir.join(output_filename)
        # vmin, vmax = vrange
        args = [
            "-o", out_path,
            "chr11:{}-{}".format(start, end),
            "-p", "-t", "hic", self.hic_path,
        ]
        result = run_klot(*args)
        print(out_path)

    @pytest.mark.parametrize("crange", [(77390001, 78600000)])
    def test_hicplot_2d(self, crange, tmpdir, run_klot, output_filename):
        start, end = crange
        out_path = tmpdir.join(output_filename)
        # vmin, vmax = vrange
        args = [
            "-o", out_path,
            "chr11:{}-{}".format(start, end),
            "-p", "-t", "hic2d", self.hic_path,
        ]
        result = run_klot(*args)
        print(out_path)

@pytest.mark.plotting
@pytest.mark.klot
class TestScorePlots:
    def setup_method(self, method):
        self.bigwig_path = self.bigwig_path = kaic.example_data["chip_bigwig"]
        self.bedgraph_path = kaic.example_data["chip_bedgraph"]
        self.peak_path = kaic.example_data["chip_peak_bed"]

    @pytest.mark.parametrize("crange", [(77390001, 78600000)])
    def test_region(self, crange, tmpdir, run_klot, output_filename):
        start, end = crange
        out_path = tmpdir.join(output_filename)
        # vmin, vmax = vrange
        args = [
            "-o", out_path,
            "chr11:{}-{}".format(start, end),
            "-p", "-t", "region", self.peak_path,
        ]
        result = run_klot(*args)
        print(out_path)

    @pytest.mark.parametrize("crange", [(77390001, 78600000)])
    def test_bigwig_bedgraph(self, crange, tmpdir, run_klot, output_filename):
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
        print(out_path)
