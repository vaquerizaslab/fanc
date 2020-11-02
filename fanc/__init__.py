"""

FAN-C
=====

Provides
    1. Classes for working with Hi-C data
    2. Classes for working with tabular data

"""
import logging
import os

# show only warnings for numexpr
logging.getLogger("numexpr").setLevel(logging.WARNING)

from genomic_regions import GenomicRegion
from .architecture.domains import InsulationScore, InsulationScores, DirectionalityIndex, \
    DirectionalityIndexes, Boundaries
from .architecture.compartments import ABCompartmentMatrix
from .architecture.comparisons import FoldChangeMatrix, DifferenceMatrix, DifferenceRegions, \
    FoldChangeRegions, DifferenceScores, FoldChangeScores
from .architecture.aggregate import AggregateMatrix, aggregate_boundaries, aggregate_loops
from .config import config
from .general import FileBased
from .map import *
from .hic import Hic
from .matrix import Edge, RegionMatrix
from .pairs import ReadPairs
from .peaks import RaoPeakInfo, RaoPeakCaller, RaoPeakFilter
from .regions import Genome, Chromosome
from .registry import class_id_dict
from .version import __version__
from .tools.load import load

# configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

example_data = dict(
    hic="test/data/test_network/rao2014.chr11_77400000_78600000.hic",
    chip_bigwig="test/data/test_plotting/CTCF_ChIP_FE_chr11_77-80Mb_mouse_embryo_fibroblasts.bigwig",
    chip_bedgraph="test/data/test_plotting/CTCF_ChIP_FE_chr11_77-80Mb_mouse_embryo_fibroblasts.bedgraph.gz",
    chip_peak_bed="test/data/test_plotting/CTCF_ChIP_FE_chr11_77-80Mb_mouse_embryo_fibroblasts.peaks.bed.gz",
    gene_gtf="test/data/test_plotting/genes_mm10_chr11_77-80Mb.gtf.gz",
)
_basepath = os.path.abspath(os.path.dirname(__file__))
example_data = {k: os.path.join(_basepath, v) for k, v in example_data.items()}
