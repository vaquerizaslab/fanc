"""

Kai-C
====

Provides
    1. Classes for working with Hi-C data
    2. Classes for working with tabular data

"""
import logging
import warnings
import os

from genomic_regions import GenomicRegion, load as gr_load
from .architecture.domains import InsulationScore, InsulationScores, DirectionalityIndex
from .architecture.compartments import ABCompartmentMatrix
from .architecture.comparisons import FoldChangeMatrix, DifferenceMatrix
from .architecture.aggregate import AggregateMatrix, aggregate_boundaries, aggregate_loops
from .config import config
from .general import FileBased
from .map import *
from .hic import Hic
from .matrix import Edge, RegionMatrix
from .pairs import ReadPairs as Pairs
from .peaks import RaoPeakInfo
from .regions import Genome, Chromosome
from .registry import class_id_dict
from .version import __version__
import tables

# configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def load(file_name, *args, **kwargs):
    mode = kwargs.pop('mode', 'r')
    file_name = os.path.expanduser(file_name)

    try:
        logger.debug("Trying FileBased classes")

        f = tables.open_file(file_name, mode='r')
        try:
            classid = f.get_node('/', 'meta_information').meta_node.attrs['_classid']
            classid = classid.decode() if isinstance(classid, bytes) else classid
        finally:
            f.close()
        cls_ = class_id_dict[classid]
        logger.debug("Detected {}".format(cls_))
        return cls_(file_name=file_name, mode=mode, *args, **kwargs)
    except (tables.HDF5ExtError, AttributeError, KeyError):
        pass

    from kaic.compatibility.juicer import JuicerHic, is_juicer
    if is_juicer(file_name):
        return JuicerHic(file_name, *args, **kwargs)

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from kaic.compatibility.cooler import is_cooler, CoolerHic
        if is_cooler(file_name):
            return CoolerHic(file_name, *args, **kwargs)
    except (ImportError, OSError):
        pass

    return gr_load(file_name, *args, **kwargs)


example_data = dict(
    hic="test/data/test_network/rao2014.chr11_77400000_78600000.hic",
    chip_bigwig="test/data/test_plotting/CTCF_ChIP_FE_chr11_77-80Mb_mouse_embryo_fibroblasts.bigwig",
    chip_bedgraph="test/data/test_plotting/CTCF_ChIP_FE_chr11_77-80Mb_mouse_embryo_fibroblasts.bedgraph.gz",
    chip_peak_bed="test/data/test_plotting/CTCF_ChIP_FE_chr11_77-80Mb_mouse_embryo_fibroblasts.peaks.bed.gz",
    gene_gtf="test/data/test_plotting/genes_mm10_chr11_77-80Mb.gtf.gz",
)
_basepath = os.path.abspath(os.path.dirname(__file__))
example_data = {k: os.path.join(_basepath, v) for k, v in example_data.items()}
