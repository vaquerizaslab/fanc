"""

Kai-C
====

Provides
    1. Classes for working with Hi-C data
    2. Classes for working with tabular data

"""
import logging
import os

from genomic_regions import GenomicRegion, load as gr_load
from kaic.architecture.genome_architecture import GenomicTrack
from kaic.architecture.hic_architecture import DirectionalityIndex, InsulationIndex, PossibleContacts, \
    ExpectedContacts, RegionContactAverage, FoldChangeMatrix, ObservedExpectedRatio, ABDomains, \
    ABDomainMatrix, MetaArray, MetaHeatmap, VectorDifference, VectorArchitecturalRegionFeature, \
    MultiVectorArchitecturalRegionFeature, cumulative_matrix
from kaic.config import config
from .general import FileBased
from .hic import Hic
from .matrix import Edge
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
    import os
    file_name = os.path.expanduser(file_name)
    try:
        f = FileBased(file_name, mode='r')
        mode = kwargs.pop('mode', 'r')
        classid = None
        try:
            classid = f.meta._classid
            classid = classid.decode() if isinstance(classid, bytes) else classid
            f.close()
            cls_ = class_id_dict[classid]
            logger.debug("Detected {}".format(cls_))
            return cls_(file_name=file_name, mode=mode, *args, **kwargs)
        except AttributeError:
            pass

            # others
            detectables = (
                ('insulation_index', InsulationIndex),
                ('directionality_index', DirectionalityIndex),
                ('contact_average', RegionContactAverage),
                ('expected_contacts', FoldChangeMatrix),
                ('distance_decay', ExpectedContacts),
                ('observed_expected', ObservedExpectedRatio),
                ('ab_domains', ABDomains),
                ('ab_domain_matrix', ABDomainMatrix),
                ('possible_contacts', PossibleContacts),
                ('meta_matrix', MetaArray),
                ('meta_heatmap', MetaHeatmap),
                ('tracks', GenomicTrack),
                ('fragments', Pairs),
                ('vector_diff', VectorDifference),
                ('region_data', VectorArchitecturalRegionFeature),
                ('array_region_data', MultiVectorArchitecturalRegionFeature),
            )

            for name, cls in detectables:
                try:
                    f.file.get_node('/' + name)
                    f.close()
                    return cls(file_name, mode=mode, *args, **kwargs)
                except tables.NoSuchNodeError:
                    pass

            f.close()
            raise ValueError("File ({}) does not have a '_classid' meta attribute. This might be fixed by loading the "
                             "class once explicitly with the appropriate class in append mode. "
                             "It was also impossible to auto-detect the file type from the file "
                             "structure.".format(file_name))
        except KeyError:
            raise ValueError("classid attribute ({}) does not have a registered class.".format(classid))
    except tables.HDF5ExtError:
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
