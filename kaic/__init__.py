"""

Kai-C
====

Provides
    1. Classes for working with Hi-C data
    2. Classes for working with tabular data

"""
from .version import __version__

from kaic.config import config
from kaic.data.genomic import Hic, Node, Edge, Genome, Chromosome, Bed, AccessOptimisedHic, load_hic, GenomicRegion, \
    BigWig
from kaic.data.general import FileBased
from kaic.data.registry import class_id_dict
from kaic.construct.seq import Reads, FragmentMappedReadPairs
from kaic.construct.seq import FragmentMappedReadPairs as Pairs  # alias for FragmentMappedReadPairs
from kaic.architecture.hic_architecture import DirectionalityIndex, InsulationIndex, PossibleContacts, \
    ExpectedContacts, RegionContactAverage, FoldChangeMatrix, ObservedExpectedRatio, ABDomains, \
    ABDomainMatrix, MetaArray, MetaHeatmap, VectorDifference, VectorArchitecturalRegionFeature, \
    MultiVectorArchitecturalRegionFeature
from kaic.data.network import RaoPeakInfo
from kaic.architecture.genome_architecture import GenomicTrack
import tables
import logging

# configure logging
logger = logging.getLogger(__name__)
try:  # Python 2.7+
    from logging import NullHandler
except ImportError:
    class NullHandler(logging.Handler):
        def emit(self, record):
            pass
logger.addHandler(NullHandler())


def load(file_name, mode='a', tmpdir=None):
    import os
    file_name = os.path.expanduser(file_name)
    try:
        f = FileBased(file_name, mode='r')
        classid = None
        try:
            classid = f.meta._classid
            f.close()
            cls_ = class_id_dict[classid]
            logger.debug("Detected {}".format(cls_))
            return cls_(file_name=file_name, mode=mode, tmpdir=tmpdir)
        except AttributeError:
            # try to detect from file structure

            # Hi-C
            try:
                n = f.file.get_node('/edges')
                from kaic.data.general import MaskedTable
                hic_class = None
                if isinstance(n, MaskedTable):
                    hic_class = Hic
                elif isinstance(n, tables.group.Group):
                    hic_class = AccessOptimisedHic

                if hic_class is not None:
                    f.close()
                    return hic_class(file_name, mode=mode, tmpdir=tmpdir)
            except tables.NoSuchNodeError:
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
                ('fragments', FragmentMappedReadPairs),
                ('reads', Reads),
                ('vector_diff', VectorDifference),
                ('region_data', VectorArchitecturalRegionFeature),
                ('array_region_data', MultiVectorArchitecturalRegionFeature),
            )

            for name, cls in detectables:
                try:
                    f.file.get_node('/' + name)
                    f.close()
                    return cls(file_name, mode=mode, tmpdir=tmpdir)
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
        # try some well-known file types

        # SAM/BAM
        import pysam
        try:
            sb = pysam.AlignmentFile(file_name, 'rb')
            if mode != 'rb':
                sb.close()
                sb = pysam.AlignmentFile(file_name, mode)
            return sb
        except (ValueError, IOError):
            pass

        import pybedtools
        f = Bed(file_name)
        try:
            ft = f.file_type
            if ft != 'empty':
                return f
        except (IndexError, pybedtools.MalformedBedLineError):
            pass

        try:
            import pyBigWig
            f = pyBigWig.open(file_name, 'r')
            if mode != 'r':
                f.close()
                f = pyBigWig.open(file_name, mode)

            return BigWig(f)
        except (ImportError, RuntimeError):
            raise ValueError("File type not recognised ({}).".format(file_name))
