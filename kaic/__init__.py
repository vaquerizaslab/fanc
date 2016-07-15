"""

KaiC
====

Provides
    1. Classes for working with Hi-C data
    2. Classes for working with tabular data

"""

from kaic.data.genomic import Hic, Node, Edge, Genome, Chromosome, Bed, AccessOptimisedHic, load_hic, GenomicRegion
from kaic.data.general import Table, FileBased
from kaic.data.registry import class_id_dict
from kaic.construct.seq import Reads, FragmentMappedReadPairs
from kaic.architecture.hic_architecture import DirectionalityIndex, InsulationIndex, PossibleContacts, \
    ExpectedContacts, load_array
from kaic.architecture.genome_architecture import GenomicTrack
import tables
import logging
logging.basicConfig(level=logging.INFO)


def load(file_name, mode='a', tmpdir=None):
    try:
        f = FileBased(file_name, mode='r')
        classid = None
        try:
            print 'here'
            classid = f.meta._classid
            print 'there'
            f.close()
            cls_ = class_id_dict[classid]
            logging.info("Detected {}".format(cls_))
            print cls_, file_name, mode, tmpdir
            return cls_(file_name=file_name, mode=mode, tmpdir=tmpdir)
        except AttributeError:
            raise ValueError("File ({}) does not have a '_classid' meta attribute. This might be fixed by loading the "
                             "class once explicitly with the appropriate class in append mode.".format(file_name))
        except KeyError:
            raise ValueError("classid attribute ({}) does not have a registered class.".format(classid))
    except tables.HDF5ExtError:
        # try some well-known file types
        import pybedtools
        f = pybedtools.BedTool(file_name)
        try:
            _ = f.file_type
        except IndexError:
            raise ValueError("File type not recognised.")

        return GenomicTrack.from_gtf(file_name=None, gtf_file=file_name)


def sample_hic():
    hic = Hic()

    # add some nodes (12 to be exact)
    nodes = []
    for i in range(1, 5000, 1000):
        nodes.append(Node(chromosome="chr1", start=i, end=i+1000-1))
    for i in range(1, 3000, 1000):
        nodes.append(Node(chromosome="chr2", start=i, end=i+1000-1))
    for i in range(1, 2000, 500):
        nodes.append(Node(chromosome="chr3", start=i, end=i+1000-1))
    hic.add_nodes(nodes)

    # add some edges with increasing weight for testing
    edges = []
    weight = 1
    for i in range(0, len(nodes)):
        for j in range(i, len(nodes)):
            edges.append(Edge(source=i, sink=j, weight=weight))
            weight += 1

    hic.add_edges(edges)

    return hic


def sample_hic_big():
    hic = Hic()

    # add some nodes (120 to be exact)
    nodes = []
    for i in range(1, 50000, 1000):
        nodes.append(Node(chromosome="chr1", start=i, end=i + 1000 - 1))
    for i in range(1, 30000, 1000):
        nodes.append(Node(chromosome="chr2", start=i, end=i + 1000 - 1))
    for i in range(1, 20000, 500):
        nodes.append(Node(chromosome="chr3", start=i, end=i + 1000 - 1))
    hic.add_nodes(nodes)

    # add some edges with increasing weight for testing
    edges = []
    weight = 1
    for i in range(0, len(nodes)):
        for j in range(i, len(nodes)):
            edges.append(Edge(source=i, sink=j, weight=weight))
            weight += 1

    hic.add_edges(edges)

    return hic


def sample_fa_hic():
    hic = AccessOptimisedHic()

    # add some nodes (120 to be exact)
    nodes = []
    for i in range(1, 50000, 1000):
        nodes.append(Node(chromosome="chr1", start=i, end=i + 1000 - 1))
    for i in range(1, 30000, 1000):
        nodes.append(Node(chromosome="chr2", start=i, end=i + 1000 - 1))
    for i in range(1, 20000, 500):
        nodes.append(Node(chromosome="chr3", start=i, end=i + 1000 - 1))
    hic.add_nodes(nodes)

    # add some edges with increasing weight for testing
    edges = []
    weight = 1
    for i in range(0, len(nodes)):
        for j in range(i, len(nodes)):
            edges.append(Edge(source=i, sink=j, weight=weight))
            weight += 1

    hic.add_edges(edges)

    return hic