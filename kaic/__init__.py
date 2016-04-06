"""

KaiC
====

Provides
    1. Classes for working with Hi-C data
    2. Classes for working with tabular data

"""

from kaic.data.genomic import Hic, Node, Edge, Genome, Chromosome, Bed, AccessOptimisedHic, load_hic
from kaic.data.general import Table 
from kaic.construct.seq import Reads, FragmentMappedReadPairs
import logging
logging.basicConfig(level=logging.INFO)


def sample_hic():
    hic = Hic()

    # add some nodes (12 to be exact)
    nodes = []
    for i in range(1,5000,1000):
        nodes.append(Node(chromosome="chr1",start=i,end=i+1000-1))
    for i in range(1,3000,1000):
        nodes.append(Node(chromosome="chr2",start=i,end=i+1000-1))
    for i in range(1,2000,500):
        nodes.append(Node(chromosome="chr3",start=i,end=i+1000-1))
    hic.add_nodes(nodes)

    # add some edges with increasing weight for testing
    edges = []
    weight = 1
    for i in range(0, len(nodes)):
        for j in range(i, len(nodes)):
            edges.append(Edge(source=i,sink=j,weight=weight))
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