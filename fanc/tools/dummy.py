from fanc.hic import LegacyHic, Hic
from fanc.matrix import Edge
from genomic_regions import GenomicRegion


def sample_hic(file_name=None, tmpdir=None):
    hic = Hic(file_name=file_name, tmpdir=tmpdir, mode='w')

    # add some nodes (12 to be exact)
    nodes = []
    for i in range(1, 5000, 1000):
        nodes.append(GenomicRegion(chromosome="chr1", start=i, end=i+1000-1))
    for i in range(1, 3000, 1000):
        nodes.append(GenomicRegion(chromosome="chr2", start=i, end=i+1000-1))
    for i in range(1, 2000, 500):
        nodes.append(GenomicRegion(chromosome="chr3", start=i, end=i+1000-1))
    hic.add_regions(nodes)

    # add some edges with increasing weight for testing
    edges = []
    weight = 1
    for i in range(0, len(nodes)):
        for j in range(i, len(nodes)):
            edges.append(Edge(source=i, sink=j, weight=weight))
            weight += 1

    hic.add_edges(edges)

    return hic


def sample_hic_big(file_name=None, tmpdir=None):
    hic = Hic(file_name=file_name, tmpdir=tmpdir, mode='w')

    # add some nodes (120 to be exact)
    nodes = []
    for i in range(1, 50000, 1000):
        nodes.append(GenomicRegion(chromosome="chr1", start=i, end=i + 1000 - 1))
    for i in range(1, 30000, 1000):
        nodes.append(GenomicRegion(chromosome="chr2", start=i, end=i + 1000 - 1))
    for i in range(1, 20000, 500):
        nodes.append(GenomicRegion(chromosome="chr3", start=i, end=i + 1000 - 1))
    hic.add_regions(nodes)

    # add some edges with increasing weight for testing
    edges = []
    for i in range(0, len(nodes)):
        for j in range(i, len(nodes)):
            weight = i + len(nodes) - (j - i)
            edges.append(Edge(source=i, sink=j, weight=weight))

    hic.add_edges(edges)

    return hic


def sample_fa_hic(file_name=None, zero_indices=set(), tmpdir=None):
    hic = Hic(file_name=file_name, tmpdir=tmpdir, mode='w')

    # add some nodes (120 to be exact)
    nodes = []
    for i in range(1, 50000, 1000):
        nodes.append(GenomicRegion(chromosome="chr1", start=i, end=i + 1000 - 1))
    for i in range(1, 30000, 1000):
        nodes.append(GenomicRegion(chromosome="chr2", start=i, end=i + 1000 - 1))
    for i in range(1, 20000, 500):
        nodes.append(GenomicRegion(chromosome="chr3", start=i, end=i + 1000 - 1))
    hic.add_regions(nodes)

    # add some edges with increasing weight for testing
    edges = []
    weight = 1
    for i in range(0, len(nodes)):
        for j in range(i, len(nodes)):
            if i not in zero_indices and j not in zero_indices:
                edges.append(Edge(source=i, sink=j, weight=weight))
            weight += 1

    hic.add_edges(edges)

    return hic


def sample_homogenous_hic(file_name=None, fill_value=1.0, zero_indices=set(), tmpdir=None):
    hic = Hic(file_name=file_name, tmpdir=tmpdir, mode='w')

    # add some nodes (120 to be exact)
    nodes = []
    for i in range(1, 50000, 1000):
        nodes.append(GenomicRegion(chromosome="chr1", start=i, end=i + 1000 - 1))
    for i in range(1, 30000, 1000):
        nodes.append(GenomicRegion(chromosome="chr2", start=i, end=i + 1000 - 1))
    for i in range(1, 20000, 500):
        nodes.append(GenomicRegion(chromosome="chr3", start=i, end=i + 1000 - 1))
    hic.add_regions(nodes)

    # add some edges with increasing weight for testing
    edges = []
    for i in range(0, len(nodes)):
        for j in range(i, len(nodes)):
            if i not in zero_indices and j not in zero_indices:
                edges.append(Edge(source=i, sink=j, weight=fill_value))

    hic.add_edges(edges)

    return hic


def sample_hic_matrix1(file_name=None, tmpdir=None):
    #     0 1 2 3 4 5 6 7 8 9
    #   #####################
    # 0 # 0 1 0 2 0 3 0 4 0 5
    # 1 #   6 0 7 0 8 0 9 0 1
    # 2 #     2 3 4 0 5 0 0 6
    # 3 #       7 0 8 9 0 1 0
    # 4 #         0 2 3 0 0 4
    # 5 #           5 6 7 8 9
    # 6 #             1 0 0 0
    # 7 #               2 3 0
    # 8 #                 0 4
    # 9 #                   5
    nodes = [
        GenomicRegion('chr1', 1, 1000),
        GenomicRegion('chr1', 1001, 2000),
        GenomicRegion('chr1', 2001, 3000),
        GenomicRegion('chr1', 3001, 4000),
        GenomicRegion('chr1', 4001, 5000),
        GenomicRegion('chr1', 5001, 6000),
        GenomicRegion('chr1', 6001, 7000),
        GenomicRegion('chr1', 7001, 8000),
        GenomicRegion('chr1', 8001, 9000),
        GenomicRegion('chr1', 9001, 10000)
    ]

    edges = [
        Edge(source=0, sink=1, weight=1), Edge(source=3, sink=5, weight=8),
        Edge(source=0, sink=3, weight=2), Edge(source=3, sink=6, weight=9),
        Edge(source=0, sink=5, weight=3), Edge(source=3, sink=8, weight=1),
        Edge(source=0, sink=7, weight=4), Edge(source=4, sink=5, weight=2),
        Edge(source=0, sink=9, weight=5), Edge(source=4, sink=6, weight=3),
        Edge(source=1, sink=1, weight=6), Edge(source=4, sink=9, weight=4),
        Edge(source=1, sink=3, weight=7), Edge(source=5, sink=5, weight=5),
        Edge(source=1, sink=5, weight=8), Edge(source=5, sink=6, weight=6),
        Edge(source=1, sink=7, weight=9), Edge(source=5, sink=7, weight=7),
        Edge(source=1, sink=9, weight=1), Edge(source=5, sink=8, weight=8),
        Edge(source=2, sink=2, weight=2), Edge(source=5, sink=9, weight=9),
        Edge(source=2, sink=3, weight=3), Edge(source=6, sink=6, weight=1),
        Edge(source=2, sink=4, weight=4), Edge(source=7, sink=7, weight=2),
        Edge(source=2, sink=6, weight=5), Edge(source=7, sink=8, weight=3),
        Edge(source=2, sink=9, weight=6), Edge(source=8, sink=9, weight=4),
        Edge(source=3, sink=3, weight=7), Edge(source=9, sink=9, weight=5)
    ]

    hic = Hic(file_name=file_name, tmpdir=tmpdir)
    hic.add_regions(nodes)
    hic.add_edges(edges)

    return hic


def sample_hic_matrix2(file_name=None, tmpdir=None):
    #     0 1 2 3 4 5 6 7 8 9
    #   #####################
    # 0 # 0 1 0 2 - 3 0 - 0 5
    # 1 #   6 0 7 - 8 0 - 0 1
    # 2 #     2 3 - 0 5 - 0 6
    # 3 #       7 - 8 9 - 1 0
    # 4 #         - - - - - -
    # 5 #           5 6 - 8 9
    # 6 #             1 - 0 0
    # 7 #               - - -
    # 8 #                 0 4
    # 9 #                   5
    nodes = [
        GenomicRegion('chr1', 1, 1000),
        GenomicRegion('chr1', 1001, 2000),
        GenomicRegion('chr1', 2001, 3000),
        GenomicRegion('chr1', 3001, 4000),
        GenomicRegion('chr1', 4001, 5000),
        GenomicRegion('chr1', 5001, 6000),
        GenomicRegion('chr1', 6001, 7000),
        GenomicRegion('chr1', 7001, 8000),
        GenomicRegion('chr1', 8001, 9000),
        GenomicRegion('chr1', 9001, 10000)
    ]

    edges = [
        Edge(source=0, sink=1, weight=1), Edge(source=3, sink=5, weight=8),
        Edge(source=0, sink=3, weight=2), Edge(source=3, sink=6, weight=9),
        Edge(source=0, sink=5, weight=3), Edge(source=3, sink=8, weight=1),
        Edge(source=0, sink=9, weight=5),
        Edge(source=1, sink=1, weight=6),
        Edge(source=1, sink=3, weight=7), Edge(source=5, sink=5, weight=5),
        Edge(source=1, sink=5, weight=8), Edge(source=5, sink=6, weight=6),
        Edge(source=1, sink=9, weight=1), Edge(source=5, sink=8, weight=8),
        Edge(source=2, sink=2, weight=2), Edge(source=5, sink=9, weight=9),
        Edge(source=2, sink=3, weight=3), Edge(source=6, sink=6, weight=1),
        Edge(source=2, sink=6, weight=5),
        Edge(source=2, sink=9, weight=6), Edge(source=8, sink=9, weight=4),
        Edge(source=3, sink=3, weight=7), Edge(source=9, sink=9, weight=5)
    ]

    hic = Hic(file_name=file_name, tmpdir=tmpdir)
    hic.add_regions(nodes)
    hic.add_edges(edges)

    return hic
