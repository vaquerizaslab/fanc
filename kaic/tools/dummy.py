from kaic.data.genomic import AccessOptimisedHic, Hic, Node, Edge


def sample_hic(file_name=None, tmpdir=None):
    hic = Hic(file_name=file_name, tmpdir=tmpdir, mode='w')

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


def sample_hic_big(file_name=None, tmpdir=None):
    hic = Hic(file_name=file_name, tmpdir=tmpdir, mode='w')

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


def sample_fa_hic(file_name=None, tmpdir=None):
    hic = AccessOptimisedHic(file_name=file_name, tmpdir=tmpdir, mode='w')

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
