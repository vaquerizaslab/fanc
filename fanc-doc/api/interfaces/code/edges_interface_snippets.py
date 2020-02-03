import fanc

# start snippet check
hic = fanc.load("examples/output/hic/binned/fanc_example_500kb.hic")
isinstance(hic, fanc.matrix.RegionPairsContainer)  # True if interface supported
# end snippet check

# start snippet create
# create a few regions
region1 = fanc.GenomicRegion(chromosome='chr1', start=1, end=1000, ix=0)
region2 = fanc.GenomicRegion(chromosome='chr1', start=1001, end=2000, ix=1)
region3 = fanc.GenomicRegion(chromosome='chr2', start=1, end=1000, ix=2)

# connect regions with edges
edge1 = fanc.Edge(region1, region2, weight=10)
edge2 = fanc.Edge(region1, region3, weight=1, foo='test')
edge3 = fanc.Edge(region2, region3, weight=2, bar=[1, 2, 3, 4])
# end snippet create

# start snippet nodes
region1 = edge1.source_node
region2 = edge1.sink_node
# end snippet nodes

# start snippet index
edge1.source  # 0
edge1.sink  # 1
# end snippet index

# start snippet region lookup
source_region = hic.regions[edge1.source]
sink_region = hic.regions[edge1.sink]
# end snippet region lookup

# start snippet fast region access
regions = list(hic.regions)
for edge in hic.edges:
    region1 = regions[edge.source]
    region2 = regions[edge.sink]
# end snippet fast region access

# start snippet only index
edge_ix = fanc.Edge(0, 1, weight=0.2)
# end snippet only index

# start snippet weight bias
edge1.weight  # returns "raw" edge weight multiplied by edge.bias
edge1.bias  # return the "correction factor" that is applied to weight
# end snippet weight bias

# start snippet weight example
edge = fanc.Edge(0, 3, weight=10)
print(edge)  # 0--3; bias: 1.0; weight: 10.0
print(edge.weight)  # 10.0
# end snippet weight example

# start snippet modify bias
edge.bias = 0.5
print(edge.weight)  # 5.0
# end snippet modify bias

# start snippet attributes
edge.foo = 1.
edge.bar = 'test'
edge.baz = [1, 2, 3, 4]
# end snippet attributes


#
# Edge iterators
#

# start snippet edge iterator property
edge_counter = 0
for edge in hic.edges:
    edge_counter += 1
# end snippet edge iterator property

# start snippet edge length
len(hic.edges)
# end snippet edge length

# start snippet edge iterator function
edge_counter = 0
for edge in hic.edges():
    edge_counter += 1
# end snippet edge iterator function

# start snippet edge iterator intra
edge_counter = 0
for edge in hic.edges(inter_chromosomal=False):
    edge_counter += 1
# end snippet edge iterator intra

# start snippet edge iterator chr19
edge_counter = 0
for edge in hic.edges('chr19'):
    edge_counter += 1
# end snippet edge iterator chr19

# start snippet edge iterator 2D selectors
# only return intra-chromosomal edges on chromosome 19
edge_counter = 0
for edge in hic.edges(('chr19', 'chr19')):
    edge_counter += 1

# only return inter-chromosomal edges between chromosome 18 and 19
edge_counter = 0
for edge in hic.edges(('chr18', 'chr19')):
    edge_counter += 1
# end snippet edge iterator 2D selectors

# start snippet edge iterator human
edge_counter = 0
for edge in hic.edges(('chr19:1mb-15mb', 'chr19:30.5mb-45000000')):
    edge_counter += 1
# end snippet edge iterator human

# start snippet edge iterator region
edge_counter = 0
region = fanc.GenomicRegion(chromosome='chr19', start=6000000, end=18000000)
for edge in hic.edges((region, region)):
    edge_counter += 1
# end snippet edge iterator region

# start snippet edge iterator norm
valid_pairs = 0
for edge in hic.edges(norm=False):
    valid_pairs += edge.weight  # here, the edge weight is an unnormalised integer
# end snippet edge iterator norm

# start snippet edge iterator lazy
weights = []
for edge in hic.edges(lazy=True):
    weights.append(edge.weight)
# end snippet edge iterator lazy

# start snippet edge iterator wrong lazy
# THIS IS WRONG!
edges = []
for edge in hic.edges(lazy=True):
    edges.append(edge)

print(edges[0].source, edges[0].sink, edges[0].weight)
# (159, 248, 0.002386864930511163)
print(edges[1].source, edges[1].sink, edges[1].weight)
# (159, 248, 0.002386864930511163)
print(edges[2].source, edges[2].sink, edges[2].weight)
# (159, 248, 0.002386864930511163)
# ...
# end snippet edge iterator wrong lazy

# start snippet edge iterator right lazy
edge_data = []
for edge in hic.edges(lazy=True):
    edge_data.append((edge.source, edge.sink, edge.weight))

print(edge_data[0][0], edge_data[0][1], edge_data[0][2])
# 84 85 0.08227859035794333
print(edge_data[1][0], edge_data[1][1], edge_data[1][2])
# 48 56 0.012760965147510347
print(edge_data[2][0], edge_data[2][1], edge_data[2][2])
# 10 153 0.0056418570804748725
# ...
# end snippet edge iterator right lazy


# start snippet edge iterator modify lazy
for edge in hic.edges(lazy=True):
    edge.weight *= 2
    edge.update()
# end snippet edge iterator modify lazy


# start snippet regions and edges
row_regions, col_regions, edges_iter = hic.regions_and_edges(('chr18', 'chr19'))
for edge in edges_iter:
    row_region = row_regions[edge.source]
    col_region = row_regions[edge.sink]
    # ...
# end snippet regions and edges


# start snippet edge data
weights = list(hic.edge_data('weight', ('chr19', 'chr19')))

# is identical to
weights = []
for weight in hic.edge_data('weight', ('chr19', 'chr19')):
    weights.append(weight)

# is identical to
weights = []
# for edge in hic.edges(('chr19', 'chr19'), lazy=True):
#     weights.append(edge.weight)
# end snippet edge data


# start snippet mappable
m = hic.mappable()
# or for a specific region
m_sub = hic.mappable('chr19')
# end snippet mappable