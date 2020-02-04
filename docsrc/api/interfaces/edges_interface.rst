.. _edges_interface:

====================
RegionPairsContainer
====================

This interface provides common properties and functions to data based on pairs of regions.
A typical example in this regard would be pairs of ligated fragments in a Hi-C library, as
represented within FAN-C by the :class:`~fanc.pairs.ReadPairs` class. But also matrix-based
data, such as in :class:`~fanc.matrix.Hic` can be interpreted as scores between pairs of
binned genomic regions, hence it also supports the :class:`~fanc.matrix.RegionPairsContainer`
interface. After loading a dataset using :func:`~fanc.tools.load.load`, you can check for
support of the :class:`~fanc.matrix.RegionPairsContainer` interface with:

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet check
    :end-before: end snippet check

The current list of FAN-C classes supporting :class:`~fanc.matrix.RegionPairsContainer` is:
:class:`~fanc.pairs.ReadPairs`,
:class:`~fanc.compatibility.cooler.CoolerHic`,
:class:`~fanc.compatibility.juicer.JuicerHic`,
:class:`~fanc.hic.Hic`,
:class:`~fanc.architecture.compartments.ABCompartmentMatrix`,
:class:`~fanc.architecture.comparisons.DifferenceMatrix`,
:class:`~fanc.architecture.comparisons.FoldChangeMatrix`,
:class:`~fanc.peaks.PeakInfo`,
and
:class:`~fanc.peaks.RaoPeakInfo`.

**************
The Edge class
**************

In :class:`~fanc.matrix.RegionPairsContainer` objects, the basic data type is
:class:`~fanc.matrix.Edge`. It "connects" two :class:`~genomic_regions.GenomicRegion`
objects, called "nodes" within the context of the edge, and supports arbitrary
property annotations, which typically include an edge "weight". Additionally,
regions are typically associated with an index :code:`ix`, which refers to their
position in the region list of the :class:`genomic_regions.RegionBased` container.

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet create
    :end-before: end snippet create

The underlying regions can be accessed using the :code:`source_node` and
:code:`sink_node` properties:

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet nodes
    :end-before: end snippet nodes

Accessing regions in this way can be computationally inefficient (depending on the
internal structure of the container), which is why the recommended way of accessing
regions connected by an edge is through their region index. This is a number describing
each region's location in the list of regions in that container.

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet index
    :end-before: end snippet index

the corresponding regions can then be looked up in the region list:

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet region lookup
    :end-before: end snippet region lookup

of, if processing a large number of edges this way, it is more efficient to obtain
the list of regions in advance:

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet fast region access
    :end-before: end snippet fast region access

It is possible to create an edge using only region indexes, without directly
linking it to a region object:

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet only index
    :end-before: end snippet only index

A call to :code:`source_node` or :code:`sink_node` will then raise a :code:`ValueError`.

:class:`~fanc.matrix.Hic` object edges and many other objects have a "weight"
property, describing, for example, their (normalised) ligation frequency or contact
probability. This properties value is internally multiplied by the value of
:code:`edge.bias` for correction on the fly:

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet weight bias
    :end-before: end snippet weight bias

If, for example, an :class:`~fanc.matrix.Edge` is created like this:

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet weight example
    :end-before: end snippet weight example

its weight will be as assigned at object creation time. We can modify that value
indirectly by changing the bias:

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet modify bias
    :end-before: end snippet modify bias

Other properties are unaffected by the bias value, and can be assigned arbitrarily:

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet attributes
    :end-before: end snippet attributes

Note, however, that not all data types support saving those arbitrary properties in the
container, and that most of the time you are working with copies of data stored in the
container, which remains unmodified even if the edge object is changed. There are exemptions
from this, which will be discussed below.


******************
The edges function
******************

:class:`~fanc.matrix.RegionPairsContainer` compatible objects are built on lists of
regions, which can be accessed and queried using the :func:`~genomic_regions.RegionBased.regions`
function (see :ref:`genomic_regions`), and lists of edges. This section shows how these
edge lists can be queried using the :func:`~fanc.matrix.RegionPairsContainer.edges`
function.

In its simplest form, :func:`~fanc.matrix.RegionPairsContainer.edges` can simply be used as
a property and returns an iterator over all edges in the object. We can use this, for example,
to count the edges in the object (not the sum of weights):

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet edge iterator property
    :end-before: end snippet edge iterator property

We can also do this much more simply (and efficiently) by using the built-in :code:`len` function:

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet edge length
    :end-before: end snippet edge length

The real power of :func:`~fanc.matrix.RegionPairsContainer.edges`, however, lies in its use as
a function:

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet edge iterator function
    :end-before: end snippet edge iterator function

Note the :code:`()`. This works exactly as the above command without the function call, but now
we can introduce additional arguments. For example, to only iterate over intra-chromosomal
edges, simply do:

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet edge iterator intra
    :end-before: end snippet edge iterator intra

To only return edges connected to chromosome 19, do:

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet edge iterator chr19
    :end-before: end snippet edge iterator chr19

Importantly, this returns all edges where either source, or sink, or both connected nodes are
on chromosome 19, including, for example, inter-chromosomal edges between chromosome 18 and 19.
To only return edges within chromosome 19 (source and sink), you can combine this with
the :code:`inter_chromosomal=False` parameter. However, it is generally more efficient to use
2-dimensional selectors:

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet edge iterator 2D selectors
    :end-before: end snippet edge iterator 2D selectors

Selectors support arbitrary region definitions and human-readable abbreviations:

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet edge iterator human
    :end-before: end snippet edge iterator human

Of course, selectors also directly support :class:`~genomic_regions.GenomicRegion` objects:

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet edge iterator region
    :end-before: end snippet edge iterator region

When dealing with normalised data (such as balanced Hi-C matrices), the returned edge weights
are already normalised. You can disable the normalisation on the fly using the
:code:`norm=False` argument:

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet edge iterator norm
    :end-before: end snippet edge iterator norm

As shown above, this can be used to count the number of valid pairs in the object, for example.

.. _lazy_evaluation:

~~~~~~~~~~~~~~~
Lazy evaluation
~~~~~~~~~~~~~~~

Hi-C datasets are often very large, with hundreds of millions, even billions of valid pairs
in the matrix. The process of creating a unique :class:`~fanc.Matrix.Edge` object for every
matrix entry can thus cumulatively take a significant amount of time. For this reason, FAN-C
offers *lazy* evaluation of edge properties. When enabled, edge data is only
read when it is requested. This, for example, avoids reading from file when it is not absolutely
necessary, and saves on time during object creation and population. Edge iterators support
lazy evaluation through the :code:`lazy` argument:

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet edge iterator lazy
    :end-before: end snippet edge iterator lazy

This is significantly faster than the default :code:`lazy=False` iteration. However, using
lazy iterators it is easy to encounter confusing situations. Because data is only provided
when explicitly requested from an edge, the following code does not work as intended:

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet edge iterator wrong lazy
    :end-before: end snippet edge iterator wrong lazy

All edges in the list are identical! This is because lazy iterations use only a single
instance of the :class:`~fanc.matrix.LazyEdge` object, which simply points to the
data location in the edge in the current iteration. This pointer is replaced for the
following edge in the next iteration, but the actual object remains the same. The result
is a list of the same :class:`~fanc.matrix.LazyEdge` object pointing to the same edge
data.

Here is the correct way of obtaining data using lazy iterators:

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet edge iterator right lazy
    :end-before: end snippet edge iterator right lazy

The example accesses the edge data in the loop and stores it independently of the
:class:`~fanc.matrix.LazyEdge` object. Lazy iterators can greatly speed up your analyses,
but **be very careful** working with them!

Another useful feature of lazy iterators is that they support data modification for
native FAN-C objects. For example, you double the edge weight of each edge in the
object like this:

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet edge iterator modify lazy
    :end-before: end snippet edge iterator modify lazy

This only works for files opened in append mode (:code:`mode='a'`), and will throw an
error for files opened in read-only mode. The :func:`~fanc.matrix.LazyEdge.update` function
ensures data is written to file after modifying it.

.. warning::

    Modifying edges this way can come with unwanted consequences, and is highly discouraged.
    Potential issues include a bias vector that no longer matches, possible duplicate edges
    with unknown effect on analysis functions, messed up mappability calculations and more.
    We always recommend to create an object from scratch with the properties you want instead
    of modifying an existing one.


***************
Other functions
***************

The :func:`~fanc.matrix.RegionPairsContainer.edges` iterator is the most versatile and useful
function in the :class:`~fanc.matrix.RegionPairsContainer` interface, but by far not the only one.
We will briefly list the most useful other functions with a brief description and examples here.

~~~~~~~~~~~~~~~~~
regions_and_edges
~~~~~~~~~~~~~~~~~

The :func:`~fanc.matrix.RegionPairsContainer.regions_and_edges` function returns three objects:
a list of regions corresponding to the selected matrix rows; a list of regions corresponding to
the selected matrix columns; and an edge iterator over the matrix subset spanned by the selector.

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet regions and edges
    :end-before: end snippet regions and edges

This is simpler and saves computation time over having separate function calls to
:func:`~fanc.regions.RegionBasedWithBins.regions` and :func:`~fanc.matrix.RegionPairsContainer.edges`.

~~~~~~~~~
edge_data
~~~~~~~~~

The function :func:`~fanc.matrix.RegionPairsContainer.edge_data` iterates over only a specific
edge attribute for a selection of edges:

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet edge data
    :end-before: end snippet edge data


~~~~~~~~
mappable
~~~~~~~~

:func:`~fanc.matrix.RegionPairsContainer.mappable` returns a boolean vector where each entry
corresponds to a region in the object. If the vector entry is :code:`True`, the regions has
at least one edge connected to it, i.e. a non-zero matrix entry, otherwise the entry is :code:`False`.

.. literalinclude:: code/edges_interface_snippets.py
    :language: python
    :start-after: start snippet mappable
    :end-before: end snippet mappable

The next interface builds on the :func:`~fanc.matrix.RegionPairsContainer` structure to populate
matrices with edge weights: :ref:`matrix_interface`.