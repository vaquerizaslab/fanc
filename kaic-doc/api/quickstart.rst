=============================
Quickstart with the Kai-C API
=============================

After you have installed Kai-C (see :ref:`kaic_installation`), you can import the module
from a Python console or script:

.. code::

    import kaic

The following tutorials will assume that you have loaded the :code:`kaic` module in this
manner.

****************
Loading datasets
****************

Any analysis typically begins with loading datasets into your workspace. Kai-C tries to make this
as simple as possible with the :func:`~kaic.__init__.load` function.
If you already have processed Hi-C files, either from the :code:`kaic` command line
application (see :ref:`kaic-auto` or :ref:`kaic-modular`), or from a compatible Hi-C application
(:code:`.cool` or :code:`.mcool` from `Cooler <https://github.com/mirnylab/cooler>`_ or :code:`.hic` from
`Juicer <https://github.com/aidenlab/juicer>`_), simply load them into your workspace
using :func:`~kaic.__init__.load` - no need to specify the type of file you are loading:

.. code::

     data = kaic.load("/path/to/file.hic")

When dealing with multi-resolution Hi-C files such as :code:`.mcool` from Cooler or :code:`.hic` from
Juicer, you can load a specific resolution using the :code:`@` notation:

.. code::

     data = kaic.load("/path/to/file.mcool@25000")

:func:`~kaic.__init__.load` is not limited to Hi-C files, but also works on any other file
produced with Kai-C, such as :class:`~kaic.pairs.RegionPairs` files, :class:`~kaic.regions.Genome`,
analysis results like :class:`~kaic.architecture.comparisons.FoldChangeMatrix` and generally most
other Kai-C files.

:func:`~kaic.__init__.load` even works on most of the common file formats for genomic
datasets, such as BED, GFF, BigWig, Tabix, BEDPE and more. Try it out on your dataset of choice and
chances are :func:`~kaic.__init__.load` can handle it. And if it does not, consider raising
an issue on `Github <https://github.com/vaquerizaslab/kaic/issues>`_ to ask for support.

Internally, :func:`~kaic.__init__.load` finds a suitable class for the type of data
in the supplied file, and opens the file using that class. For example, the result of

.. code::

    hic = kaic.Hic("output/hic/binned/kaic_example_1mb.hic")

is equivalent to

.. code::

    hic = kaic.load("output/hic/binned/kaic_example_1mb.hic")

with the big advantage that you don't need to worry about remembering class names or
their location within the Kai-C module hierarchy. In both cases, the type of the
returned object is :code:`kaic.hic.Hic`:

.. code::

    # check the type fo the object
    type(hic)  # kaic.hic.Hic

Here are a few more examples:

.. code::

    cool = kaic.load("test.cool")
    type(cool)  # kaic.compatibility.cooler.CoolerHic

    juicer = kaic.load("test_juicer.hic")
    type(juicer)  # kaic.compatibility.juicer.JuicerHic

    fragments = kaic.load("hg19_chr18_19_re_fragments.bed")
    type(fragments)  # genomic_regions.regions.Bed

    bam = kaic.load("test.bam")
    type(bam)  # pysam.libcalignmentfile.AlignmentFile

    ins = kaic.load("architecture/domains/kaic_example_100kb.insulation")
    type(ins)  # kaic.architecture.domains.InsulationScores

    # and many other data types


*****************
Common interfaces
*****************

.. note::

    Before we start introducing the neat little interface functions that Kai-C provides for
    interacting with genomic regions, region pairs and matrix data, we want to briefly discuss
    the terminology used by Kai-C:

    Chromosome conformation capture data describes associations between pairs of genomic regions.
    In the literature, you will see them called "bins", or "loci", Kai-C generally uses the term
    "regions" to describe stretches of DNA defined by a set of coordinates. Similarly, associations
    between genomic regions, such as produced by Hi-C, have been called "interactions", "contacts",
    or "proximity". While these names are certainly appropriate in some situations, they can also
    be misleading by implying more than the ligation frequency measured in an experiment. Some studies
    and tools therefore use other terminology, such as "pixels". Kai-C uses a term from network
    analysis: "**edges**". In networks, edges connect "**nodes**" and can be assigned a **weight**.
    In a Hi-C "network", nodes would correspond to regions, and the (normalised) ligation frequency
    (contact intensity, interaction frequency, proximity, ...) is assigned as a weight on an "edge"
    connecting a pair of regions. Two regions in Kai-C are considered unconnected if their edge
    weight is 0. Instead of calling the regions connected by an edge "region1" and "region2" or
    something similar, they are called "source" and "sink" in Kai-C, again borrowing from network terminology.
    While the edge in Hi-C matrices has no directionality (we cannot say that region A interacts with
    region B but not vice versa), the convention in Kai-C is that the region index of the source
    region is smaller that that of the sink region.
    We believe this is a sufficiently neutral terminology for Kai-C data.


Being able to load different data types with the same command is the first convenient
feature of :func:`~kaic.__init__.load`. The second is that most objects returned by
:func:`~kaic.__init__.load` and supported by Kai-C in general share common interfaces
that unify and greatly simply handling of different datasets. We will summarise them
here to provide you with a rough overview, and the next sections will discuss each
interface in detail.

- :class:`~genomic_regions.RegionBased`: Kai-C builds heavily on the
  :code:`genomic_regions` `package <https://github.com/vaquerizaslab/genomic_regions>`_,
  which provides a unified, powerful interface for region-based data. Its functionality
  includes iterators over regions that return the same type of object
  (:class:`~genomic_regions.GenomicRegion`) regardless of data origin.
  You can query regions in a specific interval, bin scores
  associated with regions, and much more.
  This interface supports the following file types: BED, GFF, BigWig, Tabix, BEDPE, and most
  Kai-C files. Find out more in :ref:`genomic_regions`.

- :class:`~kaic.RegionPairsContainer`: This interface extends
  :class:`~genomic_regions.RegionBased` by adding properties and functions for pairs of
  genomic regions and their relationships ("edges"). You can use the powerful
  :func:`kaic.RegionpairsContainer.edges` function to iterate over all edges in an object,
  query only a subset, normalise edge weights on the fly, hide and reveal filtered
  edges, and more. This interface supports Kai-C files like :class:`~kaic.ReadPairs` and
  :class:`~kaic.Hic`, Cooler single and multi-resolution files and Juicer files.
  To find out everything about this interface, go to :ref:`edges_interface`.

- :class:`kaic.RegionMatrixContainer`: In many cases, such as Hi-C data, edges and their
  associated scores ("weights") can be represented as a matrix. This interface extends
  :class:`~kaic.RegionPairsContainer` with matrix-specific properties and functions,
  including a versatile :func:`~RegionMatrixContainer.matrix` function for retrieving
  whole-genome or subset matrices with support for on-the-fly normalisation and O/E
  transformation, masking of unmappable regions and more. This interface supports
  Kai-C files like :class:`~kaic.Hic` and  :class:`~kaic.ComparisonMatrix`, Cooler
  single and multi-resolution files and Juicer files. Go to :ref:`matrix_interface`
  for all the details.