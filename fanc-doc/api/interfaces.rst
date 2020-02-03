.. _common_interfaces:

#################
Common interfaces
#################

.. note::

    Before we start introducing the neat little interface functions that FAN-C provides for
    interacting with genomic regions, region pairs and matrix data, we want to briefly discuss
    the terminology used by FAN-C:

    Chromosome conformation capture data describes associations between pairs of genomic regions.
    In the literature, you will see them called "bins", or "loci", FAN-C generally uses the term
    "regions" to describe stretches of DNA defined by a set of coordinates. Similarly, associations
    between genomic regions, such as produced by Hi-C, have been called "interactions", "contacts",
    or "proximity". While these names are certainly appropriate in some situations, they can also
    be misleading by implying more than the ligation frequency measured in an experiment. Some studies
    and tools therefore use other terminology, such as "pixels". FAN-C uses a term from network
    analysis: "**edges**". In networks, edges connect "**nodes**" and can be assigned a **weight**.
    In a Hi-C "network", nodes would correspond to regions, and the (normalised) ligation frequency
    (contact intensity, interaction frequency, proximity, ...) is assigned as a weight on an "edge"
    connecting a pair of regions. Two regions in FAN-C are considered unconnected if their edge
    weight is 0. Instead of calling the regions connected by an edge "region1" and "region2" or
    something similar, they are called "source" and "sink" in FAN-C, again borrowing from network terminology.
    While the edge in Hi-C matrices has no directionality (we cannot say that region A interacts with
    region B but not vice versa), the convention in FAN-C is that the region index of the source
    region is smaller that that of the sink region.
    We believe this is a sufficiently neutral terminology for FAN-C data.


Being able to load different data types with the same command is the first convenient
feature of :func:`~fanc.__init__.load`. The second is that most objects returned by
:func:`~fanc.__init__.load` and supported by FAN-C in general share common interfaces
that unify and greatly simply handling of different datasets. We will summarise them
here to provide you with a rough overview, and the next sections will discuss each
interface in detail.

- :class:`~genomic_regions.RegionBased`: FAN-C builds heavily on the
  :code:`genomic_regions` `package <https://github.com/vaquerizaslab/genomic_regions>`_,
  which provides a unified, powerful interface for region-based data. Its functionality
  includes iterators over regions that return the same type of object
  (:class:`~genomic_regions.GenomicRegion`) regardless of data origin.
  You can query regions in a specific interval, bin scores
  associated with regions, and much more.
  This interface supports the following file types: BED, GFF, BigWig, Tabix, BEDPE, and most
  FAN-C files. Find out more in :ref:`genomic_regions`.

- :class:`~fanc.RegionPairsContainer`: This interface extends
  :class:`~genomic_regions.RegionBased` by adding properties and functions for pairs of
  genomic regions and their relationships ("edges"). You can use the powerful
  :func:`fanc.RegionpairsContainer.edges` function to iterate over all edges in an object,
  query only a subset, normalise edge weights on the fly, hide and reveal filtered
  edges, and more. This interface supports FAN-C files like :class:`~fanc.ReadPairs` and
  :class:`~fanc.Hic`, Cooler single and multi-resolution files and Juicer files.
  To find out everything about this interface, go to :ref:`edges_interface`.

- :class:`fanc.RegionMatrixContainer`: In many cases, such as Hi-C data, edges and their
  associated scores ("weights") can be represented as a matrix. This interface extends
  :class:`~fanc.RegionPairsContainer` with matrix-specific properties and functions,
  including a versatile :func:`~RegionMatrixContainer.matrix` function for retrieving
  whole-genome or subset matrices with support for on-the-fly normalisation and O/E
  transformation, masking of unmappable regions and more. This interface supports
  FAN-C files like :class:`~fanc.Hic` and  :class:`~fanc.ComparisonMatrix`, Cooler
  single and multi-resolution files and Juicer files. Go to :ref:`matrix_interface`
  for all the details.


.. toctree::
   :maxdepth: 1

   Genomic regions (RegionBased interface) <interfaces/regions_interface>
   Edges (contacts, pixels, ...) (RegionPairsContainer interface) <interfaces/edges_interface>
   Matrices (RegionMatrixContainer interface) <interfaces/matrix_interface>
