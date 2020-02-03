.. _genomic_regions:

===========
RegionBased
===========

FAN-C builds extensively on the :code:`genomic_regions` package, which provides a unified
interface for most types of region-based genomic data. We highly recommend
`reading the documentation <https://vaquerizaslab.github.io/genomic_regions/>`_
of that package before going into the details of FAN-C, as many of the concepts discussed
therein are central to the handling of data in FAN-C.

You can check whether a FAN-C object supports the :class:`~genomic_regions.RegionBased`
interface with

.. code::

   import genomic_regions as gr
   isinstance(o, gr.RegionBased)  # True for objects supporting the regions interface

The current list of FAN-C objects supporting the :class:`~genomic_regions.RegionBased`
interface is:
:class:`~fanc.architecture.domains.InsulationScore`,
:class:`~fanc.architecture.domains.DirectionalityIndex`,
:class:`~fanc.architecture.domains.Boundaries`,
:class:`~fanc.architecture.domains.InsulationScores`,
:class:`~fanc.architecture.domains.DirectionalityIndexes`,
:class:`~fanc.architecture.comparisons.FoldChangeScores`,
:class:`~fanc.architecture.comparisons.DifferenceScores`,
:class:`~fanc.architecture.comparisons.DifferenceRegions`,
:class:`~fanc.architecture.comparisons.FoldChangeRegions`,
:class:`~fanc.compatibility.cooler.CoolerHic`,
:class:`~fanc.compatibility.juicer.JuicerHic`,
:class:`~fanc.hic.Hic`,
:class:`~fanc.architecture.compartments.ABCompartmentMatrix`,
:class:`~fanc.architecture.comparisons.DifferenceMatrix`,
:class:`~fanc.architecture.comparisons.FoldChangeMatrix`,
:class:`~fanc.peaks.PeakInfo`,
and
:class:`~fanc.peaks.RaoPeakInfo`.

Any object built on that foundation supports, for example, region iterators:

.. code::

   for region in hic.regions:
       print(region)
       print(region.chromosome, region.start, region.end, region.strand)
       print(region.is_forward)
       print(region.center)
       # ...

Range queries:

.. code::

   for region in hic.regions('chr1:3mb-12mb'):
       print(region.chromosome)  # chr1
       # ...

and many more convenient features. The object type returned by all of those queries
is :class:`~genomic_regions.GenomicRegion`, which has many convenient functions to
deal with region properties and operations.

.. code:: python

    len(region)  # returns the size of the region in base pairs
    region.center  # returns the base (or fraction of base) at the center of the region
    region.five_prime  # returns the starting base at the 5' end of the region
    region.three_prime  # returns the starting base at the 3' end of the region
    region.is_forward()  # True if strand is '+' or '+1'
    region.is_reverse()  # True if strand is '-' or '-1'
    region.attributes  # return all attribute names in this region object
    region.copy()  # return a shallow copy of this region
    region.to_string()  # return a region identifier string describing the region

    region = gr.as_region('chr12:12.5Mb-18Mb')
    region.overlaps('chr12:11Mb-13Mb')  # True
    region.overlaps('chr12:11Mb-11.5Mb')  # False
    region.overlaps('chr1:11Mb-13Mb')  # False

Refer to the
`genomic_regions documentation <https://github.com/vaquerizaslab/genomic_regions>`_ for
all the details.

Similarly to the :code:`regions` interface for handling collections of genomic regions,
FAN-C implements interfaces for working with pairs of genomic regions (:code:`edges`)
and matrix operations (:code:`matrix`). These work in exactly the same way for FAN-C,
Cooler, and Juicer files. Hence, all of these are directly compatible with FAN-C architectural
functions such as the insulation score or AB compartment analyses, ...

These interfaces will be introduced in the following sections, starting with :ref:`edges_interface`.