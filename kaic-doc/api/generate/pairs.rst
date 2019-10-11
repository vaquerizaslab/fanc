.. _api_pairs:

===================================
Generating and filtering read pairs
===================================

Once we have mapped and sorted BAM files (see :ref:`api_map`) we can load the paired reads
into a :class:`~kaic.pairs.ReadPairs` object.

*******************
Obtaining fragments
*******************

:class:`~kaic.pairs.ReadPairs` map reads to restriction fragments, so first of all we are
going to need a list of those. The :func:`kaic.regions.genome_regions` function has been
made for that purpose. If you provide it with a FASTA file and a restriction enzyme name,
it will build the list of fragments for you:

.. literalinclude:: code/generate_example_code.py
    :language: python
    :start-after: start snippet genome
    :end-before: end snippet genome

You can also give it a :class:`~genomic_regions.RegionBased` compatible file with fragments
directly (BED, GFF and similar work). Note that the fragments should cover the entire genome,
so don't just use a list of restriction sites!

The returned object is of class :class:`~genomic_regions.RegionBased`. See
:ref:`the main article <genomic_regions>`
if you want to learn how to interact with objects of that class.


************************************
Setting up filters at the read level
************************************

Filters can be used to remove problematic reads, such as those with poor quality
(:class:`~kaic.pairs.QualityFilter`) or with multiple mapping locations throughout
the genome (:class:`~kaic.pairs.UniquenessFilter`). These filters (of class
:class:`~kaic.pairs.ReadFilter`) have to be added already at the read pair import step,
as the necessary information is otherwise lost in later processing stages:

.. literalinclude:: code/generate_example_code.py
    :language: python
    :start-after: start snippet read filters
    :end-before: end snippet read filters

.. warning::

    Both filters are specific to reads mapped with ``bowtie2``, and there are versions of
    these filters specific for BWA-mapped reads (:class:`~kaic.pairs.BwaMemQualityFilter`
    and :class:`~kaic.pairs.BwaMemUniquenessFilter`). Make sure you use the appropriate
    filter for your chosen mapper, as the quality and uniqueness definitions used by each
    application differ substantially!

The :class:`~kaic.pairs.QualityFilter` removes all read pairs where one or both reads have
a MAPQ value lower than ``30`` (in this case). The :class:`~kaic.pairs.UniquenessFilter`
has two settings. If ``strict=True`` it will remove any read pair where either of the two
reads has the ``XS`` tag, which indicates the presence of a secondary alignment. If
``strict=False``, it further requires that the secondary alignment score is higher or equal
to the primary alignment score (``AS``). The latter setting is therefore more permissive.

Now we have set up the necessary filters, we can start importing read pairs and match them
to restriction fragments.


**********************************
Importing and accessing read pairs
**********************************

For importing read pairs from SAM files we can use :func:`~kaic.pairs.generate_pairs_split`.
It requires the two paired-end SAM (or BAM) files and fragment definitions. In addition, we
provide it with the read filters we created above, so that only reads meeting our
requirements are imported into the object. If we do not pass a file name with ``output_file``,
the resulting object is created in memory, which is useful for testing, but highly inadvisable
for actual Hi-C libraries. Here, we request four threads to be used for the loading in parallel:

.. literalinclude:: code/generate_example_code.py
    :language: python
    :start-after: start snippet import pairs
    :end-before: end snippet import pairs

The resulting object is of the type :class:`~kaic.pairs.ReadPairs`, which implements
:ref:`edges_interface`. It contains all valid read pairs matched to their respective
restriction fragments.

Let's look at some basic information:

.. literalinclude:: code/generate_example_code.py
    :language: python
    :start-after: start snippet pairs info
    :end-before: end snippet pairs info

As you can see, :class:`~kaic.pairs.ReadPairs` is a container for
:class:`~kaic.pairs.FragmentReadPair` objects. Each object contains information about its
two associated reads in the :py:attr:`~kaic.pairs.FragmentReadPair.left` and
:py:attr:`~kaic.pairs.FragmentReadPair.right` properties, respectively, which return a
:class:`~kaic.pairs.FragmentRead` object.

:class:`~kaic.pairs.FragmentRead` objects contain information about the read's mapping
location (``pair.right.position`` and ``pair.right.strand``), and the fragment it maps to
(``pair.right.fragment``). For convenience, you can also directly calculate the distance of
the read to the end of the fragment (which is the nearest restriction site) using
``pair.right.re_distance()``. The ``pair.right.fragment`` attribute is of type
:class:`~genomic_regions.regions.GenomicRegion`, which has lots of useful properties which
you can find in the `genomic_regions documentation <https://github.com/vaquerizaslab/genomic_regions>`_.
Incidentally, the above pair is a self-ligated" fragment, which will be filtered out
in post-processing.

Typically you would not access each pair individually, but iterate over all pairs, or pair
subsets with the :func:`kaic.pairs.ReadPairs.pairs` function:

.. literalinclude:: code/generate_example_code.py
    :language: python
    :start-after: start snippet pairs iter
    :end-before: end snippet pairs iter

The :func:`kaic.pairs.ReadPairs.pairs` function supports :ref:`lazy evaluation <lazy_evaluation>`
with the ``lazy=True`` argument, but make sure to read the caveats in the main article!


**************************
Read/fragment pair filters
**************************

In addition to the quality and uniqueness filters, which had to be used during read pair import,
Kai-C offers a number of :class:`~kaic.pairs.FragmentReadPairFilter` to remove problematic pairs.
First of all, however, we'll introduce one little function that you can use to reset all filters
(and restore the original pairs) in case you found your filtering t be too harsh or if you want
a fresh start, which is called :func:`~kaic.pairs.ReadPairs.reset_filters`:

.. literalinclude:: code/generate_example_code.py
    :language: python
    :start-after: start snippet pairs reset
    :end-before: end snippet pairs reset

We will illustrate how filters are used with the :class:`~kaic.pairs.SelfLigationFilter`.

.. literalinclude:: code/generate_example_code.py
    :language: python
    :start-after: start snippet pairs filter example
    :end-before: end snippet pairs filter example

In general, we first create a filter instance. The ``mask`` parameter allows us to specify a
filter name. We can also pass it a more elaborate :class:`~kaic.general.Mask` object with
properties for name and description (see commented code). Name and description of the mask will
be stored in the ``pairs`` object upon filtering. The filter is run simply by calling
:func:`~kaic.pairs.ReadPairs.filter` with the corresponding filter instance.

If you want to run multiple filters, it is more efficient to "queue" filters and then execute
them all in one go:

.. literalinclude:: code/generate_example_code.py
    :language: python
    :start-after: start snippet pairs filter queue
    :end-before: end snippet pairs filter queue

Once the :func:`~kaic.pairs.ReadPairs.filter` command has completed, the filtered pairs will
appear to have been deleted from the object.

.. literalinclude:: code/generate_example_code.py
    :language: python
    :start-after: start snippet pairs filter masked
    :end-before: end snippet pairs filter masked

In reality, filtering only masks ("hides") edges, so we have the opportunity to reset filters
or selectively disable them when iterating:

.. literalinclude:: code/generate_example_code.py
    :language: python
    :start-after: start snippet pairs filter exclude
    :end-before: end snippet pairs filter exclude

To obtain a dictionary with the filter statistics, use :func:`~kaic.pairs.ReadPairs.filter_statistics`

.. literalinclude:: code/generate_example_code.py
    :language: python
    :start-after: start snippet pairs filter stats
    :end-before: end snippet pairs filter stats


~~~~~~~~~~~~
Filter types
~~~~~~~~~~~~

Here are the filters available in Kai-C for read/fragment pairs:

- :class:`~kaic.pairs.SelfLigationFilter` removes all pairs where both reads map to the same fragment
- :class:`~kaic.pairs.ReDistanceFilter` can filter for an expected DNA (not restriction) fragment
  size in Hi-C libraries that result from fragmentation (typically sonication). It sums up the
  distance of both reads to their nearest restriction site. If that sum exceeds the cutoff set at
  filter creation, it will be marked as invalid
- :class:`~kaic.pairs.PCRDuplicateFilter` will find pairs mapping to the exact same genomic locations
  (both reads), and will only retain one copy of the exact duplicates.
- :class:`~kaic.pairs.InwardPairsFilter` and :class:`~kaic.pairs.OutwardPairsFilter` are removing
  pairs where both reads map within a specific distance on the same chromosome and are oriented
  towards or away from each other (in terms of strand), respectively. This is designed to remove
  ligation products that have likely arisen from uncut restriction sites or that are unligated.

Next, we will convert the filtered pairs to a :class:`~kaic.hic.Hic` object.
