.. _fanc-helpers:

####################
Helper tools in fanc
####################

FAN-C provides little helper tools that make working with Hi-C and associated data
somewhat easier. These are not strictly necessary for matrix generation and analysis,
but can often speed your analysis up or simply make it a little more convenient.


.. _fanc_from_txt:

========================================
fanc from-txt: import Hic from text file
========================================

You can easily import Hi-C matrices from a compatible text file format, such as that from
HiC-Pro, with ``fanc to-txt``.

.. argparse::
   :module: fanc.commands.fanc_commands
   :func: from_txt_parser
   :prog: fanc from-txt
   :nodescription:
   :nodefault:

The command requires two input files:

1. A sparse matrix with the tab-separated format ``<bin1><tab><bin2><tab><weight>``:

  .. code::

     1	1	40.385642
     1	828	5.272852
     1	1264	5.205258
     ...

2. A regions file in BED format ``<chromosome><tab><start><tab><end>[<tab><bin ID>]``:

  .. code::

     chr1	0	1000000	1
     chr1	1000000	2000000	2
     chr1	2000000	3000000	3
     chr1	3000000	4000000	4
     ...

  The ``<bin ID>`` field is optional, but if provided it must correspond to the bins used
  in the matrix file. If not provided, bin indices will be 0-based!

The FAN-C example data contains some HiC-Pro example files that you can try this out on:

.. code::

   fanc from-txt hicpro/dixon_2M_1000000_iced.matrix hicpro/dixon_2M_1000000_abs.bed hicpro/dixon_2M_1000000_iced.hic


==========================================
fanc dump: export Hic objects to text file
==========================================

You can easily export FAN-C Hic objects to a txt file using ``fanc dump``.

.. argparse::
   :module: fanc.commands.fanc_commands
   :func: dump_parser
   :prog: fanc dump
   :nodescription:
   :nodefault:

If you only pass the Hic object the ``fanc dump``, it will write all Hi-C contacts to
the command line in a tab-delimited format with the columns: chromosome1, start1, end1,
chromosome2, start2, end2, weight (number of contacts). If you add a file path as
second argument, the data will be written to that file. If you instead pass the Hic file
and two output files, the first output file will have the matrix entries in sparse notation,
and the second file will have the Hic regions/bins. You can use ``-S`` to export a full
matrix instead of a sparse one, but be warned that these can be extremely large.

If you are only interested in a specific sub-matrix, use the ``-s`` or ``--subset`` argument
of the for <chromosome>:[<start>-<end>] to export all contacts made by this particular
region across the whole genome. Use <chr>[:<start>-<end>]--<chr>[:<start><end>] to export
all contacts made between two regions. E.g. use chr1--chr1 to export the chromosome 1
sub-matrix.


=============================================
fanc subset: create Hic objects by subsetting
=============================================

It is sometimes useful to work with smaller Hi-C objects, for example for speed reasons
or to focus the analysis on a particular genomic region of interest. The ``fanc subset``
command makes it possible to create a Hic object that only contains regions and contacts
between a user-specified genomic regions from an existing Hic object.

.. argparse::
   :module: fanc.commands.fanc_commands
   :func: subset_parser
   :prog: fanc subset
   :nodescription:
   :nodefault:


=======================================
fanc downsample: downsample Hic objects
=======================================

Often Hi-C matrices have differing numbers of valid pairs, which can be a confounding factor
in many analyses. Differences can stem from varying sequencing depths, different library
qualities, or other experimental and computational factors. ``fanc downsample`` is a utility
that downsamples Hic objects to a specific number of valid pairs.

.. argparse::
   :module: fanc.commands.fanc_commands
   :func: downsample_parser
   :prog: fanc downsample
   :nodescription:
   :nodefault:

By default, the sampling is done without replacement. This requires a fairly large amount
of system memory. If you are having trouble with memory usage, use sampling with
replacement (``--with-replacement``).

.. note::

    Sampling is done on uncorrected matrix values, so you may want to apply matrix
    balancing using ``fanc hic -k`` afterwards.


==========================================
fanc fragments: in silico genome digestion
==========================================

The ``fanc pairs`` and ``fanc auto`` commands accept FASTA files as ``--genome`` argument,
and ``fanc`` conveniently calculates the restriction fragments for you using the
restriction enzyme name specified with ``--restriction-enzyme``. However, the in silico
digestion can be time-consuming, and if you are processing multiple similar Hi-C libraries,
you can use the ``fanc fragments`` utility to generate restriction fragments up front,
and use the resulting BED file as input for the ``--genome`` argument.

If you supply an integer as the second positional argument instead of a restriction enzyme
name, ``fanc fragments`` will perform binning rather than in silico digestion and return
a BED file with equally sized regions.

.. argparse::
   :module: fanc.commands.fanc_commands
   :func: fragments_parser
   :prog: fanc fragments
   :nodescription:
   :nodefault:


=====================================
fanc sort-sam: sort SAM files by name
=====================================

The ``fanc pairs`` command expects SAM/BAM files as input that have been sorted by name
(``fanc auto`` automatically sorts files). You can use ``samtools sort -n`` to sort files,
but ``fanc sam-sort`` will also do the sorting for you. it automatically chooses the fastest
sorting implementation available and also provides the option to work in a temporary folder,
which can speed the sorting up if you are working on a network volume.


.. argparse::
   :module: fanc.commands.fanc_commands
   :func: sort_sam_parser
   :prog: fanc sort-sam
   :nodescription:
   :nodefault:
