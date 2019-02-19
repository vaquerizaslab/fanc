.. _kaic-helpers:

####################
Helper tools in kaic
####################

Kai-C provides little helper tools that make working with Hi-C and associated data
somewhat easier. These are not strictly necessary for matrix generation and analysis,
but can often speed your analysis up or simply make it a little more convenient.


==========================================
kaic dump: export Hic objects to text file
==========================================

You can easily export Kai-C Hic objects to a txt file using ``kaic dump``.

.. argparse::
   :module: kaic.commands.kaic_commands
   :func: dump_parser
   :prog: kaic dump
   :nodescription:

If you only pass the Hic object the ``kaic dump``, it will write all Hi-C contacts to
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
kaic subset: create Hic objects by subsetting
=============================================

It is sometimes useful to work with smaller Hi-C objects, for example for speed reasons
or to focus the analysis on a particular genomic region of interest. The ``kaic subset``
command makes it possible to create a Hic object that only contains regions and contacts
between a user-specified genomic regions from an existing Hic object.

.. argparse::
   :module: kaic.commands.kaic_commands
   :func: subset_parser
   :prog: kaic subset
   :nodescription:


=======================================
kaic downsample: downsample Hic objects
=======================================

Often Hi-C matrices have differing numbers of valid pairs, which can be a confounding factor
in many analyses. Differences can stem from varying sequencing depths, different library
qualities, or other experimental and computational factors. ``kaic downsample`` is a utility
that downsamples Hic objects to a specific number of valid pairs.

.. argparse::
   :module: kaic.commands.kaic_commands
   :func: downsample_parser
   :prog: kaic downsample
   :nodescription:

By default, the sampling is done without replacement. This requires a fairly large amount
of system memory. If you are having trouble with memory usage, use sampling with
replacement (``--with-replacement``). Note that the samplig is done on uncorrected matrix
values, so you may want to apply matrix balancing using ``kaic hic -k`` afterwards.


==========================================
kaic fragments: in silico genome digestion
==========================================

The ``kaic pairs`` and ``kaic auto`` commands accept FASTA files as ``--genome`` argument,
and ``kaic`` conveniently calculates the restriction fragments for you using the
restriction enzyme name specified with ``--restriction-enzyme``. However, the in silico
digestion can be time-consuming, and if you are processing multiple similar Hi-C libraries,
you can use the ``kaic fragments`` utility to generate restriction fragments up front,
and use the resulting BED file as input for the ``--genome`` argument.

If you supply an integer as the second positional argument instead of a restriction enzyme
name, ``kaic fragments`` will perform binning rather than in silico digestion and return
a BED file with equally sized regions.

.. argparse::
   :module: kaic.commands.kaic_commands
   :func: fragments_parser
   :prog: kaic fragments
   :nodescription:


=====================================
kaic sort-sam: sort SAM files by name
=====================================

The ``kaic pairs`` command expects SAM/BAM files as input that have been sorted by name
(``kaic auto`` automatically sorts files). You can use ``samtools sort -n`` to sort files,
but ``kaic sam-sort`` will also do the sorting for you. it automatically chooses the fastest
sorting implementation available and also provides the option to work in a temporary folder,
which can speed the sorting up if you are working on a network volume.


.. argparse::
   :module: kaic.commands.kaic_commands
   :func: sort_sam_parser
   :prog: kaic sort-sam
   :nodescription:
