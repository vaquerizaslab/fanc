.. _kaic-comparisons:

############################
Matrix and score comparisons
############################

Kai-C provides a central utility named ``kaic compare`` to compare Hi-C matrices
and measures derived from them between two conditions, such as different cell types,
treatments, etc.

.. argparse::
   :module: kaic.commands.kaic_commands
   :func: compare_parser
   :prog: kaic compare
   :nodescription:


****************
Compare matrices
****************

When you provide ``kaic compare`` with two matrix files (e.g. Hic), these need to have
the same binning (more specifically, they must have identical regions). The comparison is
then done on a per-pixel basis. By default, the fold-change between pixels in the two
matrices is calculated, but you can also calculate the difference using ``-c difference``
(you can also use "fc" and "diff" as short-hand).

The output is a Kai-C matrix object (that can be used like a Hic object in most cases)
where entries are the fold-change (or difference) of pixels in both matrices. The
comparison is matrix1/matrix2 (or matrix1 - matrix2). For fold-change comparisons, it
might be useful to log-transform pixel values using ``-l``. Also, if pixel entries are
zero in one matrix, the fold-change could return NaN values. You can simply ignore
comparisons where one or both pixels are zero using ``-Z``. Or you can ignore infinite
values resulting from a comparison using ``-I``.

For a fair comparison, matrix entries are scaled to the total number of valid pairs on
each chromosome. This corrects for different coverage of the two matrices. If you have
chosen a normalisation that takes care of this directly, you can save some computation
time by specifying ``-S`` to omit the calculation of a scaling factor.


**************
Compare scores
**************

You can also use ``kaic compare`` to compare scores in any compatible region-based file,
such as BED, GFF, BigWig, or any region-based Kai-C output (e.g. insulation scores).
By default, the output will be a Kai-C :class:`~kaic.RegionsTable` object, but you can
change this behavior using the ``-o`` option, which you can use to output BED, GFF, or
BigWig files.


*************************
Compare multi-score files
*************************

The output of ``kaic insulation`` and ``kaic directionality`` for multiple window sizes
is a Kai-C multi-score file, which contains the insulation/directionality scores for
multiple windows sizes at once. You can directly compare all scores against each other
by providing two multi-score files as input to ``kaic compare`` - every window size will
be compared independently and the comparisons output in yet another multi-score file.
