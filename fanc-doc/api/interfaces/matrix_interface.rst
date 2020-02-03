.. _matrix_interface:

=====================
RegionMatrixContainer
=====================

This interface simplifies and unifies working with matrix data in the context of genomic
region pairs, such as you would find in a Hi-C matrix. It builds on the
:class:`~fanc.matrix.RegionPairsContainer` (see previous section :ref:`edges_interface`),
which dealt with scores and other attributes between genomic regions, and extends it by
adding functions for representing the scores in a numeric matrix.

After loading a dataset using :func:`~fanc.tools.load.load`, you can check for
support of the :class:`~fanc.matrix.RegionMatrixContainer` interface with:

.. literalinclude:: code/matrix_interface_snippets.py
    :language: python
    :start-after: start snippet check
    :end-before: end snippet check

The current list of FAN-C classes supporting :class:`~fanc.matrix.RegionMatrixContainer` is:
:class:`~fanc.compatibility.cooler.CoolerHic`,
:class:`~fanc.compatibility.juicer.JuicerHic`,
:class:`~fanc.hic.Hic`,
:class:`~fanc.architecture.compartments.ABCompartmentMatrix`,
:class:`~fanc.architecture.comparisons.DifferenceMatrix`,
:class:`~fanc.architecture.comparisons.FoldChangeMatrix`,
:class:`~fanc.peaks.PeakInfo`,
and
:class:`~fanc.peaks.RaoPeakInfo`.


*******************
The matrix function
*******************

To obtain the whole-genome, normalised matrix from an object, use the
:func:`~fanc.matrix.RegionMatrix.matrix` function:

.. literalinclude:: code/matrix_interface_snippets.py
    :language: python
    :start-after: start snippet matrix whole genome
    :end-before: end snippet matrix whole genome

Of course, the :func:`~fanc.matrix.RegionMatrixContainer.matrix` function supports matrix subsets:

.. literalinclude:: code/matrix_interface_snippets.py
    :language: python
    :start-after: start snippet matrix subset
    :end-before: end snippet matrix subset

When using tuples as keys, the first entry will select the rows, and the second entry
the columns of the matrix:

.. literalinclude:: code/matrix_interface_snippets.py
    :language: python
    :start-after: start snippet matrix row col
    :end-before: end snippet matrix row col

The returned object is of type :class:`~fanc.matrix.RegionMatrix`, which is a subclass
of Numpy's masked :class:`~numpy.ma.array` with added perks for genomic region handling.

A :class:`~fanc.matrix.RegionMatrix` can be used like any other numpy matrix,
for example calculating marginals by summing up values in rows or columns:

.. literalinclude:: code/matrix_interface_snippets.py
    :language: python
    :start-after: start snippet matrix marginals
    :end-before: end snippet matrix marginals

(this Hi-C object is normalised on a per-chromosome basis, so each marginal will be
close to 1)

Rows and columns in a matrix can be masked, i.e. their entries are considered invalid and
are ignored for most downstream analysis to prevent artifacts. By default, FAN-C masks
regions that have no edges (after filtering), typically due to mappability issues.
You can turn off masking using the :code:`mask=False` parameter:

.. literalinclude:: code/matrix_interface_snippets.py
    :language: python
    :start-after: start snippet matrix no mask
    :end-before: end snippet matrix no mask

However, we recommend working with masked matrices to ensure no unwanted edges are
part of your analyses.

:class:`~fanc.matrix.RegionMatrix` objects also keep track of the regions corresponding to
columns and rows of a matrix:

.. literalinclude:: code/matrix_interface_snippets.py
    :language: python
    :start-after: start snippet matrix regions
    :end-before: end snippet matrix regions

You can subset a :class:`~fanc.matrix.RegionMatrix` using indexes or region intervals:

.. literalinclude:: code/matrix_interface_snippets.py
    :language: python
    :start-after: start snippet matrix region matrix subset
    :end-before: end snippet matrix region matrix subset

Note that region interval definitions are always interpreted as 1-based, inclusive, and any
overlapping region is returned (in the above example the region :code:`chr19:150001-200000`
has a 1 base overlap with the requested interval).

:func:`~fanc.matrix.RegionMatrixContainer.matrix` supports all arguments also available for
:func:`~fanc.matrix.RegionPairsContainer.edges`, but it is not necessary to use lazy loading.
You can, for example, output an uncorrected matrix with

.. literalinclude:: code/matrix_interface_snippets.py
    :language: python
    :start-after: start snippet matrix no norm
    :end-before: end snippet matrix no norm

In addition, there are several parameters specific to
:func:`~fanc.matrix.RegionMatrixContainer.matrix`. Most notably, you can use the
:code:`oe=True` parameter to return an observed/expected (O/E) matrix:

.. literalinclude:: code/matrix_interface_snippets.py
    :language: python
    :start-after: start snippet matrix oe
    :end-before: end snippet matrix oe

Internally, :code:`oe=True` uses
:class:`~fanc.matrix.RegionMatrixContainer.expected_values` to calculate the expected
(average) weight of all edges connecting regions at a certain distance. The matrix
is then simply divided by the expected matrix. You may want to log2-transform the
matrix for a symmetric scale of values:

.. literalinclude:: code/matrix_interface_snippets.py
    :language: python
    :start-after: start snippet matrix log oe
    :end-before: end snippet matrix log oe

