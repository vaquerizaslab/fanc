.. _api_comparisons:


############################
Matrix and score comparisons
############################

To follow this tutorial, download the FAN-C example data, for example through our
`Keeper library <https://keeper.mpdl.mpg.de/d/147906745b634c779ed3/>`_, and set up
your Python session like this:

.. literalinclude:: code/comparisons_example_code.py
    :language: python
    :start-after: start snippet comparisons setup
    :end-before: end snippet comparisons setup

Also have a look at the command line documentation at :ref:`fanc-comparisons` for the
command-line approach to comparisons in FAN-C.

Comparisons between datasets are important to identify and highlight differences. FAN-C
has utilities to create comparison matrices and tracks by calculating their fold-change
or difference.


******************
Matrix comparisons
******************

To compare matrices in FAN-C, you can use subclasses of
:class:`~fanc.architecture.comparisons.ComparisonMatrix`. There are built-in classes for
fold-change matrices (:class:`~fanc.architecture.comparisons.FoldChangeMatrix`) and
difference matrices (:class:`~fanc.architecture.comparisons.DifferenceMatrix`), and their
usage is straightforward using the
:func:`~fanc.architecture.comparisons.ComparisonMatrix.from_matrices` function:

.. literalinclude:: code/comparisons_example_code.py
    :language: python
    :start-after: start snippet comparisons compare
    :end-before: end snippet comparisons compare

We enable ``log_cmp`` for the fold-change matrix, so comparison values are log2-transformed
and become symmetrical around 0.

Note that we could also have used ``scale=False`` in this case, to omit scaling to sequencing
depth, but if you are unsure whether your matrices are comparable without scaling it is
best to leave the setting at its default.

By default, infinite values resulting from a comparison (such as NaN from division by zero)
are omitted from the output matrix. You can keep them by disabling ``ignore_infinite``.
If you want to omit comparisons among pixels that are 0 completely, use ``ignore_zeros``.

We can show the result of the comparison using FAN-C plotting functions:

.. literalinclude:: code/comparisons_example_code.py
    :language: python
    :start-after: start snippet comparisons plot
    :end-before: end snippet comparisons plot

.. image:: images/comparisons_matrices.png

As you can see, each :class:`~fanc.architecture.comparisons.ComparisonMatrix` acts just like
a regular matrix object.


.. _fanc_compare_custom:

~~~~~~~~~~~~~~~~~~
Custom comparisons
~~~~~~~~~~~~~~~~~~

If you require a custom comparison beyond the builtin difference and fold-change, you can
easily achieve that by subclassing :class:`~fanc.architecture.comparisons.ComparisonMatrix` and
implementing a custom :func:`~fanc.architecture.comparisons.ComparisonMatrix.compare` function.
For example, if you want to create a comparison matrix where a pixels value is ``1`` if
the value in matrix 1 is larger than that in matrix 2, ``-1`` if it is smaller, and ``0`` if
they are identical, you can use:

.. literalinclude:: code/comparisons_example_code.py
    :language: python
    :start-after: start snippet comparisons custom
    :end-before: end snippet comparisons custom

Setting ``_classid`` enables loading by :func:`~fanc.load`.


*************************
Score / track comparisons
*************************

Although not necessarily Hi-C related in every case, FAN-C can also be used to compare any
kind of genomic track with scores associated with regions. File types like BED, GFF, BigWig
and many more (see :ref:`genomic_regions`) can be loaded using :func:`~fanc.load` and then
compared using :func:`~fanc.architecture.comparisons.ComparisonRegions.from_regions`. FAN-C
has built-in classes for fold-change (:class:`~fanc.architecture.comparisons.FoldChangeRegions`)
and difference (:class:`~fanc.architecture.comparisons.DifferenceRegions`):

.. literalinclude:: code/comparisons_example_code.py
    :language: python
    :start-after: start snippet regions compare
    :end-before: end snippet regions compare

We can plot it like any other region-based FAN-C object:

.. literalinclude:: code/comparisons_example_code.py
    :language: python
    :start-after: start snippet regions plot
    :end-before: end snippet regions plot

.. image:: images/comparisons_regions.png

This outputs a :class:`~fanc.regions.RegionsTable`, but you can export to file using the
:func:`~fanc.regions.RegionsTable.to_bed`, :func:`~fanc.regions.RegionsTable.to_gff`, and
:func:`~fanc.regions.RegionsTable.to_bigwig` functions.

Use ``log=True`` to log2-transform comparison values after the comparison, for example for
fold-change comparisons. You can change the attribute that is being compared using the
``attribute`` parameter, which defaults to ``"score"``. Similarly, if you want the comparison
to be saved under a different attribute name, you can specify that using ``score_field``.

~~~~~~~~~~~~~~~~~~
Custom comparisons
~~~~~~~~~~~~~~~~~~

For custom region-based comparisons, you can subclass
:func:`~fanc.architecture.comparisons.ComparisonRegions` and override
:func:`~fanc.architecture.comparisons.ComparisonRegions.compare` in the same way you would with
:class:`~fanc.architecture.comparisons.ComparisonMatrix` (see :ref:`above <fanc_compare_custom>`).


*********************************
Parameter-based score comparisons
*********************************

In :ref:`api_tads` we demonstrated how you can use
:class:`~fanc.architecture.domains.RegionScoreParameterTable` objects to store parameter-based scores,
such as window sizes in :class:`~fanc.architecture.domains.InsulationScores`. Such
:class:`~fanc.architecture.domains.RegionScoreParameterTable` objects can also be compared
using FAN-C - in this case, a separate comparison is run for each parameter. The result is a
:class:`~fanc.architecture.comparisons.ComparisonScores` object, which is based on
:class:`~fanc.architecture.domains.RegionScoreParameterTable` and can be used as such. The
comparison is done with :func:`~fanc.architecture.comparisons.ComparisonScores.from_scores`.

.. literalinclude:: code/comparisons_example_code.py
    :language: python
    :start-after: start snippet scores compare
    :end-before: end snippet scores compare

We can plot it like any other parameter-based FAN-C object:

.. literalinclude:: code/comparisons_example_code.py
    :language: python
    :start-after: start snippet scores plot
    :end-before: end snippet scores plot

.. image:: images/comparisons_scores.png
