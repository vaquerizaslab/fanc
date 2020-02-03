.. _api_oe:

#############################
Expected and O/E calculations
#############################

The following steps assume that you ran the ``fanc auto`` command in :ref:`example-fanc-auto`.
Additionally, we set up the Python session like this:

.. literalinclude:: code/oe_example_code.py
    :language: python
    :start-after: start snippet oe setup
    :end-before: end snippet oe setup

:class:`~fanc.matrix.RegionMatrixContainer` objects (see :ref:`here <matrix_interface>`) have a builtin
function to calculate expected values from existing matrix data called
:func:`~fanc.matrix.RegionMatrixContainer.expected_values`. This function calculates and returns
intra-chromosomal, intra-chromosomal per chromosome, and inter-chromosomal expected values.

.. literalinclude:: code/oe_example_code.py
    :language: bash
    :start-after: start snippet oe basic
    :end-before: end snippet oe basic

Here, ``intra_expected`` is a list of average (/expected) contact values, where the position of
the value in the list corresponds to the separation between genomic regions in bins.
``intra_expected_chromosome`` is a dictionary with chromosome names as keys, and an expected
value list as value calculated on a per-chromosome basis. ``inter_expected`` is a single, average
inter-chromosomal contact value.

When running the :func:`~fanc.matrix.RegionMatrixContainer.expected_values` function on a matrix
opened in read-only mode, you will receive a warning that the values cannot be save to the object.
This means, the expected values will be recalculated every time they are needed in downstream
analyses. You can ensure the values are saved to the object by opening the matrix in "append" mode:

.. literalinclude:: code/oe_example_code.py
    :language: python
    :start-after: start snippet oe append
    :end-before: end snippet oe append

The expected values are typically plotted on a log-log scale, as illustrated here using chromosome 19:

.. literalinclude:: code/oe_example_code.py
    :language: python
    :start-after: start snippet oe ddplot
    :end-before: end snippet oe ddplot

.. image:: images/oe_500kb.png

Note: as Hi-C matrices are normalised on a per-chromosome basis in Kai-C by default, it would be misleading
to plot the overall normalised intra-chromosomal expected values, or to use them for downstream analysis.
We can, however, also calculate the unnormalised expected values easily enough:

.. literalinclude:: code/oe_example_code.py
    :language: python
    :start-after: start snippet oe nonorm
    :end-before: end snippet oe nonorm

.. image:: images/oe_500kb_nonorm.png

Expected values rarely need to be calculated explicitly in Kai-C analysis functions, but will be calculated
(or retrieved) on demand whenever necessary. To obtain observed/expected matrices, for example, please
refer to :ref:`matrix_interface`.