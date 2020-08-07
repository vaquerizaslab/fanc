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

If you want to try the tutorial with an equivalent Cooler file, load the Hi-C file like
this instead:

.. literalinclude:: code/oe_example_code.py
    :language: python
    :start-after: start snippet alternative cooler
    :end-before: end snippet alternative cooler

or like this if you want to work with a Juicer file built from the same data:

.. literalinclude:: code/oe_example_code.py
    :language: python
    :start-after: start snippet alternative juicer
    :end-before: end snippet alternative juicer

Note that there may be minor differences in the results due to the "zooming" and balancing
applied by the different tools.

:class:`~fanc.matrix.RegionMatrixContainer` objects (see :ref:`here <matrix_interface>`) have a builtin
function to calculate expected values from existing matrix data called
:func:`~fanc.matrix.RegionMatrixContainer.expected_values`. This function calculates and returns
intra-chromosomal, intra-chromosomal per chromosome, and inter-chromosomal expected values.

.. literalinclude:: code/oe_example_code.py
    :language: python
    :start-after: start snippet oe basic
    :end-before: end snippet oe basic

Here, ``intra_expected`` is a list of average (/expected) contact values, where the position of
the value in the list corresponds to the separation between genomic regions in bins.
``intra_expected_chromosome`` is a dictionary with chromosome names as keys, and an expected
value list as value calculated on a per-chromosome basis. ``inter_expected`` is a single, average
inter-chromosomal contact value.

The expected values are typically plotted on a log-log scale, as illustrated here using chromosome 19:

.. literalinclude:: code/oe_example_code.py
    :language: python
    :start-after: start snippet oe ddplot
    :end-before: end snippet oe ddplot

.. image:: images/oe_500kb.png

FAN-C also has a built-in function for plotting the expected values,
:func:`~fanc.plotting.distance_decay_plot`. Additional named arguments
are passed on to ``ax.plot``, for example to change the line color.
The function returns a ``matplotlib`` axes object, which can then be further customised:

.. literalinclude:: code/oe_example_code.py
    :language: python
    :start-after: start snippet oe ddbuiltin
    :end-before: end snippet oe ddbuiltin

.. image:: images/oe_500kb_builtin.png

To compare the expected values of multiple samples, just provide multiple Hic objects:

.. literalinclude:: code/oe_example_code.py
    :language: python
    :start-after: start snippet oe multi
    :end-before: end snippet oe multi

.. image:: images/oe_500kb_multi.png

Note: as Hi-C matrices are normalised on a per-chromosome basis in FAN-C by default, it would be misleading
to plot the overall normalised intra-chromosomal expected values, or to use them for downstream analysis.
We can, however, also calculate the unnormalised expected values easily enough.

.. literalinclude:: code/oe_example_code.py
    :language: python
    :start-after: start snippet oe nonorm
    :end-before: end snippet oe nonorm

.. image:: images/oe_500kb_nonorm.png

If you are simply interested in plotting the unnormalised values, you can use

.. literalinclude:: code/oe_example_code.py
    :language: python
    :start-after: start snippet oe builtinnonorm
    :end-before: end snippet oe builtinnonorm

.. image:: images/oe_500kb_builtinnonorm.png

Expected values rarely need to be calculated explicitly in FAN-C analysis functions, but will be calculated
(or retrieved) on demand whenever necessary. To obtain observed/expected matrices, for example, please
refer to :ref:`matrix_interface`.
