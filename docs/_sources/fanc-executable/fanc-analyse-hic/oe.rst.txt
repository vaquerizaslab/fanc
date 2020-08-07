.. _fanc-oe:

###############
Expected values
###############

.. note::

    The following examples use the matrix files in FAN-C format. If you want to try the first few
    commands using Juicer ``.hic`` files, replace ``output/hic/binned/fanc_example_500kb.hic``
    with ``architecture/other-hic/fanc_example.juicer.hic@500kb``. If you want to work with
    Cooler files in this tutorial, use ``architecture/other-hic/fanc_example.mcool@500kb``.
    The results will be minimally different due to the "zooming" and balancing applied by
    each package.

The contact intensity in a Hi-C matrix gets progressively weaker the further apart
two loci are. The expected values follow a distinctive profile with distance for
Hi-C matrices, which can be approximated by a power law and forms an almost straight
line in a log-log plot.

To calculate the expected values of any FAN-C compatible matrix, you can use the
``fanc expected`` command:

.. argparse::
   :module: fanc.commands.fanc_commands
   :func: expected_parser
   :prog: fanc expected
   :nodescription:
   :nodefault:

*******
Example
*******

The following example calculates and plots the expected values for a 500kb resolution
Hi-C matrix of chromosome 19.

.. literalinclude:: code/oe_example_code
    :language: bash
    :start-after: start snippet expected basic
    :end-before: end snippet expected basic

The resulting plot (from ``-p``) looks like this:

.. image:: images/fanc_example_500kb_expected.png

The actual expected values are stored in ``architecture/expected/fanc_example_500kb_expected.txt``:

.. code::

    distance	Matrix_0
    0	0.24442297400748084
    500000	0.07759323503191953
    1000000	0.03699383283713825
    1500000	0.02452933204893787
    2000000	0.017725227895561607
    2500000	0.014272302693312262
    3000000	0.011708011997703627
    3500000	0.010125456912234796
    ...


*******
Options
*******

The expected values are stored in the matrix. If you are running any command that relies on
the expected values again, it will be retrieved rather than recalculated. Use ``--recalculate``
to force a re-calculation of expected values, for whatever reason.

It may be interesting to plot the expected values of unnormalised matrices, to see any ranges
where contacts are more or less abundant before normalisation. Use ``-N`` to plot the unnormalised
expected values.


*************************
Comparing expected values
*************************

When you are providing more than one matrix as input to ``fanc expected``, the expected values
for all matrices will be written to file and plotted if using the ``-p`` option:

.. literalinclude:: code/oe_example_code
    :language: bash
    :start-after: start snippet expected multi
    :end-before: end snippet expected multi

.. image:: images/expected_multi.png


************
O/E matrices
************

Using ``fancplot``, we can visualise the observed/expected Hi-C matrix, which normalised each matrix
value to its given expected value at that distance. Here, we are showing a log2-transformed
O/E matrix:

.. literalinclude:: code/oe_example_code
    :language: bash
    :start-after: start snippet expected fancplot
    :end-before: end snippet expected fancplot

.. image:: images/fanc_example_500kb_chr18_oe.png

