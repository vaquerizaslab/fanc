.. _api_compartments:


###############
AB compartments
###############

The following steps assume that you ran the ``fanc auto`` command in :ref:`example-fanc-auto`.
Additionally, we set up the Python session like this:

.. literalinclude:: code/ab_example_code.py
    :language: python
    :start-after: start snippet ab setup
    :end-before: end snippet ab setup

AB correlation matrices can very easily be obtained from Hi-C files using the
:func:`~fanc.architecture.compartments.ABCompartmentMatrix.from_hic` function:

.. literalinclude:: code/ab_example_code.py
    :language: python
    :start-after: start snippet ab matrix
    :end-before: end snippet ab matrix

The ``ab`` object acts like any FAN-C matrix (see :ref:`matrix_interface`), which means you
can query and subset the data any way you like. For example, to get the correlation matrix
of chromosome 18:

.. literalinclude:: code/ab_example_code.py
    :language: python
    :start-after: start snippet ab subset
    :end-before: end snippet ab subset

And to visualise the matrix:

.. literalinclude:: code/ab_example_code.py
    :language: python
    :start-after: start snippet ab klot-correlation
    :end-before: end snippet ab klot-correlation

.. image:: images/ab_1mb_correlation.png


############
Eigenvectors
############

The AB correlation matrix eigenvector (EV) is used to determine if a region is in the active or
the inactive compartment. It's calculation is very straightforward:

.. literalinclude:: code/ab_example_code.py
    :language: python
    :start-after: start snippet ab ev
    :end-before: end snippet ab ev

:func:`~fanc.architecture.compartments.ABCompartmentMatrix.eigenvector` returns a numpy
array with one entry per region in the AB correlation matrix (you can retrieve a matching
list of regions with :func:`~fanc.architecture.compartments.ABCompartmentMatrix.regions`).
You can also retrieve only the EV entries for a specific region using the ``sub_region``
argument, but note that the calculation is always performed on the entire genome first
to avoid biases from subsetting.

.. warning::

   Positive EV entries do not automatically mean a region is in the A compartment. In fact,
   if positive or negative entries are representing the A compartment is dependent on
   the implementation of PCA on the platform you are using. Therefore we strongly recommend
   using additional biological information to determine the correspondence between EV entry
   sign and compartment.

   One option implemented in FAN-C is to use GC content as a proxy for activity, as GC-rich
   regions have been shown to be associated with the active compartment. FAN-C implements the
   use of a genomic FASTA file, to calculate GC content and then choose the EV sign so that
   positive entries correspond to A, and negative entries to the B compartment.

   .. literalinclude:: code/ab_example_code.py
       :language: python
       :start-after: start snippet ab gc-ev
       :end-before: end snippet ab gc-ev

