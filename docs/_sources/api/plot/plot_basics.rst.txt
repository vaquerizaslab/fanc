.. _api_plot_basics:

===================
Plotting API basics
===================

FAN-C plotting functions are imported in a Python script or console using

.. literalinclude:: code/matrix_plot_examples.py
    :language: python
    :start-after: start snippet fancplot import
    :end-before: end snippet fancplot import

This will provide access to the FAN-C plot types and custom colormaps for matrix visualisation.

For the following sections we will load a :class:`~fanc.hic.Hic` file for plotting:

.. literalinclude:: code/matrix_plot_examples.py
    :language: python
    :start-after: start snippet fancplot load hic
    :end-before: end snippet fancplot load hic


***********
BasePlotter
***********

Each plot type is based on the :class:`~fanc.plotting.base_plotter.BasePlotter` class.
It controls the formatting of the X axis in genomic coordinates (kb, mb, etc), and provides
a :func:`~fanc.plotting.base_plotter.BasePlotter.plot` function, common to all FAN-C plot types,
which accepts genomic interval definitions in the form of strings (``chr18:6mb-8mb``) or
:class:`~genomic_regions.GenomicRegion` objects.

As an example of basic functionality, we will be looking at the :class:`~fanc.plotting.TriangularMatrixPlot`:

.. literalinclude:: code/matrix_plot_examples.py
    :language: python
    :start-after: start snippet fancplot triangular string
    :end-before: end snippet fancplot triangular string

.. image:: images/plot_triangular.png

The first line sets up the triangular matrix plot parameters. Except for ``vmax`` we have kept the defaults.
The second line is where the actual plotting to a genomic region happens and the plot components are
assembled. The final line simply opens an interactive plotting window.

To save the plot to file, you can use ``hp.save('/path/to/file.png')``.


You can control the formatting of the X axis using a number of different parameters.
``draw_ticks`` can be set to ``False`` to remove the major (and minor) tick marks at
major genomic locations. ``draw_major_ticks`` and ``draw_minor_ticks`` control the drawing of
major and minor ticks, respectively, in the same manner. ``draw_minor_ticks`` is ``False``
by default. To display a small legend in the bottom right of the plot that shows the
distance between major ticks, between minor ticks, and the entire plotting range, set
``draw_tick_legend`` to True. You can remove the chromosome label at the first location by
setting ``draw_chromosome_label`` to ``False``. You can invert the X axis (and the plot)
by setting ``invert_x`` to ``True``, or remove the X axis entirely by setting ``draw_x_axis``
to ``False``.


.. _api_plot_matplotlib:

*******************************
Using dedicated Matplotlib Axes
*******************************

For the highest level of control over your plot, you can build your figures on top of
`Matplotlib <https://matplotlib.org/>`_. This enables you to modify each aspect of the
plot after it has been generated. To start, import matplotlib:

.. literalinclude:: code/matrix_plot_examples.py
    :language: python
    :start-after: start snippet fancplot matplotlib
    :end-before: end snippet fancplot matplotlib

And use the ``ax`` argument of the FAN-C plot to supply it with a custom axis.

.. literalinclude:: code/matrix_plot_examples.py
    :language: python
    :start-after: start snippet fancplot plt axes
    :end-before: end snippet fancplot plt axes

.. image:: images/plot_triangular_plt.png

Matplotlib offers an incredible amount of customisation and plotting options, which are too
numerous to cover here, but we encourage you to study the `Matplotlib documentation <https://matplotlib.org/>`_
for getting your plots to look perfect!
