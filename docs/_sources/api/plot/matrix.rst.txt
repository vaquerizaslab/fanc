.. _api_plot_matrix:

============
Matrix plots
============

After importing ``fanc.plotting`` as ``fancplot``, we have access to a range of Hi-C (and related)
matrix plots. We already used the :class:`~fanc.plotting.hic_plotter.TriangularHicPlot` as an example in
:ref:`api_plot_basics`. In total there are four types of matrix visualisation built in to FAN-C:

- :class:`~fanc.plotting.hic_plotter.SquareHicPlot`, which shows a square matrix, with colors typically
  representing contact strength between the two genomic regions on the X, and Y axes, respectively.
- :class:`~fanc.plotting.hic_plotter.TriangularHicPlot`, which rotates the matrix by 45 degrees, so
  instead of a square, it now becomes a triangle. This can be particularly useful if features close to
  the matrix diagonal, such as TADs, are of interest, and if additional genomic tracks, such as CTCF
  binding, should be shown underneath
- :class:`~fanc.plotting.hic_plotter.SplitMatrixPlot`, which is similar to the
  :class:`~fanc.plotting.hic_plotter.SquareHicPlot`, but the values above and below the diagonal
  are from different matrices. This is a great way for direct matrix feature comparisons in the same
  plot.
- :class:`~fanc.plotting.plotter.MirrorMatrixPlot`, which can be used to show two matrix plots,
  particularly of the type :class:`~fanc.plotting.hic_plotter.TriangularHicPlot`, mirrored at the
  X axis


All of these are built on the base class :class:`~:fanc.plotting.base_plotter.BasePlotterMatrix`, so
they share of lot of options and functionality. Let's look at each of them in turn.


*************
Square matrix
*************

To generate a basic square matrix plot, run

.. literalinclude:: code/matrix_plot_examples.py
    :language: python
    :start-after: start snippet fancplot square string
    :end-before: end snippet fancplot square string

.. image:: images/plot_square.png

Currently, the diagonal is about the only thing visible in this example, so let's adjust the saturation
of the plot with ``vmax`` (``vmin`` can be used to set a lower bound)

To generate a basic square matrix plot, run

.. literalinclude:: code/matrix_plot_examples.py
    :language: python
    :start-after: start snippet fancplot square vmax
    :end-before: end snippet fancplot square vmax

.. image:: images/plot_square_vmax.png

FAN-C uses the custom "germany" colormap as a default, built using the
`colormap scripts <https://github.com/bids/colormap>`_ from Stéfan van der Walt and Nathaniel Smith,
which have also been used to generate the default matplotlib colormaps. It is "perceptually uniform",
and thus well-suited to represent contact intensities. We are also including a colormap called
"white_red", for a more classic Hi-C matrix look, but equally perceptually uniform:

.. literalinclude:: code/matrix_plot_examples.py
    :language: python
    :start-after: start snippet fancplot square cmap
    :end-before: end snippet fancplot square cmap

.. image:: images/plot_square_cmap.png

Any `Matplotlib-compatible colormap <https://matplotlib.org/3.1.1/gallery/color/colormap_reference.html>`_
can be used to change the look of your matrix. You can remove the colorbar, for example if you want to
add your own on a different axis, using ``show_colorbar=False``. Alternatively, if you have set up your
:ref:`axes manually <api_plot_matplotlib>`, you can move the colorbar to an axis of your choice using
``cax=<your axis>`` in the constructor.

The exponential decay of contact intensities with distance between any two loci can sometimes make
it difficult to visualise features close to the diagonal and far fro it simultaneously. Sometimes
showing colors on a log scale can improve this:

.. literalinclude:: code/matrix_plot_examples.py
    :language: python
    :start-after: start snippet fancplot square log
    :end-before: end snippet fancplot square log

.. image:: images/plot_square_log.png

You can show uncorrected values using ``matrix_norm=False`` (note the color scale change):

.. literalinclude:: code/matrix_plot_examples.py
    :language: python
    :start-after: start snippet fancplot square uncorrected
    :end-before: end snippet fancplot square uncorrected

.. image:: images/plot_square_uncorrected.png

We can also show log2 observed/expected values using ``oe=True`` and ``log=True`` When not setting the
colorbar limit explicitly, you can also make the colorbar symmetrical around 0 using ``colorbar_symmetry=0``.
We will also use a divergent colormap:

.. literalinclude:: code/matrix_plot_examples.py
    :language: python
    :start-after: start snippet fancplot square oe
    :end-before: end snippet fancplot square oe

.. image:: images/plot_square_oe.png


*****************
Triangular matrix
*****************

Any of the above customisations also work for the triangular version of the matrix plot,
:class:`~fanc.plotting.hic_plotter.TriangularHicPlot`.

.. literalinclude:: code/matrix_plot_examples.py
    :language: python
    :start-after: start snippet fancplot triangular string
    :end-before: end snippet fancplot triangular string

.. image:: images/plot_triangular.png

Additionally, we can control the height, at which the triangle is cut off, and therefore
the maximum distance at which contacts are shown between two regions, with ``max_dist``.
This value is expressed in base pairs, but supports strings in genomic format:

.. literalinclude:: code/matrix_plot_examples.py
    :language: python
    :start-after: start snippet fancplot triangular maxdist
    :end-before: end snippet fancplot triangular maxdist

.. image:: images/plot_triangular_maxdist.png


************
Split matrix
************

The :class:`~fanc.plotting.hic_plotter.SplitMatrixPlot` requires two Hi-C objects, one for
above, and one for below the diagonal. Unless you explicitly set ``scale_matrices`` to ``False``,
the plot will first scale both matrices to the same sequencing depth, in order to make them
more comparable. Depending on your normalisation, this may not be necessary.


*************
Mirrored plot
*************

To showcase the :class:`~fanc.plotting.plotter.MirrorMatrixPlot`, we can plot the Hi-C matrix,
and its O/E transformation at the same time:

.. literalinclude:: code/matrix_plot_examples.py
    :language: python
    :start-after: start snippet fancplot mirror
    :end-before: end snippet fancplot mirror

.. image:: images/plot_mirror.png

Opposed to Split matrix, it is also possible to combine matrices of different resolution, or
from different Hi-C experiments, as long as they are from the same organism.