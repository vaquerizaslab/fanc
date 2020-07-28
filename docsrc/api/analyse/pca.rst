.. _api_pca:

###
PCA
###

To follow this tutorial, download the FAN-C example data, for example through our
`Keeper library <https://keeper.mpdl.mpg.de/d/147906745b634c779ed3/>`_. Then set up
your Python session like this, loading some of our previously published
`Low-C datasets <https://doi.org/10.1038/s41467-018-06961-0>`_:

.. literalinclude:: code/pca_example_code.py
    :language: python
    :start-after: start snippet pca setup
    :end-before: end snippet pca setup

PCA is one way in FAN-C to compare different Hi-C matrices to each other. A matrix of
pixels vs matrices is assembled that contains the (normalised) contact strength of each
matrix for the respective pixel (=region pair). PCA is then run on this matrix, and the
resulting eigenvectors can be plotted to examine the variability between datasets.

In FAN-C, simply use :func:`~fanc.architecture.comparisons.hic_pca` for this purpose,
as shown here for chromosome 19:

.. literalinclude:: code/pca_example_code.py
    :language: python
    :start-after: start snippet pca run
    :end-before: end snippet pca run


We can plot the result using :func:`~fanc.plotting.pca_plot`:

.. literalinclude:: code/pca_example_code.py
    :language: python
    :start-after: start snippet pca plot
    :end-before: end snippet pca plot

.. image:: images/pca_default.png

We can easily change the colors and markers, for example by colouring all samples with
MboI and HindIII differently, and assigning different markers to samples with more or
less than 1M cells:

.. literalinclude:: code/pca_example_code.py
    :language: python
    :start-after: start snippet pca adjust
    :end-before: end snippet pca adjust

.. image:: images/pca_adjust.png

Sometimes the first eigenvector captures the library sequencing depth, so you may want to
plot the second and third EVs instead using ``eigenvectors=(1,2)`` (in this case it does
not seem to be particularly informative):

.. literalinclude:: code/pca_example_code.py
    :language: python
    :start-after: start snippet pca ev
    :end-before: end snippet pca ev

.. image:: images/pca_ev.png

This kind of analysis can be tricky, and selecting informative pixels from the matrix can
be key to getting a robust and intuitive PCA result. In the above example, we are using
several parameters to select informative pixels for the PCA. First, we are only using
pixels that are non-zero in all samples with ``ignore_zeros=True``. Second, we are sorting
the pixels, listing the ones with the largest variance across samples first, using
``strategy='variance'``. Finally, we are selecting the top 100k pixels (those with the
largest variance) first with ``sample_size=100000``.

When you are analysing matrices of higher resolution, pixels far away from the diagonal
might be dominated by noise. ``ignore_zeros`` removes most of the noisy pixels, but
additionally you might want to set a ``max_distance`` to only select pixels corresponding
to regions closer than this value. Similarly, if you want to exclude contacts close to
the diagonal, use ``min_distance``.

For more options have a look at the API reference for :func:`~fanc.plotting.pca_plot`.
