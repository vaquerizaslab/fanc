.. _fanc-plotting:

============
Plotting API
============

Here is an overview of the plot types and functions available in the plotting API.

Index of plot classes
'''''''''''''''''''''

.. currentmodule:: fanc.plotting

.. autosummary::

    GenomicFigure
    TriangularMatrixPlot
    SquareMatrixPlot
    SplitMatrixPlot
    LinePlot
    BigWigPlot
    GenePlot
    GenomicFeaturePlot
    HicComparisonPlot2D
    HicSlicePlot
    HicPeakPlot
    VerticalSplitPlot
    GenomicVectorArrayPlot
    GenomicRegionsPlot
    GenomicFeatureScorePlot
    FeatureLayerPlot
    GenomicDataFramePlot
    HighlightAnnotation
    SymmetricNorm
    LimitGroup
    MirrorMatrixPlot
    distance_decay_plot
    aggregate_plot
    saddle_plot
    pca_plot

Description and examples
''''''''''''''''''''''''

.. automodule:: fanc.plotting
    :members:
        TriangularMatrixPlot, SquareMatrixPlot, SplitMatrixPlot, HicComparisonPlot2D,
        HicSlicePlot, HicPeakPlot, VerticalSplitPlot, GenomicVectorArrayPlot,
        GenomicFeaturePlot, GenomicRegionsPlot, GenomicFeatureScorePlot,
        GenomicFigure, BigWigPlot, LinePlot,
        GenePlot, FeatureLayerPlot, GenomicDataFramePlot, HighlightAnnotation,
        SymmetricNorm, LimitGroup, MirrorMatrixPlot, distance_decay_plot, saddle_plot,
        aggregate_plot, pca_plot
    :special-members: __init__

    .. data:: example_data

        dict with user-specific paths to example data included in fanc.

    .. data:: style_ticks_whitegrid

        Seaborn style that can be passed to ``axes_style`` attribute
        of plots. Draws a grid of lines on the plot.
