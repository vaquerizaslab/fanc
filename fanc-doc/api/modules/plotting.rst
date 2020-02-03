.. _kaic-plotting:

============
Plotting API
============

Here is an overview of the plot types and functions available in the plotting API.

Index of plot classes
'''''''''''''''''''''

.. currentmodule:: kaic.plotting

.. autosummary::

    GenomicFigure
    HicPlot
    HicPlot2D
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

Description and examples
''''''''''''''''''''''''

.. automodule:: kaic.plotting
    :members:
        HicPlot, HicPlot2D, HicComparisonPlot2D,
        HicSlicePlot, HicPeakPlot, VerticalSplitPlot, GenomicVectorArrayPlot,
        GenomicFeaturePlot, GenomicRegionsPlot, GenomicFeatureScorePlot,
        GenomicFigure, BigWigPlot, LinePlot,
        GenePlot, FeatureLayerPlot, GenomicDataFramePlot, HighlightAnnotation,
        SymmetricNorm, LimitGroup
    :special-members: __init__

    .. data:: example_data

        dict with user-specific paths to example data included in kaic.

    .. data:: style_ticks_whitegrid

        Seaborn style that can be passed to ``axes_style`` attribute
        of plots. Draws a grid of lines on the plot.
