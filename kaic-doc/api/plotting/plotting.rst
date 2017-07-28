Plotting API
============

Index of plot classes
'''''''''''''''''''''

.. currentmodule:: kaic.plotting

.. autosummary::

    GenomicFigure
    HicPlot
    HicPlot2D
    BigWigPlot
    GenePlot
    GenomicFeaturePlot
    HicComparisonPlot2D
    HicSideBySidePlot2D
    HicSlicePlot
    HicPeakPlot
    VerticalSplitPlot
    GenomicVectorArrayPlot
    GenomicRegionsPlot
    GenomicFeatureScorePlot
    GenomicMatrixPlot
    GenomicTrackPlot
    FeatureLayerPlot
    GenomicDataFramePlot
    HighlightAnnotation
    SymmetricNorm
    LimitGroup

Description and examples
''''''''''''''''''''''''

.. automodule:: kaic.plotting
    :members:
        HicPlot, HicPlot2D, HicComparisonPlot2D, HicSideBySidePlot2D,
        HicSlicePlot, HicPeakPlot, VerticalSplitPlot, GenomicVectorArrayPlot,
        GenomicFeaturePlot, GenomicRegionsPlot, GenomicFeatureScorePlot,
        GenomicMatrixPlot, GenomicFigure, GenomicTrackPlot, BigWigPlot,
        GenePlot, FeatureLayerPlot, GenomicDataFramePlot, HighlightAnnotation,
        SymmetricNorm, LimitGroup
    :special-members: __init__

    .. data:: example_data

        dict with user-specific paths to example data included in kaic.

    .. data:: style_ticks_whitegrid

        Seaborn style that can be passed to ``axes_style`` attribute
        of plots. Draws a grid of lines on the plot.
