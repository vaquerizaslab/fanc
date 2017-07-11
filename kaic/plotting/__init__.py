"""
Provide plotting functions for genomic data types.

The basic idea is that figures can be composed of multiple panels which are
arranged vertically and share a common x-axis representing genomic coordinates.

Each panel is created separately and then combined into a single figure.
For example, when analyzing Hi-C data it is often interesting to correlate
features in the Hi-C map with ChIP-seq tracks. In that case, one would first
create a :class:`~kaic.plotting.hic_plotter.HicPlot` object, which visualizes
Hi-C data, and then a :class:`~kaic.plotting.plotter.BigWigPlot` object, which
can plot bigwig files that are used during ChIP-seq analysis. Finally, the two
objects are used to create a :class:`~kaic.plotting.plotter.GenomicFigure`.

Examples:

.. code:: python

    import kaic.plotting as kplot

    # Create Hic plot
    hplot = kplot.HicPlot("path/to/data.hic")
    # Create ChIP bigwig plot
    bplot = kplot.BigWigPlot("path/to/chip_data.bigwig")

    # The plots are used to generate a figure by passing them as a list
    # to the GenomicFigure constructor.
    gfig = kplot.GenomicFigure([hplot, bplot])

    # Plot a specific region of the genome
    fig, axes = gfig.plot("2:3000000-3500000")

The plot() function returns standard matplotlib Figure and Axes instances
which can be further adjusted using standard matplotlib methods.
"""

from kaic.plotting.hic_plotter import HicPlot, HicPlot2D, HicComparisonPlot2D, HicSideBySidePlot2D, \
    HicSlicePlot, HicPeakPlot
from kaic.plotting.plotter import VerticalSplitPlot, GenomicVectorArrayPlot, GenomicFeaturePlot, GenomicRegionsPlot, \
    GenomicFeatureScorePlot, GenomicMatrixPlot, GenomicFigure, GenomicTrackPlot, BigWigPlot, GenePlot, \
    FeatureLayerPlot, GenomicDataFramePlot, VerticalLineAnnotation
from kaic.plotting.helpers import append_axes, absolute_wspace_hspace, SymmetricNorm, \
                                  style_ticks_whitegrid, LimitGroup

from kaic.plotting.colormaps import *
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_style("ticks")

plt.register_cmap(name='viridis', cmap=viridis)
plt.register_cmap(name='plasma', cmap=plasma)
plt.register_cmap(name='inferno', cmap=inferno)
plt.register_cmap(name='magma', cmap=magma)
plt.register_cmap(name='RdBuWhitespace_r', cmap=fc_cmap)
plt.register_cmap(name='germany', cmap=germany_cmap)
plt.register_cmap(name='white_red', cmap=white_red)
plt.register_cmap(name='white_red_r', cmap=white_red_r)
