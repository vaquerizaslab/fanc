from kaic.plotting.hic_plotter import HicPlot, HicPlot2D, HicComparisonPlot2D, HicSideBySidePlot2D
from kaic.plotting.plotter import VerticalSplitPlot, GenomicVectorArrayPlot, GenomicFeaturePlot, GenomicRegionsPlot, \
    GenomicFeatureScorePlot, GenomicMatrixPlot, GenomicFigure, GenomicTrackPlot

from kaic.plotting.colormaps import *
import seaborn as sns
import logging

sns.set_style("ticks")
logging.basicConfig(level=logging.INFO)

sns.plt.register_cmap(name='viridis', cmap=viridis)
sns.plt.register_cmap(name='plasma', cmap=plasma)
sns.plt.register_cmap(name='inferno', cmap=inferno)
sns.plt.register_cmap(name='magma', cmap=magma)
sns.plt.register_cmap(name='RdBuWhitespace_r', cmap=fc_cmap)
