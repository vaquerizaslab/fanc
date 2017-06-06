from kaic.plotting.hic_plotter import HicPlot, HicPlot2D, HicComparisonPlot2D, HicSideBySidePlot2D, \
    HicSlicePlot, HicPeakPlot
from kaic.plotting.plotter import VerticalSplitPlot, GenomicVectorArrayPlot, GenomicFeaturePlot, GenomicRegionsPlot, \
    GenomicFeatureScorePlot, GenomicMatrixPlot, GenomicFigure, GenomicTrackPlot, BigWigPlot, GenePlot, \
    FeatureLayerPlot, GenomicDataFramePlot, VerticalLineAnnotation
from kaic.plotting.helpers import append_axes, absolute_wspace_hspace, SymmetricNorm, \
                                  style_ticks_whitegrid

from kaic.plotting.colormaps import *
import seaborn as sns

sns.set_style("ticks")

sns.plt.register_cmap(name='viridis', cmap=viridis)
sns.plt.register_cmap(name='plasma', cmap=plasma)
sns.plt.register_cmap(name='inferno', cmap=inferno)
sns.plt.register_cmap(name='magma', cmap=magma)
sns.plt.register_cmap(name='RdBuWhitespace_r', cmap=fc_cmap)
sns.plt.register_cmap(name='germany', cmap=germany_cmap)
sns.plt.register_cmap(name='white_red', cmap=white_red)
sns.plt.register_cmap(name='white_red_r', cmap=white_red_r)
