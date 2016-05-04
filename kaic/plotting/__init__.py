from kaic.plotting.plot_genomic_data import hic_correlation_plot, hic_matrix_plot, hic_contact_plot_linear, \
                                            hic_matrix_diff_plot, hic_matrix_ratio_plot, hic_ma_plot, \
                                            hic_marginals_plot, hic_triangle_plot
from kaic.plotting.plot_statistics import plot_mask_statistics, hic_ligation_structure_biases_plot, \
                                          pairs_re_distance_plot
from kaic.plotting.colormaps import *
import seaborn as sns
import logging

sns.set_style("ticks")
style_ticks_whitegrid = {
    'axes.axisbelow': True,
    'axes.edgecolor': '.15',
    'axes.facecolor': 'white',
    'axes.grid': True,
    'axes.labelcolor': '.15',
    'axes.linewidth': 1.25,
    'figure.facecolor': 'white',
    'font.family': ['sans-serif'],
    'grid.color': '.8',
    'grid.linestyle': '-',
    'image.cmap': 'Greys',
    'legend.frameon': False,
    'legend.numpoints': 1,
    'legend.scatterpoints': 1,
    'lines.solid_capstyle': 'round',
    'text.color': '.15',
    'xtick.color': '.15',
    'xtick.direction': 'out',
    'xtick.major.size': 6,
    'xtick.minor.size': 3,
    'ytick.color': '.15',
    'ytick.direction': 'out',
    'ytick.major.size': 6,
    'ytick.minor.size': 3}


logging.basicConfig(level=logging.INFO)

sns.plt.register_cmap(name='viridis', cmap=viridis)
sns.plt.register_cmap(name='plasma', cmap=plasma)
sns.plt.register_cmap(name='inferno', cmap=inferno)
sns.plt.register_cmap(name='magma', cmap=magma)
