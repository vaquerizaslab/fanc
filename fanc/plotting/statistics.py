from ..pairs import ReadPairs
from .base_plotter import GenomeCoordFormatter
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import itertools


def summary_statistics_plot(stats, ax=None, **kwargs):
    """
    Barplot with filter statistics for FAN-C objects.

    :param stats: dictionary of the form {filter_name: int, ...}
    :param ax: (optional) axis in which to plot the bars.
               Will use :code:`plt.gca()` if not specified
    :param kwargs: Keyword options passed to Seaborn :code:`barplot`
    :return: Matplotlib axis
    """
    if ax is None:
        ax = plt.gca()

    palette = kwargs.pop('palette', 'colorblind')
    labels = []
    values = []

    if 'total' in stats:
        labels.append('total')
        values.append(stats['total'])

    if 'valid' in stats:
        labels.append('valid')
        values.append(stats['valid'])

    for key, value in sorted(stats.items()):
        if key in ('total', 'valid'):
            continue
        labels.append(key)
        values.append(value)

    sns.barplot(x=np.array(labels), y=np.array(values), ax=ax,
                palette=palette, **kwargs)
    sns.despine()

    ax.set_xticklabels([label.get_text()
                        for label in ax.get_xticklabels()],
                       rotation=90)

    ax.set_ylabel("Number of read pairs")
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
    ax.figure.tight_layout()

    return ax


def ligation_bias_plot(pairs, ax=None, log=False, **kwargs):
    """
    Plot the ligation error structure of a dataset.

    :param pairs: Read pairs mapped to genomic regions (:class:`~fanc.pairs.ReadPairs`)
    :param ax: (optional) axis in which to plot the bars.
               Will use :code:`plt.gca()` if not specified
    :param log: log2-transform ratios if True
    :param kwargs: Additional arguments to pass
                   to :code:`plot`
    """
    if ax is None:
        ax = plt.gca()

    if isinstance(pairs, ReadPairs):
        x, inward_ratios, outward_ratios, bins_sizes = pairs.get_ligation_structure_biases()
    else:
        x, inward_ratios, outward_ratios, bins_sizes = pairs

    if log:
        inward_ratios = np.log2(inward_ratios) + 1
        outward_ratios = np.log2(outward_ratios) + 1

    with sns.axes_style("white", {
            "legend.frameon": True,
            "xtick.major.size": 4,
            "xtick.minor.size": 2,
            "ytick.major.size": 4,
            "ytick.minor.size": 2,
            "axes.linewidth": 0.5
    }):
        ax.set_title("Error structure by distance")
        ax.plot(x, inward_ratios, 'b', label="inward/same strand", **kwargs)
        ax.plot(x, outward_ratios, 'r', label="outward/same strand", **kwargs)
        ax.set_xscale('log')
        if log:
            ax.axhline(y=0, color='black', ls='dashed', lw=0.8)
            ax.set_ylim((-3, 3))
        else:
            ax.axhline(y=0.5, color='black', ls='dashed', lw=0.8)
            ax.set_ylim((0, 3))
        ax.set_xlabel('Gap size between fragments')
        ax.set_ylabel('Read count ratio')
        ax.legend(loc='upper right')
        sns.despine()

    ax.figure.tight_layout()
    return ax


def restriction_site_distance_plot(pairs, ax=None, max_percentile=95,
                                   sample=100000, max_distance=None,
                                   **kwargs):
    """
    Plot the distribution of read pair restriction site distances.

    The sum of distances to the nearest restriction site of each mate
    is equivalent to the insert size of the sequenced DNA molecule. As
    such, this plot can serve as a quality control for DNA fragmentation
    prior to sequencing. It can also be used to derive a cutoff for the
    :class:`~fanc.pairs.ReDistanceFilter`, which excludes mate pairs with
    very large insert sizes.

    :param pairs: :class:`~fanc.pairs.ReadPairs`
    :param ax: (optional) matplotlib axis
    :param max_percentile: Percentile of values up to which the distribution
                           is plotted. If this is set to 100, the distribution
                           will be squeezed in a small section on
                           the x axis if there are extremely large values
                           in the distribution. The default (95) makes sure
                           the distribution is properly visible.

    :param sample: If this is an integer, only <sample> random number of mate
                   pairs are plotted in the distribution to save time. If this
                   is set to None, all mate pairs are plotted. By default, a
                   sample of 100000 is shown.
    :param max_distance: If this is different from None, distances larger
                         than <max_distance> are ignored.
    :param kwargs: Keyword arguments passed to :code:`seaborn.distplot`
    :return: ax
    """
    if ax is None:
        ax = plt.gca()

    color = kwargs.pop('color', '#FBAFE4')

    distances = []
    for i, pair in enumerate(pairs.pairs(lazy=True)):
        d = pair.left.re_distance() + pair.right.re_distance()

        if max_distance is None or d <= max_distance:
            distances.append(d)
        if sample is not None and i >= sample:
            break

    distances = np.array(distances)
    if max_percentile is not None:
        median_insert, high = np.nanpercentile(distances, [50, max_percentile])
        distances = distances[distances < high]
        ax.set_xlim((0, high))
    else:
        median_insert = np.nanmedian(distances)
        ax.set_xlim(left=0)

    sns.distplot(distances, ax=ax, color=color, **kwargs)
    ax.axvline(x=median_insert, color='grey', linestyle='--')
    ax.set_xlabel("Sum of restriction site distances\n(insert size)")
    ax.set_ylabel("Density")
    sns.despine(ax=ax)
    ax.figure.tight_layout()

    return ax


def marginals_plot(matrix, chromosome, ax=None, lower=None, rel_cutoff=0.1, color='#3E896D',
                   **kwargs):
    """
    Plot Hi-C marginals vector.

    Marginals are the sum of values in each column of the matrix.
    This plot can be used to determine sensible low coverage thresholds
    for the :class:`~fanc.hic.LowCoverageFilter`.

    :param matrix: :class:`~fanc.matrix.RegionMatrixContainer`
    :param chromosome: Name of a chromosome to plot marginals for.
    :param ax: Matplotlib axis
    :param lower: Absolute lower cutoff for drawing threshold line
    :param rel_cutoff: Relative lower cutoff for drawing thresold line
    :param color: Color of line in plot
    :param kwargs: Keyword arguments passed to
                   :func:`~fanc.matrix.RegionMatrixContainer.matrix`
    :return: ax
    """
    if ax is None:
        ax = plt.gca()

    bin_size = matrix.bin_size

    kwargs['key'] = (chromosome, chromosome)
    marginals = matrix.marginals(**kwargs)

    if lower is None and rel_cutoff is not None:
        median = np.nanmedian(marginals[marginals > 0])
        lower = rel_cutoff * median

    distances = np.arange(0, bin_size * len(marginals), bin_size)
    ax.plot(distances, marginals, color=color)

    if lower is not None:
        ax.axhline(lower, color='grey', linestyle='--')

    ax.xaxis.set_major_formatter(GenomeCoordFormatter(chromosome,
                                                      minor_div=5,
                                                      display_chromosome=True,
                                                      display_scale=False))
    ax.set_ylabel("Sum of counts")
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
    ax.set_ylim(bottom=0)
    sns.despine(ax=ax)
    ax.figure.tight_layout()

    return ax


def distance_decay_plot(*matrices, ax=None, chromosome=None, **kwargs):
    if ax is None:
        ax = plt.gca()

    for matrix in matrices:
        ex, ex_chromosome, ex_inter = matrix.expected_values()

        if chromosome is not None:
            ex = ex_chromosome[chromosome]

        bin_size = matrix.bin_size
        distances = np.arange(0, bin_size * len(ex), bin_size)
        ax.plot(distances, ex, **kwargs)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_ylabel('Expected value')
    ax.set_xlabel('Genomic separation')
    ax.xaxis.set_major_formatter(GenomeCoordFormatter(chromosome if chromosome is not None else "All",
                                                      minor_div=5,
                                                      display_chromosome=False,
                                                      display_scale=False))
    ax.figure.tight_layout()
    return ax


def pca_plot(pca_res, variance=None, eigenvectors=(0, 1),
             markers=None, colors=None, names=None, ax=None):
    if markers is None:
        markers = ('^', 'o', '*', 's', 'D', 'v', 'd', 'H', 'p', '>')
    if colors is None:
        colors = ('red', 'blue', 'green', 'purple', 'yellow', 'black',
                  'orange', 'pink', 'cyan', 'lawngreen')
    markers = itertools.cycle(markers)
    colors = itertools.cycle(colors)

    xlabel = 'PC1'
    if variance is not None:
        xlabel += ' (%d%%)' % int(variance[eigenvectors[0]]*100)

    ylabel = 'PC2'
    if variance is not None:
        ylabel += ' (%d%%)' % int(variance[eigenvectors[1]]*100)

    if ax is None:
        if names is not None:
            ax = plt.subplot(121)
        else:
            ax = plt.subplot(111)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.set_title('PCA on %d samples' % pca_res.shape[0])

    for i in range(pca_res.shape[0]):
        name = names[i] if names is not None else None
        ax.plot(pca_res[i, eigenvectors[0]], pca_res[i, eigenvectors[1]],
                marker=next(markers), color=next(colors), label=name)

    if names is not None:
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    return ax.figure, ax

