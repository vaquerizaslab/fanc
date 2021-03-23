from ..pairs import ReadPairs
from .base_plotter import GenomeCoordFormatter
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
import seaborn as sns
import numpy as np
import itertools
from sklearn.decomposition import PCA


__all__ = ['summary_statistics_plot', 'ligation_bias_plot', 'restriction_site_distance_plot',
           'marginals_plot', 'distance_decay_plot', 'aggregate_plot', 'saddle_plot',
           'pca_plot']


def summary_statistics_plot(stats, ax=None, exclude=None, include=None, **kwargs):
    """
    Barplot with filter statistics for FAN-C objects.

    :param stats: dictionary of the form {filter_name: int, ...}
    :param ax: (optional) axis in which to plot the bars.
               Will use :code:`plt.gca()` if not specified
    :param exclude: Exclude this list of mask keys from the plot
    :param kwargs: Keyword options passed to Seaborn :code:`barplot`
    :return: Matplotlib axis
    """
    if ax is None:
        ax = plt.gca()

    if exclude is None:
        exclude = set()
    else:
        exclude = set(exclude)

    if include is None:
        include = set(stats.keys())
    else:
        include = set(include)

    palette = kwargs.pop('palette', 'colorblind')
    labels = []
    values = []

    if 'total' in stats and 'total' not in exclude:
        labels.append('total')
        values.append(stats['total'])

    if 'valid' in stats and 'valid' not in exclude:
        labels.append('valid')
        values.append(stats['valid'])

    for key, value in sorted(stats.items()):
        if key in ('total', 'valid') or key in exclude:
            continue
        if key not in include:
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


def distance_decay_plot(*matrices, ax=None, chromosome=None, labels=None, tight=True,
                        norm=True, **kwargs):
    """
    An distance decay (expected values) plot.

    :param matrices: Hi-C objects to be plotted
    :param ax: Optional matplotlib ax object for the plot.
               If not specified, will use ``plt.gca()``
    :param chromosome: Optional, but recommended for the default
                       chromosome-normalised matrices. The name of a chromosome to plot
    :param labels: Optional labels when providing multiple objects.
                   Will replace default labels in the legend. Must be the
                   same number as matrix objects
    :param tight: If True, uses tight figure layout. Disable for grid-based plotting.
    :param kwargs: Parameters passed on to ``ax.plot``
    :return: ax
    """
    if labels is None:
        labels = ['Matrix {}'.format(i) for i in range(len(matrices))]
    elif len(labels) != len(matrices):
        raise ValueError("Number of matrices ({}) must be equal "
                         "to number of labels ({})".format(len(matrices), len(labels)))

    if ax is None:
        ax = plt.gca()

    for i, matrix in enumerate(matrices):
        ex, ex_chromosome, ex_inter = matrix.expected_values(norm=norm)

        if chromosome is not None:
            ex = ex_chromosome[chromosome]

        bin_size = matrix.bin_size
        distances = np.arange(0, bin_size * len(ex), bin_size)

        ax.plot(distances, ex, label=labels[i], **kwargs)

    ax.set_xscale('log')
    ax.set_yscale('log')

    if len(matrices) > 1:
        ax.legend()

    ax.set_ylabel('Expected contact strength')
    ax.set_xlabel('Genomic distance')
    ax.xaxis.set_major_formatter(GenomeCoordFormatter(chromosome if chromosome is not None else "All",
                                                      minor_div=5,
                                                      display_chromosome=False,
                                                      display_scale=False))
    if tight:
        ax.figure.tight_layout()
    return ax


def pca_plot(pca_res, variance=None, eigenvectors=(0, 1),
             markers=None, colors=None, names=None, ax=None):
    """
    Plot the results of a Hi-C PCA analysis from :func:`~fanc.architecture.comparisons.hic_pca`.

    :param pca_res: The PCA result from :func:`~fanc.architecture.comparisons.hic_pca`
    :param variance: A vector specifying the explained variance of each EV in the PCA or
                     the PCA object from :func:`~fanc.architecture.comparisons.hic_pca`.
                     Optional, used to display the explained variance along the axes.
    :param eigenvectors: Tuple of length two specifying which eigenvectors (EVs) to plot.
                         0-based, (0, 1) by default for the first to EVs.
    :param markers: List of marker definitions from matplotlib (e.g. ["o", "*", "s"]).
                    Must be same length as number of samples in PCA.
    :param colors: List of colour definitions from matplotlib.
                   Must be same length as number of samples in PCA.
    :param names: Sample names for plot legend. Must be same length as number of
                  samples in PCA.
    :param ax: Optional matplotlib axes object to plot into. Otherwise uses ``plt.gca()``
    :return: figure, ax
    """
    if markers is None:
        markers = ('^', 'o', '*', 's', 'D', 'v', 'd', 'H', 'p', '>')
    if colors is None:
        colors = ('red', 'blue', 'green', 'purple', 'yellow', 'black',
                  'orange', 'pink', 'cyan', 'lawngreen')
    markers = itertools.cycle(markers)
    colors = itertools.cycle(colors)

    xlabel = 'PC{}'.format(eigenvectors[0] + 1)
    ylabel = 'PC{}'.format(eigenvectors[1] + 1)
    if variance is not None:
        if isinstance(variance, PCA):
            variance = variance.explained_variance_ratio_
        xlabel += ' (%d%%)' % int(variance[eigenvectors[0]]*100)
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


def aggregate_plot(aggregate_matrix, labels=None, vmin=None, vmax=None,
                   oe=False, log=False, colormap='bwr', ax=None, cax=None,
                   relative_label_locations=(0, 0.5, 1), plot_colorbar=True,
                   lower_triangular=False):
    if ax is None:
        ax = plt.gca()

    m = aggregate_matrix.matrix()

    if labels is None:
        labels = ['', '', '']

    if vmin is None:
        vmin = np.nanmin(m)
    if vmax is None:
        vmax = np.nanmax(m)

    if oe and log:
        abs_max = max(abs(vmin), abs(vmax))
        vmin, vmax = -1 * abs_max, abs_max

    im = ax.imshow(m, cmap=colormap, vmin=vmin, vmax=vmax, interpolation='nearest')
    if plot_colorbar:
        plt.colorbar(im, cax=cax)

    pixels = m.shape[0]
    ticks = [loc * pixels if loc != 1 else pixels - 1 for loc in relative_label_locations]

    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)
    
    # ensure same axis limits
    # need invert x here or plot will not be intuitive
    # due to imshow origin location in the top left corner
    ax.set_ylim(ax.get_xlim()[::-1])
    
    # for lower triangular just flip the matrix both ways
    if lower_triangular:
        ax.set_xlim(ax.get_xlim()[::-1])
        ax.set_ylim(ax.get_ylim()[::-1])

    return ax


def saddle_plot(ab_enrichment_matrix, cutoffs, colormap='RdBu_r',
                vmin=-0.75, vmax=0.75, only_gc=False, fig=None,
                axes=None, margin=1):

    if fig is None and axes is None:
        fig = plt.figure(figsize=(5, 5), dpi=300)

    if axes is None:
        gs = grd.GridSpec(5, 5,
                          height_ratios=[margin, 5, 1, 1, margin],
                          width_ratios=[margin, 5, 1, 1, margin])
        heatmap_ax = plt.subplot(gs[1, 1])
        barplot_ax = plt.subplot(gs[3, 1])
        cax = plt.subplot(gs[1, 3])
    else:
        heatmap_ax, barplot_ax, cax = axes

    im = None
    if heatmap_ax is not None:
        im = heatmap_ax.imshow(ab_enrichment_matrix, cmap=colormap, vmin=vmin, vmax=vmax,
                               interpolation='nearest', aspect='auto')
        heatmap_ax.set_xticks([0, ab_enrichment_matrix.shape[1] - 1])
        heatmap_ax.set_xticklabels(['active', 'inactive'])
        xlabels = heatmap_ax.get_xticklabels()
        xlabels[0].set_horizontalalignment('left')
        xlabels[1].set_horizontalalignment('right')

        heatmap_ax.set_yticks([0, ab_enrichment_matrix.shape[1] - 1])
        heatmap_ax.set_yticklabels(['active', 'inactive'], rotation=90)
        ylabels = heatmap_ax.get_yticklabels()
        ylabels[0].set_verticalalignment('bottom')
        ylabels[1].set_verticalalignment('top')

        heatmap_ax.set_ylim(heatmap_ax.get_xlim())

    if cax is not None and im is not None:
        cb = plt.colorbar(im, cax=cax)
        cb.set_ticks([vmin, 0, vmax])
        cb.set_label("log O/E")

    if barplot_ax is not None:
        pos = np.arange(ab_enrichment_matrix.shape[1])
        barplot_ax.bar(pos, cutoffs, color='grey', width=1)
        if not only_gc:
            extent = max(abs(cutoffs[0]), abs(cutoffs[-1]))
            barplot_ax.set_yticks([-1 * extent, 0, extent])
        else:
            barplot_ax.set_yticks([cutoffs[0], cutoffs[int(len(cutoffs) / 2)], cutoffs[1]])
        barplot_ax.set_xlim(heatmap_ax.get_xlim())
        barplot_ax.get_xaxis().set_visible(False)
        barplot_ax.spines['right'].set_visible(False)
        barplot_ax.spines['top'].set_visible(False)
        barplot_ax.spines['bottom'].set_visible(False)
        barplot_ax.yaxis.set_ticks_position('left')
        barplot_ax.xaxis.set_ticks_position('none')
        barplot_ax.set_ylabel("EV percentile\ncutoffs")

    return fig, [heatmap_ax, barplot_ax, cax]
