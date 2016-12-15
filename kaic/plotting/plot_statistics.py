import seaborn as sns
import numpy as np
import itertools
import tables as t
from kaic.plotting.plot_genomic_data import _prepare_backend, _plot_figure


def plot_mask_statistics(maskable, masked_table, output=None, ignore_zero=True):
    # get statistics
    stats = maskable.mask_statistics(masked_table)

    # calculate total
    if isinstance(masked_table, t.Group):
        total = 0
        for table in masked_table:
            total += table._original_len()
    else:
        total = masked_table._original_len()

    labels = ['total', 'unmasked']
    values = [total, stats['unmasked']]

    for key, value in sorted(stats.items()):
        if not key == 'unmasked':
            if not ignore_zero or value > 0:
                labels.append(key)
                values.append(value)

    if output is not None:
        old_backend = sns.plt.get_backend()
        sns.plt.switch_backend('pdf')
        sns.plt.ioff()

    barplot = sns.barplot(x=np.array(labels), y=np.array(values), palette="muted")
    sns.despine()

    if output is not None:
        barplot.figure.savefig(output)
        sns.plt.close(barplot.figure)
        sns.plt.ion()
        sns.plt.switch_backend(old_backend)
    else:
        sns.plt.show()


def hic_ligation_structure_biases_plot(pairs, output=None, log=False, *args, **kwargs):
    """
    Plot the ligation error structure of a dataset.

    :param pairs: Read pairs mapped to genomic regions (:class:`~FragmentMappedReadPairs`)
    :param output: Path to pdf file to save this plot.
    :param *args **kwargs: Additional arguments to pass
                           to :met:`~FragmentMappedReadPairs.get_ligation_structure_biases`
    """
    x, inward_ratios, outward_ratios, bins_sizes = pairs.get_ligation_structure_biases(*args, **kwargs)
    if log:
        inward_ratios = np.log2(inward_ratios) + 1
        outward_ratios = np.log2(outward_ratios) + 1
    old_backend = _prepare_backend(output)
    with sns.axes_style("white", {
            "legend.frameon": True,
            "xtick.major.size": 4,
            "xtick.minor.size": 2,
            "ytick.major.size": 4,
            "ytick.minor.size": 2,
            "axes.linewidth": 0.5
    }):
        fig = sns.plt.figure()
        fig.suptitle("Error structure by distance")
        sns.plt.plot(x, (inward_ratios), 'b', label="inward/same strand")
        sns.plt.plot(x, (outward_ratios), 'r', label="outward/same strand")
        sns.plt.xscale('log')
        if log:
            sns.plt.axhline(y=0, color='black', ls='dashed', lw=0.8)
            sns.plt.ylim(-3, 3)
        else:
            sns.plt.axhline(y=0.5, color='black', ls='dashed', lw=0.8)
            sns.plt.ylim(0, 3)
        sns.plt.xlabel('Gap size between fragments')
        sns.plt.ylabel('Read count ratio')
        sns.plt.legend(loc='upper right')
        sns.despine()
        if output is None:
            sns.plt.show()
        else:
            fig.savefig(output)
            sns.plt.close(fig)
            sns.plt.ion()
            sns.plt.switch_backend(old_backend)


def pairs_re_distance_plot(pairs, output=None, limit=10000, max_distance=None):
    distances = []
    for i, pair in enumerate(pairs.pairs(lazy=True)):
        d1 = pair.left.re_distance()
        d2 = pair.right.re_distance()
        if max_distance is None or d1 <= max_distance:
            distances.append(d1)
        if max_distance is None or d2 <= max_distance:
            distances.append(d2)
        if limit is not None and i >= limit:
            break

    old_backend = _prepare_backend(output)
    dplot = sns.distplot(distances)
    dplot.set_xlim(left=0)
    _plot_figure(dplot.figure, output, old_backend)


def mapq_hist_plot(reads, output=None, include_masked=False):
    reads = reads.reads(lazy=True, include_masked=include_masked)
    mapqs = [r.mapq for r in reads]
    old_backend = _prepare_backend(output)
    mqplot = sns.distplot(mapqs, norm_hist=False, kde=False, bins=np.arange(min(mapqs), max(mapqs)+1.5)-0.5)
    mqplot.set_xlim(left=-1, right=max(mapqs)+2)
    _plot_figure(mqplot.figure, output, old_backend)


def pca_plot(pca_res, pca_info=None, markers=None, colors=None, names=None):
    if markers is None:
        markers = ('^', 'o', '*', 's', 'D', 'v', 'd', 'H', 'p', '>')
    if colors is None:
        colors = ('red', 'blue', 'green', 'purple', 'yellow', 'black', 'orange', 'pink', 'cyan', 'lawngreen')
    markers = itertools.cycle(markers)
    colors = itertools.cycle(colors)

    xlabel = 'PC1'
    if pca_info is not None:
        xlabel += ' (%d%%)' % int(pca_info.explained_variance_ratio_[0]*100)

    ylabel = 'PC2'
    if pca_info is not None:
        ylabel += ' (%d%%)' % int(pca_info.explained_variance_ratio_[1]*100)

    if names is not None:
        ax_main = sns.plt.subplot(121)
    else:
        ax_main = sns.plt.subplot(111)
    ax_main.set_xlabel(xlabel)
    ax_main.set_ylabel(ylabel)

    ax_main.set_title('PCA on %d samples' % pca_res.shape[0])

    for i in range(pca_res.shape[0]):
        name = names[i] if names is not None else None
        ax_main.plot(pca_res[i, 0], pca_res[i, 1], marker=next(markers), color=next(colors), label=name)

    if names is not None:
        ax_main.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    return ax_main.figure, ax_main

