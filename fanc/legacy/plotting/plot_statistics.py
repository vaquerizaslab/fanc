import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import itertools
import tables as t
from fanc.plotting.plot_genomic_data import _prepare_backend, _plot_figure


def statistics_plot(stats, ax=None):
    if ax is None:
        ax = plt.gca()
    labels = []
    values = []

    if 'total' in stats:
        labels.append('total')
        values.append(stats['total'])

    if 'unmasked' in stats:
        labels.append('valid')
        values.append(stats['unmasked'])

    for key, value in sorted(stats.items()):
        if key in ('total', 'unmasked'):
            continue
        labels.append(key)
        values.append(value)

    barplot = sns.barplot(x=np.array(labels), y=np.array(values), palette="muted", ax=ax)
    sns.despine()
    return barplot


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
        old_backend = plt.get_backend()
        plt.switch_backend('pdf')
        plt.ioff()

    barplot = sns.barplot(x=np.array(labels), y=np.array(values), palette="muted")
    sns.despine()

    if output is not None:
        barplot.figure.savefig(output)
        plt.close(barplot.figure)
        plt.ion()
        plt.switch_backend(old_backend)
    else:
        plt.show()


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
        fig = plt.figure()
        fig.suptitle("Error structure by distance")
        plt.plot(x, (inward_ratios), 'b', label="inward/same strand")
        plt.plot(x, (outward_ratios), 'r', label="outward/same strand")
        plt.xscale('log')
        if log:
            plt.axhline(y=0, color='black', ls='dashed', lw=0.8)
            plt.ylim(-3, 3)
        else:
            plt.axhline(y=0.5, color='black', ls='dashed', lw=0.8)
            plt.ylim(0, 3)
        plt.xlabel('Gap size between fragments')
        plt.ylabel('Read count ratio')
        plt.legend(loc='upper right')
        sns.despine()
        if output is None:
            plt.show()
        else:
            fig.savefig(output)
            plt.close(fig)
            plt.ion()
            plt.switch_backend(old_backend)


def pairs_re_distance_plot(pairs, output=None, limit=10000, max_distance=None):
    distances = []
    for i, pair in enumerate(pairs.pairs(lazy=True)):
        d1 = pair.left.re_distance()
        d2 = pair.right.re_distance()

        d = d1 + d2

        if max_distance is None or d <= max_distance:
            distances.append(d)
        if limit is not None and i >= limit:
            break

    old_backend = _prepare_backend(output)
    dplot = sns.distplot(distances)
    dplot.set_xlim(left=10)
    dplot.set_xscale('log')
    _plot_figure(dplot.figure, output, old_backend)


def mapq_hist_plot(reads, output=None, include_masked=False):
    reads = reads.reads(lazy=True, include_masked=include_masked)
    mapqs = [r.mapq for r in reads]
    old_backend = _prepare_backend(output)
    mqplot = sns.distplot(mapqs, norm_hist=False, kde=False, bins=np.arange(min(mapqs), max(mapqs)+1.5)-0.5)
    mqplot.set_xlim(left=-1, right=max(mapqs)+2)
    _plot_figure(mqplot.figure, output, old_backend)


