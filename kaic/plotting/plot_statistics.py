import seaborn as sns
import numpy as np
from kaic.plotting.plot_genomic_data import _prepare_backend


def plot_mask_statistics(maskable, masked_table, output=None, ignore_zero=True):
    # get statistics
    stats = maskable.mask_statistics(masked_table)

    # calculate total
    total = sum(stats.values())

    labels = ['total', 'unmasked']
    values = [total, stats['unmasked']]

    for item in stats.iteritems():
        if not item[0] == 'unmasked':
            if not ignore_zero or item[1] > 0:
                labels.append(item[0])
                values.append(item[1])

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


def hic_ligation_error_structure_plot(pairs, output=None, data_points=None, skip_self_ligations=True):
    """
    Plot the ligation error structure of a dataset.

    :param pairs: Read pairs mapped to genomic regions (:class:`~FragmentMappedReadPairs`)
    :param output: Path to pdf file to save this plot.
    :param data_points: Number of data points to average per point
                        in the plot. If None (default), this will
                        be determined on a best-guess basis.
    :param skip_self_ligations: If True (default), will not consider
                                self-ligated fragments for assessing
                                the error rates.
    """
    x, inward_ratios, outward_ratios = pairs.get_error_structure(data_points, skip_self_ligations)
    old_backend = _prepare_backend(output)
    fig = sns.plt.figure()
    fig.suptitle("Error structure by distance")
    sns.plt.plot(x, inward_ratios, 'b', label="inward/same strand")
    sns.plt.plot(x, outward_ratios, 'r', label="outward/same strand")
    sns.plt.xscale('log')
    sns.plt.axhline(y=0.5, color='black', ls='dashed')
    sns.plt.ylim(0, 3)
    sns.plt.xlabel('Gap size between fragments')
    sns.plt.ylabel('Read count ratio')
    sns.plt.legend(loc='upper right')
    if output is None:
        sns.plt.show()
    else:
        fig.savefig(output)
        sns.plt.close(fig)
        sns.plt.ion()
        sns.plt.switch_backend(old_backend)
