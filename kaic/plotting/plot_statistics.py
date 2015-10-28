import seaborn as sns
import numpy as np


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


