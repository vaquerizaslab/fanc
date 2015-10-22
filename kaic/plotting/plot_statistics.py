import seaborn as sns


def plot_mask_statistics(maskable, masked_table, output=None):
    # get statistics
    stats = maskable.mask_statistics(masked_table)

    # calculate total
    total = sum(stats.values())

    labels = ['total', 'unmasked']
    values = [total, stats['unmasked']]

    for item in stats.iteritems():
        if not item[0] == 'unmasked':
            labels.append(item[0])
            values.append(item[1])

    barplot = sns.barplot(x=labels, y=values, palette="muted")
    sns.despine()
    if output is not None:
        barplot.figure.savefig(output)
    else:
        sns.plt.show()


