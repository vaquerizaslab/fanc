import seaborn as sns
import pandas
import numpy as np
import kaic.data.general as general
import kaic.data.genomic as genomic
import kaic.plotting.colormaps as cmaps

sns.plt.register_cmap(name='viridis', cmap=cmaps.viridis)
sns.plt.register_cmap(name='plasma', cmap=cmaps.plasma)
sns.plt.register_cmap(name='inferno', cmap=cmaps.inferno)
sns.plt.register_cmap(name='magma', cmap=cmaps.magma)


def hic_contact_plot_linear(hic, regions, output=None, window_size=1000000):
    if isinstance(regions, general.Table):
        new_regions = []
        for region in regions:
            region_string = "%s:%d-%d" % (region.chrom,
                                          region.start,
                                          region.end)
            new_regions.append(region_string)
        regions = new_regions

    contact_list = []
    half_window = int(window_size/2)
    bin_size = hic.bin_size()
    for i, region_string in enumerate(regions):
        feature_region = genomic.GenomicRegion.from_string(region_string)
        center = feature_region.start + int((feature_region.end-feature_region.start)/2)
        center_region = genomic.GenomicRegion(chromosome=feature_region.chromosome,
                                              start=center, end=center)

        center_node = hic.get_node(center_region)

        left_region = genomic.GenomicRegion(chromosome=feature_region.chromosome,
                                            start=max(1, center_node.start-half_window),
                                            end=center_node.start)

        right_region = genomic.GenomicRegion(chromosome=feature_region.chromosome,
                                             start=center_node.end+1,
                                             end=center_node.end+half_window)

        print left_region
        print right_region

        hic_left = hic[center_region, left_region][0]
        for j in xrange(0, len(hic_left)):
            j_r = len(hic_left)-j-1
            label = -1*bin_size*j
            val = hic_left[j_r]
            contact_list.append([label, val, str(i), 'data'])

        hic_right = hic[center_region, right_region][0]
        for j in xrange(0, len(hic_right)):
            label = bin_size*(j+1)
            val = hic_right[j]
            contact_list.append([label, val, str(i), 'data'])

    df = pandas.DataFrame(contact_list, columns=["distance", "contacts", "region", "type"])

    if output is not None:
        old_backend = sns.plt.get_backend()
        sns.plt.switch_backend('pdf')
        sns.plt.ioff()

    tsplot = sns.tsplot(data=df, time="distance", unit="region", condition="type", value="contacts",
                        estimator=np.median, err_style="unit_traces", err_palette="Reds")

    if output is not None:
        tsplot.figure.savefig(output)
        sns.plt.close(tsplot.figure)
        sns.plt.ion()
        sns.plt.switch_backend(old_backend)
    else:
        sns.plt.show()

    return df


def hic_matrix_plot(hic, output=None, key=slice(None, None, None),
                    lower_percentile=25.0, upper_percentile=98.0,
                    lower=None, upper=None, colormap='viridis'):
    hm = hic[key, key]

    if lower is None or upper is None:
        percentiles = np.percentile(hm, [lower_percentile, upper_percentile])
        if lower is None:
            lower = percentiles[0]
        if upper is None:
            upper = percentiles[1]

    if output is not None:
        old_backend = sns.plt.get_backend()
        sns.plt.switch_backend('pdf')
        sns.plt.ioff()

    heatmap = sns.heatmap(hm, vmin=lower, vmax=upper, cmap=colormap,
                          square=True, xticklabels=False, yticklabels=False)

    if output is not None:
        heatmap.figure.savefig(output)
        sns.plt.close(heatmap.figure)
        sns.plt.ion()
        sns.plt.switch_backend(old_backend)
    else:
        sns.plt.show()


def _correlation_df(hic1, hic2, include_zeros=False, in_percent=False):
    chromosomes = hic1.chromosomes()

    corr_matrix = np.zeros(shape=(len(chromosomes), len(chromosomes)))

    for chr_i in xrange(0, len(chromosomes)):
        chromosome1 = chromosomes[chr_i]
        for chr_j in xrange(chr_i, len(chromosomes)):
            chromosome2 = chromosomes[chr_j]

            m1 = hic1[chromosome1, chromosome2]
            m2 = hic2[chromosome1, chromosome2]

            contacts1 = []
            contacts2 = []
            for i in xrange(0, m1.shape[0]):
                for j in xrange(i, m1.shape[1]):
                    if not include_zeros and m1[i, j] == 0 and m2[i, j] == 0:
                        continue
                    contacts1.append(m1[i, j])
                    contacts2.append(m2[i, j])

            corr = np.corrcoef(contacts1, contacts2)[0, 1]
            if in_percent:
                corr = int(round(corr*100))
            corr_matrix[chr_i, chr_j] = corr
            corr_matrix[chr_j, chr_i] = corr

    return pandas.DataFrame(data=corr_matrix, index=chromosomes, columns=chromosomes)


def hic_correlation_plot(hic1, hic2, output=None, include_zeros=False, colormap="viridis", size=10):
    corr_df = _correlation_df(hic1, hic2, include_zeros=include_zeros, in_percent=True)

    if output is not None:
        old_backend = sns.plt.get_backend()
        sns.plt.switch_backend('pdf')
        sns.plt.ioff()

    sns.plt.figure(figsize=(size, size))
    heatmap = sns.heatmap(corr_df, vmin=5, vmax=95, cmap=colormap, square=True, annot=True, fmt=".0f")

    if output is not None:
        heatmap.figure.savefig(output)
        sns.plt.close(heatmap.figure)
        sns.plt.ion()
        sns.plt.switch_backend(old_backend)
    else:
        sns.plt.show()
