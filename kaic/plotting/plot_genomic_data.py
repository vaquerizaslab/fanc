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
                                            start=max(1, center_node.start-half_window-bin_size),
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
    tsplot = sns.tsplot(data=df, time="distance", unit="region", condition="type", value="contacts",
                        estimator=np.median, err_style="unit_traces", err_palette="Reds")

    if output is not None:
        tsplot.figure.savefig(output)
    else:
        sns.plt.show()

    return df


def hic_matrix_plot(hic, output=None, key=slice(None, None, None), zrange=(5, 40), colormap='viridis'):
    hm = hic[key, key]

    heatmap = sns.heatmap(hm, vmin=zrange[0], vmax=zrange[1], cmap=colormap)

    if output is not None:
        heatmap.figure.savefig(output)
    else:
        sns.plt.show()
