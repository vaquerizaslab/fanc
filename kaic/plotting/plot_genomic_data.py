import seaborn as sns
import kaic.data.general as general
import kaic.data.genomic as genomic


def hic_contact_plot_linear(hic, regions, window_size=1000000):
    if isinstance(regions, general.Table):
        new_regions = []
        for region in regions:
            region_string = "%s:%d-%d" % (region.chrom,
                                          region.start,
                                          region.end)
            new_regions.append(region_string)
        regions = new_regions

    half_window = int(window_size/2)
    bin_size = hic.bin_size()
    for region_string in regions:
        nodes = hic.get_nodes(region_string)
        start_node = nodes[0]
        end_node = nodes[-1]

        left_region = genomic.GenomicRegion(chromosome=start_node.chromosome,
                                            start=max(1, start_node.start-half_window),
                                            end=start_node.start-1)

        right_region = genomic.GenomicRegion(chromosome=end_node.chromosome,
                                             start=end_node.end+1,
                                             end=end_node.end+half_window)


def hic_matrix_plot(hic, output=None, key=slice(None, None, None), zrange=[5, 40], colormap='afmhot_r'):
    hm = hic[key, key]

    heatmap = sns.heatmap(hm, vmin=zrange[0], vmax=zrange[1], cmap=colormap)

    if output is not None:
        heatmap.figure.savefig(output)
    else:
        sns.plt.show()
