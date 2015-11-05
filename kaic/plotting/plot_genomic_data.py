from __future__ import division
import matplotlib as mp
import seaborn as sns
import pandas
import numpy as np
import kaic.data.general as general
import kaic.data.genomic as genomic
from mpl_toolkits.mplot3d import Axes3D
# import os.path


def _prepare_backend(output):
    if output is not None:
        old_backend = sns.plt.get_backend()
        # extension = os.path.splitext(output)[1]
        sns.plt.switch_backend('pdf')
        sns.plt.ioff()
        return old_backend
    return None

def _plot_figure(figure, output, old_backend):
    if output is not None:
        figure.figure.savefig(output)
        sns.plt.close(figure.figure)
        sns.plt.ion()
        sns.plt.switch_backend(old_backend)
    else:
        sns.plt.show()


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


def _matrix_plot(hm, output=None, lower_percentile=25.0, upper_percentile=98.0,
                 lower=None, upper=None, colormap='viridis'):

    if lower is None or upper is None:
        percentiles = np.percentile(hm, [lower_percentile, upper_percentile])
        if lower is None:
            lower = percentiles[0]
        if upper is None:
            upper = percentiles[1]

    old_backend = None
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


def hic_matrix_plot(hic, output=None, key=slice(0, None, None),
                    lower_percentile=25.0, upper_percentile=98.0,
                    lower=None, upper=None, colormap='viridis'):
    hm = hic[key, key]

    _matrix_plot(hm, output=output, lower_percentile=lower_percentile,
                 upper_percentile=upper_percentile, lower=lower,
                 upper=upper, colormap=colormap)


def hic_matrix_diff_plot(hic1, hic2, output=None, key=slice(0, None, None),
                         lower_percentile=25.0, upper_percentile=98.0,
                         lower=None, upper=None, colormap='viridis'):
    hm1 = hic1[key, key]
    hm2 = hic2[key, key]
    hm = hm1 - hm2

    _matrix_plot(hm, output=output, lower_percentile=lower_percentile,
                 upper_percentile=upper_percentile, lower=lower,
                 upper=upper, colormap=colormap)


def hic_matrix_ratio_plot(hic1, hic2, output=None, key=slice(0, None, None),
                          lower=-2, upper=2, colormap='RdBu_r', log=True):
    hm1 = hic1[key, key]
    hm2 = hic2[key, key]
    if log:
        hm = np.log2(hm1/hm2)
    else:
        hm = hm1/hm2

    _matrix_plot(hm, output=output, lower=lower,
                 upper=upper, colormap=colormap)


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


def hic_ma_plot(hic1, hic2, output=None, highlights=None, key=slice(0, None, None),
                colormap='RdBu_r', log_x=True, log_y=True, inter_chromosomal=True,
                size_factor_correct=True, plot_3d=False):

    hm1 = hic1[key, key]
    hm2 = hic2[key, key]

    # comparison correction factor
    if size_factor_correct:
        ccf = np.sum(hm1)/np.sum(hm2)
        hm2 *= ccf

    # generate a highlights dictionary
    highlights_dict = None
    if highlights is not None:
        highlights_dict = dict()
        for highlight in highlights:
            if highlight.chromosome not in highlights_dict:
                highlights_dict[highlight.chromosome] = []
            highlights_dict[highlight.chromosome].append(highlight)

    main_region_tuples = []
    main_contacts_1 = []
    main_contacts_2 = []
    main_fold_changes = []
    main_colors = []
    highlighted_region_tuples = []
    highlighted_contacts_1 = []
    highlighted_contacts_2 = []
    highlighted_fold_changes = []
    highlighted_colors = []
    if len(hm1.shape) == 2:
        for i in xrange(0, hm1.shape[0]):
            row_region = hm1.row_regions[i]
            for j in xrange(i, hm1.shape[1]):
                col_region = hm1.col_regions[j]
                if hm1[i, j] == 0 or hm2[i, j] == 0:
                    continue

                if not inter_chromosomal and row_region.chromosome != col_region.chromosome:
                    continue

                fold_change = hm1[i, j]/(hm2[i, j])
                if log_y:
                    fold_change = np.log2(fold_change)

                c = 'lightgray'
                if highlights_dict is not None:
                    overlaps_row = False
                    overlaps_col = False
                    if row_region.chromosome in highlights_dict:
                        for highlight in highlights_dict[row_region.chromosome]:
                            if row_region.overlaps(highlight):
                                overlaps_row = True
                            if col_region.overlaps(highlight):
                                overlaps_col = True
                        if row_region.chromosome != col_region.chromosome:
                            for highlight in highlights_dict[col_region.chromosome]:
                                if row_region.overlaps(highlight):
                                    overlaps_row = True
                                if col_region.overlaps(highlight):
                                    overlaps_col = True

                    if overlaps_row and overlaps_col:
                        c = 'orangered'
                    elif overlaps_col or overlaps_row:
                        c = 'khaki'

                if c == 'lightgray':
                    main_contacts_1.append(hm1[i, j])
                    main_contacts_2.append(hm2[i, j])
                    main_fold_changes.append(fold_change)
                    main_region_tuples.append((i, j))
                    main_colors.append(c)
                else:
                    highlighted_contacts_1.append(hm1[i, j])
                    highlighted_contacts_2.append(hm2[i, j])
                    highlighted_fold_changes.append(fold_change)
                    highlighted_region_tuples.append((i, j))
                    highlighted_colors.append(c)

    old_backend = _prepare_backend(output)
    scatter = sns.plt.figure()

    if not plot_3d:
        ax = scatter.add_subplot(111)
        mean_contacts = (np.array(main_contacts_1)+np.array(main_contacts_2))/2
        if log_x:
            mean_contacts = np.log10(mean_contacts)
        ax.scatter(mean_contacts, main_fold_changes, c=main_colors)
        if len(highlighted_region_tuples) > 0:
            h_mean_contacts = (np.array(highlighted_contacts_1)+np.array(highlighted_contacts_2))/2
            if log_x:
                h_mean_contacts = np.log10(h_mean_contacts)
            ax.scatter(h_mean_contacts,
                       highlighted_fold_changes, c=highlighted_colors)
    else:
        ax = scatter.add_subplot(111, projection='3d')
        if log_x:
            ax.scatter(np.log10(main_contacts_1), np.log10(main_contacts_2),
                       main_fold_changes, c=main_colors, alpha=0.3)
        else:
            ax.scatter(main_contacts_1, main_contacts_2, main_fold_changes,
                       c=main_colors, alpha=0.3)
        if len(highlighted_region_tuples) > 0:
            if log_x:
                ax.scatter(np.log10(highlighted_contacts_1), np.log10(highlighted_contacts_2),
                           highlighted_fold_changes, c=highlighted_colors)
            else:
                ax.scatter(highlighted_contacts_1, highlighted_contacts_2,
                           highlighted_fold_changes, c=highlighted_colors)

    _plot_figure(scatter, output, old_backend)

    return None


def hic_marginals_plot(hic, output=None):
    marginals = hic.marginals()

    old_backend = _prepare_backend(output)
    scatter = sns.plt.plot(marginals)
    _plot_figure(scatter, output, old_backend)


def hic_triangle_plot(hic, key=slice(0, None, None), output=None, colormap='viridis', max_height=None, axes=None,
                      n_x_labels=10, vmin=1, vmax=50):
    hm = hic[key, key]
    print hm.shape
    n = hm.shape[0]

    # mask areas we don't want to plot
    # lower triangle
    mask_lower = np.tril_indices(n, k=-1)
    hm[mask_lower] = np.nan
    if max_height is not None:
        # upper right corner
        mask_upper = np.triu_indices(n, k=max_height)
        hm[mask_upper] = np.nan
    triangle = np.ma.masked_array(hm, np.isnan(hm))

    # prepare an array of tuples that will be used to rotate triangle
    A = np.array([(y, x) for x in range(n, -1, -1) for y in range(n + 1)])
    # rotation matrix 45 degrees
    t = np.array([[0.707, 0.707], [-0.707, 0.707]])
    # "rotate" A
    A = np.dot(A, t)
    # transform A into x and y values
    X = A[:, 1].reshape((n + 1, n + 1))
    Y = A[:, 0].reshape((n + 1, n + 1))

    # flip triangle (because that pcolormesh works opposite than imshow)
    flip_triangle = np.flipud(triangle)

    with sns.axes_style("ticks"):
        new_plot = False
        if axes is None:
            new_plot = True
            fig = sns.plt.figure(figsize=(6, 3.5))
            axes = fig.add_subplot(111, frame_on=True)

        # normalize colors
        cmap = mp.cm.get_cmap(colormap)
        norm = mp.colors.BoundaryNorm(np.linspace(vmin, vmax, 50), cmap.N)

        # create plot
        caxes = sns.plt.pcolormesh(X, Y, flip_triangle, axes=axes, cmap=cmap, norm=norm)

        # re-calculate and reset axis limits
        max_x = max(A[:, 1])
        max_y = 0
        for i in xrange(flip_triangle.shape[0]):
            for j in xrange(flip_triangle.shape[1]):
                if flip_triangle[i, j] is not np.ma.masked:
                    max_y = max(max_y, Y[i, j]+2)
        axes.set_ylim((-1, max_y))
        axes.set_xlim((0, max_x))

        # set ticks to genomic regions
        all_x_ticks = np.linspace(0, max_x, n+1)
        last_label = "%s: %d" % (hm.row_regions[-1].chromosome, hm.row_regions[-1].end)
        all_x_labels = ["%s: %d" % (region.chromosome, region.start) for region in hm.row_regions] + [last_label]
        labels_dist = int(round(len(all_x_labels)/n_x_labels))
        axes.set_xticks(all_x_ticks[0:len(all_x_ticks):labels_dist])
        axes.set_xticklabels(all_x_labels[0:len(all_x_ticks):labels_dist], rotation=45, ha='right')

        # remove y ticks
        axes.set_yticks([])

        # Hide the left, right and top spines
        sns.despine(left=True)
        # hide background patch
        axes.patch.set_visible(False)

        # Only show ticks on the left and bottom spines
        axes.xaxis.set_ticks_position('bottom')

        # make figure margins accommodate labels
        sns.plt.tight_layout()
        if new_plot and output is None:
            sns.plt.show()

