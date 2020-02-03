from __future__ import division
from fanc.config import config
import seaborn as sns
import pandas
import numpy as np
from future.utils import string_types
from genomic_regions import GenomicRegion
from ..hic import LowCoverageFilter
import matplotlib.pyplot as plt


def _prepare_backend(output):
    if output is not None:
        old_backend = plt.get_backend()
        # extension = os.path.splitext(output)[1]
        sns.set_style("ticks")
        plt.switch_backend('pdf')
        plt.ioff()
        return old_backend
    return None


def _plot_figure(figure, output, old_backend):
    sns.despine()
    if output is not None:
        figure.savefig(output)
        plt.close(figure)
        plt.ion()
        plt.switch_backend(old_backend)
    else:
        plt.show()


def hic_contact_plot_linear(hic, regions, output=None, window_size=1000000):
    contact_list = []
    half_window = int(window_size/2)
    bin_size = hic.bin_size
    for i, feature_region in enumerate(regions):
        if isinstance(feature_region, string_types):
            feature_region = GenomicRegion.from_string(feature_region)
        center = feature_region.start + int((feature_region.end-feature_region.start)/2)
        center_region = GenomicRegion(chromosome=feature_region.chromosome,
                                      start=center, end=center)

        center_node = hic.get_node(center_region)

        left_region = GenomicRegion(chromosome=feature_region.chromosome,
                                    start=max(1, center_node.start-half_window),
                                    end=center_node.start)

        right_region = GenomicRegion(chromosome=feature_region.chromosome,
                                     start=center_node.end+1,
                                     end=center_node.end+half_window)

        hic_left = hic[center_region, left_region][0]
        for j in range(0, len(hic_left)):
            j_r = len(hic_left)-j-1
            label = -1*bin_size*j
            val = hic_left[j_r]
            contact_list.append([label, val, str(i), 'data'])

        hic_right = hic[center_region, right_region][0]
        for j in range(0, len(hic_right)):
            label = bin_size*(j+1)
            val = hic_right[j]
            contact_list.append([label, val, str(i), 'data'])

    df = pandas.DataFrame(contact_list, columns=["distance", "contacts", "region", "type"])

    if output is not None:
        old_backend = plt.get_backend()
        plt.switch_backend('pdf')
        plt.ioff()

    tsplot = sns.tsplot(data=df, time="distance", unit="region", condition="type", value="contacts",
                        estimator=np.median, err_style="unit_traces", err_palette="Reds")

    if output is not None:
        tsplot.figure.savefig(output)
        plt.close(tsplot.figure)
        plt.ion()
        plt.switch_backend(old_backend)
    else:
        plt.show()

    return df


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
        for i in range(0, hm1.shape[0]):
            row_region = hm1.row_regions[i]
            for j in range(i, hm1.shape[1]):
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
    scatter = plt.figure()

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


def hic_marginals_plot(hic, output=None, lower=None, upper=None, rel_cutoff=0.1):
    marginals = hic.marginals()

    if lower is None or upper is None:
        f = LowCoverageFilter(hic)
        lower_calc, upper_calc = f.calculate_cutoffs(rel_cutoff)
        if lower is None:
            lower = lower_calc
        if upper is None:
            upper = upper_calc

    old_backend = _prepare_backend(output)
    fig, ax = plt.subplots()
    ax.plot(marginals)
    ax.axhline(upper, color='r', linestyle=':')
    ax.axhline(lower, color='r', linestyle=':')
    _plot_figure(fig, output, old_backend)
