import numpy as np
from bisect import bisect_left


def vector_enrichment_profile(matrix, vector, mappable=None, per_chromosome=True,
                              percentiles=(20.0, 40.0, 60.0, 80.0, 100.0),
                              symmetric_at=None, exclude_chromosomes=()):
    if len(exclude_chromosomes) > 0:
        chromosome_bins = matrix.chromosome_bins
        exclude_vector = []
        for chromosome in matrix.chromosomes():
            if chromosome not in exclude_chromosomes:
                b = chromosome_bins[chromosome]
                for v in vector[b[0]:b[1]]:
                    exclude_vector.append(v)
                # exclude_vector += vector[b[0]:b[1]]
    else:
        exclude_vector = vector
    exclude_vector = np.array(exclude_vector)

    if symmetric_at is not None:
        lv = exclude_vector[exclude_vector <= symmetric_at]
        gv = exclude_vector[exclude_vector > symmetric_at]
        lv_cutoffs = np.nanpercentile(lv, percentiles)
        gv_cutoffs = np.nanpercentile(gv, percentiles)
        bin_cutoffs = np.concatenate((lv_cutoffs, gv_cutoffs))
    else:
        bin_cutoffs = np.nanpercentile(exclude_vector, percentiles)

    bins = []
    for value in vector:
        bins.append(bisect_left(bin_cutoffs, value))

    s = len(bin_cutoffs)
    m = np.zeros((s, s))
    c = np.zeros((s, s))

    if mappable is None:
        mappable = matrix.mappable()

    if per_chromosome:
        for chromosome in matrix.chromosomes():
            if chromosome in exclude_chromosomes:
                continue
            oem = matrix.matrix((chromosome, chromosome), oe=True, oe_per_chromosome=True)
            for i, row_region in enumerate(oem.row_regions):
                if not mappable[row_region.ix]:
                    continue
                i_bin = s - bins[row_region.ix] - 1
                for j, col_region in enumerate(oem.col_regions):
                    if not mappable[col_region.ix]:
                        continue
                    j_bin = s - bins[col_region.ix] - 1
                    value = oem[i, j]

                    m[i_bin, j_bin] += value
                    c[i_bin, j_bin] += 1
                    m[j_bin, i_bin] += value
                    c[j_bin, i_bin] += 1
    else:
        oem = matrix.matrix(oe=True)
        for i in range(oem.shape):
            if not mappable[i]:
                continue
            i_bin = s - bins[i] - 1
            for j in range(i, oem.shape):
                if not mappable[j]:
                    continue
                j_bin = s - bins[j] - 1
                value = oem[i, j]

                m[i_bin, j_bin] += value
                c[i_bin, j_bin] += 1
                m[j_bin, i_bin] += value
                c[j_bin, i_bin] += 1

    m /= c
    # m[c == 0] = 0
    rev_cutoffs = bin_cutoffs[::-1]
    for i in range(len(rev_cutoffs) - 1, 0, -1):
        if np.isclose(rev_cutoffs[i - 1], rev_cutoffs[i]):
            m[:, i - 1] = m[:, i]
            m[i - 1, :] = m[i, :]

    return np.log2(m), rev_cutoffs

