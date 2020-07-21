import numpy as np
from bisect import bisect_left
from collections import defaultdict
import intervaltree
from genomic_regions import GenomicRegion
from future.utils import string_types


def vector_enrichment_profile(matrix, vector, mappable=None, per_chromosome=True,
                              percentiles=(20.0, 40.0, 60.0, 80.0, 100.0),
                              symmetric_at=None, exclude_chromosomes=(),
                              intra_chromosomal=True, inter_chromosomal=False,
                              collapse_identical_breakpoints=False):
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

    if collapse_identical_breakpoints:
        new_bin_cutoffs = []
        for i in range(1, len(bin_cutoffs)):
            if not np.isclose(bin_cutoffs[i], bin_cutoffs[i-1]):
                new_bin_cutoffs.append(bin_cutoffs[i])
        bin_cutoffs = new_bin_cutoffs

    if isinstance(intra_chromosomal, bool):
        if intra_chromosomal:
            intra_chromosomal = 0
        else:
            intra_chromosomal = -1

    bins = []
    for value in vector:
        bins.append(bisect_left(bin_cutoffs, value))

    s = len(bin_cutoffs)
    m = np.zeros((s, s))
    c = np.zeros((s, s))

    bin_size = matrix.bin_size
    if mappable is None:
        mappable = matrix.mappable()

    if per_chromosome:
        chromosomes = matrix.chromosomes()
        for chr1_ix in range(len(chromosomes)):
            chromosome1 = chromosomes[chr1_ix]
            if chromosome1 in exclude_chromosomes:
                continue

            for chr2_ix in range(chr1_ix, len(chromosomes)):
                chromosome2 = chromosomes[chr2_ix]
                if chromosome2 in exclude_chromosomes:
                    continue

                is_intra_chromosomal = chromosome1 == chromosome2

                if not inter_chromosomal and not is_intra_chromosomal:
                    continue

                oem = matrix.matrix((chromosome1, chromosome2),
                                    oe=True, oe_per_chromosome=True)
                for i, row_region in enumerate(oem.row_regions):
                    if not mappable[row_region.ix]:
                        continue
                    i_bin = s - bins[row_region.ix] - 1
                    for j, col_region in enumerate(oem.col_regions):
                        if not mappable[col_region.ix]:
                            continue

                        if is_intra_chromosomal and abs(j - i) * bin_size < intra_chromosomal:
                            continue

                        j_bin = s - bins[col_region.ix] - 1
                        value = oem[i, j]

                        m[i_bin, j_bin] += value
                        c[i_bin, j_bin] += 1
                        m[j_bin, i_bin] += value
                        c[j_bin, i_bin] += 1
    else:
        oem = matrix.matrix(oe=True)
        for i in range(oem.shape[0]):
            if not mappable[i]:
                continue
            i_bin = s - bins[i] - 1
            for j in range(i, oem.shape[1]):
                if not mappable[j]:
                    continue

                is_intra_chromosomal = oem.row_regions[i].chromosome == oem.col_regions[j].chromosome

                if is_intra_chromosomal:
                    if abs(j - i) * bin_size < intra_chromosomal:
                        continue
                elif not inter_chromosomal:
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


class RegionScoreMatrix(np.ma.MaskedArray):
    def __new__(cls, input_matrix, row_regions=None, columns=None,
                mask=None, *args, **kwargs):
        obj = np.asarray(input_matrix).view(cls, *args, **kwargs)
        obj._row_region_trees = None
        obj.row_regions = None
        obj.columns = columns
        obj.set_row_regions(row_regions)
        obj.mask = mask
        return obj

    def _interval_tree_regions(self, regions):
        intervals = defaultdict(list)
        for i, region in enumerate(regions):
            interval = intervaltree.Interval(region.start - 1, region.end,
                                             data=i)
            intervals[region.chromosome].append(interval)

        interval_trees = {chromosome: intervaltree.IntervalTree(intervals)
                          for chromosome, intervals in intervals.items()}
        return interval_trees

    def set_row_regions(self, regions):
        self.row_regions = regions
        if regions is not None:
            self._row_region_trees = self._interval_tree_regions(regions)
        else:
            self._row_region_trees = None

    def __array_finalize__(self, obj):
        if isinstance(self, np.ma.core.MaskedArray):
            np.ma.MaskedArray.__array_finalize__(self, obj)

        if obj is None:
            return

        self.set_row_regions(getattr(obj, 'row_regions', None))
        self.columns = getattr(obj, 'columns', None)

    def __setitem__(self, key, item):
        self._setitem = True
        try:
            if isinstance(self, np.ma.core.MaskedArray):
                np.ma.MaskedArray.__setitem__(self, key, item)
            else:
                np.ndarray.__setitem__(self, key, item)
        finally:
            self._setitem = False

    def __getitem__(self, index):
        self._getitem = True

        try:
            all_cols = slice(0, len(self.columns), 1)
        except (AttributeError, TypeError):
            try:
                all_cols = slice(0, self.shape[1], 1)
            except IndexError:
                all_cols = None

        # convert string types into region indexes
        if isinstance(index, tuple):
            if len(index) == 2:
                row_key = self._convert_key(
                    index[0],
                    self._row_region_trees if hasattr(self, '_row_region_trees') else None
                )
                col_key = index[1]
                index = (row_key, col_key)
            elif len(index) == 1:
                row_key = self._convert_key(index[0],
                                            self._row_region_trees
                                            if hasattr(self, '_row_region_trees') else None)

                col_key = all_cols
                index = (row_key, )
            else:
                col_key = all_cols
                row_key = index
                index = row_key
        else:
            row_key = self._convert_key(index,
                                        self._row_region_trees
                                        if hasattr(self, '_row_region_trees') else None)

            col_key = all_cols
            index = row_key

        try:
            if isinstance(self, np.ma.core.MaskedArray):
                out = np.ma.MaskedArray.__getitem__(self, index)
            else:
                out = np.ndarray.__getitem__(self, index)
        finally:
            self._getitem = False

        if not isinstance(out, np.ndarray):
            return out

        # get regions
        try:
            row_regions = self.row_regions[row_key]
        except TypeError:
            row_regions = None

        try:
            columns = self.columns[col_key]
        except TypeError:
            columns = None

        if isinstance(row_regions, GenomicRegion):
            out.row_regions = [row_regions]
        else:
            out.row_regions = row_regions

        if not isinstance(columns, list):
            out.columns = [columns]
        else:
            out.columns = columns

        return out

    def __getslice__(self, start, stop):
        return self.__getitem__(slice(start, stop))

    def _convert_key(self, key, region_trees):
        if isinstance(key, string_types):
            key = GenomicRegion.from_string(key)

        if isinstance(key, GenomicRegion):
            start = None
            stop = None
            try:
                key_start = 0 if key.start is None else max(0, key.start - 1)
                key_end = key.end
                for interval in region_trees[key.chromosome][key_start:key_end]:
                    i = interval.data
                    start = min(i, start) if start is not None else i
                    stop = max(i + 1, stop) if stop is not None else i + 1
            except KeyError:
                raise ValueError("Requested chromosome {} was not "
                                 "found in this matrix.".format(key.chromosome))

            if start is None or stop is None:
                raise ValueError("Requested region {} was not found in this matrix.".format(key))

            return slice(start, stop, 1)
        return key