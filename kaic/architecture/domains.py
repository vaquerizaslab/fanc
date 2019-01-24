from genomic_regions import GenomicRegion
from .maxima_callers import MaximaCallerDelta
from ..regions import RegionsTable
from ..tools.matrix import apply_sliding_func, trim_stats
from scipy.stats.mstats import gmean
import numpy as np
import tables
import itertools
import warnings
import logging


logger = logging.getLogger(__name__)


class RegionScoreTable(RegionsTable):
    def __init__(*args, **kwargs):
        additional_fields = kwargs.get('additional_fields', None)
        if additional_fields is None:
            additional_fields = {}
        additional_fields['score'] = tables.Float32Col(pos=0)
        kwargs['additional_fields'] = additional_fields

        RegionsTable.__init__(*args, **kwargs)


class InsulationScore(RegionScoreTable):
    def __init__(*args, **kwargs):
        RegionScoreTable.__init__(*args, **kwargs)

    @classmethod
    def from_hic(cls, hic, window_size_in_bp, window_offset=0,
                 file_name=None, tmpdir=None, impute_missing=False,
                 na_threshold=0.5, normalise=True, normalisation_window=None,
                 trim_mean_proportion=0.0, geometric_mean=False,
                 subtract_mean=False, log=True):
        window_size = hic.distance_to_bins(window_size_in_bp)
        insulation_score_regions = cls(file_name=file_name, tmpdir=tmpdir)
        insulation_score_regions.add_regions(hic.regions, preserve_attributes=False)

        if geometric_mean:
            avg_stat = gmean
        else:
            avg_stat = np.nanmean

        chromosome_bins = hic.chromosome_bins
        window_offset += 1

        mappable = hic.mappable()

        if impute_missing:
            intra_expected, intra_expected_chromosome, _ = hic.expected_values()
        else:
            intra_expected, intra_expected_chromosome = None, None

        def _pair_to_bins(source, sink):
            # correct index by chromosome and window offset
            i = source - chromosome_start + window_offset
            j = sink - chromosome_start - window_offset

            start = max(i, j - window_size + 1)
            stop = min(j + 1, i + window_size)
            for ii_bin in range(start, stop):
                yield ii_bin

        ii_list = []
        chromosomes = hic.chromosomes()
        for chromosome in chromosomes:
            chromosome_start, chromosome_stop = chromosome_bins[chromosome]
            values_by_chromosome = [0 for _ in range(chromosome_start, chromosome_stop)]

            # add each edge weight to every insulation window that contains it
            for edge in hic.edges((chromosome, chromosome), lazy=True):
                for ii_bin in _pair_to_bins(edge.source, edge.sink):
                    weight = getattr(edge, hic._default_score_field)
                    values_by_chromosome[ii_bin] += weight

            # add imputed values, if requested
            if impute_missing:
                expected = intra_expected_chromosome[chromosome]

                covered = set()
                for ix in range(chromosome_start, chromosome_stop):
                    if not mappable[ix]:
                        sink = ix
                        for source in range(max(chromosome_start, ix - window_size - window_offset - 2), ix):
                            if (source, sink) in covered:
                                continue
                            covered.add((source, sink))
                            weight = expected[sink - source]
                            for ii_bin in _pair_to_bins(source, sink):
                                values_by_chromosome[ii_bin] += weight

                        source = ix
                        for sink in range(ix, min(chromosome_stop, ix + window_offset + window_size + 2)):
                            if (source, sink) in covered:
                                continue
                            covered.add((source, sink))
                            weight = expected[sink - source]
                            for ii_bin in _pair_to_bins(source, sink):
                                values_by_chromosome[ii_bin] += weight

                for k in range(len(values_by_chromosome)):
                    if (k - window_offset < window_size - 1
                            or k + window_offset > len(values_by_chromosome) - window_size):
                        values_by_chromosome[k] = np.nan
                    else:
                        values_by_chromosome[k] /= window_size ** 2
            # count unmappable bins in every window
            else:
                unmappable_horizontal = [0 for _ in range(chromosome_start, chromosome_stop)]
                unmappable_vertical = [0 for _ in range(chromosome_start, chromosome_stop)]

                for ix in range(chromosome_start, chromosome_stop):
                    if not mappable[ix]:
                        # horizontal
                        start_bin = ix - chromosome_start + window_offset
                        for ii_bin in range(start_bin, start_bin + window_size):
                            if 0 <= ii_bin < len(unmappable_horizontal):
                                unmappable_horizontal[ii_bin] += 1
                        # vertical
                        start_bin = ix - chromosome_start - window_offset
                        for ii_bin in range(start_bin - window_size + 1, start_bin + 1):
                            if 0 <= ii_bin < len(unmappable_vertical):
                                unmappable_vertical[ii_bin] += 1

                for k in range(len(values_by_chromosome)):
                    na_vertical = unmappable_vertical[k] * window_size
                    na_horizontal = unmappable_horizontal[k] * window_size
                    na_overlap = unmappable_horizontal[k] * unmappable_vertical[k]
                    na_total = na_vertical + na_horizontal - na_overlap

                    # take into account nan values when adding zeros
                    if ((na_total > (window_size ** 2 * na_threshold))
                            or k - window_offset < window_size - 1
                            or k + window_offset > len(values_by_chromosome) - window_size):
                        values_by_chromosome[k] = np.nan
                    else:
                        values_by_chromosome[k] /= (window_size ** 2 - na_total)

            ii_by_chromosome = values_by_chromosome

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)

                ii_by_chromosome = np.array(ii_by_chromosome)
                if normalise:
                    logger.debug("Normalising insulation index")
                    if normalisation_window is not None:
                        logger.debug("Sliding window average")
                        mean_ins = apply_sliding_func(ii_by_chromosome, normalisation_window,
                                                      func=lambda x: trim_stats(x[np.isfinite(x)],
                                                                                trim_mean_proportion,
                                                                                stat=avg_stat))
                    else:
                        logger.debug("Whole chromosome mean")
                        mean_ins = trim_stats(ii_by_chromosome[np.isfinite(ii_by_chromosome)],
                                              trim_mean_proportion, stat=avg_stat)

                    if not subtract_mean:
                        logger.debug("Dividing by mean")
                        ii_by_chromosome = ii_by_chromosome / mean_ins
                    else:
                        logger.debug("Subtracting mean")
                        ii_by_chromosome = ii_by_chromosome - mean_ins
            ii_list.append(ii_by_chromosome)

        ins_matrix = np.array(list(itertools.chain.from_iterable(ii_list)))

        if log:
            logger.debug("Log-transforming insulation index")
            ins_matrix = np.log2(ins_matrix)

        insulation_score_regions.region_data('score', ins_matrix)
        return insulation_score_regions

    def boundaries(self, min_score=None, delta_window=3, log=False,
                   sub_bin_precision=False, call_maxima=False):
        """
        Call insulation boundaries based on minima in an insulation vector of this object.

        :param min_score: Minimum difference between minimum and the closest maximum
                          in the insulation vector for a region to be considered a
                          boundary
        :param delta_window: Window size in bins to control smoothing of the delta function
                             used to calculate the derivative of the insulation index.
                             Calculation takes into account d bins upstream and d
                             bins downstream for a total window size of 2*d + 1 bins.
        :param log: Log2-transform insulation index before boundary calls
        :param sub_bin_precision: Call boundaries with sub bin precision, by taking
                                  into account the precise zero transition of the delta vector.
        :param call_maxima: Call maxima instead of minima as boundaries
        :return: list of :class:`~kaic.data.genomic.GenomicRegion`
        """
        index = self.region_data('score')
        if log:
            index = np.log2(index)
        peaks = MaximaCallerDelta(index, window_size=delta_window, sub_bin_precision=sub_bin_precision)
        if call_maxima:
            minima, scores = peaks.get_maxima()
        else:
            minima, scores = peaks.get_minima()
        regions = list(self.regions)

        boundaries = []
        for i, ix in enumerate(minima):
            if min_score is not None and scores[i] < min_score:
                continue
            if sub_bin_precision:
                region = regions[int(ix)]
                fraction = ix % 1
                shift = int((fraction - .5)*(region.end - region.start))
                b = GenomicRegion(chromosome=region.chromosome,
                                  start=region.start + shift,
                                  end=region.end + shift,
                                  score=scores[i])
            else:
                region = regions[ix]
                b = GenomicRegion(chromosome=region.chromosome,
                                  start=region.start,
                                  end=region.end,
                                  score=scores[i])
            boundaries.append(b)

        logger.info("Found {} boundaries".format(len(boundaries)))
        return boundaries
