from genomic_regions import GenomicRegion
from .maxima_callers import MaximaCallerDelta
from .helpers import RegionScoreMatrix
from ..regions import RegionsTable
from ..tools.matrix import apply_sliding_func, trim_stats
from ..tools.files import write_bed, write_bigwig, write_gff
from scipy.stats.mstats import gmean
import numpy as np
import tables
import itertools
import warnings
import logging


logger = logging.getLogger(__name__)


class RegionScoreTable(RegionsTable):

    _classid = 'REGIONSCORETABLE'

    def __init__(*args, **kwargs):
        additional_fields = kwargs.get('additional_fields', None)
        if additional_fields is None:
            additional_fields = {}
        additional_fields['score'] = tables.Float32Col(pos=0)
        kwargs['additional_fields'] = additional_fields

        RegionsTable.__init__(*args, **kwargs)

    def scores(self):
        return self.region_data('score')


class RegionMultiScoreTable(RegionsTable):

    _classid = 'REGIONMULTISCORETABLE'

    def __init__(self, score_fields=None, *args, **kwargs):
        if score_fields is not None:
            additional_fields = kwargs.get('additional_fields', None)
            if additional_fields is None:
                additional_fields = {}
            for i, score_field in enumerate(score_fields):
                additional_fields[score_field] = tables.Float32Col(pos=i)
            kwargs['additional_fields'] = additional_fields

        RegionsTable.__init__(self, *args, **kwargs)


class RegionScoreParameterTable(RegionMultiScoreTable):

    _classid = 'REGIONSCOREPARAMETERTABLE'

    def __init__(self, parameter_values=None, parameter_prefix='score_', *args, **kwargs):
        if parameter_values is None:
            RegionMultiScoreTable.__init__(self, None, *args, **kwargs)
            prefix = self.meta.parameter_prefix
            self._parameter_prefix = prefix.decode('utf-8') if isinstance(prefix, bytes) else prefix
            self._parameters = self.meta.parameter_values
            self._score_fields = self._score_field_converter()
        else:
            self._parameter_prefix = parameter_prefix
            self._parameters = parameter_values
            self._score_fields = self._score_field_converter(parameter_values)
            RegionMultiScoreTable.__init__(self, self._score_fields, *args, **kwargs)
            self.meta['parameter_values'] = parameter_values
            self.meta['parameter_prefix'] = parameter_prefix

    def _score_field_converter(self, parameters=None):
        if parameters is None:
            parameters = self._parameters
        return [self._parameter_prefix + str(v) for v in parameters]

    def scores(self, parameter, scores=None):
        score_field = self._score_field_converter([parameter])[0]
        return self.region_data(score_field, scores)

    def score_regions(self, parameter):
        score_field = self._score_field_converter([parameter])[0]
        regions = []
        for region in self.regions:
            r = GenomicRegion(chromosome=region.chromosome,
                              start=region.start, end=region.end,
                              strand=region.strand, score=getattr(region, score_field))
            regions.append(r)
        scores = RegionScoreTable()
        scores.add_regions(regions)
        return scores

    def score_matrix(self, region=None, parameters=None):
        scores = []
        regions = list(self.regions(region))
        columns = self._score_field_converter(parameters)
        for r in regions:
            row = [getattr(r, name) for name in columns]
            scores.append(row)
        return RegionScoreMatrix(np.array(scores), row_regions=regions, columns=columns)

    def _to_file(self, file_name, parameter, subset=None, _write_function=write_bed):
        score_field = self._score_field_converter([parameter])[0]

        regions = []
        for region in self.regions(subset):
            score = getattr(region, score_field)
            r = GenomicRegion(chromosome=region.chromosome,
                              start=region.start, end=region.end,
                              score=score)
            regions.append(r)
        _write_function(file_name, regions)
        return file_name

    def to_bed(self, file_name, parameter, subset=None):
        return self._to_file(file_name, parameter, subset=subset, _write_function=write_bed)

    def to_gff(self, file_name, parameter, subset=None):
        return self._to_file(file_name, parameter, subset=subset, _write_function=write_gff)

    def to_bigwig(self, file_name, parameter, subset=None):
        return self._to_file(file_name, parameter, subset=subset, _write_function=write_bigwig)


class InsulationScore(RegionScoreTable):

    _classid = 'INSULATIONSCORE'

    def __init__(*args, **kwargs):
        RegionScoreTable.__init__(*args, **kwargs)

    @classmethod
    def from_hic(cls, hic, window_size_in_bp, window_offset=0,
                 file_name=None, tmpdir=None, impute_missing=False,
                 na_threshold=0.5, normalise=True, normalisation_window=None,
                 trim_mean_proportion=0.0, geometric_mean=False,
                 subtract_mean=False, log=True):
        window_size = int(hic.distance_to_bins(window_size_in_bp) / 2)
        insulation_score_regions = cls(file_name=file_name, mode='w', tmpdir=tmpdir)
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
                        logger.debug("Dividing by mean {}".format(mean_ins))
                        ii_by_chromosome = ii_by_chromosome / mean_ins
                    else:
                        logger.debug("Subtracting mean {}".format(mean_ins))
                        ii_by_chromosome = ii_by_chromosome - mean_ins
            ii_list.append(ii_by_chromosome)

        ins_matrix = np.array(list(itertools.chain.from_iterable(ii_list)))

        if log:
            logger.debug("Log-transforming insulation index")
            ins_matrix = np.log2(ins_matrix)

        insulation_score_regions.region_data('score', ins_matrix)
        return insulation_score_regions

    def boundaries(self, *args, **kwargs):
        return Boundaries.from_insulation_score(self, *args, **kwargs)


class InsulationScores(RegionScoreParameterTable):

    _classid = 'INSULATIONSCORES'

    def __init__(*args, **kwargs):
        RegionScoreParameterTable.__init__(*args, **kwargs)

    @property
    def window_sizes(self):
        return self._parameters

    @classmethod
    def from_hic(cls, hic, window_sizes, file_name=None, tmpdir=None,
                 **kwargs):
        insulation_scores = cls(parameter_prefix='insulation_',
                                parameter_values=list(window_sizes),
                                file_name=file_name, mode='w',
                                tmpdir=tmpdir)

        for i, window_size in enumerate(window_sizes):
            logger.info("Calculating insulation score with window size {}".format(window_size))
            ii = InsulationScore.from_hic(hic, window_size, **kwargs)
            if i == 0:
                insulation_scores.add_regions(ii.regions, preserve_attributes=False)

            scores = list(ii.scores())
            insulation_scores.scores(window_size, scores)
            ii.close()

        return insulation_scores


class DirectionalityIndex(RegionScoreTable):

    _classid = 'DIRECTIONALITYINDEX'

    def __init__(*args, **kwargs):
        RegionScoreTable.__init__(*args, **kwargs)

    @classmethod
    def from_hic(cls, hic, window_size=2000000, weight_field=None,
                 file_name=None, tmpdir=None, **kwargs):
        bin_window_size = hic.distance_to_bins(window_size)
        directionality_index_regions = cls(file_name=file_name, mode='w', tmpdir=tmpdir)
        directionality_index_regions.add_regions(hic.regions, preserve_attributes=False)

        weight_field = hic._default_score_field if weight_field is None else weight_field

        n_bins = len(hic.regions)

        # find the distances of each region to the chromosome ends in bins
        boundary_dist = np.zeros(n_bins, dtype=int)
        last_chromosome = None
        last_chromosome_index = 0
        for i, region in enumerate(hic.regions(lazy=True)):
            chromosome = region.chromosome
            if last_chromosome is not None and chromosome != last_chromosome:
                chromosome_length = i - last_chromosome_index
                for j in range(chromosome_length):
                    boundary_dist[last_chromosome_index + j] = min(j, i - last_chromosome_index - 1 - j)
                last_chromosome_index = i
            last_chromosome = chromosome
        chromosome_length = n_bins - last_chromosome_index
        for j in range(chromosome_length):
            boundary_dist[last_chromosome_index + j] = min(j, n_bins - last_chromosome_index - 1 - j)

        left_sums = np.zeros(n_bins)
        right_sums = np.zeros(n_bins)
        directionality_index = np.zeros(n_bins)
        kwargs['inter_chromosomal'] = False
        kwargs['lazy'] = True
        for edge in hic.edges(**kwargs):
            source = edge.source
            sink = edge.sink
            weight = getattr(edge, weight_field)
            if source == sink:
                continue
            if sink - source <= bin_window_size:
                if boundary_dist[sink] >= sink-source:
                    left_sums[sink] += weight
                if boundary_dist[source] >= sink-source:
                    right_sums[source] += weight

        for i in range(n_bins):
            A = left_sums[i]
            B = right_sums[i]
            E = (A+B)/2
            if E != 0 and B-A != 0:
                directionality_index[i] = ((B-A)/abs(B-A)) * ((((A-E)**2)/E) + (((B-E)**2)/E))

        directionality_index_regions.region_data('score', directionality_index)

        return directionality_index_regions


class DirectionalityIndexes(RegionScoreParameterTable):

    _classid = 'DIRECTIONALITYINDEXES'

    def __init__(*args, **kwargs):
        RegionScoreParameterTable.__init__(*args, **kwargs)

    @property
    def window_sizes(self):
        return self._parameters

    @classmethod
    def from_hic(cls, hic, window_sizes, file_name=None, tmpdir=None,
                 **kwargs):
        directionality_indexes = cls(parameter_prefix='directionality_',
                                     parameter_values=list(window_sizes),
                                     file_name=file_name, mode='w',
                                     tmpdir=tmpdir)

        for i, window_size in enumerate(window_sizes):
            logger.debug("Window size {}".format(window_size))
            di = DirectionalityIndex.from_hic(hic, window_size, **kwargs)
            if i == 0:
                directionality_indexes.add_regions(di.regions, preserve_attributes=False)

            scores = list(di.scores())
            directionality_indexes.scores(window_size, scores)

        return directionality_indexes


class Boundaries(RegionScoreTable):

    _classid = 'BOUNDARIES'

    def __init__(*args, **kwargs):
        RegionScoreTable.__init__(*args, **kwargs)

    @classmethod
    def from_insulation_score(cls, insulation_score, min_score=None,
                              delta_window=3, log=False,
                              sub_bin_precision=False, call_maxima=False,
                              score_field='score',
                              **kwargs):
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
        index = list(insulation_score.region_data(score_field))
        if log:
            index = np.log2(index)
        peaks = MaximaCallerDelta(index, window_size=delta_window, sub_bin_precision=sub_bin_precision)
        if call_maxima:
            minima, scores = peaks.get_maxima()
        else:
            minima, scores = peaks.get_minima()
        regions = list(insulation_score.regions)

        boundaries = []
        for i, ix in enumerate(minima):
            if min_score is not None and scores[i] < min_score:
                continue
            if sub_bin_precision:
                region = regions[int(ix)]
                fraction = ix % 1
                shift = int((fraction - .5) * (region.end - region.start))
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

        boundary_regions = cls(**kwargs)
        boundary_regions.add_regions(boundaries)
        logger.info("Found {} boundaries".format(len(boundaries)))
        return boundary_regions
