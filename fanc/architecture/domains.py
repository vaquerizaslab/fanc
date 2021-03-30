from genomic_regions import GenomicRegion, RegionBased
from .maxima_callers import MaximaCallerDelta
from .helpers import RegionScoreMatrix
from ..regions import RegionsTable
from ..tools.matrix import apply_sliding_func, trim_stats
from ..tools.files import write_bed, write_bigwig, write_gff
from ..tools.matrix import nangmean
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
        """
        Return scores as list.

        :return: list of float
        """
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
        """
        Return scores for a specific parameter size as list.

        :param parameter: Parameter scores were calculated for, such as window size
        :param scores: If provided, set scores for this parameter to the ones in this list.
        :return: list of scores
        """
        score_field = self._score_field_converter([parameter])[0]
        return self.region_data(score_field, scores)

    def score_regions(self, parameter, **kwargs):
        """
        Construct a new object with regions that have a ``score`` attribute
        which corresponds to scores calculated with this parameter.
        :param parameter: Use scores calculated with this parameter (e.g. window size)
        :param kwargs: Keyword arguments passed to :class:`~fanc.regions.RegionsTable`
        :return: :class:`~fanc.architecture.domains.RegionScoreTable`
        """
        score_field = self._score_field_converter([parameter])[0]
        regions = []
        for region in self.regions:
            r = GenomicRegion(chromosome=region.chromosome,
                              start=region.start, end=region.end,
                              strand=region.strand, score=getattr(region, score_field))
            regions.append(r)
        scores = RegionScoreTable(**kwargs)
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
        """
        Write scores to BED file.

        :param file_name: Path to output file
        :param parameter: Parameter the scores were calculated for, such as window size
        :param subset: A :class:`~genomic_regions.GenomicRegion` or region string specifying a
                       region range to be written to file, e.g. "chr19:1-1mb"
        """
        return self._to_file(file_name, parameter, subset=subset, _write_function=write_bed)

    def to_gff(self, file_name, parameter, subset=None):
        """
        Write scores to GFF file.

        :param file_name: Path to output file
        :param parameter: Parameter the scores were calculated for, such as window size
        :param subset: A :class:`~genomic_regions.GenomicRegion` or region string specifying a
                       region range to be written to file, e.g. "chr19:1-1mb"
        """
        return self._to_file(file_name, parameter, subset=subset, _write_function=write_gff)

    def to_bigwig(self, file_name, parameter, subset=None):
        """
        Write scores to BigWig file.

        :param file_name: Path to output file
        :param parameter: Parameter the scores were calculated for, such as window size
        :param subset: A :class:`~genomic_regions.GenomicRegion` or region string specifying a
                       region range to be written to file, e.g. "chr19:1-1mb"
        """
        return self._to_file(file_name, parameter, subset=subset, _write_function=write_bigwig)


class InsulationScore(RegionScoreTable):

    _classid = 'INSULATIONSCORE'

    def __init__(*args, **kwargs):
        RegionScoreTable.__init__(*args, **kwargs)

    @classmethod
    def from_hic(cls, hic, window_size_in_bp,
                 file_name=None, tmpdir=None,
                 *args, **kwargs):

        insulation_scores = InsulationScores.from_hic(hic, window_size_in_bp, *args, **kwargs)
        insulation_score_regions = cls(file_name=file_name, mode='w', tmpdir=tmpdir)
        insulation_score_regions.add_regions(hic.regions, preserve_attributes=False)
        insulation_score_regions.region_data('score', insulation_scores.scores(window_size_in_bp))

        return insulation_score_regions

    def boundaries(self, *args, **kwargs):
        return Boundaries.from_insulation_score(self, *args, **kwargs)


class InsulationScores(RegionScoreParameterTable):

    _classid = 'INSULATIONSCORES'

    def __init__(*args, **kwargs):
        RegionScoreParameterTable.__init__(*args, **kwargs)

    @property
    def window_sizes(self):
        """
        Get a list of window sizes in this object.

        :return: list of window sizes (int)
        """
        return self._parameters

    @classmethod
    def from_hic(cls, hic, window_sizes, window_offset=0,
                 file_name=None, tmpdir=None, impute_missing=False,
                 na_threshold=0.5, normalise=True, normalisation_window=None,
                 trim_mean_proportion=0.0, geometric_mean=False,
                 subtract_mean=False, log=True):
        """
        Calculate insulation scores with multiple window sizes.

        Insulation scores provide a great way to quantify the level of
        interactions the cross each genomic region. It is calculated by summing up
        (normalised) contacts in a square next to the diagonal for each genomic
        region. Therefore, low scores correspond to highly insulated regions with
        few interactions spanning them.

        :param hic: A Hi-C object
        :param window_sizes: A window size or list of window sizes used for the
                             sliding window
        :param window_offset: An offset of the sliding window in bins from the
                              diagonal
        :param file_name: Path to file where insulation scores are saved
        :param tmpdir: Optional. If ``True``, will work with file in temporary
                       directory until it is closed
        :param impute_missing: Will replace missing / masked values in matrix
                               with their expected value prior to insulation
                               score calculation
        :param na_threshold: Fraction of missing values that is tolerated
                             in a sliding window before the score is set to NaN
        :param normalise: Normalise insulation score by dividing by chromosome mean
                          or mean of a sliding window if ``normalisation_window``
                          is set.
        :param normalisation_window: If ``None`` (default), normalisation is performed
                                     by dividing insulation scores by the chromosome mean.
                                     You can set this to a number of bins to perform a more
                                     local normalisation using average values in a window of
                                     that size
        :param trim_mean_proportion: If > 0 will use a trimmed mean for normalisation
                                     trimming this fraction of scores before calculating the
                                     mean. Use this if you expect outliers in insulation
                                     scores
        :param geometric_mean: Use a geometric mean instead of arithmetic. If using
                               log-transformed, and if you intend to subtract scores
                               from different samples for comparison, this is
                               recommended
        :param subtract_mean: For normalisation, subtract mean instead of dividing by it.
        :param log: Log2-transform insulation scores after calculation. In the default
                    parameters, this makes scores roughly symmetrical around 0.
        :return: :class:`~fanc.architecture.domains.InsulationScores`
        """
        if isinstance(window_sizes, int):
            window_sizes = [window_sizes]
        window_sizes = sorted(window_sizes)

        bin_window_sizes = [int(hic.distance_to_bins(w) / 2) for w in window_sizes]

        insulation_scores = cls(parameter_prefix='insulation_',
                                parameter_values=window_sizes,
                                file_name=file_name, mode='w',
                                tmpdir=tmpdir)
        insulation_scores.add_regions(hic.regions, preserve_attributes=False)

        if geometric_mean:
            avg_stat = nangmean
        else:
            avg_stat = np.nanmean

        chromosome_bins = hic.chromosome_bins
        window_offset += 1

        mappable = hic.mappable()

        if impute_missing:
            intra_expected, intra_expected_chromosome, _ = hic.expected_values()
        else:
            intra_expected, intra_expected_chromosome = None, None

        ii_list = [[] for _ in bin_window_sizes]
        chromosomes = hic.chromosomes()
        for chr_ix, chromosome in enumerate(chromosomes):
            logger.info("{} ({}/{})".format(chromosome, chr_ix + 1, len(chromosomes)))
            chromosome_start, chromosome_stop = chromosome_bins[chromosome]
            values_by_chromosome = [[0 for _ in range(chromosome_start, chromosome_stop)] for _ in bin_window_sizes]

            # add each edge weight to every insulation window that contains it
            for edge in hic.edges((chromosome, chromosome), lazy=True):
                i = edge.source - chromosome_start + window_offset
                j = edge.sink - chromosome_start - window_offset

                for w_ix, bin_window_size in enumerate(bin_window_sizes):
                    start = max(i, j - bin_window_size + 1)
                    stop = min(j + 1, i + bin_window_size)
                    for ii_bin in range(start, stop):
                        weight = getattr(edge, hic._default_score_field)
                        values_by_chromosome[w_ix][ii_bin] += weight

            # add imputed values, if requested
            if impute_missing:
                expected = intra_expected_chromosome[chromosome]

                covered = set()
                for ix in range(chromosome_start, chromosome_stop):
                    if not mappable[ix]:
                        for w_ix, bin_window_size in enumerate(bin_window_sizes):
                            sink = ix
                            for source in range(max(chromosome_start, ix - bin_window_size - window_offset - 2), ix):
                                if (source, sink) in covered:
                                    continue
                                covered.add((source, sink))
                                weight = expected[sink - source]
                                i = source - chromosome_start + window_offset
                                j = sink - chromosome_start - window_offset

                                start = max(i, j - bin_window_size + 1)
                                stop = min(j + 1, i + bin_window_size)
                                for ii_bin in range(start, stop):
                                    values_by_chromosome[w_ix][ii_bin] += weight

                            source = ix
                            for sink in range(ix, min(chromosome_stop, ix + window_offset + bin_window_size + 2)):
                                if (source, sink) in covered:
                                    continue
                                covered.add((source, sink))
                                weight = expected[sink - source]
                                i = source - chromosome_start + window_offset
                                j = sink - chromosome_start - window_offset

                                start = max(i, j - bin_window_size + 1)
                                stop = min(j + 1, i + bin_window_size)
                                for ii_bin in range(start, stop):
                                    values_by_chromosome[w_ix][ii_bin] += weight

                for w_ix, bin_window_size in enumerate(bin_window_sizes):
                    for k in range(len(values_by_chromosome[w_ix])):
                        if (k - window_offset < bin_window_size - 1
                                or k + window_offset > len(values_by_chromosome[w_ix]) - bin_window_size):
                            values_by_chromosome[w_ix][k] = np.nan
                        else:
                            values_by_chromosome[w_ix][k] /= bin_window_size ** 2
            # count unmappable bins in every window
            else:
                unmappable_horizontal = [[0 for _ in range(chromosome_start, chromosome_stop)] for _ in bin_window_sizes]
                unmappable_vertical = [[0 for _ in range(chromosome_start, chromosome_stop)] for _ in bin_window_sizes]

                for ix in range(chromosome_start, chromosome_stop):
                    if not mappable[ix]:
                        for w_ix, bin_window_size in enumerate(bin_window_sizes):
                            # horizontal
                            start_bin = ix - chromosome_start + window_offset
                            for ii_bin in range(start_bin, start_bin + bin_window_size):
                                if 0 <= ii_bin < len(unmappable_horizontal[w_ix]):
                                    unmappable_horizontal[w_ix][ii_bin] += 1
                            # vertical
                            start_bin = ix - chromosome_start - window_offset
                            for ii_bin in range(start_bin - bin_window_size + 1, start_bin + 1):
                                if 0 <= ii_bin < len(unmappable_vertical[w_ix]):
                                    unmappable_vertical[w_ix][ii_bin] += 1

                for w_ix, bin_window_size in enumerate(bin_window_sizes):
                    for k in range(len(values_by_chromosome[w_ix])):
                        na_vertical = unmappable_vertical[w_ix][k] * bin_window_size
                        na_horizontal = unmappable_horizontal[w_ix][k] * bin_window_size
                        na_overlap = unmappable_horizontal[w_ix][k] * unmappable_vertical[w_ix][k]
                        na_total = na_vertical + na_horizontal - na_overlap

                        # take into account nan values when adding zeros
                        if ((na_total > (bin_window_size ** 2 * na_threshold))
                                or k - window_offset < bin_window_size - 1
                                or k + window_offset > len(values_by_chromosome[w_ix]) - bin_window_size):
                            values_by_chromosome[w_ix][k] = np.nan
                        else:
                            values_by_chromosome[w_ix][k] /= (bin_window_size ** 2 - na_total)

            for w_ix, bin_window_size in enumerate(bin_window_sizes):
                ii_by_chromosome = values_by_chromosome[w_ix]

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)

                    ii_by_chromosome = np.array(ii_by_chromosome)
                    if normalise:
                        logger.debug("Normalising insulation score")
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
                    ii_list[w_ix].append(ii_by_chromosome)

        for w_ix, window_size in enumerate(window_sizes):
            ins_matrix = np.array(list(itertools.chain.from_iterable(ii_list[w_ix])))

            if log:
                logger.debug("Log-transforming insulation score")
                ins_matrix = np.log2(ins_matrix)

            insulation_scores.scores(window_size, ins_matrix)

        return insulation_scores


class DirectionalityIndex(RegionScoreTable):

    _classid = 'DIRECTIONALITYINDEX'

    def __init__(*args, **kwargs):
        RegionScoreTable.__init__(*args, **kwargs)

    @classmethod
    def from_hic(cls, hic, window_size=2000000,
                 file_name=None, tmpdir=None, *args, **kwargs):
        directionality_indexes = DirectionalityIndexes.from_hic(hic, window_size, *args, **kwargs)
        directionality_index_regions = cls(file_name=file_name, mode='w', tmpdir=tmpdir)
        directionality_index_regions.add_regions(hic.regions, preserve_attributes=False)
        directionality_index_regions.region_data('score', directionality_indexes.score(window_size))

        return directionality_index_regions


class DirectionalityIndexes(RegionScoreParameterTable):

    _classid = 'DIRECTIONALITYINDEXES'

    def __init__(*args, **kwargs):
        RegionScoreParameterTable.__init__(*args, **kwargs)

    @property
    def window_sizes(self):
        """
        Get a list of window sizes in this object.

        :return: list of window sizes (int)
        """
        return self._parameters

    @classmethod
    def from_hic(cls, hic, window_sizes, weight_field=None,
                 file_name=None, tmpdir=None, **kwargs):
        """
        Compute the directionality index for multiple window sizes.

        :param hic: A compatible Hi-C object
        :param window_sizes: A list of window sizes
        :param weight_field: Internal. Key of the weight attribute for an edge in this object.
        :param file_name: Path to output file. If not provided, will work in memory.
        :param tmpdir: Optional. If ``True``, will work in temporary directory until
                       file is closed.
        :param kwargs: Keyword arguments passed on to :class:`~fanc.matrix.RegionMatrixTable.edges`
        :return: :class:`~fanc.architecture.domains.DirectionalityIndexes`
        """
        kwargs['lazy'] = True

        if isinstance(window_sizes, int):
            window_sizes = [window_sizes]
        window_sizes = list(window_sizes)
        bin_window_sizes = [int(hic.distance_to_bins(w) / 2) for w in window_sizes]

        directionality_indexes = cls(parameter_prefix='directionality_',
                                     parameter_values=list(window_sizes),
                                     file_name=file_name, mode='w',
                                     tmpdir=tmpdir)
        directionality_indexes.add_regions(hic.regions, preserve_attributes=False)

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

        left_sums = [np.zeros(n_bins) for _ in window_sizes]
        right_sums = [np.zeros(n_bins) for _ in window_sizes]
        kwargs['inter_chromosomal'] = False
        kwargs['lazy'] = True
        for edge in hic.edges(**kwargs):
            source = edge.source
            sink = edge.sink
            weight = getattr(edge, weight_field)
            if source == sink:
                continue

            for w_ix, bin_window_size in enumerate(bin_window_sizes):
                if sink - source <= bin_window_size:
                    if boundary_dist[sink] >= sink - source:
                        left_sums[w_ix][sink] += weight
                    if boundary_dist[source] >= sink - source:
                        right_sums[w_ix][source] += weight

        mappability = hic.mappable()

        directionality_index = [np.zeros(n_bins) for _ in window_sizes]
        for w_ix, bin_window_size in enumerate(bin_window_sizes):
            corr, uncorr = 0, 0
            for i in range(n_bins):
                lm = np.sum(mappability[max(0, i - bin_window_size):i])
                rm = np.sum(mappability[i:i + bin_window_size])
                if lm / bin_window_size < 0.5 or rm / bin_window_size < 0.5:
                    directionality_index[w_ix][i] = np.nan
                    continue

                Au = left_sums[w_ix][i]
                Bu = right_sums[w_ix][i]
                # correct for mappability
                A = Au + Au / lm * (bin_window_size - lm)
                B = Bu + Bu / rm * (bin_window_size - rm)

                if Au != A:
                    corr += 1
                else:
                    uncorr += 1

                E = (A + B) / 2
                if E != 0 and B - A != 0:
                    directionality_index[w_ix][i] = ((B - A) / abs(B - A)) * ((((A - E) ** 2) / E) +
                                                                              (((B - E) ** 2) / E))

        for w_ix, window_size in enumerate(window_sizes):
            directionality_indexes.scores(window_size, directionality_index[w_ix])

        return directionality_indexes


class Boundaries(RegionScoreTable):

    _classid = 'BOUNDARIES'

    def __init__(*args, **kwargs):
        RegionScoreTable.__init__(*args, **kwargs)

    @classmethod
    def from_insulation_score(cls, insulation_score, window_size=None,
                              min_score=None,
                              delta_window=3, log=False,
                              sub_bin_precision=False, call_maxima=False,
                              score_field='score',
                              **kwargs):
        """
        Call insulation boundaries based on minima in an insulation vector of this object.

        :param insulation_score: :class:`~fanc.architecture.domains.InsulationScores` or
                                 :class:`~fanc.architecture.domains.InsulationScore` object
        :param window_size: Window size in base pairs. Only necessary for
                            :class:`~fanc.architecture.domains.InsulationScores` objects
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
        :param score_field:
        :return: list of :class:`~fanc.data.genomic.GenomicRegion`
        """
        if isinstance(insulation_score, InsulationScores):
            logger.debug("Found multiple insulation scores, extracting window size {}".format(window_size))
            if window_size is None:
                raise ValueError("Must provide window size!")
            index = list(insulation_score.region_data('insulation_{}'.format(window_size)))
            regions = list(insulation_score.regions)
        elif isinstance(insulation_score, RegionBased):
            logger.debug("Found RegionBased object, extracting regions and scores")
            index = [getattr(r, score_field) for r in insulation_score.regions]
            regions = list(insulation_score.regions)
        else:
            logger.debug("Assuming region iterator")
            index = [getattr(r, score_field) for r in insulation_score]
            regions = list(insulation_score)

        if log:
            index = np.log2(index)
        peaks = MaximaCallerDelta(index, window_size=delta_window, sub_bin_precision=sub_bin_precision)
        if call_maxima:
            minima, scores = peaks.get_maxima()
        else:
            minima, scores = peaks.get_minima()
        #regions = list(insulation_score.regions)

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
