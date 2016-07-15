from __future__ import division
from kaic.architecture.architecture import TableArchitecturalFeature, calculateondemand, ArchitecturalFeature
from kaic.architecture.genome_architecture import MatrixArchitecturalRegionFeature, VectorArchitecturalRegionFeature, \
    MatrixArchitecturalRegionFeatureFilter
from kaic.architecture.maxima_callers import MaximaCallerDelta
from kaic.data.genomic import GenomicRegion, HicEdgeFilter, Edge, Hic
from collections import defaultdict
from kaic.tools.general import ranges, to_slice
from kaic.tools.matrix import apply_sliding_func
from kaic.data.general import FileGroup
import numpy as np
import tables as t
import itertools
import logging
from bisect import bisect_left
from kaic.tools.general import RareUpdateProgressBar
logging.basicConfig(level=logging.INFO)


class HicArchitecture(object):
    def __init__(self, hic):
        self.hic = hic

    @property
    def expected_contacts(self, regions=None):
        return ExpectedContacts(self.hic, smooth=False, regions=regions)

    @property
    def possible_contacts(self, regions=None):
        return PossibleContacts(self.hic, regions=regions)

    @property
    def directionality_index(self, window_sizes=2000000):
        if isinstance(window_sizes, int):
            window_sizes = (window_sizes, )
        return DirectionalityIndex(self.hic, window_sizes=window_sizes)


class HicEdgeCollection(MatrixArchitecturalRegionFeature):
    _classid = 'HICEDGECOLLECTION'

    def __init__(self, hics=None, additional_fields=None, file_name=None, mode='a', tmpdir=None,
                 only_intra_chromosomal=False):
        if not isinstance(hics, str) and hics is not None:

            if additional_fields is not None:
                if not isinstance(additional_fields, dict) and issubclass(additional_fields, t.IsDescription):
                    # IsDescription subclass case
                    additional_fields = additional_fields.columns
            else:
                additional_fields = dict()

            original_fields = additional_fields.copy()

            self.shared_base_field_names = []
            for i, hic in enumerate(hics):
                field_descriptions = hic._field_dict
                for name, description in field_descriptions.iteritems():

                    if name.startswith("_") or name in {'source', 'sink'}:
                        continue

                    for j in xrange(len(hics)):
                        new_name = "%s_%d" % (name, j)
                        if new_name in original_fields:
                            raise ValueError("%s already in hic object, please choose different name." % new_name)
                        if new_name in additional_fields:
                            continue

                        if name not in self.shared_base_field_names:
                            self.shared_base_field_names.append(name)

                        description_args = {'dflt': description.dflt,
                                            'pos': len(additional_fields)}
                        if description.__class__ == t.description.StringCol:
                            description_args['itemsize'] = description.itemsize

                        additional_fields[new_name] = description.__class__(**description_args)
            MatrixArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode,
                                                      data_fields=additional_fields,
                                                      regions=hics[0].regions, tmpdir=tmpdir)
        elif (hics is None and file_name is not None) or file_name is None:
            if hics is not None and file_name is None:
                file_name = hics
            MatrixArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir)
            self._calculated = True

        self.hics = hics
        self.only_intra_chromosomal = only_intra_chromosomal

    def _calculate(self, *args, **kwargs):
        chromosomes = self.hics[0].chromosomes()

        # step 1: combine all hic objects into one
        #         and calculate variance
        for chr_i, chromosome1 in enumerate(chromosomes):
            for chr_j in xrange(chr_i, len(chromosomes)):
                chromosome2 = chromosomes[chr_j]

                if self.only_intra_chromosomal and chromosome1 != chromosome2:
                    continue

                logging.info("Processing chromosomes %s-%s" % (chromosome1, chromosome2))

                edges = dict()
                for i, hic in enumerate(self.hics):
                    logging.info("Processing Hic %d/%d" % (i, len(self.hics)))

                    for edge in hic.edge_subset(key=(chromosome1, chromosome2), lazy=True):
                        key = (edge.source, edge.sink)
                        if key not in edges:
                            edges[key] = {}
                        for field in self.shared_base_field_names:
                            if field not in edges[key]:
                                edges[key][field] = [None] * len(self.hics)
                            edges[key][field][i] = getattr(edge, field, None)

                for key, d in edges.iteritems():
                    for field, values in d.iteritems():
                        source = key[0]
                        sink = key[1]

                        d = dict()
                        for i, value in enumerate(values):
                            d["%s_%d" % (field, i)] = value

                        self.add_edge(Edge(source=source, sink=sink, **d), flush=False)

        self.flush()


class ExpectedContacts(TableArchitecturalFeature):
    _classid = 'EXPECTEDCONTACTS'

    def __init__(self, hic, file_name=None, mode='a', tmpdir=None, smooth=True, min_reads=400,
                 regions=None, weight_column=None, _table_name='expected_contacts'):
        if isinstance(hic, str):
            file_name = hic
            hic = None

        TableArchitecturalFeature.__init__(self, _table_name,
                                           {'distance': t.Int64Col(), 'intra': t.Float32Col(),
                                            'contacts': t.Float32Col(), 'pixels': t.Float32Col()},
                                           file_name=file_name, mode=mode, tmpdir=tmpdir)

        self.hic = hic
        self.smooth = smooth
        self.min_reads = min_reads
        self.regions = regions
        if weight_column is None:
            self.weight_column = self.hic.default_field
        else:
            self.weight_column = weight_column

    def _calculate(self):
        """
        Get intra- and inter-chromosomal expected contact counts.
        """

        # extract mappable regions from Hic object
        marginals = self.hic.marginals(weight_column=self.weight_column)
        regions = []
        all_regions = []
        region_ix = set()
        if self.regions is None:
            for region in self.hic.regions(lazy=False):
                all_regions.append(region)
                if marginals[region.ix] > 0:
                    regions.append(region)
                    region_ix.add(region.ix)
        else:
            if isinstance(self.regions, str) or isinstance(self.regions, GenomicRegion):
                self.regions = [self.regions]

            for region in self.regions:
                if isinstance(region, str):
                    region = GenomicRegion.from_string(region)

                for r in self.hic.subset(region, lazy=False):
                    all_regions.append(r)
                    if marginals[r.ix] > 0:
                        regions.append(r)
                        region_ix.add(r.ix)

        # count the number of regions per chromosome
        # and find the first and last region in each chromosome
        regions_by_chromosome = defaultdict(int)
        min_region_by_chromosome = dict()
        max_region_by_chromosome = dict()
        for region in regions:
            if (region.chromosome not in max_region_by_chromosome or
                        max_region_by_chromosome[region.chromosome] < region.ix):
                max_region_by_chromosome[region.chromosome] = region.ix

            if (region.chromosome not in min_region_by_chromosome or
                        min_region_by_chromosome[region.chromosome] > region.ix):
                min_region_by_chromosome[region.chromosome] = region.ix

            regions_by_chromosome[region.chromosome] += 1

        # find the largest distance between two regions
        # in the entire intra-chromosomal genome
        max_distance = 0
        for chromosome in min_region_by_chromosome:
            max_distance = max(max_distance,
                               max_region_by_chromosome[chromosome] - min_region_by_chromosome[chromosome])

        # find maximum theoretical distance, even if unmappable
        min_region_by_chromosome_t = dict()
        max_region_by_chromosome_t = dict()
        for region in all_regions:
            if (region.chromosome not in max_region_by_chromosome_t or
                        max_region_by_chromosome_t[region.chromosome] < region.ix):
                max_region_by_chromosome_t[region.chromosome] = region.ix

            if (region.chromosome not in min_region_by_chromosome_t or
                        min_region_by_chromosome_t[region.chromosome] > region.ix):
                min_region_by_chromosome_t[region.chromosome] = region.ix

        # find the largest distance between two regions
        # in the entire intra-chromosomal genome
        max_distance_t = 0
        for chromosome in min_region_by_chromosome_t:
            max_distance_t = max(max_distance_t,
                                 max_region_by_chromosome_t[chromosome] - min_region_by_chromosome_t[chromosome])

        # get the number of pixels at a given bin distance
        pixels_by_distance = [0.0] * (max_distance + 1)
        for chromosome, n in regions_by_chromosome.iteritems():
            for distance in xrange(0, n):
                pixels_by_distance[distance] += n-distance

        # build a reverse-lookup chromosome map to quickly
        # determine if an edge is intra-chromosomal
        chromosome_map = dict()
        for i, chromosome in enumerate(self.hic.chromosomes()):
            chromosome_map[chromosome] = i

        chromosomes = np.zeros(len(self.hic.regions), dtype=int)
        for i, region in enumerate(self.hic.regions(lazy=True)):
            chromosomes[i] = chromosome_map[region.chromosome]

        # get the number of reads at a given bin distance
        reads_by_distance = [0.0] * (max_distance + 1)
        inter_observed = 0
        for edge in self.hic.edges(lazy=True):
            source = edge.source
            sink = edge.sink
            # skip excluded regions
            if source not in region_ix or sink not in region_ix:
                continue
            # only intra-chromosomal distances
            if chromosomes[edge.source] == chromosomes[edge.sink]:
                reads_by_distance[edge.sink - edge.source] += getattr(edge, self.weight_column)
            else:
                inter_observed += getattr(edge, self.weight_column)

        with PossibleContacts(self.hic, regions=self.regions, weight_column=self.weight_column) as pc:
            intra_possible, inter_possible = pc.intra_possible(), pc.inter_possible()

        try:
            inter_expected = inter_observed/inter_possible
        except ZeroDivisionError:
            inter_expected = 0

        distance = np.arange(len(reads_by_distance))*self.hic.bin_size
        self.data('distance', distance)

        while len(reads_by_distance) < max_distance_t+1:
            reads_by_distance.append(reads_by_distance[-1])
            pixels_by_distance.append(pixels_by_distance[-1])

        # return here if smoothing not requested
        if not self.smooth:
            intra_expected = np.array(reads_by_distance)/np.array(pixels_by_distance)

            self.data('intra', intra_expected)
            self.data('contacts', reads_by_distance)
            self.data('pixels', pixels_by_distance)
            self._table.attrs['inter'] = inter_expected
            return

        # smoothing
        smoothed_reads_by_distance = np.zeros(len(reads_by_distance))
        smoothed_pixels_by_distance = np.zeros(len(pixels_by_distance))
        for i in xrange(len(reads_by_distance)):
            smoothed_reads = reads_by_distance[i]
            smoothed_pixels = pixels_by_distance[i]
            window_size = 0
            can_extend = True
            # smooth to a minimum number of reads per distance
            while smoothed_reads < self.min_reads and can_extend:
                window_size += 1
                can_extend = False
                # check if we can increase the window to the left
                if i - window_size >= 0:
                    smoothed_reads += reads_by_distance[i-window_size]
                    smoothed_pixels += pixels_by_distance[i-window_size]
                    can_extend = True
                # check if we can increase the window to the right
                if i + window_size < len(reads_by_distance):
                    smoothed_reads += reads_by_distance[i+window_size]
                    smoothed_pixels += pixels_by_distance[i+window_size]
                    can_extend = True
            smoothed_reads_by_distance[i] = smoothed_reads
            smoothed_pixels_by_distance[i] = smoothed_pixels

        intra_expected = smoothed_reads_by_distance/smoothed_pixels_by_distance

        self.data('intra', intra_expected)
        self.data('contacts', smoothed_reads_by_distance)
        self.data('pixels', smoothed_pixels_by_distance)
        self._table.attrs['inter'] = inter_expected

    @calculateondemand
    def intra_expected(self):
        return self[:, 'intra']

    @calculateondemand
    def inter_expected(self):
        return self._table.attrs['inter']

    @calculateondemand
    def distance(self):
        return self[:, 'distance']

    @calculateondemand
    def intra_contacts(self):
        return self[:, 'contacts']

    @calculateondemand
    def intra_pixels(self):
        return self[:, 'pixels']


class ObservedExpectedRatio(MatrixArchitecturalRegionFeature):
    _classid = 'OBSERVEDEXPECTEDRATIO'

    def __init__(self, hic, file_name=None, mode='a', tmpdir=None, regions=None,
                 weight_column='weight', _table_name='expected_contacts'):
        self.region_selection = regions

        # are we retrieving an existing object?
        if isinstance(hic, str) and file_name is None:
            file_name = hic
            hic = None
            MatrixArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir)
        else:
            if regions is None:
                regions = hic.regions
                self.region_conversion = {region.ix: region.ix for region in hic.regions}
            else:
                self.region_conversion = {region.ix: i for i, region in enumerate(hic.subset(regions))}
                regions = hic.subset(regions)
            MatrixArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                      data_fields={'ratio': t.Float32Col()}, regions=regions,
                                                      default_field='ratio', _table_name_edges=_table_name)
        self.hic = hic
        self.weight_column = weight_column

    def _calculate(self):
        with ExpectedContacts(self.hic, regions=self.region_selection, weight_column=self.weight_column) as ex:
            inter_expected = ex.inter_expected()
            intra_expected = ex.intra_expected()

            regions_dict = self.hic.regions_dict
            region_selection = self.region_selection if self.region_selection is not None else slice(0, None, None)
            for edge in self.hic.edge_subset(key=(region_selection, region_selection), lazy=True):
                source = edge.source
                new_source = self.region_conversion[source]
                sink = edge.sink
                new_sink = self.region_conversion[sink]
                weight = getattr(edge, self.weight_column)
                if regions_dict[source].chromosome == regions_dict[sink].chromosome:
                    expected = intra_expected[new_sink-new_source]
                else:
                    expected = inter_expected
                self.add_edge(Edge(new_source, new_sink, ratio=weight/expected), flush=False)
            self.flush()


class FoldChangeMatrix(MatrixArchitecturalRegionFeature):
    _classid = 'FOLDCHANGEMATRIX'

    def __init__(self, matrix1, matrix2=None, file_name=None, mode='a', tmpdir=None,
                 regions=None, scale_matrices=False, log2=True,
                 weight_column='weight', _table_name='expected_contacts'):
        self.region_selection = regions

        # are we retrieving an existing object?
        if isinstance(matrix1, str) and matrix2 is None and file_name is None:
            file_name = matrix1
            matrix1 = None
            MatrixArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                      default_field='fc', _table_name_edges=_table_name)
        else:
            if regions is None:
                regions = matrix1.regions
                self.region_conversion = {region.ix: region.ix for region in matrix1.regions}
            else:
                self.region_conversion = {region.ix: i for i, region in enumerate(matrix1.subset(regions))}
                regions = matrix1.subset(regions)
            MatrixArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                      data_fields={'fc': t.Float32Col()}, regions=regions,
                                                      default_field='fc', _table_name_edges=_table_name)
        self.matrix1 = matrix1
        self.matrix2 = matrix2
        self.weight_column = weight_column
        self.scale_matrices = scale_matrices
        self.log2 = log2

    def _calculate(self):
        if self.scale_matrices:
            scaling_factor = self.matrix1.scaling_factor(self.matrix2)
        else:
            scaling_factor = 1.

        chromosomes = self.chromosomes()

        for i in xrange(len(chromosomes)):
            chromosome1 = chromosomes[i]
            for j in xrange(i, len(chromosomes)):
                chromosome2 = chromosomes[j]

                edges1 = dict()
                for edge in self.matrix1.edge_subset(key=(chromosome1, chromosome2), lazy=True):
                    try:
                        source = self.region_conversion[edge.source]
                        sink = self.region_conversion[edge.sink]
                    except KeyError:
                        continue

                    edges1[(source, sink)] = getattr(edge, self.weight_column)

                for edge in self.matrix2.edge_subset(key=(chromosome1, chromosome2), lazy=True):
                    try:
                        source = self.region_conversion[edge.source]
                        sink = self.region_conversion[edge.sink]
                    except KeyError:
                        continue

                    if (source, sink) in edges1:
                        weight = edges1[(source, sink)] / (scaling_factor*edge.weight)
                        if self.log2:
                            weight = np.log2(weight)
                        self.add_edge([source, sink, weight], flush=False)

        self.flush()


class ABDomainMatrix(MatrixArchitecturalRegionFeature):
    _classid = 'ABDOMAINMATRIX'

    def __init__(self, hic, file_name=None, mode='a', tmpdir=None, regions=None,
                 ratio=True, weight_column='weight', per_chromosome=True, _table_name='ab_domains'):
        self.region_selection = regions

        # are we retrieving an existing object?
        if isinstance(hic, str) and file_name is None:
            file_name = hic
            hic = None
            MatrixArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir)
        else:
            if regions is None:
                regions = hic.regions
                self.region_conversion = {region.ix: region.ix for region in hic.regions}
            else:
                self.region_conversion = {region.ix: i for i, region in enumerate(hic.subset(regions))}
                regions = hic.subset(regions)
            MatrixArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                      data_fields={'correlation': t.Float32Col()}, regions=regions,
                                                      default_field='correlation', _table_name_edges=_table_name)
        self.hic = hic
        self.weight_column = weight_column
        self.ratio = ratio
        self.per_chromosome = per_chromosome

    def _calculate(self):
        oer = self.hic
        if self.ratio:
            oer = ObservedExpectedRatio(self.hic, weight_column=self.weight_column)

        if self.per_chromosome:
            chromosomes = self.chromosomes()
            for chromosome in chromosomes:
                    m = oer[chromosome, chromosome]
                    corr_m = np.corrcoef(m)
                    logging.info("Chromosome {}".format(chromosome))
                    with RareUpdateProgressBar(max_value=m.shape[0]) as pb:
                        for i, row_region in enumerate(m.row_regions):
                            for j, col_region in enumerate(m.col_regions):
                                if j < i:
                                    continue

                                source = self.region_conversion[row_region.ix]
                                sink = self.region_conversion[col_region.ix]
                                self.add_edge([source, sink, corr_m[i, j]], flush=False)
                            pb.update(i)
        else:
            m = oer[:]
            corr_m = np.corrcoef(m)
            with RareUpdateProgressBar(max_value=m.shape[0]) as pb:
                for i, row_region in enumerate(m.row_regions):
                    for j in xrange(i, len(m.row_regions)):
                        col_region = m.row_regions[j]
                        source = self.region_conversion[row_region.ix]
                        sink = self.region_conversion[col_region.ix]
                        self.add_edge([source, sink, corr_m[i, j]], flush=False)
                    pb.update(i)
        self.flush()

        if self.ratio:
            oer.close()


class ABDomains(VectorArchitecturalRegionFeature):
    _classid = 'ABDOMAINS'

    def __init__(self, data, file_name=None, mode='a', tmpdir=None,
                 per_chromosome=True, regions=None, _table_name='abdomains'):
        self.region_selection = regions

        # are we retrieving an existing object?
        if isinstance(data, str) and file_name is None:
            file_name = data
            data = None
            VectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir)
        else:
            if regions is None:
                regions = data.regions
            else:
                regions = data.subset(regions)

            fields = {'ev': t.Float32Col()}

            VectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                      data_fields=fields, regions=regions,
                                                      _table_name_data=_table_name)

        self.per_chromosome = per_chromosome
        self.data = data

    def _calculate(self):
        if isinstance(self.data, Hic):
            ab_data = ABDomainMatrix(self.data, regions=self.region_selection, per_chromosome=self.per_chromosome)
        else:
            ab_data = self.data

        ab_results = dict()
        if self.per_chromosome:
            for chromosome in self.chromosomes():
                m = ab_data[chromosome, chromosome]
                m[np.isnan(m)] = 0
                w, v = np.linalg.eig(m)
                ab_vector = v[:, 1]
                for i, region in enumerate(m.row_regions):
                    ab_results[region.ix] = ab_vector[i]
        else:
            m = ab_data[:]
            m[np.isnan(m)] = 0
            w, v = np.linalg.eig(m)
            ab_vector = v[:, 1]
            for i, region in enumerate(m.row_regions):
                ab_results[region.ix] = ab_vector[i]

        for region in self.regions(lazy=True):
            region.ev = ab_results[region.ix]
        self.flush()

    @calculateondemand
    def ab_domain_eigenvector(self):
        return self[:, 'ev']

    @calculateondemand
    def ab_regions(self):
        domains = []
        current_domain = None
        last_region = None
        for region in self.regions(lazy=False):
            domain_type = 'A' if region.ev >= 0 else 'B'

            if last_region is not None and region.chromosome != last_region.chromosome:
                current_domain = None

            if current_domain is None:
                current_domain = GenomicRegion(chromosome=region.chromosome, start=region.start, end=region.end,
                                               type=domain_type)
            else:
                if (region.ev < 0 and last_region.ev < 0) or (region.ev >= 0 and last_region.ev >= 0):
                    current_domain.end = region.end
                else:
                    domains.append(current_domain)
                    current_domain = GenomicRegion(chromosome=region.chromosome, start=region.start, end=region.end,
                                                   type=domain_type)
            last_region = region
        return domains


class PossibleContacts(TableArchitecturalFeature):
    _classid = 'POSSIBLECONTACTS'

    def __init__(self, hic, file_name=None, mode='a', tmpdir=None, regions=None,
                 weight_column='weight', _table_name='expected_contacts'):
        if isinstance(hic, str):
            file_name = hic
            hic = None

        TableArchitecturalFeature.__init__(self, _table_name,
                                           {'intra': t.Int32Col(), 'inter': t.Int32Col()},
                                           file_name=file_name, mode=mode, tmpdir=tmpdir)

        self.hic = hic
        self.regions = regions
        self.weight_column = weight_column

    def _calculate(self):
        marginals = self.hic.marginals(weight_column=self.weight_column)

        mappable = defaultdict(int)
        if self.regions is None:
            for r in self.hic.regions(lazy=True):
                if marginals[r.ix] > 0:
                    mappable[r.chromosome] += 1
        else:
            if isinstance(self.regions, str) or isinstance(self.regions, GenomicRegion):
                self.regions = [self.regions]

            for region in self.regions:
                if isinstance(region, str):
                    region = GenomicRegion.from_string(region)

                for r in self.hic.subset(region, lazy=True):
                    if marginals[r.ix] > 0:
                        mappable[r.chromosome] += 1

        # calculate possible combinations
        intra_possible = 0
        inter_possible = 0
        chromosomes = mappable.keys()
        for i in xrange(len(chromosomes)):
            chromosome1 = chromosomes[i]
            n1 = mappable[chromosome1]
            intra_possible += n1**2/2 + n1/2
            for j in xrange(i+1, len(chromosomes)):
                chromosome2 = chromosomes[j]
                n2 = mappable[chromosome2]
                inter_possible += n1*n2

        self.data('intra', [intra_possible])
        self.data('inter', [inter_possible])

    @calculateondemand
    def intra_possible(self):
        return self[0, 'intra']

    @calculateondemand
    def inter_possible(self):
        return self[0, 'inter']


class RowRegionMatrix(np.ndarray):
    def __new__(cls, input_matrix, regions=None, fields=None):
        obj = np.asarray(input_matrix).view(cls)
        obj.regions = regions
        obj.fields = fields
        obj._chromosome_index = None
        obj._region_index = None
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return

        self.regions = getattr(obj, 'regions', None)
        self.fields = getattr(obj, 'fields', None)
        self._chromosome_index = None
        self._region_index = None

    def _build_region_index(self):
        self._chromosome_index = defaultdict(list)
        self._region_index = defaultdict(list)
        for i, region in enumerate(self.regions):
            self._region_index[region.chromosome].append(i)
            self._chromosome_index[region.chromosome].append(region.end)

    def region_bins(self, region):
        if self._region_index is None or self._chromosome_index is None:
            self._build_region_index()
        if isinstance(region, str):
            region = GenomicRegion(region)
        start_ix = bisect_left(self._chromosome_index[region.chromosome], region.start)
        end_ix = bisect_left(self._chromosome_index[region.chromosome], region.end)
        start_region_ix = self._region_index[region.chromosome][start_ix]
        end_region_ix = self._region_index[region.chromosome][end_ix]
        return slice(start_region_ix, end_region_ix + 1)

    def __getitem__(self, index):
        self._getitem = True

        # convert string types into region indexes
        if isinstance(index, tuple):
            row_key = self._convert_region_key(index[0])
            col_key = self._convert_field_key(index[1])
            index = (row_key, col_key)
        else:
            row_key = self._convert_region_key(index)
            try:
                col_key = slice(0, len(self.fields), 1)
            except TypeError:
                col_key = None
            index = row_key

        try:
            out = np.ndarray.__getitem__(self, index)
        finally:
            self._getitem = False

        if not isinstance(out, np.ndarray):
            return out

        # get regions
        try:
            row_regions = self.regions[row_key]
        except TypeError:
            row_regions = None

        if not isinstance(col_key, list):
            try:
                col_fields = self.fields[col_key]
            except TypeError:
                col_fields = None
        else:
            try:
                col_fields = [self.fields[key] for key in col_key]
            except TypeError:
                col_fields = None

        out.fields = col_fields
        out.regions = row_regions

        return out

    def __getslice__(self, start, stop):
        return self.__getitem__(slice(start, stop))

    def _convert_region_key(self, key):
        if isinstance(key, str):
            key = GenomicRegion.from_string(key)
        if isinstance(key, GenomicRegion):
            return self.region_bins(key)
        return key

    def _convert_field_key(self, key):
        if isinstance(key, str):
            return self.fields.index(key)

        if isinstance(key, list):
            if len(key) == 0:
                raise ValueError("Length of supplied list is 0.")
            l = []
            for k in key:
                if isinstance(k, str):
                    k = self.fields.index(k)
                l.append(k)
            return l

        return key


class MultiVectorArchitecturalRegionFeature(VectorArchitecturalRegionFeature):
    _classid = 'MULTIVECTORARCHITECTURALREGIONFEATURE'

    def __init__(self, file_name=None, mode='a', data_fields=None,
                 regions=None, data=None, _table_name_data='region_data',
                 tmpdir=None):
        VectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, data_fields=data_fields,
                                                  regions=regions, data=data, _table_name_data=_table_name_data,
                                                  tmpdir=tmpdir)
        self._y_values = None

    def as_matrix(self, regions=None, keys=None):
        if regions is None:
            region_iter = self.regions(lazy=False)
        else:
            region_iter = self.subset(regions, lazy=False)

        if keys is None:
            keys = self.data_field_names
        elif isinstance(keys, slice) or isinstance(keys, int):
            keys = self.data_field_names[keys]

        array = []
        array_regions = []
        array_keys = [key for key in keys]
        for region in region_iter:
            array_regions.append(region)
            row = []
            for key in array_keys:
                row.append(getattr(region, key))
            array.append(row)

        return RowRegionMatrix(np.array(array), regions=array_regions, fields=array_keys)

    @property
    def y_values(self):
        return self._y_values

    @y_values.setter
    def y_values(self, values):
        if len(values) != len(self.data_field_names):
            raise ValueError("Length of y-values "
                             "({}) must be the same as length of data fields ({})".format(len(values),
                                                                                          len(self.data_field_names)))
        self._y_values = values


def load_array(file_name, mode='r', tmpdir=None):
    try:
        array = InsulationIndex(file_name, mode=mode, tmpdir=tmpdir)
        return array
    except t.FileModeError:
        pass

    try:
        array = DirectionalityIndex(file_name, mode=mode, tmpdir=tmpdir)
        return array
    except t.FileModeError:
        pass

    try:
        array = RegionContactAverage(file_name, mode=mode, tmpdir=tmpdir)
        return array
    except t.FileModeError:
        pass

    raise ValueError("Cannot recognise {} as array".format(file_name))


class DirectionalityIndex(MultiVectorArchitecturalRegionFeature):
    _classid = 'DIRECTIONALITYINDEX'

    def __init__(self, hic, file_name=None, mode='a', tmpdir=None,
                 weight_column=None, regions=None, window_sizes=(2000000,),
                 _table_name='directionality_index'):

        self.region_selection = regions

        # are we retrieving an existing object?
        if isinstance(hic, str) and file_name is None:
            file_name = hic
            MultiVectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                           _table_name_data=_table_name)
        else:
            di_fields = {}
            self.window_sizes = []
            for i, window_size in enumerate(window_sizes):
                di_fields['di_%d' % window_size] = t.Float32Col(pos=i)
                self.window_sizes.append(window_size)

            if regions is None:
                regions = hic.regions
            else:
                regions = hic.subset(regions)

            MultiVectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                           data_fields=di_fields, regions=regions,
                                                           _table_name_data=_table_name)

            self.hic = hic
            if weight_column is None:
                self.weight_column = self.hic.default_field
            else:
                self.weight_column = weight_column

        self.window_sizes = []
        for colname in self._regions.colnames:
            if colname.startswith("di_"):
                window_size = int(colname[3:])
                self.window_sizes.append(window_size)
        self.y_values = self.window_sizes

    def _get_boundary_distances(self):
        n_bins = len(self.hic.regions)
        # find distances to chromosome boundaries in bins
        boundary_dist = np.zeros(n_bins, dtype=int)
        last_chromosome = None
        last_chromosome_index = 0
        for i, region in enumerate(self.hic.regions(lazy=True)):
            chromosome = region.chromosome
            if last_chromosome is not None and chromosome != last_chromosome:
                chromosome_length = i-last_chromosome_index
                for j in xrange(chromosome_length):
                    boundary_dist[last_chromosome_index+j] = min(j, i-last_chromosome_index-1-j)
                last_chromosome_index = i
            last_chromosome = chromosome
        chromosome_length = n_bins-last_chromosome_index
        for j in xrange(chromosome_length):
            boundary_dist[last_chromosome_index+j] = min(j, n_bins-last_chromosome_index-1-j)

        return boundary_dist

    def _directionality_index(self, window_size=2000000):
        bin_window_size = self.hic.distance_to_bins(window_size)

        n_bins = len(self.hic.regions)
        boundary_dist = self._get_boundary_distances()

        if self.region_selection is not None:
            edge_iter = self.hic.edge_subset((self.region_selection, self.region_selection),
                                             only_intrachromosomal=True)
        else:
            edge_iter = self.hic.edges(lazy=True, only_intrachromosomal=True)

        left_sums = np.zeros(n_bins)
        right_sums = np.zeros(n_bins)
        directionality_index = np.zeros(n_bins)
        for edge in edge_iter:
            source = edge.source
            sink = edge.sink
            weight = getattr(edge, self.weight_column)
            if source == sink:
                continue
            if sink - source <= bin_window_size:
                if boundary_dist[sink] >= sink-source:
                    left_sums[sink] += weight
                if boundary_dist[source] >= sink-source:
                    right_sums[source] += weight

        for i in xrange(n_bins):
            A = left_sums[i]
            B = right_sums[i]
            E = (A+B)/2
            if E != 0 and B-A != 0:
                directionality_index[i] = ((B-A)/abs(B-A)) * ((((A-E)**2)/E) + (((B-E)**2)/E))

        if self.region_selection is not None:
            nodes_ix = self.hic._getitem_nodes(key=self.region_selection, as_index=True)

            if not isinstance(nodes_ix, list):
                nodes_ix = [nodes_ix]

            di_sub = []
            for node_range in ranges(nodes_ix):
                di_sub += list(directionality_index[node_range[0]:node_range[1]+1])

            return np.array(di_sub)
        return directionality_index

    def _calculate(self):
        for window_size in self.window_sizes:
            logging.info("Calculating directionality index for window size {}".format(window_size))
            directionality_index = self._directionality_index(window_size)
            self.data("di_%d" % window_size, directionality_index)

    @calculateondemand
    def directionality_index(self, window_size=None):
        if window_size is None:
            window_size = self.window_sizes[0]
        return self[:, 'di_%d' % window_size]


class InsulationIndex(MultiVectorArchitecturalRegionFeature):
    _classid = 'INSULATIONINDEX'

    def __init__(self, hic, file_name=None, mode='a', tmpdir=None,
                 regions=None, relative=False, offset=0, normalise=False, impute_missing=True,
                 window_sizes=(200000,), log=False, _normalisation_window=300,
                 subtract_mean=False, _table_name='insulation_index'):
        self.region_selection = regions

        # are we retrieving an existing object?
        if isinstance(hic, str) and file_name is None:
            file_name = hic
            hic = None
            MultiVectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                           _table_name_data=_table_name)
        else:
            if regions is None:
                regions = hic.regions
            else:
                regions = hic.subset(regions)

            ii_fields = {}
            self.window_sizes = []
            for i, window_size in enumerate(window_sizes):
                ii_fields['ii_%d' % window_size] = t.Float32Col(pos=i)
                self.window_sizes.append(window_size)

            MultiVectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                           data_fields=ii_fields, regions=regions,
                                                           _table_name_data=_table_name)

        self.window_sizes = []
        for colname in self._regions.colnames:
            if colname.startswith("ii_"):
                window_size = int(colname[3:])
                self.window_sizes.append(window_size)
        self.y_values = self.window_sizes

        self.offset = offset
        self.hic = hic
        self.relative = relative
        self.impute_missing = impute_missing
        self.normalise = normalise
        self.normalisation_window = _normalisation_window
        self.subtract_mean = subtract_mean
        self.log = log

    def _insulation_index(self, d1, d2, mask_thresh=.5, aggr_func=np.ma.mean, _mappable=None, _expected=None):
        if self.region_selection is not None:
            regions = self.hic.subset(self.region_selection)
        else:
            regions = self.hic.regions

        if _mappable is None:
            _mappable = self.hic.mappable()

        if self.impute_missing and _expected is None:
            _expected = ExpectedContacts(self.hic, smooth=True)

        nan_chromosome_counter = 0
        nan_mask_counter = 0
        nan_invalid_counter = 0
        logging.debug("Starting processing")
        skipped = 0
        last_chromosome = None
        ins_by_chromosome = []
        for i, r in enumerate(regions):
            if r.chromosome != last_chromosome:
                logging.info("Processing chromosome {}".format(r.chromosome))
                last_chromosome = r.chromosome
                ins_by_chromosome.append(list())
                unmasked = self.hic.as_matrix(key=(r.chromosome, r.chromosome),
                                              mask_missing=False, impute_missing=False)
                mask = np.zeros(unmasked.shape, dtype=bool)
                for z, ix in enumerate(xrange(r.ix, r.ix + unmasked.shape[0])):
                    if not _mappable[ix]:
                        mask[z, :] = True
                        mask[:, z] = True
                hic_matrix = np.ma.MaskedArray(unmasked, mask=mask)
                if self.impute_missing:
                    hic_matrix = self.hic._impute_missing_contacts(hic_matrix=hic_matrix,
                                                                   _expected_contacts=_expected)
                logging.info("Matrix loaded.")

            rix = len(ins_by_chromosome[-1])

            if rix < d2 or hic_matrix.shape[0] - rix <= d2 + 1:
                ins_by_chromosome[-1].append(np.nan)
                nan_chromosome_counter += 1
                continue

            if not _mappable[i]:
                ins_by_chromosome[-1].append(np.nan)
                nan_mask_counter += 1
                continue

            up_rel_slice = (slice(rix - d2, rix - d1), slice(rix - d2, rix - d1))
            down_rel_slice = (slice(rix + d1 + 1, rix + d2 + 1), slice(rix + d1 + 1, rix + d2 + 1))
            ins_slice = (slice(rix + d1 + 1, rix + d2 + 1), slice(rix - d2, rix - d1))

            if ((self.relative and np.sum(hic_matrix.mask[up_rel_slice]) > ((d2-d1)**2)*mask_thresh) or
                    (self.relative and np.sum(hic_matrix.mask[down_rel_slice]) > ((d2-d1)**2)*mask_thresh) or
                    np.sum(hic_matrix.mask[ins_slice]) > ((d2-d1)**2)*mask_thresh):
                # If too close to the edge of chromosome or
                # if more than half of the entries in this quadrant are masked (unmappable)
                # exclude it from the analysis
                skipped += 1
                ins_by_chromosome[-1].append(np.nan)
                nan_invalid_counter += 1
                continue

            if not self.relative:
                s = aggr_func(hic_matrix[ins_slice].data
                              if self.impute_missing else hic_matrix[ins_slice])
                ins_by_chromosome[-1].append(s)
            else:
                if not self.impute_missing:
                    s = (aggr_func(hic_matrix[ins_slice]) /
                         aggr_func(np.ma.dstack((hic_matrix[up_rel_slice],
                                                 hic_matrix[down_rel_slice]))))
                    ins_by_chromosome[-1].append(s)
                else:
                    s = (aggr_func(hic_matrix[ins_slice].data) /
                         aggr_func(np.ma.dstack((hic_matrix[up_rel_slice].data,
                                                 hic_matrix[down_rel_slice].data))))
                    ins_by_chromosome[-1].append(s)

        logging.info("Skipped {} regions because >{:.1%} of matrix positions were masked".format(skipped, mask_thresh))

        for i in xrange(len(ins_by_chromosome)):
            ins_by_chromosome[i] = np.array(ins_by_chromosome[i])
            if self.normalise:
                if self.normalisation_window is not None:
                    mean_ins = apply_sliding_func(ins_by_chromosome[i], self.normalisation_window, func=np.nanmean)
                else:
                    mean_ins = np.nanmean(ins_by_chromosome[i])
                if not self.subtract_mean:
                    ins_by_chromosome[i] = ins_by_chromosome[i] / mean_ins
                else:
                    ins_by_chromosome[i] = ins_by_chromosome[i] - mean_ins

        ins_matrix = np.array(list(itertools.chain.from_iterable(ins_by_chromosome)))

        logging.info("__nans__\nchromosome boundary: {}\nmasked region: {}\ninvalid balance: {}\n total: {}/{}".format(
            nan_chromosome_counter, nan_mask_counter, nan_invalid_counter, np.sum(np.isnan(ins_matrix)), len(ins_matrix)
        ))
        if self.log:
            return np.log2(ins_matrix)
        else:
            return ins_matrix

    def _calculate(self):
        mappable = self.hic.mappable()
        if self.impute_missing:
            ex = ExpectedContacts(self.hic, smooth=True)
        else:
            ex = None
        offset_bins = self.hic.distance_to_bins(self.offset)
        for window_size in self.window_sizes:
            logging.info("Calculating insulation index for window size {}".format(window_size))
            bins = self.hic.distance_to_bins(window_size)
            insulation_index = self._insulation_index(offset_bins, offset_bins+bins,
                                                      _mappable=mappable, _expected=ex)
            self.data("ii_%d" % window_size, insulation_index)
        if ex is not None:
            ex.close()

    @calculateondemand
    def insulation_index(self, window_size):
        if window_size is None:
            window_size = self.window_sizes[0]
        return self[:, 'ii_%d' % window_size]

    def boundaries(self, window_size, min_score=None, delta_window=7, log=False):
        index = self.insulation_index(window_size)
        if log:
            index = np.log2(index)
        peaks = MaximaCallerDelta(index, window_size=delta_window)
        minima, scores = peaks.get_minima()
        regions = list(self.regions)

        boundaries = []
        for i, ix in enumerate(minima):
            if min_score is not None and scores[i] < min_score:
                continue
            region = regions[ix]
            b = GenomicRegion(chromosome=region.chromosome, start=region.start, end=region.end, score=scores[i])
            boundaries.append(b)

        logging.info("Found {} boundaries for window size {}".format(len(boundaries), window_size))
        return boundaries


class RegionContactAverage(MultiVectorArchitecturalRegionFeature):
    _classid = 'REGIONCONTACTAVERAGE'

    def __init__(self, matrix, file_name=None, mode='a', tmpdir=None,
                 window_sizes=(200000,), regions=None, offset=0, padding=1, impute_missing=True,
                 _table_name='contact_average'):
        self.region_selection = regions

        # are we retrieving an existing object?
        if isinstance(matrix, str) and file_name is None:
            file_name = matrix
            matrix = None
            MultiVectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                           _table_name_data=_table_name)
        else:
            if regions is None:
                regions = matrix.regions
                self.region_conversion = {region.ix: region.ix for region in matrix.regions}
            else:
                self.region_conversion = {region.ix: i for i, region in enumerate(matrix.subset(regions))}
                regions = matrix.subset(regions)

            av_fields = {}
            self.window_sizes = []
            n = 0
            for i, window_size in enumerate(window_sizes):
                av_fields['av_%d' % window_size] = t.Float32Col(pos=n)
                av_fields['avl_%d' % window_size] = t.Float32Col(pos=n+1)
                av_fields['avr_%d' % window_size] = t.Float32Col(pos=n+2)
                n += 3
                self.window_sizes.append(window_size)

            MultiVectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                           data_fields=av_fields, regions=regions,
                                                           _table_name_data=_table_name)

        self.window_sizes = []
        for colname in self._regions.colnames:
            if colname.startswith("av_"):
                window_size = int(colname[3:])
                self.window_sizes.append(window_size)

        self.offset = offset
        self.padding = padding
        self.matrix = matrix
        self.impute_missing = impute_missing

    def _contact_average(self, window_size, offset=0, padding=1, _aggr_func=np.ma.mean):
        av_values = dict()
        for chromosome in self.chromosomes():
            matrix = self.matrix.as_matrix(key=(chromosome, chromosome), mask_missing=True,
                                           impute_missing=self.impute_missing)

            # region index
            for i, region in enumerate(matrix.row_regions):
                ix = self.region_conversion[region.ix]
                slice_left = slice(max(0, i-offset-window_size), max(0, i-offset+1))
                slice_right = slice(min(i+offset, matrix.shape[0]), min(i+offset+window_size+1, matrix.shape[0]))
                slice_vertical = slice(max(0, i-padding), min(i+padding, matrix.shape[0]))

                value_left = _aggr_func(matrix[slice_vertical, slice_left])
                value_right = _aggr_func(matrix[slice_vertical, slice_right])

                av_values[ix] = (value_left, value_right)

        av_left = np.zeros(len(self.regions))
        av_right = np.zeros(len(self.regions))
        for region in self.regions(lazy=True):
            if region.ix in av_values:
                av_left[region.ix], av_right[region.ix] = av_values[region.ix]
        return av_left, av_right

    def _calculate(self):
        offset_bins = self.matrix.distance_to_bins(self.offset)
        for window_size in self.window_sizes:
            bins = self.matrix.distance_to_bins(window_size)
            av_values_left, av_values_right = self._contact_average(bins, offset_bins)
            av_values = (av_values_left + av_values_right)/2
            self.data("av_%d" % window_size, av_values)
            self.data("avl_%d" % window_size, av_values_left)
            self.data("avr_%d" % window_size, av_values_right)

    @calculateondemand
    def average_contacts(self, window_size):
        if window_size is None:
            window_size = self.window_sizes[0]
        return self[:, 'av_%d' % window_size]


class MetaMatrixBase(ArchitecturalFeature, FileGroup):
    _classid = 'METAMATRIXBASE'

    def __init__(self, array=None, regions=None, window_width=50, data_selection=None,
                 file_name=None, mode='a', tmpdir=None,
                 _group_name='meta_base'):

        ArchitecturalFeature.__init__(self)

        if isinstance(array, str) and file_name is None:
            file_name = array
            array = None

        FileGroup.__init__(self, _group_name, file_name, mode=mode, tmpdir=tmpdir)

        try:
            self.meta_matrix = self.file.get_node(self._group, 'meta_matrix')
            self._calculated = True
        except t.NoSuchNodeError:
            self.meta_matrix = None

        self.array = array
        self.window_width = self.array.distance_to_bins(window_width)
        self._data_selection = None
        self._matrix_shape = None
        self.data_selection = data_selection
        self.regions = regions

    @property
    def data_selection(self):
        return self._data_selection

    @data_selection.setter
    def data_selection(self, selection):
        matrix_shape = [self.window_width * 2 + 1, 0]
        if selection is None:
            selection = range(0, len(self.array.data_field_names))

        if isinstance(selection, xrange):
            selection = list(selection)

        if isinstance(selection, int) or isinstance(selection, str):
            selection = [selection]

        if isinstance(selection, list) or isinstance(selection, tuple):
            new_list = []
            for i in selection:
                ix = i
                if isinstance(i, str):
                    ix = self.array.data_field_names.index(i)
                if ix <= len(self.array.data_field_names):
                    new_list.append(ix)

            matrix_shape[1] = len(new_list)
            self._data_selection = new_list
        else:
            raise ValueError("Unsupported data_selection type({})".format(type(selection)))

        self._matrix_shape = tuple(matrix_shape)

    def _sub_matrices(self):
        chromosome_regions = defaultdict(list)

        for i, region in enumerate(self.regions):
            chromosome_regions[region.chromosome].append(((region.start + region.end) / 2, region, i))

        ds = self.data_selection
        try:
            ds = to_slice(ds)
        except ValueError:
            pass

        for chromosome in self.array.chromosomes():
            matrix = self.array.as_matrix(chromosome)
            for pos, region, i in chromosome_regions[chromosome]:
                bin_range = matrix.region_bins(GenomicRegion(start=pos, end=pos, chromosome=chromosome))
                for region_ix in xrange(bin_range.start, bin_range.stop):
                    yield i, region, matrix[region_ix - self.window_width:region_ix + self.window_width + 1, ds]

    @calculateondemand
    def _calculate(self, *args, **kwargs):
        pass


class MetaArray(MetaMatrixBase):
    _classid = 'METAARRAY'

    def __init__(self, array=None, regions=None, window_width=50000, data_selection=None,
                 file_name=None, mode='a', tmpdir=None,
                 _group_name='meta_matrix'):
        MetaMatrixBase.__init__(self, array=array, regions=regions, window_width=window_width,
                                data_selection=data_selection, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                _group_name=_group_name)

    def _calculate(self):
        avg_matrix = np.zeros(self._matrix_shape)
        count_matrix = np.zeros(self._matrix_shape)
        for _, _, m_sub in self._sub_matrices():
            if m_sub.shape == self._matrix_shape:
                count_matrix += np.isnan(m_sub) == False
                m_sub[np.isnan(m_sub)] = 0
                avg_matrix += m_sub

        avg_matrix /= count_matrix

        self.meta_matrix = self.file.create_carray(self._group, 'meta_matrix', t.Float32Atom(),
                                                   tuple(reversed(self._matrix_shape)))
        self.meta_matrix[:] = avg_matrix.T

    @calculateondemand
    def matrix(self):
        return self.meta_matrix[:]

    def x(self):
        x = []
        for i in xrange(-1*self.window_width, self.window_width + 1):
            d = self.array.bins_to_distance(i)
            x.append(d)
        return x

    def y(self):
        y = []
        for i in self.data_selection:
            y.append(self.array.y_values[i])
        return y


class MetaHeatmap(MetaMatrixBase):
    _classid = 'METAHEATMAP'

    def __init__(self, array=None, regions=None, window_width=50000, data_selection=None,
                 file_name=None, mode='a', tmpdir=None,
                 _group_name='meta_heatmap'):

        if data_selection is None and array is not None:
            data_selection = array.data_field_names[0]

        if not isinstance(data_selection, str) and not isinstance(data_selection, int):
            raise ValueError("data_selection parameter must be "
                             "int, str, or None, but is {}".format(type(data_selection)))

        MetaMatrixBase.__init__(self, array=array, regions=regions, window_width=window_width,
                                data_selection=data_selection, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                _group_name=_group_name)

    def _calculate(self):
        order = []
        heatmap = []
        for i, region, m_sub in self._sub_matrices():
            order.append(i)
            if m_sub.shape == self._matrix_shape:
                heatmap.append(m_sub[:, 0])
            else:
                m = np.zeros((self._matrix_shape[0],))
                m[m == 0] = np.nan
                heatmap.append(m)

        heatmap = np.array(heatmap)[order]
        self.meta_matrix = self.file.create_carray(self._group, 'meta_matrix', t.Float32Atom(),
                                                   heatmap.shape)
        self.meta_matrix[:] = heatmap

    @calculateondemand
    def heatmap(self):
        return self.meta_matrix[:]

    def x(self):
        x = []
        for i in xrange(-1*self.window_width, self.window_width + 1):
            d = self.array.bins_to_distance(i)
            x.append(d)
        return x


class ZeroWeightFilter(MatrixArchitecturalRegionFeatureFilter):
    """
    Filter edges where every associated weight is 0.
    """
    def __init__(self, mask=None):
        """
        Initialize filter with chosen parameters.

        :param mask: Optional Mask object describing the mask
                     that is applied to filtered edges.
        """
        MatrixArchitecturalRegionFeatureFilter.__init__(self, mask=mask)

    def valid_edge(self, edge):
        """
        Check if an edge is on (or near) the diagonal of the :class:`~Hic` matrix.
        """
        i = 0
        while True:
            try:
                weight = getattr(edge, 'weight_' + str(i))
                if weight == 0:
                    return False
                i += 1
            except AttributeError:
                return True


class ExpectedObservedCollectionFilter(MatrixArchitecturalRegionFeatureFilter):
    """
    Filter edges where every associated weight is 0.
    """
    def __init__(self, collection, fold_change=1, filter_when_single_invalid=False, mask=None):
        """
        Initialize filter with chosen parameters.

        :param mask: Optional Mask object describing the mask
                     that is applied to filtered edges.
        """
        MatrixArchitecturalRegionFeatureFilter.__init__(self, mask=mask)

        # count number of matrices
        n_hic = 0
        while 'weight_' + str(n_hic) in collection.field_names:
            n_hic += 1

        self.intra_expected = dict()
        self.inter_expected = dict()
        for i in xrange(n_hic):
            print 'weight_' + str(i)
            with ExpectedContacts(collection, weight_column='weight_' + str(i)) as ex:
                self.intra_expected[i] = ex.intra_expected()
                self.inter_expected[i] = ex.inter_expected()

        self.n_hic = n_hic
        self.fold_change = fold_change
        self.filter_single = filter_when_single_invalid
        self.regions_dict = collection.regions_dict

    def valid_edge(self, edge):
        """
        Check if an edge weight is at least fold_change above
        the expected weight for this contact.
        """
        source = edge.source
        sink = edge.sink
        intra = False
        if self.regions_dict[source].chromosome == self.regions_dict[sink].chromosome:
            intra = True
        n_failed = 0
        for i in xrange(self.n_hic):
            if intra:
                expected = self.intra_expected[i][abs(sink-source)]
            else:
                expected = self.inter_expected[i]

            if getattr(edge, 'weight_' + str(i)) < self.fold_change*expected:
                if self.filter_single:
                    return False
                else:
                    n_failed += 1
        if n_failed == self.n_hic:
            return False
        return True


class BackgroundLigationCollectionFilter(MatrixArchitecturalRegionFeatureFilter):
    """
    Filter edges where every associated weight is 0.
    """
    def __init__(self, collection, fold_change=1, filter_when_single_invalid=False,
                 all_contacts=True, mask=None):
        """
        Initialize filter with chosen parameters.

        :param mask: Optional Mask object describing the mask
                     that is applied to filtered edges.
        """
        MatrixArchitecturalRegionFeatureFilter.__init__(self, mask=mask)

        regions_dict = collection.regions_dict
        # count number of matrices
        n_hic = 0
        while 'weight_' + str(n_hic) in collection.field_names:
            n_hic += 1

        inter_count = defaultdict(int)
        inter_sum = defaultdict(int)
        for edge in collection.edges(lazy=True):
            intra = regions_dict[edge.source].chromosome == regions_dict[edge.sink].chromosome
            for i in xrange(n_hic):
                if intra:
                    inter_count[i] += 1
                    inter_sum[i] += getattr(edge, 'weight_' + str(i))

        if all_contacts:
            with PossibleContacts(collection, weight_column='weight_0') as pc:
                for i in xrange(n_hic):
                    inter_count[i] = pc.inter_possible()

        self.cutoff = dict()
        for i in xrange(n_hic):
            if inter_count[i] == 0:
                self.cutoff[i] = 0
            else:
                self.cutoff[i] = fold_change*(inter_sum[i]/inter_count[i])

        self.n_hic = n_hic
        self.filter_single = filter_when_single_invalid

    def valid_edge(self, edge):
        """
        Check if an edge weight is at least fold_change above
        the expected weight for this contact.
        """
        n_failed = 0
        for i in xrange(self.n_hic):
            if getattr(edge, 'weight_' + str(i)) < self.cutoff[i]:
                if self.filter_single:
                    return False
                else:
                    n_failed += 1
        if n_failed == self.n_hic:
            return False
        return True


class MinMaxDistanceCollectionFilter(MatrixArchitecturalRegionFeatureFilter):
    """
    Filter edges where every associated weight is 0.
    """
    def __init__(self, collection, min_distance, max_distance, mask=None):
        """
        Initialize filter with chosen parameters.

        :param mask: Optional Mask object describing the mask
                     that is applied to filtered edges.
        """
        MatrixArchitecturalRegionFeatureFilter.__init__(self, mask=mask)

        self.min_distance_bins = collection.distance_to_bins(min_distance)
        self.max_distance_bins = collection.distance_to_bins(max_distance)
        self.regions_dict = collection.regions_dict

    def valid_edge(self, edge):
        """
        Check if an edge weight is at least fold_change above
        the expected weight for this contact.
        """
        source = edge.source
        sink = edge.sink
        distance_bins = sink-source

        # inter-chromosomal are valid by default
        if self.regions_dict[source].chromosome != self.regions_dict[sink].chromosome:
            return True

        if self.min_distance_bins <= distance_bins <= self.max_distance_bins:
            return True

        return False


class DiagonalCollectionFilter(MatrixArchitecturalRegionFeatureFilter):
    """
    Filter contacts in the diagonal of a :class:`~Hic` matrix.
    """
    def __init__(self, distance=0, mask=None):
        """
        Initialize filter with chosen parameters.

        :param distance: Distance from the diagonal up to which
                         contacts will be filtered
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered edges.
        """
        MatrixArchitecturalRegionFeatureFilter.__init__(self, mask=mask)
        self.distance = distance

    def valid_edge(self, edge):
        """
        Check if an edge is on (or near) the diagonal of the :class:`~Hic` matrix.
        """
        if abs(edge.source-edge.sink) <= self.distance:
            return False
        return True


class BackgroundLigationFilter(HicEdgeFilter):
    """
    Filter a :class:`~HicEdge` if it does not have a weight
    larger than  [fold_change*background ligation frequency].

    Background ligation frequency is estimated as the average
    of all non-zero inter-chromosomal contacts of this Hic object.
    """
    def __init__(self, hic, fold_change=5, all_contacts=False, mask=None):
        """
        Initialize filter with these settings.

        :param hic: The :class:`~Hic` object that this
                    filter will be called on. Needed for
                    contact count calculation.
        :param fold_change: Lowest acceptable edge weight is calculated
                            as fold_change*(inter_sum/inter_count)
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered edges.
        """
        HicEdgeFilter.__init__(self, mask=mask)

        regions_dict = hic.regions_dict

        inter_count = 0
        inter_sum = 0
        for edge in hic.edges(lazy=True):
            if regions_dict[edge.source].chromosome != regions_dict[edge.sink].chromosome:
                inter_count += 1
                inter_sum += edge.weight

        if all_contacts:
            with PossibleContacts(hic) as pc:
                inter_count = pc.inter_possible()

        if inter_count == 0:
            self.cutoff = 0
        else:
            self.cutoff = fold_change*(inter_sum/inter_count)

    def valid_edge(self, edge):
        """
        Check if an edge weight is below background ligation frequency.
        """
        if edge.weight < self.cutoff:
            return False
        return True


class ExpectedObservedEnrichmentFilter(HicEdgeFilter):
    """
    Filter a :class:`~HicEdge` if it does not have a weight
    larger than fold_change times its expected value.
    """
    def __init__(self, hic, fold_change=2, mask=None):
        """
        Initialize filter with these settings.

        :param hic: The :class:`~Hic` object that this
                    filter will be called on. Needed for
                    expected contact count calculation.
        :param fold_change: Lowest acceptable edge weight is calculated
                            as fold_change*expected
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered edges.
        """
        HicEdgeFilter.__init__(self, mask=mask)

        with ExpectedContacts(hic) as ex:
            self.intra_expected = ex.intra_expected()
            self.inter_expected = ex.inter_expected()
        self.regions_dict = hic.regions_dict
        self.fold_change = fold_change

    def valid_edge(self, edge):
        """
        Check if an edge weight is at least fold_change above
        the expected weight for this contact.
        """
        source = edge.source
        sink = edge.sink
        if self.regions_dict[source].chromosome == self.regions_dict[sink].chromosome:
            expected = self.intra_expected[abs(sink-source)]
        else:
            expected = self.inter_expected

        if edge.weight < self.fold_change*expected:
            return False
        return True