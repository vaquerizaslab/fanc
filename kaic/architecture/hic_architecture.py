from __future__ import division
from kaic.architecture.architecture import TableArchitecturalFeature, calculateondemand
from kaic.architecture.genome_architecture import MatrixArchitecturalRegionFeature, VectorArchitecturalRegionFeature, \
    MatrixArchitecturalRegionFeatureFilter
from kaic.data.genomic import GenomicRegion, HicEdgeFilter, Edge
from collections import defaultdict
import numpy as np
import tables as t
import logging
logging.basicConfig(level=logging.INFO)


class HicArchitecture(object):
    def __init__(self, hic):
        self.hic = hic

    def expected_contacts(self, regions=None):
        return ExpectedContacts(self.hic, smooth=False, regions=regions)

    def possible_contacts(self, regions=None):
        return PossibleContacts(self.hic, regions=regions)

    def directionality_index(self, window_sizes=2000000):
        if isinstance(window_sizes, int):
            window_sizes = (window_sizes, )
        return DirectionalityIndex(self.hic, window_sizes=window_sizes)


class HicEdgeCollection(MatrixArchitecturalRegionFeature):
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
                field_descriptions = hic._edges.coldescrs
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
    def __init__(self, hic, file_name=None, mode='a', tmpdir=None, smooth=True, min_reads=400,
                 regions=None, weight_column='weight', _table_name='expected_contacts'):
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
        self.weight_column = weight_column

    def _calculate(self):
        """
        Get intra- and inter-chromosomal expected contact counts.
        """

        # extract mappable regions from Hic object
        marginals = self.hic.marginals(weight_column=self.weight_column)
        regions = []
        region_ix = set()
        if self.regions is None:
            for region in self.hic.regions(lazy=False):
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
            max_distance = max(max_distance, max_region_by_chromosome[chromosome]-min_region_by_chromosome[chromosome])

        # get the number of pixels at a given bin distance
        pixels_by_distance = np.zeros(max_distance + 1)
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
        reads_by_distance = np.zeros(max_distance + 1)
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

        self.data('distance', np.arange(len(reads_by_distance))*self.hic.bin_size)

        # return here if smoothing not requested
        if not self.smooth:
            intra_expected = reads_by_distance/pixels_by_distance

            self.data('intra', intra_expected)
            self.data('contacts', reads_by_distance)
            self.data('pixels', pixels_by_distance)
            self._table.attrs['inter'] = inter_expected
            return

        # smoothing
        smoothed_reads_by_distance = np.zeros(max_distance + 1)
        smoothed_pixels_by_distance = np.zeros(max_distance + 1)
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


class PossibleContacts(TableArchitecturalFeature):
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


class DirectionalityIndex(VectorArchitecturalRegionFeature):
    def __init__(self, hic, file_name=None, mode='a', tmpdir=None,
                 window_sizes=(2000000,), _table_name='directionality_index'):

        # are we retrieving an existing object?
        if isinstance(hic, str) and file_name is None:
            file_name = hic
            hic = None
            VectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir)
        else:
            di_fields = {}
            self.window_sizes = []
            for i, window_size in enumerate(window_sizes):
                di_fields['di_%d' % window_size] = t.Float32Col(pos=i)
                self.window_sizes.append(window_size)
            VectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                      data_fields=di_fields, regions=hic.regions,
                                                      _table_name_data=_table_name)

        self.window_sizes = []
        for colname in self._regions.colnames:
            if colname.startswith("di_"):
                window_size = int(colname[3:])
                self.window_sizes.append(window_size)

        self.hic = hic

    def _get_boundary_distances(self):
        n_bins = len(self.regions)
        # find distances to chromosome boundaries in bins
        boundary_dist = np.zeros(n_bins, dtype=int)
        last_chromosome = None
        last_chromosome_index = 0
        for i, region in enumerate(self.regions(lazy=True)):
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
        bin_size = self.hic.bin_size
        bin_window_size = int(window_size/bin_size)
        if window_size % bin_size > 0:
            bin_window_size += 1

        n_bins = len(self.regions())
        boundary_dist = self._get_boundary_distances()

        left_sums = np.zeros(n_bins)
        right_sums = np.zeros(n_bins)
        directionality_index = np.zeros(n_bins)
        for edge in self.hic.edges(lazy=True):
            source = edge.source
            sink = edge.sink
            weight = edge.weight
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
        return directionality_index

    def _calculate(self):
        for window_size in self.window_sizes:
            directionality_index = self._directionality_index(window_size)
            self.data("di_%d" % window_size, directionality_index)

    @calculateondemand
    def directionality_index(self, window_size=None):
        if window_size is None:
            window_size = self.window_sizes[0]
        return self[:, 'di_%d' % window_size]


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