from .config import config
from .regions import Chromosome, Genome
from .matrix import RegionMatrixTable, RegionMatrixContainer, Edge, LazyEdge
from abc import abstractmethod, ABCMeta
from future.utils import with_metaclass, string_types, viewitems
from .tools.load import load
from .tools.general import distribute_integer, RareUpdateProgressBar
from .tools.matrix import restore_sparse_rows, remove_sparse_rows
from .general import MaskFilter, MaskedTableView
from collections import defaultdict
import multiprocessing as mp
import threading
import queue
import numpy as np
import warnings
import logging
import msgpack
import copy

logger = logging.getLogger(__name__)
fanc_access_lock = threading.Lock()


def _edge_overlap_split_rao(original_edge, overlap_map):
    """
    Resolve the distribution of contacts when binning using
    Rao and Rowley et al. 2014 approach.
    """
    original_source = original_edge[0]
    original_sink = original_edge[1]
    original_weight = original_edge[2]

    new_source_nodes = overlap_map[original_source]
    new_sink_nodes = overlap_map[original_sink]

    if len(new_source_nodes) == 0:
        return []
    elif len(new_source_nodes) == 1:
        new_source_nodes = [new_source_nodes[0][0]]
    else:
        new_source_nodes = [new_source_nodes[0][0], new_source_nodes[-1][0]]

    if len(new_sink_nodes) == 0:
        return []
    elif len(new_sink_nodes) == 1:
        new_sink_nodes = [new_sink_nodes[0][0]]
    else:
        new_sink_nodes = [new_sink_nodes[0][0], new_sink_nodes[-1][0]]

    edges = {}
    for new_source in new_source_nodes:
        for new_sink in new_sink_nodes:
            if new_source <= new_sink:
                edges[(new_source, new_sink)] = 0
            else:
                edges[(new_sink, new_source)] = 0

    edges_list = []
    try:
        weights = distribute_integer(original_weight, len(edges))
        for i, key_pair in enumerate(edges):
            edges_list.append([key_pair[0], key_pair[1], weights[i]])
    except IndexError:
        pass

    return edges_list


def _get_overlap_map(old_regions, new_regions):
    # 1. organize regions in self by chromosome
    new_region_map = defaultdict(list)
    for i, new_region in enumerate(new_regions):
        new_region_map[new_region.chromosome].append([new_region.start, new_region.end, i])

    # 2. iterate over regions in hic to find overlap
    def _get_overlap(new_region, old_region):
        new_region_length = new_region[1] - new_region[0] + 1
        overlap = min(old_region[1], new_region[1]) - max(old_region[0], new_region[0]) + 1
        return max(0, overlap / new_region_length)

    old_to_new = {}
    current_chromosome = ''
    current_ix = 0
    for i, old_region in enumerate(old_regions):
        old_to_new[i] = []
        if current_chromosome != old_region.chromosome:
            current_ix = 0
            current_chromosome = old_region.chromosome

        if current_chromosome not in new_region_map:
            continue

        found_overlap = True
        while found_overlap:
            found_overlap = False
            if current_ix < len(new_region_map[current_chromosome]):
                new_region = new_region_map[current_chromosome][current_ix]
                overlap = _get_overlap(new_region, [old_region.start, old_region.end, i])
                if overlap > 0:
                    old_to_new[i].append([new_region[2], overlap])
                    current_ix += 1
                    found_overlap = True
                elif old_region.start > new_region[1]:
                    current_ix += 1
                    found_overlap = True

        current_ix -= 1

    return old_to_new


def _bin_hic_partition_worker(hic_file, qin, qout,
                              overlap_map, _edges_by_overlap_method,
                              access_lock):

    try:
        overlap_map = msgpack.loads(overlap_map, strict_map_key=False)
        while True:
            worker_input = qin.get()
            if worker_input is None:
                logger.debug("Received stop signal, worker terminating.")
                break

            partition1, partition2 = worker_input
            logger.debug("Received {}-{}".format(partition1, partition2))

            with access_lock:
                hic = None
                try:
                    hic = load(hic_file)
                    _weight_field = hic._default_score_field
                    all_edges = []
                    key = slice(partition1[0], partition1[1], 1), slice(partition2[0], partition2[1], 1)
                    for edge in hic.edges_dict(key, lazy=True):
                        try:
                            all_edges.append([edge['source'], edge['sink'], edge[_weight_field]])
                        except KeyError:
                            all_edges.append([edge['source'], edge['sink'], edge['weight']])

                finally:
                    if hic is not None:
                        hic.close()

            edges = defaultdict(int)
            for old_source, old_sink, old_weight in all_edges:
                try:
                    for new_source, new_sink, new_weight in _edges_by_overlap_method(
                            [old_source, old_sink, old_weight], overlap_map):
                        if new_weight != 0:
                            edges[(new_source, new_sink)] += new_weight
                except KeyError:
                    warnings.warn("key".format(old_source, old_sink))

            qout.put(msgpack.dumps(edges))
    except Exception as e:
        qout.put(e)


class Hic(RegionMatrixTable):
    """
    Central class for working with Hi-C data.

    This class adds functions for matrix binning and filtering to
    the base class :class:`~RegionMatrixTable`.
    """

    _classid = 'HIC'

    def __init__(self, file_name=None, mode='a', tmpdir=None,
                 partition_strategy='auto',
                 additional_region_fields=None, additional_edge_fields=None,
                 _table_name_regions='regions', _table_name_edges='edges',
                 _edge_buffer_size=config.edge_buffer_size):
        RegionMatrixTable.__init__(self, file_name=file_name,
                                   mode=mode, tmpdir=tmpdir,
                                   additional_region_fields=additional_region_fields,
                                   additional_edge_fields=additional_edge_fields,
                                   partition_strategy=partition_strategy,
                                   _table_name_regions=_table_name_regions,
                                   _table_name_edges=_table_name_edges,
                                   _edge_buffer_size=_edge_buffer_size)

    def load_from_hic(self, hic, threads=1, chromosomes=None,
                      _edges_by_overlap_method=_edge_overlap_split_rao,
                      _regions_soft_max=50000):
        """
        Load data from another :class:`~Hic` object.

        If this object has no associated regions, the regions and
        contacts of the provided object will simply be copied.

        If regions are already present, the contacts of the provided
        matrix will be binned into the regions of this object using the
        overlap method provided.

        :param hic: Another :class:`~Hic` object
        :param threads: Number of parallel processing threads. More threads also
                        means higher memory usage.
        :param _edges_by_overlap_method: A function that maps reads from
                                         one genomic region to others using
                                         a supplied overlap map. By default
                                         it uses the Rao et al. (2014) method.
                                         See :func:`~_edge_overlap_split_rao`
        :param _regions_soft_max: Maximum dimension of each processed submatrix
                                  per thread. This is a soft maximum, which may be
                                  increased as required for very large chromosomes or
                                  small bin sizes
        """
        # if we do not have any nodes in this Hi-C object...
        if len(self.regions) == 0:
            self.add_regions(hic.regions, preserve_attributes=True)
            self.add_edges(hic.edges(lazy=True))

        # if already have nodes in this HiC object...
        else:
            logger.info("Binning Hi-C contacts")

            # create region "overlap map"
            overlap_map = _get_overlap_map(hic.regions(lazy=False), self.regions(lazy=False))

            file_name = hic.file.filename
            m = mp.Manager()
            access_lock = m.Lock()
            qout = m.Queue()
            qin = m.Queue()

            if not isinstance(hic, RegionMatrixTable):
                edge_counter = 0
                with RareUpdateProgressBar(max_value=len(hic.edges), silent=config.hide_progressbars,
                                           prefix="Binning") as pb:
                    if chromosomes is None:
                        chromosomes = hic.chromosomes()
                    for i in range(len(chromosomes)):
                        for j in range(i, len(chromosomes)):
                            logger.debug("Chromosomes: {}-{}".format(chromosomes[i], chromosomes[j]))
                            edges = defaultdict(int)
                            for edge in hic.edges_dict((chromosomes[i], chromosomes[j]), lazy=True, norm=False):
                                old_source, old_sink = edge['source'], edge['sink']
                                try:
                                    old_weight = edge[hic._default_score_field]
                                except KeyError:
                                    old_weight = edge['weight']

                                for new_source, new_sink, new_weight in _edges_by_overlap_method(
                                        [old_source, old_sink, old_weight], overlap_map):
                                    if new_weight != 0:
                                        edges[(new_source, new_sink)] += new_weight

                                edge_counter += 1
                                pb.update(edge_counter)
                            logger.debug("Adding edges {}/{} ({})".format(i, j, len(edges)))
                            for (source, sink), weight in edges.items():
                                self.add_edge_simple(source, sink, weight=weight)
            else:
                pool = None
                try:
                    logger.info("Launching processes")
                    with mp.get_context("spawn").Pool(threads, _bin_hic_partition_worker,
                                                      (file_name,
                                                       qin, qout,
                                                       msgpack.dumps(overlap_map),
                                                       _edges_by_overlap_method,
                                                       access_lock)) as pool:

                        n_fragments = len(hic.regions)
                        thread_max = min(int(n_fragments/threads), _regions_soft_max)
                        partitions = [[0, 0]]
                        previous_chromosome = None
                        for i, region in enumerate(hic.regions(lazy=True)):
                            if previous_chromosome is not None and (region.chromosome != previous_chromosome
                                                                    or i == n_fragments - 1):
                                partition_size = partitions[-1][1] - partitions[-1][0]
                                current_size = i - partitions[-1][1]
                                partition_empty = thread_max - partition_size

                                if partition_size == 0 or partition_empty > current_size/2:
                                    partitions[-1][1] = i
                                else:
                                    partitions.append([partitions[-1][1], i])
                            previous_chromosome = region.chromosome
                        partitions[-1][1] = n_fragments

                        logger.info("Submitting partitions")

                        # submit intra-chromosomal first to distribute
                        # load among workers more evenly
                        n_chunks = 0
                        for partition in partitions:
                            qin.put((partition, partition))
                            n_chunks += 1

                        # then submit inter-chromosomal partitions
                        for cix1, partition1 in enumerate(partitions):
                            for cix2 in range(cix1 + 1, len(partitions)):
                                partition2 = partitions[cix2]
                                qin.put((partition1, partition2))
                                n_chunks += 1

                        self._disable_edge_indexes()
                        logger.info("Collecting results")
                        with RareUpdateProgressBar(max_value=n_chunks, prefix="Binning") as pb:
                            for i in range(n_chunks):
                                out = qout.get(block=True)
                                if isinstance(out, Exception):
                                    raise out
                                edges = msgpack.loads(out, use_list=False, strict_map_key=False)
                                for (source, sink), weight in edges.items():
                                    self.add_edge_simple(source, sink, weight=weight)
                                pb.update(i)
                finally:
                    for i in range(threads):
                        qin.put(None)

                    if pool is not None:
                        pool.terminate()
            logger.debug("Final flush")
            self.flush()

    def bin(self, bin_size, threads=1, chromosomes=None, *args, **kwargs):
        """
        Map edges in this object to equidistant bins.

        :param bin_size: Bin size in base pairs
        :param threads: Number of threads used for binning
        :return: :class:`~Hic` object
        """
        # find chromosome lengths
        logger.info("Constructing binned genome...")
        if chromosomes is None:
            chromosomes = self.chromosomes()
        chromosome_sizes = {chromosome: 0 for chromosome in chromosomes}
        for region in self.regions():
            if region.chromosome not in chromosome_sizes:
                continue
            if chromosome_sizes[region.chromosome] < region.end:
                chromosome_sizes[region.chromosome] = region.end

        chromosome_list = []
        for chromosome in chromosomes:
            chromosome_list.append(Chromosome(name=chromosome, length=self.chromosome_lengths[chromosome]))

        genome = Genome(chromosomes=chromosome_list)
        regions = genome.get_regions(bin_size)
        genome.close()

        logger.info("Binning edges...")
        if 'mode' not in kwargs:
            kwargs['mode'] = 'w'
        hic = self.__class__(*args, **kwargs)
        hic.add_regions(regions.regions(lazy=True), preserve_attributes=False)
        regions.close()
        hic.load_from_hic(self, threads=threads, chromosomes=chromosomes)

        return hic

    def bias_vector(self, vector=None):
        """
        Get or set the vector of region biases in this object.

        This internally sets the "bias" attribute of each region in the
        object.

        :param vector: a numpy array with bias values
        :return: a numpy array with bias values
        """
        if vector is not None:
            self.region_data('bias', vector)

        biases = np.array(list(self.region_data('bias')))
        return biases

    def filter_diagonal(self, distance=0, queue=False):
        """
        Convenience function that applies a :class:`~DiagonalFilter`.

        :param distance: Distance from the diagonal up to which matrix entries
                         will be filtered/removed. The default, 0, filters only the
                         diagonal itself.
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        """
        mask = self.add_mask_description('diagonal',
                                         'Mask the diagonal of the Hic matrix (up to distance %d)' % distance)
        diagonal_filter = DiagonalFilter(self, distance=distance, mask=mask)
        self.filter(diagonal_filter, queue)

    def filter_low_coverage_regions(self, rel_cutoff=None, cutoff=None, queue=False):
        """
        Convenience function that applies a :class:`~LowCoverageFilter`.

        The cutoff can be provided in two ways:
        1. As an absolute threshold. Regions with contact count below this
        absolute threshold are filtered
        2. As a fraction relative to the median contact count of all regions.

        If both is supplied, whichever threshold is lower will be selected.

        If no parameter is supplied, rel_cutoff will be chosen as 0.1.

        :param rel_cutoff: A cutoff as a fraction (0-1) of the median contact count of all
                           regions.
        :param cutoff: A cutoff in absolute contact counts (can be float) below
                       which regions are considered "low coverage"
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        """
        if cutoff is None and rel_cutoff is None:
            rel_cutoff = 0.1

        mask = self.add_mask_description('low_coverage',
                                         'Mask low coverage regions in the Hic matrix '
                                         '(absolute cutoff {:.4}, relative '
                                         'cutoff {:.1%}'.format(float(cutoff) if cutoff else 0.,
                                                                float(rel_cutoff) if rel_cutoff else 0.))

        low_coverage_filter = LowCoverageFilter(self, rel_cutoff=rel_cutoff, cutoff=cutoff, mask=mask)
        self.filter(low_coverage_filter, queue)

    def filter_statistics(self):
        stats = self.mask_statistics(self._edges)
        return stats

    def normalise(self, method='KR', **kwargs):
        if method.lower() == 'kr':
            bias_vector = kr_balancing(self, **kwargs)
        elif method.lower() == 'ice':
            bias_vector = ice_balancing(self, **kwargs)
        elif method.lower() == 'vc' or method.lower() == 'vanilla':
            bias_vector = vanilla_coverage_norm(self, **kwargs)
        elif method.lower() in {'sqrt_vc', 'sqrt-vc', 'sqrt_vanilla', 'sqrt-vanilla',
                                'vc_sqrt', 'vc-sqrt', 'vanilla-sqrt', 'vanilla_sqrt'}:
            bias_vector = sqrt_vanilla_coverage_norm(self, **kwargs)
        else:
            raise ValueError("Unknown normalisation method: {}".format(method))
        return bias_vector


class LegacyHic(Hic):

    _classid = 'ACCESSOPTIMISEDHIC'

    def __init__(self, file_name=None, mode='a', tmpdir=None,
                 partition_strategy='chromosome',
                 additional_region_fields=None, additional_edge_fields=None,
                 _table_name_regions='nodes', _table_name_edges='edges',
                 _edge_buffer_size=config.edge_buffer_size):
        RegionMatrixTable.__init__(self, file_name=file_name,
                                   mode=mode, tmpdir=tmpdir,
                                   additional_region_fields=additional_region_fields,
                                   additional_edge_fields=additional_edge_fields,
                                   partition_strategy=partition_strategy,
                                   _table_name_regions=_table_name_regions,
                                   _table_name_edges=_table_name_edges,
                                   _edge_buffer_size=_edge_buffer_size)
        self._expected_value_cache = None

    def expected_values(self, selected_chromosome=None, norm=True,
                        force=False, *args, **kwargs):
        if not force and self._expected_value_cache is not None:
            intra_expected, chromosome_intra_expected, inter_expected = self._expected_value_cache
        else:
            intra_expected, chromosome_intra_expected, \
                inter_expected = RegionMatrixContainer.expected_values(self, norm=norm)
            self._expected_value_cache = intra_expected, chromosome_intra_expected, inter_expected

        if selected_chromosome is not None:
            return chromosome_intra_expected[selected_chromosome]

        return intra_expected, chromosome_intra_expected, inter_expected

    def marginals(self, masked=True, *args, **kwargs):
        """
        Get the marginals vector of this Hic matrix.

        Sums up all contacts for each bin of the Hi-C matrix.
        Unmappable regoins will be masked in the returned vector unless
        the :code:`masked` parameter is set to :code:`False`.

        By default, corrected matrix entries are summed up.
        To get uncorrected matrix marginals use :code:`norm=False`.
        Generally, all parameters accepted by :func:`~RegionMatrixContainer.edges`
        are supported.

        :param masked: Use a numpy masked array to mask entries
                       corresponding to unmappable regions
        :param kwargs: Keyword arguments passed to :func:`~RegionPairsContainer.edges`
        """
        kwargs.setdefault('lazy', True)
        row_regions, col_regions, edges_iter = self.regions_and_matrix_entries(*args, **kwargs)
        min_ix = min(row_regions[0].ix, col_regions[0].ix)
        max_ix = max(row_regions[-1].ix, col_regions[-1].ix)

        marginals = np.zeros(max_ix - min_ix + 1)

        logger.debug("Calculating marginals...")
        for i, (source, sink, weight) in enumerate(edges_iter):
            if source <= sink:
                marginals[source] += weight
            if source < sink:
                marginals[sink] += weight

        if masked:
            mask = np.zeros(len(marginals), dtype=bool)
            for i in range(len(marginals)):
                mask[i] = marginals[i] == 0
            marginals = np.ma.masked_where(mask, marginals)

        return marginals

    def mappable(self):
        return self.marginals() > 0

    def bias_vector(self):
        bias = [row['bias'] for row in self.file.get_node('/', 'node_annot').iterrows()]
        if len(bias) == len(self.regions):
            return np.array(bias)
        else:
            return np.repeat(1., len(self.regions))


class HicEdgeFilter(with_metaclass(ABCMeta, MaskFilter)):
    """
    Abstract class that provides filtering functionality for the
    edges/contacts in a :class:`~Hic` object.

    Extends MaskFilter and overrides valid(self, row) to make
    :class:`~HicEdge` filtering more "natural".

    To create custom filters for the :class:`~Hic` object, extend this
    class and override the valid_edge(self, edge) method.
    valid_edge should return False for a specific :class:`~HicEdge` object
    if the object is supposed to be filtered/masked and True
    otherwise. See :class:`~DiagonalFilter` for an example.

    Pass a custom filter to the :func:`~Hic.filter` method in :class:`~Hic`
    to apply it.
    """

    def __init__(self, hic=None, mask=None):
        """
        Initialize HicEdgeFilter.

        :param mask: The Mask object that should be used to mask
                     filtered :class:`~HicEdge` objects. If None the default
                     Mask will be used.
        """
        super(HicEdgeFilter, self).__init__(mask)
        self._lazy_edge = LazyEdge(None)
        self._hic = hic

    @abstractmethod
    def valid_edge(self, edge):
        """
        Determine if a :class:`~HicEdge` object is valid or should
        be filtered.

        When implementing custom HicEdgeFilter this method must be
        overridden. It should return False for :class:`~HicEdge` objects that
        are to be fitered and True otherwise.

        Internally, the :class:`~Hic` object will iterate over all HicEdge
        instances to determine their validity on an individual
        basis.

        :param edge: A :class:`~HicEdge` object
        :return: True if :class:`~HicEdge` is valid, False otherwise
        """
        pass

    def set_hic_object(self, hic_object):
        """
        Set the :class:`~Hic` instance to be filtered by this
        HicEdgeFilter.

        Used internally by :class:`~Hic` instance.

        :param hic_object: :class:`~Hic` object
        """
        self._hic = hic_object

    def valid(self, row):
        """
        Map valid_edge to MaskFilter.valid(self, row).

        :param row: A pytables Table row.
        :return: The boolean value returned by valid_edge.
        """
        edge = self._hic._row_to_edge(row, lazy_edge=self._lazy_edge)
        return self.valid_edge(edge)


class DiagonalFilter(HicEdgeFilter):
    """
    Filter contacts in the diagonal of a :class:`~Hic` matrix.
    """
    def __init__(self, hic, distance=0, mask=None):
        """
        Initialize filter with chosen parameters.

        :param distance: Distance from the diagonal up to which
                         contacts will be filtered
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered edges.
        """
        HicEdgeFilter.__init__(self, mask=mask)
        self.set_hic_object(hic)
        self.distance = distance

    def valid_edge(self, edge):
        """
        Check if an edge is on (or near) the diagonal of the :class:`~Hic` matrix.
        """
        if abs(edge.source-edge.sink) <= self.distance:
            return False
        return True


class LowCoverageFilter(HicEdgeFilter):
    """
    Filter a :class:`~HicEdge` if it connects a region that
    does not have a contact count larger than a specified
    cutoff.

    If the cutoff is not provided, it is automatically
    chosen at 10% of the mean contact count of all regions.
    """
    def __init__(self, hic_object, cutoff=None, rel_cutoff=None, mask=None):
        """
        Initialize filter with these settings.

        The cutoff can be provided in two ways:
        1. As an absolute threshold. Regions with contact count below this
        absolute threshold are filtered
        2. As a fraction relative to the median contact count of all regions.

        If both is supplied, whichever threshold is lower will be selected.

        :param hic_object: The :class:`~Hic` object that this
                           filter will be called on. Needed for
                           contact count calculation.
        :param rel_cutoff: A cutoff as a fraction (0-1) of the median contact count of all
                           regions. If cutoff and rel_cutoff are None, will be set to 10%
        :param cutoff: A cutoff in absolute contact counts (can be float) below
                       which regions are considered "low coverage"
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered edges.
        """
        HicEdgeFilter.__init__(self, mask=mask)
        self.set_hic_object(hic_object)

        self._marginals = hic_object.marginals(norm=False, masked=False)
        if cutoff is None and rel_cutoff is None:
            rel_cutoff = 0.1
            logger.info("Using default 10 percent relative coverage as cutoff")

        if cutoff is not None and rel_cutoff is not None:
            cutoff = min(cutoff if cutoff else float("inf"),
                         self.calculate_cutoffs(rel_cutoff)[0] if rel_cutoff else float("inf"))
        elif rel_cutoff is not None:
            cutoff = self.calculate_cutoffs(rel_cutoff)[0]
        logger.info("Final absolute cutoff threshold is {:.4}".format(float(cutoff)))

        self._regions_to_mask = set()
        for i, contacts in enumerate(self._marginals):
            if contacts < cutoff:
                self._regions_to_mask.add(i)
        logger.info("Selected a total of {} ({:.1%}) regions to be masked".format(
            len(self._regions_to_mask), len(self._regions_to_mask)/len(hic_object.regions)))

    def calculate_cutoffs(self, fraction_threshold=0.05):
        lower = np.nanmedian(self._marginals[self._marginals > 0])*fraction_threshold
        upper = np.nanmedian(self._marginals[self._marginals > 0])+lower
        return lower, upper

    def valid_edge(self, edge):
        """
        Check if an edge falls into a low-coverage region.
        """
        if edge.source in self._regions_to_mask:
            return False
        if edge.sink in self._regions_to_mask:
            return False
        return True


def ice_balancing(hic, tolerance=1e-2, max_iterations=500, whole_matrix=True,
                  inter_chromosomal=True, intra_chromosomal=True, restore_coverage=False,
                  sqrt=True):
    """
    Apply ICE balancing to Hi-C matrices.

    Iteratively calculates and divides by the matrix margins.

    :param hic: Hi-C object
    :param tolerance: Error tolerance (marginal error)
    :param max_iterations: Maximum number of iterations to perform
                           to achieve error tolerance
    :param whole_matrix: Correct the whole matrix at once.
                         Default is to correct each chromosome individually.
    :param inter_chromosomal: Include inter-chromosomal contacts in balancing (only whole matrix)
    :param intra_chromosomal: Include intra-chromosomal contacts in balancing (only whole matrix)
    :param restore_coverage: Restore the matrix to its original coverage after balancing,
                             i.e. the sum of contacts in the matrix after balancing remains
                             (roughly) the same
    :return: bias vector
    """
    logger.info("Starting ICE matrix balancing")

    if not whole_matrix:
        bias_vectors = []
        for chromosome in hic.chromosomes():
            logger.debug("Chromosome {}".format(chromosome))
            region_converter = dict()
            bias_vector = []
            n_regions = 0
            for i, region in enumerate(hic.regions(chromosome)):
                region_converter[region.ix] = i
                bias_vector.append(1)
                if region.valid:
                    n_regions += 1
            bias_vector = np.array(bias_vector, dtype='float64')

            marginal_error = tolerance + 1
            current_iteration = 0

            total_weight = 0
            edges = []
            for e in hic.edges_dict((chromosome, chromosome), lazy=True, norm=False):
                source, sink, weight = e['source'], e['sink'], e['weight']
                edges.append([source, sink, weight])
                total_weight += weight
                if source != sink:
                    total_weight += weight

            if len(edges) > 0:
                while (marginal_error > tolerance and
                       current_iteration < max_iterations):
                    m = np.zeros(len(bias_vector), dtype='float64')
                    for source, sink, weight in edges:
                        source_sub = region_converter[source]
                        sink_sub = region_converter[sink]
                        m[source_sub] += weight
                        if source != sink:
                            m[sink_sub] += weight

                    try:
                        marginal_error = _marginal_error(m)
                    except Exception as e:
                        print(e)
                        print(edges)
                        raise

                    if sqrt:
                        m = np.sqrt(m)
                    else:
                        # multiply with constant factor so marginals are 1
                        bias_mean = np.mean(m[m != 0])
                        marginal_mean = np.sqrt(np.sum(m) / n_regions)
                        m = m * marginal_mean / bias_mean

                    bias_vector *= m

                    for i in range(len(edges)):
                        source = region_converter[edges[i][0]]
                        sink = region_converter[edges[i][1]]
                        edges[i][2] = 0 if m[sink] == 0 else edges[i][2] / m[source] / m[sink]

                    current_iteration += 1

                    if restore_coverage:
                        bias_vector = bias_vector / np.sqrt(total_weight / n_regions)

                    logger.debug("Iteration: %d, error: %lf" % (current_iteration, marginal_error))
            else:
                warnings.warn("Chromosome {} has no valid edges, skipping normalisation!".format(chromosome))

            bias_vectors.append(bias_vector)
        logger.info("Done.")
        logger.info("Adding bias vector...")
        bias_vector = np.concatenate(bias_vectors)

        logger.info("Done.")
    else:
        bias_vector = np.ones(len(hic.regions), float)
        marginal_error = tolerance + 1
        current_iteration = 0
        logger.info("Collecting edges")

        total_weight = 0
        edges = []
        for e in hic.edges(norm=False, lazy=True,
                           intra_chromosomal=intra_chromosomal,
                           inter_chromosomal=inter_chromosomal):
            source, sink, weight = e.source, e.sink, e.weight
            edges.append([source, sink, weight])
            total_weight += weight
            if source != sink:
                total_weight += weight

        logger.info("Starting iterations")
        while (marginal_error > tolerance and
               current_iteration < max_iterations):
            m = np.zeros(len(bias_vector), dtype='float64')
            for i in range(len(edges)):
                source = edges[i][0]
                sink = edges[i][1]
                m[source] += edges[i][2]
                if source != sink:
                    m[sink] += edges[i][2]

            bias_vector *= np.sqrt(m)
            marginal_error = _marginal_error(m)
            for i in range(len(edges)):
                source = edges[i][0]
                sink = edges[i][1]
                edges[i][2] = 0 if m[sink] == 0 else edges[i][2]/np.sqrt(m[source])/np.sqrt(m[sink])

            current_iteration += 1
            logger.debug("Iteration: %d, error: %lf" % (current_iteration, marginal_error))

        if restore_coverage:
            bias_vector = bias_vector / np.sqrt(total_weight / len(hic.regions))

    with np.errstate(divide='ignore'):
        bias_vector = 1/bias_vector
    bias_vector[~np.isfinite(bias_vector)] = 0

    hic.region_data('bias', bias_vector)
    return bias_vector


def _marginal_error(marginals, percentile=99.9):
    marginals = marginals[marginals != 0]
    error = np.percentile(np.abs(marginals - marginals.mean()), percentile)
    return error / marginals.mean()


def vanilla_coverage_norm(*args, **kwargs):
    """
    Apply vanilla coverage normalisation to Hi-C matrices.

    Identical to ice_balancing with max_iterations set to 1 and sqrt to False.

    :param args: see ice_balancing
    :param kwargs: ice_balancing
    :return: bias vector (numpy)
    """
    kwargs['max_iterations'] = 1
    kwargs['sqrt'] = False
    return ice_balancing(*args, **kwargs)


def sqrt_vanilla_coverage_norm(*args, **kwargs):
    """
    Apply vanilla coverage normalisation to Hi-C matrices with sqrt bias vectors.

    Identical to ice_balancing with max_iterations set to 1.

    :param args: see ice_balancing
    :param kwargs: ice_balancing
    :return: bias vector (numpy)
    """
    kwargs['max_iterations'] = 1
    return ice_balancing(*args, **kwargs)


def kr_balancing(hic, whole_matrix=True, intra_chromosomal=True, inter_chromosomal=True,
                 restore_coverage=False):

    if not whole_matrix:
        bias_vectors = []
        for chromosome in hic.chromosomes():
            m = hic.matrix((chromosome, chromosome), norm=False)
            m_corrected, bias_vector_chromosome = correct_matrix(m, restore_coverage=restore_coverage)
            bias_vectors.append(bias_vector_chromosome)
        bias_vector = np.concatenate(bias_vectors)
    else:
        logger.debug("Fetching whole genome matrix")
        m = hic.matrix(norm=False)
        cb = hic.chromosome_bins
        if not intra_chromosomal:
            for chromosome, bins in cb.items():
                m[bins[0]:bins[1], bins[0]:bins[1]] = 0
        if not inter_chromosomal:
            for chromosome, bins in cb:
                m[0:bins[0], bins[0]:bins[1]] = 0
                m[bins[1]:, bins[0]:bins[1]] = 0

        m_corrected, bias_vector = correct_matrix(m, restore_coverage=restore_coverage)

    hic.region_data('bias', bias_vector)
    return bias_vector


def correct_matrix(m, max_attempts=50, restore_coverage=False):
    # remove zero-sum rows
    removed_rows = []
    m_nonzero, ixs = remove_sparse_rows(m, cutoff=0)
    removed_rows.append(ixs)

    has_errors = True
    iterations = 0
    x = None
    while has_errors:
        has_errors = False

        try:
            x = get_bias_vector(m_nonzero)
        except ValueError as e:
            logger.debug("Matrix balancing failed (this can happen!), \
                          removing sparsest rows to try again. Error: \
                          %s" % str(e))
            m_nonzero, ixs = remove_sparse_rows(m_nonzero)
            removed_rows.append(ixs)
            has_errors = True

        iterations += 1
        if iterations > max_attempts:
            raise RuntimeError("Exceeded maximum attempts (%d)" % max_attempts)

    if restore_coverage:
        x = x*np.sqrt(np.sum(m_nonzero)/m_nonzero.shape[0])

    logger.debug("Applying bias vector")
    m_nonzero = x*m_nonzero*x[:, np.newaxis]

    logger.debug(removed_rows)
    logger.debug("Restoring {} sets ({} total) sparse rows.".format(
        len(removed_rows), sum(len(x) for x in removed_rows)))
    # restore zero rows
    m_nonzero = restore_sparse_rows(m_nonzero, removed_rows)
    x = restore_sparse_rows(x, removed_rows)

    return m_nonzero, x


def get_bias_vector(A, x0=None, tol=1e-06, delta=0.1, Delta=3, fl=0, high_precision=False, outer_limit=300):
    logger.debug("Starting matrix balancing")

    with warnings.catch_warnings():
        warnings.filterwarnings('error')

        try:
            # basic variables
            # n=size_(A,1)
            if not isinstance(A, np.ndarray):
                try:
                    if high_precision:
                        A = np.array(A, dtype=np.float128)
                    else:
                        A = np.array(A)
                except AttributeError:
                    A = np.array(A)
            n = A.shape[0]
            # e=ones_(n,1)
            try:
                if high_precision:
                    e = np.ones(n, dtype=np.float128)
                else:
                    e = np.ones(n)
            except AttributeError:
                e = np.ones(n)

            if not x0:
                try:
                    if high_precision:
                        x0 = np.ones(n, np.float128)
                    else:
                        x0 = np.ones(n)
                except AttributeError:
                    x0 = np.ones(n)
            else:
                try:
                    if high_precision:
                        x0 = np.array(x0, np.float128)
                    else:
                        x0 = np.array(x0)
                except AttributeError:
                    x0 = np.array(x0)
            res = np.array([])
            g = 0.9
            etamax = 0.1
            eta = etamax
            stop_tol = tol * 0.5
            x = x0.copy()
            rt = tol ** 2
            # v=x.dot((A * x))
            v = x*A.dot(x)
            rk = 1 - v
            # rho_km1=rk.T * rk
            rho_km1 = rk.T.dot(rk)
            rout = rho_km1
            rho_km2 = rho_km1
            rold = rout
            MVP = 0
            i = 0

            n_iterations_outer = 0
            while rout > rt:
                n_iterations_outer += 1

                if n_iterations_outer > outer_limit:
                    raise ValueError("Number of iterations has exceeded the limit (%d)." % outer_limit)

                i += 1
                k = 0
                y = e.copy()
                innertol = max(eta ** 2 * rout, rt)
                n_iterations_inner = 0
                while rho_km1 > innertol:
                    n_iterations_inner += 1

                    k += 1
                    if k == 1:
                        try:
                            Z = rk / v
                        except Warning:
                            raise ValueError("v=0; Remove zero or sparse rows")
                        p = Z.copy()
                        rho_km1 = rk.T.dot(Z)
                    else:
                        beta = rho_km1 / rho_km2
                        p = Z + beta * p
                    # w = x.*(A*(x.*p)) + v.*p;
                    w = x*A.dot(x*p) + v*p
                    alpha = rho_km1 / p.T.dot(w)
                    ap = alpha * p
                    ynew = y + ap
                    if min(ynew) <= delta:
                        if delta == 0:
                            break
                        ind = np.where(ap < 0)[0]
                        # gamma = min((delta  - y(ind))./ap(ind));
                        gamma = min((delta-y[ind])/ap[ind])
                        y = y + gamma * ap
                        break
                    if max(ynew) >= Delta:
                        ind = np.where(ynew > Delta)[0]
                        gamma = min((Delta-y[ind])/ap[ind])
                        y = y + gamma * ap
                        break
                    y = ynew.copy()
                    rk = rk - alpha * w
                    rho_km2 = rho_km1.copy()
                    Z = rk / v
                    rho_km1 = rk.T.dot(Z)

                try:
                    x = x*y
                except Warning:
                    raise ValueError("Value in x or y too small to represent numerically. Try removing sparse rows")

                v = x*A.dot(x)
                rk = 1 - v
                rho_km1 = rk.T.dot(rk)
                rout = rho_km1.copy()
                MVP = MVP + k + 1
                rat = rout / rold
                rold = rout.copy()
                res_norm = np.sqrt(rout)
                eta_o = eta
                eta = g * rat
                if g * eta_o ** 2 > 0.1:
                    eta = max(eta, g * eta_o ** 2)
                eta = max(min(eta, etamax), stop_tol / res_norm)
                if fl == 1:
                    res = np.array([[res], [res_norm]])

            logger.debug("Matrix-vector products = %d\n" % MVP)
            logger.debug("Outer iterations: %d" % n_iterations_outer)
        except Warning as e:
            logger.error(str(e))
            raise ValueError("Generic catch all warnings")

    return x
