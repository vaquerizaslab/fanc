"""
TODO

update mappable region handling
update expected value handling
update normalisation handling

"""

import logging
import os
import warnings

import numpy as np
import tables
from genomic_regions import RegionBased, GenomicRegion, as_region
import intervaltree

from .config import config
from .regions import LazyGenomicRegion, RegionsTable, RegionBasedWithBins
from .tools.general import RareUpdateProgressBar, ranges, create_col_index, range_overlap
from .general import Maskable, MaskedTable

from collections import defaultdict
from future.utils import string_types

from bisect import bisect_right

logger = logging.getLogger(__name__)


class Edge(object):
    """
    A contact / an Edge between two genomic regions.

    .. attribute:: source

        The index of the "source" genomic region. By convention,
        source <= sink.

    .. attribute:: sink

        The index of the "sink" genomic region.

    .. attribute:: weight

        The weight or contact strength of the edge. Can, for
        example, be the number of reads mapping to a contact.
    """
    def __init__(self, source, sink, **kwargs):
        """
        :param source: The index of the "source" genomic region
                       or :class:`~Node` object.
        :param sink: The index of the "sink" genomic region
                     or :class:`~Node` object.
        :param data: The weight or of the edge or a dictionary with
                     other fields
        """
        self._source = source
        self._sink = sink
        self.field_names = []

        for key, value in kwargs.items():
            setattr(self, key.decode() if isinstance(key, bytes) else key, value)
            self.field_names.append(key)

    @property
    def source(self):
        try:
            return self._source.ix
        except AttributeError:
            return self._source

    @property
    def sink(self):
        try:
            return self._sink.ix
        except AttributeError:
            return self._sink

    @property
    def source_node(self):
        if isinstance(self._source, GenomicRegion):
            return self._source
        raise RuntimeError("Source not not provided during object initialization!")

    @property
    def sink_node(self):
        if isinstance(self._sink, GenomicRegion):
            return self._sink
        raise RuntimeError("Sink not not provided during object initialization!")

    def __repr__(self):
        base_info = "{}--{}".format(self.source, self.sink)
        for field in self.field_names:
            base_info += "; {}: {}".format(field, str(getattr(self, field)))
        return base_info


class LazyEdge(Edge):
    def __init__(self, row, nodes_table=None):
        self._row = row
        self._nodes_table = nodes_table

    def __getattr__(self, item):
        try:
            return self._row[item]
        except KeyError:
            raise AttributeError("No such attribute {}".format(item))

    def set_row_field(self, key, value, auto_update=True):
        self._row[key] = value
        if auto_update:
            self._row.update()

    @property
    def source_node(self):
        if self._nodes_table is None:
            raise RuntimeError("Must set the _nodes_table attribute before calling this method!")

        if self._source_node is None:
            source_row = self._nodes_table[self.source]
            return LazyGenomicRegion(source_row)
        return self._source_node

    @property
    def sink_node(self):
        if self._nodes_table is None:
            raise RuntimeError("Must set the _nodes_table attribute before calling this method!")

        if self._sink_node is None:
            sink_row = self._nodes_table[self.sink]
            return LazyGenomicRegion(sink_row)
        return self._sink_node

    def __repr__(self):
        base_info = "{}--{}".format(self.source, self.sink)
        return base_info


def as_edge(edge):
    if isinstance(edge, Edge):
        return edge

    try:
        return Edge(**edge)
    except TypeError:
        pass

    if isinstance(edge, tuple) and len(edge) > 1 and \
            isinstance(edge[0], GenomicRegion) and isinstance(edge[1], GenomicRegion):
        try:
            source, sink = edge[0].ix, edge[1].ix
            if len(edge) > 2:
                return Edge(source, sink, weight=edge[2])
            return Edge(source, sink)
        except AttributeError:
            pass

    try:
        source, sink = edge[0], edge[1]
        try:
            weight = edge[2]
        except IndexError:
            weight = 1

        return Edge(source, sink, weight=weight)
    except (TypeError, IndexError):
        pass

    try:
        weight = getattr(edge, 'weight', None)
        if weight is not None:
            return Edge(edge.source, edge.sink, weight=weight)
        return Edge(edge.source, edge.sink)
    except AttributeError:
        pass

    raise ValueError("{} of type {} not recognised as edge "
                     "/ contact!".format(edge, type(edge)))


class RegionPairsContainer(RegionBased):

    def __init__(self):
        RegionBased.__init__(self)
        self._default_value = 1.0
        self._default_score_field = None

    def _add_edge(self, edge, *args, **kwargs):
        raise NotImplementedError("Subclass must override this function")

    def _edges_iter(self, *args, **kwargs):
        raise NotImplementedError("Subclass must implement _edges_iter "
                                  "to enable iterating over edges!")

    def _edges_subset(self, key=None, row_regions=None, col_regions=None,
                      *args, **kwargs):
        raise NotImplementedError("Subclass must implement _edges_subset "
                                  "to enable iterating over edge subsets!")

    def _edges_length(self):
        return sum(1 for _ in self.edges)

    def _edges_getitem(self, item, row_regions=None, col_regions=None,
                       *args, **kwargs):
        raise NotImplementedError("Subclass must implement _edges_getitem "
                                  "to enable getting specific edges!")

    def _key_to_regions(self, key, *args, **kwargs):
        if isinstance(key, tuple):
            if len(key) == 2:
                row_key, col_key = key
            else:
                raise ValueError("Cannot retrieve edge table rows using key {}".format(key))
        else:
            row_key = key
            col_key = None

        if isinstance(row_key, list) and isinstance(row_key[0], GenomicRegion):
            row_regions = row_key
        else:
            row_regions = self.regions(row_key, *args, **kwargs)

        if isinstance(col_key, list) and isinstance(col_key[0], GenomicRegion):
            col_regions = col_key
        else:
            col_regions = self.regions(col_key, *args, **kwargs)

        return row_regions, col_regions

    def _min_max_region_ix(self, regions):
        min_ix = len(self.regions)
        max_ix = 0
        for region in regions:
            min_ix = min(min_ix, region.ix)
            max_ix = max(max_ix, region.ix)
        return min_ix, max_ix

    def __len__(self):
        return self._edges_length()

    def __getitem__(self, item):
        return self._edges_getitem(item)

    def __iter__(self):
        return self.edges()

    def add_contact(self, contact, *args, **kwargs):
        return self.add_edge(contact, *args, **kwargs)

    def add_edge(self, edge, check_nodes_exist=True, *args, **kwargs):
        """
        Add an edge to this object.

        :param edge: :class:`~Edge`, dict with at least the
                     attributes source and sink, optionally weight,
                     or a list of length 2 (source, sink) or 3
                     (source, sink, weight).
        :param check_nodes_exist: Make sure that there are nodes
                                  that match source and sink indexes
        """
        edge = as_edge(edge)

        if check_nodes_exist:
            n_regions = len(self.regions)
            if edge.source >= n_regions or edge.sink >= n_regions:
                raise ValueError("Node index ({}/{}) exceeds number of nodes ({}) in object".format(
                    edge.source, edge.sink, n_regions
                ))

        self._add_edge(edge, *args, **kwargs)

    def add_edges(self, edges, *args, **kwargs):
        """
        Bulk-add edges from a list.

        :param edges: List (or iterator) of edges. See
                      :func:`~RegionMatrixTable.add_edge`
                      for details
        """
        for edge in edges:
            self.add_edge(edge, *args, **kwargs)

    def add_contacts(self, contacts, *args, **kwargs):
        return self.add_edges(contacts, *args, **kwargs)

    @property
    def edges(self):
        """
        Iterate over :class:`~Edge` objects.

        :return: Iterator over :class:`~Edge`
        """

        class EdgeIter(object):
            def __init__(self, this):
                self._regions_pairs = this

            def __getitem__(self, item):
                return self._regions_pairs._edges_getitem(item)

            def __iter__(self):
                return self()

            def __call__(self, key=None, *args, **kwargs):
                norm = kwargs.pop("norm", True)
                intra_chromosomal = kwargs.pop("intra_chromosomal", True)
                inter_chromosomal = kwargs.pop("inter_chromosomal", True)

                row_regions, col_regions = self._regions_pairs._key_to_regions(key)
                if isinstance(row_regions, GenomicRegion):
                    row_regions = [row_regions]
                else:
                    row_regions = list(row_regions)

                if isinstance(col_regions, GenomicRegion):
                    col_regions = [col_regions]
                else:
                    col_regions = list(col_regions)

                regions = dict()
                for rr in (row_regions, col_regions):
                    for region in rr:
                        regions[region.ix] = region

                if key is None:
                    edge_iter = self._regions_pairs._edges_iter(*args, **kwargs)
                else:
                    edge_iter = self._regions_pairs._edges_subset(key, row_regions, col_regions,
                                                                  *args, **kwargs)

                weight_field = kwargs.pop('weight_field', self._regions_pairs._default_score_field)
                bias_field = kwargs.pop('bias_field', 'bias')

                for edge in edge_iter:
                    row_region = regions[edge.source]
                    col_region = regions[edge.sink]
                    if not intra_chromosomal and row_region.chromosome == col_region.chromosome:
                        continue
                    if not inter_chromosomal and row_region.chromosome != col_region.chromosome:
                        continue

                    try:
                        if not row_region.valid or not col_region.valid:
                            continue
                    except AttributeError:
                        pass

                    try:
                        weight = getattr(edge, weight_field)
                        if norm:
                            try:
                                row_bias = getattr(regions[edge.source], bias_field, 1.0)
                                col_bias = getattr(regions[edge.sink], bias_field, 1.0)
                                bias = row_bias * col_bias
                            except AttributeError:
                                bias = 1.0
                        else:
                            bias = 1.0

                        setattr(edge, weight_field, weight * bias)
                    except (TypeError, AttributeError):
                        pass

                    yield edge

            def __len__(self):
                return self._regions_pairs._edges_length()

        return EdgeIter(self)

    def edge_subset(self, key=None, *args, **kwargs):
        """
        Get a subset of edges.

        :param key: Possible key types are:

                    Region types

                    - Node: Only the ix of this node will be used for
                      identification
                    - GenomicRegion: self-explanatory
                    - str: key is assumed to describe a genomic region
                      of the form: <chromosome>[:<start>-<end>:[<strand>]],
                      e.g.: 'chr1:1000-54232:+'

                    Node types

                    - int: node index
                    - slice: node range

                    List types

                    - list: This key type allows for a combination of all
                      of the above key types - the corresponding matrix
                      will be concatenated


                    If the key is a 2-tuple, each entry will be treated as the
                    row and column key, respectively,
                    e.g.: 'chr1:0-1000, chr4:2300-3000' will extract the Hi-C
                    map of the relevant regions between chromosomes 1 and 4.
        :return: generator (:class:`~Edge`)
        """
        return self.edges(key, *args, **kwargs)

    @classmethod
    def merge(cls, pairs, *args, **kwargs):
        if 'mode' not in kwargs:
            kwargs['mode'] = 'w'
        merged_pairs = cls(*args, **kwargs)

        pairs = [pair_object for pair_object in pairs]

        regions = list(pairs[0].regions)
        for pair_object in pairs[1:]:
            try:
                for r1, r2 in zip(regions, pair_object.regions):
                    assert r1.chromosome == r2.chromosome
                    assert r1.start == r2.start
                    assert r1.end == r2.end
            except AssertionError:
                raise ValueError("Regions in pair objects are not identical, cannot perform merge!")

        merged_pairs.add_regions(regions)

        for pair_object in pairs:
            merged_pairs.add_edges(pair_object.edges(lazy=True))

        return merged_pairs

    def regions_and_edges(self, key, *args, **kwargs):
        row_regions, col_regions = self._key_to_regions(key)
        if key is None:
            edges = self._edges_iter(*args, **kwargs)
        else:
            edges = self._edges_subset(key, row_regions, col_regions, *args, **kwargs)

        return row_regions, col_regions, edges

    def mappable(self):
        """
        Get the mappability vector of this matrix.
        """
        return np.array([True if r.valid else False for r in self.regions])


class RegionMatrixContainer(RegionPairsContainer, RegionBasedWithBins):
    def __init__(self):
        RegionPairsContainer.__init__(self)
        self._default_value = 0.0
        self._default_score_field = 'weight'

    def regions_and_matrix_entries(self, key, norm=True, oe=False, oe_per_chromosome=True,
                                   bias_field='bias', score_field=None,
                                   *args, **kwargs):
        row_regions, col_regions = self._key_to_regions(key)
        if isinstance(row_regions, GenomicRegion):
            row_regions = [row_regions]
        else:
            row_regions = list(row_regions)

        if isinstance(col_regions, GenomicRegion):
            col_regions = [col_regions]
        else:
            col_regions = list(col_regions)

        try:
            row_offset = row_regions[0].ix
            col_offset = col_regions[0].ix
        except IndexError:
            return row_regions, col_regions, []

        basic_iter = self._matrix_entries(key, row_regions, col_regions,
                                          score_field=score_field,
                                          *args, **kwargs)

        if norm:
            biases = dict()
            for regions in (row_regions, col_regions):
                for region in regions:
                    biases[region.ix] = getattr(region, bias_field, 1.0)
        else:
            biases = defaultdict(lambda: 1)

        def offset_iter(edge_iter):
            for source, sink, weight in edge_iter:
                i = source - row_offset
                j = sink - col_offset
                if i >= 0 and j >= 0:
                    yield source, sink, i, j, weight * biases[source] * biases[sink]

                l = source - col_offset
                k = sink - row_offset
                if (i, j) != (k, l) and k >= 0 and l >= 0:
                    yield source, sink, k, l, weight * biases[source] * biases[sink]

        if not oe:
            entry_iter = ((i, j, weight)
                          for source, sink, i, j, weight in offset_iter(basic_iter))

        else:
            intra_expected, chromosome_intra_expected, inter_expected = self.expected_values(norm=norm)

            entry_iter = ((i, j,
                           weight / chromosome_intra_expected[row_regions[i].chromosome][abs(source - sink)]
                           if row_regions[i].chromosome == col_regions[j].chromosome
                           else weight / inter_expected)
                          for source, sink, i, j, weight in offset_iter(basic_iter))

        return row_regions, col_regions, entry_iter

    def _matrix_entries(self, key, row_regions, col_regions,
                        score_field=None, *args, **kwargs):
        if score_field is None:
            score_field = self._default_score_field

        for edge in self._edges_subset(key, row_regions, col_regions, *args, **kwargs):
            yield (edge.source, edge.sink,
                   getattr(edge, score_field, self._default_value))

    def matrix(self, key=None, norm=True, oe=False,
               score_field=None, bias_field='bias',
               default_value=None, mask=True,
               _mappable=None):

        if score_field is None:
            score_field = self._default_score_field

        if default_value is None:
            default_value = self._default_value

        row_regions, col_regions, matrix_entries = self.regions_and_matrix_entries(key, norm=norm, oe=oe,
                                                                                   score_field=score_field,
                                                                                   bias_field=bias_field)

        m = np.full((len(row_regions), len(col_regions)), default_value)

        for source, sink, weight in matrix_entries:
            ir = source
            jr = sink
            if 0 <= ir < m.shape[0] and 0 <= jr < m.shape[1]:
                m[ir, jr] = weight

        if isinstance(key, tuple) and len(key) == 2:
            if isinstance(key[0], int) and isinstance(key[1], int):
                return m[0, 0]
            elif isinstance(key[0], int):
                m = m[0, :]
            elif isinstance(key[1], int):
                m = m[:, 0]

        return RegionMatrix(m, row_regions=row_regions, col_regions=col_regions, mask=mask)

    def __getitem__(self, item):
        return self.matrix(item)

    def possible_contacts(self):
        logger.debug("Calculating possible counts")
        regions = list(self.regions)
        chromosomes = self.chromosomes()

        cb = self.chromosome_bins
        chromosome_max_distance = defaultdict(int)
        max_distance = 0
        chromosome_subtractions = dict()
        for chromosome in chromosomes:
            start, stop = cb[chromosome]
            max_distance = max(max_distance, stop - start)
            chromosome_max_distance[chromosome] = max(chromosome_max_distance[chromosome], stop - start)
            chromosome_subtractions[chromosome] = np.zeros(stop - start,
                                                           dtype='int32')

        chromosome_mappable = defaultdict(int)
        chromosome_unmappable = defaultdict(set)
        for i, mappable in enumerate(self.mappable()):
            chromosome = regions[i].chromosome
            if not mappable:  # unmappable
                s = chromosome_subtractions[chromosome]
                o = cb[chromosome][0]
                ix = i - o
                # horizontal
                s[0: len(s) - ix] += 1
                # vertical
                for j in range(1, ix + 1):
                    if ix - j not in chromosome_unmappable[chromosome]:
                        s[j] += 1
                chromosome_unmappable[chromosome].add(ix)
            else:
                chromosome_mappable[chromosome] += 1

        inter_total = 0
        intra_total = [0] * max_distance
        chromosome_intra_total = dict()
        for chromosome, d in chromosome_max_distance.items():
            chromosome_intra_total[chromosome] = [0] * d

        for i, chromosome in enumerate(chromosomes):
            start, stop = cb[chromosome]
            count = stop - start

            # intra-chromosomal
            s = chromosome_subtractions[chromosomes[i]]
            for distance in range(0, count):
                intra_total[distance] += count - distance - s[distance]
                chromosome_intra_total[chromosome][distance] += count - distance - s[distance]

            # inter-chromosomal
            for j in range(i + 1, len(chromosomes)):
                count_mappable = chromosome_mappable[chromosomes[i]]
                count2_mappable = chromosome_mappable[chromosomes[j]]
                inter_total += count_mappable * count2_mappable

        return intra_total, chromosome_intra_total, inter_total

    def expected_values(self, selected_chromosome=None, norm=True):
        # get all the bins of the different chromosomes
        chromosome_bins = self.chromosome_bins
        chromosome_dict = defaultdict(list)

        chromosome_max_distance = defaultdict(int)
        max_distance = 0
        for chromosome, (start, stop) in chromosome_bins.items():
            max_distance = max(max_distance, stop - start)
            chromosome_max_distance[chromosome] = max(chromosome_max_distance[chromosome], stop - start)

            for i in range(start, stop):
                chromosome_dict[i] = chromosome

        chromosome_intra_sums = dict()
        chromosome_intra_expected = dict()
        for chromosome, d in chromosome_max_distance.items():
            chromosome_intra_sums[chromosome] = [0.0] * d
            chromosome_intra_expected[chromosome] = [0.0] * d

        # get the sums of edges at any given distance
        marginals = [0.0] * len(self.regions)
        inter_sums = 0.0
        intra_sums = [0.0] * max_distance
        with RareUpdateProgressBar(max_value=len(self.edges), prefix='Expected') as pb:
            for i, edge in enumerate(self.edges(lazy=True, norm=norm)):
                source, sink = edge.source, edge.sink
                weight = getattr(edge, self._default_score_field)

                source_chromosome = chromosome_dict[source]
                sink_chromosome = chromosome_dict[sink]

                marginals[source] += weight
                marginals[sink] += weight

                if sink_chromosome != source_chromosome:
                    inter_sums += weight
                else:
                    distance = sink - source
                    intra_sums[distance] += weight
                    chromosome_intra_sums[source_chromosome][distance] += weight
                pb.update(i)

        intra_total, chromosome_intra_total, inter_total = self.possible_contacts()

        # expected values
        inter_expected = 0 if inter_total == 0 else inter_sums/inter_total

        intra_expected = [0.0] * max_distance
        bin_size = self.bin_size
        distances = []
        for d in range(max_distance):
            distances.append(bin_size * d)

            # whole genome
            count = intra_total[d]
            if count > 0:
                intra_expected[d] = intra_sums[d] / count

        # chromosomes
        for chromosome in chromosome_intra_expected:
            for d in range(chromosome_max_distance[chromosome]):
                chromosome_count = chromosome_intra_total[chromosome][d]
                if chromosome_count > 0:
                    chromosome_intra_expected[chromosome][d] = chromosome_intra_sums[chromosome][d] / chromosome_count

        if selected_chromosome is not None:
            return chromosome_intra_expected[selected_chromosome]

        return intra_expected, chromosome_intra_expected, inter_expected

    def marginals(self, weight_column=None, norm=True):
        """
        Get the marginals vector of this Hic matrix.
        """
        if weight_column is None:
            weight_column = self._default_score_field

        # prepare marginals dict
        marginals = np.zeros(len(self.regions), float)

        logger.debug("Calculating marginals...")
        with RareUpdateProgressBar(max_value=len(self.edges), silent=config.hide_progressbars) as pb:
            for i, edge in enumerate(self.edges(lazy=True, norm=norm)):
                marginals[edge.source] += getattr(edge, weight_column)
                if edge.source != edge.sink:
                    marginals[edge.sink] += getattr(edge, weight_column)
                pb.update(i)

        return marginals


class RegionPairsTable(RegionPairsContainer, Maskable, RegionsTable):

    _classid = 'REGIONPAIRSTABLE'

    def __init__(self, file_name=None, mode='a', tmpdir=None,
                 additional_region_fields=None, additional_edge_fields=None,
                 partitioning_strategy='chromosome',
                 _table_name_regions='regions', _table_name_edges='edges',
                 _edge_buffer_size=1000000):

        """
        Initialize a :class:`~RegionPairsTable` object.

        :param file_name: Path to a save file
        :param mode: File mode to open underlying file
        :param additional_region_fields: Additional fields (in PyTables notation) associated with
                                         edge data, e.g. {'weight': tables.Float32Col()}
        :param _table_name_regions: (Internal) name of the HDF5 node for regions
        :param _table_name_edges: (Internal) name of the HDF5 node for edges
        :param _edge_buffer_size: (Internal) size of edge / contact buffer
        """

        # private variables
        self._edges_dirty = False
        self._mappability_dirty = False
        self._partitioning_strategy = partitioning_strategy

        file_exists = False
        if file_name is not None:
            file_name = os.path.expanduser(file_name)
            if os.path.exists(file_name):
                file_exists = True

        # initialise inherited objects
        RegionPairsContainer.__init__(self)

        if additional_region_fields is None:
            additional_region_fields = {}
        additional_region_fields['valid'] = tables.BoolCol(dflt=True)

        RegionsTable.__init__(self, file_name=file_name, _table_name_regions=_table_name_regions,
                              mode=mode, tmpdir=tmpdir, additional_fields=additional_region_fields)
        Maskable.__init__(self, self.file)

        self._edge_table_dict = dict()
        if file_exists and mode != 'w':
            # retrieve edge tables and partitions
            self._edges = self.file.get_node('/', _table_name_edges)
            self._partition_breaks = getattr(self.meta, 'partition_breaks', None)
            if self._partition_breaks is None:
                self._update_partitions()
        else:
            self._edges = self.file.create_group('/', _table_name_edges)

            basic_fields = {
                'source': tables.Int32Col(pos=0),
                'sink': tables.Int32Col(pos=1),
            }
            if additional_edge_fields is not None:
                if not isinstance(additional_edge_fields, dict) and issubclass(additional_edge_fields,
                                                                               tables.IsDescription):
                    # IsDescription subclass case
                    additional_edge_fields = additional_edge_fields.columns

                for field in additional_edge_fields.values():
                    if field._v_pos is None:
                        field._v_pos = len(additional_edge_fields)

                current = len(basic_fields)
                for key, value in sorted(additional_edge_fields.items(), key=lambda x: x[1]._v_pos):
                    if key not in basic_fields:
                        if value._v_pos is not None:
                            value._v_pos = current
                            current += 1
                        basic_fields[key] = value

            self._partition_breaks = None
            self._update_partitions()

            self._create_edge_table(0, 0, fields=basic_fields)

        # update edge table dict
        self._edge_table_dict = dict()
        for edge_table in self._edges:
            self._edge_table_dict[(edge_table.attrs['source_partition'],
                                   edge_table.attrs['sink_partition'])] = edge_table

        # update field names
        self._source_field_ix = 0
        self._sink_field_ix = 0
        self.field_names = []
        self._field_names_dict = dict()
        self._edge_field_defaults = dict()
        self._update_field_names()

        # set up edge buffer
        self._edge_buffer = defaultdict(list)
        self._edge_buffer_size = _edge_buffer_size

    def _flush_table_edge_buffer(self):
        for (source_partition, sink_partition), records in self._edge_buffer.items():
            if not (source_partition, sink_partition) in self._edge_table_dict:
                self._create_edge_table(source_partition, sink_partition)
            table = self._edge_table_dict[(source_partition, sink_partition)]
            table.append(records)
        self._edge_buffer = defaultdict(list)

    def _flush_regions(self):
        if self._regions_dirty:
            RegionsTable._flush_regions(self)
            self._update_partitions()

    def _flush_edges(self, silent=config.hide_progressbars):
        if self._edges_dirty:
            if len(self._edge_buffer) > 0:
                logger.debug("Adding buffered edges...")
                self._flush_table_edge_buffer()

            with RareUpdateProgressBar(max_value=sum(1 for _ in self._edges), silent=silent) as pb:
                for i, edge_table in enumerate(self._edges):
                    edge_table.flush(update_index=False, log_progress=False)
                    pb.update(i)

            self._enable_edge_indexes()
            for i, edge_table in enumerate(self._edges):
                edge_table.flush(update_index=True, log_progress=False)
                pb.update(i)
            self._edges_dirty = False

            self._update_mappability()

    def flush(self, silent=config.hide_progressbars):
        """
        Write data to file and flush buffers.

        :param silent: do not print flush progress
        """
        self._flush_regions()
        self._flush_edges(silent=silent)

    def _disable_edge_indexes(self):
        for edge_table in self._edge_table_dict.values():
            edge_table.cols.source.remove_index()
            edge_table.cols.sink.remove_index()
            edge_table.disable_mask_index()

    def _enable_edge_indexes(self):
        for edge_table in self._edge_table_dict.values():
            create_col_index(edge_table.cols.source)
            create_col_index(edge_table.cols.sink)
            edge_table.enable_mask_index()

    def _update_partitions(self):
        partition_breaks = []

        if self._partitioning_strategy == 'chromosome':
            previous_chromosome = None
            for i, region in enumerate(self.regions(lazy=True)):
                if region.chromosome != previous_chromosome and previous_chromosome is not None:
                    partition_breaks.append(i)
                previous_chromosome = region.chromosome
        elif isinstance(self._partitioning_strategy, int):
            n_regions = len(self.regions)
            for i in range(self._partitioning_strategy, n_regions, self._partitioning_strategy):
                partition_breaks.append(i)

        self._partition_breaks = partition_breaks
        try:
            self.meta['partition_breaks'] = partition_breaks
        except tables.FileModeError:
            pass

    def _update_field_names(self):
        """
        Set internal object variables related to edge table field names.
        """
        edge_table = self._edge_table_dict[(0, 0)]

        # update field names
        self._source_field_ix = 0
        self._sink_field_ix = 0
        self.field_names = []
        self._field_names_dict = dict()
        self._edge_field_defaults = dict()
        for i, name in enumerate(edge_table.colnames):
            if not name.startswith("_"):
                self.field_names.append(name)
            if name == 'source':
                self._source_field_ix = i
            if name == 'sink':
                self._sink_field_ix = i
            self._field_names_dict[name] = i
            self._edge_field_defaults[name] = edge_table.coldescrs[name].dflt

    def _create_edge_table(self, source_partition, sink_partition, fields=None):
        """
        Create and register an edge table for a partition combination.
        """
        if fields is None:
            fields = self._edge_table_dict[(0, 0)].coldescrs

        if (source_partition, sink_partition) in self._edge_table_dict:
            return self._edge_table_dict[(source_partition, sink_partition)]

        edge_table = MaskedTable(self._edges,
                                 'chrpair_' + str(source_partition) + '_' + str(sink_partition),
                                 fields, ignore_reserved_fields=True)
        edge_table.attrs['source_partition'] = source_partition
        edge_table.attrs['sink_partition'] = sink_partition

        # index
        create_col_index(edge_table.cols.source)
        create_col_index(edge_table.cols.sink)

        self._edge_table_dict[(source_partition, sink_partition)] = edge_table
        return edge_table

    def _get_edge_table_tuple(self, source, sink):
        if source > sink:
            source, sink = sink, source

        source_partition = self._get_partition_ix(source)
        sink_partition = self._get_partition_ix(sink)

        return source_partition, sink_partition

    def _add_edge(self, edge, row=None, replace=False):
        """
        Add an edge to an internal edge table.
        """
        if not self._edges_dirty:
            self._edges_dirty = True
            self._disable_edge_indexes()

        source, sink = edge.source, edge.sink
        if source > sink:
            source, sink = sink, source

        if row is None:
            record = [None] * len(self._field_names_dict)
            for name, ix in self._field_names_dict.items():
                try:
                    record[ix] = getattr(edge, name)
                except AttributeError:
                    record[ix] = self._edge_field_defaults[name]
            record[self._field_names_dict['source']] = source
            record[self._field_names_dict['sink']] = sink

            self._add_edge_from_tuple(record)
        else:
            row['source'] = source
            row['sink'] = sink
            for name in self.field_names:
                if not name == 'source' and not name == 'sink':
                    try:
                        value = getattr(edge, name)
                        if replace:
                            row[name] = value
                        else:
                            row[name] += value
                    except AttributeError:
                        pass
            row.update()

    def _add_edge_from_tuple(self, edge):
        if not self._edges_dirty:
            self._edges_dirty = True
            self._disable_edge_indexes()

        source = edge[self._source_field_ix]
        sink = edge[self._sink_field_ix]
        if source > sink:
            source, sink = sink, source
        source_partition, sink_partition = self._get_edge_table_tuple(source, sink)

        self._edge_buffer[(source_partition, sink_partition)].append(tuple(edge))
        if sum(len(records) for records in self._edge_buffer.values()) > self._edge_buffer_size:
            self._flush_table_edge_buffer()

    def add_edges(self, edges, *args, **kwargs):
        if self._regions_dirty:
            self._flush_regions()
        RegionPairsContainer.add_edges(self, edges, *args, **kwargs)
        self._flush_edges()

    def _get_partition_ix(self, region_ix):
        """
        Bisect the partition table to get the partition index for a region index.
        """
        return bisect_right(self._partition_breaks, region_ix)

    def _is_partition_covered(self, partition_ix, region_ix_start, region_ix_end):
        try:
            partition_end = self._partition_breaks[partition_ix]
        except IndexError:
            partition_end = len(self.regions) - 1

        if partition_ix > 0:
            partition_start = self._partition_breaks[partition_ix - 1]
        else:
            partition_start = 0

        if region_ix_start <= partition_start and region_ix_end >= partition_end:
            return True
        return False

    def _edge_subset_rows(self, key):
        row_regions, col_regions = self._key_to_regions(key, lazy=True)

        return self._edge_subset_rows_from_regions(
            row_regions, col_regions
        )

    def _edge_subset_rows_from_regions(self, row_regions, col_regions):
        row_start, row_end = self._min_max_region_ix(row_regions)
        col_start, col_end = self._min_max_region_ix(col_regions)

        row_partition_start = self._get_partition_ix(row_start)
        row_partition_end = self._get_partition_ix(row_end)
        col_partition_start = self._get_partition_ix(col_start)
        col_partition_end = self._get_partition_ix(col_end)

        partition_extracted = set()
        for a in range(row_partition_start, row_partition_end + 1):
            for b in range(col_partition_start, col_partition_end + 1):
                if b < a:
                    i, j = b, a
                else:
                    i, j = a, b

                if (i, j) in partition_extracted:
                    continue
                else:
                    partition_extracted.add((i, j))

                row_covered = self._is_partition_covered(a, row_start, row_end)
                col_covered = self._is_partition_covered(b, col_start, col_end)

                if (i, j) not in self._edge_table_dict:
                    continue

                edge_table = self._edge_table_dict[(i, j)]

                # if we need to get all regions in a table, return the whole thing
                if row_covered and col_covered:
                    for row in edge_table:
                        yield row

                # otherwise only return the subset defined by the respective indices
                else:
                    condition = "(%d < source) & (source < %d) & (% d < sink) & (sink < %d)"
                    condition1 = condition % (row_start - 1, row_end + 1, col_start - 1, col_end + 1)
                    condition2 = condition % (col_start - 1, col_end + 1, row_start - 1, row_end + 1)

                    if row_start > col_start:
                        condition1, condition2 = condition2, condition1

                    overlap = range_overlap(row_start, row_end, col_start, col_end)

                    for edge_row in edge_table.where(condition1):
                        yield edge_row

                    for edge_row in edge_table.where(condition2):
                        if overlap is not None:
                            if (overlap[0] <= edge_row['source'] <= overlap[1]) and (
                                    overlap[0] <= edge_row['sink'] <= overlap[1]):
                                continue

                        yield edge_row

    def _matrix_entries(self, key, row_regions, col_regions,
                        score_field=None, *args, **kwargs):
        if score_field is None:
            score_field = self._default_score_field

        for row in self._edge_subset_rows_from_regions(row_regions, col_regions):
            yield (row['source'], row['sink'], row[score_field])

    def _row_to_edge(self, row, lazy=False, **kwargs):
        if not lazy:
            source = row["source"]
            sink = row["sink"]
            d = dict()
            for field in self.field_names:
                if field != 'source' and field != 'sink':
                    value = row[field]
                    value = value.decode() if isinstance(value, bytes) else value
                    d[field] = value

            source_node = self.regions[source]
            sink_node = self.regions[sink]
            return Edge(source_node, sink_node, **d)
        else:
            return LazyEdge(row, self._regions)

    def _edges_subset(self, key=None, row_regions=None, col_regions=None,
                      *args, **kwargs):
        for row in self._edge_subset_rows_from_regions(row_regions, col_regions):
            yield self._row_to_edge(row, *args, **kwargs)

    def _edges_iter(self, *args, **kwargs):
        for i in range(0, len(self._partition_breaks) + 1):
            for j in range(i, len(self._partition_breaks) + 1):
                if (i, j) in self._edge_table_dict:
                    for row in self._edge_table_dict[(i, j)]:
                        yield self._row_to_edge(row, *args, **kwargs)

    def _edges_length(self):
        s = 0
        for edge_table in self._edge_table_dict.values():
            s += len(edge_table)
        return s

    def _edges_getitem(self, item, *args, **kwargs):
        result = list(self.edges(item))
        if isinstance(item, tuple) and len(item) == 2 and isinstance(item[0], int) and isinstance(item[1], int):
            return result[0]
        if isinstance(item, int):
            return result[0]
        return result

    def _update_mappability(self):
        logger.info("Updating region mappability")
        mappable = [False] * len(self.regions)
        with RareUpdateProgressBar(max_value=len(self.edges), prefix='Mappability',
                                   silent=config.hide_progressbars) as pb:
            for i, edge in enumerate(self.edges(lazy=True, norm=False)):
                try:
                    weight = getattr(edge, self._default_score_field, self._default_value)
                except TypeError:
                    weight = self._default_value

                if weight > 0:
                    mappable[edge.source] = True
                    mappable[edge.sink] = True
                pb.update(i)
        self.region_data('valid', mappable)

    def filter(self, edge_filter, queue=False, log_progress=not config.hide_progressbars):
        """
        Filter edges in this object by using a
        :class:`~HicEdgeFilter`.

        :param edge_filter: Class implementing :class:`~HicEdgeFilter`.
                            Must override valid_edge method, ideally sets mask parameter
                            during initialization.
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        :param log_progress: If true, process iterating through all edges
                             will be continuously reported.
        """
        total = 0
        filtered = 0
        if not queue:
            with RareUpdateProgressBar(max_value=len(self._edge_table_dict),
                                       silent=not log_progress) as pb:
                for i, edge_table in enumerate(self._edge_table_dict.values()):
                    stats = edge_table.filter(edge_filter, _logging=False)
                    for key, value in stats.items():
                        if key != 0:
                            filtered += stats[key]
                        total += stats[key]
                    pb.update(i)
            if log_progress:
                logger.info("Total: {}. Filtered: {}".format(total, filtered))
        else:
            for edge_table in self._edge_table_dict.values():
                edge_table.queue_filter(edge_filter)

        self._update_mappability()

    def run_queued_filters(self, log_progress=not config.hide_progressbars):
        """
        Run queued filters.

        :param log_progress: If true, process iterating through all edges
                             will be continuously reported.
        """
        total = 0
        filtered = 0
        with RareUpdateProgressBar(max_value=len(self._edge_table_dict),
                                   silent=not log_progress) as pb:
            for i, edge_table in enumerate(self._edge_table_dict.values()):
                stats = edge_table.run_queued_filters(_logging=False)
                for key, value in stats.items():
                    if key != 0:
                        filtered += stats[key]
                    total += stats[key]
                pb.update(i)
        if log_progress:
            logger.info("Total: {}. Filtered: {}".format(total, filtered))

        self._update_mappability()

    def sample(self, n, exact=False, file_name=None):
        if isinstance(n, RegionPairsContainer):
            n = len(n.edges)

        region_pairs = []
        if not exact:
            weights = []
            logger.info("Using sampling with replacement")
            for edge in self.edges(lazy=True):
                region_pairs.append((edge.source, edge.sink))
                weights.append(edge.weight)
            s = sum(weights)
            p = [w / s for w in weights]
        else:
            p = None
            logger.info("Using sampling without replacement")
            for edge in self.edges(lazy=True):
                for i in range(int(edge.weight)):
                    region_pairs.append((edge.source, edge.sink))

        new_pairs = self.__class__(file_name=file_name, mode='w')
        new_pairs.add_regions(self.regions, preserve_attributes=False)
        new_edges = defaultdict(int)
        for new_pair_ix in np.random.choice(len(region_pairs), size=n, replace=not exact, p=p):
            new_edges[region_pairs[new_pair_ix]] += 1
        new_edges = [[source, sink, weight] for (source, sink), weight in new_edges.items()]
        new_pairs.add_edges(new_edges)

        return new_pairs

    def subset(self, *regions, **kwargs):
        """
        Subset a Hic object by specifying one or more subset regions.

        :param regions: string or GenomicRegion object(s)
        :param kwargs: Supports
                       file_name: destination file name of subset Hic object;
                       tmpdir: if True works in tmp until object is closed
        :return: Hic
        """
        file_name = kwargs.get("file_name", None)
        tmpdir = kwargs.get('tmpdir', None)

        new_pairs = self.__class__(file_name=file_name, mode='w', tmpdir=tmpdir)

        ix_converter = {}
        ix = 0
        new_regions = []
        for region_string in regions:
            for region in self.regions(region_string):
                ix_converter[region.ix] = ix
                ix += 1
                new_regions.append(region)
        new_pairs.add_regions(new_regions, preserve_attributes=False)

        for i, region_string1 in enumerate(regions):
            for j in range(i, len(regions)):
                region_string2 = regions[j]
                for edge in self.edges((region_string1, region_string2), lazy=True):
                    source = ix_converter[edge.source]
                    sink = ix_converter[edge.sink]
                    new_pairs.add_edge([source, sink, edge.weight])
        new_pairs.flush()

        return new_pairs


class RegionMatrixTable(RegionMatrixContainer, RegionPairsTable):

    _classid = 'REGIONMATRIXTABLE'

    def __init__(self, file_name=None, mode='a', tmpdir=None,
                 partitioning_strategy='chromosome',
                 additional_region_fields=None, additional_edge_fields=None,
                 default_score_field='weight', default_value=0.0,
                 _table_name_regions='regions', _table_name_edges='edges',
                 _table_name_expected_values='expected_values',
                 _edge_buffer_size=1000000):

        self._default_score_field = default_score_field
        self._default_value = default_value

        if additional_edge_fields is None:
            additional_edge_fields = {}
        additional_edge_fields['weight'] = tables.Float64Col(pos=0)

        if additional_region_fields is None:
            additional_region_fields = {}
        additional_region_fields['valid'] = tables.BoolCol(dflt=True)
        additional_region_fields['bias'] = tables.Float64Col(dflt=1.0)

        RegionPairsTable.__init__(self,
                                  file_name=file_name, mode=mode, tmpdir=tmpdir,
                                  additional_region_fields=additional_region_fields,
                                  additional_edge_fields=additional_edge_fields,
                                  partitioning_strategy=partitioning_strategy,
                                  _table_name_regions=_table_name_regions,
                                  _table_name_edges=_table_name_edges,
                                  _edge_buffer_size=_edge_buffer_size)
        RegionMatrixContainer.__init__(self)

        file_exists = False
        if file_name is not None and os.path.exists(os.path.expanduser(file_name)):
            file_exists = True

        # create expected value group
        if file_exists and mode != 'w':
            try:
                self._expected_value_group = self.file.get_node('/', _table_name_expected_values)
            except tables.NoSuchNodeError:
                if mode in ['a', 'r+']:
                    self._expected_value_group = self.file.create_group('/', _table_name_expected_values)
                else:
                    self._expected_value_group = None
        else:
            self._expected_value_group = self.file.create_group('/', _table_name_expected_values)

    def _remove_expected_values(self):
        if self._expected_value_group is not None:
            try:
                self.file.remove_node(self._expected_value_group, 'corrected', recursive=True)
            except tables.NoSuchNodeError:
                pass

            try:
                self.file.remove_node(self._expected_value_group, 'uncorrected', recursive=True)
            except tables.NoSuchNodeError:
                pass

    def _flush_edges(self, silent=config.hide_progressbars):
        if self._edges_dirty:
           self._remove_expected_values()

        RegionPairsTable._flush_edges(self, silent=silent)

    def set_biases(self, biases):
        self.region_data('bias', biases)
        self._remove_expected_values()

    def expected_values(self, selected_chromosome=None, norm=True):
        group_name = 'corrected' if norm else 'uncorrected'

        if self._expected_value_group is not None:
            try:
                group = self.file.get_node(self._expected_value_group, group_name)
                if selected_chromosome is not None:
                    return self.file.get_node(group, '_' + selected_chromosome)[:]
                else:
                    intra_expected = None
                    inter_expected = None
                    chromosome_intra_expected = {}
                    for node in self.file.walk_nodes(group):
                        if isinstance(node, tables.Group):
                            continue
                        if node.name == '__intra__':
                            intra_expected = node[:]
                        elif node.name == '__inter__':
                            inter_expected = node[0]
                        else:
                            if node.name.startswith('_'):
                                chromosome = node.name[1:]
                                chromosome_intra_expected[chromosome] = node[:]

                    if intra_expected is not None and inter_expected is not None \
                            and len(chromosome_intra_expected) > 0:
                        return intra_expected, chromosome_intra_expected, inter_expected
            except tables.NoSuchNodeError:
                pass

        intra_expected, chromosome_intra_expected, inter_expected = RegionMatrixContainer.expected_values(self, norm=norm)

        # try saving to object
        try:
            try:
                group = self.file.get_node(self._expected_value_group, group_name)
            except tables.NoSuchNodeError:
                group = self.file.create_group(self._expected_value_group, group_name)

            self.file.create_array(group, '__intra__',
                                   np.array(intra_expected), "Intra-chromosomal expected values")
            self.file.create_array(group, '__inter__',
                                   np.array([inter_expected]), "Inter-chromosomal expected value")
            for chromosome, values in chromosome_intra_expected.items():
                self.file.create_array(group, '_' + chromosome,
                                       np.array(values), "Intra-chromosomal expected "
                                                         "value {}".format(chromosome))
        except tables.FileModeError:
            warnings.warn("Matrix file opened in read-only mode, "
                          "cannot save expected values to object. "
                          "Use mode 'a' to add expected values to "
                          "an existing object!")

        if selected_chromosome is not None:
            return chromosome_intra_expected[selected_chromosome]

        return intra_expected, chromosome_intra_expected, inter_expected

    @classmethod
    def merge(cls, matrices, *args, **kwargs):
        if isinstance(matrices, RegionMatrixContainer):
            matrices = [matrices]

        if 'mode' not in kwargs:
            kwargs['mode'] = 'w'
        merged_matrix = cls(*args, **kwargs)

        matrices = [matrix_object for matrix_object in matrices]  # ensure list

        regions = list(matrices[0].regions)
        for matrix_object in matrices[1:]:
            try:
                for r1, r2 in zip(regions, matrix_object.regions):
                    assert r1.chromosome == r2.chromosome
                    assert r1.start == r2.start
                    assert r1.end == r2.end
            except AssertionError:
                raise ValueError("Regions in matrix objects are not identical, cannot perform merge!")

        merged_matrix.add_regions(regions, preserve_attributes=False)

        default_field = merged_matrix._default_score_field

        chromosomes = merged_matrix.chromosomes()
        for i in range(len(chromosomes)):
            chromosome1 = chromosomes[i]
            for j in range(i, len(chromosomes)):
                chromosome2 = chromosomes[j]

                edges = defaultdict(int)
                for matrix_object in matrices:
                    for edge in matrix_object.edges((chromosome1, chromosome2), lazy=True):
                        edges[edge.source, edge.sink] += getattr(edge, default_field)

                for (source, sink), weight in edges.items():
                    e = Edge(source=source, sink=sink)
                    setattr(e, default_field, weight)
                    merged_matrix.add_edge(e)

        merged_matrix.flush()
        return merged_matrix


class RegionMatrix(np.ma.MaskedArray):
    def __new__(cls, input_matrix, col_regions=None, row_regions=None,
                mask=True, *args, **kwargs):
        obj = np.asarray(input_matrix).view(cls, *args, **kwargs)
        obj._row_region_trees = None
        obj._col_region_trees = None
        obj.col_regions = None
        obj.row_regions = None
        obj.set_col_regions(col_regions)
        obj.set_row_regions(row_regions)
        obj._do_mask = mask
        if mask:
            obj._apply_mask()
        return obj

    def _apply_mask(self):
        if self.row_regions is not None and self.col_regions is not None:
            mask = np.zeros(self.shape, dtype=bool)

            try:
                row_offset = self.row_regions[0].ix
                col_offset = self.col_regions[0].ix
            except IndexError:
                return

            for regions in (self.row_regions, self.col_regions):
                for region in regions:
                    valid = getattr(region, 'valid', True)
                    if not valid:
                        row_ix = region.ix - row_offset
                        if 0 <= row_ix < self.shape[0]:
                            mask[row_ix] = True

                        col_ix = region.ix - col_offset
                        if 0 <= col_ix <= self.shape[1]:
                            mask[:, col_ix] = True

            self.mask = mask

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

    def set_col_regions(self, regions):
        self.col_regions = regions
        if regions is not None:
            self._col_region_trees = self._interval_tree_regions(regions)
        else:
            self._col_region_trees = None

    def __array_finalize__(self, obj):
        if isinstance(self, np.ma.core.MaskedArray):
            np.ma.MaskedArray.__array_finalize__(self, obj)

        if obj is None:
            return

        self.set_row_regions(getattr(obj, 'row_regions', None))
        self.set_col_regions(getattr(obj, 'col_regions', None))
        mask = getattr(obj, '_do_mask', True)
        self._do_mask = mask
        if mask:
            self._apply_mask()

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
            all_cols = slice(0, len(self.col_regions), 1)
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
                col_key = self._convert_key(
                    index[1],
                    self._col_region_trees if hasattr(self, '_col_region_trees') else None
                )
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
            col_regions = self.col_regions[col_key]
        except TypeError:
            col_regions = None

        if isinstance(row_regions, GenomicRegion):
            out.row_regions = [row_regions]
        else:
            out.row_regions = row_regions

        if isinstance(col_regions, GenomicRegion):
            out.col_regions = [col_regions]
        else:
            out.col_regions = col_regions

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
