import logging

import numpy as np
from genomic_regions import RegionBased, GenomicRegion

from .config import config
from .regions import LazyGenomicRegion
from .tools.general import RareUpdateProgressBar

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
            base_info += "\n\t{}: {}".format(field, str(getattr(self, field)))
        return base_info + "\n"


class LazyEdge(Edge):
    def __init__(self, row, nodes_table=None, auto_update=True):
        self.reserved = {'_row', '_nodes_table', 'auto_update', '_source_node', '_sink_node'}
        self._row = row
        self._nodes_table = nodes_table
        self.auto_update = auto_update
        self._source_node = None
        self._sink_node = None

    def _set_item(self, item, value):
        self._row[item] = value
        if self.auto_update:
            self.update()

    def __getattr__(self, item):
        if item == 'reserved' or item in self.reserved:
            return object.__getattribute__(self, item)
        try:
            value = self._row[item]
            value = value.decode() if isinstance(value, bytes) else value
            return value
        except KeyError:
            raise AttributeError("Attribute not supported (%s)" % str(item))

    def __setattr__(self, key, value):
        if key == 'reserved' or key in self.reserved:
            super(LazyEdge, self).__setattr__(key, value)
        else:
            self._row[key] = value
            if self.auto_update:
                self.update()

    def update(self):
        self._row.update()

    @property
    def source(self):
        return self._row['source']

    @property
    def sink(self):
        return self._row['sink']

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

    def _add_edge(self, edge, *args, **kwargs):
        raise NotImplementedError("Subclass must override this function")

    def _edge_from_object(self, edge):
        return edge

    def _edge_from_dict(self, edge):
        source, sink = edge['source'], edge['sink']

        attributes = dict()
        for name, value in edge.items():
            if not name == 'source' and not name == 'sink':
                attributes[name] = value

        return Edge(source, sink, **attributes)

    def _edge_from_list(self, edge):
        source, sink = edge[0], edge[1]
        try:
            weight = edge[2]
        except IndexError:
            weight = 1

        return Edge(source, sink, weight=weight)

    def _edges_iter(self, *args, **kwargs):
        raise NotImplementedError("Subclass must implement _edges_iter "
                                  "to enable iterating over edges!")

    def _edges_subset(self, key=None, *args, **kwargs):
        raise NotImplementedError("Subclass must implement _edges_subset "
                                  "to enable iterating over edge subsets!")

    def _edges_length(self):
        return sum(1 for _ in self.edges)

    def _edges_getitem(self, item, *args, **kwargs):
        raise NotImplementedError("Subclass must implement _edges_getitem "
                                  "to enable getting specific edges!")

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
                raise ValueError("Node index exceeds number of nodes in object")

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
                return self._regions_pairs._edges_iter()

            def __call__(self, key=None, *args, **kwargs):
                if key is None:
                    return self._regions_pairs._edges_iter(*args, **kwargs)
                else:
                    return self._regions_pairs.edges_subset(key, *args, **kwargs)

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
        return self._edges_subset(key, *args, **kwargs)

    def mappable(self):
        """
        Get the mappability vector of this matrix.
        """
        logger.debug("Calculating mappability...")

        mappable = np.zeros(len(self.regions), dtype=bool)
        with RareUpdateProgressBar(max_value=len(self.edges), silent=config.hide_progressbars) as pb:
            for i, edge in enumerate(self.edges(lazy=True)):
                mappable[edge.source] = True
                mappable[edge.sink] = True
                pb.update(i)
        return mappable
