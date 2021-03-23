"""



"""

import logging
import os
import warnings
from bisect import bisect_right
from collections import defaultdict

import intervaltree
import numpy as np
import tables
from future.utils import string_types

from genomic_regions import RegionBased, GenomicRegion
from .config import config
from .general import Maskable, MaskedTable
from .regions import LazyGenomicRegion, RegionsTable, RegionBasedWithBins
from .tools.general import RareUpdateProgressBar, create_col_index, range_overlap, str_to_int
from .tools.load import load
import datetime

logger = logging.getLogger(__name__)


class Edge(object):
    """
    A contact / an Edge between two genomic regions.

    .. attribute:: source

        The index of the "source" genomic region. By convention,
        source <= sink.

    .. attribute:: sink

        The index of the "sink" genomic region.

    .. attribute:: bias

        Bias factor obtained via normalisation of the Hi-C matrix

    .. attribute:: source_node

        The first :class:`~fanc.GenomicRegion` in this contact

    .. attribute:: sink_node

        The second :class:`~fanc.GenomicRegion` in this contact
    """
    def __init__(self, source, sink, _weight_field='weight', **kwargs):
        """
        :param source: The index of the "source" genomic region
                       or :class:`~fanc.GenomicRegion` object.
        :param sink: The index of the "sink" genomic region
                     or :class:`~fanc.GenomicRegion` object.
        :param kwargs: Other key, value pairs to be stored as
                       :class:`~Edge` attributes
        """
        object.__setattr__(self, '_source', source)
        object.__setattr__(self, '_sink', sink)
        object.__setattr__(self, 'bias', 1.)
        object.__setattr__(self, 'expected', None)
        object.__setattr__(self, '_weight_field', _weight_field)
        object.__setattr__(self, '_weight', None)

        for key, value in kwargs.items():
            setattr(self, key.decode() if isinstance(key, bytes) else key, value)

    def __getattr__(self, item):
        if item == '_weight_field' or item != self._weight_field:
            return object.__getattribute__(self, item)

        if self.expected is None:
            return object.__getattribute__(self, '_weight') * self.bias
        else:
            return (object.__getattribute__(self, '_weight') * self.bias) / self.expected

    def __setattr__(self, key, value):
        if key == object.__getattribute__(self, '_weight_field'):
            object.__setattr__(self, '_weight', value)
        else:
            object.__setattr__(self, key, value)

    def __getitem__(self, item):
        try:
            return getattr(self, item)
        except AttributeError:
            raise KeyError("No such key: {}".format(item))

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
        raise ValueError("Source not not provided during object initialization!")

    @source_node.setter
    def source_node(self, value):
        self._source = value

    @property
    def source_region(self):
        return self.source_node

    @property
    def sink_node(self):
        if isinstance(self._sink, GenomicRegion):
            return self._sink
        raise ValueError("Sink not not provided during object initialization!")

    @sink_node.setter
    def sink_node(self, value):
        self._sink = value

    @property
    def sink_region(self):
        return self.sink_node

    def __repr__(self):
        base_info = "{}--{}".format(self.source, self.sink)
        for field in dir(self):
            if not field.startswith('_') and not field == 'source' and not field == 'sink':
                try:
                    value = getattr(self, field)
                    base_info += "; {}: {}".format(field, value)
                except ValueError:
                    pass
        return base_info


class LazyEdge(object):
    """
    An :class:`~Edge` equivalent supporting lazy loading.

    .. attribute:: source

        The index of the "source" genomic region. By convention,
        source <= sink.

    .. attribute:: sink

        The index of the "sink" genomic region.

    .. attribute:: bias

        Bias factor obtained via normalisation of the Hi-C matrix

    .. attribute:: source_node

        The first :class:`~fanc.GenomicRegion` in this contact

    .. attribute:: sink_node

        The second :class:`~fanc.GenomicRegion` in this contact
    """
    def __init__(self, row, regions_table=None, _weight_field='weight'):
        self._row = row
        self._regions_table = regions_table
        self.bias = 1.
        self.expected = None
        self._weight_field = _weight_field

    def __getattr__(self, item):
        try:
            return self._row[item]
        except KeyError:
            raise AttributeError("{} does not exist in edge".format(item))

    @property
    def weight(self):
        if self.expected is None:
            return self._row[self._weight_field] * self.bias
        else:
            return (self._row[self._weight_field] * self.bias) / self.expected

    @property
    def source_node(self):
        if self._regions_table is None:
            raise RuntimeError("Must set the _regions_table attribute before calling this method!")

        source_row = self._regions_table[self.source]
        return LazyGenomicRegion(source_row)

    @property
    def sink_node(self):
        if self._regions_table is None:
            raise RuntimeError("Must set the _regions_table attribute before calling this method!")

        sink_row = self._regions_table[self.sink]
        return LazyGenomicRegion(sink_row)

    @property
    def source_region(self):
        return self.source_node

    @property
    def sink_region(self):
        return self.sink_node

    def __repr__(self):
        return "<{}.{} for row {}>".format(self.__module__, self.__class__.__name__, self._row)


class MutableLazyEdge(LazyEdge):
    def __init__(self, row, regions_table=None, _weight_field='weight'):
        self.__dict__['_row'] = row
        self.__dict__['_regions_table'] = regions_table
        self.__dict__['bias'] = 1.
        self.__dict__['expected'] = None
        self.__dict__['_weight_field'] = _weight_field

    def __setattr__(self, key, value):
        if not key.startswith('_'):
            try:
                self._row[key] = value
                return
            except KeyError:
                pass
        object.__setattr__(self, key, value)

    def update(self):
        """
        Write changes to PyTables row to file.
        """
        self._row.update()


def as_edge(edge):
    """
    Convert input to :class:`~Edge`.

    :param edge: Can be :class:`~Edge`,
                 tuple or list of the form (source, sink, weight),
                 tuple of the form (:class:`~fanc.GenomicRegion`, :class:`~fanc.GenomicRegion`),
                 dict, or :class:`~Edge` equivalent
    :return: :class:`~Edge`
    """
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
    """
    Class representing pairs of genomic regions.

    This is the basic interface for all pair and matrix classes in this module.
    It inherits all methods from :class:`~genomic_regions.RegionBased`, and is
    therefore based on a list of genomic regions (:class:`~fanc.GenomicRegion`)
    representing the underlying genome. You can use the
    :func:`~genomic_regions.RegionBased.regions` method to access genomic regions
    in a intuitive fashion, for example:

    .. code ::

        for region in rpc.regions('chr1'):
            # do something with region
            print(region)

    For more details on region access, see the :code:`genomic_regions`
    documentation, on which this module is built.

    :class:`~RegionPairsContainer` adds methods for *pairs* of genomic regions
    on top of the :class:`~genomic_regions.RegionBased` methods for individual
    regions. In the nomenclature of this module, which borrows from network
    analysis terminology, a pair of regions is represented by an :class:`~Edge`.

    .. code ::

        # iterate over all region pairs / edges in chr1
        for edge in rpc.edges(("chr1", "chr1")):
            # do something with edge / region pair
            region1 = edge.source_region
            region2 = edge.sink_region

    for more details see the :func:`~RegionPairsContainer.edges` method help.

    This class itself is only an interface and cannot actually be used to add
    regions and region pairs. Implementations of this interface, i.e. subclasses
    such as :class:`~RegionPairsTable` must override various hidden methods
    to give them full functionality.

    * :func:`~RegionPairsContainer._add_edge` is used to save region pairs / edges
      to the object. It receives a single :class:`~Edge` as input and should
      return the index of the added edge.

    * :func:`~RegionPairsContainer._edges_iter` is required by
      :func:`~RegionPairsContainer.edges`. It is used to iterate over all
      edges in the object in no particular order. It should return a generator
      of :class:`~Edge` objects representing all region pairs in the object.

    * :func:`~RegionPairsContainer._edges_subset` is also used by
      :func:`~RegionPairsContainer.edges`. It is used to iterate over a subset of
      edges in this object. It receives as input a :code:`key` representing the requested
      subset (further described in :func:`~RegionPairsContainer.edges`), and
      two lists of :class:`~fanc.GenomicRegion` objects, :code:`row_regions` and
      :code:`col_regions` representing the two dimensions of regions selected
      by :code:`key`. It should return an iterator over :class:`~Edge` objects.

    * :func:`~RegionPairsContainer._edges_getitem` is used by
      :func:`~RegionPairsContainer.edges` for retrieval of edges by bracket notation.
      For integer input, it should return a single :class:`~Edge`, for :class:`~slice`
      input a list of :class:`~Edge` objects.

    The above methods cover all the basic :class:`~RegionPairsContainer` functionality,
    but for speed improvements you may also want to override the following method,
    which by default iterates over all edges

    * :func:`~RegionPairsContainer._edges_length` which returns the total number of
      edges in the object

    """
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
        """
        Alias for :func:`~RegionPairsContainer.add_edge`

        :param contact: :class:`~Edge`
        :param args: Positional arguments passed to
                     :func:`~RegionPairsContainer._add_edge`
        :param kwargs: Keyword arguments passed to
                       :func:`~RegionPairsContainer._add_edge`
        """
        return self.add_edge(contact, *args, **kwargs)

    def add_edge(self, edge, check_nodes_exist=True, *args, **kwargs):
        """
        Add an edge / contact between two regions to this object.

        :param edge: :class:`~Edge`, dict with at least the
                     attributes source and sink, optionally weight,
                     or a list of length 2 (source, sink) or 3
                     (source, sink, weight).
        :param check_nodes_exist: Make sure that there are nodes
                                  that match source and sink indexes
        :param args: Positional arguments passed to
                     :func:`~RegionPairsContainer._add_edge`
        :param kwargs: Keyword arguments passed to
                       :func:`~RegionPairsContainer._add_edge`
        """
        if isinstance(edge, Edge):
            self._add_edge(edge, *args, **kwargs)

            if check_nodes_exist:
                n_regions = len(self.regions)
                if edge.source >= n_regions or edge.sink >= n_regions:
                    raise ValueError("Node index ({}/{}) exceeds number of nodes ({}) in object".format(
                        edge.source, edge.sink, n_regions
                    ))
        elif isinstance(edge, list) or isinstance(edge, tuple):
            self.add_edge_from_list(edge)
        elif isinstance(edge, dict):
            self.add_edge_from_dict(edge)
        else:
            edge = as_edge(edge)
            self.add_edge_from_edge(edge)

    def add_edge_from_list(self, edge, *args, **kwargs):
        """
        Direct method to add an edge from list or tuple input.

        :param edge: List or tuple. Should be of length 2
                     (source, sink) or 3 (source, sink, weight)
        """
        return self.add_edge(as_edge(edge), *args, **kwargs)

    def add_edge_from_dict(self, edge, *args, **kwargs):
        """
        Direct method to add an edge from dict input.

        :param edge: dict with at least the keys "source"
                     and "sink". Additional keys will be loaded
                     as edge attributes
        """
        return self.add_edge(as_edge(edge), *args, **kwargs)

    def add_edge_from_edge(self, edge, *args, **kwargs):
        """
        Direct method to add an edge from :class:`~Edge` input.

        :param edge: :class:`~Edge`
        """
        return self.add_edge(as_edge(edge), *args, **kwargs)

    def add_edge_simple(self, source, sink, weight=None, *args, **kwargs):
        """
        Direct method to add an edge from :class:`~Edge` input.

        :param source: Source region index
        :param sink: Sink region index
        :param weight: Weight of the edge
        """
        return self.add_edge(Edge(source=source, sink=sink, weight=weight), *args, **kwargs)

    def add_edges(self, edges, *args, **kwargs):
        """
        Bulk-add edges from a list.

        List items can be any of the supported edge types,
        list, tuple, dict, or :class:`~Edge`. Repeatedly
        calls :func:`~RegionPairsContainer.add_edge`, so
        may be inefficient for large amounts of data.

        :param edges: List (or iterator) of edges. See
                      :func:`~RegionMatrixTable.add_edge`
                      for details
        """
        for edge in edges:
            self.add_edge(edge, *args, **kwargs)

    def add_contacts(self, contacts, *args, **kwargs):
        """
        Alias for :func:`~RegionPairsTable.add_edges`
        """
        return self.add_edges(contacts, *args, **kwargs)

    @property
    def edges(self):
        """
        Iterate over contacts / edges.

        :func:`~RegionPairsContainer.edges` is the central function of
        :class:`~RegionPairsContainer`. Here, we will use the
        :class:`~fanc.Hic` implementation for demonstration purposes,
        but the usage is exactly the same for all compatible
        objects implementing :class:`~RegionPairsContainer`, including
        :class:`~fanc.compatibility.juicer.JuicerHic` and
        :class:`~fanc.compatibility.cooler.CoolerHic`.

        .. code ::

            import fanc

            # file from FAN-C examples
            hic = fanc.load("output/hic/binned/fanc_example_1mb.hic")

        We can easily find the number of edges in the sample
        :class:`~fanc.Hic` object:

        .. code ::

            len(hic.edges)  # 8695

        When used in an iterator context, :func:`~RegionPairsContainer.edges`
        iterates over all edges in the :class:`~RegionPairsContainer`:

        .. code ::

            for edge in hic.edges:
                # do something with edge
                print(edge)
                # 42--42; bias: 5.797788472650082e-05; sink_node: chr18:42000001-43000000; source_node: chr18:42000001-43000000; weight: 0.12291311562018173
                # 24--28; bias: 6.496381719803623e-05; sink_node: chr18:28000001-29000000; source_node: chr18:24000001-25000000; weight: 0.025205961072838057
                # 5--76; bias: 0.00010230955745211447; sink_node: chr18:76000001-77000000; source_node: chr18:5000001-6000000; weight: 0.00961709840049876
                # 66--68; bias: 8.248432587969082e-05; sink_node: chr18:68000001-69000000; source_node: chr18:66000001-67000000; weight: 0.03876763316345468
                # ...

        Calling :func:`~RegionPairsContainer.edges` as a method has the
        same effect:

        .. code ::

            # note the '()'
            for edge in hic.edges():
                # do something with edge
                print(edge)
                # 42--42; bias: 5.797788472650082e-05; sink_node: chr18:42000001-43000000; source_node: chr18:42000001-43000000; weight: 0.12291311562018173
                # 24--28; bias: 6.496381719803623e-05; sink_node: chr18:28000001-29000000; source_node: chr18:24000001-25000000; weight: 0.025205961072838057
                # 5--76; bias: 0.00010230955745211447; sink_node: chr18:76000001-77000000; source_node: chr18:5000001-6000000; weight: 0.00961709840049876
                # 66--68; bias: 8.248432587969082e-05; sink_node: chr18:68000001-69000000; source_node: chr18:66000001-67000000; weight: 0.03876763316345468
                # ...

        Rather than iterate over all edges in the object, we can select only a subset.
        If the key is a string or a :class:`~fanc.GenomicRegion`, all non-zero edges connecting
        the region described by the key to any other region are returned. If the key is a
        tuple of strings or :class:`~fanc.GenomicRegion`, only edges between the two regions
        are returned.

        .. code ::

            # select all edges between chromosome 19
            # and any other region:
            for edge in hic.edges("chr19"):
                print(edge)
                # 49--106; bias: 0.00026372303696871666; sink_node: chr19:27000001-28000000; source_node: chr18:49000001-50000000; weight: 0.003692122517562033
                # 6--82; bias: 0.00021923129703834945; sink_node: chr19:3000001-4000000; source_node: chr18:6000001-7000000; weight: 0.0008769251881533978
                # 47--107; bias: 0.00012820949175399097; sink_node: chr19:28000001-29000000; source_node: chr18:47000001-48000000; weight: 0.0015385139010478917
                # 38--112; bias: 0.0001493344481069762; sink_node: chr19:33000001-34000000; source_node: chr18:38000001-39000000; weight: 0.0005973377924279048
                # ...

            # select all edges that are only on
            # chromosome 19
            for edge in hic.edges(('chr19', 'chr19')):
                print(edge)
                # 90--116; bias: 0.00021173151730025176; sink_node: chr19:37000001-38000000; source_node: chr19:11000001-12000000; weight: 0.009104455243910825
                # 135--135; bias: 0.00018003890596887822; sink_node: chr19:56000001-57000000; source_node: chr19:56000001-57000000; weight: 0.10028167062466517
                # 123--123; bias: 0.00011063368998965993; sink_node: chr19:44000001-45000000; source_node: chr19:44000001-45000000; weight: 0.1386240135570439
                # 92--93; bias: 0.00040851066434864896; sink_node: chr19:14000001-15000000; source_node: chr19:13000001-14000000; weight: 0.10090213409411629
                # ...

            # select inter-chromosomal edges
            # between chromosomes 18 and 19
            for edge in hic.edges(('chr18', 'chr19')):
                print(edge)
                # 49--106; bias: 0.00026372303696871666; sink_node: chr19:27000001-28000000; source_node: chr18:49000001-50000000; weight: 0.003692122517562033
                # 6--82; bias: 0.00021923129703834945; sink_node: chr19:3000001-4000000; source_node: chr18:6000001-7000000; weight: 0.0008769251881533978
                # 47--107; bias: 0.00012820949175399097; sink_node: chr19:28000001-29000000; source_node: chr18:47000001-48000000; weight: 0.0015385139010478917
                # 38--112; bias: 0.0001493344481069762; sink_node: chr19:33000001-34000000; source_node: chr18:38000001-39000000; weight: 0.0005973377924279048
                # ...

        By default, :func:`~RegionPairsContainer.edges` will retrieve all edge attributes,
        which can be slow when iterating over a lot of edges. This is why all file-based FAN-C
        :class:`~RegionPairsContainer` objects support lazy loading, where attributes
        are only read on demand.

        .. code ::

            for edge in hic.edges('chr18', lazy=True):
                print(edge.source, edge.sink, edge.weight, edge)
                # 42 42 0.12291311562018173 <fanc.matrix.LazyEdge for row /edges/chrpair_0_0.row (Row), pointing to row #0>
                # 24 28 0.025205961072838057 <fanc.matrix.LazyEdge for row /edges/chrpair_0_0.row (Row), pointing to row #1>
                # 5 76 0.00961709840049876 <fanc.matrix.LazyEdge for row /edges/chrpair_0_0.row (Row), pointing to row #2>
                # 66 68 0.03876763316345468 <fanc.matrix.LazyEdge for row /edges/chrpair_0_0.row (Row), pointing to row #3>
                # ...

        .. warning :: The lazy iterator reuses the :class:`~LazyEdge` object in every iteration,
                      and overwrites the :class:`~LazyEdge` attributes. Therefore **do not** use
                      lazy iterators if you need to store edge objects for later access.
                      For example, the following code works as expected
                      :code:`list(hic.edges())`, with all :class:`~Edge` objects stored in the
                      list, while this code :code:`list(hic.edges(lazy=True))`
                      will result in a list of identical :class:`~LazyEdge` objects. Always ensure
                      you do all edge processing in the loop when working with lazy iterators!

        When working with normalised contact frequencies, such as obtained through
        matrix balancing in the example above, :func:`~RegionPairsContainer.edges`
        automatically returns normalised edge weights. In addition, the :code:`bias`
        attribute will (typically) have a value different from 1.

        When you are interested in the raw contact frequency, use the :code:`norm=False`
        parameter:

        .. code ::

            for edge in hic.edges('chr18', lazy=True, norm=False):
                print(edge.source, edge.sink, edge.weight)
                # 42 42 2120.0
                # 24 28 388.0
                # 5 76 94.0
                # 66 68 470.0
                # ...

        You can also choose to omit all intra- or inter-chromosomal edges using
        :code:`intra_chromosomal=False` or :code:`inter_chromosomal=False`, respectively.

        :return: Iterator over :class:`~Edge` or equivalent.
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
                check_valid = kwargs.pop('check_valid', True)
                oe = kwargs.pop('oe', False)
                oe_per_chromosome = kwargs.pop('oe_per_chromosome', True)

                start = datetime.datetime.now()

                if norm and hasattr(self._regions_pairs, 'bias_vector'):
                    bias = self._regions_pairs.bias_vector()
                else:
                    bias = np.repeat(1., len(self._regions_pairs.regions))

                if oe:
                    if not hasattr(self._regions_pairs, 'expected_values'):
                        raise ValueError("Cannot perform O/E transformation because this object does not "
                                         "support the expected_values function!")
                    expected_genome, expected_intra, expected_inter = self._regions_pairs.expected_values(norm=norm)
                    if not oe_per_chromosome and expected_genome is None:
                        raise ValueError("Expected values were not calculated for the whole genome. "
                                         "Please note that this option "
                                         "is only supported for FAN-C Hi-C files. If this is a FAN-C file, "
                                         "please report this as a bug.")
                else:
                    expected_genome, expected_intra, expected_inter = None, None, None

                valid = [getattr(r, 'valid', True) for r in self._regions_pairs.regions(lazy=True)]

                # getting regions
                row_regions, col_regions = self._regions_pairs._key_to_regions(key)
                if isinstance(row_regions, GenomicRegion):
                    row_regions = [row_regions]
                if isinstance(col_regions, GenomicRegion):
                    col_regions = [col_regions]

                row_regions_by_chromosome = defaultdict(list)
                for r in row_regions:
                    row_regions_by_chromosome[r.chromosome].append(r)

                col_regions_by_chromosome = defaultdict(list)
                for r in col_regions:
                    col_regions_by_chromosome[r.chromosome].append(r)

                d = datetime.datetime.now() - start
                # print("Startup: {}".format(d.total_seconds()))

                chromosome_pairs = set()
                for row_chromosome, row_chromosome_regions in row_regions_by_chromosome.items():
                    for col_chromosome, col_chromosome_regions in col_regions_by_chromosome.items():
                        if (col_chromosome, row_chromosome) in chromosome_pairs:
                            continue
                        chromosome_pairs.add((row_chromosome, col_chromosome))

                        if row_chromosome == col_chromosome:
                            if not intra_chromosomal:
                                continue
                            if oe:
                                ex = expected_intra[row_chromosome] if oe_per_chromosome else expected_genome
                            else:
                                ex = np.repeat(None, len(self._regions_pairs.regions))
                        else:
                            if not inter_chromosomal and row_chromosome != col_chromosome:
                                continue
                            ex = np.repeat(expected_inter, len(self._regions_pairs.regions))

                        for edge in self._regions_pairs._edges_subset(
                                (GenomicRegion(row_chromosome,
                                               start=row_chromosome_regions[0].start,
                                               end=row_chromosome_regions[-1].end),
                                 GenomicRegion(col_chromosome,
                                               start=col_chromosome_regions[0].start,
                                               end=col_chromosome_regions[-1].end)),
                                row_chromosome_regions, col_chromosome_regions,
                                *args, **kwargs):
                            source, sink = edge.source, edge.sink
                            if check_valid and (not valid[source] or not valid[sink]):
                                continue
                            edge.bias = bias[source] * bias[sink]
                            edge.expected = ex[abs(sink - source)]
                            yield edge

            def __len__(self):
                return self._regions_pairs._edges_length()

        return EdgeIter(self)

    def edges_dict(self, *args, **kwargs):
        """
        Edges iterator with access by bracket notation.

        This iterator **always** returns unnormalised edges.

        :return: dict or dict-like iterator
        """
        kwargs['norm'] = False
        return self.edges(*args, **kwargs)

    def edge_subset(self, key=None, *args, **kwargs):
        """
        Get a subset of edges.

        This is an alias for :func:`~RegionPairsContainer.edges`.

        :return: generator (:class:`~Edge`)
        """
        return self.edges(key, *args, **kwargs)

    @staticmethod
    def regions_identical(pairs):
        """
        Check if the regions in all objects in the list are identical.

        :param pairs: :class:`~list` of :class:`~genomic_regions.RegionBased`
                      objects
        :return: True if chromosome, start, and end are identical between
                 all regions in the same list positions.
        """
        logger.info("Checking if regions are identical")
        regions = list(pairs[0].regions)
        for matrix_object in pairs[1:]:
            try:
                for r1, r2 in zip(regions, matrix_object.regions):
                    assert r1.chromosome == r2.chromosome
                    assert r1.start == r2.start
                    assert r1.end == r2.end
            except AssertionError:
                return False
        return True

    @classmethod
    def merge(cls, pairs, *args, **kwargs):
        """
        Merge two or more :class:`~RegionPairsContainer` objects.

        :param pairs: :class:`~list` of :class:`~RegionPairsContainer`
        :param args: Positional arguments passed to constructor of this
                     class
        :param kwargs: Keyword arguments passed to constructor of this
                       class
        """
        if 'mode' not in kwargs:
            kwargs['mode'] = 'w'
        merged_pairs = cls(*args, **kwargs)

        pairs = [pair_object for pair_object in pairs]

        if not RegionPairsContainer.regions_identical(pairs):
            raise ValueError("Regions in pair objects are not identical, cannot perform merge!")

        merged_pairs.add_regions(pairs[0].regions(lazy=True))

        for pair_object in pairs:
            merged_pairs.add_edges(pair_object.edges(lazy=True))

        return merged_pairs

    def edge_data(self, attribute, *args, **kwargs):
        """
        Iterate over specific edge attribute.

        :param attribute: Name of the attribute, e.g. "weight"
        :param args: Positional arguments passed to :func:`~RegionPairsContainer.edges`
        :param kwargs: Keyword arguments passed to :func:`~RegionPairsContainer.edges`
        :return: iterator over edge attribute
        """
        kwargs.setdefault('lazy', True)
        for edge in self.edges(*args, **kwargs):
            yield getattr(edge, attribute)

    def regions_and_edges(self, key, *args, **kwargs):
        """
        Convenient access to regions and edges selected by key.

        :param key: Edge selector, see :func:`~RegionPairsContainer.edges`
        :param args: Positional arguments passed to :func:`~RegionPairsContainer.edges`
        :param kwargs: Keyword arguments passed to :func:`~RegionPairsContainer.edges`
        :return: list of row regions, list of col regions, iterator over edges
        """
        row_regions, col_regions = self._key_to_regions(key)
        if isinstance(row_regions, GenomicRegion):
            row_regions = [row_regions]
        else:
            row_regions = list(row_regions)

        if isinstance(col_regions, GenomicRegion):
            col_regions = [col_regions]
        else:
            col_regions = list(col_regions)

        edges = self.edges((row_regions, col_regions), *args, **kwargs)

        return row_regions, col_regions, edges

    def mappable(self, region=None):
        """
        Get the mappability of regions in this object.

        A "mappable" region has at least one contact to another region
        in the genome.

        :return: :class:`~np.array` where True means mappable
                 and False unmappable
        """
        return np.array([True if getattr(r, 'valid', True) else False
                         for r in self.regions(region, lazy=True)])


class RegionMatrixContainer(RegionPairsContainer, RegionBasedWithBins):
    """
    Class representing matrices where pixels correspond to genomic region pairs.

    This is the common interface for all matrix-based classes, such as
    :class:`~fanc.Hic` or :class:`~fanc.FoldChangeMatrix`. It provides
    access to specialised matrix methods, most importantly
    :func:`~RegionMatrixContainer.matrix`, which assembles :mod:`numpy`
    arrays from the list of pairwise contacts stored in each object.

    It inherits all region methods from :class:`~genomic_regions.RegionBased`,
    and all edge/contact methods from :class:`~RegionPairsContainer`.
    You can use the same type of keys for :func:`~RegionMatrixContainer.matrix`
    that you would use for :func:`~RegionPairsContainer.edges`, and additionally
    have the option to retrieve the observed/expected matrix.

    .. code ::

        import fanc
        hic = fanc.load("output/hic/binned/fanc_example_1mb.hic")

        # get the whole-genome matrix
        m = hic.matrix()
        type(m)  # fanc.matrix.RegionMatrix
        isinstance(m, np.ndarray)  # True
        m.shape  # 139, 139

        # get just the chromosome 18 intra-chromosomal matrix
        m = hic.matrix(('chr18', 'chr18'))
        m.shape  # 79, 79

        # get all rows of the whole-genome matrix
        # corresponding to chromosome 18
        m = hic.matrix('chr18')
        m.shape  # 79, 139

        # get unnormalised chromosome 18 matrix
        m = hic.matrix(('chr18', 'chr18'), norm=False)

        # get chromosome 18 O/E matrix
        m = hic.matrix(('chr18', 'chr18'), oe=True)

        # get log2-transformed chromosome 18 O/E matrix
        m = hic.matrix(('chr18', 'chr18'), oe=True, log=True)

    """
    def __init__(self):
        RegionPairsContainer.__init__(self)
        self._default_value = 0.0
        self._default_score_field = 'weight'

    def regions_and_matrix_entries(self, key=None, score_field=None, *args, **kwargs):
        """
        Convenient access to non-zero matrix entries and associated regions.

        :param key: Edge key, see :func:`~RegionPairsContainer.edges`
        :param oe: If True, will divide observed values by their expected value
                   at the given distance. False by default
        :param oe_per_chromosome: If True (default), will do a per-chromosome O/E
                                  calculation rather than using the whole matrix
                                  to obtain expected values
        :param score_field: (optional) any edge attribute that returns a number
                            can be specified here for filling the matrix. Usually
                            this is defined by the :code:`_default_score_field`
                            attribute of the matrix class.
        :param args: Positional arguments passed to :func:`~RegionPairsContainer.edges`
        :param kwargs: Keyword arguments passed to :func:`~RegionPairsContainer.edges`
        :return: list of row regions, list of col regions, iterator over (i, j, weight) tuples
        """
        row_regions, col_regions, edges_iter = self.regions_and_edges(key, *args, **kwargs)

        try:
            row_offset = row_regions[0].ix
            col_offset = col_regions[0].ix
        except IndexError:
            return row_regions, col_regions, []

        if score_field is None:
            score_field = self._default_score_field

        def offset_iter(edge_iter):
            for edge in edge_iter:
                source, sink, weight = edge.source, edge.sink, getattr(edge, score_field, self._default_value)

                i = source - row_offset
                j = sink - col_offset
                if i >= 0 and j >= 0:
                    yield source, sink, i, j, weight

                l = source - col_offset
                k = sink - row_offset
                if (i, j) != (k, l) and k >= 0 and l >= 0:
                    yield source, sink, k, l, weight

        #edges_iter = self.edges((row_regions, col_regions), *args, **kwargs)
        entry_iter = ((i, j, weight)
                      for _, _, i, j, weight in offset_iter(edges_iter))

        return row_regions, col_regions, entry_iter

    def matrix(self, key=None,
               log=False,
               default_value=None, mask=True, log_base=2,
               *args, **kwargs):
        """
        Assemble a :class:`~RegionMatrix` from region pairs.

        :param key: Matrix selector. See :func:`~fanc.matrix.RegionPairsContainer.edges`
                    for all supported key types
        :param log: If True, log-transform the matrix entries. Also see log_base
        :param log_base: Base of the log transformation. Default: 2; only used when
                         log=True
        :param default_value: (optional) set the default value of matrix entries
                              that have no associated edge/contact
        :param mask: If False, do not mask unmappable regions
        :param args: Positional arguments passed to
                     :func:`~fanc.matrix.RegionMatrixContainer.regions_and_matrix_entries`
        :param kwargs: Keyword arguments passed to
                       :func:`~fanc.matrix.RegionMatrixContainer.regions_and_matrix_entries`
        :return: :class:`~fanc.matrix.RegionMatrix`
        """

        if default_value is None:
            default_value = self._default_value

        if kwargs.get('oe', False):
            default_value = 1.0

        kwargs['lazy'] = True
        row_regions, col_regions, matrix_entries = self.regions_and_matrix_entries(key,
                                                                                   *args,
                                                                                   **kwargs)

        m = np.full((len(row_regions), len(col_regions)), default_value)

        for source, sink, weight in matrix_entries:
            ir = source
            jr = sink
            if 0 <= ir < m.shape[0] and 0 <= jr < m.shape[1]:
                m[ir, jr] = weight

        if log:
            m = np.log(m) / np.log(log_base)
            m[~np.isfinite(m)] = default_value

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
        """
        Calculate the possible number of contacts in the genome.

        This calculates the number of potential region pairs in
        a genome for any possible separation distance, taking into
        account the existence of unmappable regions.

        It will calculate one number for inter-chromosomal pairs,
        return a list with the number of possible pairs where the
        list index corresponds to the number of bins separating two regions,
        and a dictionary of lists for each chromosome.

        :return: possible intra-chromosomal pairs,
                 possible intra-chromosomal pairs by chromosome,
                 possible inter-chromosomal pairs
        """
        logger.debug("Calculating possible counts")

        logger.debug("Setup for possible counts")
        chromosomes = self.chromosomes()
        mappability = self.mappable()
        cb = self.chromosome_bins

        max_distance = 0
        chromosome_intra_total = dict()
        chromosome_mappable_counts = dict()
        for chromosome in chromosomes:
            logger.debug("Possible counts for {}".format(chromosome))
            chromosome_start_bin, chromosome_end_bin = cb[chromosome]
            mappable_chromosome = np.array(mappability[chromosome_start_bin:chromosome_end_bin])
            max_distance = max(max_distance, chromosome_end_bin - chromosome_start_bin)

            possible_by_distance = np.zeros(chromosome_end_bin - chromosome_start_bin)
            sub = np.ones(chromosome_end_bin - chromosome_start_bin)
            d = len(mappable_chromosome)
            for i, mappable in enumerate(mappable_chromosome):
                if mappable:
                    possible_by_distance[:(d - i)] += np.ones(d - i)
                else:
                    # subtract vertical
                    sub[i] = 0
                    possible_by_distance[:(i+1)] -= sub[:(i+1)][::-1]

            chromosome_intra_total[chromosome] = possible_by_distance
            chromosome_mappable_counts[chromosome] = np.sum(mappable_chromosome)

        possible_by_distance_whole_matrix = defaultdict(int)
        for possible_by_distance in chromosome_intra_total.values():
            for distance, count in enumerate(possible_by_distance):
                possible_by_distance_whole_matrix[distance] += count

        intra_total = [possible_by_distance_whole_matrix[i]
                       for i in range(max_distance)]

        inter_total = 0
        for i in range(len(chromosomes)):
            chromosome1 = chromosomes[i]
            for j in range(i + 1, len(chromosomes)):
                chromosome2 = chromosomes[j]
                inter_total += chromosome_mappable_counts[chromosome1] * chromosome_mappable_counts[chromosome2]

        return intra_total, chromosome_intra_total, inter_total

    def expected_values_and_marginals(self, selected_chromosome=None, norm=True,
                                      *args, **kwargs):
        """
        Calculate the expected values for genomic contacts at all distances
        and the whole matrix marginals.

        This calculates the expected values between genomic regions
        separated by a specific distance. Expected values are calculated
        as the average weight of edges between region pairs with the same
        genomic separation, taking into account unmappable regions.

        It will return a tuple with three values: a list of genome-wide
        intra-chromosomal expected values (list index corresponds to number
        of separating bins), a dict with chromosome names as keys and
        intra-chromosomal expected values specific to each chromosome, and
        a float for inter-chromosomal expected value.

        :param selected_chromosome: (optional) Chromosome name. If provided,
                                    will only return expected values for this
                                    chromosome.
        :param norm: If False, will calculate the expected values on the
                     unnormalised matrix.
        :param args: Not used in this context
        :param kwargs: Not used in this context
        :return: list of intra-chromosomal expected values,
                 dict of intra-chromosomal expected values by chromosome,
                 inter-chromosomal expected value
        """
        weight_field = getattr(self, '_default_score_field', None)
        default_value = getattr(self, '_default_value', 1.)

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
        valid = [False] * len(self.regions)
        inter_sums = 0.0
        intra_sums = [0.0] * max_distance
        with RareUpdateProgressBar(max_value=len(self.edges), prefix='Expected') as pb:
            for i, edge in enumerate(self.edges(lazy=True, norm=norm, check_valid=False)):
                source, sink = edge.source, edge.sink
                try:
                    weight = getattr(edge, weight_field)
                except AttributeError:
                    weight = default_value

                source_chromosome = chromosome_dict[source]
                sink_chromosome = chromosome_dict[sink]

                marginals[source] += weight
                marginals[sink] += weight
                if weight != self._default_value:
                    valid[source] = True
                    valid[sink] = True

                if sink_chromosome != source_chromosome:
                    inter_sums += weight
                else:
                    distance = sink - source
                    intra_sums[distance] += weight
                    chromosome_intra_sums[source_chromosome][distance] += weight
                pb.update(i)

        intra_total, chromosome_intra_total, inter_total = self.possible_contacts()

        # expected values
        inter_expected = 0 if inter_total == 0 else inter_sums / inter_total

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
                    chromosome_intra_expected[chromosome][d] = chromosome_intra_sums[chromosome][
                                                                   d] / chromosome_count

        if selected_chromosome is not None:
            return chromosome_intra_expected[selected_chromosome], marginals, valid

        return intra_expected, chromosome_intra_expected, inter_expected, marginals, valid

    def expected_values(self, selected_chromosome=None, norm=True, *args, **kwargs):
        """
        Calculate the expected values for genomic contacts at all distances.

        This calculates the expected values between genomic regions
        separated by a specific distance. Expected values are calculated
        as the average weight of edges between region pairs with the same
        genomic separation, taking into account unmappable regions.

        It will return a tuple with three values: a list of genome-wide
        intra-chromosomal expected values (list index corresponds to number
        of separating bins), a dict with chromosome names as keys and
        intra-chromosomal expected values specific to each chromosome, and
        a float for inter-chromosomal expected value.

        :param selected_chromosome: (optional) Chromosome name. If provided,
                                    will only return expected values for this
                                    chromosome.
        :param norm: If False, will calculate the expected values on the
                     unnormalised matrix.
        :param args: Not used in this context
        :param kwargs: Not used in this context
        :return: list of intra-chromosomal expected values,
                 dict of intra-chromosomal expected values by chromosome,
                 inter-chromosomal expected value

        """
        result = self.expected_values_and_marginals(selected_chromosome=selected_chromosome,
                                                    norm=norm, *args, **kwargs)
        return result[:-2]

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
            for r in row_regions + col_regions:
                mask[r.ix - min_ix] = not r.valid
            marginals = np.ma.masked_where(mask, marginals)

        return marginals

    def scaling_factor(self, matrix, weight_column=None):
        """
        Compute the scaling factor to another matrix.

        Calculates the ratio between the number of contacts in
        this Hic object to the number of contacts in another
        Hic object.

        :param matrix: A :class:`~Hic` object
        :param weight_column: Name of the column to calculate the scaling factor on
        :return: float
        """
        if weight_column is None:
            weight_column = self._default_score_field

        logger.info("Calculating scaling factor...")
        m1_sum = 0
        for v1 in self.edge_data(weight_column, lazy=True):
            if np.isfinite(v1):
                m1_sum += v1

        m2_sum = 0
        for v2 in matrix.edge_data(weight_column, lazy=True):
            if np.isfinite(v2):
                m2_sum += v2

        scaling_factor = m1_sum / m2_sum
        logger.debug("Scaling factor: {}/{} = {}".format(m1_sum, m2_sum, scaling_factor))
        return scaling_factor

    def deepcopy(self, target_class=None, bias=True, **kwargs):
        cls = self.__class__ if target_class is None else target_class

        copy = cls(**kwargs)
        copy.add_regions(self.regions(lazy=True))

        total = len(self.edges)
        with RareUpdateProgressBar(max_value=total) as pb:
            for i, edge in enumerate(self.edges(lazy=True, norm=False, oe=False)):
                copy.add_edge_simple(edge.source, edge.sink, edge.weight)
                pb.update(i)
        copy.flush(update_mappability=False)

        needs_update = True
        try:
            if bias and hasattr(self, 'bias_vector') and hasattr(copy, 'bias_vector'):
                copy.bias_vector(self.bias_vector())
                needs_update = False
            copy.flush(update_mappability=False)
        except IndexError:
            warnings.warn("Could not copy index vector. Hic version may be too old. "
                          "Please run matrix balancing again after deepcopy!")

        if needs_update:
            copy._update_mappability()

        return copy


class TableBuffer(object):
    def __init__(self, matrix, buffer_size='3G', large_distance=3, large_fraction=0.5):
        self._matrix = matrix
        self._buffer_size = buffer_size
        self._template_row = None
        self._buffer = dict()
        self._counter = dict()
        self._buffer_size = buffer_size
        self._small_buffer_size = 0
        self._large_buffer_size = 0
        self._large_distance = large_distance
        self._large_fraction = large_fraction
        self._colnames = None
        self._colindices = None
        self._small_buffer_size = None
        self._large_buffer_size = None
        self._is_initialised = False
        self._source_field = None
        self._sink_field = None
        self._weight_field = None

    def initialise_buffers(self):
        logger.debug("Initialising edge buffer!")
        matrix_table = self._matrix._edge_table(0, 0)
        self._colnames = matrix_table.colnames
        self._colindices = {name: ix for ix, name in enumerate(self._colnames)}
        self._source_field = self._colindices['source']
        self._sink_field = self._colindices['sink']
        if self._matrix._default_score_field is not None:
            self._weight_field = self._colindices[self._matrix._default_score_field]
        else:
            try:
                self._weight_field = self._colnames.index('weight')
            except ValueError:
                pass

        dtypes = [(name, matrix_table.coldtypes[name]) for name in self._colnames]
        self._template_row = np.empty(1, dtype=dtypes)[0]

        for i, name in enumerate(self._colnames):
            self._template_row[i] = matrix_table.coldflts[name]

        n_partitions = int(len(self._matrix._partition_breaks) ** 2 / 2 +
                           len(self._matrix._partition_breaks))
        n_large_partitions = max(1, self._large_distance * len(self._matrix._partition_breaks))
        buffer_size_bytes = str_to_int(self._buffer_size)

        self._large_buffer_size = int(buffer_size_bytes * self._large_fraction /
                                      n_large_partitions / self._template_row.nbytes)
        self._small_buffer_size = int(buffer_size_bytes * (1 - self._large_fraction) /
                                      max(1, n_partitions - n_large_partitions) /
                                      self._template_row.nbytes)

        logger.debug("Partitions: {} ({})".format(n_partitions, n_large_partitions))
        logger.debug("Buffer sizes ({}): {}/{}".format(self._buffer_size, self._large_buffer_size,
                                                       self._small_buffer_size))

        self._is_initialised = True

    def _reset_buffer_table(self, partition):
        self._counter[partition] = 0
        if partition[1] - partition[0] < self._large_distance:
            self._buffer[partition] = np.repeat(self._template_row, max(500, self._large_buffer_size))
        else:
            self._buffer[partition] = np.repeat(self._template_row, max(500, self._small_buffer_size))

    def flush(self, partition=None):
        if not self._is_initialised:
            return

        if partition is None:
            logger.debug("Flushing all buffers")
            partitions = list(self._buffer.keys())
        else:
            logger.debug("Flushing buffer partition {}".format(partition))
            partitions = [partition]

        with RareUpdateProgressBar(max_value=len(partitions),
                                   silent=config.hide_progressbars or partition is not None,
                                   prefix="Buffers") as pb:
            for i, partition in enumerate(partitions):
                edge_table = self._matrix._edge_table(partition[0], partition[1], create_index=False)
                buffer_table = self._buffer[partition]
                ix = self._counter[partition]
                flush = False
                if ix == buffer_table.shape[0]:
                    edge_table.append(buffer_table)
                    flush = True
                elif ix > 0:
                    edge_table.append(buffer_table[:ix])
                    flush = True
                if flush:
                    edge_table.flush(update_index=False)
                del self._buffer[partition]
                del self._counter[partition]
                pb.update(i)

    def _current_buffer_row(self, partition):
        if not self._matrix._edges_dirty:
            logger.debug("Disabling edge indexes")
            self._matrix._edges_dirty = True
            self._matrix._disable_edge_indexes()

        try:
            ix = self._counter[partition]
        except KeyError:
            if not self._is_initialised:
                self.initialise_buffers()
            self._reset_buffer_table(partition)
            ix = 0

        try:
            row = self._buffer[partition][ix]
        except IndexError:
            self.flush(partition=partition)
            self._reset_buffer_table(partition)
            ix = 0
            row = self._buffer[partition][ix]

        self._counter[partition] += 1

        return row

    def add(self, edge):
        if isinstance(edge, Edge):
            self.add_edge(edge)
        elif isinstance(edge, dict):
            self.add_dict(edge)
        elif isinstance(edge, list) or isinstance(edge, tuple) or isinstance(edge, np.ndarray):
            self.add_list(edge)
        else:
            raise ValueError("Edge format ({}) not supported!".format(type(edge)))

    def add_edge(self, edge, partition=None):
        if partition is None:
            partition = self._matrix._get_edge_table_tuple(edge.source, edge.sink)

        row = self._current_buffer_row(partition)
        for i, name in enumerate(self._colnames):
            try:
                row[i] = getattr(edge, name)
            except AttributeError:
                pass

    def add_list(self, edge_list, partition=None):
        if partition is None:
            partition = self._matrix._get_edge_table_tuple(edge_list[self._colindices['source']],
                                                           edge_list[self._colindices['sink']])

        row = self._current_buffer_row(partition)
        for i, value in enumerate(edge_list):
            row[i] = value

    def add_dict(self, edge_dict, partition=None):
        if partition is None:
            partition = self._matrix._get_edge_table_tuple(edge_dict['source'],
                                                           edge_dict['sink'])

        row = self._current_buffer_row(partition)
        for name in edge_dict.keys():
            try:
                ix = self._colindices[name]
                row[ix] = edge_dict[name]
            except KeyError:
                continue

    def add_weight(self, source, sink, weight=None, partition=None):
        if partition is None:
            partition = self._matrix._get_edge_table_tuple(source, sink)

        row = self._current_buffer_row(partition)
        row[self._source_field] = source
        row[self._sink_field] = sink
        if weight is not None:
            row[self._weight_field] = weight


class RegionPairsTable(RegionPairsContainer, Maskable, RegionsTable):
    """
    HDF5 implementation of the :class:`~RegionPairsContainer` interface.
    """

    _classid = 'REGIONPAIRSTABLE'

    def __init__(self, file_name=None, mode='a', tmpdir=None,
                 additional_region_fields=None, additional_edge_fields=None,
                 partition_strategy='auto',
                 _table_name_regions='regions', _table_name_edges='edges',
                 _edge_buffer_size=config.edge_buffer_size, _edge_table_prefix='chrpair_'):
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
        self._partition_strategy = partition_strategy
        self._edge_table_prefix = _edge_table_prefix

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

        if file_exists and mode != 'w':
            # retrieve edge tables and partitions
            self._edges = self.file.get_node('/', _table_name_edges)
            self._partition_breaks = getattr(self.meta, 'partition_breaks', None)
            self._partition_strategy = getattr(self.meta, 'partition_strategy', 'chromosome')
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

            self._edge_table(0, 0, fields=basic_fields)

        # update field names
        self._source_field_ix = 0
        self._sink_field_ix = 0
        self.field_names = []
        self._field_names_dict = dict()
        self._edge_field_defaults = dict()
        self._update_field_names()

        # set up edge buffer
        self._edge_buffer = TableBuffer(self, buffer_size=_edge_buffer_size)

    def _edge_table(self, source_partition, sink_partition, fields=None, create_if_missing=True, create_index=True):
        """
        Create and register an edge table for a partition combination.
        """
        edge_table_name = self._edge_table_prefix + str(source_partition) + '_' + str(sink_partition)
        try:
            return getattr(self._edges, edge_table_name)
        except tables.NoSuchNodeError:
            if not create_if_missing:
                raise ValueError("The edge table {}/{} cannot be found!".format(source_partition,
                                                                                sink_partition))

        if fields is None:
            fields = self._edge_table(0, 0).coldescrs

        edge_table = MaskedTable(self._edges,
                                 edge_table_name,
                                 fields, ignore_reserved_fields=True,
                                 expectedrows=10000000,
                                 create_mask_index=create_index)
        edge_table.attrs['source_partition'] = source_partition
        edge_table.attrs['sink_partition'] = sink_partition

        # index
        if create_index:
            create_col_index(edge_table.cols.source)
            create_col_index(edge_table.cols.sink)

        return edge_table

    def _has_edge_table(self, source_partition, sink_partition):
        edge_table_name = self._edge_table_prefix + str(source_partition) + '_' + str(sink_partition)
        try:
            getattr(self._edges, edge_table_name)
        except tables.NoSuchNodeError:
            return False
        return True

    def _iter_edge_tables(self):
        if self._partition_breaks is None:
            return

        for source_partition in range(len(self._partition_breaks) + 1):
            for sink_partition in range(source_partition, len(self._partition_breaks) + 1):
                try:
                    yield (source_partition, sink_partition), self._edge_table(source_partition,
                                                                               sink_partition,
                                                                               create_if_missing=False)
                except ValueError:
                    pass

    def _flush_regions(self):
        if self._regions_dirty:
            RegionsTable._flush_regions(self)
            self._update_partitions()

    def _flush_edges(self, silent=config.hide_progressbars, update_mappability=True):
        if self._edges_dirty:
            logger.debug("Flushing edge buffer")
            self._edge_buffer.flush()

            logger.debug("Flushing all edge tables and updating index")
            for _, edge_table in self._iter_edge_tables():
                edge_table.flush(update_index=True, log_progress=False)
            logger.debug("Done updating index")

            self._enable_edge_indexes()
            self._edges_dirty = False

            if update_mappability:
                self._update_mappability()

    def flush(self, silent=config.hide_progressbars, update_mappability=True):
        """
        Write data to file and flush buffers.

        :param silent: do not print flush progress
        :param update_mappability: After writing data, update mappability and expected values
        """
        self._flush_regions()
        self._flush_edges(silent=silent, update_mappability=update_mappability)

    def _disable_edge_indexes(self):
        logger.debug("Disabling edge indexes")
        for ix, edge_table in self._iter_edge_tables():
            logger.debug("Disabling {}".format(ix))
            if edge_table.cols.source.is_indexed:
                edge_table.cols.source.remove_index()
            if edge_table.cols.sink.is_indexed:
                edge_table.cols.sink.remove_index()
            edge_table.disable_mask_index()
            edge_table.flush()

    def _enable_edge_indexes(self):
        logger.debug("Enabling edge indexes")
        for _, edge_table in self._iter_edge_tables():
            if not edge_table.cols.source.is_indexed:
                create_col_index(edge_table.cols.source)
            if not edge_table.cols.sink.is_indexed:
                create_col_index(edge_table.cols.sink)
            edge_table.enable_mask_index()
            edge_table.flush()

    def _update_partitions(self):
        logger.debug("Updating partitions!")
        n_regions = len(self.regions)
        if n_regions == 0:
            return

        logger.debug("Partition strategy: {}".format(self._partition_strategy))
        logger.debug("Regions: {}".format(n_regions))

        partition_breaks = []
        if self._partition_strategy == 'auto':
            size = max(1000, int(n_regions / 100))
            self._partition_strategy = size

        if self._partition_strategy == 'chromosome':
            previous_chromosome = None
            for i, region in enumerate(self.regions(lazy=True)):
                if region.chromosome != previous_chromosome and previous_chromosome is not None:
                    partition_breaks.append(i)
                previous_chromosome = region.chromosome
        elif isinstance(self._partition_strategy, int) or isinstance(self._partition_strategy, np.int64):
            for i in range(self._partition_strategy, n_regions, int(self._partition_strategy)):
                partition_breaks.append(i)
        elif (isinstance(self._partition_strategy, list) or
              isinstance(self._partition_strategy, tuple)):
            partition_breaks = self._partition_strategy
        else:
            raise ValueError("{} is not a valid partition strategy!".format(self._partition_strategy))

        logger.debug("Partition breaks: {}".format(partition_breaks))
        self._partition_breaks = partition_breaks
        try:
            self.meta['partition_strategy'] = self._partition_strategy
            self.meta['partition_breaks'] = partition_breaks
        except tables.FileModeError:
            pass

    def _update_field_names(self):
        """
        Set internal object variables related to edge table field names.
        """
        edge_table = self._edge_table(0, 0)

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
        if row is None:
            self._edge_buffer.add_edge(edge)
        else:
            row['source'] = edge.source
            row['sink'] = edge.sink
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

    def add_edge_from_dict(self, edge, *args, **kwargs):
        self._edge_buffer.add_dict(edge, **kwargs)

    def add_edge_from_list(self, edge, *args, **kwargs):
        if len(edge) < 3:
            self._edge_buffer.add_weight(edge[0], edge[1], **kwargs)
        else:
            self._edge_buffer.add_weight(edge[0], edge[1], weight=edge[2], **kwargs)

    def add_edge_simple(self, source, sink, weight=None, *args, **kwargs):
        self._edge_buffer.add_weight(source, sink, weight=weight, **kwargs)

    def add_edge_from_edge(self, edge, *args, **kwargs):
        self._edge_buffer.add_edge(edge, **kwargs)

    def _add_edge_from_tuple(self, edge):
        self._edge_buffer.add_list(edge)

    def add_edges(self, edges, flush=True, *args, **kwargs):
        if self._regions_dirty:
            self._flush_regions()

        for edge in edges:
            self.add_edge(edge)

        if flush:
            self._edge_buffer.flush()
            self._enable_edge_indexes()
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

    def _edge_subset_rows(self, key=None, *args, **kwargs):
        row_regions, col_regions = self._key_to_regions(key, lazy=False)

        return self._edge_subset_rows_from_regions(
            row_regions, col_regions, *args, **kwargs
        )

    def _edge_subset_rows_from_regions(self, row_regions, col_regions, excluded_filters=0,
                                       *args, **kwargs):
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

                try:
                    edge_table = self._edge_table(i, j, create_if_missing=False)
                except ValueError:
                    continue

                # if we need to get all regions in a table, return the whole thing
                if row_covered and col_covered:
                    for row in edge_table.iterrows(excluded_filters=excluded_filters,
                                                   maskable=self):
                        yield row

                # otherwise only return the subset defined by the respective indices
                else:
                    condition = "(%d < source) & (source < %d) & (% d < sink) & (sink < %d)"
                    condition1 = condition % (row_start - 1, row_end + 1, col_start - 1, col_end + 1)
                    condition2 = condition % (col_start - 1, col_end + 1, row_start - 1, row_end + 1)

                    if row_start > col_start:
                        condition1, condition2 = condition2, condition1

                    overlap = range_overlap(row_start, row_end, col_start, col_end)

                    for edge_row in edge_table.where(condition1, excluded_filters=excluded_filters,
                                                     maskable=self):
                        yield edge_row

                    for edge_row in edge_table.where(condition2, excluded_filters=excluded_filters,
                                                     maskable=self):
                        if overlap is not None:
                            if (overlap[0] <= edge_row['source'] <= overlap[1]) and (
                                    overlap[0] <= edge_row['sink'] <= overlap[1]):
                                continue

                        yield edge_row

    def _matrix_entries(self, key, row_regions, col_regions,
                        score_field=None, *args, **kwargs):
        if score_field is None:
            score_field = self._default_score_field

        for row in self._edge_subset_rows_from_regions(row_regions, col_regions, *args, **kwargs):
            yield (row['source'], row['sink'], row[score_field])

    def _row_to_edge(self, row, lazy_edge=None, **kwargs):
        if lazy_edge is None:
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

        lazy_edge._row = row
        return lazy_edge

    def _edges_subset(self, key=None, row_regions=None, col_regions=None,
                      lazy=False, lazy_edge=None, weight_field='weight',
                      writable=False, *args, **kwargs):
        if lazy and lazy_edge is None:
            if writable:
                lazy_edge = MutableLazyEdge(None, self._regions, _weight_field=weight_field)
            else:
                lazy_edge = LazyEdge(None, self._regions, _weight_field=weight_field)
        else:
            lazy_edge = None

        excluded_filters = kwargs.get('excluded_filters', 0)

        for row in self._edge_subset_rows_from_regions(row_regions, col_regions,
                                                       excluded_filters=excluded_filters):
            yield self._row_to_edge(row, lazy_edge=lazy_edge, **kwargs)

    def _edges_iter(self, lazy=False, lazy_edge=None, weight_field='weight', *args, **kwargs):
        if lazy and lazy_edge is None:
            lazy_edge = LazyEdge(None, self._regions, _weight_field=weight_field)
        else:
            lazy_edge = None

        excluded_filters = kwargs.get('excluded_filters', 0)

        for (i, j), edge_table in self._iter_edge_tables():
            for row in edge_table.iterrows(maskable=self, excluded_filters=excluded_filters):
                yield self._row_to_edge(row, lazy_edge=lazy_edge, **kwargs)

    def edges_dict(self, *args, **kwargs):
        return self._edge_subset_rows(*args, **kwargs)

    def _edges_length(self):
        s = 0
        for _, edge_table in self._iter_edge_tables():
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
        weight_field = getattr(self, '_default_score_field', None)
        default_value = getattr(self, '_default_value', 1.)

        with RareUpdateProgressBar(max_value=len(self.edges), prefix='Mappability',
                                   silent=config.hide_progressbars) as pb:
            for i, edge in enumerate(self.edges_dict(lazy=True)):
                try:
                    if weight_field is not None:
                        weight = edge[weight_field]
                    else:
                        weight = default_value
                except KeyError:
                    weight = default_value

                if weight != 0:
                    mappable[edge['source']] = True
                    mappable[edge['sink']] = True
                pb.update(i)

        self.region_data('valid', mappable)

    def filter(self, edge_filter, queue=False, log_progress=not config.hide_progressbars):
        """
        Filter edges in this object by using a
        :class:`~fanc.general.MaskFilter`.

        :param edge_filter: Class implementing :class:`~fanc.general.MaskFilter`.
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      :func:`~RegionPairsTable.run_queued_filters`
        :param log_progress: If true, process iterating through all edges
                             will be continuously reported.
        """
        total = 0
        filtered = 0
        if not queue:
            with RareUpdateProgressBar(max_value=sum(1 for _ in self._edges),
                                       silent=not log_progress,
                                       prefix="Filter") as pb:
                for i, (_, edge_table) in enumerate(self._iter_edge_tables()):
                    stats = edge_table.filter(edge_filter, _logging=False)
                    for key, value in stats.items():
                        if key != 0:
                            filtered += stats[key]
                        total += stats[key]
                    pb.update(i)
            if log_progress:
                logger.info("Total: {}. Filtered: {}".format(total, filtered))
            self._update_mappability()
        else:
            self._queued_filters.append(edge_filter)

    def run_queued_filters(self, log_progress=not config.hide_progressbars):
        """
        Run queued filters.

        :param log_progress: If true, process iterating through all edges
                             will be continuously reported.
        """
        total = 0
        filtered = 0
        with RareUpdateProgressBar(max_value=sum(1 for _ in self._edges),
                                   silent=not log_progress,
                                   prefix="Filter") as pb:
            for i, (_, edge_table) in enumerate(self._iter_edge_tables()):
                for f in self._queued_filters:
                    edge_table.queue_filter(f)

                stats = edge_table.run_queued_filters(_logging=False)
                for key, value in stats.items():
                    if key != 0:
                        filtered += stats[key]
                    total += stats[key]
                pb.update(i)
        if log_progress:
            logger.info("Total: {}. Filtered: {}".format(total, filtered))

        self._queued_filters = []
        self._update_mappability()

    def reset_filters(self, log_progress=not config.hide_progressbars):
        with RareUpdateProgressBar(max_value=sum(1 for _ in self._edges),
                                   silent=not log_progress,
                                   prefix="Reset") as pb:
            for i, (_, edge_table) in enumerate(self._iter_edge_tables()):
                edge_table.reset_all_masks(silent=True)
                pb.update(i)
        self._update_mappability()

    def downsample(self, n, file_name=None):
        """
        Sample edges from this object.

        Sampling is always done on uncorrected Hi-C matrices.

        :param n: Sample size or reference object. If n < 1 will be interpreted as
                  a fraction of total reads in this object.
        :param file_name: Output file name for down-sampled object.
        :return: :class:`~RegionPairsTable`
        """
        logger.info("Collecting valid pairs")
        total = int(sum(e.weight for e in self.edges(lazy=True, norm=False)))

        if isinstance(n, string_types) and os.path.exists(os.path.expanduser(n)):
            with load(n) as ref:
                n = sum(e.weight for e in ref.edges(lazy=True, norm=False))
        elif isinstance(n, RegionPairsContainer):
            logger.info("Using reference Hi-C object to downsample")
            n = sum(e.weight for e in n.edges(lazy=True, norm=False))
        else:
            n = float(n)
            if n < 1:
                logger.info("Using fraction of valid pairs to downsample")
                n = int(n*total)
            else:
                logger.info("Using specific number to downsample")
                n = int(n)
        logger.info("Final n: {}/{}".format(n, total))

        logger.info("Determining random sample")
        choice = np.random.choice(int(total), size=int(n), replace=False)
        choice = sorted(choice)

        logger.info("Adding sampled pairs to new object...")
        new_pairs = self.__class__(file_name=file_name, mode='w')
        new_pairs.add_regions(self.regions, preserve_attributes=False)

        with RareUpdateProgressBar(max_value=total) as pb:
            choice_counter = 0
            pairs_counter = 0
            try:
                for edge in self.edges(lazy=True, norm=False):
                    new_weight = 0
                    for i in range(int(edge.weight)):
                        while pairs_counter == choice[choice_counter]:
                            new_weight += 1
                            choice_counter += 1

                        pairs_counter += 1

                        pb.update(pairs_counter)

                    if new_weight > 0:
                        new_pairs.add_edge_simple(edge.source, edge.sink, new_weight)
            except IndexError:
                pass

        new_pairs.flush()
        return new_pairs

    @classmethod
    def merge_region_pairs_tables(cls, pairs, check_regions_identical=True,
                                  *args, **kwargs):
        try:
            for pair in pairs:
                assert isinstance(pair, RegionPairsTable)

            # check partitions are identical
            breaks = [p._partition_breaks for p in pairs]
            for i in range(1, len(breaks)):
                assert np.array_equal(breaks[0], breaks[i])

            if check_regions_identical:
                if not RegionPairsContainer.regions_identical(pairs):
                    raise ValueError("Regions in pair objects are not identical, cannot perform merge!")
                else:
                    logger.info("Regions identical")

            logger.info("Creating merged pairs object")
            kwargs['mode'] = 'w'
            kwargs['partition_strategy'] = breaks[0]
            new_pairs = cls(*args, **kwargs)

            new_pairs.add_regions(pairs[0].regions(lazy=True))
            new_pairs._disable_edge_indexes()

            # create edge tables
            partition_pairs = set()
            for pair in pairs:
                for (source_partition, sink_partition), _ in pair._iter_edge_tables():
                    new_pairs._edge_table(source_partition, sink_partition)
                    partition_pairs.add((source_partition, sink_partition))

            logger.info("Starting fast pair merge")
            with RareUpdateProgressBar(max_value=len(partition_pairs), prefix='Merge') as pb:
                for i, (source_partition, sink_partition) in enumerate(partition_pairs):
                    edge_table = new_pairs._edge_table(source_partition, sink_partition)
                    fields = edge_table.colnames
                    new_row = edge_table.row
                    for pair in pairs:
                        if not pair._has_edge_table(source_partition, sink_partition):
                            continue
                        for row in pair._edge_table(source_partition, sink_partition).iterrows():
                            for field in fields:
                                new_row[field] = row[field]
                            new_row.append()
                    edge_table.flush()
                    pb.update(i)
            new_pairs._edges_dirty = True

            new_pairs.flush()
        except (AttributeError, AssertionError):
            raise ValueError("Partitioning is not identical, cannot "
                             "perform region pairs table merge")
        return new_pairs

    @classmethod
    def merge(cls, pairs, *args, **kwargs):
        """
        Merge two or more :class:`~RegionPairsTable` objects.

        :param pairs: list of :class:`~RegionPairsTable`
        :return: merged :class:`~RegionPairsTable`
        """
        pairs = [pair for pair in pairs]
        if not RegionPairsContainer.regions_identical(pairs):
            raise ValueError("Regions in pair objects are not identical, "
                             "cannot perform merge!")

        try:
            return cls.merge_region_pairs_tables(pairs, check_regions_identical=False, **kwargs)
        except ValueError:
            logger.info("Pair objects not compatible with fast merge, "
                        "performing regular merge")

        return RegionPairsContainer.merge(cls, pairs, *args, **kwargs)

    def subset(self, *regions, **kwargs):
        """
        Subset a Hic object by specifying one or more subset regions.

        :param regions: string or GenomicRegion object(s)
        :param kwargs: Supports
                       file_name: destination file name of subset Hic object;
                       tmpdir: if True works in tmp until object is closed
                       additional parameters are passed to
                       :func:`~RegionMatrixTable.edges`
        :return: :class:`~fanc.Hic`
        """
        logger.debug("Subsetting matrix")
        file_name = kwargs.get("file_name", None)
        tmpdir = kwargs.get('tmpdir', None)
        norm = kwargs.get('norm', False)

        new_pairs = self.__class__(file_name=file_name, mode='w', tmpdir=tmpdir)

        logger.debug("Adding subset regions")
        ix_converter = {}
        ix = 0
        new_regions = []
        for region_string in regions:
            for region in self.regions(region_string):
                ix_converter[region.ix] = ix
                ix += 1
                new_regions.append(region)
        new_pairs.add_regions(new_regions, preserve_attributes=False)
        bias_vector = [r.bias for r in new_regions]

        logger.debug("Adding subset edges")
        for i, region_string1 in enumerate(regions):
            for j in range(i, len(regions)):
                region_string2 = regions[j]
                for edge in self.edges((region_string1, region_string2), lazy=True, norm=norm, **kwargs):
                    source = ix_converter[edge.source]
                    sink = ix_converter[edge.sink]
                    new_pairs.add_edge([source, sink, edge.weight])
        new_pairs.flush(update_mappability=False)

        logger.debug("Adding subset bias vector")
        new_pairs.bias_vector(bias_vector)

        return new_pairs


class RegionMatrixTable(RegionMatrixContainer, RegionPairsTable):
    """
    HDF5 implementation of the :class:`~RegionMatrixContainer` interface.
    """

    _classid = 'REGIONMATRIXTABLE'

    def __init__(self, file_name=None, mode='a', tmpdir=None,
                 partition_strategy='auto',
                 additional_region_fields=None, additional_edge_fields=None,
                 default_score_field='weight', default_value=0.0,
                 _table_name_regions='regions', _table_name_edges='edges',
                 _table_name_expected_values='expected_values',
                 _edge_buffer_size=config.edge_buffer_size):

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
                                  partition_strategy=partition_strategy,
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

    def _flush_edges(self, **kwargs):
        if self._edges_dirty:
           self._remove_expected_values()

        RegionPairsTable._flush_edges(self, **kwargs)

    def set_biases(self, biases):
        self.region_data('bias', biases)
        self._remove_expected_values()

    def expected_values_and_marginals(self, selected_chromosome=None, norm=True,
                                      force=False, *args, **kwargs):
        group_name = 'corrected' if norm else 'uncorrected'

        if not force and self._expected_value_group is not None:
            try:
                group = self.file.get_node(self._expected_value_group, group_name)
                if selected_chromosome is not None:
                    return self.file.get_node(group, '_' + selected_chromosome)[:]
                else:
                    intra_expected = None
                    inter_expected = None
                    marginals = None
                    chromosome_intra_expected = {}
                    for node in self.file.walk_nodes(group):
                        if isinstance(node, tables.Group):
                            continue
                        if node.name == '__intra__':
                            intra_expected = node[:]
                        elif node.name == '__marginals__':
                            marginals = node[:]
                        elif node.name == '__inter__':
                            inter_expected = node[0]
                        else:
                            if node.name.startswith('_'):
                                chromosome = node.name[1:]
                                chromosome_intra_expected[chromosome] = node[:]
                    valid = self.region_data('valid')
                    if intra_expected is not None and inter_expected is not None and marginals is not None \
                            and len(chromosome_intra_expected) > 0:
                        return intra_expected, chromosome_intra_expected, inter_expected, marginals, valid
            except tables.NoSuchNodeError:
                pass

        (intra_expected, chromosome_intra_expected,
         inter_expected, marginals, valid) = RegionMatrixContainer.expected_values_and_marginals(self, norm=norm, *args,
                                                                                                 **kwargs)

        # try saving to object
        logger.debug("Attempting to save expected values and marginals to file")
        if hasattr(self, '_expected_value_group') and self._expected_value_group is not None:
            try:
                try:
                    logger.debug("Removing old expected value vectors")
                    self.file.remove_node(self._expected_value_group, group_name, recursive=True)
                except tables.NoSuchNodeError:
                    pass

                logger.debug("Creating expected value group")
                group = self.file.create_group(self._expected_value_group, group_name)

                logger.debug("Saving intra-chromosomal expected values")
                self.file.create_array(group, '__intra__',
                                       np.array(intra_expected), "Intra-chromosomal expected values")
                logger.debug("Saving inter-chromosomal expected values")
                self.file.create_array(group, '__inter__',
                                       np.array([inter_expected]), "Inter-chromosomal expected value")
                logger.debug("Saving marginals")
                self.file.create_array(group, '__marginals__',
                                       np.array(marginals), "Marginals")

                for chromosome, values in chromosome_intra_expected.items():
                    logger.debug("Saving intra-chromosomal expected values {}".format(chromosome))
                    self.file.create_array(group, '_' + chromosome,
                                           np.array(values), "Intra-chromosomal expected "
                                                             "value {}".format(chromosome))
            except tables.FileModeError:
                warnings.warn("Matrix file opened in read-only mode, "
                              "cannot save expected values to object. "
                              "Run 'fanc expected <matrix_file>' on the "
                              "command line or in Python "
                              "use mode 'a' to add expected values to "
                              "an existing object. The results of the current "
                              "computation are not affected if you don't "
                              "do this, but it will speed things up in the future.")

        try:
            self.region_data('valid', valid)
        except (OSError, KeyError):  # ignore older Hic versions and read-only files
            pass

        if selected_chromosome is not None:
            return chromosome_intra_expected[selected_chromosome], marginals

        return intra_expected, chromosome_intra_expected, inter_expected, marginals, valid

    def _update_mappability(self):
        _ = self.expected_values_and_marginals(force=True)

    def region_data(self, key, value=None):
        data = RegionPairsTable.region_data(self, key, value)

        if key == 'bias' and value is not None:
            logger.debug("Recalculating mappability and expected values after bias vector change!")
            self._update_mappability()

        return data

    @classmethod
    def merge_region_matrix_tables(cls, matrices, check_regions_identical=True,
                                   *args, **kwargs):
        try:
            for matrix in matrices:
                assert isinstance(matrix, RegionMatrixTable)

            # check partitions are identical
            breaks = [m._partition_breaks for m in matrices]
            for i in range(1, len(breaks)):
                assert np.array_equal(breaks[0], breaks[i])

            if check_regions_identical:
                if not RegionPairsContainer.regions_identical(matrices):
                    raise ValueError("Regions in matrix objects are not "
                                     "identical, cannot perform merge!")
                else:
                    logger.info("Regions identical")
        except (AttributeError, AssertionError) as e:
            logger.error(e)
            raise ValueError("Partitioning is not identical, cannot "
                             "perform region pairs table merge")

        logger.info("Creating merged matrix object")
        kwargs['mode'] = 'w'
        kwargs['partition_strategy'] = breaks[0]

        new_matrix = cls(*args, **kwargs)

        logger.info("Adding regions to merged matrix")
        new_matrix.add_regions(matrices[0].regions(lazy=True))

        logger.info("Preparing internal file structure")
        # create edge tables
        partition_pairs = set()
        for matrix in matrices:
            for (source_partition, sink_partition), _ in matrix._iter_edge_tables():
                new_matrix._edge_table(source_partition, sink_partition)
                partition_pairs.add((source_partition, sink_partition))

        new_matrix._disable_edge_indexes()

        default_field = getattr(new_matrix, '_default_score_field', 'weight')
        logger.info("Starting fast matrix merge")
        with RareUpdateProgressBar(max_value=len(partition_pairs), prefix="Merge") as pb:
            for i, (source_partition, sink_partition) in enumerate(partition_pairs):
                edges = defaultdict(int)
                for matrix in matrices:
                    if not matrix._has_edge_table(source_partition, sink_partition):
                        continue

                    for row in matrix._edge_table(source_partition, sink_partition).iterrows():
                        edges[(row['source'], row['sink'])] += row[default_field]

                edge_table = new_matrix._edge_table(source_partition, sink_partition)
                new_row = edge_table.row
                for (source, sink), weight in edges.items():
                    new_row['source'] = source
                    new_row['sink'] = sink
                    new_row[default_field] = weight
                    new_row.append()
                edge_table.flush()
                pb.update(i)
        logger.info("Done merging matrices")
        new_matrix._edges_dirty = True

        logger.debug("Flushing changes to file")
        new_matrix.flush()
        logger.debug("Done flushing changes to file")

        return new_matrix

    @classmethod
    def merge(cls, matrices, *args, **kwargs):
        """
        Merge multiple :class:`~RegionMatrixContainer` objects.

        Merging is done by adding the weight of edges in each object.

        :param matrices: list of :class:`~RegionMatrixContainer`
        :return: merged :class:`~RegionMatrixContainer`
        """
        matrices = [matrix for matrix in matrices]
        if not RegionPairsContainer.regions_identical(matrices):
            raise ValueError("Regions in matrix objects are not identical, "
                             "cannot perform merge!")

        try:
            return cls.merge_region_matrix_tables(matrices, check_regions_identical=False,
                                                  *args, **kwargs)
        except ValueError:
            logger.info("Pair objects not compatible with fast merge, "
                        "performing regular merge")

        logger.info("Adding {} regions to merged matrix".format(len(matrices[0].regions)))
        kwargs['mode'] = 'w'
        merged_matrix = cls(*args, **kwargs)
        merged_matrix.add_regions(matrices[0].regions(lazy=True), preserve_attributes=False)

        default_field = getattr(merged_matrix, '_default_score_field', 'weight')

        logger.info("Adding edges to merged matrix")
        chromosomes = merged_matrix.chromosomes()
        lc = len(chromosomes)
        with RareUpdateProgressBar(max_value=int((lc ** 2 + lc) / 2), prefix="Merge") as pb:
            chromosome_pair_ix = 0
            for i in range(len(chromosomes)):
                chromosome1 = chromosomes[i]
                for j in range(i, len(chromosomes)):
                    chromosome2 = chromosomes[j]
                    logger.debug("Adding edges from {} vs {}".format(chromosome1, chromosome2))
                    chromosome_pair_ix += 1

                    edges = defaultdict(int)
                    for matrix_object in matrices:
                        for edge in matrix_object.edges((chromosome1, chromosome2),
                                                        lazy=True, norm=False):
                            edges[edge.source, edge.sink] += getattr(edge, default_field)

                    merged_matrix.add_edges(edges, flush=False)
                    pb.update(chromosome_pair_ix)

        merged_matrix.flush()
        return merged_matrix


class RegionMatrix(np.ma.MaskedArray):
    """
    Subclass of :class:`~np.ma.masked_array` with genomic region support.

    Objects of this type are returned by :class:`~RegionMatrixContainer.matrix`.
    :class:`~RegionMatrix` supports subsetting by :class:`~fanc.GenomicRegion`
    and region strings of the form :code:`<chromosome>[:<start>-<end>]`.

    .. code::

        import fanc
        hic = fanc.load("output/hic/binned/fanc_example_1mb.hic")

        m = hic.matrix(('chr18', 'chr18'))
        type(m)  # fanc.matrix.RegionMatrix

        m_sub = m['chr18:1-5mb', 'chr18:1-10mb']
        type(m_sub)  # fanc.matrix.RegionMatrix
        m.shape  # 5, 10
        m_sub.row_regions  # [chr18:1-1000000, chr18:1000001-2000000,
                           #  chr18:2000001-3000000, chr18:3000001-4000000,
                           #  chr18:4000001-5000000]

    If the associated row or col regions have a :code:`False` :code:`valid`
    attribute, the rows/cols of the ::class:`~RegionMatrix` will be masked.

    .. attribute:: row_regions

        A list of regions matching the first matrix dimension

    .. attribute:: col_regions

        A list of regions matching the second matrix dimension
    """
    def __new__(cls, input_matrix, row_regions=None, col_regions=None,
                mask=True, *args, **kwargs):
        obj = np.asarray(input_matrix).view(cls, *args, **kwargs)
        obj._row_region_trees = None
        obj._col_region_trees = None
        obj.col_regions = None
        obj.row_regions = None
        obj.set_row_regions(row_regions)
        obj.set_col_regions(col_regions)
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
                        if 0 <= col_ix < self.shape[1]:
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

        try:
            if isinstance(row_regions, GenomicRegion):
                out.row_regions = [row_regions]
            else:
                out.row_regions = row_regions

            if isinstance(col_regions, GenomicRegion):
                out.col_regions = [col_regions]
            else:
                out.col_regions = col_regions
        except AttributeError:
            pass

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
