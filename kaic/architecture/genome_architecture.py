from __future__ import division
from kaic.data.genomic import RegionMatrixTable, RegionsTable, GenomicRegion
from kaic.architecture.architecture import ArchitecturalFeature, calculateondemand, _get_pytables_data_type
from kaic.data.general import MaskFilter
import tables as t
import numpy as np
import itertools as it
from collections import defaultdict
from abc import abstractmethod, ABCMeta
from future.utils import with_metaclass
import logging
logger = logging.getLogger(__name__)


class MatrixArchitecturalRegionFeature(RegionMatrixTable, ArchitecturalFeature):
    """
    Process and store matrix-based, genomic region-associated data.

    Behaves similarly to :class:`~Hic`, as they both inherit from
    :class:`~RegionMatrixTable`.
    """

    _classid = 'MATRIXARCHITECTURALREGIONFEATURE'

    def __init__(self, file_name=None, mode='a', data_fields=None, default_field=None,
                 regions=None, edges=None, _table_name_regions='region_data',
                 _table_name_edges='edges', tmpdir=None):
        RegionMatrixTable.__init__(self, file_name=file_name, additional_fields=data_fields,
                                   mode=mode, tmpdir=tmpdir, default_field=default_field,
                                   _table_name_nodes=_table_name_regions,
                                   _table_name_edges=_table_name_edges)
        ArchitecturalFeature.__init__(self)

        if len(self._edges) > 0:
            self._calculated = True

        if regions is not None:
            self.add_regions(regions)

        # process data
        if edges is not None:
            self.add_edges(edges)

    @calculateondemand
    def as_matrix(self, key=slice(0, None, None), values_from=None, mask_missing=False, impute_missing=False):
        """
        See :class:`~RegionMatrixTable`
        """
        return RegionMatrixTable.as_matrix(self, key=key, values_from=values_from,
                                           mask_missing=mask_missing, impute_missing=impute_missing)

    @calculateondemand
    def _get_nodes_from_key(self, key, as_index=False):
        return RegionMatrixTable._get_nodes_from_key(self, key, as_index=as_index)

    @calculateondemand
    def _get_matrix(self, row_ranges=None, col_ranges=None, weight_column=None):
        return RegionMatrixTable._get_matrix(self, row_ranges=row_ranges, col_ranges=col_ranges,
                                             weight_column=weight_column)

    @calculateondemand
    def _getitem_nodes(self, key, as_index=False):
        return RegionMatrixTable._getitem_nodes(self, key, as_index=as_index)

    @calculateondemand
    def as_data_frame(self, key, weight_column=None):
        """
        See :class:`~RegionMatrixTable`
        """
        return RegionMatrixTable.as_data_frame(self, key, weight_column=weight_column)

    @calculateondemand
    def get_node(self, key):
        """
        See :class:`~RegionMatrixTable`
        """
        return RegionMatrixTable.get_node(self, key)

    @calculateondemand
    def get_edge(self, ix, lazy=False):
        """
        See :class:`~RegionMatrixTable`
        """
        return RegionMatrixTable.get_edge(self, ix, lazy=lazy)

    @calculateondemand
    def _nodes_iter(self):
        return RegionMatrixTable._nodes_iter(self)

    @calculateondemand
    def edges_sorted(self, sortby, *args, **kwargs):
        """
        See :class:`~RegionMatrixTable`
        """
        return RegionMatrixTable.edges_sorted(self, sortby, *args, **kwargs)

    @calculateondemand
    def _edges_iter(self):
        return RegionMatrixTable._edges_iter(self)

    @abstractmethod
    def _calculate(self, *args, **kwargs):
        raise NotImplementedError("This method must be overridden in subclass!")

    @calculateondemand
    def filter(self, edge_filter, queue=False, log_progress=False):
        """
        Filter edges in this object by using a
        :class:`~MatrixArchitecturalRegionFeatureFilter`.

        :param edge_filter: Class implementing :class:`~MatrixArchitecturalRegionFeatureFilter`.
                            Must override valid_edge method, ideally sets mask parameter
                            during initialization.
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        :param log_progress: If true, process iterating through all edges
                             will be continuously reported.
        """
        edge_filter.set_matrix_object(self)
        if not queue:
            self._edges.filter(edge_filter, _logging=log_progress)
        else:
            self._edges.queue_filter(edge_filter)

    def run_queued_filters(self, log_progress=False):
        """
        Run queued filters.

        :param log_progress: If true, process iterating through all edges
                             will be continuously reported.
        """
        self._edges.run_queued_filters(_logging=log_progress)


class MatrixArchitecturalRegionFeatureFilter(with_metaclass(ABCMeta, MaskFilter)):
    """
    Abstract class that provides filtering functionality for the
    edges/contacts in a :class:`~MatrixArchitecturalRegionFeature` object.

    Extends MaskFilter and overrides valid(self, row) to make
    :class:`~Edge` filtering more "natural".

    To create custom filters for the :class:`~MatrixArchitecturalRegionFeature`
    object, extend this class and override the valid_edge(self, edge) method.
    valid_edge should return False for a specific :class:`~Edge` object
    if the object is supposed to be filtered/masked and True
    otherwise.
    """

    def __init__(self, mask=None):
        """
        Initialize MatrixArchitecturalRegionFeatureFilter.

        :param mask: The Mask object that should be used to mask
                     filtered :class:`~Edge` objects. If None the default
                     Mask will be used.
        """
        super(MatrixArchitecturalRegionFeatureFilter, self).__init__(mask)
        self._hic = None

    @abstractmethod
    def valid_edge(self, edge):
        """
        Determine if an :class:`~Edge` object is valid or should
        be filtered.

        When implementing custom MatrixArchitecturalRegionFeatureFilter this
        method must be overridden. It should return False for :class:`~HicEdge`
        objects that are to be fitered and True otherwise.

        Internally, the :class:`~MatrixArchitecturalRegionFeature` object
        will iterate over all Edge instances to determine their validity on
        an individual basis.

        :param edge: A :class:`~Edge` object
        :return: True if :class:`~Edge` is valid, False otherwise
        """
        pass

    def set_matrix_object(self, matrix_object):
        """
        Set the :class:`~MatrixArchitecturalRegionFeature` instance to
        be filtered by this MatrixArchitecturalRegionFeatureFilter.

        Used internally by :class:`~MatrixArchitecturalRegionFeature`
        instance.

        :param matrix_object: :class:`~MatrixArchitecturalRegionFeature`
                              object
        """
        self._matrix = matrix_object

    def valid(self, row):
        """
        Map valid_edge to MaskFilter.valid(self, row).

        :param row: A pytables Table row.
        :return: The boolean value returned by valid_edge.
        """
        edge = self._matrix._row_to_edge(row, lazy=True)
        return self.valid_edge(edge)


class VectorArchitecturalRegionFeature(RegionsTable, ArchitecturalFeature):
    """
    Process and store vector data associated with genomic regions.
    """

    _classid = 'VECTORARCHITECTURALREGIONFEATURE'

    def __init__(self, file_name=None, mode='a', data_fields=None,
                 regions=None, data=None, _table_name_data='region_data',
                 tmpdir=None):
        """

        :param file_name: Path to file to save/retrieve feature information
        :param mode: File mode ('r' = read-only, 'w' = (over)write, 'a' = append)
        :param data_fields: Dictionary of class with PyTables column definitions.
        :param regions: Any iterable with objects accepted by
                        :func:`~kaic.data.genomic.GenomicRegions.add_region`.
        :param data: Dict of iterables with vector data to import
        :param tmpdir: Path to temporary directory.
        """
        RegionsTable.__init__(self, file_name=file_name, mode=mode,
                              additional_fields=data_fields, tmpdir=tmpdir,
                              _table_name_regions=_table_name_data)
        ArchitecturalFeature.__init__(self)

        if len(self._regions) > 0:
            self._calculated = True

        if regions is not None:
            self.add_regions(regions)

        # process data
        if data is not None:
            self.add_data(data)

    def flush(self):
        self._regions.flush()

    @property
    def data_field_names(self):
        names = []
        for name in self._regions.colnames:
            if name not in ("ix", "chromosome", "start", "end", "strand"):
                names.append(name)
        return names

    @classmethod
    def from_regions_and_data(cls, regions, data, file_name=None, mode='a', tmpdir=None, data_name='data'):
        if not isinstance(data, dict):
            # assume this is a vector
            data = {
                data_name: data
            }

        data_fields = dict()
        for data_name, vector in data.items():
            string_size = 0
            for value in vector:
                table_type = _get_pytables_data_type(value)
                if table_type != t.StringCol:
                    data_fields[data_name] = table_type(pos=len(data_fields))
                    break
                else:
                    string_size = max(string_size, len(value))
            if string_size > 0:
                data_fields[data_name] = t.StringCol(string_size, pos=len(data_fields))

        self = cls(file_name=file_name, mode=mode, data_fields=data_fields,
                   regions=regions, data=data, tmpdir=tmpdir)
        return self

    def add_data(self, data, name="data"):
        """
        Add vector-data to this object. If there is exsting data in this
        object with the same name, it will be replaced

        :param data: Either an iterable with the same number of items as
                     regions in this object, or a dictionary of iterables
                     if multiple objects should be imported
        :param name: (optional) name of the data set if data is a single
                     iterable
        """

        if not isinstance(data, dict):
            # assume this is a vector
            data = {
                name: data
            }

        for data_name, vector in data.items():
            self.data(data_name, vector)

    @calculateondemand
    def __setitem__(self, key, value):
        if isinstance(key, tuple):
            n_rows = sum(1 for _ in self._get_rows(key[0], lazy=True,
                                                   auto_update=False, force_iterable=True))
            rows = self._get_rows(key[0], lazy=True, auto_update=False, force_iterable=True)
            column_selectors = key[1]
        else:
            n_rows = sum(1 for _ in self._get_rows(slice(0, None, None), lazy=True,
                                                   auto_update=False, force_iterable=True))
            rows = self._get_rows(slice(0, None, None), lazy=True, auto_update=False, force_iterable=True)
            column_selectors = key

        if not isinstance(column_selectors, list):
            column_selectors = [column_selectors]

        colnames = []
        for column_selector in column_selectors:
            if isinstance(column_selector, int) or isinstance(column_selector, slice):
                colnames.append(self._regions.colnames[column_selector])
            elif isinstance(column_selector, str):
                colnames.append(column_selector)

        if isinstance(value, list):
            if len(value) != n_rows:
                raise ValueError("Number of elements in selection does not "
                                 "match number of elements to be replaced!")
            for i, row in enumerate(rows):
                value_row = value[i]
                for j, colname in enumerate(colnames):
                    try:
                        v = getattr(value_row, colname)
                    except AttributeError:
                        if len(colnames) == 1:
                            v = value_row
                        else:
                            try:
                                v = value_row[j]
                            except TypeError:
                                raise ValueError("Bad value format")
                            except KeyError:
                                v = value[colname]
                    setattr(row, colname, v)
                    row.update()
        else:
            if n_rows != 1:
                raise ValueError("Can only replace selection with elements in a list")

            for row in rows:
                setattr(row, colnames[0], value)
                row.update()
        self.flush()

    def __getitem__(self, item):
        if isinstance(item, tuple):
            return self._get_columns(item[1], regions=self._get_rows(item[0]))
        else:
            return self._get_rows(item)

    @calculateondemand
    def _get_rows(self, item, lazy=False, auto_update=True, force_iterable=False):
        if isinstance(item, int):
            if force_iterable:
                return (self._row_to_region(row, lazy=lazy, auto_update=auto_update)
                        for row in self._regions.iterrows(item, item+1, 1))
            return self._row_to_region(self._regions[item], lazy=lazy, auto_update=auto_update)

        if isinstance(item, slice):
            return (self._row_to_region(row, lazy=lazy, auto_update=auto_update)
                    for row in self._regions.iterrows(item.start, item.stop, item.step))

        if isinstance(item, str):
            item = GenomicRegion.from_string(item)

        if isinstance(item, GenomicRegion):
            return self.subset(item, lazy=lazy, auto_update=auto_update)

    @calculateondemand
    def _get_columns(self, item, regions=None):
        if regions is None:
            regions = self._get_rows(slice(0, None, None), lazy=True)

        is_list = True
        if isinstance(item, int):
            colnames = [self._regions.colnames[item]]
            is_list = False
        elif isinstance(item, slice):
            colnames = self._regions.colnames[item]
        elif isinstance(item, str):
            colnames = [item]
            is_list = False
        elif isinstance(item, list):
            colnames = item
        else:
            raise KeyError("Unrecognised key type (%s)" % str(type(item)))

        if not isinstance(regions, GenomicRegion):
            results_dict = defaultdict(list)
            for region in regions:
                for name in colnames:
                    results_dict[name].append(getattr(region, name))
        else:
            results_dict = dict()
            for name in colnames:
                results_dict[name] = getattr(regions, name)

        if is_list:
            return results_dict
        return results_dict[colnames[0]]

    @abstractmethod
    def _calculate(self, *args, **kwargs):
        raise NotImplementedError("This method must be overridden in subclass!")


class BasicRegionTable(VectorArchitecturalRegionFeature):
    """
    Convenience class to store region-based meta information.
    Sample usage:

    .. code: python
       # 10 regions
       regions = [GenomicRegion(chromosome='chr1', start=1+i, end=i+1000) for i in xrange(0, 10000, 1000)]

       # string data
       s = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']
       # int data
       i = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

       t = BasicRegionTable(regions, fields=('mystr', 'myint'), types=(str, int), data={'mystr': s, 'myint': i})
    """
    _classid = 'VECTORARCHITECTURALREGIONFEATURE'

    def __init__(self, regions, fields=None, types=None, data=None,
                 file_name=None, mode='a', tmpdir=None,
                 _string_size=100, _group_name='region_table'):
        """
        Initialize a BasicRegionTable object.
        :param regions: Iterable with :class:`~GenomicRegion` objects.
        :param fields: Either dictionary of class with PyTables column definitions or
                       list of field names. If the latter, must provide a list with
                       types of equal length.
        :param types: List of Python or PyTables data types for table columns.
        :param data: Dict of iterables with vector data to import
        :param file_name: Path to file to save/retrieve feature information
        :param mode: File mode ('r' = read-only, 'w' = (over)write, 'a' = append)
        :param tmpdir: Path to temporary directory.
        """
        if isinstance(regions, str):
            if file_name is None:
                file_name = regions
                regions = None
            else:
                raise ValueError("fields cannot be string unless file_name is None")

        if regions is None and file_name is not None:
            VectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode,
                                                      _table_name_data=_group_name, tmpdir=tmpdir)
        else:
            pt_fields = {}
            if fields is not None:
                if isinstance(fields, dict):
                    for field, field_type in fields.items():
                        pt_fields[field] = _get_pytables_data_type(field_type)
                else:
                    if types is None or not len(fields) == len(types):
                        raise ValueError("fields (%d) must be the same length as types (%d)" % (len(fields),
                                                                                                len(types)))
                    for i, field in enumerate(fields):
                        pt_fields[field] = _get_pytables_data_type(types[i])

            data_fields = {}
            for data_name, table_type in pt_fields.items():
                if table_type != t.StringCol:
                    data_fields[data_name] = table_type(pos=len(data_fields))
                else:
                    data_fields[data_name] = table_type(_string_size, pos=len(data_fields))

            VectorArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode, data_fields=data_fields,
                                                      regions=regions, data=data, _table_name_data=_group_name,
                                                      tmpdir=tmpdir)

    def _calculate(self, *args, **kwargs):
        pass


def _get_typed_array(input_iterable, nan_strings, count=-1):
    try:
        return np.fromiter((0 if x in nan_strings else x for x in input_iterable), int, count)
    except ValueError:
        pass
    try:
        return np.fromiter((np.nan if x in nan_strings else x for x in input_iterable), float, count)
    except ValueError:
        pass
    return np.fromiter(input_iterable, str, count)


def _is_simple_type(data_type):
    if data_type in {int, float, bool, str, long, unicode}:
        return True
    return False


class GenomicTrack(BasicRegionTable):
    """
    Update of legacy GenomicTrack class backed by PyTables table.
    """
    _classid = 'GENOMICTRACK'

    def __init__(self, file_name=None, title=None, data_dict=None, regions=None, _table_name_tracks='tracks',
                 mode='a', tmpdir=None):
        """
        Initialize a genomic track vector.

        :param file_name: Storage location of the genomic track HDF5 file
        :param data_dict: Dictionary containing data tracks as numpy arrays.
                          The arrays must have as many elements in the first
                          dimension as there are regions.
        :param title: The overall title of the track.
        :param regions: An iterable of (:class: `~kaic.data.genomic.GenomicRegion`)
                        or String elemnts that describe regions.
        """
        matrix_data = {}
        fields = {}
        if data_dict is not None:
            for field, values in data_dict.items():
                data_type = type(values[0])
                if _is_simple_type(data_type):
                    fields[field] = type(values[0])
                else:
                    matrix_data[field] = values

        for key in matrix_data.keys():
            del data_dict[key]

        self._matrix_tracks = set()
        BasicRegionTable.__init__(self, regions=regions, fields=fields, data=data_dict, file_name=file_name,
                                  _group_name=_table_name_tracks, mode=mode, tmpdir=tmpdir)

        for key, values in matrix_data.items():
            self.data(key, values)

        if title is not None:
            self.title = title

        for node in self.file.iter_nodes(self._group):
            if node.name != 'regions' and node.name not in self._matrix_tracks:
                self._matrix_tracks.add(node.name)

    @property
    def title(self):
        try:
            return self._regions.attrs['title']
        except KeyError:
            return None

    @title.setter
    def title(self, title):
        self._regions.attrs['title'] = title

    @property
    def _tracks(self):
        return self.data_field_names

    @classmethod
    def from_gtf(cls, file_name, gtf_file, store_attrs=None, nan_strings=(".", ""), sort=False):
        """
        Import a GTF file as GenomicTrack.

        :param file_name: Storage location of the genomic track HDF5 file
        :param gtf_file: Location of GTF file_name
        :param store_attrs: List or listlike
                            Only store attributes in the list
        :param nan_strings: These characters will be considered NaN for parsing.
                            Will become 0 for int arrays, np.nan for float arrays
                            and left as is for string arrays.
        :param sort: If True, sort regions in GTF file by genomic location
        """
        import pybedtools as pbt
        gtf = pbt.BedTool(gtf_file)
        n = len(gtf)
        regions = []
        values = {}
        if not sort:
            gtf_iter = gtf
        else:
            gtf_iter = gtf.sort()
        for i, f in enumerate(gtf_iter):
            regions.append(GenomicRegion(chromosome=f.chrom, start=f.start, end=f.end, strand=f.strand))
            # If input is a GTF file, also store the type and source fields
            if f.file_type in ("gff", "gtf"):
                f.attrs["source"] = f.fields[1]
                f.attrs["feature"] = f.fields[2]
            # Check if there is a new attribute that hasn't occured before
            for k in f.attrs.keys():
                if k not in values and (not store_attrs or k in store_attrs):
                    if i > 0:
                        # Fill up values for this attribute with nan
                        values[k] = [nan_strings[0]]*i
                    else:
                        values[k] = []
            for k in values.keys():
                values[k].append(f.attrs.get(k, nan_strings[0]))
        for k, v in values.items():
            values[k] = _get_typed_array(v, nan_strings=nan_strings, count=n)
        return cls(file_name=file_name, data_dict=values, regions=regions)

    def to_bedgraph(self, prefix, tracks=None, skip_nan=True):
        """
        Exoprt data as BEDgraph file.

        :param prefix: Prefix of file - entire path in the form /path/to/prefix
        :param tracks: Name of track/field to save
        :param skip_nan: If True, NaN values will not be written to file
        """
        if tracks is None:
            tracks = self._tracks
        elif not isinstance(tracks, list) and not isinstance(tracks, tuple):
            tracks = [tracks]
        for track in tracks:
            logger.info("Writing track {}".format(track))
            with open("{}{}.bedgraph".format(prefix, track), "w") as f:
                for r, v in it.izip(self.regions, self[track]):
                    if skip_nan and np.isnan(v):
                        continue
                    f.write("{}\t{}\t{}\t{}\n".format(r.chromosome, r.start - 1, r.end, v))

    def __getitem__(self, item):
        if isinstance(item, str):
            if item in self._tracks:
                return self[:, item]
            if item in self._matrix_tracks:
                return self.data(item)[:]
        elif isinstance(item, int) or isinstance(item, slice):
            tracks = dict()
            for name in self._tracks:
                tracks[name] = self[item, name]
            for name in self._matrix_tracks:
                tracks[name] = np.array(getattr(self._group, name)[item])
            return tracks
        else:
            return BasicRegionTable.__getitem__(self, item)

    @property
    def tracks(self):
        """
        Get a dictionary with all tracks/vectors in the data set.
        :return: dict
        """
        return {t_name: self[:, t_name] for t_name in self._tracks}

    def data(self, key, value=None):
        """
        Retrieve or add vector-data OR matrix data to this object. If there is exsting data in this
        object with the same name, it will be replaced

        :param key: Name of the data column
        :param value: vector with region-based data (one entry per region)
        """
        # Hack for matrix-based data
        if value is not None and not _is_simple_type(type(value[0])):
            a = np.array(value)
            if a.shape[0] != len(self.regions):
                raise ValueError("First dimension of values must have as many elements "
                                 "({}) as there are regions ({})".format(value.shape, len(self._regions)))
            self.file.create_array(self._group, key, value)
            self._matrix_tracks.add(key)

        if key in self._matrix_tracks:
            return getattr(self._group, key)[:]

        return RegionsTable.data(self, key, value=value)
