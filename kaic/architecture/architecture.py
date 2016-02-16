from kaic.data.genomic import RegionsTable, GenomicRegion
from abc import abstractmethod, ABCMeta
import tables as t
from collections import defaultdict

def _get_pytables_data_type(value):
    if isinstance(value, float):
        return t.Float32Col
    if isinstance(value, int):
        return t.Int32Col
    if isinstance(value, bool):
        return t.BoolCol
    if isinstance(value, str):
        return t.StringCol


class ArchitecturalFeature(object):
    __metaclass__ = ABCMeta

    @classmethod
    @abstractmethod
    def from_file(cls, file_name, *args, **kwargs):
        pass


class VectorArchitecturalRegionFeature(RegionsTable):
    def __init__(self, file_name=None, mode='a', data_fields=None,
                 regions=None, data=None, _table_name_regions='regions',
                 tmpdir=None):
        RegionsTable.__init__(self, regions=regions, file_name=file_name, mode=mode,
                              additional_fields=data_fields, tmpdir=tmpdir,
                              _table_name_regions=_table_name_regions)

        # process data
        if data is not None:
            self.add_data(data)

    @classmethod
    def from_regions_and_data(cls, regions, data, file_name=None, mode='a', tmpdir=None, data_name='data'):
        if not isinstance(data, dict):
            # assume this is a vector
            data = {
                data_name: data
            }

        data_fields = dict()
        for data_name, vector in data.iteritems():
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

        for data_name, vector in data.iteritems():
            self.data(data_name, vector)

    def __getitem__(self, item):
        if isinstance(item, tuple):
            return self._get_columns(item[1], regions=self._get_rows(item[0]))
        else:
            return self._get_rows(item)

    def _get_rows(self, item, lazy=False, auto_update=True):
        if isinstance(item, int):
            return self._row_to_region(self._regions[item], lazy=lazy, auto_update=auto_update)

        if isinstance(item, slice):
            return (self._row_to_region(row, lazy=lazy, auto_update=auto_update)
                    for row in self._regions.iterrows(item.start, item.stop, item.step))

        if isinstance(item, str):
            item = GenomicRegion.from_string(item)

        if isinstance(item, GenomicRegion):
            return self.subset(item, lazy=lazy, auto_update=auto_update)

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
                    results_dict[name].append(getattr(region, name, None))
        else:
            results_dict = dict()
            for name in colnames:
                results_dict[name] = getattr(regions, name, None)

        if is_list:
            return results_dict
        return results_dict[colnames[0]]



class HicArchitecture(object):
    def __init__(self, hic):
        self.hic = hic

    def distance_decay(self, per_chromosome=False):
        pass
