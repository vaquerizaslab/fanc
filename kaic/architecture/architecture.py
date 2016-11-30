from kaic.data.general import FileGroup
from abc import abstractmethod
import tables as t
from collections import defaultdict
import inspect


def _get_pytables_data_type(value):
    """
    Convert a value or class to a PyTables data type.
    :param value: Primitive value
    :return: PyTables Col object
    """

    if inspect.isclass(value):
        compare_method = issubclass
        if compare_method(value, t.Col):
            return value
    else:
        compare_method = isinstance
        if compare_method(value, t.Col):
            return value.__class__

    if compare_method(value, float):
        return t.Float32Col
    if compare_method(value, int):
        return t.Int32Col
    if compare_method(value, bool):
        return t.BoolCol
    if compare_method(value, str):
        return t.StringCol


def calculateondemand(func):
    """
    Decorator: run object calculations before executing method.

    :param func: An :class:`~ArchitecturalFeature` object's method
    :return: The original func's return value
    """
    def inner(self, *args, **kwargs):
        if 'force' in kwargs and isinstance(kwargs['force'], bool):
            force = kwargs['force']
            del kwargs['force']
        else:
            force = False
        if force or not self._calculated:
            self.calculate()

        return func(self, *args, **kwargs)
    return inner


class ArchitecturalFeature(object):
    """
    Base class for architectural features.

    Subclasses must override _calculate method to provide full functionality.
    _calculate should save calculated results in the object.

    Subclass methods decorated with :func:`calculateondemand` will trigger
    the execution of :func:`ArchitecturalFeature::calculate` on first run.
    If :func:`ArchitecturalFeature::calculate` has already been run before,
    the method is called directly (accessing saved data).
    """
    def __init__(self):
        self._calculated = False

    def calculate(self, *args, **kwargs):
        """
        Run _calculate of this object and set _calculated status to True.
        """
        self._calculated = True
        self._calculate(*args, **kwargs)

    def _calculate(self, *args, **kwargs):
        pass


class TableArchitecturalFeature(FileGroup, ArchitecturalFeature):
    """
    Calculate and store tabular data.

    Provides column name access to data, vector assignment, and other
    useful features.

    Examples:

    .. code:: python

        # data assignment
        taf.data('int_data', [1, 2, 3, 4, 5])
        taf.data('text_data', ['a', 'b', 'c', 'd', 'e'])

        # data retrieval
        taf.data('int_data')  # [1, 2, 3, 4, 5]
        taf.data('text_data')  # ['a', 'b', 'c', 'd', 'e']

        taf[1:4, 'a']  # [2, 3, 4]
    """

    _classid = 'TABLEARCHITECTURALFEATURE'

    def __init__(self, group, fields=None, file_name=None, mode='a',
                 tmpdir=None, _table_name='table_architecture'):
        if isinstance(fields, str):
            if file_name is None:
                file_name = fields
            else:
                raise ValueError("Must provide fields if creating new table!")

        FileGroup.__init__(self, group, file_name=file_name, mode=mode, tmpdir=tmpdir)
        ArchitecturalFeature.__init__(self)

        try:
            self._table = getattr(self._group, _table_name)
            # there is data in the table, no need to recalculate
            if len(self._table) > 0:
                self._calculated = True
        except t.NoSuchNodeError:
            self._table = t.Table(self._group, _table_name, fields)

    def flush(self):
        self._table.flush()

    @classmethod
    def from_data(cls, group, data, file_name=None, mode='a', tmpdir=None, data_name='data'):
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

        self = cls(group, fields=data_fields, file_name=file_name, mode=mode, tmpdir=tmpdir)
        self.add_data(data)
        return self

    def data(self, key, value=None):
        """
        Retrieve or add vector-data to this object. If there is exsting data in this
        object with the same name, it will be replaced. If the value vector is longer
        than the existing table, the table will be expanded to fit.

        :param key: Name of the data column
        :param value: vector with region-based data (one entry per region)
        """
        if key not in self._table.colnames:
            raise KeyError("%s is unknown region attribute" % key)

        if value is not None:
            for i, row in enumerate(self._table):
                row[key] = value[i]
                row.update()
            self._table.flush()

            row = self._table.row
            for i in range(len(self._table), len(value)):
                row[key] = value[i]
                row.append()
            self._table.flush()

        return (row[key] for row in self._table)

    def add_data(self, data, data_name="data"):
        """
        Add vector-data to this object. If there is exsting data in this
        object with the same name, it will be replaced. If the value vector is longer
        than the existing table, the table will be expanded to fit.

        :param data: Either an iterable with the same number of items as
                     regions in this object, or a dictionary of iterables
                     if multiple objects should be imported
        :param data_name: (optional) name of the data set if data is a single
                          iterable
        """

        if not isinstance(data, dict):
            # assume this is a vector
            data = {
                data_name: data
            }

        for data_name, vector in data.items():
            self.data(data_name, vector)

    @calculateondemand
    def __setitem__(self, key, value):
        if isinstance(key, tuple):
            n_rows = sum(1 for _ in self._get_rows(key[0], lazy=True, as_row=True))
            rows = self._get_rows(key[0], lazy=True, as_row=True)
            column_selectors = key[1]
        else:
            n_rows = sum(1 for _ in self._get_rows(slice(0, None, None), lazy=True, as_row=True))
            rows = self._get_rows(slice(0, None, None), lazy=True, as_row=True)
            column_selectors = key

        if not isinstance(column_selectors, list):
            column_selectors = [column_selectors]

        colnames = []
        for column_selector in column_selectors:
            if isinstance(column_selector, int) or isinstance(column_selector, slice):
                colnames.append(self._table.colnames[column_selector])
            elif isinstance(column_selector, str):
                colnames.append(column_selector)

        if isinstance(value, list):
            print(len(value), n_rows)
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
                    row[colname] = v
                    row.update()

        else:
            if n_rows != 1:
                raise ValueError("Can only replace selection with elements in a list")

            for row in rows:
                row[colnames[0]] = value
                row.update()
        self.flush()

    def __getitem__(self, item):
        if isinstance(item, tuple):
            return self._get_columns(item[1], rows=self._get_rows(item[0]))
        else:
            return self._get_rows(item)

    def _row_to_entry(self, row, lazy=False, auto_update=True):
        if lazy:
            return row
        d = dict()
        for colname in self._table.colnames:
            d[colname] = row[colname]
        return d

    @calculateondemand
    def _get_rows(self, item, lazy=False, auto_update=True, as_row=False):
        if isinstance(item, int):
            if as_row:
                iterator = self._table.iterrows(item, item+1, 1)
                return iterator
            return self._row_to_entry(self._table[item], lazy=lazy, auto_update=auto_update)

        if isinstance(item, slice):
            iterator = self._table.iterrows(item.start, item.stop, item.step)
            if as_row:
                return iterator
            return (self._row_to_entry(row, lazy=lazy, auto_update=auto_update)
                    for row in iterator)

    @calculateondemand
    def _get_columns(self, item, rows=None):
        if rows is None:
            rows = self._get_rows(slice(0, None, None), lazy=True)

        is_list = True
        if isinstance(item, int):
            colnames = [self._table.colnames[item]]
            is_list = False
        elif isinstance(item, slice):
            colnames = self._table.colnames[item]
        elif isinstance(item, str):
            colnames = [item]
            is_list = False
        elif isinstance(item, list):
            colnames = item
        else:
            raise KeyError("Unrecognised key type (%s)" % str(type(item)))

        if not isinstance(rows, dict):
            results_dict = defaultdict(list)
            for row in rows:
                for name in colnames:
                    results_dict[name].append(row[name])
        else:
            results_dict = dict()
            for name in colnames:
                results_dict[name] = rows[name]

        if is_list:
            return results_dict
        return results_dict[colnames[0]]

    @abstractmethod
    def _calculate(self, *args, **kwargs):
        raise NotImplementedError("This method must be overridden in subclass!")


class BasicTable(TableArchitecturalFeature):
    """
    Convenience class to store row-based information.
    """
    _classid = 'BASICTABLE'

    def __init__(self, fields, types=None, file_name=None, mode='a',
                 _string_size=100, _group_name='basic_table'):
        if isinstance(fields, str):
            if file_name is None:
                file_name = fields
                fields = None
            else:
                raise ValueError("fields cannot be string unless file_name is None")

        if fields is None:
            TableArchitecturalFeature.__init__(self, _group_name, file_name, mode=mode)
        else:
            pt_fields = {}
            if isinstance(fields, dict):
                for field, field_type in fields.items():
                    pt_fields[field] = _get_pytables_data_type(field_type)
            else:
                if types is None or len(fields) != len(types):
                    raise ValueError("fields (%d) must be the same length as types (%d)" % (len(fields), len(types)))
                for i, field in enumerate(fields):
                    pt_fields[field] = _get_pytables_data_type(types[i])

            data_fields = {}
            for data_name, table_type in pt_fields.items():
                if table_type != t.StringCol:
                    data_fields[data_name] = table_type(pos=len(data_fields))
                else:
                    data_fields[data_name] = table_type(_string_size, pos=len(data_fields))

            TableArchitecturalFeature.__init__(self, _group_name, fields=data_fields,
                                               file_name=file_name, mode=mode)

    def _calculate(self, *args, **kwargs):
        pass

