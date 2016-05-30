"""
Provide basic convenience data types.

Most classes are based on pytables and hdf5 dictionaries,
allowing on-disk storage and therefore processing of large
files. Other features include indexing and querying.
"""

from __future__ import division
import tables as t
import kaic.fixes.pytables_nrowsinbuf_inheritance_fix
from kaic.tools.files import create_or_open_pytables_file, is_hdf5_file
import numpy as np
import warnings
from kaic.tools.general import RareUpdateProgressBar
import os
import time
import logging
from tables.exceptions import NoSuchNodeError
from abc import ABCMeta, abstractmethod
from datetime import datetime
from kaic.tools.lru import lru_cache
import shutil
import binascii
from collections import defaultdict
logging.basicConfig(level=logging.INFO)
_filter = t.Filters(complib="blosc", complevel=2, shuffle=True)


def _convert_to_tables_type(col_type, pos=None):
    """
    Convert Python and numpy data types to pytables Col class.
    """
    
    _typemap = {
        str: t.StringCol, # @UndefinedVariable
        np.string_: t.StringCol, # @UndefinedVariable
        np.string0: t.StringCol, # @UndefinedVariable
        np.str_: t.StringCol, # @UndefinedVariable
        int: t.Int32Col, # @UndefinedVariable
        np.int8: t.Int32Col, # @UndefinedVariable
        np.int16: t.Int32Col, # @UndefinedVariable
        np.int32: t.Int32Col, # @UndefinedVariable
        np.int64: t.Int64Col, # @UndefinedVariable
        np.uint8: t.Int32Col, # @UndefinedVariable
        np.uint16: t.Int32Col, # @UndefinedVariable
        np.uint32: t.Int32Col, # @UndefinedVariable
        np.uint64: t.Int64Col, # @UndefinedVariable
        float: t.Float32Col, # @UndefinedVariable
        np.float16: t.Float32Col, # @UndefinedVariable
        np.float32: t.Float32Col, # @UndefinedVariable
        np.float64: t.Float64Col, # @UndefinedVariable
        np.bool: t.BoolCol, # @UndefinedVariable
        np.bool8: t.BoolCol, # @UndefinedVariable
        np.bool_: t.BoolCol # @UndefinedVariable
    }
    
    if isinstance(col_type,t.Col):
        return col_type
    if col_type in _typemap:
        if (col_type is str or col_type is np.string_ or
            col_type is np.string0 or col_type is np.str_):
            return _typemap[col_type](255,pos=pos)
        return _typemap[col_type](pos=pos)
    raise ValueError("Unknown column type " + str(col_type))


def _structured_array_to_table_type(a, rownames=None, colnames=None, return_type=None):
    """
    Convert a record/structured numpy array to TableRow, TableCol, or TableArray.
    
    Args:
        a (numpy record array): array to be converted
        rownames (list): a list of strings to be used as rownames
        colnames (list): a list of strings to be used as colnames
        
    Returns:
        TableRow: if record array has exactly one row
        TableCol: if record array has exactly one column
        TableArray: if record array has at least 2 rows and columns
                    or if a is None
    """
    
    # is row or single value
    if len(a) == 1:
        if len(a[0]) == 1:
            return a[0][0]
        if return_type is None:
            return TableRow(a[0], rowname=rownames[0], colnames=colnames)
        else:
            return return_type.from_row(a[0])

    # is column
    if len(a) > 1 and len(a[0]) == 1:
        l = []
        for row in a:
            l.append(row[0])
        return TableCol(np.array(l), colname=colnames[0], rownames=rownames)

    # is matrix or nothing
    if return_type is None:
        return TableArray(a, colnames=colnames, rownames=rownames)
    else:
        l = []
        for row in a:
            l.append(return_type.from_row(row))
        return l
    

def _to_list_and_names(a):
    """
    Convert all sorts of input data into suitable Table input. 
    
    Args:
        a (array, list, dictionary): data that can be expressed as
                                     a table with column and row names
    
    Returns:
        list: a list containing the rows of a table as elements, either
              as lists or dictionaries
        list: a list of columns names
        list: a list of column types as python data types, such as str, int
    """
    
    l = []
    colnames = None
    coltypes = None

    # structured/record array
    try:
        a.dtype.fields.keys()
        try:
            colnames = a.colnames
        except AttributeError:
            colnames = list(a.dtype.names)
        
        for i in range(0,len(a)):
            row = []
            for j in range(0,len(a[i])):
                row.append(a[i][j])
            if i == 0:
                coltypes =[]
                for j in range(0,len(a[i])):
                    coltypes.append(type(a[i][j]))
            l.append(row)
        return l, colnames, coltypes
    except AttributeError:
        pass

    # 2D array or list
    try:
        for i in range(0,len(a)):
            row = []
            for j in range(0,len(a[i])):
                row.append(a[i][j])
            l.append(row)
            if i == 0:
                colnames = [str(x) for x in range(0,len(a[i]))]
                coltypes = []
                for j in range(0,len(a[i])):
                    coltypes.append(type(a[i][j]))
        return l, colnames, coltypes
    except KeyError:
        # dictionary
        for i in range(0,len(a)):
            l.append(a[i])

        return l, None, coltypes
    except TypeError:
        pass

    # 1D array or list
    try:
        row = []
        for j in range(0,len(a)):
            row.append(a[j])
            if j == 0:
                colnames = [str(x) for x in range(0,len(a))]
        l.append(row)
        return l, colnames, coltypes
    except TypeError:
        pass
    except KeyError:
        # dictionary
        l.append(a)
        return l, None, coltypes

    l = [[a]]
    colnames = ['0']
    return l, colnames, coltypes


def _file_to_data(file_name, sep="\t", has_header=None, types=None):
    """
    Convert the contents of a file to Table-readable data.
    
    Args:
        sep (str): field separator in text file, defaults to "\t"
        has_header (bool): does the file contain a header line?
        types (list): list of Python data types to be used for
                      Table columns. If None (default), data type
                      for every column will be str
    Returns:
        list: a list containing the rows of a table lists
        list: a list of columns names
        list: a list of row names
    """
    
    previous = []
    data = []
    rownames = []
    header = []
    with open(file_name,'r') as f:
        # header(?)
        
        if has_header:
            header_line = f.readline().rstrip()
            header = str.split(header_line, sep)

        i = 0
        line = f.readline().rstrip()
        while line != '':
            fields = str.split(line, sep)
            
            if i == 1 and (has_header is None) and (len(previous) == len(fields)-1):
                header = previous
            elif i > 0:

                if (header is not None) and (len(previous) == len(header) + 1):
                    rownames.append(previous[0])
                    previous = previous[1:]
    
                if types is not None:
                    for i in range(0, len(previous)):
                        previous[i] = types[i](previous[i])
    
                data.append(previous)
                
            previous = fields
            line = f.readline().rstrip()
            i += 1
                
        if (header is not None) and (len(previous) == len(header) + 1):
            rownames.append(previous[0])
            previous = previous[1:]

        if types is not None:
            for i in range(0, len(previous)):
                previous[i] = types[i](previous[i])
        data.append(previous)

    return data, header, rownames


class FileBased(object):
    def __init__(self, file_name=None, mode='a', tmpdir=None):
        # open file or keep in memory
        if hasattr(self, 'file'):
            if not isinstance(self.file, t.file.File):
                raise ValueError("'file' attribute already exists, but is no pytables File")
            else:
                return
        self.file = None
        self.tmp_file = None
        self.file_name = file_name
        self.tmp_file_name = None
        self._mode = mode
        if tmpdir is None:
            self.tmp_file_name = None
            self._init_file(file_name, mode)
        else:
            logging.info("Working in temporary directory...")
            self.tmp_file_name = os.path.join(tmpdir, self._generate_tmp_file_name())
            logging.info("Temporary output file: {}".format(self.tmp_file_name))
            if mode in ['w', 'x', 'w-']:
                pass
            elif mode in ['r+', 'r']:
                shutil.copyfile(file_name, self.tmp_file_name)
            elif mode in ['a'] and os.path.isfile(file_name):
                shutil.copyfile(file_name, self.tmp_file_name)
            self._init_file(self.tmp_file_name, mode)
        self.closed = False
    
    def close(self):
        self.file.close()
        self.closed = True

    def finalize(self):
        if self.closed:
            if self.tmp_file_name:
                if self._mode not in ['r']:
                    logging.info("Moving temporary output file to destination {}".format(self.file_name))
                    shutil.copyfile(self.tmp_file_name, self.file_name)
        else:
            raise IOError('The file has to be closed before copying the tmp file. Use close()')

    def cleanup(self):
        if self.closed:
            os.remove(self.tmp_file_name)
        else:
            raise IOError('The file has to be closed before deleting the tmp file. Use close()')

    def _generate_tmp_file_name(self):
        rand_str = binascii.b2a_hex(os.urandom(15))
        return "tmp_{}.h5".format(rand_str)

    def _init_file(self, file_name, mode):
        if file_name is None:
            self.file = create_or_open_pytables_file()
        elif type(file_name) == str:
            self.file = create_or_open_pytables_file(file_name, mode=mode)
        elif isinstance(file_name, t.file.File):
            self.file = file_name
        else:
            raise TypeError("file_name is not a recognisable type")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        if exc_type is None:
            if self.tmp_file_name:
                self.finalize()
                self.cleanup()
            return True
        else:
            if self.tmp_file_name:
                self.cleanup()
            return False

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, mode):
        if mode not in ['r', 'w', 'x', 'w-', 'r+', 'a']:
            raise ValueError('Unknown mode {}'.format(mode))
        if mode in ['r', 'r+'] and not os.path.isfile(self.file_name):
            raise OSError('The file {} does not exists'.format(self.file_name))
        if mode in ['x', 'w-'] and os.path.isfile(self.file_name):
            raise OSError('The file {} already exists'.format(self.file_name))
        self._mode = mode


class FileGroup(FileBased):
    def __init__(self, group, file_name=None, mode='a', tmpdir=None):
        FileBased.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir)

        try:
            group_node = self.file.get_node("/" + group)
            if not isinstance(group_node, t.group.Group):
                raise TypeError("%s is not a group, but %s" % (group, str(type(group_node))))
            self._group = group_node
        except t.NoSuchNodeError:
            self._group = self.file.create_group('/', group)


class TableRow(tuple):
    """
    In-memory row of a Table.
    
    Provides convenient access to fields using the Table column name
    or field index. Extends a Python tuple, so all expected tuple
    functionality is there.
    """
    
    def __new__(cls, t_tuple, rowname=None, colnames=None):
        obj = super(TableRow, cls).__new__(cls, tuple(t_tuple))
        obj.rowname = rowname
        obj.colnames = colnames
        obj._colnames_dict = {colnames[i]:i for i in range(0, len(colnames))}
        return obj

    # BUG FIX in python
    def __getslice__(self, start, stop):
        return self.__getitem__(slice(start, stop))

    def _get_cols(self, key):
        return self.__getitem__(key)

    def __getitem__(self, key):

        # do processing of "special" keys
        if type(key) is str:
            key = self._colnames_dict[key]

        try:
            res = []
            cn = []

            l = len(key)
            for i in range(0, l):
                if type(key[i]) is str:
                    key[i] =  self._colnames_dict[key[i]]
                res.append(super(TableRow, self).__getitem__(key[i]))
                cn.append(self.colnames[key[i]])
            return TableRow(res, rowname=self.rowname, colnames=cn)
        except TypeError:
            pass

        # get original value from pre-processed key
        cn = self.colnames[key]
        res = super(TableRow, self).__getitem__(key)

        if type(res) is tuple:
            return TableRow(res, rowname=self.rowname, colnames=cn)

        return res

    def __getattr__(self, name):
        if name in self._colnames_dict:
            return self.__getitem__(name)
        raise AttributeError


class TableCol(np.ndarray):
    """
    In-memory column of a Table.
    
    Provides convenient access to fields using the Table row name
    or field index. Extends a numpy array, so can be used in numpy
    calculations.
    """
    
    def __new__(cls, array, colname='0', rownames=None):
        obj = array.view(cls)
        obj.colname = colname
        obj.rownames = rownames
        obj._rownames_dict = {rownames[i]:i for i in range(0,len(rownames))}
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        # do any attribute setting here
        # e.g.
        # self.foo = getattr(obj, 'foo', None)
        self.colname = getattr(obj, 'colname', None)
        self.rownames = getattr(obj, 'rownames', None)

    # BUG FIX in python
    def __getslice__(self, start, stop):
        return self.__getitem__(slice(start, stop))

    def __getitem__(self, key):
        # do processing of "special" keys
        if type(key) is str:
            key = self._rownames_dict[key]

        try:
            res = []
            rn = []

            l = len(key)
            for i in range(0, l):
                if type(key[i]) is str:
                    key[i] =  self._rownames_dict[key[i]]
                res.append(super(TableCol, self).__getitem__(key[i]))
                rn.append(self.rownames[key[i]])
            return TableCol(np.array(res), colname=self.colname, rownames=rn)
        except TypeError:
            pass

        # get original value from pre-processed key
        rn = self.rownames[key]
        res = super(TableCol, self).__getitem__(key)

        if isinstance(res, np.ndarray):
            return TableCol(res, colname=self.colname, rownames=rn)

        return res

    def __getattr__(self, name):
        if name in self._rownames_dict:
            return self.__getitem__(name)
        raise AttributeError


class TableArray(np.ndarray):
    """
    Simplified Table version without pytables backing.
    
    Provides equivalent access to table data as Table class,
    but does not have advanced features such as indexing,
    querying, or saving to file.
    """
    
    def __new__(cls, array, colnames=None, rownames=None):
        obj = array.view(cls)
        # may set additional attributes (row names, etc)
        # e.g.
        # obj.foo = foo
        obj.colnames=colnames
        if colnames is not None:
            obj._colnames_dict = {colnames[i]:i for i in range(0, len(colnames))}
        obj.rownames=rownames
        if rownames is not None:
            obj._rownames_dict = {rownames[i]:i for i in range(0, len(rownames))}

        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        # do any attribute setting here
        # e.g.
        # self.foo = getattr(obj, 'foo', None)
        self.rownames = getattr(obj, 'rownames', None)
        self.colnames = getattr(obj, 'colnames', None)

    def __getitem__(self, key):
        if type(key) is tuple:
            if len(key) == 2:
                if type(key[0]) is slice and key[0].start is None and key[0].step is None and key[0].stop is None:
                    return self._get_cols(key[1])
                if type(key[1]) is slice and key[1].start is None and key[1].step is None and key[1].stop is None:
                    return self._get_rows(key[0])

                t = self._get_rows(key[0])

                return t._get_cols(key[1])

            else:
                raise KeyError("Unrecognized key " + str(key))

        return self._get_rows(key)

    def __getattr__(self, name):
        if name in self._colnames_dict:
            return self._get_cols(name)
        raise AttributeError

    def __getslice__(self, start, stop):
        return self.__getitem__(slice(start, stop))

    def _get_rows(self, key):
        rn = []

        # prepare structured numpy array
        # get column dtypes (always get all columns in _get_rows)
        dtype = self.dtype

        if type(key) is str:
            key = self._rownames_dict[key]

        try:
            res = []

            l = len(key)
            for i in range(0, l):
                if type(key[i]) is str:
                    key[i] =  self._rownames_dict[key[i]]
                res.append(super(TableArray, self).__getitem__(key[i]))
                rn.append(self.rownames[key[i]])

            if len(res) == 0 or len(dtype) > 0:
                x = np.zeros((len(res),), dtype=dtype)
                x[:] = res
            else:
                x = np.zeros((len(res), len(self.colnames)), dtype=dtype)

            return _structured_array_to_table_type(x, colnames=self.colnames, rownames=rn)

        except TypeError:
            pass

        # regular index or slice
        rows = super(TableArray, self).__getitem__(key)

        # multiple rows
        l = []
        try:

            n_rows = rows.shape[0]

            if n_rows > 1:
                for row in rows:
                    l.append(row)
            else:
                rows = rows[0]
                raise IndexError("Only one row")
        # single row
        except IndexError:
            l.append(rows)

        a = np.zeros((len(l),), dtype=dtype)
        a[:] = l

        rn = self.rownames[key]
        return _structured_array_to_table_type(a, colnames=self.colnames, rownames=rn)

    def _get_cols(self, key):
        cn = []

        dtype = self.dtype

        if type(key) is str:
            key = self._colnames_dict[key]

        # 1:5
        if type(key) is slice:
            start = 0 if key.start is None else key.start
            stop = len(self.colnames) if key.stop is None else (min(self.dim()[1],key.stop))
            step = 1 if key.step is None else key.step

            key = [i for i in range(start,stop,step)]

        try:
            dtypes = []
            l = len(key)
            for i in range(0, l):
                if type(key[i]) is str:
                    key[i] =  self._colnames_dict[key[i]]

                name = self.colnames[key[i]]
                cn.append(name)

                try:
                    dtypes.append((name, dtype[key[i]]))
                except KeyError:
                    dtypes = dtype

            a = np.zeros((len(self),), dtype=dtypes)
            a[:] = [tuple([super(TableArray, self).__getitem__(i)[k] for k in key]) for i in range(0, len(self))]

            return _structured_array_to_table_type(a, colnames=cn, rownames=self.rownames)

        except TypeError:
            pass

        try:
            dtypes = dtype[key]
        except KeyError:
            dtypes = dtype

        b = [(super(TableArray, self).__getitem__(i)[key],) for i in range(0, len(self))]
        name = self.colnames[key]
        mydt = np.dtype([(name,str(dtypes))])
        a = np.array(b,dtype=mydt)

        return _structured_array_to_table_type(a, colnames=name, rownames=self.rownames)

    def dim(self):
        """
        Return the number of rows and columns as a tuple
        """
        
        if len(self) > 0:
            return len(self), len(super(TableArray, self).__getitem__(0))
        return len(self), 0


class TableObject(object):
    __metaclass__ = ABCMeta
    
    @classmethod
    @abstractmethod
    def from_row(cls, row):
        pass
    
    def __getitem__(self, key):
        try:
            return self.__getattribute__(key)
        except AttributeError:
            raise IndexError("No item " + str(key) + " in object")


class Table(object):
    """
    Table class for saving ad accessing tabular data.
    
    Modeled a bit after panda's DataFrame (i.e. R-like functionality),
    but is inherently pytables-backed and has custom access to Table
    data through []-selectors.
    
    Also provides advanced querying capabilities through pytables.
    """
    
    def __init__(self, data=None, file_name=None,
                       ncols=0, nrows=0,
                       colnames=None, rownames=None,
                       col_types=None, default_type=str,
                       table_name='table', return_type=None):
        """
        Create a Table object.
        
        Can create an empty Table, load a previously saved Table
        object from file, or fill a table with data supplied as
        constructor argument(s).
        
        Table will be backed by an hdf5 dictionary, either on file
        or in memory.
        
        Creates a default pytables table with a standard '_rowname'
        field for keeping track of row names. Column names will be
        stored directly in the colnames object of the pytables table.
        
        Args:
            data (str, list): if type is string, will try to load
                              existing Table from hdf5 file or if 
                              that fails, will try to read data 
                              from tab-separated values file
                              else will treat argument as input
                              for the table, typically a list of
                              row data
            file_name (str):  file name as str to save the table
                              to - will create hdf5 dictionary at
                              this location
            ncols (int):      number of columns to create in Table
                              - will generally be inferred from
                              input data, only supply when creating
                              empty Table
            nrows (int):      number of rows to create in Table
                              - will generally be inferred from
                              input data, only supply when creating
                              empty Table
            colnames (list):  list of column names for Table
            rownames (list):  list of row names for Table
            col_types (list): list of Python data types to be used
                              as data types for Table columns; if
                              None will default to default_type
                              for each column
            default_type (type): the default column type as Python
                                 data type (str, int, float, bool)
            table_name (str): name to be used internally for storing
                              pytables table
        """
                
        # parse potential unnamed argument
        if data is not None:
            if type(data) is str:
                data = os.path.expanduser(data)
                # if file does not exist or file is hdf5 dict
                # treat it as an hdf5 file name
                if not os.path.isfile(data) or is_hdf5_file(data):
                    file_name = data
                    data = None
                # else it is treated as a data file...
                else:
                    data, tmp_colnames, tmp_rownames = _file_to_data(data, has_header=True, types=col_types)
                    if colnames is None:
                        colnames = tmp_colnames
                    if rownames is None:
                        rownames = tmp_rownames
            else:
                data, c, ty = _to_list_and_names(data)
                if colnames is None:
                    colnames = c
                if col_types is None:
                    col_types = ty
        
        # check return type
        if return_type is not None and not issubclass(return_type, TableObject):
            raise ValueError("return_type must inherit TableObject")
        else:
            self._row_type = return_type
        
        # open file or keep in memory
        if hasattr(self, 'file'):
            if not isinstance(self.file, t.file.File):
                raise ValueError("'file' attribute already exists, but is no pytables File")
        else:
            if file_name is None:
                self.file = create_or_open_pytables_file()
            elif type(file_name) == str:
                self.file = create_or_open_pytables_file(file_name)
            else:
                self.file = file_name
            
        # set a few sensible defaults
        self._rowname_field = '_rowname'

        # create or retrieve table
        if not table_name in self.file.root:

            # set a few sensible defaults
            self._rowname_field = '_rowname'

            # columns
            if colnames is not None:
                ncols = max(ncols, len(colnames))
            else:
                colnames = []

            for i in range(len(colnames), ncols):
                colnames.append(str(i))

            if col_types is None:
                col_types = []

            for i in range(len(col_types), ncols):
                col_types.append(default_type)

            # rows
            if rownames is not None:
                nrows = max(nrows, len(rownames))
            else:
                rownames = []

            if data is not None:
                nrows = max(nrows, len(data))

            for i in range(len(rownames), nrows):
                rownames.append(str(i))

            self._table_name = table_name

            columns_dict = dict()
            columns_dict[self._rowname_field] = t.StringCol(50, pos=0)
            
            for i, colname in enumerate(colnames):
                columns_dict[colname] = _convert_to_tables_type(col_types[i], pos=i+1)
                        
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                self._table = self.file.create_table("/", table_name, columns_dict)
        else:
            self._table = self.file.get_node('/' + table_name)

        # clean up
        self._table.flush()
        
        # data is file
        if type(data) is str:
            data = _file_to_data(data, types=col_types)

        # add any potential data
        if data is not None:
            logging.info("Adding data")
            self.append(data)
        elif nrows > 0:
            dt = self._table[0:0].dtype
            dtypes = []
            for name in dt.names:
                if name is not self._rowname_field:
                    dtypes.append((name,dt[name]))
            data = np.zeros((nrows,), dtype=dtypes)
            self.append(data)
        self.set_rownames(rownames)
        
    def _flush(self):
        self._table.flush()

    def save_as(self, file_name, table_name='table'):
        """
        Save current table at different location.
        
        Creates a copy of the hdf5 file backing this table
        and saves it at the specified location. Will update
        local variables so that the Table continues to 
        behave as expected, but backed by the new file.
        
        Args:
            file_name (str):  path of the target file
            table_name (str): name to be used internally for storing
                              pytables table
        """
        
        self.file.copy_file(file_name)
        self.file.close()
        self.file = create_or_open_pytables_file(file_name)
        self._table = self.file.get_node('/' + table_name)
        
    def export(self, file_name, sep="\t", include_colnames=True, include_rownames=True):
        """
        Export table as plain-text file.
        
        By default, writes Table data to the provided location
        with column names and row names.
        
        Args:
            file_name (str):         location to export the Table to
            sep (str):               field separator (default '\t')
            include_colnames (bool): write column names in first line
                                     of file
            include_rownames (bool): write row names as first field
                                     in each line
        """
        
        rn = self.rownames
        cn = self.colnames
        with open(file_name,'w') as o:
            if include_colnames:
                for j in range(0,len(cn)):
                    o.write(cn[j])
                    if j < len(cn)-1:
                        o.write(sep)
                    else:
                        o.write("\n")
            
            for i in range(0,len(self)):
                row = self[i]
                if include_rownames:
                    o.write(rn[i] + sep)
                for j in range(0,len(row)):
                    o.write(row[j])
                    if j < len(row)-1:
                        o.write(sep)
                    else:
                        o.write("\n")

    def close(self):
        """
        Close the underlying hdf5 dictionary file.
        """
        self.file.close()
    
    @property
    def colnames(self):
        """
        Return the column names of this Table as a list.
        """
        return self._table.colnames[1:]
        
    @property
    def rownames(self):
        """
        Return the row names of this Table as a list.
        """
        return self._table[:]['_rowname']
    
    def set_rownames(self, names):
        """
        Set the row names of this Table to this list.
        
        Args:
            names (list): a list of names to be used as rownames.
                          Length must match Table length.
        
        Raises:
            ValueError:  if names are not same length as Table
        """
        if names is not None:
            if not len(names) == len(self):
                raise ValueError("names must be the same length as table")
            i = 0
            for row in self._table.iterrows():
                row[self._rowname_field] = names[i]
                row.update()
                i += 1

            self._table.flush()

    def _append_row_list(self, l, flush=True, rowname=None):
        """
        Append a list representing a row to the Table.
        
        Args:
            l (list):      a list the same length as table columns,
                           representing a row, OR a list with one 
                           additional field representing a row with
                           the row name in the first field
            flush (bool):  write data to disk after appending to 
                           Table (default True, but False combined
                           with manual flush is useful in bulk import)
            rowname (str): row name of the imported row. If None
                           defaults to the current length of the
                           Table (or the first field of the row,
                           see above)
        """
        row = self._table.row

        try:

            if len(l) == len(self._table.colnames)-1:
                if rowname is not None:
                    row[self._rowname_field] = rowname
                else:
                    row[self._rowname_field] = len(self)
                for i in range(1,len(self._table.colnames)):
                    row[self._table.colnames[i]] = l[i-1]
            else:
                for i in range(0,len(self._table.colnames)):
                    row[self._table.colnames[i]] = l[i]

        except (TypeError, KeyError, IndexError, ValueError):
            raise TypeError("l is not a list")

        row.append()

        if flush:
            self._table.flush()

    def _append_row_dict(self, d, flush=True, rowname=None):
        """
        Append a dict representing a row to the Table.
        
        Args:
            d (dict):      a dict where keys are column names
            flush (bool):  write data to disk after appending to 
                           Table (default True, but False combined
                           with manual flush is useful in bulk import)
            rowname (str): row name of the imported row. If None
                           defaults to the current length of the
                           Table
        """
        
        row = self._table.row

        try:
            if not self._rowname_field in d:
                if rowname is not None:
                    d[self._rowname_field] = rowname
                else:
                    d[self._rowname_field] = len(self)
            for name in d:
                row[name] = d[name]
        except (TypeError, KeyError, AttributeError, IndexError, ValueError):
            raise TypeError("d is not a dictionary")

        row.append()

        if flush:
            self._table.flush()

    def append(self, data, rownames=None):
        """
        Append data to the table.
        
        Args:
        data (...): Data can have several formats:
            (a) list of lists (or with list-like indexing)
                each entry in the list represents one row
                with fields in the same order as Table
                columns.
                Rows may include row name as first field.
            (b) list of dicts
                each entry in the list is a dict, keys in
                the dict are column names. To import row
                names, add '_rowname' entry
            (c) TSV file
                A tab-separated value (TSV) file without
                header line. May include row names in 
                first field
            (d) list
                A single list representing a row in the
                Table
            (e) dict
                A single dict with keys corresponding to
                table column names. To import row
                names, add '_rowname' entry
        rownames (list): optional list of rownames, overwritten by
            any rownames that may be in the data
        """
        
        # data is file
        if type(data) is str:
            data = _file_to_data(data)
        
        original_length = len(self)
        # data is a list of lists?
        try:
            i = 0
            for l in data:
                if rownames is not None:
                    self._append_row_list(l, flush=False, rowname=rownames[i])
                else:
                    self._append_row_list(l, flush=False, rowname=str(original_length+i))
                i += 1
            self._table.flush()
            return
        except TypeError:
            pass

        # data is a list of dicts?
        try:
            i = 0
            for d in data:
                if rownames is not None:
                    self._append_row_dict(d, flush=False, rowname=rownames[i])
                else:
                    self._append_row_dict(d, flush=False,rowname=str(original_length+i))
                i += 1
            self._table.flush()
            return
        except TypeError:
            pass
        
        try:
            rowname = rownames[0]
        except TypeError:
            rowname = rownames
            
        # data is a list?
        try:
            self._append_row_list(data, flush=True, rowname=rowname)
            return
        except TypeError:
            pass

        # data is a dictionary
        try:
            self._append_row_dict(data, flush=True, rowname=rowname)
            return
        except TypeError:
            raise ValueError("Data type unsupported")

    def __len__(self):
        return len(self._table)

    def dim(self):
        """
        Return the number of rows and columns as a tuple
        """
        return len(self), len(self._table.colnames)-1

    def __repr__(self, *args, **kwargs):
        # find out maximum column width
        show = 9
        m = 3 if show < len(self) else 0

        max_col_width = [max(len(str(x)),m) for x in self._table.colnames]
        max_col_width[0] = 0
        for i in range(0,min(show,len(self))):
            for j in range(0,len(max_col_width)):
                strlen = min(20,len(str(self._table[i][j])))
                max_col_width[j] = max(max_col_width[j], strlen)

        r = ' ' * max_col_width[0] + ' '
        for j in range(1,len(self._table.colnames)):
            s = str(self._table.colnames[j])
            r += s.rjust(max_col_width[j]) + " "
        r += "\n"

        for i in range(0,min(show,len(self))):
            for j in range(0,len(max_col_width)):
                s = str(self._table[i][j])
                if len(s) > 17:
                    s = s[0:17] + '...'
                r += s.rjust(max_col_width[j]) + " "
            if i < len(self)-1:
                r += "\n"

        r += ' ' * max_col_width[0] + ' '
        if show < len(self):
            for j in range(1,len(self._table.colnames)):
                r += '...'.rjust(max_col_width[j]) + ' '
            r += '\n%d rows' % len(self)

        return r

    def __del__(self):
        self.file.close()

    def _get_rows(self, key, use_row_type=True):

        rn = []
        cn = []

        # prepare structured numpy array
        # get column dtypes (always get all columns in _get_rows)
        table_dtypes = self._table[0:0].dtype
        dtypes = []
        for i in range(1,len(self._table.colnames)):
            name = self._table.colnames[i]
            cn.append(name)
            dt = str(table_dtypes[name])
            dtypes.append((name,dt))

        l = []
        # "rowname" or [] or ["rowname1", "rowname2", ...] or [1, 3, ...]
        if type(key) is str or type(key) is list:

            # "rowname"
            if type(key) is str:
                # query pytables table
                condition = "%s == '%s'" % (self._rowname_field,key)
                rows = [x.fetch_all_fields() for x in self._table.where(condition)]
                for row in rows:
                    rn.append(row[0])
                    l.append(tuple(row)[1:])

            # [] or ["rowname1", "rowname2", ...] or [1, 3, ...]
            else:
                # []
                if len(key) == 0:
                    pass
                # ["rowname1", "rowname2", ...]
                elif type(key[0]) is str:
                    # concatenate rows with name in list
                    rows = []
                    for k in key:
                        rows = rows + [x.fetch_all_fields() for x in self._table.where("%s == '%s'" % (self._rowname_field,k))]
                    for row in rows:
                        rn.append(row[0])
                        l.append(tuple(row)[1:])
                # [1, 3, ...]
                else:
                    for i in key:
                        row = self._table[i]
                        rn.append(row[0])
                        l.append(tuple(row)[1:])

        # regular index or slice
        else:
            rows = self._table[key]

            # multiple rows
            try:
                n_rows = rows.shape[0]

                if n_rows > 1:
                    for row in rows:
                        rn.append(row[0])
                        l.append(tuple(row)[1:])
                else:
                    rows = rows[0]
                    raise IndexError("Only one row")
            # single row
            except IndexError:
                rn.append(rows[0])
                l.append(tuple(rows)[1:])

        a = np.zeros((len(l),), dtype=dtypes)
        a[:] = l
        
        if use_row_type:
            return _structured_array_to_table_type(a, rownames=rn, colnames=cn, return_type=self._row_type)
        return _structured_array_to_table_type(a, rownames=rn, colnames=cn)

    def _get_cols(self, key):
        cn = []
        rn = []

        # 1 - convert to string
        if type(key) is int:
            key = self._table.colnames[key+1]

        # 'A'
        if type(key) is str:
            dt = self._table[0:0].dtype[key]
            b = [(x[key],) for x in self._table.iterrows()]
            mydt = np.dtype([(key,str(dt))])
            a = np.array(b,dtype=mydt)

            cn.append(key)
            rn = [x[self._rowname_field] for x in self._table.iterrows()]

            return _structured_array_to_table_type(a, rownames=rn, colnames=cn)

        # 1:5
        if type(key) is slice:
            start = 0 if key.start is None else key.start
            stop = len(self._table.colnames)-1 if key.stop is None else (min(self.dim()[1],key.stop))
            step = 1 if key.step is None else key.step

            key = [self._table.colnames[i+1] for i in range(start,stop,step)]

        try:
            dtypes = []
            keys = []
            for k in key:
                if type(k) is int:
                    k = self._table.colnames[k+1]
                keys.append(k)
                dt = self._table[0:0].dtype[k]
                dtypes.append((k,dt))

            cn = keys
            a = np.zeros((len(self),), dtype=dtypes)
            a[:] = [tuple([x[k] for k in keys]) for x in self._table.iterrows()]
            rn = [x[self._rowname_field] for x in self._table.iterrows()]

        except IndexError:
            raise IndexError("Cannot retrieve column with value " + str(key))

        return _structured_array_to_table_type(a, rownames=rn, colnames=cn)

    def __getitem__(self, key):

        if type(key) is tuple:
            if len(key) == 2:
                if type(key[0]) is slice and key[0].start is None and key[0].step is None and key[0].stop is None:
                    return self._get_cols(key[1])
                if type(key[1]) is slice and key[1].start is None and key[1].step is None and key[1].stop is None:
                    return self._get_rows(key[0])

                t = self._get_rows(key[0], use_row_type=False)

                return t._get_cols(key[1])

            else:
                raise KeyError("Unrecognized key " + str(key))

        return self._get_rows(key)
    
    def _get_item_as_table_object(self, key, cls):
        res = self._get_rows(key)
                
        if isinstance(res, TableRow):
            return cls.from_row(res)
        elif isinstance(res, TableArray):
            l = []
            for row in res:
                l.append(cls.from_row(row))
        return res

    def __getattr__(self, name):
        if name in self.colnames:
            return self._get_cols(name)
        return self.__getattribute__(name)

    def __iter__(self):
        low = 0
        high = len(self)-1
        table = self
        class Iter:
            def __init__(self, low, high):
                self.current = low
                self.high = high

            def __iter__(self):
                return self

            def next(self):
                if self.current > self.high:
                    raise StopIteration
                self.current += 1
                return table._get_rows(self.current-1)

        return Iter(low,high)

    def as_array(self):
        """
        Return table as numpy record array.
        """
        return self[:]

    def where(self, query):
        """
        Query Table and return sub-Table.
        
        Use pytables-style syntax
        (http://www.pytables.org/usersguide/condition_syntax.html)
        to query the Table.
        
        Args:
            query: A pytables-style query string, e.g.
                   "(chromosome == 'chr1') & (start > 25000)"
        
        Return:
            TableArray: A TableArray with the subset
                        of rows matching the query
        """
        
        rn = []
        cn = []

        # prepare structured numpy array
        # get column dtypes (always get all columns in _get_rows)
        table_dtypes = self._table[0:0].dtype
        dtypes = []
        for i in range(1,len(self._table.colnames)):
            name = self._table.colnames[i]
            cn.append(name)
            dt = str(table_dtypes[name])
            dtypes.append((name,dt))
        
        l = []
        rows = [x.fetch_all_fields() for x in self._table.where(query)]
        for row in rows:
            rn.append(row[0])
            l.append(tuple(row)[1:])
        
        a = np.zeros((len(l),), dtype=dtypes)
        a[:] = l

        return _structured_array_to_table_type(a, rownames=rn, colnames=cn, return_type=self._row_type)


class Mask(object):
    """
    Class providing Mask details.
    """
    
    def __init__(self, ix, name, description=''):
        self.ix = ix
        self.name = name
        self.description = description

    def __repr__(self):
        return "%d. %s: %s" % (self.ix, self.name, self.description)


class Maskable(FileBased):
    """
    Class that adds masking functionality to tables.
    
    It can be desirable to hide a subset of data in 
    a table rather than deleting it outright. This is 
    especially useful in pytables tables, as they are
    not designed for frequent deletions and may 
    suffer from performance problems as a consequence.
    The process of hiding data is referred to as 
    'masking'.
    
    This class creates a pytables Table that contains
    meta-information about masks, i.e. mask IDs, mask
    names, and mask descriptions. Since it is backed
    by pytables, Mask descriptions will be saved to 
    the corresponding hdf5 file automatically. This 
    class additionally provides tools to create and 
    retrieve Mask descriptions.
    
    Note that this class only provides meta-
    information for Masks. Masking can technically 
    performed without Maskable, but then the reasons
    for masking will not be saved.
    
    The actual masking is performed using a MaskFilter
    on a MaskedTable object using the 
    MaskedTable.filter(MaskFilter) function. 
    """
    
    class MaskDescription(t.IsDescription):
        ix = t.Int16Col(pos=0)
        name = t.StringCol(50, pos=1)
        description = t.StringCol(255, pos=2)

    def __init__(self, data=None, file_name=None, table_name="mask", tmpdir=None):
        """
        Enable recording of masking in pytables-backed object.
        
        This constructor is built with considerable flexibility
        in the ways it can be called. The 'data' parameter
        serves here as a multi-purpose container that simplifies
        calling the constructor with just one argument.
        
        Args:
            data (...):         multi-purpose container, simplifies one-argument
                                constructor calls. Can stand in for a pytables 
                                Table or File, or a path.
                (None):         An existing mask description table is expected in 
                                self._mask. If none is found, the constructor will 
                                look for a pytables File in self.file and attempt 
                                to create the mask description table in there.
                (str):          An exisiting pytable File at this location will be
                                opened in append mode or one will be created.
                (tables.file.File):
                                A mask description table will be created in this
                                pytables File.
                (tables.table.Table):
                                This pytables Table will be used as mask 
                                description table (if it has the necessary fields
                                ix, name, description).
            table_name (str):   name of the mask table that is created in
                                pytables file, does not usually need to be 
                                modified
        """
        # parse potential unnamed argument
        if data is not None:
            # data is file name
            if type(data) is str or isinstance(data, t.file.File):
                if file_name is None:
                    file_name = data
                    data = None
            elif type(data) == t.table.Table:
                self._set_mask_table(data)
        else:
            if hasattr(self, '_mask'):
                self._set_mask_table(self._mask)

        FileBased.__init__(self, file_name, tmpdir=tmpdir)
                
        if (not hasattr(self, '_mask') or self._mask is None) and self.file is not None:
            try:
                self._mask = self.file.get_node('/' + table_name)
            except NoSuchNodeError:
                self._mask = self.file.create_table("/", table_name, Maskable.MaskDescription)
                row = self._mask.row
                row['ix'] = 0
                row['name'] = 'default'
                row['description'] = 'Default mask'
                row.append()
                self._mask.flush()
        
    def _set_mask_table(self, table):
        if type(table) == t.table.Table:
            if ('ix' not in table.colnames or
                    'name' not in table.colnames or
                    'description' not in table.colnames):
                raise ValueError("Object already has a mask table, \
                                  but it does not have all the necessary \
                                  columns (ix, name, description)")
            self._mask = table
        else:
            raise ValueError("Table is not a MaskDescription table")

    def add_mask_description(self, name, description):
        """
        Add a mask description to the _mask table and return its ID.
        

        Args:
            name (str):        name of the mask
            description (str): description of the mask
            
        Returns:
            int:    id of the mask
        """
        ix = len(self._mask)
        row = self._mask.row
        row['ix'] = ix
        row['name'] = name
        row['description'] = description
        row.append()
        self._mask.flush()
        self.clear_caches()
        
        return Mask(ix, name, description)

    def clear_caches(self):
        self.get_mask.cache_clear()
        self.get_masks.cache_clear()

    def masks(self):
        this = self

        class MaskIter:
            def __init__(self):
                self.iter = iter(this._mask)

            def __iter__(self):
                return self

            def next(self):
                row = self.iter.next()
                return Maskable._row_to_mask(row)

        return MaskIter()

    @staticmethod
    def _row_to_mask(row):
        return Mask(name=row['name'], ix=row['ix'], description=row['description'])

    def get_binary_mask_from_masks(self, masks):
        o = []
        for m in masks:
            if type(m) == type('str'):
                o.append(2**self.get_mask(m).ix)
            elif isinstance(m, MaskFilter):
                o.append(2**m.mask_ix)
            elif isinstance(m, Mask):
                o.append(2**m.ix)
            elif type(m) == type(1):
                o.append(2**m)
            else:
                raise ValueError('Can only get binary mask from mask names, indexes and MaskFilter instances')
        return sum(o)

    @lru_cache(maxsize=1000)
    def get_mask(self, key):
        """
        Search _mask table for key and return Mask.
        
        Args:
            key (str): search by mask name
            key (int): search by mask ID
            
        Returns:
            Mask
        """
        
        if type(key) == int:
            for row in self._mask.where("ix == %d" % key):
                return Maskable._row_to_mask(row)
        else:
            for row in self._mask.where("name == '%s'" % str(key)):
                return Maskable._row_to_mask(row)
        return KeyError("Unrecognised key type")
    
    @lru_cache(maxsize=1000)
    def get_masks(self, ix):
        """
        Extract mask IDs encoded in parameter and return masks.
        
        IDs are powers of 2, so a single int field in the table can hold
        multiple masks by simply adding up the IDs. Similar principle to
        UNIX chmod (although that uses base 8)
        
        Args:
            ix (int): integer that is the sum of powers of 2. Note that this value
                      is not necessarily itself a power of 2.
            
        Returns:
            list (Mask): list of Masks extracted from ix
        """
        
        key_copy = ix

        last_ix = int(self._mask[-1][0])
        
        masks = []
        while last_ix > 0:
            if key_copy - 2**last_ix >= 0:
                mask = self.get_mask(last_ix)
                if mask is not None:
                    masks.append(mask)
                key_copy -= 2**last_ix
            last_ix -= 1
        if key_copy > 0:
            masks.append(self.get_mask(0))
        
        return list(reversed(masks))

    def mask_statistics(self, table, include_unmasked=True):
        masks = defaultdict(int)
        if include_unmasked:
            masks['unmasked'] = 0

        stats = defaultdict(int)

        for row in table._iter_visible_and_masked():
            stats[row[table._mask_field]] += 1

        for mask_bit, count in stats.iteritems():
            row_masks = self.get_masks(mask_bit)

            found_masks = False
            for mask in row_masks:
                found_masks = True
                masks[mask.name] += count

            if not found_masks and include_unmasked:
                masks['unmasked'] += count
        return masks


class MaskedTableView(object):
    def __init__(self, masked_table, it=None, excluded_masks=0):
        self.masked_table = masked_table
        self.iter = it if it else self.masked_table._iter_visible_and_masked()
        self.excluded_mask_ix = excluded_masks

    def __iter__(self):
        return self

    def next(self):
        row = self.iter.next()
        # bit-shift magic! Go @alexis!
        # a is a subset of b if and only if a | b == b.
        # If this condition is satisfied for each byte, return TRUE. Otherwise return FALSE
        while row[self.masked_table._mask_field] | self.excluded_mask_ix != self.excluded_mask_ix:
            row = self.iter.next()
        return row


class MaskedTable(t.Table):
    """
    Wrapper that adds masking functionality to a pytables table. 
    
            
    MaskedTable is simply a wrapper around a
    pytables Table with certain columns. It
    adds masking functionality to table rows.
        
    It can be desirable to hide a subset of data in 
    a table rather than deleting it outright. This is 
    especially useful in pytables tables, as they are
    not designed for frequent deletions and may 
    suffer from performance problems as a consequence.
    
    This class extends the functionality of classes
    that work on the basis of a pytables table. It
    allows the masking of rows and provides table iterators
    and len() functions that ignore masked rows, and 
    therefore give the illusion of deleted rows.
    Masked rows, however, can be inspected separately
    at any point.
    
    The minimum requirement for a table to be masked
    is the inclusion of a 'mask' column. For full
    iterator and len() functionality, the table also
    needs to include a 'mask_ix' column. Checks for
    this are performed when specifying the table to
    be masked.
    """
    
    _c_classid = 'MASKEDTABLE'
    
    def __init__(self, parentnode, name, description=None,
                 title="", filters=None, expectedrows=None,
                 chunkshape=None, byteorder=None, _log=False,
                 mask_field='_mask', mask_index_field='_mask_ix'):
        """
        Pytables Table extension to provide masking functionality.
        
        Args:
            table (tables.Table):
                pytables Table with at least a 'mask' column.
                'mask_ix' column required for full indexing
                functionality.
        """
        
        # set instance variables
        self._queued_filters = []
        self._mask_field = mask_field
        self._mask_index_field = mask_index_field
        
        if description is not None:
            # try converting description to dict
            try:
                description = description.columns
            except AttributeError:
                pass
            
            # fill in fields required for masking
            if isinstance(description, dict):
                masked_description = description.copy()
            else:
                raise ValueError("Unrecognised description type (%s)" % str(type(description)))
            
            # check that reserved keys are not used
            if masked_description.has_key(mask_field):
                raise ValueError("%s field is reserved in MaskedTable!" % mask_field)
            if masked_description.has_key(mask_index_field):
                raise ValueError("%s field is reserved in MaskedTable!" % mask_index_field)
            
            # add mask fields to description
            masked_description[mask_field] = t.Int32Col()
            masked_description[mask_index_field] = t.Int64Col()
        else:
            masked_description = None
                
        t.Table.__init__(self, parentnode, name,
                         description=masked_description, title=title,
                         filters=filters,
                         expectedrows=expectedrows,
                         _log=_log,
                         chunkshape=chunkshape,
                         byteorder=byteorder)
        
        mask_ix_col = getattr(self.cols, self._mask_index_field)
        if not mask_ix_col.is_indexed:
            mask_ix_col.create_index()

    def flush(self, update_index=False):
        """
        Flush buffered rows.
        
        Also updates the mask index, if requested.
        """
        self._flush(update_index)
    
    def _flush(self, update_index=False):
        # commit any previous changes
        t.Table.flush(self)

        if update_index:
            self._update_ix()
            # force flush of index if
            # autoindex is disabled
            if not self.autoindex:
                self.flush_rows_to_index()
            # commit index changes
            t.Table.flush(self)

    def iterrows(self, start=None, stop=None, step=None, excluded_masks=0):
        it = t.Table.iterrows(self, start, stop, step)
        return MaskedTableView(self, it, excluded_masks=excluded_masks)

    def itersorted(self, sortby, checkCSI=False,
                   start=None, stop=None, step=None, excluded_masks=0):
        it = t.Table.itersorted(self, sortby, checkCSI=checkCSI,
                                start=start, stop=stop, step=step)
        return MaskedTableView(self, it, excluded_masks=excluded_masks)

    def _iter_visible_and_masked(self):
        """
        Return an iterator over all rows, including masked ones.
        """
        return t.Table.iterrows(self)

    def masked_rows(self):
        """
        Return an iterator over masked rows.
        """
        this = self
        it = self._iter_visible_and_masked()

        class MaskedRows(MaskedTableView):
            def __init__(self, masked_table, it):
                super(MaskedRows, self).__init__(masked_table, it)

            def next(self):
                row = self.iter.next()
                while row[this._mask_field] == 0:
                    row = self.iter.next()
                return row

            def __getitem__(self, key):
                if isinstance(key, int):
                    if key >= 0:
                        key = -1*key - 1
                        res = [x.fetch_all_fields() for x in
                               super(MaskedTable, this).where("%s == %d" % (this._mask_index_field, key))]
                        if len(res) == 1:
                            return res[0]
                        if len(res) == 0:
                            raise IndexError("Index %d out of bounds" % key)
                        raise RuntimeError("Duplicate row for key %d" % key)
                    else:
                        l = this._visible_len()
                        return self[l+key]
                else:
                    raise KeyError('Cannot retrieve row with key ' + str(key))

        return MaskedRows(self, it)

    def __getitem__(self, key):
        return self._get_visible_item(key)
    
    def _get_visible_item(self, key):
        if type(key) == slice:
            res = []
            # set sensible defaults
            start = key.start
            if start is None:
                start = 0
            stop = key.stop
            if stop is None:
                stop = len(self)
            step = key.step
            if step is None:
                step = 1
                  
            for i in range(start, stop, step):
                res.append(self._get_visible_item(i))
            return res
        else:
            try:
                # this could be a numpy int
                try:
                    key = key.item()
                except AttributeError:
                    pass
                
                key = int(key)
                
                if key >= 0:
                    res = [x.fetch_all_fields() for x in self.where("%s == %d" % (self._mask_index_field,key))]
                    if len(res) == 1:
                        return res[0]
                    if len(res) == 0:
                        raise IndexError("Index %d out of bounds" % key)
                    raise RuntimeError("Duplicate row for key %d" % key)
                else:
                    l = self._visible_len()
                    if l+key >= 0:
                        return self._get_visible_item(l+key)
                    raise KeyError('Cannot retrieve row with key %s' % str(key))
            except ValueError as e:
                raise KeyError('Cannot retrieve row with key %s (%s)' % (str(key), str(e)))
    
    def _original_getitem(self, key):
        return t.Table.__getitem__(self, key)
    
    def __len__(self):
        """
        Return the 'perceived' length of the masked table.
          
        If the table has masked rows, these will not be counted.
        """
          
        return self._visible_len()
    
    def _visible_len(self):
        if 'masked_length' not in self.attrs or self.attrs['masked_length'] == -1:
            return sum(1 for _ in iter(self.where("%s >= 0" % self._mask_index_field)))
        return int(self.attrs['masked_length'])
    
    def _original_len(self):
        return t.Table.__len__(self)
    
    # new index update method
    def _update_ix(self):
        """
        Update the row indexes of the Table.
        
        Should be run after masking. Will assign auto-
        incrementing integers (from 0) to the mask index
        field of each row in the table if it is not 
        masked, -1 otherwise.
        """
        
        logging.info("Updating mask indices")

        l = self._original_len()
        with RareUpdateProgressBar(max_value=l) as pb:
            ix = 0
            masked_ix = -1
            for i, row in enumerate(self._iter_visible_and_masked()):
                if row[self._mask_field] > 0:
                    row[self._mask_index_field] = masked_ix
                    masked_ix -= 1
                else:
                    row[self._mask_index_field] = ix
                    ix += 1
                row.update()
                pb.update(i)
            self.attrs['masked_length'] = ix
    
    @lru_cache(maxsize=1000)
    def _get_masks(self, binary_mask):
        def bits(n):
            while n:
                b = n & (~n+1)
                yield b
                n ^= b
        return list(bits(binary_mask))

    def _row_masks(self, row):
        return self._get_masks(row[self._mask_field])

    def _has_mask(self, row, mask):
        return mask in self._row_masks(row)

    def filter(self, mask_filter, _logging=False):
        """
        Run a MaskFilter on this table.
        
        This functions calls the MaskFilter.valid function on
        every row and masks them if the function returns False.
        After running the filter, the table index is updated
        to match only unmasked rows.

        :param mask_filter: A :class:`~MaskFilter` object
        :param _logging: Print progress to stderr
        """

        total = 0
        ix = 0
        mask_ix = -1

        # progress bar
        l = self._original_len()
        pb = RareUpdateProgressBar(max_value=l)
        if _logging:
            pb.start()

        # statistics
        stats = defaultdict(int)

        for i, row in enumerate(self._iter_visible_and_masked()):
            total += 1

            if not mask_filter.valid(row):
                row[self._mask_field] += 2**mask_filter.mask_ix

            stats[row[self._mask_field]] += 1

            # update index
            if row[self._mask_field] > 0:
                row[self._mask_index_field] = mask_ix
                mask_ix -= 1
            else:
                row[self._mask_index_field] = ix
                ix += 1
            row.update()

            if _logging:
                pb.update(i)
        self.attrs['masked_length'] = ix

        if _logging:
            pb.finish()
            logging.info("Total: %d. Filtered: %d" % (total, -1*(mask_ix-1)))
                    
        self.flush(update_index=False)

        return stats

    def queue_filter(self, filter_definition):
        """
        Add a MaskFilter to filter queue.
        
        Queued filters can be run at a later time using
        the run_queued_filters function.
        """
        self._queued_filters.append(filter_definition)
        
    def run_queued_filters(self, _logging=False):
        """
        Run queued MaskFilters.
        
        MaskFilters can be queued using the
        queue_filter function.

        :param _logging: If True, prints log to stderr
        """
        ix = 0
        mask_ix = -1
        total = 0

        # progress bar
        l = self._original_len()
        pb = RareUpdateProgressBar(max_value=l)
        if _logging:
            pb.start()

        for i, row in enumerate(self._iter_visible_and_masked()):
            total += 1
            for f in self._queued_filters:
                if not f.valid(row):
                    row[self._mask_field] += 2**f.mask_ix
                    
            # update index
            if row[self._mask_field] > 0:
                row[self._mask_index_field] = mask_ix
                mask_ix -= 1
            else:
                row[self._mask_index_field] = ix
                ix += 1
            row.update()
            
            if _logging:
                pb.update(i)

        self.attrs['masked_length'] = ix

        if _logging:
            pb.finish()
            logging.info("Total: %d. Filtered: %d" % (total, mask_ix-1))

        self.flush(update_index=False)

    def where(self, condition, condvars=None,
              start=None, stop=None, step=None):
        condition = "(" + condition + ") & (%s >= 0)" % self._mask_index_field
        return super(MaskedTable, self).where(condition, condvars=None,
                                              start=None, stop=None, step=None)


class MaskFilter(object):
    """
    Abstract class that defines a filter for MaskedTable.
    """

    __metaclass__ = ABCMeta
    
    def __init__(self, mask = None, mask_ix=0, mask_name='default', mask_description="Default mask."):
        """
        Create a MaskFilter.
        
        Sets values for mask index, name, and description
        to be used by MaskedTable.filter function.
        
        Args:
            mask (Mask, int):    if Mask, the attributes mask_ix,
                                 mask_name, and mask_description
                                 will be copied here.
                                 if int, mask_name and mask_description
                                 will be set to defaults, while
                                 mask_ix will be set to this value.
            mask_ix (int):       Index of the mask to be applied.
                                 Defaults to 0
            mask_name (str):     Name of the mask to be applied.
                                 Defaults to 'default'
            mask_description (str):
                                 Description of the mask to be applied.
                                 Defaults to 'Default mask'
        """
        
        if mask is not None:
            if isinstance(mask, Mask):
                self.mask_ix = mask.ix
                self.maks_name = mask.name
                self.mask_description = mask.description
            else:
                self.mask_ix = mask
                self.maks_name = mask_name
                self.mask_description = mask_description
        else:
            self.mask_ix = mask_ix
            self.maks_name = mask_name
            self.mask_description = mask_description
            
        if not type(self.mask_ix) == int:
            raise ValueError("mask_ix must be an integer!")
    
    @abstractmethod
    def valid(self, row):
        """
        Test if a row is valid according to this filter.
        
        Args:
            row (tables.tableextension.Row):
                A pytables Table Row
        
        Returns:
            bool: True if row is valid, False otherwise
        """
        pass
        

class Meta(object):
    def __init__(self, date=None, message='', name='', value=None, level=0, category=''):
        self.name = name
        self.value = value
        if date is None:
            self.date = datetime.now()
        elif isinstance(date,datetime):
            self.date = date
        else:
            self.date = datetime.fromtimestamp(date)
        
        self.category = category
        self.message = message
        self.level = level
        
    def __repr__(self):
        return "%s %s: %s" % (str(self.date), self.name, self.message)


class MetaContainer(FileBased):
    """
    Class that provides recording of meta-information.
    """

    class MetaDescription(t.IsDescription):
        date = t.Float32Col(pos=0)
        level = t.Int16Col(pos=1)
        name = t.StringCol(50, pos=2)
        category = t.StringCol(50, pos=3)
        message = t.StringCol(255, pos=4)
        value = t.Int32Col(pos=5, dflt=-1)
    
    def __init__(self, data=None, file_name=None,
                 group_name='meta', logger_name="meta", tmpdir=None):
        """
        Enable recording of meta-information in pytables-backed object.
        
        This constructor is built with considerable flexibility
        in the ways it can be called. The 'data' parameter
        serves here as a multi-purpose container that simplifies
        calling the constructor with just one argument.
        
        Args:
            data (...):         multi-purpose container, simplifies one-argument
                                constructor calls. Can stand in for a pytables 
                                Table or File, or a path.
                (None):         An existing meta description table is expected in 
                                self._meta. If none is found, the constructor will 
                                look for a pytables File in self.file and attempt 
                                to create the meta description table in there.
                (str):          An exisiting pytable File at this location will be
                                opened in append mode or one will be created.
                (tables.file.File):
                                A meta description table will be created in this
                                pytables File.
                (tables.table.Table):
                                This pytables Table will be used as meta 
                                description table (if it has the necessary fields
                                date, category, description, name, and value).
            table_name (str):   name of the meta table that is created in
                                pytables file, does not usually need to be 
                                modified
        """
        
        # parse potential unnamed argument
        if data is not None:
            # data is file name
            if type(data) is str or isinstance(data, t.file.File):                
                if file_name is None:
                    file_name = data
                    data = None
        
        if file_name is not None and isinstance(file_name, str):
            file_name = os.path.expanduser(file_name)
        
        FileBased.__init__(self, file_name, tmpdir=tmpdir)
        
        try:
            self._meta = self.file.get_node("/" + group_name + "/main")
            self._meta_values = self.file.get_node("/" + group_name + "/values")
        except NoSuchNodeError:
            # create reads group
            group = self.file.create_group("/", group_name, 'Meta info group',
                                           filters=t.Filters(complib="blosc",
                                                             complevel=2, shuffle=True))
            self._meta = self.file.create_table(group, "main",
                                                MetaContainer.MetaDescription)
            self._meta_values = self.file.create_vlarray(group, "values", t.ObjectAtom(),
                                                         expectedrows=100)
        
        self._meta_logger = logging.getLogger(logger_name)
        self._meta_level_debug = 0
        self._meta_level_info = 1
        self._meta_level_warn = 2
        self._meta_level_error = 3
    
    def log_info(self, message, save=True):
        if save:
            self.add_meta(message, level=self._meta_level_info)
        self._meta_logger.info(message)
    
    def log_debug(self, message, save=True):
        if save:
            self.add_meta(message, level=self._meta_level_debug)
        self._meta_logger.debug(message)
    
    def log_warn(self, message, save=True):
        if save:
            self.add_meta(message, level=self._meta_level_warn)
        self._meta_logger.warn(message)
    
    def log_error(self, message, save=True):
        if save:
            self.add_meta(message, level=self._meta_level_error)
        self._meta_logger.error(message)
    
    def add_meta(self, message='', name='', value=None, level=0, category=''):
        row = self._meta.row
        
        row['date'] = time.time()
        row['name'] = name
        row['category'] = category
        row['message'] = message
        row['level'] = level
        if value is not None:
            value_ix = len(self._meta_values)
            row['value'] = value_ix
            self._meta_values.append(value)
            self._meta_values.flush()
        row.append()
        self._meta.flush()
    
    def get_meta(self, key):

        if isinstance(key, int):
            row = self._meta[key]
            value_ix = row["value"]
            if value_ix >= 0:
                value = self._meta_values[value_ix]
            else:
                value = None
            return Meta(date=row["date"], message=row["message"],
                        name=row["name"], value=value, level=row["level"],
                        category=row["category"])
        if isinstance(key, str):
            metas = []
            for row in self._meta.where("name == '%s'" % key):
                value_ix = row["value"]
                value = self._meta_values[value_ix]
                meta = Meta(date=row["date"], message=row["message"],
                            name=row["name"], value=value, level=row["level"],
                            category=row["category"])
                metas.append(meta)
            return metas
        
    def meta_history(self, n=20):
        """
        Return a list of recent meta_info.
        """
        
        n = min(len(self._meta)-1, n)
        history = []
        while n > 0:
            history.append(self.get_meta(n))
            n -= 1
        
        return history
