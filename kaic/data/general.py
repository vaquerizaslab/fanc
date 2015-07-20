'''
Provide basic convenience data types.

Most classes are based on pytables and hdf5 dictionaries,
allowing on-disk storage and therefore processing of large
files. Other features include indexing and querying.
'''

import tables as t
from kaic.tools.files import create_or_open_pytables_file, is_hdf5_file,\
    random_name
from __builtin__ import isinstance
import random # @UnusedImport
import string # @UnusedImport
import numpy as np
import warnings
import os.path

import logging
from tables.exceptions import NoSuchNodeError
logging.basicConfig(level=logging.INFO)



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

def _structured_array_to_table_type(a, rownames=None, colnames=None):
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
        return TableRow(a[0], rowname=rownames[0], colnames=colnames)

    # is column
    if len(a) > 1 and len(a[0]) == 1:
        l = []
        for row in a:
            l.append(row[0])
        return TableCol(np.array(l), colname=colnames[0], rownames=rownames)

    # is matrix or nothing
    return TableArray(a, colnames=colnames, rownames=rownames)


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
        colnames = a.dtype.fields.keys()
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
                    for i in range(0,len(previous)):
                        previous[i] = types[i](previous[i])
    
                data.append(previous)
                
            previous = fields
            line = f.readline().rstrip()
            i += 1
                
        if (header is not None) and (len(previous) == len(header) + 1):
            rownames.append(previous[0])
            previous = previous[1:]

        if types is not None:
            for i in range(0,len(previous)):
                previous[i] = types[i](previous[i])
        data.append(previous)

    return data, header, rownames

class TableRow(tuple):
    """
    In-memory row of a Table.
    
    Provides convenient access to fields using the Table column name
    or field index. Extends a Python tuple, so all expected tuple
    functionality is there.
    """
    
    def __new__(cls, t, rowname=None, colnames=None):
        obj = super(TableRow, cls).__new__(cls, tuple(t))
        obj.rowname = rowname
        obj.colnames = colnames
        obj._colnames_dict = {colnames[i]:i for i in range(0,len(colnames))}
        return obj

    def __init__(self, t, rowname='0', colnames=None):
        pass

    # BUG FIX in python
    def __getslice__(self, start, stop) :
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
    def __getslice__(self, start, stop) :
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
            obj._colnames_dict = {colnames[i]:i for i in range(0,len(colnames))}
        obj.rownames=rownames
        if rownames is not None:
            obj._rownames_dict = {rownames[i]:i for i in range(0,len(rownames))}

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

    def __getslice__(self, start, stop) :
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
                    dtypes.append((name,dtype[key[i]]))
                except KeyError:
                    dtypes = dtype


            a = np.zeros((len(self),), dtype=dtypes)
            a[:] = [tuple([super(TableArray, self).__getitem__(i)[k] for k in key]) for i in range(0,len(self))]

            return _structured_array_to_table_type(a, colnames=cn, rownames=self.rownames)

        except TypeError:
            pass


        try:
            dtypes = dtype[key]
        except KeyError:
            dtypes=dtype


        b = [(super(TableArray, self).__getitem__(i)[key],) for i in range(0,len(self))]
        name = self.colnames[key]
        mydt = np.dtype([(name,str(dtypes))])
        a = np.array(b,dtype=mydt)

        return _structured_array_to_table_type(a, colnames=name, rownames=self.rownames)



    def dim(self):
        """
        Return the number of rows and columns as a tuple
        """
        
        if len(self) > 0:
            return (len(self), len(super(TableArray, self).__getitem__(0)))
        return (len(self),0)


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
                       table_name='table'):
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

        # open file or keep in memory
        if file_name == None:
            file_name = random_name()
            self.file = create_or_open_pytables_file(file_name, inMemory=True)
        else:
            self.file = create_or_open_pytables_file(file_name, inMemory=False)


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

            #rows
            if rownames is not None:
                nrows = max(nrows, len(rownames))
            else:
                rownames = []

            if data is not None:
                nrows = max(nrows, len(data))

            for i in range(len(rownames), nrows):
                rownames.append(str(i))

            self._table_name = table_name

            columns_dict = {}
            columns_dict[self._rowname_field] = t.StringCol(50,pos=0) # @UndefinedVariable

            for i in range(0,len(colnames)):
                columns_dict[colnames[i]] = _convert_to_tables_type(col_types[i], pos=i+1)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                self._table = self.file.create_table("/", table_name, columns_dict)
        else:
            self._table = self.file.get_node('/' + table_name)

        # clean up
        self._table.flush()

        # add any potential data
        if data is not None:
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
        except (TypeError,KeyError,AttributeError,IndexError, ValueError), e:
            print e
            raise TypeError("d is not a dictionary")

        row.append()

        if flush:
            self._table.flush()


    def append(self, data, rownames=None):
        """
        Append data to the table.
        
        Args:
        data (...):      Data can have several formats:
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
        return (len(self), len(self._table.colnames)-1)
    

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

    def _get_rows(self, key):

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

                t = self._get_rows(key[0])

                return t._get_cols(key[1])

            else:
                raise KeyError("Unrecognized key " + str(key))

        return self._get_rows(key)

    def __getattr__(self, name):
        if name in self.colnames:
            return self._get_cols(name)
        raise AttributeError


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

        return _structured_array_to_table_type(a, rownames=rn, colnames=cn)



class Mask(object):
    """
    Class providing Mask details.
    """
    
    def __init__(self, ix, name, description=''):
        self.ix = ix
        self.name = name
        self.description = description

class Maskable(object):
    """
    
    """
    
    class MaskDescription(t.IsDescription):
        ix = t.Int16Col(pos=0)
        name = t.StringCol(50,pos=1)
        description = t.StringCol(255,pos=2)
    
    
    def __init__(self, data=None, table_name="mask"):
        super(Maskable, self).__init__()
        
        if not hasattr(self, 'file') or self.file is None:

            if data is None or not type(data) == str:
                self.file = create_or_open_pytables_file(random_name(), inMemory=True)
            else:
                self.file = create_or_open_pytables_file(data, inMemory=False)
        
        else:
            if not type(self.file) == t.file.File:
                raise ValueError("Object has file attribute, but it is not a pytables File object")
        
        
        if hasattr(self, '_mask') and type(self._mask) == t.table.Table:
            logging.info("Found existing _mask attribute")
            return
        
        if hasattr(self, '_mask'):
            raise ValueError("Object already has a _mask attribute, but it is not from a Maskable object")
        
        try:
            self._mask = self.file.get_node('/' + table_name)
        except NoSuchNodeError:
            self._mask = self.file.create_table("/", table_name, Maskable.MaskDescription)
            row = self._mask.row
            row['ix'] = 1
            row['name'] = 'default'
            row['description'] = 'Default mask'
            row.append()
            self._mask.flush()
        
    def add_mask_description(self, name, description):
        ix = 1
        if len(self._mask) > 0:
            ix = self._mask[-1][0]*2
        row = self._mask.row
        row['ix'] = ix
        row['name'] = name
        row['description'] = description
        row.append()
        self._mask.flush()
        
        return int(ix)
    
    def get_mask(self, key):
        if type(key) == int:
            condition = "ix == %d" % key
            res = [x.fetch_all_fields() for x in self._mask.where(condition)]
        else:
            res = [x.fetch_all_fields() for x in self._mask.where("name == %s" % str(key))]

        if len(res) == 1:
            return Mask(res[0][0], res[0][1], res[0][2])
        return None
    
    def get_masks(self, key):
        key_copy = key

        last = int(self._mask[-1][0])

        masks = []
        while last > 1:
            if key_copy - last >= 0:
                mask = self.get_mask(last)
                if mask is not None:
                    masks.append(mask)
                key_copy -= last
            last = last/2
        if key_copy > 0:
            masks.append(self.get_mask(0))
        
        return masks
    
    
class MetaContainer(object):
    class MetaDescription(t.IsDescription):
        date = t.Time32Col(pos=0)
        category = t.StringCol(50,pos=1)
        name = t.StringCol(255,pos=2)
        value = t.StringCol(255, pos=3)
        
    def __init__(self, data=None, table_name="meta"):
        super(MetaContainer, self).__init__()
        
        if not hasattr(self, 'file') or self.file is None:

            if data is None or not type(data) == str:
                self.file = create_or_open_pytables_file(random_name(), inMemory=True)
            else:
                self.file = create_or_open_pytables_file(data, inMemory=False)
        
        else:
            if not type(self.file) == t.file.File:
                raise ValueError("Object has file attribute, but it is not a pytables File object")
        
        
        if hasattr(self, '_meta') and type(self._meta) == t.table.Table:
            logging.info("Found existing _meta attribute")
            return
        
        if hasattr(self, '_meta'):
            raise ValueError("Object already has a _meta attribute, but it is not from a MetaContainer object")
        
        try:
            self._meta = self.file.get_node('/' + table_name)
        except NoSuchNodeError:
            self._meta = self.file.create_table("/", table_name, MetaContainer.MetaDescription)
            self._meta.flush()
            
    