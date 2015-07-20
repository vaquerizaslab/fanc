'''
Created on Jun 8, 2015

@author: kkruse1
'''

import tables as t
from kaic.tools.files import create_or_open_pytables_file, is_hdf5_file,\
    random_name
from __builtin__ import isinstance
import random # @UnusedImport
import string # @UnusedImport
from collections import OrderedDict
import numpy as np
import warnings
import os.path

import logging
from tables.exceptions import NoSuchNodeError
logging.basicConfig(level=logging.INFO)



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
_inv_typemap = {v: k for k, v in _typemap.items()}

def _convert_to_tables_type(col_type, pos=None):
    if isinstance(col_type,t.Col):
        return col_type
    if col_type in _typemap:
        if col_type is str or col_type is np.string_:
            return _typemap[col_type](255,pos=pos)
        return _typemap[col_type](pos=pos)
    raise ValueError("Unknown column type " + str(col_type))

def _structured_array_to_table_type(a, rownames=None, colnames=None):

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
    data = []
    rownames = []
    header = None
    with open(file_name,'r') as f:
        # header(?)
        if has_header:
            header_line = f.readline().rstrip()
            header = str.split(header_line, sep)

        line = f.readline().rstrip()
        while line != '':
            fields = str.split(line, sep)
            if (header is not None) and (len(fields) == len(header) + 1):
                rownames.append(fields[0])
                fields = fields[1:]

            if types is not None:
                for i in range(0,len(fields)):
                    fields[i] = types[i](fields[i])

            data.append(fields)
            line = f.readline().rstrip()

    return data, header, rownames

class TableRow(tuple):
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
        if len(self) > 0:
            return (len(self), len(super(TableArray, self).__getitem__(0)))
        return (len(self),0)


class Table(object):

    def __init__(self, data=None, file_name=None,
                       ncols=0, nrows=0,
                       col_names=None, row_names=None,
                       col_types=None, default_type=str,
                       table_name='table'):

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
                    data, tmp_col_names, tmp_row_names = _file_to_data(data, has_header=True, types=col_types)
                    if col_names is None:
                        col_names = tmp_col_names
                    if row_names is None:
                        row_names = tmp_row_names
            else:
                data, c, ty = _to_list_and_names(data)
                if col_names is None:
                    col_names = c
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
            if col_names is not None:
                ncols = max(ncols, len(col_names))
            else:
                col_names = []

            for i in range(len(col_names), ncols):
                col_names.append(str(i))

            if col_types is None:
                col_types = []

            for i in range(len(col_types), ncols):
                col_types.append(default_type)

            #rows
            if row_names is not None:
                nrows = max(nrows, len(row_names))
            else:
                row_names = []

            if data is not None:
                nrows = max(nrows, len(data))

            for i in range(len(row_names), nrows):
                row_names.append(str(i))

            self._table_name = table_name

            columns_dict = OrderedDict()
            columns_dict[self._rowname_field] = t.StringCol(50,pos=0) # @UndefinedVariable

            for i in range(0,len(col_names)):
                columns_dict[col_names[i]] = _convert_to_tables_type(col_types[i], pos=i+1)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                self._table = self.file.create_table("/", table_name, columns_dict)
        else:
            self._table = self.file.get_node('/' + table_name)

        # clean up
        self._table.flush()

        # add potential data
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
        self.set_rownames(row_names)




    @classmethod
    def from_structured_array(cls, a, name=None, file_name=None, table_name='table', row_names=None):
        try:
            col_names = a.dtype.names
        except AttributeError:
            raise ValueError("Input must be a structured array")

        # "ghetto" conversion of numpy file types to python
        try:
            col_types = [type(a.dtype.fields[x][0].type(0).item()) for x in col_names]
        except TypeError:
            col_types = [type(a.dtype.type(0).item())]

        if row_names is True:
            try:
                row_names = [x[0] for x in a]
            except IndexError:
                row_names = a[0]

        t = cls(col_names=col_names, col_types=col_types, file_name=file_name, name=name, table_name=table_name)
        t.append(a, row_names)
        return t


    def save(self,file_name, table_name='table'):
        self.file.copy_file(file_name)
        self.file.close()
        self.file = create_or_open_pytables_file(file_name)
        self._table = self.file.get_node('/' + table_name)
        
    def export(self, file_name, sep="\t", include_colnames=True, include_rownames=True):
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
        self.file.close()
    
    @property
    def colnames(self):
        return self._table.colnames[1:]
        
    @property
    def rownames(self):
        return self._table[:]['_rowname']
    
    def set_rownames(self, names):
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
        row = self._table.row

        try:
            if not self._rowname_field in d:
                d[self._rowname_field] = len(self)
            for name in d:
                row[name] = d[name]
        except (TypeError,KeyError,AttributeError,IndexError, ValueError):
            raise TypeError("d is not a dictionary")

        row.append()

        if flush:
            self._table.flush()


    def append(self, data, rownames=None):
        # data is file
        if type(data) is str:
            self.append_file(data)
            return

        # data is a list of lists?
        try:
            i = 0
            for l in data:
                if rownames is not None:
                    self._append_row_list(l, flush=False, rowname=rownames[i])
                else:
                    self._append_row_list(l, flush=False)
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
                    self._append_row_dict(d, flush=False)
                i += 1
            self._table.flush()
            return
        except TypeError:
            pass

        # data is a list?
        try:
            self._append_row_list(data, flush=True, rowname=rownames)
            return
        except TypeError:
            pass

        # data is a dictionary
        try:
            self._append_row_dict(data, flush=True, rowname=rownames)
            return
        except TypeError:
            raise ValueError("Data type unsupported")



    def __len__(self):
        return len(self._table)

    def dim(self):
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


    def as_array(self, include_rownames=False):

        return self[:]


    def where(self, query):
        
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
    def __init__(self, ix, name, description=''):
        self.ix = ix
        self.name = name
        self.description = description

class Maskable(object):
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
            
    