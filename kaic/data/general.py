'''
Created on Jun 8, 2015

@author: kkruse1
'''

import tables as t
from kaic.tools.files import create_or_open_pytables_file
from __builtin__ import isinstance
import random
import string
from collections import OrderedDict
import numpy as np
import warnings

_typemap = {
    str: t.StringCol, # @UndefinedVariable
    int: t.Int32Col, # @UndefinedVariable
    float: t.Float32Col, # @UndefinedVariable
    bool: t.BoolCol # @UndefinedVariable
}
_inv_typemap = {v: k for k, v in _typemap.items()}

def _convert_to_tables_type(col_type, pos=None):
    if isinstance(col_type,t.Col):
        return col_type
    if col_type in _typemap:
        if col_type is str:
            return _typemap[col_type](255,pos=pos)
        return _typemap[col_type](pos=pos)
    raise ValueError("Unknown column type " + str(col_type))
    




        
class Table(object):
    
    def __init__(self, file_name=None, name=None,
                       ncols=0, nrows=0,
                       col_names=None, row_names=None,
                       col_types=None, default_type=str,
                       table_name='table', data=None):
        
        # open file if exists
        # open file or keep in memory
        in_memory  = False
        if file_name == None:
            in_memory = True
        
        if in_memory:
            rs = ''.join(random.SystemRandom().choice(string.uppercase + string.digits) for _ in xrange(6))  # @UndefinedVariable
            file_name = rs
        
        self.file = create_or_open_pytables_file(file_name, inMemory=in_memory)
        
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
            name = name if name else self._table.attrs.name
        
        # clean up
        self._table.attrs.name = name if name else file_name
        self.name = self._table.attrs.name
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
        self.rownames(row_names)
        
        
        
        
        
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
    
    def save(self,file_name):
        t = Table.from_structured_array(self._table[:], name=self._table.attrs.name, file_name=file_name, table_name=self._table_name, row_names=True)
        self.file.close()
        self.file = t.file
        self._table = t._table
        
    
    def rownames(self, names=None):
        if names is not None:
            if not len(names) == len(self):
                raise ValueError("names must be the same length as table")
            i = 0
            for row in self._table.iterrows():
                row[self._rowname_field] = names[i]
                row.update()
                i += 1
    
            self._table.flush()
        return self._table[:]['_rowname']
            
    def _id2ix(self, key):
        try:
            ix = self._table.colnames.index(key)
        except ValueError:
            raise KeyError("%s is not a valid column name" % key)
        
        return ix
    
    def _ix2id(self, ix):
        try:
            key = self._table.colnames[ix]
        except IndexError:
            raise KeyError("%d is not a column index (max %d)" % (ix, len(self._table.colnames)))
        return key
    
    def _append_row_list(self, l, flush=True, rowname=None):
        row = self._table.row
                
        try:
            
            if len(l) != len(self._table.colnames):
                if rowname is not None:
                    row[self._rowname_field] = rowname
                else:
                    row[self._rowname_field] = len(self)
                for i in range(1,len(self._table.colnames)):
                    row[self._table.colnames[i]] = l[i-1]
            else:
                for i in range(0,len(self._table.colnames)):
                    row[self._table.colnames[i]] = l[i]
        
        except (TypeError, KeyError, IndexError, ValueError), e:
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
        
    
#     def name(self, name=None):
#         if name is not None:
#             self._table.attrs.name = name
#         return self._table.attrs.name
    
    def __len__(self):
        return len(self._table)
    
    def dim(self):
        return (len(self), len(self._table.colnames))
    
    def __repr__(self, *args, **kwargs):
        # find out maximum column width
        show = 9
        m = 3 if show < len(self) else 0
        
        max_col_width = [max(len(str(x)),m) for x in self._table.colnames]
        max_col_width[0] = 0
        for i in range(0,min(show,len(self))):
            for j in range(0,len(max_col_width)):
                max_col_width[j] = max(max_col_width[j], len(str(self._table[i][j])))
        
        r = 'Name: %s\n' % self._table.attrs.name
        r += ' ' * max_col_width[0] + ' '
        for j in range(1,len(self._table.colnames)):
            r += str(self._table.colnames[j]).rjust(max_col_width[j]) + " "
        r += "\n"
            
        
        for i in range(0,min(show,len(self))):
            for j in range(0,len(max_col_width)):
                r += str(self._table[i][j]).rjust(max_col_width[j]) + " "
            if i < len(self)-1:
                r += "\n"
        
        r += ' ' * max_col_width[0] + ' '
        if show < len(self):
            for j in range(1,len(self._table.colnames)):
                r += '...'.rjust(max_col_width[j]) + ' '
            r += '\n%d rows' % len(self)
            
        
        return r
    
    
    def _get_rows(self, key, as_table=None):
        if type(key) is str or type(key) is list:
            table_dtypes = self._table[0:0].dtype
            dtypes = []
            for name in self._table.colnames:
                dt = str(table_dtypes[name])
                dtypes.append((name,dt))
            
            if type(key) is str:
                condition = "%s == '%s'" % (self._rowname_field,key)
                l = [x.fetch_all_fields() for x in self._table.where(condition)]
            else:
                
                l = []
                if len(key) == 0:
                    return Table()
                elif type(key[0]) is str:
                    for k in key:
                        l = l + [x.fetch_all_fields() for x in self._table.where("%s == '%s'" % (self._rowname_field,k))]
                else:
                    for i in key:
                        l.append(self._table[i])

            a = np.zeros((len(l),), dtype=dtypes)
            a[:] = l
        else:
            a = self._table[key]
        
        if as_table is True:
            return Table.from_structured_array(a,name=self.name,row_names=True)
        if as_table is False:
            return a
        
        
        if type(a) is np.void:
            return a
        if len(a) == 1:
            if len(a[0]) == 1:
                return a[0][0]
            return a[0]
            
        return Table.from_structured_array(a,name=self.name,row_names=True)
    
    def _get_cols(self, key):
        rowname_dt = self._table[0:0].dtype[self._rowname_field]
        
        if type(key) is int:
            key = self._table.colnames[key+1]
                
        if type(key) is str:
            dt = self._table[0:0].dtype[key]
            b = [(x[self._rowname_field],x[key],) for x in self._table.iterrows()]
            mydt = np.dtype([(self._rowname_field,str(rowname_dt)),(key,str(dt))])
            a = np.array(b,dtype=mydt)

            if len(a) == 1:
                if len(a[0]) == 2:
                    return a[0][1]
                return a[0]
            return Table.from_structured_array(a,name=self.name,row_names=True)
        
        if type(key) is slice:
            start = 0 if key.start is None else key.start
            stop = len(self._table.colnames)-1 if key.stop is None else (min(self.dim()[1],key.stop))
            step = 1 if key.step is None else key.step
            
            key = [self._table.colnames[i+1] for i in range(start,stop,step)]
            
        try:
            dtypes = [(self._rowname_field,rowname_dt)]
            keys = [self._rowname_field]
            for k in key:
                if type(k) is int:
                    k = self._table.colnames[k+1]
                keys.append(k)
                dt = self._table[0:0].dtype[k]
                dtypes.append((k,dt))
            a = np.zeros((len(self),), dtype=dtypes)
            a[:] = [tuple([x[k] for k in keys]) for x in self._table.iterrows()]
            
            if len(a) == 1:
                if len(a[0]) == 2:
                    return a[0][1]
                return a[0]
        except IndexError:
            raise IndexError("Cannot retrieve column with value " + str(key))
            
        return Table.from_structured_array(a,name=self.name,row_names=True)
            
    
    
    def __getitem__(self, key):
        
        if type(key) is tuple:
            t = self._get_rows(key[0], as_table=True)
            return t._get_cols(key[1])
        
        return self._get_rows(key)
        
    
    
    def as_array(self, include_rownames=False):
        if include_rownames:
            return self._table[:]
        
        dtypes = []
        keys = self._table.colnames[1:]
        for k in keys:
            dt = self._table[0:0].dtype[k]
            dtypes.append((k,dt))
        if len(dtypes) == 0:
            a = np.zeros(0)
        elif len(dtypes) > 1:
            a = np.zeros((len(self),), dtype=dtypes)
            a[:] = [tuple([x[k] for k in keys]) for x in self._table.iterrows()]
        else:
            a = np.array([x[keys[0]] for x in self._table.iterrows()])
        
        return a
        
        
        
        
        
        
        
        
        
    #def select(self, query):