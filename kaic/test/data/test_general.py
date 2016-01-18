import tables as t
from tables.registry import class_id_dict
import numpy as np
import pytest
from kaic.data.general import Table, _to_list_and_names, TableRow, TableCol,\
    TableArray, _convert_to_tables_type, _structured_array_to_table_type,\
    _file_to_data, Mask, Maskable, MaskedTable, MaskFilter, MetaContainer,\
    FileBased
from __builtin__ import classmethod
import os
from kaic.tools.files import create_or_open_pytables_file

class TestSupport:
    
    def test_list_conversion(self):
        # structured array
        x = np.zeros((5,),dtype=('i4,f4,a10'))
        x[:] = [(1,2.,'Hello'),(2,3.,"World"),(3,4.,'this'),(4,5.,"is"),(5,6.,'me')]
        l, n, _ = _to_list_and_names(x)
        assert l[0] == [1,2.,'Hello']
        assert l[4] == [5,6.,'me']
        assert n == ['f0','f1','f2']
        #assert t == [int, float, str]
        
        
        # 2D array
        x = np.zeros((5,3))
        l, n, _ = _to_list_and_names(x)
        assert l[0] == [0.,0.,0.]
        assert l[4] == [0.,0.,0.]
        assert n == ['0','1','2']
        
        # 2D list
        x = [(1,2.,'Hello'),(2,3.,"World"),(3,4.,'this'),(4,5.,"is"),(5,6.,'me')]
        l, n, _ = _to_list_and_names(x)
        assert l[0] == [1,2.,'Hello']
        assert l[4] == [5,6.,'me']
        assert n == ['0','1','2']
        
        # 1D array
        x = np.array([1,2,3])
        l, n, _ = _to_list_and_names(x)
        assert l[0] == [1,2,3]
        assert n == ['0','1','2']
        
        # 1D list
        x = [1,2,3]
        l, n, _ = _to_list_and_names(x)
        assert l[0] == [1,2,3]
        assert n == ['0','1','2']
        
        # scalar
        x = 1
        l, n, _ = _to_list_and_names(x)
        assert l[0] == [1]
        assert n == ['0']
        
    def test_typemap(self):
        assert _convert_to_tables_type(str) == t.StringCol(255)
        assert _convert_to_tables_type(int) == t.Int32Col()
        assert _convert_to_tables_type(float) == t.Float32Col()
        assert _convert_to_tables_type(bool) == t.BoolCol()
        assert _convert_to_tables_type(np.string_) == t.StringCol(255)
        assert _convert_to_tables_type(np.int32) == t.Int32Col()
        assert _convert_to_tables_type(np.float32) == t.Float32Col()
        assert _convert_to_tables_type(np.bool_) == t.BoolCol()
        
        assert _convert_to_tables_type(t.BoolCol()) == t.BoolCol()
        with pytest.raises(ValueError):
            _convert_to_tables_type(Table)
    
    def test_table_conversion(self):
        x = np.zeros((5,),dtype=('i4,f4,a10'))
        x[:] = [(1,2.,'Hello'),(2,3.,"World"),(3,4.,'this'),(4,5.,"is"),(5,6.,'me')]
        a = _structured_array_to_table_type(x, rownames=['a','b','c','d','e'], colnames=['A','B','C'])
        assert type(a) == TableArray
        
        x = np.zeros((1,),dtype=('i4,f4,a10'))
        x[:] = [(1,2.,'Hello')]
        a = _structured_array_to_table_type(x, rownames=['a'], colnames=['A','B','C'])
        assert type(a) == TableRow
        
        x = np.zeros((5,),dtype=[('a','i4')])
        x[:] = [(1,),(2,),(3,),(4,),(5,)]
        a = _structured_array_to_table_type(x, rownames=['a','b','c','d','e'], colnames=['A'])
        assert type(a) == TableCol
        
    def test_file_read(self):
        this_dir = os.path.dirname(os.path.realpath(__file__))
        tsv1_file = this_dir + "/test_general/tsv.test.txt"
        tsv2_file = this_dir + "/test_general/tsv.test2.txt"
        
        data, colnames, rownames = _file_to_data(tsv1_file)
        assert np.array_equal(data[0], ['A','B','C'])
        assert np.array_equal(data[1], ['a','1','1.'])
        assert np.array_equal(data[2], ['b','2','2.'])
        assert np.array_equal(colnames, [])
        assert np.array_equal(rownames, [])
        
        data, colnames, rownames = _file_to_data(tsv1_file, has_header=True)
        assert np.array_equal(data[0], ['a','1','1.'])
        assert np.array_equal(data[1], ['b','2','2.'])
        assert np.array_equal(colnames, ['A','B','C'])
        assert np.array_equal(rownames, [])
        
        data, colnames, rownames = _file_to_data(tsv1_file, has_header=True, types=[str,int,float])
        assert np.array_equal(data[0], ['a',1,1.])
        assert np.array_equal(data[1], ['b',2,2.])
        assert np.array_equal(colnames, ['A','B','C'])
        assert np.array_equal(rownames, [])
        
        # auto-detect header
        data, colnames, rownames = _file_to_data(tsv2_file)
        assert np.array_equal(data[0], ['a','1','1.'])
        assert np.array_equal(data[1], ['b','2','2.'])
        assert np.array_equal(colnames, ['A','B','C'])
        assert np.array_equal(rownames, ['a','b'])
        
        data, colnames, rownames = _file_to_data(tsv2_file, has_header=True)
        assert np.array_equal(data[0], ['a','1','1.'])
        assert np.array_equal(data[1], ['b','2','2.'])
        assert np.array_equal(colnames, ['A','B','C'])
        assert np.array_equal(rownames, ['a','b'])
        
        data, colnames, rownames = _file_to_data(tsv2_file, has_header=True, types=[str,int,float])
        assert np.array_equal(data[0], ['a',1,1.])
        assert np.array_equal(data[1], ['b',2,2.])
        assert np.array_equal(colnames, ['A','B','C'])
        assert np.array_equal(rownames, ['a','b'])
        
        
    
class TestTableRow:
    
    @classmethod
    def setup_method(self, method):
        data = (1, 2., 'Hello')
        self.row = TableRow(data, rowname='a', colnames=['A','B','C'])
        
    def test_get_item(self):
        # index
        assert self.row[0] == 1
        assert self.row[1] == 2.
        assert self.row[2] == 'Hello'
        
        # name
        assert self.row['A'] == 1
        assert self.row['B'] == 2.
        assert self.row['C'] == 'Hello'
        
        # index list
        x = self.row[[0,2]]
        assert x == (1,'Hello')
        assert np.array_equal(x.colnames, ['A','C'])
        
        # name list
        x = self.row[['A','C']]
        assert x == (1,'Hello')
        assert np.array_equal(x.colnames, ['A','C'])
        
        # mixed list
        x = self.row[[0,'C']]
        assert x == (1,'Hello')
        assert np.array_equal(x.colnames, ['A','C'])
        
        # slice
        x = self.row[0:2]
        assert x == (1,2.)
        assert np.array_equal(x.colnames, ['A','B'])
        
        x = self.row[0:1]
        assert x == (1,)
        assert np.array_equal(x.colnames, ['A'])
        
    def test_get_attr(self):
        assert self.row.A == 1
        assert self.row.B == 2.
        assert self.row.C == 'Hello'
        
        
class TestTableCol:
    
    @classmethod
    def setup_method(self, method):
        data = np.array([1, 2, 3])
        self.col = TableCol(data, colname='a', rownames=['A','B','C'])
        
    def test_get_item(self):
        # index
        assert self.col[0] == 1
        assert self.col[1] == 2
        assert self.col[2] == 3
        
        # name
        assert self.col['A'] == 1
        assert self.col['B'] == 2
        assert self.col['C'] == 3
        
        # index list
        x = self.col[[0,2]]
        assert np.array_equal(x, [1,3])
        assert np.array_equal(x.rownames, ['A','C'])
        
        # name list
        x = self.col[['A','C']]
        assert np.array_equal(x, [1,3])
        assert np.array_equal(x.rownames, ['A','C'])
        
        # mixed list
        x = self.col[[0,'C']]
        assert np.array_equal(x, [1,3])
        assert np.array_equal(x.rownames, ['A','C'])
        
        # slice
        x = self.col[0:2]
        assert np.array_equal(x, [1,2])
        assert np.array_equal(x.rownames, ['A','B'])
        
        x = self.col[0:1]
        assert np.array_equal(x, [1])
        assert np.array_equal(x.rownames, ['A'])
        
        
class TestTableArray:
    
    @classmethod
    def setup_method(self, method):
        x = np.zeros((5,),dtype=('i4,f4,a10'))
        x[:] = [(1,2.,'Hello'),(2,3.,"World"),(3,4.,'this'),(4,5.,"is"),(5,6.,'me')]
        
        
        self.table = TableArray(x, rownames=['a','b','c','d','e'], colnames=['A','B','C'])
    
    def test_row_selection(self):
        
        # single rows by ix
        assert self.table._get_rows(0) == (1,2.,'Hello')
        assert self.table._get_rows(1) == (2,3.,'World')
        assert self.table._get_rows(2) == (3,4.,'this')
        assert self.table._get_rows(3) == (4,5.,'is')
        assert self.table._get_rows(4) == (5,6.,'me')
        
        # single rows by name
        assert self.table._get_rows('a') == (1,2.,'Hello')
        assert self.table._get_rows('b') == (2,3.,'World')
        assert self.table._get_rows('c') == (3,4.,'this')
        assert self.table._get_rows('d') == (4,5.,'is')
        assert self.table._get_rows('e') == (5,6.,'me')
        
        # single rows by list
        assert self.table._get_rows([0]) == (1,2.,'Hello')
        assert self.table._get_rows([1]) == (2,3.,'World')
        assert self.table._get_rows([2]) == (3,4.,'this')
        assert self.table._get_rows([3]) == (4,5.,'is')
        assert self.table._get_rows([4]) == (5,6.,'me')
        
        # multiple rows by list
        x = self.table._get_rows([0,1])
        assert tuple(x[0]) == (1,2.,'Hello')
        assert tuple(x[1]) == (2,3.,'World')
        assert np.array_equal(x.colnames, ['A', 'B', 'C'])
        assert np.array_equal(x.rownames, ['a','b'])
        
        x = self.table._get_rows([4,2])
        assert tuple(x[0]) == (5,6.,'me')
        assert tuple(x[1]) == (3,4.,'this')
        assert np.array_equal(x.colnames, ['A', 'B', 'C'])
        assert np.array_equal(x.rownames, ['e','c'])
        
        # slice
        x = self.table._get_rows(slice(0,2,1))
        assert tuple(x[0]) == (1,2.,'Hello')
        assert tuple(x[1]) == (2,3.,'World')
        assert np.array_equal(x.colnames, ['A', 'B', 'C'])
        assert np.array_equal(x.rownames, ['a','b'])
        
        x = self.table._get_rows(slice(0,1,1))
        assert x == (1,2.,'Hello')
        assert np.array_equal(x.colnames, ['A', 'B', 'C'])
        assert x.rowname == 'a'
        
    def test_col_selection(self):
        
        # single cols
        assert np.array_equal(self.table._get_cols(0),(1,2,3,4,5))
        assert np.array_equal(self.table._get_cols('A'),(1,2,3,4,5))
        
        assert np.array_equal(self.table._get_cols(1),(2.,3.,4.,5.,6.))
        assert np.array_equal(self.table._get_cols('B'),(2.,3.,4.,5.,6.))
        
        assert np.array_equal(self.table._get_cols(2),('Hello','World','this','is','me'))
        assert np.array_equal(self.table._get_cols('C'),('Hello','World','this','is','me'))
        
        x = self.table._get_cols(0)
        assert np.array_equal(x.colname, 'A')
        assert np.array_equal(x.rownames, ['a','b','c','d','e'])
        
        
        
        
        # single rows by list
        assert np.array_equal(self.table._get_cols([0]),(1,2,3,4,5))
        assert np.array_equal(self.table._get_cols([1]),(2.,3.,4.,5.,6.))
        assert np.array_equal(self.table._get_cols([2]),('Hello','World','this','is','me'))
        
        
        # multiple rows by list
        # multiple cols
        x = self.table._get_cols([1,2])
        assert x[0][0] == 2.
        assert x[4][1] == 'me'
        assert len(x[0]) == 2
        assert np.array_equal(x.colnames, ['B', 'C'])
        assert np.array_equal(x.rownames, ['a','b','c','d','e'])
        
        # (reverse)
        x = self.table._get_cols([2,1])
        assert x[0][1] == 2.
        assert x[4][0] == 'me'
        assert len(x[0]) == 2
        assert np.array_equal(x.colnames, ['C','B'])
        assert np.array_equal(x.rownames, ['a','b','c','d','e'])
        
        x = self.table._get_cols(['B','C'])
        assert x[0][0] == 2.
        assert x[4][1] == 'me'
        assert len(x[0]) == 2
        assert np.array_equal(x.colnames, ['B', 'C'])
        assert np.array_equal(x.rownames, ['a','b','c','d','e'])
        
        
        # slice
        x = self.table._get_cols(slice(1,3,1))
        assert x[0][0] == 2.
        assert x[4][1] == 'me'
        assert len(x[0]) == 2
        assert np.array_equal(x.colnames, ['B', 'C'])
        assert np.array_equal(x.rownames, ['a','b','c','d','e'])

        

class TestTable:
    
    
    @classmethod
    def setup_method(self, method):
        x = np.zeros((5,),dtype=('i4,f4,a10'))
        x[:] = [(1,2.,'Hello'),(2,3.,"World"),(3,4.,'this'),(4,5.,"is"),(5,6.,'me')]
        
        
        self.table = Table(x, rownames=['a','b','c','d','e'], colnames=['A','B','C'])

    
    def test_intialize(self):
        # from tsv
        # no rownames
        current_dir = os.path.dirname(os.path.realpath(__file__))
        table1 = Table(current_dir + "/test_general/tsv.test.txt")
        assert table1.dim()[0] == 2
        assert table1.dim()[1] == 3
        assert np.array_equal(table1.colnames, ['A','B','C'])
        assert np.array_equal(table1.rownames, ['0','1'])
        assert table1[0] == ('a','1','1.')
        assert table1[1] == ('b','2','2.')
        assert table1[0] != ('a',1,1.)
        # with rownames
        table2 = Table(current_dir + "/test_general/tsv.test2.txt")
        assert table2.dim()[0] == 2
        assert table2.dim()[1] == 3
        assert np.array_equal(table2.colnames, ['A','B','C'])
        assert np.array_equal(table2.rownames, ['a','b'])
        assert table2[0] == ('a','1','1.')
        assert table2[1] == ('b','2','2.')
        # with rownames and col types
        table3 = Table(current_dir + "/test_general/tsv.test2.txt", col_types=[str,int,float])
        assert table3.dim()[0] == 2
        assert table3.dim()[1] == 3
        assert np.array_equal(table3.colnames, ['A','B','C'])
        assert np.array_equal(table3.rownames, ['a','b'])
        assert table3[0] == ('a',1,1.)
        assert table3[1] == ('b',2,2.)
        
        # from record array
        x = np.zeros((5,),dtype=('i4,f4,a10'))
        x[:] = [(1,2.,'Hello'),(2,3.,"World"),(3,4.,'this'),(4,5.,"is"),(5,6.,'me')]
        table4 = Table(x, rownames=['a','b','c','d','e'], colnames=['A','B','C'])
        assert table4.dim()[0] == 5
        assert table4.dim()[1] == 3
        assert np.array_equal(table4.colnames, ['A','B','C'])
        assert np.array_equal(table4.rownames, ['a','b','c','d','e'])
    
    def test_save_and_load(self, tmpdir):
        current_dir = os.path.dirname(os.path.realpath(__file__))
        h5_file_name = str(tmpdir) + "/test.h5"
        table2 = Table(current_dir + "/test_general/tsv.test.txt")
        
        table2.save_as(h5_file_name)
        
        assert table2.dim()[0] == 2
        assert table2.dim()[1] == 3
        assert np.array_equal(table2.colnames, ['A','B','C'])
        assert np.array_equal(table2.rownames, ['0','1'])
        assert table2[0] == ('a','1','1.')
        assert table2[1] == ('b','2','2.')
        
        assert os.path.exists(h5_file_name)
        
    def test_create_and_load(self, tmpdir):
        current_dir = os.path.dirname(os.path.realpath(__file__))
        h5_file_name = str(tmpdir) + "/test.h5"
        table1 = Table(current_dir + "/test_general/tsv.test.txt", file_name=h5_file_name)
        table1.close()
        
        table2 = Table(h5_file_name)
        assert table2.dim()[0] == 2
        assert table2.dim()[1] == 3
        assert np.array_equal(table2.colnames, ['A','B','C'])
        assert np.array_equal(table2.rownames, ['0','1'])
        assert table2[0] == ('a','1','1.')
        assert table2[1] == ('b','2','2.')
    
    def test_read_and_write_tsv(self, tmpdir):
        current_dir = os.path.dirname(os.path.realpath(__file__))
        table1 = Table(current_dir + "/test_general/tsv.test2.txt")
        
        # defaults
        default_file_name = str(tmpdir) + "/default_test.txt"
        table1.export(default_file_name)
        table2 = Table(default_file_name)
        assert table2.dim()[0] == 2
        assert table2.dim()[1] == 3
        assert np.array_equal(table2.colnames, ['A','B','C'])
        assert np.array_equal(table2.rownames, ['a','b'])
        assert table2[0] == ('a','1','1.')
        assert table2[1] == ('b','2','2.')
        
        table3 = Table(default_file_name)
        assert np.array_equal(table3.colnames, ['A','B','C'])
        assert np.array_equal(table3.rownames, ['a','b'])
        assert table3[0] == ('a','1','1.')
        assert table3[1] == ('b','2','2.')
        
        # no rownames
        no_rownames_file_name = str(tmpdir) + "/nr_test.txt"
        table1.export(no_rownames_file_name, include_rownames=False)
        table3 = Table(no_rownames_file_name)
        assert table3.dim()[0] == 2
        assert table3.dim()[1] == 3
        assert np.array_equal(table3.colnames, ['A','B','C'])
        assert np.array_equal(table3.rownames, ['0','1'])
        assert table3[0] == ('a','1','1.')
        assert table3[1] == ('b','2','2.')
        
    def test_append_row_list(self):
        # default rowname
        self.table._append_row_list([6,7.,'bla'])
        assert self.table[5] == (6,7.0,'bla')
        assert self.table.rownames[5] == '5'
        
        self.table._append_row_list(['f',6,7.,'bla'])
        assert self.table[6] == (6,7.0,'bla')
        assert self.table.rownames[6] == 'f'
        
        self.table._append_row_list([6,7.,'bla'], rowname='g')
        assert self.table[7] == (6,7.0,'bla')
        assert self.table.rownames[7] == 'g'
        
    def test_append_row_dict(self):
        self.table._append_row_dict({'A': 6, 'B': 7., 'C': 'bla'})
        assert self.table[5] == (6,7.0,'bla')
        assert self.table.rownames[5] == '5'
        
        self.table._append_row_dict({'A': 6, 'B': 7., 'C': 'bla', '_rowname': 'f'})
        assert self.table[6] == (6,7.0,'bla')
        assert self.table.rownames[6] == 'f'
        
    def test_append(self):
        # append a row by list:
        self.table.append([6,7.,'bla'], rownames='f')
        assert self.table[5] == (6,7.0,'bla')
        assert self.table.rownames[5] == 'f'
        
        # append a row by dict
        self.table._append_row_dict({'A': 6, 'B': 7., 'C': 'bla'})
        assert self.table[6] == (6,7.0,'bla')
        assert self.table.rownames[6] == '6'
        
        # append multiple rows by list
        self.table.append([[6,7.,'foo'],[7,8.,'bar']], rownames=['g','h'])
        assert self.table[7] == (6,7.0,'foo')
        assert self.table.rownames[7] == 'g'
        assert self.table[8] == (7,8.0,'bar')
        assert self.table.rownames[8] == 'h'
        
        # append multiple dicts by list
        self.table.append([{'A': 6, 'B': 7., 'C': 'foo'},{'A': 7, 'B': 8., 'C': 'bar'}])
        assert self.table[9] == (6,7.0,'foo')
        assert self.table.rownames[9] == '9'
        assert self.table[10] == (7,8.0,'bar')
        assert self.table.rownames[10] == '10'
        
    
    def test_row_selection(self):
        
        # single rows
        assert self.table._get_rows(0) == (1,2.,'Hello')
        assert self.table._get_rows('a') == (1,2.,'Hello')
        x = self.table._get_rows(0)
        assert np.array_equal(x.colnames, ['A','B','C'])
        assert x.rowname == 'a'
        
        # multiple distinct rows by list
        x = self.table._get_rows([1,3])
        assert tuple(x[0]) == (2,3.,"World")
        assert tuple(x[1]) == (4,5.,"is")
        assert np.array_equal(x.colnames, ['A','B','C'])
        assert np.array_equal(x.rownames, ['b','d'])
        
        x = self.table._get_rows(['b','d'])
        assert tuple(x[0]) == (2,3.,"World")
        assert tuple(x[1]) == (4,5.,"is")
        assert np.array_equal(x.colnames, ['A','B','C'])
        assert np.array_equal(x.rownames, ['b','d'])
        
        # in reverse order
        x = self.table._get_rows([3,1])
        assert tuple(x[0]) == (4,5.,"is")
        assert tuple(x[1]) == (2,3.,"World")
        assert np.array_equal(x.colnames, ['A','B','C'])
        assert np.array_equal(x.rownames, ['d','b'])
        
        x = self.table._get_rows(['d','b'])
        assert tuple(x[0]) == (4,5.,"is")
        assert tuple(x[1]) == (2,3.,"World")
        assert np.array_equal(x.colnames, ['A','B','C'])
        assert np.array_equal(x.rownames, ['d','b'])
        
        # slices
        x = self.table._get_rows(slice(0,2,1))
        assert tuple(x[0]) == (1,2.,'Hello')
        assert tuple(x[1]) == (2,3.,"World")
        
        x = self.table._get_rows(slice(0,1,1))
        assert x == (1,2.,'Hello')

        x = self.table._get_rows(slice(0,10,1))
        assert tuple(x[0]) == (1,2.,'Hello')
        assert tuple(x[4]) == (5,6.,"me")
        
        print x
        
    def test_col_selection(self):
        
        # single cols
        assert np.array_equal(self.table._get_cols(0),(1,2,3,4,5))
        assert np.array_equal(self.table._get_cols('A'),(1,2,3,4,5))
        
        assert np.array_equal(self.table._get_cols(1),(2.,3.,4.,5.,6.))
        assert np.array_equal(self.table._get_cols('B'),(2.,3.,4.,5.,6.))
        
        assert np.array_equal(self.table._get_cols(2),('Hello','World','this','is','me'))
        assert np.array_equal(self.table._get_cols('C'),('Hello','World','this','is','me'))
        
        x = self.table._get_cols(0)
        assert np.array_equal(x.colname, 'A')
        assert np.array_equal(x.rownames, ['a','b','c','d','e'])
        
        
        # multiple cols
        x = self.table._get_cols([1,2])
        assert x[0][0] == 2.
        assert x[4][1] == 'me'
        assert len(x[0]) == 2
        assert np.array_equal(x.colnames, ['B', 'C'])
        assert np.array_equal(x.rownames, ['a','b','c','d','e'])
        
        # (reverse)
        x = self.table._get_cols([2,1])
        assert x[0][1] == 2.
        assert x[4][0] == 'me'
        assert len(x[0]) == 2
        assert np.array_equal(x.colnames, ['C','B'])
        assert np.array_equal(x.rownames, ['a','b','c','d','e'])
        
        x = self.table._get_cols(['B','C'])
        assert x[0][0] == 2.
        assert x[4][1] == 'me'
        assert len(x[0]) == 2
        assert np.array_equal(x.colnames, ['B', 'C'])
        assert np.array_equal(x.rownames, ['a','b','c','d','e'])
        
        # slice
        x = self.table._get_cols(slice(1,3,1))
        assert x[0][0] == 2.
        assert x[4][1] == 'me'
        assert len(x[0]) == 2
        assert np.array_equal(x.colnames, ['B', 'C'])
        assert np.array_equal(x.rownames, ['a','b','c','d','e'])

        
        
    def test_selection(self):
        
        # single rows
        assert self.table[0] == (1,2.,'Hello')
        assert self.table['a'] == (1,2.,'Hello')
        x = self.table[0]
        assert np.array_equal(x.colnames, ['A','B','C'])
        assert x.rowname == 'a'
        
        # multiple distinct rows by list
        x = self.table[[1,3]]
        assert tuple(x[0]) == (2,3.,"World")
        assert tuple(x[1]) == (4,5.,"is")
        assert np.array_equal(x.colnames, ['A','B','C'])
        assert np.array_equal(x.rownames, ['b','d'])
        
        x = self.table[['b','d']]
        assert tuple(x[0]) == (2,3.,"World")
        assert tuple(x[1]) == (4,5.,"is")
        assert np.array_equal(x.colnames, ['A','B','C'])
        assert np.array_equal(x.rownames, ['b','d'])
        
        # in reverse order
        x = self.table[[3,1]]
        assert tuple(x[0]) == (4,5.,"is")
        assert tuple(x[1]) == (2,3.,"World")
        assert np.array_equal(x.colnames, ['A','B','C'])
        assert np.array_equal(x.rownames, ['d','b'])
        
        x = self.table[['d','b']]
        assert tuple(x[0]) == (4,5.,"is")
        assert tuple(x[1]) == (2,3.,"World")
        assert np.array_equal(x.colnames, ['A','B','C'])
        assert np.array_equal(x.rownames, ['d','b'])
        
        # slices
        x = self.table[0:2]
        assert tuple(x[0]) == (1,2.,'Hello')
        assert tuple(x[1]) == (2,3.,"World")
        
        x = self.table[0:1]
        assert x == (1,2.,'Hello')
        
        
        
        # single cols
        assert np.array_equal(self.table[:,0],(1,2,3,4,5))
        assert np.array_equal(self.table[:,'A'],(1,2,3,4,5))
        
        assert np.array_equal(self.table[:,1],(2.,3.,4.,5.,6.))
        assert np.array_equal(self.table[:,'B'],(2.,3.,4.,5.,6.))
        
        assert np.array_equal(self.table[:,2],('Hello','World','this','is','me'))
        assert np.array_equal(self.table[:,'C'],('Hello','World','this','is','me'))
        
        x = self.table[:,0]
        assert np.array_equal(x.colname, 'A')
        assert np.array_equal(x.rownames, ['a','b','c','d','e'])
        
        
        # multiple cols
        x = self.table[:,[1,2]]
        assert x[0][0] == 2.
        assert x[4][1] == 'me'
        assert len(x[0]) == 2
        assert np.array_equal(x.colnames, ['B', 'C'])
        assert np.array_equal(x.rownames, ['a','b','c','d','e'])
        
        # (reverse)
        x = self.table[:,[2,1]]
        assert x[0][1] == 2.
        assert x[4][0] == 'me'
        assert len(x[0]) == 2
        assert np.array_equal(x.colnames, ['C','B'])
        assert np.array_equal(x.rownames, ['a','b','c','d','e'])
        
        x = self.table[:,['B','C']]
        assert x[0][0] == 2.
        assert x[4][1] == 'me'
        assert len(x[0]) == 2
        assert np.array_equal(x.colnames, ['B', 'C'])
        assert np.array_equal(x.rownames, ['a','b','c','d','e'])
        
        # slice
        x = self.table[:,1:3]
        assert x[0][0] == 2.
        assert x[4][1] == 'me'
        assert len(x[0]) == 2
        assert np.array_equal(x.colnames, ['B', 'C'])
        assert np.array_equal(x.rownames, ['a','b','c','d','e'])
        
        
        
        
        # select single value
        assert self.table[0,0] == 1
        assert self.table[0,1] == 2.
        assert self.table[0,2] == 'Hello'
        
        # select row subset
        assert self.table[0,1:3] == (2., 'Hello')
        
        # select column subset
        assert np.array_equal(self.table[1:4,'A'], [2,3,4]) 
        
        # select array subset
        x = self.table[['a','c'],['B','C']]
        assert tuple(x[0]) == (2., 'Hello')
        assert tuple(x[1]) == (4., 'this')
        
        assert x[0,0] == 2.
        
        y = x[:,:]
        assert tuple(y[0]) == (2., 'Hello')
        assert tuple(y[1]) == (4., 'this')
        
    
    def test_where(self):
        # simple query
        x = self.table.where("A < 3")
        assert tuple(x[0]) == (1,2.,'Hello')
        assert tuple(x[1]) == (2,3.,"World")
        
        # complex query
        x = self.table.where("(A < 4) & (B <= 3.0)")
        assert tuple(x[0]) == (1,2.,'Hello')
        assert tuple(x[1]) == (2,3.,"World")
        assert len(x) == 2
        
        # complex query
        x = self.table.where('(A < 4) & (B <= 3.0) & (C == "Hello")')
        assert tuple(x) == (1,2.,'Hello')
        
    
    def test_attributes(self):
        assert np.array_equal(self.table.A,(1,2,3,4,5))
        assert np.array_equal(self.table.B,(2.,3.,4.,5.,6.))
        assert np.array_equal(self.table.C,('Hello','World','this','is','me'))


class TestMask:
    
    def test_instantiate(self):
        mask = Mask(0,'test','test description')
        assert mask.ix == 0
        assert mask.name == 'test'
        assert mask.description == 'test description'
        
        mask = Mask(0,'test')
        assert mask.ix == 0
        assert mask.name == 'test'
        assert mask.description == ''

class TestMaskable:
    
    def test_instantiate(self, tmpdir):
        # no args
        maskable1 = Maskable()
        assert isinstance(maskable1._mask, t.table.Table)
        
        # string args
        maskable2 = Maskable(str(tmpdir) + "/test1.h5")
        assert isinstance(maskable2._mask, t.table.Table)
        
        # file args
        h5_file = create_or_open_pytables_file(str(tmpdir) + "/test2.h5")
        maskable3 = Maskable(h5_file)
        assert isinstance(maskable3._mask, t.table.Table)
        
        # table args
        h5_file2 = create_or_open_pytables_file(str(tmpdir) + "/test3.h5")
        table = h5_file2.create_table("/", 'mask', Maskable.MaskDescription)
        maskable4 = Maskable(table)
        assert isinstance(maskable4._mask, t.table.Table)
        
        # inherited
        class MaskableContainerTest1(Maskable):
            def __init__(self):
                Maskable.__init__(self)
        
        mc1 = MaskableContainerTest1()
        assert isinstance(mc1._mask, t.table.Table)
        
        class MaskableContainerTest2(Maskable):
            def __init__(self, h5_file):
                self.file = h5_file
                Maskable.__init__(self)
        
        h5_file3 = create_or_open_pytables_file(str(tmpdir) + "/test4.h5")
        mc2 = MaskableContainerTest2(h5_file3)
        assert isinstance(mc2._mask, t.table.Table)
    
    def test_get_default_masks(self):
        maskable = Maskable()
        default = maskable.get_mask(0)
        assert default.ix == 0
        assert default.name == 'default'
        assert default.description == 'Default mask'
        
        default = maskable.get_mask('default')
        assert default.ix == 0
        assert default.name == 'default'
        assert default.description == 'Default mask'
    
    def test_add_description(self):
        maskable = Maskable()
        
        mask = maskable.add_mask_description('test', 'description')
        maskret = maskable.get_mask(1)
        assert maskret.ix == 1
        assert maskret.name == 'test'
        assert maskret.description == 'description'
        
        assert maskret.ix == mask.ix
        assert maskret.name == mask.name
        assert maskret.description == mask.description
        
    def test_get_masks_by_ix(self):
        
        maskable = Maskable()
        
        maskable.add_mask_description('one', '1')
        maskable.add_mask_description('two', '2')
        maskable.add_mask_description('three', '3')
        maskable.add_mask_description('four', '4')
        
        ix = 31
        masks = maskable.get_masks(ix)
        assert masks[0].name == 'default'
        assert masks[1].name == 'one'
        assert masks[2].name == 'two'
        assert masks[3].name == 'three'
        assert masks[4].name == 'four'
        
        ix = 29
        masks = maskable.get_masks(ix)
        assert masks[0].name == 'default'
        assert masks[1].name == 'two'
        assert masks[2].name == 'three'
        assert masks[3].name == 'four'
        
        ix = 16
        masks = maskable.get_masks(ix)
        assert masks[0].name == 'four'
        
        ix = -1
        masks = maskable.get_masks(ix)
        assert len(masks) == 0
        
class TestMaskedTable:
    class ExampleFilter(MaskFilter):
        def __init__(self, cutoff=25):
            super(TestMaskedTable.ExampleFilter, self).__init__()
            self.cutoff = cutoff
            
        def valid(self, test):
            if test['b'] < self.cutoff:
                return False
            return True
            
    @classmethod
    def setup_method(self, method):
        f = create_or_open_pytables_file()
        test_description = {
            'a': t.StringCol(50,pos=0),
            'b': t.Int32Col(pos=1),
            'c': t.Float32Col(pos=2)
        }
        self.table = MaskedTable(f.get_node("/"), "test", test_description)
        
        row = self.table.row
        for i in range(0,50):
            row['a'] = "test_%d" % i 
            row['b'] = 0 + i
            row['c'] = 0.0 + i
            row.append()
        
        self.table.flush(update_index=True)
        
        self.filtered_table = MaskedTable(f.get_node("/"), "test_filter", test_description)
        
        row = self.filtered_table.row
        for i in range(0,50):
            row['a'] = "test_%d" % i 
            row['b'] = 0 + i
            row['c'] = 0.0 + i
            row.append()
        
        self.filtered_table.flush(update_index=True)
        self.filtered_table.filter(TestMaskedTable.ExampleFilter())
        
    def test_initialize(self):
        assert len(self.table) == 50
        for i in range(0,50):
            assert self.table[i][4] == i
    
    def test_len(self):
        assert len(self.filtered_table) == 25
    
    def test_select(self):
        # single positive index        
        assert self.filtered_table[0][1] == 25
        # single negative index
        assert self.filtered_table[-1][1] == 49
        # slice
        x = self.filtered_table[1:3]
        assert np.array_equal(tuple(x[0]), ('test_26',26,26.0,0,1))
        assert np.array_equal(tuple(x[1]), ('test_27',27,27.0,0,2))
    
    def test_filter(self):
        self.table.filter(TestMaskedTable.ExampleFilter())
        
        i = 0
        masked_i = -1
        for row in self.table:
            if row[self.table._mask_field] > 0:
                assert row[self.table._mask_index_field] == masked_i
                assert row[self.table._mask_field] == 1
                masked_i -= 1
            else:
                assert row[self.table._mask_index_field] == i
                i += 1
            
    def test_masked(self):
        assert self.filtered_table.masked_rows()[0][1] == 0
        assert self.filtered_table.masked_rows()[-1][1] == 24
        
        masked_ix = -1
        for row in self.filtered_table.masked_rows():
            assert row[self.table._mask_index_field] == masked_ix
            masked_ix -= 1

    def test_exclude_filters(self):
        t = self.filtered_table
        assert len(list(t.iterrows(excluded_masks=1))) == 50

class TestMeta:
    def setup_method(self, method):
        self.meta = MetaContainer()

    def test_add_and_get(self):
        self.meta.add_meta("test", "test_name", 34, 0, "testing")
        meta_info = self.meta.get_meta(0)

    def test_history(self):
        self.meta.log_info("test1")
        self.meta.log_warn("test2")
        self.meta.log_error("test3")
        self.meta.log_debug("test4")

        history = self.meta.meta_history()


class RegisteredTable(t.Table):
    # Class identifier. Enough for registration,
    # which is taken care of by MetaNode 
    _c_classid = 'REGISTEREDTABLE'
    
    def __init__(self, parentnode, name,
                 description=None, title="", filters=None,
                 expectedrows=None, chunkshape=None,
                 byteorder=None, _log=True):

        
        # normal Table init
        t.Table.__init__(self, parentnode, name, description,
                         title, filters, expectedrows, chunkshape,
                         byteorder, _log)


class TestFileBased:
    def test_read_only(self, tmpdir):
        f = FileBased(str(tmpdir) + "/test.file")
        f.file.create_table("/", "test1", {'a': t.Int32Col()})
        f.close()

        r = FileBased(str(tmpdir) + "/test.file", mode='r')
        with pytest.raises(t.FileModeError):
            r.file.create_table("/", "test2", {'b': t.Int32Col()})
        r.close()

    def test_tmp(self, tmpdir):
        filename = str(tmpdir) + "/test.file"
        f = FileBased(file_name=filename, mode='a', tmpdir='/tmp')
        assert os.path.isfile(filename) == False
        assert os.path.isfile(f.tmp_file_name) == True
        f.close()
        f.finalize()
        assert os.path.isfile(filename) == True
        f.cleanup()
        assert os.path.isfile(f.tmp_file_name) == False

    def test_tmp_with(self, tmpdir):
        filename = str(tmpdir) + "/test.file"
        with FileBased(file_name=filename, mode='a', tmpdir='/tmp') as f:
            assert os.path.isfile(filename) == False
            assert os.path.isfile(f.tmp_file_name) == True
        assert os.path.isfile(filename) == True
        assert os.path.isfile(f.tmp_file_name) == False

    def test_tmp_with_exception(self, tmpdir):
        filename = str(tmpdir) + "/test.file"
        with pytest.raises(Exception):
            with FileBased(file_name=filename, mode='a', tmpdir='/tmp') as f:
                assert os.path.isfile(filename) == False
                assert os.path.isfile(f.tmp_file_name) == True
                try:
                    raise Exception
                except:
                    assert os.path.isfile(filename) == False
                    assert os.path.isfile(f.tmp_file_name) == False

    def test_tmp_with_existing(self, tmpdir):
        filename = str(tmpdir) + "/test.file"
        f = FileBased(str(tmpdir) + "/test.file")
        f.file.create_table("/", "test1", {'a': t.Int32Col()})
        f.close()
        assert os.path.isfile(filename) == True
        with FileBased(file_name=filename, mode='a', tmpdir='/tmp') as f:
            f.file.create_table("/", "test2", {'b': t.Int32Col()})
            assert os.path.isfile(f.tmp_file_name) == True
        assert os.path.isfile(filename) == True
        assert os.path.isfile(f.tmp_file_name) == False
        with FileBased(file_name=filename, mode='r') as f:
            assert 'test1' in f.file.root
            assert 'test2' in f.file.root

class TestPytablesInheritance:
    class MinTable(t.Table):

        def __init__(self):
    
            f = t.open_file('bla', 'a', driver="H5FD_CORE",driver_core_backing_store=0)
    
            description = {
                'c1': t.Int32Col(),
                'c2': t.Int16Col()
            }
    
            t.Table.__init__(self, f.get_node('/'), 'test',
                             description=description)
    
            self._enable_index()
    
            row = self.row
            row['c1'] = 0
            row.append()
            self.flush()
    
        def _enable_index(self):
            self.cols.c1.create_index()
            
    
    
    def test_table_wrapper_with_index(self):
        m = TestPytablesInheritance.MinTable()
        # next line fails in unpatched code
        [x.fetch_all_fields() for x in m.where('c1 == 0')]

    def test_no_wrapper_with_index(self):
        f = t.open_file('bla2', 'a', driver="H5FD_CORE",driver_core_backing_store=0)
        table = t.Table(f.get_node('/'),'test',{ 'c1': t.Int32Col(), 'c2': t.Int16Col() },title='test')
        table.row['c1'] = 0
        table.row.append()
        table.flush()
        table.cols.c1.create_index()
        [x.fetch_all_fields() for x in table.where('c1 == 0')]
        
    def test_registered_table(self, tmpdir):
        h5_file = create_or_open_pytables_file(str(tmpdir) + "/test.h5")
        table = RegisteredTable(h5_file.get_node('/'), 'test', { 'c1': t.Int32Col(), 'c2': t.Int16Col() }, 'test')
        print class_id_dict
        table.row['c1'] = 0
        table.row.append()
        table.flush()
        table.close()
        h5_file.close()

        h5_file_ret = create_or_open_pytables_file(str(tmpdir) + "/test.h5")
        table = h5_file_ret.get_node('/','test')
        assert type(table) == RegisteredTable
        #assert 0
