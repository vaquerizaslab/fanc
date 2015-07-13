
import numpy as np
from kaic.data.general import Table, _to_list_and_names, TableRow, TableCol,\
    TableArray
from __builtin__ import classmethod
import os

class TestSupport:
    
    def test_list_conversion(self):
        # structured array
        x = np.zeros((5,),dtype=('i4,f4,a10'))
        x[:] = [(1,2.,'Hello'),(2,3.,"World"),(3,4.,'this'),(4,5.,"is"),(5,6.,'me')]
        l, n, t = _to_list_and_names(x)
        assert l[0] == [1,2.,'Hello']
        assert l[4] == [5,6.,'me']
        assert n == ['f0','f1','f2']
        #assert t == [int, float, str]
        
        
        # 2D array
        x = np.zeros((5,3))
        l, n, t = _to_list_and_names(x)
        assert l[0] == [0.,0.,0.]
        assert l[4] == [0.,0.,0.]
        assert n == ['0','1','2']
        
        # 2D list
        x = [(1,2.,'Hello'),(2,3.,"World"),(3,4.,'this'),(4,5.,"is"),(5,6.,'me')]
        l, n, t = _to_list_and_names(x)
        assert l[0] == [1,2.,'Hello']
        assert l[4] == [5,6.,'me']
        assert n == ['0','1','2']
        
        # 1D array
        x = np.array([1,2,3])
        l, n, t = _to_list_and_names(x)
        assert l[0] == [1,2,3]
        assert n == ['0','1','2']
        
        # 1D list
        x = [1,2,3]
        l, n, t = _to_list_and_names(x)
        assert l[0] == [1,2,3]
        assert n == ['0','1','2']
        
        # scalar
        x = 1
        l, n, t = _to_list_and_names(x)
        assert l[0] == [1]
        assert n == ['0']


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
        
        
        self.table = Table(x, row_names=['a','b','c','d','e'], col_names=['A','B','C'])

    
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
        table4 = Table(x, row_names=['a','b','c','d','e'], col_names=['A','B','C'])
        assert table4.dim()[0] == 5
        assert table4.dim()[1] == 3
        assert np.array_equal(table4.colnames, ['A','B','C'])
        assert np.array_equal(table4.rownames, ['a','b','c','d','e'])
    
    def test_save_and_load(self, tmpdir):
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
        