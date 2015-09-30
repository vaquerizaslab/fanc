'''
Created on Sep 29, 2015

@author: kkruse1
'''

import tables as t
from kaic.data.general import FileBased
import random
import string

class DummyObject(object):
    def __init__(self, x=None, y=None, z=None):
        self.x = x
        self.y = y
        self.z = z
    
class DummyObjectFromRow(object):
    def __init__(self, row):
        self.row = row
    
    @property
    def x(self, x=None):
        if x is not None:
            self.row['x'] = x
            self.row.update()
        return self.row['x']
    
    @property
    def y(self, y=None):
        if y is not None:
            self.row['y'] = y
            self.row.update()
        return self.row['y']
    
    @property
    def z(self, z=None):
        if z is not None:
            self.row['z'] = z
            self.row.update()
        return self.row['z']
    
    
        

class BasicTable(FileBased):
    class BasicDummyDescription(t.IsDescription):
        x = t.Int64Col(pos=0)
        y = t.Float32Col(pos=1)
        z = t.StringCol(100, pos=2)
        
    def __init__(self):
        FileBased.__init__(self, None)
        
        self._table = t.Table(self.file.root, 'test',
                              BasicTable.BasicDummyDescription,
                              expectedrows=10000)
        
        row = self._table.row
        for i in xrange(0,10000):
            x = random.randint(0, 9000000000000000000)
            y = random.uniform(0, 1000000)
            z = ''.join(random.choice(string.lowercase) for _ in range(100))
            
            row['x'] = x
            row['y'] = y
            row['z'] = z
            row.append()
        
        self._table.flush()
    

class StandardIteratorObject(BasicTable):
    def __init__(self):
        super(StandardIteratorObject, self).__init__()
        
    def __iter__(self):
        this = self
        class StandardIter:
            def __init__(self):
                self.iter = iter(this._table)
                
            def __iter__(self):
                return self
            
            def next(self):
                return self.iter.next()
            
        return StandardIter()


class ConversionIteratorObject(BasicTable):
    def __init__(self):`
        super(ConversionIteratorObject, self).__init__()
        
    def __iter__(self):
        this = self
        class ConversionIter:
            def __init__(self):
                self.iter = iter(this._table)
                
            def __iter__(self):
                return self
            
            def next(self):
                row = self.iter.next()
                return DummyObject(x=row['x'], y=row['y'], z=row['z'])
            
        return ConversionIter()
    
class FromRowIteratorObject(BasicTable):
    def __init__(self):
        super(FromRowIteratorObject, self).__init__()
        
    def __iter__(self):
        this = self
        class FromRowIter:
            def __init__(self):
                self.iter = iter(this._table)
                
            def __iter__(self):
                return self
            
            def next(self):
                row = self.iter.next()
                return DummyObjectFromRow(row)
            
        return FromRowIter()



def test_iterators(test_subject):
    res = ["%d %.5f %s" % (row.x, row.y, row.z) for row in test_subject]
    return res

def test_iterators_row(test_subject):
    res = ["%d %.5f %s" % (row['x'], row['y'], row['z']) for row in test_subject]
    return res


