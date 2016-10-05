import tables as t
import numpy as np
import pytest
from kaic.data.general import Mask, Maskable, MaskedTable, MaskFilter, FileBased
from __builtin__ import classmethod
import os
from kaic.tools.files import create_or_open_pytables_file


class TestMask:
    
    def test_instantiate(self):
        mask = Mask(0, 'test', 'test description')
        assert mask.ix == 0
        assert mask.name == 'test'
        assert mask.description == 'test description'
        
        mask = Mask(0, 'test')
        assert mask.ix == 0
        assert mask.name == 'test'
        assert mask.description == ''


class TestMaskable:
    
    def test_instantiate(self, tmpdir):
        # no args
        maskable1 = Maskable()
        assert isinstance(maskable1._mask, t.table.Table)
        maskable1.close()
        
        # string args
        maskable2 = Maskable(str(tmpdir) + "/test1.h5")
        assert isinstance(maskable2._mask, t.table.Table)
        maskable2.close()
        
        # file args
        h5_file = create_or_open_pytables_file(str(tmpdir) + "/test2.h5")
        maskable3 = Maskable(h5_file)
        assert isinstance(maskable3._mask, t.table.Table)
        maskable3.close()
        h5_file.close()
        
        # table args
        h5_file2 = create_or_open_pytables_file(str(tmpdir) + "/test3.h5")
        table = h5_file2.create_table("/", 'mask', Maskable.MaskDescription)
        maskable4 = Maskable(table)
        assert isinstance(maskable4._mask, t.table.Table)
        maskable4.close()
        h5_file2.close()
        
        # inherited
        class MaskableContainerTest1(Maskable):
            def __init__(self):
                Maskable.__init__(self)
        
        mc1 = MaskableContainerTest1()
        assert isinstance(mc1._mask, t.table.Table)
        mc1.close()
        
        class MaskableContainerTest2(Maskable):
            def __init__(self, h5_file):
                self.file = h5_file
                Maskable.__init__(self)
        
        h5_file3 = create_or_open_pytables_file(str(tmpdir) + "/test4.h5")
        mc2 = MaskableContainerTest2(h5_file3)
        assert isinstance(mc2._mask, t.table.Table)
        mc2.close()
    
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
        maskable.close()
    
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
        maskable.close()
        
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
        maskable.close()
        
class TestMaskedTable:
    class ExampleFilter(MaskFilter):
        def __init__(self, cutoff=25):
            super(TestMaskedTable.ExampleFilter, self).__init__()
            self.cutoff = cutoff
            
        def valid(self, test):
            if test['b'] < self.cutoff:
                return False
            return True
            
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
        self.file = f

    def teardown_method(self, method):
        self.file.close()
        
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
            self.file = f
    
        def _enable_index(self):
            self.cols.c1.create_index()
            
    
    
    def test_table_wrapper_with_index(self):
        m = TestPytablesInheritance.MinTable()
        # next line fails in unpatched code
        [x.fetch_all_fields() for x in m.where('c1 == 0')]
        m.file.close()

    def test_no_wrapper_with_index(self):
        f = t.open_file('bla2', 'a', driver="H5FD_CORE",driver_core_backing_store=0)
        table = t.Table(f.get_node('/'),'test',{ 'c1': t.Int32Col(), 'c2': t.Int16Col() },title='test')
        table.row['c1'] = 0
        table.row.append()
        table.flush()
        table.cols.c1.create_index()
        [x.fetch_all_fields() for x in table.where('c1 == 0')]
        f.close()
        
    def test_registered_table(self, tmpdir):
        h5_file = create_or_open_pytables_file(str(tmpdir) + "/test.h5")
        table = RegisteredTable(h5_file.get_node('/'), 'test', { 'c1': t.Int32Col(), 'c2': t.Int16Col() }, 'test')
        table.row['c1'] = 0
        table.row.append()
        table.flush()
        table.close()
        h5_file.close()

        h5_file_ret = create_or_open_pytables_file(str(tmpdir) + "/test.h5")
        table = h5_file_ret.get_node('/','test')
        assert type(table) == RegisteredTable
        h5_file_ret.close()


class TestMetaInformation:

    def test_create(self):
        with FileBased() as f:
            f.meta.test = 'test'
            assert f.meta.test == 'test'

            f.meta['test2'] = 1
            assert f.meta['test2'] == 1

            assert f.meta._classid == 'FILEBASED'

            with pytest.raises(AttributeError):
                _ = f.meta.foo

            with pytest.raises(KeyError):
                _ = f.meta['foo']

    def test_load(self, tmpdir):
        with FileBased(file_name=str(tmpdir) + '/test.h5') as f:
            f.meta.test = 'test'
            f.meta.test2 = 1

        with FileBased(file_name=str(tmpdir) + '/test.h5', mode='r') as f:
            assert f.meta.test == 'test'
            assert f.meta.test2 == 1

            with pytest.raises(t.FileModeError):
                f.meta.test = 'foo'

            with pytest.raises(t.FileModeError):
                f.meta['test'] = 'foo'
