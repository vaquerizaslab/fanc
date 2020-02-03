from fanc.architecture.architecture import TableArchitecturalFeature, calculateondemand
import tables as t
import types


class TAF(TableArchitecturalFeature):
    """
    This only exists so we can instantiate a VectorArchitecturalRegionFeature
    for testing.
    """
    def __init__(self, group, fields, file_name=None, mode='a', tmpdir=None):
        TableArchitecturalFeature.__init__(self, group, fields,
                                           file_name=file_name, mode=mode, tmpdir=tmpdir)

    def _calculate(self, *args, **kwargs):
        self.data('a', [1, 2, 3, 4, 5])
        self.data('b', ['a', 'b', 'c', 'd', 'e'])

    @calculateondemand
    def a(self, key=None):
        if key is None:
            return self[:, 'a']
        return self[key, 'a']


class TestTableArchitecturalFeature:
    def setup_method(self, method):
        self.taf = TAF('test', {'a': t.Int32Col(), 'b': t.StringCol(10)})

    def teardown_method(self, method):
        self.taf.close()

    def test_get_rows(self):
        assert isinstance(self.taf[0], dict)
        assert self.taf[0]['a'] == 1
        assert self.taf[0]['b'] == 'a'

        results = self.taf[1:3]
        assert isinstance(results, types.GeneratorType)
        for i, r in enumerate(results):
            assert r['a'] == i+2
            assert r['b'] == 'abcd'[i+1]

    def test_get_columns(self):
        # let's test single ones
        results = self.taf[0, 'a']
        assert results == 1
        results = self.taf[1, 'b']
        assert results == 'b'

        # int
        results = self.taf[0, 0]
        assert results == 1
        results = self.taf[1, 1]
        assert results == 'b'

        # slice
        results = self.taf[0, 0:2]
        assert results == {'a': 1, 'b': 'a'}
        results = self.taf[1, 0:2]
        assert results == {'a': 2, 'b': 'b'}

    def test_calculateondemand_decorator(self):
        assert self.taf.a() == [1, 2, 3, 4, 5]
        assert self.taf.a(1) == 2

    def test_set_item(self):
        assert self.taf[0, 'a'] == 1
        self.taf[0, 'a'] = 9
        assert self.taf[0, 'a'] == 9

        assert self.taf[:, 'a'] == [9, 2, 3, 4, 5]
        self.taf[:, 'a'] = [5, 6, 7, 8, 9]
        assert self.taf[:, 'a'] == [5, 6, 7, 8, 9]

        assert self.taf[1:4, 'b'] == ['b', 'c', 'd']
        self.taf[1:4, 'b'] = ['f', 'ghi', 'h']
        assert self.taf[1:4, 'b'] == ['f', 'ghi', 'h']
