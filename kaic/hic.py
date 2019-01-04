from .matrix import RegionMatrixTable


class Hic(RegionMatrixTable):
    def __init__(self, file_name=None, mode='a', tmpdir=None,
                 partitioning_strategy='chromosome',
                 _table_name_regions='regions', _table_name_edges='edges',
                 _edge_buffer_size=1000000):
        RegionMatrixTable.__init__(self, file_name=file_name,
                                   mode=mode, tmpdir=tmpdir,
                                   partitioning_strategy=partitioning_strategy,
                                   _table_name_regions=_table_name_regions,
                                   _table_name_edges=_table_name_edges,
                                   _edge_buffer_size=_edge_buffer_size)
