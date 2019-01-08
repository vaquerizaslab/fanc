import logging

import h5py
import numpy as np
import pandas
from cooler import Cooler
from cooler.io import parse_cooler_uri, create
from genomic_regions import GenomicRegion

from ..matrix import RegionMatrixContainer, Edge

logger = logging.getLogger(__name__)


def to_cooler(hic, path):
    """
    Export Hi-C data as cooler file. Only contacts that have not been
    filtered are exported.
    https://github.com/mirnylab/cooler/
    If input Hi-C matrix is uncorrected, the uncorrected matrix is stored.
    If it is corrected, the uncorrected matrix is stored and the bias vector.
    Cooler always calculates corrected matrix on-the-fly from the uncorrected
    matrix and the bias vector.
    :param hic: Hi-C file in any compatible (RegionMatrixContainer) format
    :param path: Output path for cooler file
    """

    contact_dtype = [("source", np.int_), ("sink", np.int_), ("weight", np.float_)]
    bias = hic.region_data('bias')

    logger.info("Loading contacts")
    contact_array = np.fromiter((edge.source, edge.sink, edge.weight
                                 for edge in hic.edges(lazy=True, norm=False)),
                                dtype=contact_dtype, count=len(hic.edges))
    logger.info("Sorting contacts")
    order = np.argsort(contact_array, order=("source", "sink"))
    counts = np.rint(contact_array["weight"]).astype(np.int_)
    contact_dict = {
        "bin1_id": contact_array["source"][order],
        "bin2_id": contact_array["sink"][order],
        "count": counts[order],
    }
    region_dicts = [{"chrom": r.chromosome, "start": r.start - 1, "end": r.end} for r in hic.regions()]
    region_df = pandas.DataFrame(region_dicts)
    logger.info("Writing cooler")

    create(cool_uri=path, bins=region_df, pixels=contact_dict)

    cool_path, group_path = parse_cooler_uri(path)
    logger.info("Writing bias vector")
    # Copied this section from
    # https://github.com/mirnylab/cooler/blob/356a89f6a62e2565f42ff13ec103352f20d251be/cooler/cli/balance.py#L195
    with h5py.File(cool_path, 'r+') as h5:
        grp = h5[group_path]
        # add the bias column to the file
        h5opts = dict(compression='gzip', compression_opts=6)
        grp['bins'].create_dataset("weight", data=bias, **h5opts)


class LazyCoolerEdge(Edge):
    def __init__(self, series, c):
        self._series = series
        self._c = c

    @property
    def source(self):
        return int(self._series.bin1_id)

    @property
    def sink(self):
        return int(self._series.bin2_id)

    @property
    def weight(self):
        return float(self._series['count'])

    @property
    def source_node(self):
        return self._c.regions[self.source]

    @property
    def sink_node(self):
        return self._c.regions[self.sink]

    @property
    def field_names(self):
        return ['weight']


class CoolerMatrix(RegionMatrixContainer, Cooler):
    def __init__(self, *args, **kwargs):
        Cooler.__init__(self, *args, **kwargs)
        RegionMatrixContainer.__init__(self)

    def _series_to_region(self, series, ix=None):
        index = set(series.index)
        index.remove('chrom')
        try:
            index.remove('ix')
        except KeyError:
            pass

        kwargs = {name: series[name] for name in index}
        kwargs['chromosome'] = series.chrom
        if ix is not None:
            kwargs['ix'] = ix
        return GenomicRegion(**kwargs)

    def _region_iter(self, *args, **kwargs):
        for df in self.bins():
            yield self._series_to_region(df.iloc[0], ix=df.index[0])

    def _region_subset(self, region, *args, **kwargs):
        query = "{}".format(region.chromosome)
        if region.start is not None and region.end is not None:
            query += ':{}-{}'.format(region.start, region.end)

        df = self.bins().fetch(query)
        for index, row in df.iterrows():
            yield self._series_to_region(row, ix=index)

    def _get_regions(self, item, *args, **kwargs):
        regions = []
        df = self.bins()[item]
        for index, row in df.iterrows():
            regions.append(self._series_to_region(row, ix=index))

        if isinstance(item, int):
            return regions[0]
        return regions

    def _region_len(self):
        return len(self.bins())

    def chromosomes(self):
        cs = []
        for index, row in self.chroms()[:].iterrows():
            cs.append(row['name'])
        return cs

    @property
    def chromosome_lengths(self):
        cl = {}
        for index, row in self.chroms()[:].iterrows():
            cl[row['name']] = row['length']
        return cl

    def _series_to_edge(self, series, lazy=True):
        if not lazy:
            index = set(series.index)
            index.remove('bin1_id')
            index.remove('bin2_id')
            index.remove('count')

            kwargs = {name: series[name] for name in index}
            kwargs['source'] = int(series.bin1_id)
            kwargs['sink'] = int(series.bin2_id)
            kwargs['weight'] = float(series.count)
            kwargs['source_node'] = self.regions[series.bin1_id]
            kwargs['sink_node'] = self.regions[series.bin2_id]

            return Edge(**kwargs)
        else:
            return LazyCoolerEdge(series, self)

    def _edges_iter(self, *args, **kwargs):
        for df in self.bins():
            yield self._series_to_edge(df.iloc[0])

    def _edges_subset(self, key=None, row_regions=None, col_regions=None, *args, **kwargs):
        row_start, row_end = self._min_max_region_ix(row_regions)
        col_start, col_end = self._min_max_region_ix(col_regions)

        df = Cooler.matrix(self, as_pixels=True, balance=False)[row_start:row_end+1, col_start:col_end+1]

        for index, row in df.iterrows():
            yield self._series_to_edge(row, *args, **kwargs)

    def _edges_getitem(self, item, *args, **kwargs):
        edges = []
        df = self.pixels()[item]
        for index, row in df.iterrows():
            edges.append(self._series_to_edge(row))

        if isinstance(item, int):
            return edges[0]
        return edges

    def _edges_length(self):
        return len(self.pixels())
