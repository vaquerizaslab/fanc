from cooler import Cooler
from ..matrix import RegionMatrixContainer, Edge
from genomic_regions import GenomicRegion


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


class CoolerRegions(RegionMatrixContainer, Cooler):
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

    def _edges_subset(self, key=None, *args, **kwargs):
        row_regions, col_regions = self._key_to_regions(key)

        row_start, row_end = self._min_max_region_ix(row_regions)
        col_start, col_end = self._min_max_region_ix(col_regions)

        df = Cooler.matrix(self, as_pixels=True, balance=False)[row_start:row_end+1, col_start:col_end+1]

        for index, row in df.iterrows():
            yield self._series_to_edge(row, *args, **kwargs)

    def _edges_getitem(self, item, *args, **kwargs):
        edges = []
        df = self.pixels()[item]
        for index, row in df.iterrows():
            print(row)
            edges.append(self._series_to_edge(row))

        if isinstance(item, int):
            return edges[0]
        return edges

    def _edges_length(self):
        return len(self.pixels())
