import logging

import os
import h5py
import numpy as np
import pandas
import cooler
from genomic_regions import GenomicRegion
from ..tools.sambam import natural_cmp
from functools import cmp_to_key
import tempfile


from ..matrix import RegionMatrixContainer, Edge

logger = logging.getLogger(__name__)


def is_cooler(file_name):
    try:
        if "@" in file_name:
            hic_file, at_resolution = file_name.split("@")
            file_name = hic_file + '::resolutions/{}'.format(at_resolution)
        cooler.Cooler(file_name)
        return True
    except KeyError:
        return False


def to_cooler(hic, path, balance=True, multires=True,
              resolutions=None, n_zooms=10, threads=1,
              chunksize=100000, max_resolution=5000000,
              natural_order=True, chromosomes=None,
              **kwargs):
    """
    Export Hi-C data as Cooler file.

    Only contacts that have not been
    filtered are exported. https://github.com/mirnylab/cooler/

    Single resolution files:
    If input Hi-C matrix is uncorrected, the uncorrected matrix is stored.
    If it is corrected, the uncorrected matrix is stored along with bias vector.
    Cooler always calculates corrected matrix on-the-fly from the uncorrected
    matrix and the bias vector.

    Multi-resolution files (default):


    :param hic: Hi-C file in any compatible (RegionMatrixContainer) format
    :param path: Output path for cooler file
    :param balance: Include bias vector in cooler output (single res) or perform
                    iterative correction (multi res)
    :param multires: Generate a multi-resolution cooler file
    :param resolutions: Resolutions in bp (int) for multi-resolution cooler output
    :param chunksize: Number of pixels processed at a time in cooler
    :param kwargs: Additional arguments passed to cooler.iterative_correction
    """
    base_resolution = hic.bin_size

    tmp_files = []
    try:
        if multires:
            if resolutions is None:
                resolutions = [base_resolution * 2 ** i for i in range(n_zooms)
                               if base_resolution * 2 ** i < max_resolution]
            else:
                for r in resolutions:
                    if r % base_resolution != 0:
                        raise ValueError("Resolution {} must be a multiple of "
                                         "base resolution {}!".format(r, base_resolution))

            single_path = tempfile.NamedTemporaryFile(delete=False, suffix='.cool').name
            tmp_files.append(single_path)
            multi_path = path
        else:
            single_path = path
            multi_path = None

        natural_key = cmp_to_key(natural_cmp)
        if chromosomes is None:
            chromosomes = hic.chromosomes()
            if natural_order:
                chromosomes = sorted(chromosomes, key=lambda x: natural_key(x.encode('utf-8')))

        logger.info("Loading genomic regions")
        ix_converter = dict()
        regions = []
        region_order = []
        new_region_index = 0
        for chromosome in chromosomes:
            for region in hic.regions(chromosome, lazy=True):
                regions.append((region.chromosome,
                                region.start - 1,
                                region.end))
                ix_converter[region.ix] = new_region_index
                region_order.append(region.ix)
                new_region_index += 1
        region_df = pandas.DataFrame(regions, columns=['chrom', 'start', 'end'])

        def pixel_iter():
            for chri in range(len(chromosomes)):
                chromosome1 = chromosomes[chri]
                for chrj in range(chri, len(chromosomes)):
                    chromosome2 = chromosomes[chrj]

                    logger.info("{} - {}".format(chromosome1, chromosome2))

                    def chromosome_pixel_iter():
                        for edge in hic.edges((chromosome1, chromosome2), norm=False, lazy=True):
                            source, sink = ix_converter[edge.source], ix_converter[edge.sink]
                            if sink < source:
                                source, sink = sink, source
                            yield source, sink, edge.weight

                    pixels = np.fromiter(chromosome_pixel_iter(),
                                         dtype=[("bin1_id", np.int_),
                                                ("bin2_id", np.int_),
                                                ("count", np.float_)])
                    pixels = np.sort(pixels, order=("bin1_id", "bin2_id"))
                    yield pandas.DataFrame(pixels)

        logger.info("Writing cooler")
        cooler.create_cooler(cool_uri=single_path, bins=region_df, pixels=pixel_iter(), ordered=False)

        cool_path, group_path = cooler.util.parse_cooler_uri(single_path)

        if not multires:
            if balance:
                logger.info("Writing bias vector from FAN-C matrix")
                bias = hic.bias_vector()[np.array(region_order)]

                # Copied this section from
                # https://github.com/mirnylab/cooler/blob/356a89f6a62e2565f42ff13ec103352f20d251be/cooler/cli/balance.py#L195
                with h5py.File(cool_path, 'r+') as h5:
                    grp = h5[group_path]
                    # add the bias column to the file
                    h5opts = dict(compression='gzip', compression_opts=6)
                    grp['bins'].create_dataset("weight", data=bias, **h5opts)
            return CoolerHic(single_path)
        else:
            cooler.zoomify_cooler(single_path, multi_path, resolutions, chunksize, nproc=threads)
            if balance:
                logger.info("Balancing zoom resolutions...")
                for resolution in resolutions:
                    uri = multi_path + "::resolutions/" + str(resolution)
                    cool_path, group_path = cooler.util.parse_cooler_uri(uri)
                    cool = cooler.Cooler(uri)
                    bias, stats = cooler.balance_cooler(cool, chunksize=chunksize, **kwargs)
                    with h5py.File(cool_path, 'r+') as h5:
                        grp = h5[group_path]
                        # add the bias column to the file
                        h5opts = dict(compression='gzip', compression_opts=6)
                        grp['bins'].create_dataset("weight", data=bias, **h5opts)
                        grp['bins']['weight'].attrs.update(stats)
            return CoolerHic(multi_path + '::resolutions/{}'.format(base_resolution))
    finally:
        for tmp_file in tmp_files:
            os.remove(tmp_file)


class LazyCoolerEdge(Edge):
    def __init__(self, series, c):
        self._series = series
        self._c = c
        self._weight_field = 'weight'
        self._bias = 1.

    def __getattr__(self, item):
        if item != self._weight_field:
            return object.__getattribute__(self, item)
        return object.__getattribute__(self, item) * self._bias

    def __getitem__(self, item):
        try:
            return getattr(self, item)
        except AttributeError:
            raise KeyError("No such key: {}".format(item))

    @property
    def bias(self):
        return self._bias

    @bias.setter
    def bias(self, b):
        self._bias = b

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


class LazyCoolerRegion(GenomicRegion):
    def __init__(self, series, ix=None):
        self._series = series
        self.ix = ix

    def __getattr__(self, item):
        return getattr(self._series, item)

    @property
    def chromosome(self):
        try:
            return self._series.chrom
        except Exception as e:
            print(e)
            raise

    @property
    def start(self):
        return self._series.start + 1

    @property
    def bias(self):
        return self._series.weight

    @property
    def strand(self):
        try:
            return self._series.strand
        except AttributeError:
            return 1


class CoolerHic(RegionMatrixContainer, cooler.Cooler):
    def __init__(self, *args, **kwargs):
        largs = list(args)
        if "@" in args[0]:
            hic_file, at_resolution = args[0].split("@")
            largs[0] = hic_file + '::resolutions/{}'.format(at_resolution)
        cooler.Cooler.__init__(self, *largs, **kwargs)
        RegionMatrixContainer.__init__(self)
        self._mappability = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return True

    def _series_to_region(self, series, ix=None, lazy_region=None):
        if lazy_region is not None:
            lazy_region._series = series
            lazy_region.ix = ix
            return lazy_region

        index = set(series.index)
        index.remove('chrom')
        try:
            index.remove('ix')
        except KeyError:
            pass

        kwargs = {name: series[name] for name in index}
        kwargs['chromosome'] = series.chrom
        try:
            kwargs['bias'] = series.weight
        except AttributeError:
            kwargs['bias'] = 1.
        kwargs['start'] = series.start + 1
        if ix is not None:
            kwargs['ix'] = ix
        return GenomicRegion(**kwargs)

    def _region_iter(self, lazy=False, *args, **kwargs):
        lazy_region = LazyCoolerRegion(None) if lazy else None

        for df in self.bins():
            yield self._series_to_region(df.iloc[0], ix=df.index[0], lazy_region=lazy_region)

    def _region_subset(self, region, lazy=False, *args, **kwargs):
        lazy_region = LazyCoolerRegion(None) if lazy else None

        query = "{}".format(region.chromosome)
        if region.start is not None and region.end is not None:
            query += ':{}-{}'.format(region.start - 1, region.end)

        df = self.bins().fetch(query)
        for index, row in df.iterrows():
            yield self._series_to_region(row, ix=index, lazy_region=lazy_region)

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

    def _chromosome_bins(self, *args, **kwargs):
        return RegionMatrixContainer._chromosome_bins(self, lazy=True)

    def _series_to_edge(self, series, lazy_edge=None):
        if lazy_edge is None:
            index = set(series.index)
            index.remove('bin1_id')
            index.remove('bin2_id')
            index.remove('count')

            kwargs = {name: series[name] for name in index}
            kwargs['source'] = int(series['bin1_id'])
            kwargs['sink'] = int(series['bin2_id'])
            kwargs['weight'] = float(series['count'])
            kwargs['source_node'] = self.regions[kwargs['source']]
            kwargs['sink_node'] = self.regions[kwargs['sink']]

            return Edge(**kwargs)
        else:
            lazy_edge._series = series
            return lazy_edge

    def _edges_iter(self, lazy=False, *args, **kwargs):
        lazy_edge = LazyCoolerEdge(None, self) if lazy else None
        selector = cooler.Cooler.matrix(self, as_pixels=True, balance=False)
        for df in selector:
            for index, row in df.iterrows():
                yield self._series_to_edge(row, lazy_edge=lazy_edge)

    def _edges_subset(self, key=None, row_regions=None, col_regions=None,
                      lazy=False, *args, **kwargs):
        lazy_edge = LazyCoolerEdge(None, self) if lazy else None
        row_start, row_end = self._min_max_region_ix(row_regions)
        col_start, col_end = self._min_max_region_ix(col_regions)

        df = cooler.Cooler.matrix(self, as_pixels=True, balance=False)[row_start:row_end+1, col_start:col_end+1]

        for index, row in df.iterrows():
            yield self._series_to_edge(row, lazy_edge=lazy_edge, *args, **kwargs)

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

    def mappable(self, region=None):
        """
        Get the mappability vector of this matrix.
        """
        if self._mappability is not None:
            return self._mappability

        mappable = [False] * len(self.regions(region))
        for edge in self.edges(region, lazy=True):
            mappable[edge.source] = True
            mappable[edge.sink] = True
        self._mappability = mappable
        return mappable
