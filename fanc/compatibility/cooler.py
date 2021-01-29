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
import warnings
from ..tools.files import tmp_file_name
from ..tools.general import str_to_int
import shutil


from ..matrix import RegionMatrixContainer, Edge

logger = logging.getLogger(__name__)


def is_cooler(file_name):
    try:
        if "@" in file_name:
            fields = file_name.split("@")
            if len(fields) != 2:
                return False
            hic_file, at_resolution = fields
            file_name = hic_file + '::resolutions/{}'.format(str_to_int(at_resolution))
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
                    if len(pixels) > 0:
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


class LazyCoolerEdge(object):
    def __init__(self, series, c):
        self._series = series
        self._c = c
        self._weight_field = 'weight'
        self.bias = 1.
        self.expected = None

    def __getattr__(self, item):
        return getattr(self._series, item)

    @property
    def weight(self):
        if self.expected is None:
            return float(self._series.count) * self.bias
        else:
            return (float(self._series.count) * self.bias) / self.expected

    def __getitem__(self, item):
        try:
            return getattr(self, item)
        except AttributeError:
            raise KeyError("No such key: {}".format(item))

    @property
    def source(self):
        return int(self._series.bin1_id)

    @property
    def sink(self):
        return int(self._series.bin2_id)

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
    def attributes(self):
        try:
            base_attributes = self._series.index.tolist()
        except AttributeError:
            base_attributes = [a for a in dir(self._series) if not a.startswith('_')]

        return list(set(base_attributes + ['chromosome', 'bias', 'strand', 'start', 'end', 'weight']))

    @property
    def chromosome(self):
        return self._series.chrom

    @property
    def start(self):
        return self._series.start + 1

    @property
    def weight(self):
        try:
            return self._series.weight
        except AttributeError:
            return 1.0

    @property
    def bias(self):
        return self.weight

    @property
    def strand(self):
        try:
            return self._series.strand
        except AttributeError:
            return 1


class CoolerHic(RegionMatrixContainer, cooler.Cooler):
    def __init__(self, *args, **kwargs):
        mode = kwargs.pop("mode", 'r')
        tmpdir = kwargs.pop('tmpdir', None)

        if mode != 'r':
            warnings.warn("Mode {} not compatible with CoolerHic. "
                          "File will be opened in read-only mode and "
                          "changes will not be saved to file!")
        largs = list(args)
        at_resolution = None

        hic_file = largs[0]
        if "@" in largs[0]:
            hic_file, at_resolution = largs[0].split("@")

        if tmpdir is None or (isinstance(tmpdir, bool) and not tmpdir):
            self.tmp_file_name = None
        else:
            logger.info("Working in temporary directory...")
            if isinstance(tmpdir, bool):
                tmpdir = tempfile.gettempdir()
            else:
                tmpdir = os.path.expanduser(tmpdir)
            self.tmp_file_name = tmp_file_name(tmpdir, prefix='tmp_fanc', extension='.mcool')
            logger.info("Temporary file: {}".format(self.tmp_file_name))
            shutil.copyfile(hic_file, self.tmp_file_name)
            hic_file = self.tmp_file_name

        if at_resolution is None:
            largs[0] = hic_file
        else:
            largs[0] = hic_file + '::resolutions/{}'.format(str_to_int(at_resolution))

        cooler.Cooler.__init__(self, *largs, **kwargs)
        RegionMatrixContainer.__init__(self)
        self._mappability = None
        self._expected_value_and_marginals_cache = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return exc_type is None

    def close(self, remove_tmp=True):
        """
        Close this Juicer file and run exit operations.

        If file was opened with tmpdir in read-only mode:
        close file and delete temporary copy.

        :param remove_tmp: If False, does not delete temporary copy of file.
        """
        if self.tmp_file_name is not None and remove_tmp:
            os.remove(self.tmp_file_name)

    def _tuple_to_region(self, t, lazy_region=None):
        if lazy_region is not None:
            lazy_region._series = t
            lazy_region.ix = t.Index
            return lazy_region

        index = set(a for a in dir(t) if not a.startswith('_'))
        index.remove('chrom')
        index.remove('Index')
        try:
            index.remove('ix')
        except KeyError:
            pass

        kwargs = {name: getattr(t, name) for name in index}
        kwargs['chromosome'] = t.chrom
        try:
            kwargs['bias'] = t.weight
        except AttributeError:
            kwargs['bias'] = 1.
        kwargs['start'] = t.start + 1
        kwargs['ix'] = t.Index
        return GenomicRegion(**kwargs)

    def _region_iter(self, lazy=False, *args, **kwargs):
        lazy_region = LazyCoolerRegion(None) if lazy else None

        for t in self.bins()[:].itertuples():
            yield self._tuple_to_region(t, lazy_region=lazy_region)

    def _region_subset(self, region, lazy=False, *args, **kwargs):
        lazy_region = LazyCoolerRegion(None) if lazy else None

        query = "{}".format(region.chromosome)
        if region.start is not None and region.end is not None:
            query += ':{}-{}'.format(region.start - 1, region.end)

        try:
            df = self.bins().fetch(query)
        except ValueError:  # region might be beyond chromosome bounds
            start = region.start if region.start > 0 else 1
            chromosome_end = self.chromosome_lengths[region.chromosome]
            end = region.end if region.end <= chromosome_end else chromosome_end
            query = '{}:{}-{}'.format(region.chromosome, start - 1, end)
            df = self.bins().fetch(query)

        for t in df.itertuples():
            yield self._tuple_to_region(t, lazy_region=lazy_region)

    def _get_regions(self, item, *args, **kwargs):
        regions = []
        df = self.bins()[item]
        for t in df.itertuples():
            regions.append(self._tuple_to_region(t))

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

    def _tuple_to_edge(self, t, lazy_edge=None):
        if lazy_edge is None:
            index = set(a for a in dir(t) if not a.startswith('_'))
            index.remove('Index')
            index.remove('bin1_id')
            index.remove('bin2_id')
            index.remove('count')

            kwargs = {name: getattr(t, name) for name in index}
            kwargs['ix'] = int(t.Index)
            kwargs['source'] = int(t.bin1_id)
            kwargs['sink'] = int(t.bin2_id)
            kwargs['weight'] = float(t.count)
            kwargs['source_node'] = self.regions[kwargs['source']]
            kwargs['sink_node'] = self.regions[kwargs['sink']]

            return Edge(**kwargs)
        else:
            lazy_edge._series = t
            return lazy_edge

    def _edges_iter(self, lazy=False, *args, **kwargs):
        lazy_edge = LazyCoolerEdge(None, self) if lazy else None

        chromosomes = self.chromosomes()
        for i in range(len(chromosomes)):
            for j in range(i, len(chromosomes)):
                df = cooler.Cooler.matrix(self, as_pixels=True, balance=False).fetch(chromosomes[i],
                                                                                     chromosomes[j])
                for t in df.itertuples():
                    yield self._tuple_to_edge(t, lazy_edge=lazy_edge)

    def _edges_subset(self, key=None, row_regions=None, col_regions=None,
                      lazy=False, *args, **kwargs):
        lazy_edge = LazyCoolerEdge(None, self) if lazy else None

        row_start, row_end = self._min_max_region_ix(row_regions)
        col_start, col_end = self._min_max_region_ix(col_regions)

        dfs = []
        if col_start > row_end:  # no overlap, upper triangular
            df = cooler.Cooler.matrix(self, as_pixels=True, balance=False)[row_start:row_end + 1, col_start:col_end + 1]
            dfs.append(df)
        elif row_start > col_end:  # no overlap, lower triangular
            df = cooler.Cooler.matrix(self, as_pixels=True, balance=False)[col_start:col_end + 1, row_start:row_end + 1]
            dfs.append(df)
        else:  # overlapping
            overlap_start, overlap_end = max(row_start, col_start), min(row_end, col_end)
            df = cooler.Cooler.matrix(self, as_pixels=True,
                                      balance=False)[overlap_start:overlap_end + 1, overlap_start:overlap_end + 1]
            dfs.append(df)

            if col_start < overlap_start:  # lower triangular
                df = cooler.Cooler.matrix(self, as_pixels=True,
                                          balance=False)[col_start:overlap_start, row_start:row_end + 1]
                dfs.append(df)

            if col_end > overlap_end:  # upper triangular
                df = cooler.Cooler.matrix(self, as_pixels=True,
                                          balance=False)[row_start:row_end + 1, overlap_end + 1:col_end + 1]
                dfs.append(df)

            if row_start < overlap_start:  # upper triangular
                df = cooler.Cooler.matrix(self, as_pixels=True,
                                          balance=False)[row_start: overlap_start, overlap_start:overlap_end + 1]
                dfs.append(df)

            if row_end > overlap_end:  # lower triangular
                df = cooler.Cooler.matrix(self, as_pixels=True,
                                          balance=False)[overlap_start:overlap_end + 1, overlap_end + 1: row_end + 1]
                dfs.append(df)

        for i, df in enumerate(dfs):
            for t in df.itertuples():
                yield self._tuple_to_edge(t, lazy_edge=lazy_edge)

    def _edges_getitem(self, item, *args, **kwargs):
        edges = []
        df = self.pixels()[item]
        for t in df.itertuples():
            edges.append(self._tuple_to_edge(t))

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

        mappable = [False] * len(list(self.regions(region)))
        for edge in self.edges(region, lazy=True):
            mappable[edge.source] = True
            mappable[edge.sink] = True
        self._mappability = mappable
        return mappable

    def bias_vector(self):
        return np.array([r.bias if np.isfinite(r.bias) else 0 for r in self.regions(lazy=True)])

    def expected_values_and_marginals(self, selected_chromosome=None, norm=True,
                                      *args, **kwargs):
        if self._expected_value_and_marginals_cache is not None:
            logger.debug("Using cached expected values")
            intra_expected, chromosome_intra_expected, \
            inter_expected, marginals, valid = self._expected_value_and_marginals_cache
        else:
            intra_expected, chromosome_intra_expected, \
            inter_expected, marginals, valid = RegionMatrixContainer.expected_values_and_marginals(self, norm=norm,
                                                                                                   *args, **kwargs)
            self._expected_value_and_marginals_cache = intra_expected, chromosome_intra_expected, \
                                                       inter_expected, marginals, valid

        if selected_chromosome is not None:
            return chromosome_intra_expected[selected_chromosome], marginals, valid

        return intra_expected, chromosome_intra_expected, inter_expected, marginals, valid
