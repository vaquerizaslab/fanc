import logging
import struct
import zlib

import numpy as np
from genomic_regions import GenomicRegion

from ..regions import Genome
from ..hic import Hic
from ..pairs import ReadPairs
from ..matrix import RegionMatrixContainer, Edge
from ..config import config
from ..tools.files import tmp_file_name
from ..tools.general import str_to_int

import os
import subprocess
import warnings
import tempfile
import gzip
import shutil

from collections import defaultdict

logger = logging.getLogger(__name__)


def is_juicer(file_name):
    try:
        if "@" in file_name and not os.path.exists(file_name):
            fields = file_name.split("@")
            if len(fields) == 2:
                hic_file, _ = fields
            elif len(fields) == 3:
                hic_file, _, _ = fields
            else:
                raise ValueError("Too many fields for Juicer '@' notation!")
            file_name = hic_file
        with open(file_name, 'rb') as req:
            try:
                magic_string = struct.unpack('<3s', req.read(3))[0]
            except struct.error:
                return False

            if magic_string != b"HIC":
                return False
            else:
                return True
    except KeyError:
        return False


def convert_juicer_to_hic(juicer_file, genome_file, resolution, juicer_tools_jar_path=None,
                          norm='NONE', output_file=None, inter_chromosomal=True,
                          chromosomes=None):
    if juicer_tools_jar_path is None:
        juicer_tools_jar_path = config.juicer_tools_jar_path
    if juicer_tools_jar_path is None:
        raise ValueError("Just must provide the juicer tools jar path or set it in fanc.conf!")

    hic = Hic(file_name=output_file, mode='w')

    # regions
    logger.info("Building genome {}".format(genome_file))
    genome = Genome.from_string(genome_file)
    regions = genome.get_regions(resolution, chromosomes=chromosomes)
    hic.add_regions(regions)
    genome.close()
    regions.close()

    if chromosomes is None:
        chromosomes = hic.chromosomes()
    chromosome_bins = hic.chromosome_bins

    logger.info("Extracting edges from hic file")
    nan_counter = 0
    for i in range(len(chromosomes)):
        chromosome1 = chromosomes[i]
        offset1 = chromosome_bins[chromosome1][0]

        for j in range(i, len(chromosomes)):
            if i != j and not inter_chromosomal:
                continue
            chromosome2 = chromosomes[j]
            offset2 = chromosome_bins[chromosome2][0]

            logger.info("{} -- {}".format(chromosome1, chromosome2))

            juicer_command = ['java', '-jar', juicer_tools_jar_path,
                              'dump', 'observed', norm, juicer_file,
                              chromosome1, chromosome2, 'BP', str(resolution)]

            juicer_process = subprocess.Popen(juicer_command, stdout=subprocess.PIPE)

            for line in juicer_process.stdout:
                fields = line.rstrip().split()

                try:
                    start, end, weight = int(fields[0]), int(fields[1]), float(fields[2])
                except ValueError:
                    continue

                start_ix = int(start / resolution) + offset1
                end_ix = int(end / resolution) + offset2
                if start_ix > end_ix:
                    start_ix, end_ix = end_ix, start_ix

                if not np.isfinite(weight):
                    nan_counter += 1
                    continue

                hic.add_edge([start_ix, end_ix, weight], check_nodes_exist=False)

    hic.flush()

    if nan_counter > 0:
        logger.warning("{} contacts could not be imported, "
                       "because they had non-finite values.".format(nan_counter))

    return hic


def to_juicer(pairs, juicer_file, juicer_tools_jar_path=None,
              resolutions=None, fragment_map=False, tmp=False,
              verbose=False):
    tmpdir = None if not tmp else tempfile.mkdtemp()
    try:
        if juicer_tools_jar_path is None:
            juicer_tools_jar_path = config.juicer_tools_jar_path
        if juicer_tools_jar_path is None:
            raise ValueError("Just must provide the juicer tools jar path or set it in fanc.conf!")

        # check that juicer tools can be executed
        if not os.path.exists(juicer_tools_jar_path):
            raise ValueError("{} does not exist or is unreadable".format(juicer_tools_jar_path))

        help_command = ['java', '-jar', juicer_tools_jar_path, '-h']
        logger.debug("Testing juicer: {}".format(" ".join(help_command)))
        res = subprocess.call(help_command)
        if res != 0:
            raise RuntimeError("juicer had nonzero exit status ({})! "
                               "It looks like your .jar path is either wrong ".format(res))

        if isinstance(pairs, ReadPairs):
            pairs = [pairs]

        chromosomes = pairs[0].chromosomes()

        with tempfile.NamedTemporaryFile(suffix='.chrom.sizes', mode='w', delete=False) as chrom_sizes_file:
            logger.info("Writing chromosome sizes")
            for chromosome, size in pairs[0].chromosome_lengths.items():
                chrom_sizes_file.write("{}\t{}\n".format(chromosome, size))
            chrom_sizes_file.flush()

            # make restriction site file
            with tempfile.NamedTemporaryFile(prefix="re_sites_", suffix='.txt', mode='w', dir=tmpdir) as re_file:
                logger.info("Making restriction map")
                for chromosome in chromosomes:
                    re_sites = []
                    for region in pairs[0].regions(chromosome, lazy=True):
                        re_sites.append(str(region.end))
                    re_file.write("{} {}\n".format(chromosome, " ".join(re_sites)))
                re_file.flush()

                logger.info("Writing pairs")
                with tempfile.NamedTemporaryFile(prefix='reads_', suffix='.txt.gz', dir=tmpdir) as f:
                    f.flush()
                    fgz = gzip.GzipFile(mode='wb', fileobj=f)

                    for chri, chromosome1 in enumerate(chromosomes):
                        for chrj in range(chri, len(chromosomes)):
                            chromosome2 = chromosomes[chrj]
                            for p in pairs:
                                for pair in p.pairs((chromosome1, chromosome2), lazy=True):
                                    if pair.left.fragment.chromosome <= pair.right.fragment.chromosome:
                                        fgz.write("{} {} {} {} {} {} {} {}\n".format(
                                            1 if pair.left.strand == -1 else 0,
                                            pair.left.fragment.chromosome,
                                            pair.left.position,
                                            pair.left.fragment.ix,
                                            1 if pair.right.strand == -1 else 0,
                                            pair.right.fragment.chromosome,
                                            pair.right.position,
                                            pair.right.fragment.ix
                                        ).encode('utf-8'))
                                    else:
                                        fgz.write("{} {} {} {} {} {} {} {}\n".format(
                                            1 if pair.right.strand == -1 else 0,
                                            pair.right.fragment.chromosome,
                                            pair.right.position,
                                            pair.right.fragment.ix,
                                            1 if pair.left.strand == -1 else 0,
                                            pair.left.fragment.chromosome,
                                            pair.left.position,
                                            pair.left.fragment.ix
                                        ).encode('utf-8'))
                    fgz.close()

                    f.flush()

                    pre_command = ['java', '-jar', juicer_tools_jar_path, 'pre']
                    if fragment_map:
                        pre_command += ['-f', re_file.name]
                    if resolutions is not None:
                        pre_command += ['-r', ",".join([str(r) for r in resolutions])]
                    if tmp:
                        pre_command += ['-t', tmpdir]
                    if verbose:
                        pre_command.append('-v')
                    pre_command += [f.name, juicer_file, chrom_sizes_file.name]

                    logger.info("Making Juicer file: {}".format(" ".join(pre_command)))
                    res = subprocess.call(pre_command)
                    if res != 0:
                        raise RuntimeError("juicer pre had nonzero exit status ({})!".format(res))
        return JuicerHic(juicer_file)
    finally:
        if tmpdir is not None:
            shutil.rmtree(tmpdir)


def _read_cstr(f):
    """
    Copyright (c) 2016 Aiden Lab

    :param f: open binary file
    :return: str
    """
    buf = ""
    while True:
        b = f.read(1)
        b = b.decode('utf-8', 'backslashreplace')
        if b is None or b == '\0':
            return str(buf)
        else:
            buf = buf + b


class LazyJuicerEdge(object):
    def __init__(self, source, sink, matrix, **kwargs):
        self._matrix = matrix
        self._weight_field = 'weight'
        self.source = source
        self.sink = sink
        self.bias = 1.
        self.expected = None

    def __getattribute__(self, item):
        if item == '_weight_field' or item != self._weight_field:
            return object.__getattribute__(self, item)

        if self.expected is None:
            return object.__getattribute__(self, item) * self.bias
        else:
            return (object.__getattribute__(self, item) * self.bias) / self.expected

    def __getitem__(self, item):
        try:
            return getattr(self, item)
        except AttributeError:
            raise KeyError("No such key: {}".format(item))

    @property
    def source_node(self):
        return self._matrix.region_by_ix(self.source)

    @property
    def sink_node(self):
        return self._matrix.region_by_ix(self.sink)

    @property
    def field_names(self):
        return ['weight']


class JuicerHic(RegionMatrixContainer):
    def __init__(self, hic_file, resolution=None, mode='r', tmpdir=None, norm='KR'):
        RegionMatrixContainer.__init__(self)
        if '@' in hic_file:
            fields = hic_file.split("@")
            if len(fields) == 2:
                hic_file, at_resolution = fields
            elif len(fields) == 3:
                hic_file, at_resolution, norm = fields
            else:
                raise ValueError("Too many fields for Juicer '@' notation!")

            if resolution is not None and str_to_int(at_resolution) != resolution:
                raise ValueError("Conflicting resolution specifications: "
                                 "{} and {}".format(at_resolution, resolution))
            resolution = str_to_int(at_resolution)

        if tmpdir is None or (isinstance(tmpdir, bool) and not tmpdir):
            self.tmp_file_name = None
            self._hic_file = hic_file
        else:
            logger.info("Working in temporary directory...")
            if isinstance(tmpdir, bool):
                tmpdir = tempfile.gettempdir()
            else:
                tmpdir = os.path.expanduser(tmpdir)
            self.tmp_file_name = tmp_file_name(tmpdir, prefix='tmp_fanc', extension='_juicer.hic')
            logger.info("Temporary file: {}".format(self.tmp_file_name))
            shutil.copyfile(hic_file, self.tmp_file_name)
            self._hic_file = self.tmp_file_name

        bp_resolutions, _ = self.resolutions()
        if resolution is None:
            resolution = bp_resolutions[0]
            warnings.warn("No resolution chosen for Juicer Hic - using {}bp. "
                          "Specify a custom resolution using <.hic file>@<resolution>".format(resolution))
        if resolution not in bp_resolutions:
            raise ValueError("Resolution {} not supported ({})".format(resolution, bp_resolutions))

        if mode != 'r':
            warnings.warn("Mode {} not compatible with JuicerHic. "
                          "File will be opened in read-only mode and "
                          "changes will not be saved to file!")

        self._resolution = resolution
        self._normalisation = norm
        self._unit = 'BP'
        self._mappability = None

        if not is_juicer(hic_file):
            raise ValueError("File {} does not seem to be a .hic "
                             "file produced with juicer!".format(hic_file))

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

    @property
    def version(self):
        with open(self._hic_file, 'rb') as req:
            req.read(4)  # jump to version location
            return struct.unpack('<i', req.read(4))[0]

    def _master_index(self):
        with open(self._hic_file, 'rb') as req:
            req.read(8)  # jump to master index location
            return struct.unpack('<q', req.read(8))[0]

    @staticmethod
    def _skip_to_attributes(req):
        req.seek(0)

        # skip magic, version, master
        req.read(16)

        # skip genome
        while req.read(1).decode("utf-8") != '\0':
            pass

    @staticmethod
    def _skip_to_chromosome_lengths(req):
        JuicerHic._skip_to_attributes(req)

        # skip attributes
        n_attributes = struct.unpack('<i', req.read(4))[0]
        for _ in range(0, n_attributes):
            # skip key
            while req.read(1).decode("utf-8", "backslashreplace") != '\0':
                pass
            # skip value
            while req.read(1).decode("utf-8", "backslashreplace") != '\0':
                pass

    @staticmethod
    def _skip_to_resolutions(req):
        JuicerHic._skip_to_chromosome_lengths(req)

        n_chromosomes = struct.unpack('<i', req.read(4))[0]
        for _ in range(0, n_chromosomes):
            while req.read(1).decode("utf-8", "backslashreplace") != '\0':
                pass
            req.read(4)

    @staticmethod
    def _skip_to_chromosome_sites(req):
        JuicerHic._skip_to_resolutions(req)

        n_resolutions = struct.unpack('<i', req.read(4))[0]
        for _ in range(0, n_resolutions):
            req.read(4)

        n_fragment_resolutions = struct.unpack('<i', req.read(4))[0]
        for _ in range(0, n_fragment_resolutions):
            req.read(4)

    @property
    def juicer_attributes(self):
        with open(self._hic_file, 'rb') as req:
            JuicerHic._skip_to_attributes(req)
            attributes = {}
            n_attributes = struct.unpack('<i', req.read(4))[0]
            for _ in range(0, n_attributes):
                key = _read_cstr(req)
                value = _read_cstr(req)
                attributes[key] = value
        return attributes

    @property
    def chromosome_lengths(self):
        with open(self._hic_file, 'rb') as req:
            JuicerHic._skip_to_chromosome_lengths(req)

            chromosome_lengths = {}
            n_chromosomes = struct.unpack('<i', req.read(4))[0]
            for _ in range(0, n_chromosomes):
                name = _read_cstr(req)
                length = struct.unpack('<i', req.read(4))[0]
                chromosome_lengths[name] = length

        return chromosome_lengths

    def _chromosome_bins(self, *args, **kwargs):
        chromosome_lengths = self.chromosome_lengths
        bin_size = self.bin_size
        current_chromosome_start_bin = 0
        chromosome_bins = {}
        for chromosome in self.chromosomes():
            chromosome_length = chromosome_lengths[chromosome]
            end = int(np.ceil(chromosome_length / bin_size))
            chromosome_bins[chromosome] = [current_chromosome_start_bin, current_chromosome_start_bin + end]
            current_chromosome_start_bin += end
        return chromosome_bins

    def _all_chromosomes(self):
        with open(self._hic_file, 'rb') as req:
            JuicerHic._skip_to_chromosome_lengths(req)

            chromosomes = []
            n_chromosomes = struct.unpack('<i', req.read(4))[0]
            for _ in range(0, n_chromosomes):
                name = _read_cstr(req)
                req.read(4)
                chromosomes.append(name)

        return chromosomes

    def chromosomes(self):
        chromosomes = []
        for chromosome in self._all_chromosomes():
            if chromosome.lower() != 'all':
                chromosomes.append(chromosome)

        return chromosomes

    def resolutions(self):
        with open(self._hic_file, 'rb') as req:
            JuicerHic._skip_to_resolutions(req)

            resolutions = []
            n_resolutions = struct.unpack('<i', req.read(4))[0]
            for _ in range(0, n_resolutions):
                resolution = struct.unpack('<i', req.read(4))[0]
                resolutions.append(resolution)

            fragment_resolutions = []
            n_fragment_resolutions = struct.unpack('<i', req.read(4))[0]
            for _ in range(0, n_fragment_resolutions):
                resolution = struct.unpack('<i', req.read(4))[0]
                fragment_resolutions.append(resolution)

            return resolutions, fragment_resolutions

    @staticmethod
    def _skip_to_footer(req):
        req.seek(0)
        req.read(8)  # jump to master index location
        master_index = struct.unpack('<q', req.read(8))[0]

        # jump to footer location
        req.seek(master_index)
        req.read(4)  # skip number of bytes

    @staticmethod
    def _skip_to_expected_values(req):
        JuicerHic._skip_to_footer(req)

        n_entries = struct.unpack('<i', req.read(4))[0]
        for _ in range(n_entries):
            while req.read(1).decode("utf-8", "backslashreplace") != '\0':
                pass
            req.read(12)

    @staticmethod
    def _skip_to_normalised_expected_values(req):
        JuicerHic._skip_to_expected_values(req)

        n_vectors = struct.unpack('<i', req.read(4))[0]
        for _ in range(n_vectors):
            while req.read(1).decode("utf-8", "backslashreplace") != '\0':
                pass
            req.read(4)

            n_values = struct.unpack('<i', req.read(4))[0]
            for j in range(n_values):
                req.read(8)

            n_scaling_factors = struct.unpack('<i', req.read(4))[0]
            for _ in range(n_scaling_factors):
                req.read(12)

    @staticmethod
    def _skip_to_normalisation_vectors(req):
        JuicerHic._skip_to_normalised_expected_values(req)

        n_vectors = struct.unpack('<i', req.read(4))[0]
        for _ in range(n_vectors):
            while req.read(1).decode("utf-8", "backslashreplace") != '\0':
                pass
            while req.read(1).decode("utf-8", "backslashreplace") != '\0':
                pass
            req.read(4)

            n_values = struct.unpack('<i', req.read(4))[0]
            for j in range(n_values):
                req.read(8)

            n_scaling_factors = struct.unpack('<i', req.read(4))[0]
            for _ in range(n_scaling_factors):
                req.read(12)

    def _matrix_positions(self):
        """
        Copyright (c) 2016 Aiden Lab
        """

        with open(self._hic_file, 'rb') as req:
            JuicerHic._skip_to_footer(req)

            chromosome_pair_positions = {}
            n_entries = struct.unpack('<i', req.read(4))[0]
            for _ in range(n_entries):
                key = tuple(_read_cstr(req).split('_'))
                file_position = struct.unpack('<q', req.read(8))[0]
                req.read(4)  # skip size in bytes
                chromosome_pair_positions[key] = file_position
            return chromosome_pair_positions

    @staticmethod
    def _expected_value_vectors_from_pos(req, normalisation=None, unit='BP'):
        expected_values = defaultdict(list)
        scaling_factors = defaultdict(dict)

        n_vectors = struct.unpack('<i', req.read(4))[0]
        for _ in range(n_vectors):
            if normalisation:
                entry_normalisation = _read_cstr(req)
                entry_unit = _read_cstr(req)
            else:
                entry_unit = _read_cstr(req)
                entry_normalisation = 'NONE'

            bin_size = struct.unpack('<i', req.read(4))[0]

            ev = []
            n_values = struct.unpack('<i', req.read(4))[0]

            for j in range(n_values):
                v = struct.unpack('<d', req.read(8))[0]
                ev.append(v)

            if entry_unit == unit and (normalisation is None or entry_normalisation == normalisation):
                expected_values[bin_size] = ev

            sf = dict()
            n_scaling_factors = struct.unpack('<i', req.read(4))[0]
            for _ in range(n_scaling_factors):
                chromosome_index = struct.unpack('<i', req.read(4))[0]
                f = struct.unpack('<d', req.read(8))[0]
                sf[chromosome_index] = f

            if normalisation is None or entry_normalisation == normalisation:
                scaling_factors[bin_size] = sf

        return expected_values, scaling_factors

    def expected_value_vector(self, chromosome, normalisation=None, resolution=None):
        if normalisation is None:
            normalisation = self._normalisation

        if resolution is None:
            resolution = self._resolution

        chromosome_ix = self._all_chromosomes().index(chromosome)

        if normalisation == 'NONE':
            vectors, scaling_factors = self.expected_value_vectors()
        else:
            vectors, scaling_factors = self.normalised_expected_value_vectors(normalisation)

        try:
            sf = scaling_factors[resolution][chromosome_ix]
        except KeyError:
            warnings.warn("Cannot find an expected value scaling factor "
                          "for {}, setting to 0.".format(chromosome))
            sf = 0.0
        return np.array(vectors[resolution]) / sf

    def expected_value_vectors(self):
        with open(self._hic_file, 'rb') as req:
            JuicerHic._skip_to_expected_values(req)

            return JuicerHic._expected_value_vectors_from_pos(req)

    def normalised_expected_value_vectors(self, normalisation=None):
        if normalisation is None:
            normalisation = self._normalisation

        with open(self._hic_file, 'rb') as req:
            JuicerHic._skip_to_normalised_expected_values(req)

            return JuicerHic._expected_value_vectors_from_pos(req, normalisation=normalisation)

    def expected_values(self, selected_chromosome=None, norm=True, *args, **kwargs):
        def _fill_norm_vector(v, n, fill):
            if len(v) < n:
                v = np.append(v, [fill for _ in range(n - len(v))])
            return v

        cb = self.chromosome_bins
        if selected_chromosome is not None:
            if norm:
                expected_chromosome = self.expected_value_vector(selected_chromosome, self._normalisation)
            else:
                expected_chromosome = self.expected_value_vector(selected_chromosome, 'NONE')

            chromosome_size_in_bins = cb[selected_chromosome][1] - cb[selected_chromosome][0]
            expected_chromosome = _fill_norm_vector(expected_chromosome, chromosome_size_in_bins, 0)

            return expected_chromosome

        intra_expected = dict()
        for chromosome in self.chromosomes():
            if norm:
                expected_chromosome = self.expected_value_vector(chromosome, self._normalisation)
            else:
                expected_chromosome = self.expected_value_vector(chromosome, 'NONE')

            chromosome_size_in_bins = cb[chromosome][1] - cb[chromosome][0]
            expected_chromosome = _fill_norm_vector(expected_chromosome, chromosome_size_in_bins, 0)

            intra_expected[chromosome] = expected_chromosome

        return None, intra_expected, None

    def normalisation_vector(self, chromosome, normalisation=None, resolution=None, unit=None):
        if resolution is None:
            resolution = self._resolution

        if normalisation is None:
            normalisation = self._normalisation

        if normalisation == 'NONE':
            bins = int(np.ceil(self.chromosome_lengths[chromosome]/resolution))
            return [1.0] * bins

        if unit is None:
            unit = self._unit

        chromosomes = self.chromosomes()
        chromosome_index = chromosomes.index(chromosome) + 1

        existing_normalisations = set()
        with open(self._hic_file, 'rb') as req:
            JuicerHic._skip_to_normalisation_vectors(req)

            n_entries = struct.unpack('<i', req.read(4))[0]
            for _ in range(n_entries):
                entry_normalisation = _read_cstr(req)
                entry_chromosome_index = struct.unpack('<i', req.read(4))[0]
                entry_unit = _read_cstr(req)
                entry_resolution = struct.unpack('<i', req.read(4))[0]
                file_position = struct.unpack('<q', req.read(8))[0]
                req.read(4)  # skip size in bytes
                existing_normalisations.add(entry_normalisation)

                if (entry_chromosome_index == chromosome_index and
                        entry_normalisation == normalisation and
                        entry_resolution == resolution and
                        entry_unit == unit):
                    req.seek(file_position)
                    vector = []
                    n_values = struct.unpack('<i', req.read(4))[0]
                    for _ in range(n_values):
                        v = struct.unpack('<d', req.read(8))[0]
                        vector.append(v)

                    return vector
        
        if normalisation not in existing_normalisations:
            raise ValueError("Cannot find normalisation '{}'".format(normalisation))
        
        warnings.warn("Cannot find normalisation vector for "
                      "chromosome: {chromosome}, normalisation: {norm}, "
                      "resolution: {res}, unit: {unit}. This could "
                      "indicate that {norm} normalisation did not "
                      "work for this chromosome. Will return NaN instead.".format(
                          chromosome=chromosome, norm=normalisation, 
                          res=resolution, unit=unit
                          )
                      )
        c_bins = self.chromosome_bins[chromosome]
        return np.repeat(np.nan, c_bins[1] - c_bins[0])

    def region_by_ix(self, ix):
        chromosome_lengths = self.chromosome_lengths

        offset_ix = 0
        remaining_ix = ix
        chromosomes = self.chromosomes()
        current_chromosome = None
        for chromosome in chromosomes:
            current_chromosome = chromosome
            chromosome_length = chromosome_lengths[chromosome]
            if chromosome.lower() == 'all':
                continue

            ixs = int(np.ceil(chromosome_length / self._resolution))
            if remaining_ix > ixs:
                offset_ix += ixs
                remaining_ix -= ixs
            else:
                break

        region_ix = offset_ix + remaining_ix
        start = remaining_ix * self._resolution + 1
        return GenomicRegion(chromosome=current_chromosome, start=start,
                             end=min(start + self._resolution - 1,
                                     chromosome_lengths[current_chromosome]),
                             ix=region_ix)

    def _chromosome_ix_offset(self, target_chromosome):
        chromosome_lengths = self.chromosome_lengths
        if target_chromosome not in chromosome_lengths:
            raise ValueError("Chromosome {} not in matrix.".format(target_chromosome))

        offset_ix = 0
        chromosomes = self.chromosomes()
        for chromosome in chromosomes:
            chromosome_length = chromosome_lengths[chromosome]
            if chromosome.lower() == 'all':
                continue

            if target_chromosome == chromosome:
                return offset_ix

            ixs = int(np.ceil(chromosome_length / self._resolution))
            offset_ix += ixs

    def _region_start(self, region):
        region = self._convert_region(region)
        offset_ix = self._chromosome_ix_offset(region.chromosome)
        region_start = region.start if region.start is not None else 1
        ix = int((region_start - 1) / self._resolution)
        start = self._resolution * ix + 1
        return offset_ix + ix, start

    def _region_iter(self, *args, **kwargs):
        chromosome_lengths = self.chromosome_lengths

        chromosomes = self.chromosomes()
        for chromosome in chromosomes:
            chromosome_length = chromosome_lengths[chromosome]
            if chromosome.lower() == 'all':
                continue

            offset_ix = self._chromosome_ix_offset(chromosome)

            try:
                norm = self.normalisation_vector(chromosome)
            except ValueError:
                warnings.warn("Cannot find bias vector ({}, {} resolution) for chromosome {}."
                              "Continuing by masking the corresponding regions. If this is unexpected, "
                              "try choosing another normalisation method!".format(self._normalisation,
                                                                                  self._resolution, chromosome))
                norm = None

            for i, start in enumerate(range(1, chromosome_length, self._resolution)):
                if norm is None or np.isnan(norm[i]) or norm[i] == 0:
                    valid = False
                    bias = 1.0
                else:
                    valid = True
                    try:
                        bias = 1/norm[i]
                    except IndexError:
                        bias = 1.0

                if np.isnan(bias):
                    bias = 1.0

                end = min(start + self._resolution - 1, chromosome_length)
                region = GenomicRegion(chromosome=chromosome, start=start,
                                       end=end, bias=bias, valid=valid,
                                       ix=int(offset_ix + i))
                yield region

    def _region_subset(self, region, *args, **kwargs):
        subset_ix, subset_start = self._region_start(region)

        cl = self.chromosome_lengths[region.chromosome]
        norm = self.normalisation_vector(region.chromosome)
        for i, start in enumerate(range(subset_start, min(cl, region.end), self._resolution)):
            end = min(start + self._resolution - 1, cl, region.end)
            bias_ix = int(start / self._resolution)

            if bias_ix >= len(norm):
                break
            
            if np.isnan(norm[bias_ix]) or norm[bias_ix] == 0:
                valid = False
                bias = 1.0
            else:
                valid = True
                try:
                    bias = 1/norm[bias_ix]
                except IndexError:
                    bias = 1.0

            if np.isnan(bias):
                bias = 1.0

            r = GenomicRegion(chromosome=region.chromosome, start=int(start),
                              end=int(end), bias=bias, valid=valid,
                              ix=int(subset_ix + i))
            yield r

    def _region_len(self):
        length = 0
        for chromosome, chromosome_length in self.chromosome_lengths.items():
            if chromosome.lower() == 'all':
                continue

            ixs = int(np.ceil(chromosome_length / self._resolution))
            length += ixs
        return length

    def _read_block(self, req, file_position, block_size_in_bytes):
        req.seek(file_position)
        block_compressed = req.read(block_size_in_bytes)
        block = zlib.decompress(block_compressed)

        n_records = struct.unpack('<i', block[0:4])[0]
        if self.version < 7:
            for i in range(n_records):
                x = struct.unpack('<i', block[(12 * i + 4):(12 * i + 8)])[0]
                y = struct.unpack('<i', block[(12 * i + 8):(12 * i + 12)])[0]
                weight = struct.unpack('<f', block[(12 * i + 12):(12 * i + 16)])[0]
                yield x, y, weight
        else:
            x_offset = struct.unpack('<i', block[4:8])[0]
            y_offset = struct.unpack('<i', block[8:12])[0]
            use_short = not struct.unpack('<b', block[12:13])[0] == 0
            block_type = struct.unpack('<b', block[13:14])[0]
            index = 0

            if block_type == 1:
                row_count = struct.unpack('<h', block[14:16])[0]
                temp = 16
                for i in range(row_count):
                    y_raw = struct.unpack('<h', block[temp:(temp + 2)])[0]
                    temp += 2
                    y = y_raw + y_offset
                    col_count = struct.unpack('<h', block[temp:(temp + 2)])[0]
                    temp += 2
                    for j in range(col_count):
                        x_raw = struct.unpack('<h', block[temp:(temp + 2)])[0]
                        temp += 2
                        x = x_offset + x_raw
                        if not use_short:
                            weight = struct.unpack('<h', block[temp:(temp + 2)])[0]
                            temp += 2
                        else:
                            weight = struct.unpack('<f', block[temp:(temp + 4)])[0]
                            temp += 4
                        yield x, y, weight

                        index += 1
            elif block_type == 2:
                temp = 14
                n_points = struct.unpack('<i', block[temp:(temp + 4)])[0]
                temp += 4
                w = struct.unpack('<h', block[temp:(temp + 2)])[0]
                temp += 2
                for i in range(n_points):
                    row = int(i / w)
                    col = i - row * w
                    x = int(x_offset + col)
                    y = int(y_offset + row)
                    if not use_short:
                        weight = struct.unpack('<h', block[temp:(temp + 2)])[0]
                        temp += 2
                        if weight != -32768:
                            yield x, y, weight
                            index += 1
                    else:
                        weight = struct.unpack('<f', block[temp:(temp + 4)])[0]
                        temp += 4
                        if weight != 0x7fc00000:
                            yield x, y, weight
                            index = index + 1

    def _read_matrix(self, region1, region2):
        region1 = self._convert_region(region1)
        region2 = self._convert_region(region2)

        chromosomes = self._all_chromosomes()
        chromosome1_ix = chromosomes.index(region1.chromosome)
        chromosome2_ix = chromosomes.index(region2.chromosome)

        if chromosome1_ix > chromosome2_ix:
            region1, region2 = region2, region1
            chromosome1_ix, chromosome2_ix = chromosome2_ix, chromosome1_ix

        region1_chromosome_offset = self._chromosome_ix_offset(region1.chromosome)
        region2_chromosome_offset = self._chromosome_ix_offset(region2.chromosome)

        try:
            matrix_file_position = self._matrix_positions()[(str(chromosome1_ix), str(chromosome2_ix))]
        except KeyError:
            return

        with open(self._hic_file, 'rb') as req:
            req.seek(matrix_file_position)
            req.read(8)  # skip chromosome index

            block_bin_count = None
            block_column_count = None
            block_map = dict()
            n_resolutions = struct.unpack('<i', req.read(4))[0]
            for i in range(n_resolutions):
                unit = _read_cstr(req)
                req.read(20)  # skip reserved but unused fields

                bin_size = struct.unpack('<i', req.read(4))[0]
                if unit == self._unit and bin_size == self._resolution:
                    block_bin_count = struct.unpack('<i', req.read(4))[0]
                    block_column_count = struct.unpack('<i', req.read(4))[0]

                    n_blocks = struct.unpack('<i', req.read(4))[0]
                    for b in range(n_blocks):
                        block_number = struct.unpack('<i', req.read(4))[0]
                        file_position = struct.unpack('<q', req.read(8))[0]
                        block_size_in_bytes = struct.unpack('<i', req.read(4))[0]
                        block_map[block_number] = (file_position, block_size_in_bytes)
                else:
                    req.read(8)

                    n_blocks = struct.unpack('<i', req.read(4))[0]
                    for b in range(n_blocks):
                        req.read(16)

            if block_bin_count is None or block_column_count is None:
                raise ValueError("Matrix data for {} {} not found!".format(self._resolution, self._unit))

            region1_bins = int(region1.start / self._resolution), int(region1.end / self._resolution) + 1
            region2_bins = int(region2.start / self._resolution), int(region2.end / self._resolution) + 1

            col1, col2 = int(region1_bins[0] / block_bin_count), int(region1_bins[1] / block_bin_count)
            row1, row2 = int(region2_bins[0] / block_bin_count), int(region2_bins[1] / block_bin_count)

            blocks = set()
            for r in range(row1, row2 + 1):
                for c in range(col1, col2 + 1):
                    block_number = r * block_column_count + c
                    blocks.add(block_number)

            if region1.chromosome == region2.chromosome:
                for r in range(col1, col2 + 1):
                    for c in range(row1, row2 + 1):
                        block_number = r * block_column_count + c
                        blocks.add(block_number)

            for block_number in blocks:
                try:
                    file_position, block_size_in_bytes = block_map[block_number]

                    for x, y, weight in self._read_block(req, file_position, block_size_in_bytes):
                        if x < y:
                            if region1_bins[0] <= x < region1_bins[1] - 1 and region2_bins[0] <= y < region2_bins[1] - 1:
                                yield x + region1_chromosome_offset, y + region2_chromosome_offset, weight
                            elif region1.chromosome == region2.chromosome:
                                if region1_bins[0] <= y < region1_bins[1] - 1 and region2_bins[0] <= x < region2_bins[1] - 1:
                                    yield x + region1_chromosome_offset, y + region2_chromosome_offset, weight
                        else:
                            if region1_bins[0] <= x < region1_bins[1] - 1 and region2_bins[0] <= y < region2_bins[1] - 1:
                                yield y + region2_chromosome_offset, x + region1_chromosome_offset, weight

                except KeyError:
                    logger.debug("Could not find block {}".format(block_number))

    def _edges_subset(self, key=None, row_regions=None, col_regions=None,
                      lazy=False, *args, **kwargs):

        if row_regions[0].chromosome != row_regions[-1].chromosome:
            raise ValueError("Cannot subset rows across multiple chromosomes!")

        if col_regions[0].chromosome != col_regions[-1].chromosome:
            raise ValueError("Cannot subset columns across multiple chromosomes!")

        regions_by_ix = {}
        for region in row_regions + col_regions:
            regions_by_ix[region.ix] = region

        row_span = GenomicRegion(chromosome=row_regions[0].chromosome,
                                 start=row_regions[0].start,
                                 end=row_regions[-1].end)

        col_span = GenomicRegion(chromosome=col_regions[0].chromosome,
                                 start=col_regions[0].start,
                                 end=col_regions[-1].end)

        if not lazy:
            for x, y, weight in self._read_matrix(row_span, col_span):
                if x > y:
                    x, y = y, x
                yield Edge(source=regions_by_ix[x],
                           sink=regions_by_ix[y],
                           weight=weight)
        else:
            edge = LazyJuicerEdge(source=0, sink=0, weight=1.0, matrix=self)
            for x, y, weight in self._read_matrix(row_span, col_span):
                if x > y:
                    x, y = y, x
                edge.source, edge.sink, edge.weight = x, y, weight
                yield edge

    def _edges_iter(self, *args, **kwargs):
        chromosomes = self.chromosomes()
        for ix1 in range(len(chromosomes)):
            chromosome1 = chromosomes[ix1]
            for ix2 in range(ix1, len(chromosomes)):
                chromosome2 = chromosomes[ix2]

                key = chromosome1, chromosome2
                row_regions, col_regions = self._key_to_regions(key)

                if isinstance(row_regions, GenomicRegion):
                    row_regions = [row_regions]
                else:
                    row_regions = list(row_regions)

                if isinstance(col_regions, GenomicRegion):
                    col_regions = [col_regions]
                else:
                    col_regions = list(col_regions)

                for edge in self._edges_subset(key=key, row_regions=row_regions,
                                               col_regions=col_regions, *args, **kwargs):
                    yield edge

    @property
    def bin_size(self):
        return self._resolution

    def mappable(self, region=None):
        """
        Get the mappability vector of this matrix.
        """
        if region is None:
            if self._mappability is not None:
                return self._mappability

            mappable = []
            for chromosome, chromosome_length in self.chromosome_lengths.items():
                if chromosome.lower() == 'all':
                    continue

                norm = self.normalisation_vector(chromosome)
                for i, _ in enumerate(range(1, chromosome_length, self._resolution)):
                    if np.isnan(norm[i]):
                        mappable.append(False)
                    else:
                        mappable.append(True)

            self._mappability = mappable

            return mappable
        else:
            return np.array([r.valid for r in self.regions(region, lazy=True)])

    def bias_vector(self):
        return np.array([r.bias for r in self.regions(lazy=True)])
