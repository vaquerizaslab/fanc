from __future__ import division

import os
import copy
import gzip
import logging
import multiprocessing as mp
from queue import Empty
import threading
import uuid
from abc import abstractmethod, ABCMeta
from bisect import bisect_right
from builtins import object
from collections import defaultdict
from timeit import default_timer as timer

import msgpack
import numpy as np
import pysam
import tables as t
from future.utils import with_metaclass, viewitems

from genomic_regions import GenomicRegion
from .regions import genome_regions
from .matrix import Edge, RegionPairsTable
from .hic import Hic
from .config import config
from .general import MaskFilter, MaskedTable, Mask
from .tools.general import RareUpdateProgressBar, add_dict, find_alignment_match_positions, WorkerMonitor
from .tools.sambam import natural_cmp

logger = logging.getLogger(__name__)


def generate_pairs(sam1_file, sam2_file, regions,
                   restriction_enzyme=None, read_filters=(),
                   output_file=None, check_sorted=True,
                   threads=1, batch_size=10000000):
    regions = genome_regions(regions, restriction_enzyme=restriction_enzyme)

    sb = SamBamReadPairGenerator(sam1_file, sam2_file, check_sorted=check_sorted)
    for f in read_filters:
        sb.add_filter(f)

    pairs = ReadPairs(file_name=output_file, mode='w')

    pairs.add_regions(regions, preserve_attributes=False)
    pairs.add_read_pairs(sb, threads=threads, batch_size=batch_size)

    return pairs


class Monitor(WorkerMonitor):
    def __init__(self, value=0):
        WorkerMonitor.__init__(self, value=value)
        self.generating_pairs_lock = threading.Lock()

        with self.generating_pairs_lock:
            self.generating_pairs = True

    def set_generating_pairs(self, value):
        with self.generating_pairs_lock:
            self.generating_pairs = value

    def is_generating_pairs(self):
        with self.generating_pairs_lock:
            return self.generating_pairs


def _fragment_info_worker(monitor, input_queue, output_queue, fi, fe):
    worker_uuid = uuid.uuid4()
    logger.debug("Starting fragment info worker {}".format(worker_uuid))

    while True:
        # wait for input
        monitor.set_worker_idle(worker_uuid)
        logger.debug("Worker {} waiting for input".format(worker_uuid))
        read_pairs = input_queue.get(True)
        monitor.set_worker_busy(worker_uuid)
        logger.debug('Worker {} reveived input!'.format(worker_uuid))
        read_pairs = msgpack.loads(read_pairs)

        fragment_infos = []
        skipped_counter = 0
        for (chrom1, pos1, flag1), (chrom2, pos2, flag2) in read_pairs:
            chrom1 = chrom1.decode() if isinstance(chrom1, bytes) else chrom1
            chrom2 = chrom2.decode() if isinstance(chrom2, bytes) else chrom2

            try:
                pos_ix1 = bisect_right(fe[chrom1], pos1)
                pos_ix2 = bisect_right(fe[chrom2], pos2)

                f_ix1, f_chromosome_ix1, f_start1, f_end1 = fi[chrom1][pos_ix1]
                f_ix2, f_chromosome_ix2, f_start2, f_end2 = fi[chrom2][pos_ix2]

                r_strand1 = -1 if flag1 & 16 else 1
                r_strand2 = -1 if flag2 & 16 else 1

                fragment_infos.append(
                    ((pos1, r_strand1, f_ix1, f_chromosome_ix1, f_start1, f_end1),
                     (pos2, r_strand2, f_ix2, f_chromosome_ix2, f_start2, f_end2))
                )
            except (KeyError, IndexError):
                skipped_counter += 1
        logger.debug("Worker {} skipped {} pairs".format(worker_uuid, skipped_counter))
        output_queue.put(msgpack.dumps(fragment_infos))
        del read_pairs


def _read_pairs_worker(read_pairs, input_queue, monitor, batch_size=100000):
    logger.debug("Starting read pairs worker")
    try:
        read_pairs_batch = []
        for read1, read2 in read_pairs:
            read_pairs_batch.append((
                (read1.reference_name, read1.pos, read1.flag),
                (read2.reference_name, read2.pos, read2.flag)
            ))
            if len(read_pairs_batch) >= batch_size:
                logger.debug("Submitting read pair batch ({}) to input queue".format(batch_size))
                input_queue.put(msgpack.dumps(read_pairs_batch))
                read_pairs_batch = []
                monitor.increment()
        if len(read_pairs_batch) > 0:
            logger.debug("Submitting read pair batch ({}) to input queue".format(batch_size))
            input_queue.put(msgpack.dumps(read_pairs_batch))
            monitor.increment()
    finally:
        monitor.set_generating_pairs(False)
    logger.debug("Terminating read pairs worker")


class MinimalRead(object):
    def __init__(self, chromosome, position, strand):
        self.chromosome = chromosome
        self.reference_name = chromosome
        self.position = position
        self.strand = strand
        self.flag = 0 if strand == '+' or strand == 1 else -1
        self.pos = position


class ReadPairGenerator(object):
    def __init__(self):
        self.filters = []
        self._filter_stats = defaultdict(int)
        self._total_pairs = 0
        self._valid_pairs = 0

        unmapped_filter = UnmappedFilter(mask=Mask('unmapped', 'Mask unmapped reads', ix=len(self.filters)))
        self.add_filter(unmapped_filter)
        self._unmapped_filter_ix = len(self.filters) - 1

    def _iter_read_pairs(self, *args, **kwargs):
        raise NotImplementedError("Class must override iter_read_pairs")

    def add_filter(self, read_filter):
        if not isinstance(read_filter, ReadFilter):
            raise ValueError("argument must be an instance of class ReadFilter!")
        self.filters.append(read_filter)

    def stats(self):
        filter_names = []
        for i, f in enumerate(self.filters):
            if hasattr(f, 'mask_name') and f.mask_name is not None:
                filter_names.append(f.mask_name)
            else:
                filter_names.append('filter_{}'.format(i))

        stats = dict()
        for i, count in self._filter_stats.items():
            stats[filter_names[i]] = count
        stats['unmasked'] = self._valid_pairs
        stats['total'] = self._total_pairs
        return stats

    def __iter__(self):
        self._filter_stats = defaultdict(int)
        self._total_pairs = 0
        self._valid_pairs = 0
        for (read1, read2) in self._iter_read_pairs():
            valid_reads = True
            for i, f in enumerate(self.filters):
                if not f.valid_read(read1) or not f.valid_read(read2):
                    self._filter_stats[i] += 1
                    valid_reads = False
            if valid_reads:
                yield (read1, read2)
                self._valid_pairs += 1
            self._total_pairs += 1


class TxtReadPairGenerator(ReadPairGenerator):
    def __init__(self, valid_pairs_file, sep=None,
                 chr1_field=1, pos1_field=2, strand1_field=3,
                 chr2_field=4, pos2_field=5, strand2_field=6):
        ReadPairGenerator.__init__(self)
        self._file_name = valid_pairs_file
        self.sep = sep
        self.chr1_field = chr1_field
        self.pos1_field = pos1_field
        self.strand1_field = strand1_field
        self.chr2_field = chr2_field
        self.pos2_field = pos2_field
        self.strand2_field = strand2_field

        if self._file_name.endswith('.gz') or self._file_name.endswith('gzip'):
            self._open_file = gzip.open
        else:
            self._open_file = open

        if not os.path.exists(valid_pairs_file):
            raise ValueError("File {} does not exist!".format(valid_pairs_file))

        # check that this file is valid:
        with self._open_file(valid_pairs_file, 'rt') as f:
            for line in f:
                line = line.rstrip()
                if line.startswith('#') or line == '':
                    continue

                fields = line.split(sep)

                # find max field
                max_field_ix = 0
                for field_number in {chr1_field, chr2_field,
                                     pos1_field, pos2_field,
                                     strand1_field, strand2_field}:
                    if field_number is not None:
                        max_field_ix = max(max_field_ix, field_number)

                if len(fields) < max_field_ix + 1:
                    raise ValueError("Not enough fields ({}) in file {}".format(len(fields),
                                                                                valid_pairs_file))

                # ensure we can transform pos fields to int
                int(fields[pos1_field])
                int(fields[pos2_field])
                # ensure strand fields are valid if specified
                if strand1_field is not None and fields[strand1_field] not in {'+', '-', '.', '-1', '1', '+1'}:
                    raise ValueError("Cannot read strand1 field!")
                if strand2_field is not None and fields[strand2_field] not in {'+', '-', '.', '-1', '1', '+1'}:
                    raise ValueError("Cannot read strand2 field!")

                break

    def _iter_read_pairs(self, *args, **kwargs):
        with self._open_file(self._file_name, 'rt') as f:
            for line in f:
                line = line.rstrip()
                if line == '' or line.startswith('#'):
                    continue
                fields = line.split(self.sep)

                strand1 = fields[self.strand1_field] if self.strand1_field is not None else '+'
                strand2 = fields[self.strand2_field] if self.strand2_field is not None else '+'

                read1 = MinimalRead(chromosome=fields[self.chr1_field],
                                    position=int(fields[self.pos1_field]),
                                    strand=strand1)
                read2 = MinimalRead(chromosome=fields[self.chr2_field],
                                    position=int(fields[self.pos2_field]),
                                    strand=strand2)
                yield (read1, read2)


class HicProPairGenerator(TxtReadPairGenerator):
    def __init__(self, file_name):
        TxtReadPairGenerator.__init__(self, file_name, sep="\t",
                                      chr1_field=1, pos1_field=2, strand1_field=3,
                                      chr2_field=4, pos2_field=5, strand2_field=6)


class FourDNucleomePairGenerator(TxtReadPairGenerator):
    def __init__(self, pairs_file):
        if pairs_file.endswith('.gz') or pairs_file.endswith('gzip'):
            open_file = gzip.open
        else:
            open_file = open

        columns = dict()
        with open_file(pairs_file, 'rt') as f:
            for line_ix, line in enumerate(f):
                if line_ix == 0 and not line.startswith("## pairs format"):
                    raise ValueError("Not a 4D nucleome pairs format file."
                                     "Missing '## pairs format X.X' header line.")

                line = line.rstrip()
                if not line.startswith('#'):
                    raise ValueError("Pairs file does not contain a "
                                     "'#columns' entry in the header")

                if line.startswith('#columns:'):
                    _, columns_field = line.split(':')
                    for i, name in columns_field.split():
                        columns[name] = i

        TxtReadPairGenerator.__init__(self, pairs_file, sep=None,
                                      chr1_field=columns['chr1'],
                                      pos1_field=columns['pos1'],
                                      strand1_field=columns['strand1'] if 'strand1' in columns else None,
                                      chr2_field=columns['chr2'],
                                      pos2_field=columns['pos2'],
                                      strand2_field=columns['strand2'] if 'strand2' in columns else None,
                                      )


class SamBamReadPairGenerator(ReadPairGenerator):
    def __init__(self, sam_file1, sam_file2, check_sorted=True):
        ReadPairGenerator.__init__(self)
        self.sam_file1 = sam_file1
        self.sam_file2 = sam_file2
        self._check_sorted = check_sorted
        if not os.path.exists(self.sam_file1):
            raise ValueError("File {} does not exist!".format(self.sam_file1))
        if not os.path.exists(self.sam_file2):
            raise ValueError("File {} does not exist!".format(self.sam_file2))

    def _iter_read_pairs(self, *args, **kwargs):
        max_dist_same_locus = kwargs.get('max_dist_same_locus', 100)
        logger.info("Starting to generate read pairs from SAM")

        def _all_reads(iterator, last_read=None):
            reads = []
            if last_read is not None:
                if last_read.is_unmapped:
                    self._filter_stats[self._unmapped_filter_ix] += 1
                else:
                    reads.append(last_read)

            next_read = None
            try:
                next_read = next(iterator)
                while len(reads) == 0 or natural_cmp(next_read.qname.encode(), reads[0].qname.encode()) == 0:
                    if not next_read.is_unmapped:
                        reads.append(next_read)
                    else:
                        self._filter_stats[self._unmapped_filter_ix] += 1
                    next_read = next(iterator)
            except StopIteration:
                if len(reads) == 0:
                    raise
            return reads[0].qname.encode(), reads, next_read

        def _find_pair(reads1, reads2):
            """
            :return: read1, read2, is_chimeric
            """
            if len(reads1) == len(reads2) == 1:
                return reads1[0], reads2[0], False
            elif (len(reads1) == 1 and len(reads2) == 2) or (len(reads2) == 1 and len(reads1) == 2):
                if len(reads2) > len(reads1):
                    reads1, reads2 = reads2, reads1

                read2 = reads2[0]
                match_pos2 = find_alignment_match_positions(read2, longest=True)[0]
                if match_pos2 is None:
                    return None, None, True

                read1 = None
                same_locus = False
                for read in reads1:
                    if read.reference_id != read2.reference_id:
                        read1 = read
                    else:
                        match_pos1 = find_alignment_match_positions(read, longest=True)[0]
                        if match_pos1 is None:
                            return None, None, True

                        if min(abs(match_pos1[0] - match_pos2[1]),
                               abs(match_pos1[1] - match_pos2[0])) > max_dist_same_locus:
                            read1 = read
                        else:
                            same_locus = True

                if same_locus:
                    return read1, read2, True
            return None, None, False

        with pysam.AlignmentFile(self.sam_file1) as sam1:
            with pysam.AlignmentFile(self.sam_file2) as sam2:
                normal_pairs = 0
                chimeric_pairs = 0
                abnormal_pairs = 0

                sam1_iter = iter(sam1)
                sam2_iter = iter(sam2)
                try:
                    qname1, reads1, next_read1 = _all_reads(sam1_iter)
                    qname2, reads2, next_read2 = _all_reads(sam2_iter)
                    while True:
                        check1 = False
                        check2 = False
                        previous_qname1 = qname1
                        previous_qname2 = qname2

                        cmp = natural_cmp(qname1, qname2)
                        if cmp == 0:  # read name identical
                            read1, read2, is_chimeric = _find_pair(reads1, reads2)
                            if read1 is not None and read2 is not None:
                                yield (read1, read2)
                                if is_chimeric:
                                    chimeric_pairs += 1
                                else:
                                    normal_pairs += 1
                            else:
                                abnormal_pairs += 1
                            qname1, reads1, next_read1 = _all_reads(sam1_iter, last_read=next_read1)
                            qname2, reads2, next_read2 = _all_reads(sam2_iter, last_read=next_read2)
                            check1, check2 = True, True
                        elif cmp < 0:  # first pointer behind
                            qname1, reads1, next_read1 = _all_reads(sam1_iter, last_read=next_read1)
                            check1 = True
                        else:  # second pointer behind
                            qname2, reads2, next_read2 = _all_reads(sam2_iter, last_read=next_read2)
                            check2 = True

                        # check that the files are sorted
                        if self._check_sorted:
                            if check1 and natural_cmp(previous_qname1, qname1) > 0:
                                raise ValueError("First SAM file is not sorted by "
                                                 "read name (samtools sort -n)! Read names:"
                                                 "{} and {}".format(previous_qname1, qname1))
                            if check2 and natural_cmp(previous_qname2, qname2) > 0:
                                raise ValueError("Second SAM file is not sorted by "
                                                 "read name (samtools sort -n)! Read names:"
                                                 "{} and {}".format(previous_qname2, qname2))
                except StopIteration:
                    logger.info("Done generating read pairs.")
                    logger.info("Normal pairs: {}".format(normal_pairs))
                    logger.info("Chimeric pairs: {}".format(chimeric_pairs))
                    logger.info("Abnormal pairs: {}".format(abnormal_pairs))

    def filter_quality(self, cutoff=30):
        """
        Convenience function that registers a QualityFilter.
        The actual algorithm and rationale used for filtering will depend on the
        internal _mapper attribute.

        :param cutoff: Minimum mapping quality (mapq) a read must have to pass
                       the filter
        """
        mask = Mask('mapq', 'Mask read pairs with a mapping quality lower than {}'.format(cutoff), ix=len(self.filters))
        quality_filter = QualityFilter(cutoff, mask)
        self.add_filter(quality_filter)

    def filter_unmapped(self):
        """
        Convenience function that registers an UnmappedFilter.
        """
        mask = Mask('unmapped', 'Mask read pairs that are unmapped', ix=len(self.filters))
        unmapped_filter = UnmappedFilter(mask)
        self.add_filter(unmapped_filter)

    def filter_non_unique(self, strict=True):
        """
        Convenience function that registers a UniquenessFilter.
        The actual algorithm and rationale used for filtering will depend on the
        internal _mapper attribute.

        :param strict: If True will filter if XS tag is present. If False,
                       will filter only when XS tag is not 0. This is applied if
                       alignments are from bowtie2.
        """
        mask = Mask('uniqueness', 'Mask reads that do not map uniquely (according to XS tag)', ix=len(self.filters))
        uniqueness_filter = UniquenessFilter(strict, mask)
        self.add_filter(uniqueness_filter)


class FragmentRead(object):
    """
    Class representing a fragment-mapped read.

    .. attribute:: fragment

        A :class:`~kaic.data.genomic.GenomicRegion` delineated by
        restriction sites.

    .. attribute:: position

        The position of this read in base-pairs (1-based) from the
        start of the chromosome it maps to.

    .. attribute:: strand

        The strand this read maps to (-1 or +1).
    """

    def __init__(self, fragment=None, position=None, strand=0, qname_ix=None):
        """
        Initialize this :class:`~FragmentRead` object.

        :param fragment: A :class:`~kaic.data.genomic.GenomicRegion` delineated by
                         restriction sites.
        :param position: The position of this read in base-pairs (1-based) from the
                         start of the chromosome it maps to.
        :param strand: The strand this read maps to (-1 or +1).
        """
        self.fragment = fragment
        self.position = position
        self.strand = strand
        self.qname_ix = qname_ix

    def re_distance(self):
        return min(abs(self.position - self.fragment.start),
                   abs(self.position - self.fragment.end))

    def __repr__(self):
        return "%s: %d-(%d[%d])-%d" % (self.fragment.chromosome,
                                       self.fragment.start,
                                       self.position,
                                       self.strand,
                                       self.fragment.end)


class LazyFragmentRead(FragmentRead):
    def __init__(self, row, pairs, side="left"):
        self.row = row
        self.pairs = pairs
        self.side = side

    @property
    def position(self):
        return self.row[self.side + "_read_position"]

    @property
    def strand(self):
        return self.row[self.side + "_read_strand"]

    @property
    def qname_ix(self):
        return self.row[self.side + "_read_qname_ix"]

    @property
    def fragment(self):
        return LazyFragment(self.row, self.pairs, side=self.side)


class LazyFragment(GenomicRegion):
    def __init__(self, row, pairs, ix=None, side="left"):
        self.row = row
        self.pairs = pairs
        self.side = side
        self.static_ix = ix

    @property
    def chromosome(self):
        return self.pairs._ix_to_chromosome[self.row[self.side + "_fragment_chromosome"]]

    @property
    def start(self):
        return self.row[self.side + "_fragment_start"]

    @property
    def end(self):
        return self.row[self.side + "_fragment_end"]

    @property
    def strand(self):
        return 1

    @property
    def ix(self):
        if self.static_ix is None:
            return self.row['source'] if self.side == 'left' else self.row['sink']
        return self.static_ix


class ReadPairs(RegionPairsTable):
    _classid = 'READPAIRS'

    def __init__(self, file_name=None, mode='a',
                 _group_name='fragment_map',
                 _table_name_fragments='fragments',
                 _table_name_pairs='pairs',
                 tmpdir=None):
        """
        Initialize empty FragmentMappedReadPairs object.

        :param file_name: Path to a file that will be created to save
                          this object or path to an existing HDF5 file
                          representing a FragmentMappedReadPairs object.
        :param group_name: Internal, name for hdf5 group that info for
                           this object will be saved under
        :param table_name_fragments: Internal, name of the HDF5 node
                                     that will house the region/fragment
                                     data
        """
        RegionPairsTable.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                  additional_edge_fields={
                                      'ix': t.Int32Col(pos=0),
                                      'left_read_position': t.Int64Col(pos=1),
                                      'left_read_strand': t.Int8Col(pos=2),
                                      'left_fragment_start': t.Int64Col(pos=3),
                                      'left_fragment_end': t.Int64Col(pos=4),
                                      'left_fragment_chromosome': t.Int32Col(pos=5),
                                      'right_read_position': t.Int64Col(pos=6),
                                      'right_read_strand': t.Int8Col(pos=7),
                                      'right_fragment_start': t.Int64Col(pos=8),
                                      'right_fragment_end': t.Int64Col(pos=9),
                                      'right_fragment_chromosome': t.Int32Col(pos=10)
                                  },
                                  _table_name_regions=_table_name_fragments,
                                  _table_name_edges=_table_name_pairs)

        self._pairs = self._edges
        self._pair_count = sum(edge_table._original_len()
                               for _, edge_table in self._iter_edge_tables())
        self._ix_to_chromosome = dict()
        self._chromosome_to_ix = dict()
        self._update_references()

    def _update_references(self):
        """
        Update internal chromosome index dictionaries.
        """
        chromosomes = []
        for region in self.regions(lazy=True):
            if len(chromosomes) == 0 or chromosomes[-1] != region.chromosome:
                chromosomes.append(region.chromosome)

        for i, chromosome in enumerate(chromosomes):
            self._ix_to_chromosome[i] = chromosome
            self._chromosome_to_ix[chromosome] = i

    def _flush_regions(self):
        if self._regions_dirty:
            self._regions.flush()
            self._update_references()
            self._regions_dirty = False

    def flush(self, update_mappability=True, silent=config.hide_progressbars):
        RegionPairsTable.flush(self, silent=silent)

    def _read_fragment_info(self, read):
        chromosome = read.reference_name
        fragment_info = None
        for row in self._regions.where(
                        "(start <= %d) & (end >= %d) & (chromosome == b'%s')" % (read.pos, read.pos, chromosome)):
            fragment_info = [row['ix'], self._chromosome_to_ix[chromosome], row['start'], row['end']]

        if fragment_info is None:
            raise ValueError("No matching region can be found for {}".format(read))

        return fragment_info

    def _read_pair_fragment_info(self, read_pair):
        read1, read2 = read_pair
        f_ix1, f_chromosome_ix1, f_start1, f_end1 = self._read_fragment_info(read1)
        f_ix2, f_chromosome_ix2, f_start2, f_end2 = self._read_fragment_info(read2)
        r_strand1 = -1 if read1.flag & 16 else 1
        r_strand2 = -1 if read2.flag & 16 else 1
        return ((read1.pos, r_strand1, f_ix1, f_chromosome_ix1, f_start1, f_end1),
                (read2.pos, r_strand2, f_ix2, f_chromosome_ix2, f_start2, f_end2))

    def _read_pairs_fragment_info(self, read_pairs, threads=4, batch_size=1000000, timeout=180):
        fragment_infos = defaultdict(list)
        fragment_ends = defaultdict(list)
        for region in self.regions(lazy=True):
            chromosome = region.chromosome
            fragment_infos[chromosome].append((region.ix, self._chromosome_to_ix[chromosome],
                                               region.start, region.end))
            fragment_ends[chromosome].append(region.end)

        worker_pool = None
        t_pairs = None
        try:
            monitor = Monitor()
            input_queue = mp.Queue(maxsize=2*threads)
            output_queue = mp.Queue(maxsize=2*threads)

            monitor.set_generating_pairs(True)
            t_pairs = threading.Thread(target=_read_pairs_worker, args=(read_pairs, input_queue,
                                                                        monitor, batch_size))
            t_pairs.daemon = True
            t_pairs.start()

            worker_pool = mp.Pool(threads, _fragment_info_worker,
                                  (monitor, input_queue, output_queue, fragment_infos, fragment_ends))

            output_counter = 0
            while output_counter < monitor.value() or not monitor.workers_idle() or monitor.is_generating_pairs():
                try:
                    read_pair_infos = output_queue.get(block=True, timeout=timeout)

                    for read1_info, read2_info in msgpack.loads(read_pair_infos):
                        yield read1_info, read2_info
                    output_counter += 1
                    del read_pair_infos
                except Empty:
                    logger.debug("Reached SAM pair generator timeout. This could mean that no "
                                 "valid read pairs were found after filtering. "
                                 "Check filter settings!")
        finally:
            if worker_pool is not None:
                worker_pool.terminate()
            if t_pairs is not None:
                t_pairs.join()

    def _add_infos(self, fi1, fi2):
        r_pos1, r_strand1, f_ix1, f_chromosome_ix1, f_start1, f_end1 = fi1
        r_pos2, r_strand2, f_ix2, f_chromosome_ix2, f_start2, f_end2 = fi2

        edge = Edge(ix=self._pair_count,
                    source=f_ix1, sink=f_ix2,
                    left_read_position=r_pos1, right_read_position=r_pos2,
                    left_read_strand=r_strand1, right_read_strand=r_strand2,
                    left_fragment_start=f_start1, right_fragment_start=f_start2,
                    left_fragment_end=f_end1, right_fragment_end=f_end2,
                    left_fragment_chromosome=f_chromosome_ix1,
                    right_fragment_chromosome=f_chromosome_ix2)

        self._add_pair(edge)

    def _fast_add_infos(self, fi1, fi2, default_edge):
        if fi1[2] > fi2[2]:
            r_pos1, r_strand1, f_ix1, f_chromosome_ix1, f_start1, f_end1 = fi2
            r_pos2, r_strand2, f_ix2, f_chromosome_ix2, f_start2, f_end2 = fi1
        else:
            r_pos1, r_strand1, f_ix1, f_chromosome_ix1, f_start1, f_end1 = fi1
            r_pos2, r_strand2, f_ix2, f_chromosome_ix2, f_start2, f_end2 = fi2

        edge = copy.copy(default_edge)
        edge[self._field_names_dict['ix']] = self._pair_count
        edge[self._field_names_dict['source']] = f_ix1
        edge[self._field_names_dict['sink']] = f_ix2
        edge[self._field_names_dict['left_read_position']] = r_pos1
        edge[self._field_names_dict['right_read_position']] = r_pos2
        edge[self._field_names_dict['left_read_strand']] = r_strand1
        edge[self._field_names_dict['right_read_strand']] = r_strand2
        edge[self._field_names_dict['left_fragment_start']] = f_start1
        edge[self._field_names_dict['right_fragment_start']] = f_start2
        edge[self._field_names_dict['left_fragment_end']] = f_end1
        edge[self._field_names_dict['right_fragment_end']] = f_end2
        edge[self._field_names_dict['left_fragment_chromosome']] = f_chromosome_ix1
        edge[self._field_names_dict['right_fragment_chromosome']] = f_chromosome_ix2

        self._add_edge_from_tuple(tuple(edge))

        self._pair_count += 1

    def _default_edge_list(self):
        record = [None] * len(self._field_names_dict)
        for name, ix in self._field_names_dict.items():
            record[ix] = self._edge_field_defaults[name]
        return record

    def add_read_pair(self, read_pair, flush=True):
        fi1, fi2 = self._read_pair_fragment_info(read_pair)
        self._add_infos(fi1, fi2)
        if flush:
            self.flush()

    def _flush_fragment_info_buffer(self):
        for (source_partition, sink_partition), edges in self._edge_buffer.items():
            edge_table = self._edge_table(source_partition, sink_partition)
            row = edge_table.row

            for fi1, fi2 in edges:
                if fi1[2] > fi2[2]:
                    fi1, fi2 = fi2, fi1

                row['ix'] = self._pair_count
                row['source'] = fi1[2]
                row['sink'] = fi2[2]
                row['left_read_position'] = fi1[0]
                row['right_read_position'] = fi2[0]
                row['left_read_strand'] = fi1[1]
                row['right_read_strand'] = fi2[1]
                row['left_fragment_start'] = fi1[4]
                row['right_fragment_start'] = fi2[4]
                row['left_fragment_end'] = fi1[5]
                row['right_fragment_end'] = fi2[5]
                row['left_fragment_chromosome'] = fi1[3]
                row['right_fragment_chromosome'] = fi2[3]
                row.append()
                self._pair_count += 1

            edge_table.flush()
        self._edge_buffer = defaultdict(list)

    def add_read_pairs(self, read_pairs, batch_size=1000000, threads=1):
        self._edges_dirty = True
        self._disable_edge_indexes()

        start_time = timer()
        chunk_start_time = timer()
        pairs_counter = 0
        for fi1, fi2 in self._read_pairs_fragment_info(read_pairs, batch_size=batch_size, threads=threads):
            source_partition, sink_partition = self._get_edge_table_tuple(fi1[2], fi2[2])
            self._edge_buffer[(source_partition, sink_partition)].append([fi1, fi2])

            pairs_counter += 1
            if pairs_counter % self._edge_buffer_size == 0:
                self._flush_fragment_info_buffer()
                end_time = timer()
                logger.debug("Wrote {} pairs in {}s (current 1M chunk: {}s)".format(
                    pairs_counter, end_time - start_time, end_time - chunk_start_time
                ))
                chunk_start_time = timer()
        self._flush_fragment_info_buffer()
        end_time = timer()
        logger.debug("Wrote {} pairs in {}s".format(
            pairs_counter, end_time - start_time
        ))

        logger.info('Done saving read pairs.')

        if isinstance(read_pairs, ReadPairGenerator):
            stats = read_pairs.stats()
            if 'read_filter_stats' not in self.meta:
                self.meta.read_filter_stats = stats
            else:
                self.meta.read_filter_stats = add_dict(self.meta.read_filter_stats, stats)

        self.flush()
        logger.info("Done adding pairs.")

    def _add_pair(self, pair):
        self.add_edge(pair, check_nodes_exist=False, replace=True)
        self._pair_count += 1

    def _pair_from_row(self, row, lazy=False):
        """
        Convert a PyTables row to a FragmentReadPair
        """
        if lazy:
            left_read = LazyFragmentRead(row, self, side="left")
            right_read = LazyFragmentRead(row, self, side="right")
        else:
            fragment1 = GenomicRegion(start=row['left_fragment_start'],
                                      end=row['left_fragment_end'],
                                      chromosome=self._ix_to_chromosome[row['left_fragment_chromosome']],
                                      ix=row['source'])
            fragment2 = GenomicRegion(start=row['right_fragment_start'],
                                      end=row['right_fragment_end'],
                                      chromosome=self._ix_to_chromosome[row['right_fragment_chromosome']],
                                      ix=row['sink'])

            left_read = FragmentRead(fragment1, position=row['left_read_position'],
                                     strand=row['left_read_strand'])
            right_read = FragmentRead(fragment2, position=row['right_read_position'],
                                      strand=row['right_read_strand'])

        return FragmentReadPair(left_read=left_read, right_read=right_read, ix=row['ix'])

    def get_ligation_structure_biases(self, sampling=None, skip_self_ligations=True):

        """
        Compute the ligation biases (inward and outward to same-strand) of this data set.

        :param sampling: Approximate number of data points to average per point
                         in the plot. If None (default), this will
                         be determined on a best-guess basis.
        :param skip_self_ligations: If True (default), will not consider
                                    self-ligated fragments for assessing
                                    the error rates.
        """
        l = len(self)
        type_same = 0
        type_inward = 1
        type_outward = 2

        def _init_gaps_and_types():
            same_count = 0
            inward_count = 0
            outward_count = 0
            same_fragment_count = 0
            inter_chrm_count = 0
            gaps = []
            types = []

            with RareUpdateProgressBar(max_value=len(self), silent=config.hide_progressbars) as pb:
                for i, pair in enumerate(self):
                    if pair.is_same_fragment():
                        same_fragment_count += 1
                        if skip_self_ligations:
                            continue
                    if pair.is_same_chromosome():
                        gap_size = pair.get_gap_size()
                        if gap_size > 0:
                            gaps.append(gap_size)
                            if pair.is_outward_pair():
                                types.append(type_outward)
                                outward_count += 1
                            elif pair.is_inward_pair():
                                types.append(type_inward)
                                inward_count += 1
                            else:
                                types.append(0)
                                same_count += 1
                    else:
                        inter_chrm_count += 1
                    pb.update(i)

            logger.info("Pairs: %d" % l)
            logger.info("Inter-chromosomal: {}".format(inter_chrm_count))
            logger.info("Same fragment: {}".format(same_fragment_count))
            logger.info("Same: {}".format(same_count))
            logger.info("Inward: {}".format(inward_count))
            logger.info("Outward: {}".format(outward_count))
            return gaps, types

        def _sort_data(gaps, types):
            points = zip(gaps, types)
            sorted_points = sorted(points)
            return zip(*sorted_points)

        def _guess_sampling(sampling):
            if sampling is None:
                sampling = max(100, int(l * 0.0025))
            logger.info("Number of data points averaged per point in plot: {}".format(sampling))
            return sampling

        def _calculate_ratios(gaps, types, sampling):
            x = []
            inward_ratios = []
            outward_ratios = []
            bin_sizes = []
            counter = 0
            same_counter = 0
            mids = 0
            outwards = 0
            inwards = 0
            for typ, gap in zip(types, gaps):
                mids += gap
                if typ == type_same:
                    same_counter += 1
                elif typ == type_inward:
                    inwards += 1
                else:
                    outwards += 1
                counter += 1
                if same_counter > sampling:
                    x.append(int(mids / counter))
                    inward_ratios.append(inwards / same_counter)
                    outward_ratios.append(outwards / same_counter)
                    bin_sizes.append(counter)
                    same_counter = 0
                    counter = 0
                    mids = 0
                    outwards = 0
                    inwards = 0
            return list(map(np.array, [x, inward_ratios, outward_ratios, bin_sizes]))

        gaps, types = _init_gaps_and_types()
        # sort data
        gaps, types = _sort_data(gaps, types)
        # best guess for number of data points
        sampling = _guess_sampling(sampling)
        # calculate ratios
        return _calculate_ratios(gaps, types, sampling)

    @staticmethod
    def _auto_dist(dists, ratios, sample_sizes, p=0.05, expected_ratio=0.5):
        """
        Function that attempts to infer sane distances for filtering inward
        and outward read pairs

        :param dists: List of distances in bp.
        :param ratios: List of ratios
        """

        def x_prop(p_obs, p_exp, n):
            obs = p_obs * n
            exp = p_exp * n
            p = (obs + exp) / (n * 2)
            return abs((p_exp - p_obs) / np.sqrt(p * (1 - p) * (2 / n)))

        ratios = np.clip(ratios, 0.0, 1.0)
        z_scores = np.array([x_prop(r, expected_ratio, b) for r, b in zip(ratios, sample_sizes)])
        which_valid = z_scores < 1.96
        which_valid_indices = np.argwhere(which_valid).flatten()
        if len(which_valid_indices) > 0:
            return int(dists[which_valid_indices[0]])
        return None

    def filter(self, pair_filter, queue=False, log_progress=not config.hide_progressbars):
        pair_filter.set_pairs_object(self)

        total = 0
        filtered = 0
        if not queue:
            with RareUpdateProgressBar(max_value=sum(1 for _ in self._edges),
                                       silent=not log_progress) as pb:
                for i, (_, edge_table) in enumerate(self._iter_edge_tables()):
                    stats = edge_table.filter(pair_filter, _logging=False)
                    for key, value in stats.items():
                        if key != 0:
                            filtered += stats[key]
                        total += stats[key]
                    pb.update(i)
            if log_progress:
                logger.info("Total: {}. Filtered: {}".format(total, filtered))
        else:
            for _, edge_table in self._iter_edge_tables():
                edge_table.queue_filter(pair_filter)

    def run_queued_filters(self, log_progress=not config.hide_progressbars):
        """
        Run queued filters.

        :param log_progress: If true, process iterating through all edges
                             will be continuously reported.
        """
        total = 0
        filtered = 0
        with RareUpdateProgressBar(max_value=sum(1 for _ in self._edges),
                                   silent=not log_progress) as pb:
            for i, (_, edge_table) in enumerate(self._iter_edge_tables()):
                stats = edge_table.run_queued_filters(_logging=False)
                for key, value in stats.items():
                    if key != 0:
                        filtered += stats[key]
                    total += stats[key]
                pb.update(i)
        if log_progress:
            logger.info("Total: {}. Filtered: {}".format(total, filtered))

    def filter_pcr_duplicates(self, threshold=3, queue=False):
        """
        Convenience function that applies an :class:`~PCRDuplicateFilter`.

        :param threshold: If distance between two alignments is smaller or equal the threshold, the alignments
                          are considered to be starting at the same position
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        """
        mask = self.add_mask_description('pcr_duplicate', 'Mask read pairs that are considered PCR duplicates')
        pcr_duplicate_filter = AOPCRDuplicateFilter(pairs=self, threshold=threshold, mask=mask)
        self.filter(pcr_duplicate_filter, queue)

    def filter_inward(self, minimum_distance=None, queue=False, *args, **kwargs):
        """
        Convenience function that applies an :class:`~InwardPairsFilter`.

        :param minimum_distance: Minimum distance inward-facing read
                                 pairs must have to pass this filter
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        :param *args **kwargs: Additional arguments to pass
                               to :met:`~FragmentMappedReadPairs.get_ligation_structure_biases`
        """
        if minimum_distance is None:
            dists, inward_ratios, _, bins_sizes = self.get_ligation_structure_biases(*args, **kwargs)
            minimum_distance = self._auto_dist(dists, inward_ratios, bins_sizes)
        if minimum_distance:
            mask = self.add_mask_description('inward',
                                             'Mask read pairs that are inward facing and < {}bp apart'
                                             .format(minimum_distance))
            logger.info("Filtering out inward facing read pairs < {} bp apart".format(minimum_distance))
            inward_filter = InwardPairsFilter(minimum_distance=minimum_distance, mask=mask)
            self.filter(inward_filter, queue)
        else:
            raise Exception('Could not automatically detect a sane distance threshold for filtering inward reads')

    def filter_outward(self, minimum_distance=None, queue=False, *args, **kwargs):
        """
        Convenience function that applies an :class:`~OutwardPairsFilter`.

        :param minimum_distance: Minimum distance outward-facing read
                                 pairs must have to pass this filter
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        :param *args **kwargs: Additional arguments to pass
                               to :met:`~FragmentMappedReadPairs.get_ligation_structure_biases`
        """
        if minimum_distance is None:
            dists, _, outward_ratios, bins_sizes = self.get_ligation_structure_biases(*args, **kwargs)
            minimum_distance = self._auto_dist(dists, outward_ratios, bins_sizes)
        if minimum_distance:
            mask = self.add_mask_description('outward',
                                             'Mask read pairs that are outward facing and < {}bp apart'
                                             .format(minimum_distance))
            logger.info("Filtering out outward facing read pairs < {} bp apart".format(minimum_distance))
            outward_filter = OutwardPairsFilter(minimum_distance=minimum_distance, mask=mask)
            self.filter(outward_filter, queue)
        else:
            raise Exception('Could not automatically detect a sane distance threshold for filtering outward reads')

    def filter_ligation_products(self, inward_threshold=None, outward_threshold=None, queue=False, *args, **kwargs):
        """
        Convenience function that applies an :class:`~OutwardPairsFilter` and an :class:`~InwardPairsFilter`.

        :param inward_threshold: Minimum distance inward-facing read
                                 pairs must have to pass this filter. If None, will be infered from the data
        :param outward_threshold: Minimum distance outward-facing read
                                 pairs must have to pass this filter. If None, will be infered from the data
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        :param *args **kwargs: Additional arguments to pass
                               to :met:`~FragmentMappedReadPairs.get_ligation_structure_biases`
        """
        self.filter_inward(inward_threshold, queue=queue, *args, **kwargs)
        self.filter_outward(outward_threshold, queue=queue, *args, **kwargs)

    def filter_re_dist(self, maximum_distance, queue=False):
        """
        Convenience function that applies an :class:`~ReDistanceFilter`.

        :param maximum_distance: Maximum distance a read can have to the
                                 nearest region border (=restriction site)
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        """
        mask = self.add_mask_description('re-dist',
                                         'Mask read pairs where a read is >%dbp away from nearest RE site' % maximum_distance)
        re_filter = ReDistanceFilter(maximum_distance=maximum_distance, mask=mask)
        self.filter(re_filter, queue)

    def filter_self_ligated(self, queue=False):
        """
        Convenience function that applies an :class:`~SelfLigationFilter`.

        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        """
        mask = self.add_mask_description('self_ligated',
                                         'Mask read pairs the represet a self-ligated fragment')
        self_ligation_filter = SelfLigationFilter(mask=mask)
        self.filter(self_ligation_filter, queue)

    def __iter__(self):
        """
        Iterate over unfiltered fragment-mapped read pairs.
        """
        return self.pairs(lazy=False)

    def _edge_row_iter(self, intra_chromosomal=True, inter_chromosomal=True, excluded_filters=()):
        """
        Yield rows in edge tables, ordered by partition.
        """
        excluded_masks = self.get_binary_mask_from_masks(excluded_filters)

        for (ix1, ix2), edge_table in self._iter_edge_tables():
            if not intra_chromosomal and ix1 == ix2:
                continue
            if not inter_chromosomal and ix1 != ix2:
                continue

            for row in edge_table.iterrows(excluded_masks=excluded_masks):
                yield row

    def pairs(self, lazy=False, excluded_filters=()):
        for row in self._edge_row_iter(excluded_filters=excluded_filters):
            yield self._pair_from_row(row, lazy=lazy)

    def get_edge(self, item, *row_conversion_args, **row_conversion_kwargs):
        """
        Get an edge by index.

        :param row_conversion_args: Arguments passed to :func:`RegionPairs._row_to_edge`
        :param row_conversion_args: Keyword arguments passed to :func:`RegionPairs._row_to_edge`
        :return: :class:`~Edge`
        """
        if item < 0:
            item += len(self)

        l = 0
        for _, edge_table in self._iter_edge_tables():
            if l <= item < l + len(edge_table):
                res = edge_table[item - l]
                return self._row_to_edge(res, *row_conversion_args, **row_conversion_kwargs)
            l += len(edge_table)
        raise IndexError("index out of range (%d)" % item)

    def __getitem__(self, item):
        if isinstance(item, int):
            edge = self.get_edge(item)
            fragment1 = GenomicRegion(start=edge.left_fragment_start,
                                      end=edge.left_fragment_end,
                                      chromosome=self._ix_to_chromosome[edge.left_fragment_chromosome],
                                      ix=edge.source)
            fragment2 = GenomicRegion(start=edge.right_fragment_start,
                                      end=edge.right_fragment_end,
                                      chromosome=self._ix_to_chromosome[edge.right_fragment_chromosome],
                                      ix=edge.sink)

            left_read = FragmentRead(fragment1, position=edge.left_read_position,
                                     strand=edge.left_read_strand)
            right_read = FragmentRead(fragment2, position=edge.right_read_position,
                                      strand=edge.right_read_strand)

            return FragmentReadPair(left_read=left_read, right_read=right_read, ix=edge.ix)
        else:
            pairs = []
            for row in self.edges.get_row_range(item):
                pairs.append(self._pair_from_row(row, lazy=False))
            return pairs

    def __len__(self):
        l = 0
        for _, edge_table in self._iter_edge_tables():
            l += len(edge_table)
        return l

    def to_hic(self, file_name=None, tmpdir=None, _hic_class=Hic):
        hic = _hic_class(file_name=file_name, mode='w', tmpdir=tmpdir)
        hic.add_regions(self.regions(), preserve_attributes=False)

        hic._disable_edge_indexes()

        l = len(self)
        pairs_counter = 0
        with RareUpdateProgressBar(max_value=l, silent=config.hide_progressbars) as pb:
            for _, pairs_edge_table in self._iter_edge_tables():

                partition_edge_buffer = defaultdict(dict)
                for row in pairs_edge_table:
                    key = (row['source'], row['sink'])
                    source_partition = self._get_partition_ix(key[0])
                    sink_partition = self._get_partition_ix(key[1])
                    if key not in partition_edge_buffer[(source_partition, sink_partition)]:
                        partition_edge_buffer[(source_partition, sink_partition)][key] = 0
                    partition_edge_buffer[(source_partition, sink_partition)][key] += 1
                    pb.update(pairs_counter)
                    pairs_counter += 1

                for hic_partition_key, edge_buffer in viewitems(partition_edge_buffer):
                    hic_edge_table = hic._edge_table(hic_partition_key[0], hic_partition_key[1])
                    row = hic_edge_table.row

                    for (source, sink), weight in viewitems(edge_buffer):
                        row['source'] = source
                        row['sink'] = sink
                        row[hic._default_score_field] = float(weight)
                        row.append()
                    hic_edge_table.flush(update_index=False)
        hic.flush()

        hic._enable_edge_indexes()

        return hic

    def pairs_by_chromosomes(self, chromosome1, chromosome2, lazy=False):
        chromosome_bins = self.chromosome_bins
        if chromosome1 not in chromosome_bins or chromosome2 not in chromosome_bins:
            raise ValueError("Chromosomes {}/{} not in object".format(chromosome1, chromosome2))
        source_partition = self._get_partition_ix(chromosome_bins[chromosome1][0])
        sink_partition = self._get_partition_ix(chromosome_bins[chromosome2][0])
        if source_partition > sink_partition:
            source_partition, sink_partition = sink_partition, source_partition
        for row in self._edge_table_dict[(source_partition, sink_partition)]:
            yield self._pair_from_row(row, lazy=lazy)

    def filter_statistics(self):
        try:
            read_stats = self.meta.read_filter_stats
        except AttributeError:
            read_stats = dict()

        pair_stats = self.mask_statistics(self._pairs)
        if 'unmasked' in pair_stats:
            read_stats['unmasked'] = pair_stats['unmasked']
        pair_stats.update(read_stats)
        return pair_stats


class FragmentReadPair(object):
    """
    Container for two paired :class:`~FragmentRead` objects.
    """

    def __init__(self, left_read, right_read, ix=None):
        self.left = left_read
        self.right = right_read
        self.ix = ix

    def is_same_chromosome(self):
        """
        Check if both reads are mapped to the same chromosome.

        :return: True is reads map to the same chromosome, False
                 otherwise
        """
        return self.left.fragment.chromosome == self.right.fragment.chromosome

    def is_inward_pair(self):
        """
        Check if reads form an "inward-facing" pair.

        A pair is considered inward-facing if the left read maps
        to the plus the right read to the minus strand and both
        reads map to the same chromosome.

        :return: True is reads are inward-facing, False otherwise
        """
        if not self.is_same_chromosome():
            return False

        if self.left.strand == 1 and self.right.strand == -1:
            return True
        return False

    def is_outward_pair(self):
        """
        Check if reads form an "outward-facing" pair.

        A pair is considered outward-facing if the left read maps
        to the minus the right read to the plus strand and both
        reads map to the same chromosome.

        :return: True is reads are outward-facing, False otherwise
        """
        if not self.is_same_chromosome():
            return False

        if self.left.strand == -1 and self.right.strand == 1:
            return True
        return False

    def is_same_pair(self):
        """
        Check if reads face in the same direction.

        :return: True if reads map to the same fragment,
                 False otherwise.
        """
        if not self.is_same_chromosome():
            return False

        if self.left.strand == self.right.strand:
            return True
        return False

    def is_same_fragment(self):
        """
        Check if reads map to the same fragment.

        :return: True if reads map to the same fragment,
                 False otherwise.
        """
        if not self.is_same_chromosome():
            return False

        return self.left.fragment.start == self.right.fragment.start

    def get_gap_size(self):
        """
        Get the gap size in base pairs between the fragments these
        reads map to.

        :return: 0 if reads map to the same fragment or neighboring
                 fragments, the distance between fragments if they
                 are on the same chromosome, None otherwise
        """
        if not self.is_same_chromosome():
            return None

        if self.is_same_fragment():
            return 0

        gap = self.right.fragment.start - self.left.fragment.end

        if gap == 1:  # neighboring fragments
            return 0

        return gap

    def __getitem__(self, key):
        if key == 0:
            return self.left
        if key == 1:
            return self.right
        raise KeyError("Can only access read [0] and read [1]")

    def __repr__(self):
        left_repr = self.left.__repr__()
        right_repr = self.right.__repr__()
        return "{} -- {}".format(left_repr, right_repr)


class ReadFilter(MaskFilter):
    """
    Abstract class that provides filtering functionality for the
    Reads object.

    Extends MaskFilter and overrides valid(self, read) to make
    Read filtering more "natural".

    To create custom filters for the Reads object, extend this
    class and override the valid_read(self, read) method.
    valid_read should return False for a specific Read object
    if the object is supposed to be filtered/masked and True
    otherwise. See :class:`~QualityFilter` for an example.

    Pass a custom filter to the filter method in :class:`~Reads`
    to aplly it.
    """

    def __init__(self, mask=None):
        """
        Initialize ReadFilter.

        :param mask: The Mask object that should be used to mask
                     filtered Read objects. If None the default
                     Mask will be used.
        """
        super(ReadFilter, self).__init__(mask)
        self._reads = None

    def valid_read(self, read):
        """
        Determine if a Read object is valid or should be filtered.

        When implementing custom ReadFilters this method must be
        overridden. It should return False for Read objects that
        are to be fitered and True otherwise.

        Internally, the Reads object will iterate over all Read
        instances to determine their validity on an individual
        basis.

        :param read: A :class:`~Read` object
        :return: True if Read is valid, False otherwise
        """
        raise NotImplementedError("ReadFilters must implement valid_read function")

    def set_reads_object(self, reads_object):
        """
        Set the Reads instance to be filtered by this ReadFilter.

        Used internally by Reads instance.

        :param reads_object: Reads object
        """
        self._reads = reads_object

    def valid(self, row):
        """
        Map valid_read to MaskFilter.valid(self, row).

        :param row: A pytables Table row.
        :return: The boolean value returned by valid_read.
        """
        read = self._reads._row2read(row, lazy=True)
        return self.valid_read(read)


class QualityFilter(ReadFilter):
    """
    Filter mapped reads based on mapping quality.
    """

    def __init__(self, cutoff=30, mask=None):
        """
        :param cutoff: Lowest mapping quality that is still
                       considered acceptable for a mapped read.
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered reads.
        """
        super(QualityFilter, self).__init__(mask)
        self.cutoff = cutoff

    def valid_read(self, read):
        """
        Check if a read has a mapq >= cutoff.
        """
        return read.mapq >= self.cutoff


class ContaminantFilter(ReadFilter):
    """
    Filter mapped reads based on mapping quality.
    """

    def __init__(self, contaminant_reads, mask=None):
        """
        :param contaminant_reads: A SAM/BAM file representing a
                                  contaminant
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered reads.
        """
        super(ContaminantFilter, self).__init__(mask)

        self.contaminant_names = set()
        with pysam.AlignmentFile(contaminant_reads) as contaminant:
            for read in contaminant:
                self.contaminant_names.add(read.qname)

    def valid_read(self, read):
        """
        Check if a read also maps to a contaminant
        """
        if read.qname in self.contaminant_names:
            return False
        return True


class BwaMemQualityFilter(ReadFilter):
    """
    Filters `bwa mem` generated alignements base on the alignment score
    (normalized by the length of the alignment).
    """
    def __init__(self, cutoff=0.90, mask=None):
        """
        :param cutoff: Ratio of the alignment score to the maximum score
                       possible for an alignment that long
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered reads.
        """
        super(BwaMemQualityFilter, self).__init__(mask)
        self.cutoff = cutoff

    def valid_read(self, read):
        """
        Check if a read has a high alignment score.
        """
        if read.alen:
            return float(read.get_tag('AS')) / read.alen >= self.cutoff
        return False


class UniquenessFilter(ReadFilter):
    """
    Filter reads that do not map uniquely to the reference sequence.
    """

    def __init__(self, strict=True, mask=None):
        """
        :param strict: If True, valid_read checks only for the
                       presence of an XS tag. If False, the value
                       of an XS tag also has to be different from 0.
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered reads.
        """
        self.strict = strict
        super(UniquenessFilter, self).__init__(mask)

    def valid_read(self, read):
        """
        Check if a read has an XS tag.

        If strict is disabled checks if a read has an XS tag and
        the value of the XS tag id different from 0.
        """
        tags = read.tags
        if self.strict:
            for tag in tags:
                if tag[0] == 'XS':
                    return False
        else:
            tags = {tag[0]: tag[1] for tag in tags}

            if tags['AS'] == tags['XS']:
                return False
        return True


class BwaMemUniquenessFilter(ReadFilter):
    """
    Filters `bwa mem` generated alignements based on whether they are unique or not.
    The presence of a non-zero XS tag does not mean a read is a multi-mapping one.
    Instead, we make sure that the ratio XS/AS is inferior to a certain threshold.
    """
    def __init__(self, strict=False, mask=None):
        """
        :param strict: If True, valid_read checks only for the
                       presence of an XA tag. If False, the edit
                       distance (NM) of an alternative alignment has to be
                       the same or better as the original NM.
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered reads.
        """
        super(BwaMemUniquenessFilter, self).__init__(mask)
        if strict:
            self.valid_read = self._valid_read_strict
        else:
            self.valid_read = self._valid_read

    def _valid_read(self, read):
        try:
            xa = read.get_tag('XA')
        except KeyError:
            return True

        try:
            nm = read.get_tag('NM')
        except KeyError:
            return False

        for alt in xa.split(';'):
            if alt == '':
                continue
            _, _, _, nm_alt = alt.split(',')
            if int(nm_alt) <= nm:
                return False
        return True

    def _valid_read_strict(self, read):
        if read.has_tag('XA'):
            return False
        return True

    def valid_read(self, read):
        pass


class UnmappedFilter(ReadFilter):
    """
    Filter reads that do not map to the reference sequence.
    """
    def __init__(self, mask=None):
        """
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered reads.
        """
        super(UnmappedFilter, self).__init__(mask)

    def valid_read(self, read):
        """
        Check if the the flag bit 4 is set.
        """
        if read.flag & 4:
            return False
        return True


class FragmentMappedReadPairFilter(with_metaclass(ABCMeta, MaskFilter)):
    """
    Abstract class that provides filtering functionality for the
    :class:`~FragmentReadPair` object.

    Extends MaskFilter and overrides valid(self, read) to make
    :class:`~FragmentReadPair` filtering more "natural".

    To create custom filters for the :class:`~FragmentMappedReadPairs`
    object, extend this
    class and override the :func:`~FragmentMappedReadPairFilter.valid_pair` method.
    valid_pair should return False for a specific :class:`~FragmentReadPair` object
    if the object is supposed to be filtered/masked and True
    otherwise. See :class:`~InwardPairsFilter` for an example.

    Pass a custom filter to the filter method in :class:`~FragmentMappedReadPairs`
    to apply it.
    """

    def __init__(self, mask=None):
        super(FragmentMappedReadPairFilter, self).__init__(mask)
        self.pairs = None

    def set_pairs_object(self, pairs):
        self.pairs = pairs

    @abstractmethod
    def valid_pair(self, fr_pair):
        pass

    def valid(self, row):
        """
        Map validity check of rows to pairs.
        """
        pair = self.pairs._pair_from_row(row, lazy=True)
        return self.valid_pair(pair)


class InwardPairsFilter(FragmentMappedReadPairFilter):
    """
    Filter inward-facing read pairs at a distance less
    than a specified cutoff.
    """

    def __init__(self, minimum_distance=10000, mask=None):
        """
        Initialize filter with filter settings.

        :param minimum_distance: Minimum distance below which
                                 reads are invalidated
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered reads.
        """
        super(InwardPairsFilter, self).__init__(mask=mask)
        self.minimum_distance = minimum_distance

    def valid_pair(self, pair):
        """
        Check if a pair is inward-facing and <minimum_distance apart.
        """
        if not pair.is_inward_pair():
            return True

        if pair.get_gap_size() > self.minimum_distance:
            return True
        return False


class PCRDuplicateFilter(FragmentMappedReadPairFilter):
    """
    Masks alignments that are suspected to be PCR duplicates.
    In order to be considered duplicates, two pairs need to have identical
    start positions of their respective left alignments AND of their right alignments.
    """

    def __init__(self, pairs, threshold=3, mask=None):
        """
        Initialize filter with filter settings.

        :param pairs: The :class:`~FragmentMappedReadPairs` instance that the filter will be
                      applied to
        :param threshold: If distance between two alignments is smaller or equal the threshold,
                          the alignments are considered to be starting at the same position
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered reads.
        """
        FragmentMappedReadPairFilter.__init__(self, mask=mask)
        self.threshold = threshold
        self.pairs = pairs
        self.duplicates_set = set()
        self.duplicate_stats = defaultdict(int)
        self._mark_duplicates(self.pairs._pairs)

        n_dups = len(self.duplicates_set)
        percent_dups = 1. * n_dups / self.pairs._pairs._original_len()
        logger.info("PCR duplicate stats: " +
                    "{} ({:.1%}) of pairs marked as duplicate. ".format(n_dups, percent_dups) +
                    " (multiplicity:occurances) " +
                    " ".join("{}:{}".format(k, v) for k, v in self.duplicate_stats.items()))

    def _mark_duplicates(self, edge_table):
        # In order for sorted iteration to work, column needs to be indexed
        try:
            edge_table.cols.left_read_position.create_csindex()
            index_existed = False
        except ValueError:  # Index already exists
            index_existed = True

        # Using itersorted from Table class, since MaskedTable.itersorted only yields unmasked entries
        all_iter = super(MaskedTable, edge_table).itersorted(sortby="left_read_position")
        cur_pos = {}
        cur_duplicates = {}
        for p in all_iter:
            pair = self.pairs._pair_from_row(p, lazy=True)
            chrm = (pair.left.fragment.chromosome, pair.right.fragment.chromosome)
            if cur_pos.get(chrm) is None:
                cur_pos[chrm] = (pair.left.position, pair.right.position)
                cur_duplicates[chrm] = 1
                continue
            if (abs(pair.left.position - cur_pos[chrm][0]) <= self.threshold and
                    abs(pair.right.position - cur_pos[chrm][1]) <= self.threshold):
                self.duplicates_set.add(pair.ix)
                cur_duplicates[chrm] += 1
                continue
            if cur_duplicates[chrm] > 1:
                self.duplicate_stats[cur_duplicates[chrm]] += 1
            cur_pos[chrm] = (pair.left.position, pair.right.position)
            cur_duplicates[chrm] = 1

        if not index_existed:
            edge_table.cols.left_read_position.remove_index()

    def valid_pair(self, pair):
        """
        Check if a pair is duplicated.
        """
        if pair.ix in self.duplicates_set:
            return False
        return True


class AOPCRDuplicateFilter(PCRDuplicateFilter):
    """
    Masks alignments that are suspected to be PCR duplicates.
    In order to be considered duplicates, two pairs need to have identical
    start positions of their respective left alignments AND of their right alignments.
    """

    def __init__(self, pairs, threshold=3, mask=None):
        """
        Initialize filter with filter settings.

        :param pairs: The :class:`~FragmentMappedReadPairs` instance that the filter will be
                      applied to
        :param threshold: If distance between two alignments is smaller or equal the threshold,
                          the alignments are considered to be starting at the same position
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered reads.
        """
        FragmentMappedReadPairFilter.__init__(self, mask=mask)
        self.threshold = threshold
        self.pairs = pairs
        self.duplicates_set = set()
        self.duplicate_stats = defaultdict(int)
        original_len = 0
        for _, edge_table in self.pairs._iter_edge_tables():
            original_len += edge_table._original_len()
            self._mark_duplicates(edge_table)

        n_dups = len(self.duplicates_set)
        percent_dups = 1. * n_dups / original_len
        logger.info("PCR duplicate stats: " +
                    "{} ({:.1%}) of pairs marked as duplicate. ".format(n_dups, percent_dups) +
                    " (multiplicity:occurances) " +
                    " ".join("{}:{}".format(k, v) for k, v in self.duplicate_stats.items()))


class OutwardPairsFilter(FragmentMappedReadPairFilter):
    """
    Filter outward-facing read pairs at a distance less
    than a specified cutoff.
    """

    def __init__(self, minimum_distance=10000, mask=None):
        """
        Initialize filter with filter settings.

        :param minimum_distance: Minimum distance below which
                                 outward-facing reads are invalidated
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered reads.
        """
        super(OutwardPairsFilter, self).__init__(mask=mask)
        self.minimum_distance = minimum_distance

    def valid_pair(self, pair):
        if not pair.is_outward_pair():
            return True

        if pair.get_gap_size() > self.minimum_distance:
            return True
        return False


class ReDistanceFilter(FragmentMappedReadPairFilter):
    """
    Filters read pairs where one or both reads are more than
    maximum_distance away from the nearest restriction site.
    """

    def __init__(self, maximum_distance=10000, mask=None):
        super(ReDistanceFilter, self).__init__(mask=mask)
        self.maximum_distance = maximum_distance

    def valid_pair(self, pair):
        """
        Check if any read is >maximum_distance away from RE site.
        """
        d1 = min(abs(pair.left.position - pair.left.fragment.start),
                 abs(pair.left.position - pair.left.fragment.end))
        d2 = min(abs(pair.right.position - pair.right.fragment.start),
                 abs(pair.right.position - pair.right.fragment.end))

        if d1 + d2 > self.maximum_distance:
            return False

        return True


class SelfLigationFilter(FragmentMappedReadPairFilter):
    """
    Filters read pairs where one or both reads are more than
    maximum_distance away from the nearest restriction site.
    """

    def __init__(self, mask=None):
        super(SelfLigationFilter, self).__init__(mask=mask)

    def valid_pair(self, pair):
        """
        Check if any read is >maximum_distance away from RE site.
        """
        if pair.is_same_fragment():
            return False
        return True
