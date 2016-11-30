"""
Module for working with sequencing data.

This module forms the first part of the Hi-C pipeline after externally
mapping FASTQ sequencing reads to a reference genome with a read aligner
such as Bowtie 2.

Facilities for loading and filtering reads from SAM/BAM files are
provided in the :class:`~Reads` class.

To pair reads (e.g. from paired-end sequencing) the
:class:`~FragmentMappedReadPairs` class can be used.

For example, let's assume we have two SAM files of reads mapped to a
reference sequence using Bowtie 2. A sample process of loading, pairing,
and filtering could look like this:

.. code:: python

    # load both sides of reads separately
    reads1 = Reads("/path/to/first_paired_end_reads")
    reads2 = Reads("/path/to/second_paired_end_reads")

    # filter lists of reads
    for reads in reads1, reads2:
        reads.filter_unmapped(queue=True)
        reads.filter_non_unique(strict=False, queue=True)
        reads.filter_quality(30, queue=True)
        reads.run_queued_filters()

    # load genome object
    genome = kaic.data.genomic.Genome.from_folder("/path/to/fasta_folder/")
    # extract genomic regions delineated
    # by restriction fragments
    fragments = genome.get_regions('HindIII')

    # pair up reads and assign to fragments
    pairs = FragmentMappedReadPairs(file_name="/path/to/save_file")
    pairs.load(reads1, reads2, fragments)

    # investigate error structure
    pairs.plot_ligation_error()
    # ... choose cutoffs

    # filter pairs
    pairs.filter_inward(10000)
    pairs.filter_outward(25000)
    pairs.filter_re_dist(500)

    # wrap up
    pairs.close()

"""

from __future__ import division
import tables as t
import pysam
from kaic.tools.general import RareUpdateProgressBar
from kaic.tools.files import is_sambam_file
from kaic.data.general import Maskable, MaskFilter, MaskedTable, FileBased
import os
from tables.exceptions import NoSuchNodeError
from abc import abstractmethod, ABCMeta
from bisect import bisect_right
from kaic.tools.general import bit_flags_from_int, CachedIterator
from kaic.tools.lru import lru_cache
from kaic.data.genomic import RegionsTable, GenomicRegion
import msgpack as pickle
import numpy as np
import hashlib
from functools import partial
from collections import defaultdict
from future.utils import with_metaclass
import logging
logger = logging.getLogger(__name__)


class Reads(Maskable, FileBased):
    """
    Load and filter mapped reads from a SAM/BAM file.

    .. code:: python

        # build object from SAM file
        reads = Reads(sambam_file="path/to/sam_file",
                      file_name="path/to/save_file")
        # apply built-in filters
        reads.filter_quality(30)
        reads.filter_non_unique()
        # ... continue working with object
        reads.close()

        # load previously saved Reads object from file
        reads = Reads("path/to/reads_file")
        # ... continue working with object
        reads.close()

    This class can load mapped reads from SAM or BAM files (internally uses
    pysam to process files). It will read and save all read attributes as
    provided in the file. Reads are automatically sorted by name using
    samtools.

    (Read)Filters can be applied to selectively filter out reads that match
    certain criteria. There are two built-in filters that can be accessed
    through their respective methods:

    QualityFilter (filter_quality): Filters reads that have a mapping quality
    lower than the specified cutoff

    UniquenessFilter (filter_non_unique): Filters reads that are not uniquely
    mapping, i.e. that have the XS tag

    It is also possible to create custom filters by inheriting the ReadFilter
    class and overriding the abstract valid_read method. See ReadFilter class
    for details.

    Important note: Read mapping positions (reference coordinates) are
    1-based. Pysam by default maps coordinaes to a 0-based system, this will
    be overridden.
    """

    _classid = 'READS'

    class ReadsDefinition(t.IsDescription):
        ix = t.Int32Col(pos=0)
        qname = t.Int32Col(pos=1, dflt=-1)
        flag = t.Int32Col(pos=2)
        ref = t.Int32Col(pos=3)
        pos = t.Int64Col(pos=4)
        mapq = t.Int32Col(pos=5)
        cigar = t.Int32Col(pos=6, dflt=-1)
        rnext = t.Int32Col(pos=7)
        pnext = t.Int32Col(pos=8)
        tlen = t.Int32Col(pos=9)
        seq = t.Int32Col(pos=10, dflt=-1)
        qual = t.Int32Col(pos=11, dflt=-1)
        tags = t.Int32Col(pos=12, dflt=-1)
        qname_ix = t.Float64Col(pos=13, dflt=-1)

    def __init__(self, sambam_file=None, file_name=None, mode='a',
                 _group_name='reads', mapper=None, tmpdir=None):
        """
        Create Reads object and optionally load SAM file.

        :param sambam_file: A SAM or BAM file. Tested with Bowtie2 output.
                            If provided will load SAM file into the object
                            at the end of initialization. This is the
                            recommended way of loading from file, since this
                            will automatically determine the space necessary
                            for saving read names and sequences.
                            Note: If the path does not lead to a SAM/BAM file
                            but to an HDF5 dict this method will attempt to
                            load the dict as a Reads object. Convenient for
                            loading previously created Reads objects.
        :param file_name: Either the path to an existing, previously created
                          Reads object file - in this case the object will
                          be loaded. Or the name or a non-existing file that
                          will be used to save the Reads object. If not
                          provided, all operations are performed in memory.
        :param mode: File mode ('r' = read-only, 'w' = (over)write, 'a' = append)
        :param _group_name: (internal) Name for the HDF5 group that will house
                            the Reads object's tables.
        :param mapper: Mapper that was used to align the reads. If None, will
                       try to autodetect from SAM header. Current valid mapper
                       values are ['bowtie2', 'bwa']. If other, default algorithms
                       will be used for filters.
        :return: Reads
        """

        if (sambam_file is not None and
                file_name is None and not
                is_sambam_file(sambam_file)):
            file_name = sambam_file
            sambam_file = None

        FileBased.__init__(self, file_name, mode=mode, tmpdir=tmpdir)
        Maskable.__init__(self, self.file)

        self._row_counter = {
            'reads': 0,
            'tags': 0,
            'cigar': 0,
            'qname': 0,
            'qual': 0,
            'seq': 0
        }

        # try to retrieve existing tables
        # Reads group
        try:
            self._file_group = self.file.get_node('/' + _group_name)
        except NoSuchNodeError:
            # create reads group
            self._file_group = self.file.create_group("/", _group_name, 'Reads group',
                                                      filters=t.Filters(complib="blosc",
                                                                        complevel=2, shuffle=True))
        # Main table
        try:
            self._reads = self._file_group.main
            self._row_counter['reads'] = len(self._reads)
        except NoSuchNodeError:
            self._reads = MaskedTable(self._file_group, 'main', Reads.ReadsDefinition)

        # Header attribute
        try:
            self._header = self._reads._v_attrs.header
        except AttributeError:
            logger.warn("No header attributes found in existing table")
            self._header = None

        # Reference names
        try:
            self._ref = self._reads._v_attrs.ref
        except AttributeError:
            logger.warn("No ref attributes found in existing table")
            self._ref = None

        # Qname table
        try:
            self._qname = self._file_group.qname
            self._row_counter['qname'] = len(self._qname)
        except NoSuchNodeError:
            self._qname = None

        # Cigar table
        try:
            self._cigar = self._file_group.cigar
            self._row_counter['cigar'] = len(self._cigar)
        except NoSuchNodeError:
            self._cigar = None

        # Seq table
        try:
            self._seq = self._file_group.seq
            self._row_counter['seq'] = len(self._seq)
        except NoSuchNodeError:
            self._seq = None

        # Qual table
        try:
            self._qual = self._file_group.qual
            self._row_counter['qual'] = len(self._qual)
        except NoSuchNodeError:
            self._qual = None

        # Tags table
        try:
            self._tags = self._file_group.tags
            self._row_counter['tags'] = len(self._tags)
        except NoSuchNodeError:
            self._tags = None

        # mapper
        self.mapper = mapper

        # load reads
        if sambam_file and is_sambam_file(sambam_file):
            self.load(sambam_file, ignore_duplicates=True, mapper=mapper)

    @staticmethod
    def sambam_size(sambam):
        """
        Determine number of reads in SAM/BAM file.

        :param sambam: SAM/BAM file. Can be a string (path to file) or a
                       pysam AlignmentFile.
        :return: (int) read count
        """
        if type(sambam) == str:
            sambam = pysam.AlignmentFile(sambam, 'rb')

        count = sum(1 for _ in iter(sambam))
        sambam.close()
        return count

    @property
    def mapper(self):
        return self._mapper

    @mapper.setter
    def mapper(self, mapper=None):
        if mapper:
            self._mapper = mapper
            logger.info('Mapper was explicitly set to {}'.format(self._mapper))
        else:
            try:
                self._mapper = self.header['PG'][0]['ID']
                logger.info('Mapper was detected to be {}'.format(self._mapper))
            except (KeyError, AttributeError, TypeError):
                self._mapper = None
                logger.warn('Could not auto-detect mapping program from SAM header')

    def load(self, sambam, ignore_duplicates=True,
             store_qname=True, store_cigar=True,
             store_seq=True, store_tags=True,
             store_qual=True, sample_size=None, mapper=None):
        """
        Load mapped reads from SAM/BAM file.

        Will sort mapped reads by name if not explicitly deactivated.
        By default, successive reads with duplicate names will be ignored.

        :param sambam: A string that describes the path to the SAM/BAM file
                       to be loaded or a pysam AlignmentFile.
        :param ignore_duplicates: If True (default) will not load reads that
                                  have the same name as the last loaded read.
                                  Will load all reads if False.
        """

        logger.info("Loading SAM/BAM file")
        # get file names
        try:
            file_name = sambam.filename
        except AttributeError:
            file_name = sambam

        sambam = pysam.AlignmentFile(file_name, 'rb')

        # count number of reads
        logger.info("Counting number of reads...")
        n_reads = sum(1 for _ in sambam)
        sambam.close()
        sambam = pysam.AlignmentFile(file_name, 'rb')
        logger.info("Done.")

        logger.info("Estimating field sizes")
        qname_length, seq_length, cigar_length, tags_length = Reads.determine_field_sizes(file_name, sample_size,
                                                                                          store_qname=True,
                                                                                          store_cigar=True,
                                                                                          store_seq=True,
                                                                                          store_tags=True,
                                                                                          store_qual=True)
        if sample_size is not None:
            qname_length *= 2
            seq_length *= 2
            cigar_length *= 2
            tags_length *= 2

        # create string tables if they do not yet exist
        if self._qname is None and store_qname:
            self._qname = self.file.create_earray(self._file_group, 'qname',
                                                  t.StringAtom(itemsize=qname_length), (0,))

        if self._cigar is None and store_cigar:
            self._cigar = self.file.create_earray(self._file_group, 'cigar',
                                                  t.StringAtom(itemsize=cigar_length), (0,))

        if self._seq is None and store_seq:
            self._seq = self.file.create_earray(self._file_group, 'seq',
                                                t.StringAtom(itemsize=seq_length), (0,))

        if self._qual is None and store_qual:
            self._qual = self.file.create_earray(self._file_group, 'qual',
                                                 t.StringAtom(itemsize=seq_length), (0,))

        if self._tags is None and store_tags:
            self._tags = self.file.create_earray(self._file_group, 'tags',
                                                 t.StringAtom(itemsize=tags_length), (0,))

        # header
        self._reads._v_attrs.header = {k: sambam.header[k]
                                       for k in sambam.header
                                       if k in ('HD', 'RG', 'PG')}
        self._header = sambam.header

        # Tool used to map reads
        self.mapper = mapper
        if self.mapper == 'bwa':
            ignore_duplicates = False
        
        # references
        self._reads._v_attrs.ref = sambam.references
        self._ref = sambam.references

        logger.info("Loading mapped reads...")
        with RareUpdateProgressBar(max_value=n_reads) as pb:
            last_name = ""
            for i, read in enumerate(sambam):
                if i % 1000000 == 0:
                    self.flush(update_index=False, update_csi=False)
                if ignore_duplicates and read.qname == last_name:
                    continue
                self.add_read(read, flush=False, store_cigar=store_cigar,
                              store_seq=store_seq, store_qual=store_qual,
                              store_qname=store_qname, store_tags=store_tags)
                last_name = read.qname
                pb.update(i)
        self.flush()

        logger.info("Done.")

    def _update_csi(self):
        if not self._reads.cols.qname_ix.is_indexed:
            logger.info("Sorting on qname_ix...")
            self._reads.cols.qname_ix.create_csindex()
            logger.info("Done.")
        elif not self._reads.cols.qname_ix.index.is_csi:
            logger.info("qname_ix sorting is stale, reindexing...")
            self._reads.cols.qname_ix.reindex()
            logger.info("Done.")

    def _is_sorted(self):
        if (self._reads.cols.qname_ix.index is None or
                not self._reads.cols.qname_ix.index.is_csi):
            return False
        return True

    @staticmethod
    def determine_field_sizes(sambam, sample_size=10000,
                              store_qname=True, store_cigar=True,
                              store_seq=True, store_tags=True,
                              store_qual=True):
        """
        Determine the sizes of relevant fields in a SAM/BAM file.

        :param sambam: A string that describes the path to the SAM/BAM file
                       to be loaded or a pysam AlignmentFile.
        :param sample_size: Number of lines to sample to determine field
                            sizes.
        :return: qname length, sequence length, cigar length, tags length
        """
        if type(sambam) == str:
            sambam = pysam.AlignmentFile(sambam, 'rb')

        qname_length = 0
        seq_length = 0
        cigar_length = 0
        tags_length = 0
        i = 0
        for r in sambam:
            i += 1
            if store_qname:
                qname_length = max(qname_length, len(r.qname))
            if store_seq or store_qual:
                seq_length = max(seq_length, len(r.seq))
            if store_cigar:
                cigar = r.cigar
                if cigar is not None:
                    cigar_dump = pickle.dumps(cigar)
                    cigar_length = max(cigar_length, len(cigar_dump))
            if store_tags:
                tags = r.tags
                if tags is not None:
                    tags_dump = pickle.dumps(tags)
                    tags_length = max(tags_length, len(tags_dump))

            if sample_size is not None and i >= sample_size:
                break

            if i % 100000 == 0:
                logger.info(i)
        sambam.close()

        return qname_length+10, seq_length+10, cigar_length+10, tags_length+10

    def add_read(self, read, flush=True,
                 store_qname=True, store_cigar=True,
                 store_seq=True, store_tags=True,
                 store_qual=True):
        """
        Add a read with all its attributes to this object.

        :param read: Generally a pysam.AlignedSegment object. However, can be
                     any object that provides the same fields (qname, flag,
                     reference_id, pos, mapq, rnext, pnext, tlen, seq, qual),
                     such as a Read object from this package.
        :param flush: If True, read information will be flushed from the buffer
                      and written to the object directly. When adding single
                      reads, this is recomended.
                      However, when adding reads in bulk, it is generally
                      faster to set this to False and call the flush() method
                      directly after the import.
        """
        reads_row = self._reads.row

        # main read info
        reads_row['flag'] = read.flag
        reads_row['ref'] = read.reference_id
        if read.pos >= 0:
            if isinstance(read, pysam.AlignedRead) or isinstance(read, pysam.AlignedSegment):
                reads_row['pos'] = read.pos+1
            else:
                reads_row['pos'] = read.pos
        else:
            reads_row['pos'] = read.pos
        reads_row['mapq'] = read.mapq
        reads_row['rnext'] = read.rnext
        reads_row['pnext'] = read.pnext
        reads_row['tlen'] = read.tlen

        # string info
        qname = read.qname
        if store_qname:
            self._qname.append([qname])
            reads_row['qname'] = self._row_counter['qname']
            self._row_counter['qname'] += 1
        reads_row['qname_ix'] = float(int(hashlib.md5(qname).hexdigest(), 16))

        cigar = read.cigar
        if store_cigar and cigar is not None:
            self._cigar.append([pickle.dumps(cigar)])
            reads_row['cigar'] = self._row_counter['cigar']
            self._row_counter['cigar'] += 1

        if store_seq:
            self._seq.append([read.seq])
            reads_row['seq'] = self._row_counter['seq']
            self._row_counter['seq'] += 1

        if store_qual:
            self._qual.append([read.qual])
            reads_row['qual'] = self._row_counter['qual']
            self._row_counter['qual'] += 1

        tags = read.tags
        if store_tags and tags is not None:
            tags_dump = pickle.dumps(tags)
            self._tags.append([tags_dump])
            reads_row['tags'] = self._row_counter['tags']
            self._row_counter['tags'] += 1

        reads_row['ix'] = self._row_counter['reads']
        reads_row.append()
        self._row_counter['reads'] += 1

        if flush:
            self.flush()

    def flush(self, update_index=True, update_csi=True):
        """
        Write the latest changes to this object to file.
        """
        self._reads.flush(update_index=update_index)
        if update_csi:
            self._update_csi()
        if self._tags is not None:
            self._tags.flush()
            self._row_counter['tags'] = len(self._tags)
        if self._cigar is not None:
            self._cigar.flush()
            self._row_counter['cigar'] = len(self._cigar)
        if self._seq is not None:
            self._seq.flush()
            self._row_counter['seq'] = len(self._seq)
        if self._qual is not None:
            self._qual.flush()
            self._row_counter['qual'] = len(self._qual)
        if self._qname is not None:
            self._qname.flush()
            self._row_counter['qname'] = len(self._qname)

    @property
    def header(self):
        """
        Retrieve the header of a loaded SAM/BAM file.

        :return: The SAM/BAM header as extracted by pysam.
        """
        return self._header

    @property
    def chromosomes(self):
        return self._ref

    def _ix2ref(self, ix):
        """
        Convert a reference_id (ix) to the reference name.
        """
        if self._ref is None:
            raise RuntimeError("Chromosome reference for left read not present")
        return self._ref[ix]

    def _row2read(self, row, lazy=False):
        """
        Convert a row from the internal _reads pytables table to Read object.
        """
        if lazy:
            return LazyRead(row, self)

        tags = None
        tags_ix = row['tags']
        if tags_ix >= 0:
            tags_str = self._tags[tags_ix]
            try:
                tags = pickle.loads(tags_str)
            except pickle.UnpackValueError:
                tags = pickle.loads(tags_str + '\x00')

        cigar = None
        cigar_ix = row['cigar']
        if cigar_ix >= 0:
            cigar_str = self._cigar[cigar_ix]
            try:
                cigar = pickle.loads(cigar_str)
            except pickle.UnpackValueError:
                cigar = pickle.loads(cigar_str + '\x00')

        qname = None
        qname_ix = row['qname']
        if qname_ix >= 0:
            qname = self._qname[qname_ix]

        qual = None
        qual_ix = row['qual']
        if qual_ix >= 0:
            qual = self._qual[qual_ix]

        seq = None
        seq_ix = row['seq']
        if seq_ix >= 0:
            seq = self._seq[seq_ix]

        ref_ix = row['ref']
        ref = self._ix2ref(ref_ix)

        return Read(qname=qname, flag=row['flag'], ref=ref,
                    pos=row['pos'], mapq=row['mapq'], cigar=cigar, rnext=row['rnext'],
                    pnext=row['pnext'], tlen=row['tlen'], seq=seq, qual=qual,
                    tags=tags, reference_id=ref_ix, qname_ix=row['qname_ix'])

    def reads(self, lazy=False, sort_by_qname_ix=False, excluded_filters=()):
        """
        Iterate over _reads table and convert each result to Read.

        :param lazy: Lazily load read properties (only works inside loop!)
        :param sort_by_qname_ix: Iterate by ascending qname_ix
        :return: ReadsIter that iterates over visible reads
        """

        excluded_masks = self.get_binary_mask_from_masks(excluded_filters)

        # ensure sorting on qname_ix column
        if sort_by_qname_ix:
            if not self._is_sorted():
                try:
                    logger.info("Sorting qnames...")
                    self._update_csi()
                except t.exceptions.FileModeError:
                    raise RuntimeError("This object is not sorted by qname_ix! "
                                       "Cannot sort manually, because file is in read-only mode.")
            it = self._reads.itersorted('qname_ix', excluded_masks=excluded_masks)
        else:
            it = self._reads.iterrows(excluded_masks=excluded_masks)

        return (self._row2read(row, lazy=lazy) for row in it)

    def __iter__(self):
        return self.reads()

    def __getitem__(self, key):
        """
        Get a Read by index or slice.
        """
        if isinstance(key, int):
            return self._row2read(self._reads[key])
        if isinstance(key, slice):
            reads = []
            for row in self._reads[key]:
                reads.append(self._row2read(row))
            return reads
        raise KeyError("Key %s not supported" % str(key))

    def __len__(self):
        """
        Return number of (unmasked) reads in object.
        """
        return len(self._reads)

    def where(self, query):
        """
        Search through reads using queries.

        .. code:: python

            Example:

            result = reads.where("(mapq < 30) & (flag == 0)")

        :param query: A query string in pytables query format (see
                      PyTables documentation).
        :return: A list of matching Read objects.
        """
        reads = []
        for row in self._reads.where(query):
            reads.append(self._row2read(row))
        return reads

    def get_read_by_qname(self, qname):
        for read in self.reads(lazy=True):
            if read.qname == qname:
                return read
        return None

    def filter(self, read_filter, queue=False, log_progress=False):
        """
        Filter reads using a ReadFilter object.

        :param read_filter: Class implementing ReadFilter. Must override
                            valid_read method, ideally sets mask parameter
                            during initialization.
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        :param log_progress: If true, process iterating through all reads
                             will be continuously reported.
        """
        read_filter.set_reads_object(self)
        if not queue:
            self._reads.filter(read_filter, _logging=log_progress)
        else:
            self._reads.queue_filter(read_filter)

    def filter_quality(self, cutoff=30, queue=False):
        """
        Convenience function that applies a QualityFilter.
        The actual algorithm and rationale used for filtering will depend on the
        internal _mapper attribute.

        :param cutoff: Minimum mapping quality a read must have to pass
                       the filter
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        """
        mask = self.add_mask_description('mapq', 'Mask read pairs with a mapping quality lower than %d' % cutoff)
        if self._mapper == 'bwa':
            quality_filter = BwaMemQualityFilter(cutoff, mask)
        else:
            quality_filter = QualityFilter(cutoff, mask)
        self.filter(quality_filter, queue)

    def filter_unmapped(self, queue=False):
        """
        Convenience function that applies an UnmappedFilter.

        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        """
        mask = self.add_mask_description('unmapped', 'Mask read pairs that are unmapped')
        unmapped_filter = UnmappedFilter(mask)
        self.filter(unmapped_filter, queue)
            
    def filter_non_unique(self, strict=True, cutoff=3, queue=False):
        """
        Convenience function that applies a UniquenessFilter.
        The actual algorithm and rationale used for filtering will depend on the
        internal _mapper attribute.

        :param cutoff: Minimum mapq value if the alignments are from bwa-mem.
        :param strict: If True will filter if XS tag is present. If False,
                       will filter only when XS tag is not 0. This is applied if
                       alignments are from bowtie2.
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        """
        if self._mapper == 'bwa':
            mask = self.add_mask_description('uniqueness', 'Mask reads that do not map uniquely (mapq <= {})'.format(cutoff))
            uniqueness_filter = BwaMemUniquenessFilter(cutoff, mask)
        else:
            mask = self.add_mask_description('uniqueness', 'Mask reads that do not map uniquely (according to XS tag)')
            uniqueness_filter = UniquenessFilter(strict, mask)
        self.filter(uniqueness_filter, queue)

    def run_queued_filters(self, log_progress=False):
        """
        Run queued filters.

        :param log_progress: If true, process iterating through all reads
                             will be continuously reported.
        """
        self._reads.run_queued_filters(_logging=log_progress)

    def filtered_reads(self):
        """
        Iterate over filtered reads.

        Reads are never deleted when using a filter, only masked.
        This function iterates over all masked reads.
        """
        this = self

        class MaskedReadsIter:
            def __init__(self):
                self.iter = this._reads.masked_rows()

            def __iter__(self):
                return self

            def next(self):
                row = self.iter.next()
                read = this._row2read(row)

                masks = this.get_masks(row[this._reads._mask_field])

                return MaskedRead(qname=read.qname, flag=read.flag, ref=read.ref,
                                  pos=read.pos, mapq=read.mapq, cigar=read.cigar, rnext=read.rnext,
                                  pnext=read.pnext, tlen=read.tlen, seq=read.seq, qual=read.qual,
                                  tags=read.tags, reference_id=read.reference_id, qname_ix=read.qname_ix,
                                  masks=masks)
        return MaskedReadsIter()


class Read(object):
    """
    Object representing a mapped sequencing read.

    This is the main class used by a Reads object to represent mapped
    read information.
    """

    def __init__(self, qname=None, flag=0, ref="",
                 pos=0, mapq=0, cigar=None, rnext=0,
                 pnext=0, tlen=0, seq=None, qual=None,
                 tags=None, reference_id=None, qname_ix=None):
        """
        Initialize a Read with specific attributes.

        :param qname: String identifier (generally unique) of the Read
        :param flag: Base 2 bit flag. Meaning can differ depending on
                     mapper output.
        :param ref: Reference sequence name
        :param pos: Mapping position of read along reference sequence
                    (1-based)
        :param mapq: Mapping quality
        :param cigar: CIGAR string
        :param rnext: Reference sequence name for the next mapped read.
        :param pnext: Mapping position of the next read
        :param tlen: Signed observed template length
        :param seq: Sequence of the segment/mapped read
        :param qual: Quality string
        :param tags: Dictionary of tag-value pairs that provide
                     additional information about the alignment
        :param reference_id: ID (integer) of the reference sequence
        :param qname_ix: qname converted into a unique float representation
        """

        self.qname = qname
        self.flag = flag
        self.ref = ref
        self.pos = pos
        self.mapq = mapq
        self.cigar = cigar
        self.rnext = rnext
        self.pnext = pnext
        self.tlen = tlen
        self.seq = seq
        self.qual = qual
        self.tags = tags
        self.reference_id = reference_id
        self.qname_ix = qname_ix

    @property
    def strand(self):
        """
        Get the strand that this read is aligned to.

        :return: -1 of if aligned to the negative strand,
                 +1 if alinged to the positive strand
        """
        bit_flags = bit_flags_from_int(self.flag)
        if 4 in bit_flags:
            return -1
        return 1

    @property
    def alen(self):
        """
        Return the length of the aligned portion of the read
        """
        valids = [0]
        return sum([i[1] for i in self.cigar if i[0] in valids])

    def get_tag(self, key):
        """
        Return the value of a tag. None if does not exist

        :param key: Key/name of alignment tag
        :return: Value of tag of none if tag not present
        """
        for tag in self.tags:
            if tag[0] == key:
                return tag[1]
        return None

    def __getitem__(self, key):
        """
        Retrieve attribute with bracket notation for convenience.
        """
        try:
            value = self.__getattribute__(key)
            return value
        except:
            raise KeyError("Read does not have %s attribute" % str(key))

    def __repr__(self):
        return "%s, ref: %s, pos: %d" % (self.qname, self.ref, self.pos)


class MaskedRead(Read):
    """
    Read that is masked for some reason.

    Provides all the functionality of Read plus access to the
    masks that are hiding this read in the Reads object.
    """
    def __init__(self, qname=None, flag=0, ref="",
                 pos=0, mapq=0, cigar=None, rnext=0,
                 pnext=0, tlen=0, seq=None, qual=None,
                 tags=None, reference_id=None,
                 qname_ix=None, masks=None):
        """
        Initialize a MaskedRead with specific attributes.

        :param qname: String identifier (generally unique) of the Read
        :param flag: Base 2 bit flag. Meaning can differ depending on
                     mapper output.
        :param ref: Reference sequence name
        :param pos: Mapping position of read along reference sequence
                    (1-based)
        :param mapq: Mapping quality
        :param cigar: CIGAR string
        :param rnext: Reference sequence name for the next mapped read.
        :param pnext: Mapping position of the next read
        :param tlen: Signed observed template length
        :param seq: Sequence of the segment/mapped read
        :param qual: Quality string
        :param tags: Dictionary of tag-value pairs that provide
                     additional information about the alignment
        :param reference_id: ID (integer) of the reference sequence
        :param masks: A list of Mask objects
        """
        super(MaskedRead, self).__init__(qname=qname, flag=flag, ref=ref,
                                         pos=pos, mapq=mapq, cigar=cigar, rnext=rnext,
                                         pnext=pnext, tlen=tlen, seq=seq, qual=qual,
                                         tags=tags, reference_id=reference_id, qname_ix=qname_ix)
        self.masks = masks

    def __repr__(self):
        representation = super(MaskedRead, self).__repr__()
        if self.masks is not None:
            mask_names = []
            for mask in self.masks:
                mask_names.append(mask.name)
            return "%s (%s)" % (representation, ", ".join(mask_names))
        return representation


class LazyRead(Read):
    def __init__(self, row, parent):
        self.row = row
        self.parent = parent

    @property
    def flag(self):
        return self.row["flag"]

    @property
    def ref(self):
        return self.parent._ix2ref(self.row['ref'])

    @property
    def pos(self):
        return self.row["pos"]

    @property
    def mapq(self):
        return self.row["mapq"]

    @property
    def rnext(self):
        return self.row["rnext"]

    @property
    def pnext(self):
        return self.row["pnext"]

    @property
    def tlen(self):
        return self.row["tlen"]

    @property
    def cigar(self):
        ix = self.row["cigar"]
        if ix >= 0:
            cigar = self.parent._cigar[ix]
            try:
                return pickle.loads(cigar)
            except pickle.UnpackValueError:
                return pickle.loads(cigar + '\x00')
        return None

    @property
    def qname(self):
        ix = self.row["qname"]
        if ix >= 0:
            return self.parent._qname[ix]
        return None

    @property
    def seq(self):
        ix = self.row["seq"]
        if ix >= 0:
            return self.parent._seq[ix]
        return None

    @property
    def qual(self):
        ix = self.row["qual"]
        if ix >= 0:
            return self.parent._qual[ix]
        return None

    @property
    def tags(self):
        ix = self.row["tags"]
        if ix >= 0:
            tags = self.parent._tags[ix]
            try:
                return pickle.loads(tags)
            except pickle.UnpackValueError:
                return pickle.loads(tags + '\x00')
        return None

    @property
    def reference_id(self):
        return self.row['ref']

    @property
    def qname_ix(self):
        return self.row['qname_ix']


#
# Filters Reads
#
class ReadFilter(with_metaclass(ABCMeta, MaskFilter)):
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

    @abstractmethod
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
        pass

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
        :param contaminant_reads: A :class:`~Reads` object representing a
                                  contaminant
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered reads.
        """
        super(ContaminantFilter, self).__init__(mask)

        self.contaminant_names = set()
        for read in contaminant_reads.reads(lazy=True):
            self.contaminant_names.add("%.0f" % read.qname_ix)

    def valid_read(self, read):
        """
        Check if a read has a mapq >= cutoff.
        """
        if "%.0f" % read.qname_ix in self.contaminant_names:
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
        xs_tag = read.get_tag('XS')
        if xs_tag is not None and (self.strict or xs_tag != 0):
            return False
        return True


class BwaMemUniquenessFilter(ReadFilter):
    """
    Filters `bwa mem` generated alignements based on whether they are unique or not.
    The presence of a non-zero XS tag does not mean a read is a multi-mapping one.
    Instead, we make sure that the ratio XS/AS is inferior to a certain threshold.
    """
    def __init__(self, cutoff=3, mask=None):
        """
        :param cutoff: MAPQ score above which an alignment should be in order to be considered
                       unique. In practice, a cutoff of 3 ensures that no other alignment is
                       equally good.
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered reads.
        """
        super(BwaMemUniquenessFilter, self).__init__(mask)
        self.cutoff = cutoff

    def valid_read(self, read):
        """
        Check if a read has a high alignment score.
        """
        return read.mapq > self.cutoff


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
        Check if the flag attribute contains 2**2.
        """
        if 2 in bit_flags_from_int(read.flag):
            return False
        return True


class PairLoader(with_metaclass(ABCMeta, object)):

    class Aln(object):
        def __init__(self, row):
            self.qname_ix = row.qname_ix
            self.pos = row.pos
            self.strand = row.strand
            self.ref = row.ref
            self.flag = row.flag
            self.cigar = row.cigar

    def __init__(self, pairs, ignore_duplicates=True, _in_memory_index=True, excluded_filters=()):
        self._pairs = pairs
        self._reads1 = None
        self._reads2 = None
        self.fragment_ends = None
        self.fragment_infos = None
        self.add_read_single = self._pairs.add_read_single
        self.add_read_pair = self._pairs.add_read_pair
        self._in_memory_index = _in_memory_index
        self.ignore_duplicates = ignore_duplicates
        self.excluded_filters = excluded_filters

    def load(self, reads1, reads2, regions=None):
        self._reads1 = reads1
        self._reads2 = reads2
        reads1, reads2, regions, add_read_single, add_read_pair = self.setup(reads1, reads2, regions)
        self.load_pairs_from_reads(reads1, reads2, regions, add_read_single, add_read_pair)
        self._pairs.flush(update_index=True)

    def setup(self, reads1, reads2, regions):
        if regions is not None:
            logger.info("Adding regions...")
            self._pairs.add_regions(regions)
            logger.info("Done.")

        # generate index for fragments
        logger.info("Generating region index...")
        fragment_infos = None
        fragment_ends = None
        if self._in_memory_index:
            fragment_infos = {}
            fragment_ends = {}
            for region in self._pairs.regions():
                if region.chromosome not in fragment_ends:
                    fragment_ends[region.chromosome] = []
                    fragment_infos[region.chromosome] = []
                fragment_infos[region.chromosome].append((region.ix,
                                                          self._pairs._chromosome_to_ix[region.chromosome],
                                                          region.start, region.end))
                fragment_ends[region.chromosome].append(region.end)

        if isinstance(self._reads1, str):
            logger.info("Loading reads 1")
            reads1 = Reads(sambam_file=self._reads1)

        if isinstance(self._reads2, str):
            logger.info("Loading reads 2")
            reads2 = Reads(sambam_file=self._reads2)

        self.fragment_ends = fragment_ends
        self.fragment_infos = fragment_infos

        add_read_single = partial(self._pairs.add_read_single, flush=False,
                                  _fragment_ends=self.fragment_ends,
                                  _fragment_infos=self.fragment_infos)
        add_read_pair = partial(self._pairs.add_read_pair, flush=False,
                                _fragment_ends=self.fragment_ends,
                                _fragment_infos=self.fragment_infos)

        return reads1, reads2, regions, add_read_single, add_read_pair

    @abstractmethod
    def load_pairs_from_reads(self, reads1, reads2, regions, add_read_single, add_read_pair):
        pass


class Bowtie2PairLoader(PairLoader):
    def __init__(self, pairs, ignore_duplicates=True, _in_memory_index=True, excluded_filters=None):
        super(Bowtie2PairLoader, self).__init__(pairs,
                                                ignore_duplicates=ignore_duplicates,
                                                _in_memory_index=_in_memory_index,
                                                excluded_filters=excluded_filters)

    @staticmethod
    def get_next_read(iterator):
        try:
            r = iterator.next()
            return r
        except StopIteration:
            return None

    def load_pairs_from_reads(self, reads1, reads2, regions, add_read_single, add_read_pair):
        logger.info("Adding read pairs...")

        iter1 = reads1.reads(lazy=True, sort_by_qname_ix=True, excluded_filters=self.excluded_filters)
        iter2 = reads2.reads(lazy=True, sort_by_qname_ix=True, excluded_filters=self.excluded_filters)

        # add and map reads
        i = 0
        last_r1_name_ix = ''
        last_r2_name_ix = ''
        r1 = self.get_next_read(iter1)
        r2 = self.get_next_read(iter2)
        r1_count = 0
        r2_count = 0
        pair_count = 0
        single_count = 0

        total = len(reads1) + len(reads2)

        with RareUpdateProgressBar(max_value=total) as pb:
            while r1 is not None and r2 is not None:
                i += 1
                if r1.qname_ix == last_r1_name_ix:
                    if not self.ignore_duplicates:
                        raise ValueError("Duplicate left read QNAME %s" % r1.qname)
                    r1 = self.get_next_read(iter1)
                    r1_count += 1
                elif r2.qname_ix == last_r2_name_ix:
                    if not self.ignore_duplicates:
                        raise ValueError("Duplicate right read QNAME %s" % r2.qname)
                    r2 = self.get_next_read(iter2)
                    r2_count += 1
                elif abs(r1.qname_ix-r2.qname_ix) < 0.5:
                    add_read_pair(r1, r2)
                    last_r1_name_ix = r1.qname_ix
                    last_r2_name_ix = r2.qname_ix
                    r1 = self.get_next_read(iter1)
                    r2 = self.get_next_read(iter2)
                    r1_count += 1
                    r2_count += 1
                    pair_count += 1
                elif r1.qname_ix-r2.qname_ix < 0:
                    add_read_single(r1)
                    last_r1_name_ix = r1.qname_ix
                    r1 = self.get_next_read(iter1)
                    r1_count += 1
                    single_count += 1
                else:
                    add_read_single(r2)
                    last_r2_name_ix = r2.qname_ix
                    r2 = self.get_next_read(iter2)
                    r2_count += 1
                    single_count += 1

                pb.update(r1_count + r2_count)

                if i % 1000000 == 0:
                    self._pairs.flush(update_index=False)

            # add remaining unpaired reads
            while r1 is not None:
                i += 1
                if r1.qname_ix == last_r1_name_ix:
                    if not self.ignore_duplicates:
                        raise ValueError("Duplicate left read QNAME %s" % r1.qname)
                else:
                    add_read_single(r1)
                last_r1_name_ix = r1.qname_ix
                r1 = self.get_next_read(iter1)
                r1_count += 1
                single_count += 1

                pb.update(r1_count + r2_count)

                if i % 1000000 == 0:
                    self._pairs.flush(update_index=False)

            while r2 is not None:
                i += 1
                if r2.qname_ix == last_r2_name_ix:
                    if not self.ignore_duplicates:
                        raise ValueError("Duplicate right read QNAME %s" % r2.qname)
                else:
                    add_read_single(r2)
                last_r2_name_ix = r2.qname_ix
                r2 = self.get_next_read(iter2)
                r2_count += 1
                single_count += 1

                pb.update(r1_count + r2_count)

                if i % 1000000 == 0:
                    self._pairs.flush(update_index=False)

        logger.info("Left reads: %d, right reads: %d" % (r1_count, r2_count))
        logger.info("Pairs: %d. Single: %d" % (pair_count, single_count))


class BwaMemPairLoader(PairLoader):
    def __init__(self, pairs, _in_memory_index=True, excluded_filters=()):
        self.add_read_single = None
        self.add_read_pair = None
        super(BwaMemPairLoader, self).__init__(pairs,
                                               ignore_duplicates=False,
                                               _in_memory_index=_in_memory_index,
                                               excluded_filters=excluded_filters)

    @staticmethod
    def get_all_read_alns(it):
        alns = list()
        alns.append(it.next())
        if alns[0] is not None:
            name_ix = alns[0].qname_ix
            while 1:
                r = it.next()
                if r is None:
                    break
                elif r.qname_ix == name_ix:
                    alns.append(r)
                else:
                    it.prev()
                    break
        return alns

    @staticmethod
    def sort_alns(alns):
        """
        Sorts split alignments of a read according to their order on the template.
        Thanks Clemens for the inspiration!
        """
        def _get_match_part(a):
            "Find part of alignment that is covered by M (matches)"
            cigar = a.cigar if a.strand == 1 else reversed(a.cigar)
            m = []
            cur_pos = 0
            for c in cigar:
                if c[0] == 0:
                    m.extend([cur_pos, cur_pos + c[1]])
                if c[0] in (0, 1, 4, 5): # if c[0] in 'MISH'
                    cur_pos += c[1]
            return min(m), max(m)
        # Sort alignments based on their match positions
        segments = [None] * len(alns)
        for i, a in enumerate(alns):
            m = _get_match_part(a)
            segments[i] = (m[0], m[1], i)
        sorted_segments = sorted(segments, key=lambda x: x[0])
        return [alns[s[2]] for s in sorted_segments]

    def process_bwa_alns(self, head=None, tail=None):
        """
        Decide how to handle alignments
        """
        if head and tail:
            if not(len(head) == len(tail) == 1):
                head = self.sort_alns(head)
                tail = self.sort_alns(tail)
            self.add_read_pair(head[0], tail[0])
        else:
            alns = head if head is not None else tail
            if len(alns) == 1:
                self.add_read_single(alns[0])
            elif len(alns) == 2:
                self.add_read_pair(*alns)
            else:
                pass

    def load_pairs_from_reads(self, reads1, reads2, regions, add_read_single, add_read_pair):
        logger.info("Adding read pairs...")
        self.add_read_single = add_read_single
        self.add_read_pair = add_read_pair

        iter1 = CachedIterator(
            (PairLoader.Aln(r) for r in reads1.reads(lazy=True, sort_by_qname_ix=True)), 2
        )
        iter2 = CachedIterator(
            (PairLoader.Aln(r) for r in reads2.reads(lazy=True, sort_by_qname_ix=True)), 2
        )

        # add and map reads
        i = 0
        r1 = self.get_all_read_alns(iter1)
        r2 = self.get_all_read_alns(iter2)
        r1_count = 0
        r2_count = 0

        total = len(reads1) + len(reads2)

        with RareUpdateProgressBar(max_value=total) as pb:
            while r1[0] is not None and r2[0] is not None:
                i += 1
                if abs(r1[0].qname_ix-r2[0].qname_ix) < 0.5:
                    self.process_bwa_alns(r1, r2)
                    r1 = self.get_all_read_alns(iter1)
                    r2 = self.get_all_read_alns(iter2)
                    r1_count += 1
                    r2_count += 1
                elif r1[0].qname_ix-r2[0].qname_ix < 0:
                    self.process_bwa_alns(r1)
                    r1 = self.get_all_read_alns(iter1)
                    r1_count += 1
                else:
                    self.process_bwa_alns(r2)
                    r2 = self.get_all_read_alns(iter2)
                    r2_count += 1

                pb.update(r1_count + r2_count)

                if i % 1000000 == 0:
                    self._pairs.flush(update_index=False)

            # add remaining unpaired reads
            while r1[0] is not None:
                i += 1
                self.process_bwa_alns(r1)
                r1 = self.get_all_read_alns(iter1)
                r1_count += 1

                pb.update(r1_count + r2_count)

                if i % 1000000 == 0:
                    self._pairs.flush(update_index=False)

            while r2[0] is not None:
                i += 1
                self.process_bwa_alns(r2)
                r2 = self.get_all_read_alns(iter2)
                r2_count += 1

                pb.update(r1_count + r2_count)

                if i % 1000000 == 0:
                    self._pairs.flush(update_index=False)

        logger.info("Left reads: %d, right reads: %d" % (r1_count, r2_count))


class FragmentMappedReadPairs(Maskable, RegionsTable, FileBased):
    """
    Map pairs of reads to restriction fragments in a reference genome.

    Example:

    ..  code:: python

        # load, pair up, and map reads to restriction fragments
        pairs = FragmentMappedReadPairs("/path/to/save_file")
        genome = Genome.from_folder("/path/to/fasta_folder/")
        pairs.load(reads_object1, reads_object2, genome.get_regions('HindIII'))

        # filter read pairs
        pairs.filter_inward(10000)
        pairs.filter_outward(25000)
        pairs.filter_re_dist(500)

        # ... more processing
        pairs.close()

    This class provides methods to pair reads from two different lists, map them to
    a list of genomic regions (generally resctriction fragments obtained through the
    :class:`~kaic.data.genomic.Genome` object), and filter the list of pairs by
    different criteria.

    It provides several built-in pair filters:

    - :func:`~FragmentMappedReadPairs.filter_inward`: Filters pairs where the left read
      maps to the positive and the right read to the negative strand ("inward-facing"),
      and their distance is smaller than a specified cutoff.

    - :func:`~FragmentMappedReadPairs.filter_outward`: Filters pairs where the left read
      maps to the negative and the right read to the positive strand ("outward-facing"),
      and their distance is smaller than a specified cutoff.

    - :func:`~FragmentMappedReadPairs.filter_re_dist`: Filters pairs where one or both
      reads are further than a specified distance away from the nearest fragment border
      (= restriction site)

    Custom filters can be created by implementing the FragmentMappedReadPairFilter
    class and overriding the valid_pair method.
    """

    _classid = 'FRAGMENTMAPPEDREADPAIRS'

    class FragmentsMappedReadPairDescription(t.IsDescription):
        """
        Needed by PyTables to build read pairs table.
        """
        ix = t.Int32Col(pos=0)
        left_read_qname_ix = t.Float64Col(pos=1)
        left_read_position = t.Int64Col(pos=2)
        left_read_strand = t.Int8Col(pos=3)
        left_fragment = t.Int32Col(pos=4, dflt=-1)
        left_fragment_start = t.Int64Col(pos=5)
        left_fragment_end = t.Int64Col(pos=6)
        left_fragment_chromosome = t.Int32Col(pos=7)
        right_read_qname_ix = t.Float64Col(pos=8)
        right_read_position = t.Int64Col(pos=9)
        right_read_strand = t.Int8Col(pos=10)
        right_fragment = t.Int32Col(pos=11, dflt=-1)
        right_fragment_start = t.Int64Col(pos=12)
        right_fragment_end = t.Int64Col(pos=13)
        right_fragment_chromosome = t.Int32Col(pos=14)

    class FragmentsMappedReadSingleDescription(t.IsDescription):
        """
        Needed by PyTables to build single reads table.
        """
        ix = t.Int32Col(pos=0)
        read_qname_ix = t.Float64Col(pos=1)
        read_position = t.Int64Col(pos=2)
        read_strand = t.Int8Col(pos=3)
        fragment = t.Int32Col(pos=4, dflt=-1)
        fragment_start = t.Int64Col(pos=5)
        fragment_end = t.Int64Col(pos=6)
        fragment_chromosome = t.Int32Col(pos=7)

    def __init__(self, file_name=None,
                 mode='a',
                 group_name='fragment_map',
                 table_name_fragments='fragments',
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

        if file_name is not None and isinstance(file_name, str):
            file_name = os.path.expanduser(file_name)

        FileBased.__init__(self, file_name, mode=mode, tmpdir=tmpdir)
        RegionsTable.__init__(self, file_name=self.file, _table_name_regions=table_name_fragments)

        # generate tables from inherited classes
        Maskable.__init__(self, self.file)

        # try to retrieve existing table
        try:
            self._pairs = self.file.get_node('/' + group_name + '/mapped_read_pairs')
            self._single = self.file.get_node('/' + group_name + '/mapped_read_single')
            self._pair_count = self._pairs._original_len()
            self._single_count = self._single._original_len()
        # or build table from scratch
        except NoSuchNodeError:
            # create group
            group = self.file.create_group("/", group_name, 'Mapped read pairs group',
                                           filters=t.Filters(complib="blosc",
                                                             complevel=2, shuffle=True))

            self._pairs = MaskedTable(group, 'mapped_read_pairs',
                                      FragmentMappedReadPairs.FragmentsMappedReadPairDescription)

            self._single = MaskedTable(group, 'mapped_read_single',
                                       FragmentMappedReadPairs.FragmentsMappedReadSingleDescription)
            self._pair_count = 0
            self._single_count = 0

    def flush(self, update_index=True):
        RegionsTable.flush(self)
        self._pairs.flush(update_index=update_index)
        self._single.flush(update_index=update_index)

    def load(self, reads1, reads2, regions=None, ignore_duplicates=True, _in_memory_index=True, excluded_filters=()):
        """
        Load paired reads and map them to genomic regions (e.g. RE-fragments).

        This method contains all the necessary functionality to load data
        into the FragmentMappedReadPairs object in one go. If provided, it will
        initially store a list of genomic regions - generally restriction
        fragments. It will then advance through two lists of reads
        simultaneously (they MUST be sorted by qname either by samtools or
        by an equivalent sorting method), and pair and map them to
        the provided regions along the way.

        Reads that can't be paired, for example when one side has been
        filtered in an earlier processing step, are considered single and
        will be saved separately.

        :param reads1: Can be a :class:`~Reads` object, the path to a sorted
                       SAM/BAM file, or a sorted pysam AlignmentFile. This
                       constitutes the first half of mapped reads from paired-
                       end sequencing.
        :param reads2: Can be a :class:`~Reads` object, the path to a sorted
                       SAM/BAM file, or a sorted pysam AlignmentFile. This
                       constitutes the second half of mapped reads from
                       paired-end sequencing.
        :param regions: A list of genomic regions that will be used to assign
                        reads. Look at :class:`~kaic.data.genomic.RegionsTable`
                        for more details on allowed formats.
        :param ignore_duplicates: Will ignore all duplicates of previously
                                  loaded reads with the same name.
        :param _in_memory_index: If True (default), will load the genomic
                                 regions into memory for mapping. Set to
                                 False if memory is an issue, but be prepared
                                 for a very long runtime.
        """
        if reads1.mapper == 'bwa' and reads2.mapper == 'bwa':
            loader = BwaMemPairLoader(self, _in_memory_index, excluded_filters)
        else:
            loader = Bowtie2PairLoader(self, ignore_duplicates, _in_memory_index, excluded_filters)

        loader.load(reads1, reads2, regions=regions)

    def add_read_pair(self, read1, read2, flush=True,
                      _fragment_ends=None, _fragment_infos=None):
        """
        Add a pair of reads to this object.

        :param read1: Left side of a paired-end read.
        :param read2: Right side of a paired-end read.
        :param flush: Write changes to file immediately after adding pair.
        :param _fragment_ends: (Internal) Can be used to provide a dict
                               with chromosome names as keys and lists of
                               fragment end positions as values. Speeds up
                               mapping, but preferably use
                               :func:`~FragmentMappedReadPairs.load` instead!
        :param _fragment_infos: (Internal) See previous argument, but provide
                                fragment info ([index, chromosome index, start, end])
                                instead of end positions.
        """

        fragment_infos1 = self._find_fragment_info(read1.ref, read1.pos,
                                                   _fragment_ends=_fragment_ends, _fragment_infos=_fragment_infos)
        fragment_infos2 = self._find_fragment_info(read2.ref, read2.pos,
                                                   _fragment_ends=_fragment_ends, _fragment_infos=_fragment_infos)

        # both must be integer if successfully mapped
        if fragment_infos1 is not None and fragment_infos2 is not None:
            if fragment_infos1[0] <= fragment_infos2[0]:
                fragment_ix1, fragment_chromosome1, fragment_start1, fragment_end1 = fragment_infos1
                fragment_ix2, fragment_chromosome2, fragment_start2, fragment_end2 = fragment_infos2
            else:
                tmp_read = read1
                read1 = read2
                read2 = tmp_read
                fragment_ix1, fragment_chromosome1, fragment_start1, fragment_end1 = fragment_infos2
                fragment_ix2, fragment_chromosome2, fragment_start2, fragment_end2 = fragment_infos1

            try:
                read1_strand = read1.strand
            except AttributeError:
                read1_strand = None

            if read1_strand is None:
                bit_flags = bit_flags_from_int(read1.flag)
                if 4 in bit_flags:
                    read1_strand = -1
                else:
                    read1_strand = 1

            try:
                read2_strand = read2.strand
            except AttributeError:
                read2_strand = None

            if read2_strand is None:
                bit_flags = bit_flags_from_int(read2.flag)
                if 4 in bit_flags:
                    read2_strand = -1
                else:
                    read2_strand = 1

            row = self._pairs.row
            row['ix'] = self._pair_count
            row['left_read_qname_ix'] = read1.qname_ix
            row['left_read_position'] = read1.pos
            row['left_read_strand'] = read1_strand
            row['left_fragment'] = fragment_ix1
            row['left_fragment_start'] = fragment_start1
            row['left_fragment_end'] = fragment_end1
            row['left_fragment_chromosome'] = fragment_chromosome1
            row['right_read_qname_ix'] = read2.qname_ix
            row['right_read_position'] = read2.pos
            row['right_read_strand'] = read2_strand
            row['right_fragment'] = fragment_ix2
            row['right_fragment_start'] = fragment_start2
            row['right_fragment_end'] = fragment_end2
            row['right_fragment_chromosome'] = fragment_chromosome2

            row.append()
            self._pair_count += 1

            if flush:
                self._pairs.flush(update_index=True)

    def add_read_single(self, read, flush=True,
                        _fragment_ends=None, _fragment_infos=None):
        fragment_info = self._find_fragment_info(read.ref, read.pos,
                                                 _fragment_ends=_fragment_ends, _fragment_infos=_fragment_infos)

        try:
            read_strand = read.strand
        except AttributeError:
            read_strand = None

        if read_strand is None:
            bit_flags = bit_flags_from_int(read.flag)
            if 4 in bit_flags:
                read_strand = -1
            else:
                read_strand = 1

        if fragment_info is not None:
            row = self._single.row
            row['ix'] = self._single_count
            row['read_qname_ix'] = read.qname_ix
            row['read_position'] = read.pos
            row['read_strand'] = read_strand
            row['fragment'] = fragment_info[0]
            row['fragment_start'] = fragment_info[2]
            row['fragment_end'] = fragment_info[3]
            row['fragment_chromosome'] = fragment_info[1]
            row.append()

            if flush:
                self._single.flush(update_index=True)
            self._single_count += 1

    def _find_fragment_info(self, chromosome, position,
                            _fragment_ends=None, _fragment_infos=None):
        """
        Find the index of a fragment by genomic coordinate and chromosome name.
        """
        # binary search for fragment
        fragment_infos = None
        if _fragment_ends is not None and _fragment_infos is not None:
            try:
                pos_ix = bisect_right(_fragment_ends[chromosome], position)
                fragment_infos = _fragment_infos[chromosome][pos_ix]
            except KeyError:
                # potentially keep a record of unmatched chromosome names
                pass
        else:
            for row in self._regions.where(
                            "(start <= %d) & (end >= %d) & (chromosome == '%s')" % (position, position, chromosome)):
                fragment_infos = [row['ix'], self._chromosome_to_ix[chromosome], row['start'], row['end']]

        return fragment_infos

    def _pair_from_row(self, row, lazy=False):
        """
        Convert a pytables row to a FragmentReadPair
        """
        if lazy:
            left_read = LazyFragmentRead(row, self, side="left")
            right_read = LazyFragmentRead(row, self, side="right")
        else:
            fragment1 = GenomicRegion(start=row['left_fragment_start'],
                                      end=row['left_fragment_end'],
                                      chromosome=self._ix_to_chromosome[row['left_fragment_chromosome']],
                                      ix=row['left_fragment'])
            fragment2 = GenomicRegion(start=row['right_fragment_start'],
                                      end=row['right_fragment_end'],
                                      chromosome=self._ix_to_chromosome[row['right_fragment_chromosome']],
                                      ix=row['right_fragment'])

            left_read = FragmentRead(fragment1, position=row['left_read_position'],
                                     strand=row['left_read_strand'], qname_ix=row['left_read_qname_ix'])
            right_read = FragmentRead(fragment2, position=row['right_read_position'],
                                      strand=row['right_read_strand'], qname_ix=row['right_read_qname_ix'])

        return FragmentReadPair(left_read=left_read, right_read=right_read, ix=row['ix'])

    @lru_cache(maxsize=10)
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

            with RareUpdateProgressBar(max_value=len(self)) as pb:
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
                    x.append(int(mids/counter))
                    inward_ratios.append(inwards/same_counter)
                    outward_ratios.append(outwards/same_counter)
                    bin_sizes.append(counter)
                    same_counter = 0
                    counter = 0
                    mids = 0
                    outwards = 0
                    inwards = 0
                    same = 0
            return list(map(np.array, [x, inward_ratios, outward_ratios, bin_sizes]))

        gaps, types = _init_gaps_and_types()
        # sort data
        gaps, types = _sort_data(gaps, types)
        # best guess for number of data points
        sampling = _guess_sampling(sampling)
        # calculate ratios
        return _calculate_ratios(gaps, types, sampling)

    def filter(self, pair_filter, queue=False, log_progress=False):
        """
        Filter read pairs in this object by using a
        :class:`~FragmentMappedReadPairFilter`.

        :param pair_filter: Class implementing :class:`~FragmentMappedReadPairFilter`.
                            Must override valid_read method, ideally sets mask parameter
                            during initialization.
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        :param log_progress: If true, process iterating through all pairs
                             will be continuously reported.
        """
        pair_filter.set_pairs_object(self)
        if not queue:
            self._pairs.filter(pair_filter, _logging=log_progress)
        else:
            self._pairs.queue_filter(pair_filter)

    def run_queued_filters(self, log_progress=False):
        """
        Run queued filters.

        :param log_progress: If true, process iterating through all reads
                             will be continuously reported.
        """
        self._pairs.run_queued_filters(_logging=log_progress)

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
            p = (obs+exp) / (n*2)
            return abs((p_exp-p_obs) / np.sqrt(p*(1-p) * (2/n)))
        ratios = np.clip(ratios, 0.0, 1.0)
        z_scores = np.array([x_prop(r, expected_ratio, b) for r, b in zip(ratios, sample_sizes)])
        which_valid = z_scores < 1.96
        which_valid_indices = np.argwhere(which_valid).flatten()
        if len(which_valid_indices) > 0:
            return int(dists[which_valid_indices[0]])
        return None

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
        pcr_duplicate_filter = PCRDuplicateFilter(pairs=self, threshold=threshold, mask=mask)
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

    def pairs(self, lazy=False, excluded_filters=()):
        """
        Iterate over unfiltered fragment-mapped read pairs.
        """
        excluded_masks = self.get_binary_mask_from_masks(excluded_filters)
        it = self._pairs.iterrows(excluded_masks=excluded_masks)
        return (self._pair_from_row(i, lazy=lazy) for i in it)
    
    def __getitem__(self, key):
        """
        Get read pairs directly using int or slice as key.
        """
        if isinstance(key, int):
            row = self._pairs[key]
            return self._pair_from_row(row)
        elif isinstance(key, slice):
            pairs = []
            for row in self._pairs[key]:
                pairs.append(self._pair_from_row(row))
            return pairs
    
    def __len__(self):
        """
        Get the number of read pairs in this object.
        """
        return len(self._pairs)


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
        return min(abs(self.position-self.fragment.start),
                   abs(self.position-self.fragment.end))

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
            return self.row[self.side + "_fragment"]
        return self.static_ix


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
        return "%s -- %s" % (left_repr, right_repr)

    
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
        super(PCRDuplicateFilter, self).__init__(mask=mask)
        self.threshold = threshold
        self.pairs = pairs
        # In order for sorted iteration to work, column needs to be indexed
        try:
            self.pairs._pairs.cols.left_read_position.create_csindex()
            index_existed = False
        except ValueError: # Index already exists
            index_existed = True
        # Using itersorted from Table class, since MaskedTable.itersorted only yields unmasked entries
        all_iter = super(MaskedTable, self.pairs._pairs).itersorted(sortby="left_read_position")
        cur_pos = {}
        cur_duplicates = {}
        self.duplicates_set = set()
        duplicate_stats = defaultdict(int)
        for p in all_iter:
            pair = self.pairs._pair_from_row(p, lazy=True)
            chrm = (pair.left.fragment.chromosome, pair.right.fragment.chromosome)
            if cur_pos.get(chrm) is None:
                cur_pos[chrm] = (pair.left.position, pair.right.position)
                cur_duplicates[chrm] = 1
                continue
            if (abs(pair.left.position - cur_pos[chrm][0]) <= threshold and
                abs(pair.right.position - cur_pos[chrm][1]) <= threshold):
                self.duplicates_set.add(pair.ix)
                cur_duplicates[chrm] += 1
                continue
            if cur_duplicates[chrm] > 1:
                duplicate_stats[cur_duplicates[chrm]] += 1
            cur_pos[chrm] = (pair.left.position, pair.right.position)
            cur_duplicates[chrm] = 1
        if not index_existed:
            self.pairs._pairs.cols.left_read_position.remove_index()
        n_dups = len(self.duplicates_set)
        percent_dups = 1.*n_dups/self.pairs._pairs._original_len()
        logger.info("PCR duplicate stats: " +
                    "{} ({:.1%}) of pairs marked as duplicate. ".format(n_dups, percent_dups) +
                    " (multiplicity:occurances) " +
                    " ".join("{}:{}".format(k, v) for k, v in duplicate_stats.items()))

    def valid_pair(self, pair):
        """
        Check if a pair is duplicated.
        """
        if pair.ix in self.duplicates_set:
            return False
        return True


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

    def __init__(self, maximum_distance=500, mask=None):
        super(ReDistanceFilter, self).__init__(mask=mask)
        self.maximum_distance = maximum_distance
    
    def valid_pair(self, pair):
        """
        Check if any read is >maximum_distance away from RE site.
        """
        for read in [pair.left, pair.right]:
            if (read.position - read.fragment.start > self.maximum_distance
                    and read.fragment.end - read.position > self.maximum_distance):
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
