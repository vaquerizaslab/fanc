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
from kaic.tools.files import is_sambam_file
from kaic.data.general import Maskable, MetaContainer, MaskFilter, MaskedTable, FileBased
import tempfile
import os
import logging
from tables.exceptions import NoSuchNodeError
from abc import abstractmethod, ABCMeta
from bisect import bisect_right
from kaic.tools.general import bit_flags_from_int
from kaic.data.genomic import RegionsTable, GenomicRegion, LazyGenomicRegion
import subprocess
import msgpack as pickle
import numpy as np


class Reads(Maskable, MetaContainer, FileBased):
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

    def __init__(self, sambam_file=None, file_name=None,
                 qname_length=60, seq_length=200, read_only=False,
                 _group_name='reads', mapper=None):
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
        :param qname_length: If it is for some reason not wanted to load reads
                             during initialization (for example, when building
                             the Reads object from scratch using the add_read
                             method) this parameter can be used to set the
                             expected length of read names. If
                             the respective strings are longer than the field
                             size specified here, they will be truncated.
        :param seq_length: If it is for some reason not wanted to load reads
                           during initialization (for example, when building
                           the Reads object from scratch using the add_read
                           method) this parameter can be used to set the
                           expected length of sequences. If
                           the respective strings are longer than the field
                           size specified here, they will be truncated.
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

        FileBased.__init__(self, file_name, read_only=read_only)
        Maskable.__init__(self, self.file)
        MetaContainer.__init__(self, self.file)

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
            self._reads = self._file_group.get_node('main')
            self._row_counter['reads'] = len(self._reads)
        except NoSuchNodeError:
            self._reads = MaskedTable(self._file_group, 'main', Reads.ReadsDefinition)

        # Header attribute
        try:
            self._header = self._reads._v_attrs.header
        except AttributeError:
            logging.warn("No header attributes found in existing table")
            self._header = None

        # Reference names
        try:
            self._ref = self._reads._v_attrs.ref
        except AttributeError:
            logging.warn("No ref attributes found in existing table")
            self._ref = None

        # Qname table
        try:
            self._qname = self._file_group.get_node('qname')
            self._row_counter['qname'] = len(self._qname)
        except NoSuchNodeError:
            self._qname = None

        # Cigar table
        try:
            self._cigar = self._file_group.get_node('cigar')
            self._row_counter['cigar'] = len(self._cigar)
        except NoSuchNodeError:
            self._cigar = None

        # Seq table
        try:
            self._seq = self._file_group.get_node('seq')
            self._row_counter['seq'] = len(self._seq)
        except NoSuchNodeError:
            self._seq = None

        # Qual table
        try:
            self._qual = self._file_group.get_node('qual')
            self._row_counter['qual'] = len(self._qual)
        except NoSuchNodeError:
            self._qual = None

        # Tags table
        try:
            self._tags = self._file_group.get_node('tags')
            self._row_counter['tags'] = len(self._tags)
        except NoSuchNodeError:
            self._tags = None

        # load reads
        if sambam_file and is_sambam_file(sambam_file):
            self.load(sambam_file, ignore_duplicates=True)

        if mapper:
            self._mapper = mapper
        else:
            try:
                self._mapper = self.header['PG'][0]['ID']
            except (KeyError, AttributeError, TypeError):
                self._mapper = None
                logging.warn('Could not auto-detect mapping program from SAM header')

    @staticmethod
    def sambam_size(sambam):
        """
        Determine number of reads in SAM/BAM file.

        :param sambam: SAM/BAM file. Can be a string (path to file) or a
                       pysam AlignmentFile.
        :return: (int) read count
        """
        if type(sambam) == str:
            sambam = pysam.AlignmentFile(sambam, 'rb')  # @UndefinedVariable

        count = sum(1 for _ in iter(sambam))
        sambam.close()
        return count

    def close(self):
        """
        Close the file backing this object.
        """
        self.file.close()

    def load(self, sambam, ignore_duplicates=True, is_sorted=False,
             store_qname=True, store_cigar=True,
             store_seq=True, store_tags=True,
             store_qual=True, sample_fields=False):
        """
        Load mapped reads from SAM/BAM file.

        Will sort mapped reads by name if not explicitly deactivated.
        By default, successive reads with duplicate names will be ignored.

        :param sambam: A string that describes the path to the SAM/BAM file
                       to be loaded or a pysam AlignmentFile.
        :param ignore_duplicates: If True (default) will not load reads that
                                  have the same name as the last loaded read.
                                  Will load all reads if False.
        :param is_sorted: If False will sort the file using samtools (if
                          available, otherwise will attempt to sort using
                          pysam.sort). If True assumes that the file is
                          already sorted.
        """

        self.log_info("Loading SAM/BAM file")
        # get file names
        try:
            file_name = sambam.filename
        except AttributeError:
            file_name = sambam

        # sort files if required
        if not is_sorted:
            self.log_info("Sorting...")
            tmp_file = tempfile.NamedTemporaryFile(delete=False, suffix=".bam")
            tmp_file.close()
            logging.info(file_name)
            logging.info(os.path.splitext(tmp_file.name)[0])
            try:
                subprocess.call(["samtools", "sort", "-n", file_name,
                                 os.path.splitext(tmp_file.name)[0]])
            except OSError:
                pysam.sort('-n', file_name, os.path.splitext(tmp_file.name)[0])
            self.log_info("Done. Reading sorted BAM file...")
            sambam = pysam.AlignmentFile(tmp_file.name, 'rb')
            self.log_info("Done...")
        else:
            sambam = pysam.AlignmentFile(file_name, 'rb')

        if sample_fields:
            qname_length, seq_length, cigar_length, tags_length = Reads.determine_field_sizes(file_name, 10000)
            qname_length *= 2
            seq_length *= 2
            cigar_length *= 2
            tags_length *= 2
        else:
            qname_length, seq_length, cigar_length, tags_length = Reads.determine_field_sizes(file_name, None)

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

        # references
        self._reads._v_attrs.ref = sambam.references
        self._ref = sambam.references

        self.log_info("Loading mapped reads...")
        last_name = ""
        for i, read in enumerate(sambam):
            if i % 10000 == 0:
                self.log_info("%d" % i, save=False)
            if ignore_duplicates and read.qname == last_name:
                continue
            self.add_read(read, flush=False, store_cigar=store_cigar,
                          store_seq=store_seq, store_qual=store_qual,
                          store_qname=store_qname, store_tags=store_tags)
            last_name = read.qname
        self.flush()

        self.log_info("Done.")

    @staticmethod
    def determine_field_sizes(sambam, sample_size=10000):
        """
        Determine the sizes of relevant fields in a SAM/BAM file.

        :param sambam: A string that describes the path to the SAM/BAM file
                       to be loaded or a pysam AlignmentFile.
        :param sample_size: Number of lines to sample to determine field
                            sizes.
        :return: qname length, sequence length
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
            qname_length = max(qname_length, len(r.qname))
            seq_length = max(seq_length, len(r.seq))
            cigar = r.cigarstring
            if cigar is not None:
                cigar_dump = pickle.dumps(cigar)
                cigar_length = max(cigar_length, len(cigar_dump))
            tags = r.tags
            if tags is not None:
                tags_dump = pickle.dumps(tags)
                tags_length = max(tags_length, len(tags_dump))
            if sample_size is not None and i >= sample_size:
                break

            if i % 100000 == 0:
                logging.info(i)
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
        if store_qname:
            self._qname.append([read.qname])
            reads_row['qname'] = self._row_counter['qname']
            self._row_counter['qname'] += 1

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

    def flush(self):
        """
        Write the latest changes to this object to file.
        """
        self._reads.flush(update_index=True)
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
                    tags=tags, reference_id=ref_ix)

    def reads(self, lazy=False):
        """
        Iterate over _reads table and convert each result to Read.

        :param lazy: Lazily load read properties (only works inside loop!)
        :return: ReadsIter that iterates over visible reads
        """
        this = self

        class ReadsIter:
            def __init__(self):
                self.iter = iter(this._reads)

            def __iter__(self):
                return self

            def next(self):
                row = self.iter.next()
                return this._row2read(row, lazy=lazy)
        return ReadsIter()

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

    def filter_non_unique(self, strict=True, cutoff=0.5, queue=False):
        """
        Convenience function that applies a UniquenessFilter.
        The actual algorithm and rationale used for filtering will depend on the
        internal _mapper attribute.

        :param cutoff: Ratio of the secondary to the primary alignment score. Smaller
                       values mean that the next best secondary alignment is of substantially
                       lower quality than the primary one, and that the latter can be considered
                       as unique. Used only if the reads have been aligned with bwa-mem.
        :param strict: If True will filter if XS tag is present. If False,
                       will filter only when XS tag is not 0.
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        """
        mask = self.add_mask_description('uniqueness', 'Mask read pairs that do not map uniquely (according to XS tag)')
        if self._mapper == 'bwa':
            uniqueness_filter = BwaMemUniquenessFilter(cutoff, mask)
        else:
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
                                  tags=read.tags, masks=masks)
        return MaskedReadsIter()


class Read(object):
    """
    Object representing a mapped sequencing read.

    This is the main class used by a Reads object to represent mapped
    read information.
    """

    def __init__(self, qname="", flag=0, ref="",
                 pos=0, mapq=0, cigar="", rnext=0,
                 pnext=0, tlen=0, seq="", qual="",
                 tags={}, reference_id=None):
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
    def __init__(self, qname="", flag=0, ref="",
                 pos=0, mapq=0, cigar="", rnext=0,
                 pnext=0, tlen=0, seq="", qual="",
                 tags={}, reference_id=None, masks=None):
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
                                         tags=tags, reference_id=reference_id)
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


#
# Filters Reads
#
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

    __metaclass__ = ABCMeta

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
        if xs_tag is not None or (not self.strict and xs_tag == 0):
            return False
        return True


class BwaMemUniquenessFilter(ReadFilter):
    """
    Filters `bwa mem` generated alignements based on whether they are unique or not.
    The presence of a non-zero XS tag does not mean a read is a multi-mapping one.
    Instead, we make sure that the ratio XS/AS is inferior to a certain threshold.
    """
    def __init__(self, cutoff=0.5, mask=None):
        """
        :param cutoff: Ratio of the secondary to the primary alignment score. Smaller
                       values mean that the next best secondary alignment is of substantially
                       lower quality than the primary one, and that the latter can be considered
                       as unique
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered reads.
        """
        super(BwaMemUniquenessFilter, self).__init__(mask)
        self.cutoff = cutoff

    def valid_read(self, read):
        """
        Check if a read has a high alignment score.
        """
        alignment_score = read.get_tag('AS')
        nextbest_score = read.get_tag('XS')
        if alignment_score:
            return float(nextbest_score) / alignment_score <= self.cutoff
        return False


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


class FragmentMappedReadPairs(Maskable, MetaContainer, RegionsTable, FileBased):
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

    class FragmentMappedReadDescription(t.IsDescription):
        """
        Needed by PyTables to build mapped read table.
        """
        ix = t.Int32Col(pos=0)
        position = t.Int64Col(pos=2)
        strand = t.Int8Col(pos=3)

    class FragmentsMappedReadPairDescription(t.IsDescription):
        """
        Needed by PyTables to build read pairs table.
        """
        ix = t.Int32Col(pos=0)
        left_read = t.Int32Col(pos=1)
        left_fragment = t.Int32Col(pos=2, dflt=-1)
        right_read = t.Int32Col(pos=3)
        right_fragment = t.Int32Col(pos=4, dflt=-1)

    class FragmentsMappedReadSingleDescription(t.IsDescription):
        """
        Needed by PyTables to build single reads table.
        """
        ix = t.Int32Col(pos=0)
        read = t.Int32Col(pos=1)
        fragment = t.Int32Col(pos=2, dflt=-1)

    def __init__(self, file_name=None,
                 read_only=False,
                 group_name='fragment_map',
                 table_name_fragments='fragments'):
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

        FileBased.__init__(self, file_name, read_only=read_only)
        RegionsTable.__init__(self, file_name=self.file, _table_name_regions=table_name_fragments)

        # generate tables from inherited classes
        Maskable.__init__(self, self.file)
        MetaContainer.__init__(self, self.file)

        # try to retrieve existing table
        try:
            self._reads = self.file.get_node('/' + group_name + '/mapped_reads')
            self._pairs = self.file.get_node('/' + group_name + '/mapped_read_pairs')
            self._single = self.file.get_node('/' + group_name + '/mapped_read_single')
            self._read_count = len(self._reads)
            self._pair_count = self._pairs._original_len()
            self._single_count = self._single._original_len()
        # or build table from scratch
        except NoSuchNodeError:
            # create group
            group = self.file.create_group("/", group_name, 'Mapped read pairs group',
                                           filters=t.Filters(complib="blosc",
                                                             complevel=2, shuffle=True))
            # create main tables
            self._reads = t.Table(group, 'mapped_reads',
                                  FragmentMappedReadPairs.FragmentMappedReadDescription)

            self._pairs = MaskedTable(group, 'mapped_read_pairs',
                                      FragmentMappedReadPairs.FragmentsMappedReadPairDescription)

            self._single = MaskedTable(group, 'mapped_read_single',
                                       FragmentMappedReadPairs.FragmentsMappedReadSingleDescription)
            self._read_count = 0
            self._pair_count = 0
            self._single_count = 0

        try:
            self._pairs.cols.left_fragment.create_csindex()
        except ValueError:
            # Index exists, no problem!
            pass

    def load(self, reads1, reads2, regions=None, ignore_duplicates=True, _in_memory_index=True):
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
        if regions is not None:
            self.log_info("Adding regions...")
            self.add_regions(regions)
            self.log_info("Done.")

        # generate index for fragments
        fragment_ixs = None
        fragment_ends = None
        if _in_memory_index:
            fragment_ixs = {}
            fragment_ends = {}
            for region in self.regions():
                if region.chromosome not in fragment_ends:
                    fragment_ends[region.chromosome] = []
                    fragment_ixs[region.chromosome] = []
                fragment_ixs[region.chromosome].append(region.ix)
                fragment_ends[region.chromosome].append(region.end)

        if isinstance(reads1, str):
            self.log_info("Loading reads 1")
            reads1 = Reads(sambam_file=reads1)

        if isinstance(reads2, str):
            self.log_info("Loading reads 2")
            reads2 = Reads(sambam_file=reads2)

        iter1 = iter(reads1)
        iter2 = iter(reads2)

        def get_next_read(iterator):
            try:
                r = iterator.next()
                return r
            except StopIteration:
                return None

        # add and map reads
        i = 0
        last_r1_name = ''
        last_r2_name = ''
        r1 = get_next_read(iter1)
        r2 = get_next_read(iter2)
        r1_count = 0
        r2_count = 0
        while r1 is not None and r2 is not None:
            c = cmp_natural(r1.qname, r2.qname)

            i += 1
            if r1.qname == last_r1_name:
                if not ignore_duplicates:
                    raise ValueError("Duplicate left read QNAME %s" % r1.qname)
                r1 = get_next_read(iter1)
                r1_count += 1
            elif r2.qname == last_r2_name:
                if not ignore_duplicates:
                    raise ValueError("Duplicate right read QNAME %s" % r2.qname)
                r2 = get_next_read(iter2)
                r2_count += 1
            elif c == 0:
                self.add_read_pair(r1, r2, flush=False, _fragment_ends=fragment_ends, _fragment_ixs=fragment_ixs)
                last_r1_name = r1.qname
                last_r2_name = r2.qname
                r1 = get_next_read(iter1)
                r2 = get_next_read(iter2)
                r1_count += 1
                r2_count += 1
            elif c < 0:
                self.add_read_single(r1, flush=False, _fragment_ends=fragment_ends, _fragment_ixs=fragment_ixs)
                last_r1_name = r1.qname
                r1 = get_next_read(iter1)
                r1_count += 1
            else:
                self.add_read_single(r2, flush=False, _fragment_ends=fragment_ends, _fragment_ixs=fragment_ixs)
                last_r2_name = r2.qname
                r2 = get_next_read(iter2)
                r2_count += 1

            if i % 100000 == 0:
                logging.info("%d reads processed" % i)

        # add remaining unpaired reads
        while r1 is not None:
            if r1.qname == last_r1_name:
                if not ignore_duplicates:
                    raise ValueError("Duplicate left read QNAME %s" % r1.qname)
            else:
                self.add_read_single(r1, flush=False, _fragment_ends=fragment_ends, _fragment_ixs=fragment_ixs)
            last_r1_name = r1.qname
            r1 = get_next_read(iter1)
            r1_count += 1

        while r2 is not None:
            if r2.qname == last_r2_name:
                if not ignore_duplicates:
                    raise ValueError("Duplicate right read QNAME %s" % r2.qname)
            else:
                self.add_read_single(r2, flush=False, _fragment_ends=fragment_ends, _fragment_ixs=fragment_ixs)
            last_r2_name = r2.qname
            r2 = get_next_read(iter2)
            r2_count += 1

        logging.info('Counts: R1 %d R2 %d' % (r1_count, r2_count))

        self._reads.flush()
        self._pairs.flush(update_index=True)
        self._single.flush(update_index=True)

    def add_read_pair(self, read1, read2, flush=True, _fragment_ends=None, _fragment_ixs=None):
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
        :param _fragment_ixs: (Internal) See previous argument, but provide
                              fragment indices instead of end positions.
        """
        ix1 = self._add_read(read1, flush=flush)
        ix2 = self._add_read(read2, flush=flush)
        fragment_ix1 = self._find_fragment_ix(read1.ref, read1.pos,
                                              _fragment_ends=_fragment_ends, _fragment_ixs=_fragment_ixs)
        fragment_ix2 = self._find_fragment_ix(read2.ref, read2.pos,
                                              _fragment_ends=_fragment_ends, _fragment_ixs=_fragment_ixs)

        # both must be integer if successfully mapped
        if fragment_ix1 is not None and fragment_ix2 is not None:
            row = self._pairs.row
            row['ix'] = self._pair_count
            if fragment_ix1 <= fragment_ix2:
                row['left_read'] = ix1
                row['right_read'] = ix2
                row['left_fragment'] = fragment_ix1
                row['right_fragment'] = fragment_ix2
            else:
                row['left_read'] = ix2
                row['right_read'] = ix1
                row['left_fragment'] = fragment_ix2
                row['right_fragment'] = fragment_ix1
            row.append()
            self._pair_count += 1

            if flush:
                self._pairs.flush(update_index=True)

    def add_read_single(self, read, flush=True, _fragment_ends=None, _fragment_ixs=None):
        ix = self._add_read(read, flush=flush)
        fragment_ix = self._find_fragment_ix(read.ref, read.pos, _fragment_ends=_fragment_ends, _fragment_ixs=_fragment_ixs)
        if fragment_ix is not None:
            row = self._single.row
            row['ix'] = self._single_count
            row['fragment'] = fragment_ix
            row['read'] = ix
            row.append()

            if flush:
                self._single.flush(update_index=True)
            self._single_count += 1

    def _find_fragment_ix(self, chromosome, position, _fragment_ends=None, _fragment_ixs=None):
        """
        Find the index of a fragment by genomic coordinate and chromosome name.
        """
        # binary search for fragment
        fragment_ix = None
        if _fragment_ends is not None and _fragment_ixs is not None:
            try:
                pos_ix = bisect_right(_fragment_ends[chromosome], position)
                fragment_ix = _fragment_ixs[chromosome][pos_ix]
            except KeyError:
                # potentially keep a record of unmatched chromosome names
                pass
        else:
            for row in self._regions.where(
                            "(start <= %d) & (end >= %d) & (chromosome == '%s')" % (position, position, chromosome)):
                fragment_ix = row['ix']

        return fragment_ix

    def _add_read(self, read, flush=True):
        """
        Add position and strand information of a read to reads table.
        """
        ix = self._read_count
        row = self._reads.row
        row['ix'] = ix
        row['position'] = read.pos
        if hasattr(read, 'strand') and read.strand is not None:
            row['strand'] = read.strand
        else:
            bit_flags = bit_flags_from_int(read.flag)
            if 4 in bit_flags:
                row['strand'] = -1
            else:
                row['strand'] = 1
        row.append()

        self._read_count += 1

        if flush:
            self._reads.flush()

        return ix

    def _pair_from_row(self, row, lazy=False):
        """
        Convert a pytables row to a FragmentReadPair
        """
        if lazy:
            left_read = LazyFragmentRead(row, self, side="left")
            right_read = LazyFragmentRead(row, self, side="right")
            return FragmentReadPair(left_read=left_read, right_read=right_read)

        ix1 = row['left_read']
        ix2 = row['right_read']
        fragment_ix1 = row['left_fragment']
        fragment_ix2 = row['right_fragment']

        read1_row = self._reads[ix1]
        fragment1_row = self._regions[fragment_ix1]
        fragment1 = GenomicRegion(fragment1_row['start'], fragment1_row['end'], fragment1_row['chromosome'])
        left_read = FragmentRead(fragment1, position=read1_row['position'], strand=read1_row['strand'])

        read2_row = self._reads[ix2]
        fragment2_row = self._regions[fragment_ix2]
        fragment2 = GenomicRegion(fragment2_row['start'], fragment2_row['end'], fragment2_row['chromosome'])
        right_read = FragmentRead(fragment2, position=read2_row['position'], strand=read2_row['strand'])

        return FragmentReadPair(left_read=left_read, right_read=right_read)

    def get_error_structure(self, data_points=None, skip_self_ligations=True):

        """
        Compute the ligation error structure of this data set.

        :param data_points: Number of data points to average per point
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

        def _log_done(i):
            last_percent = -1
            if i % int(l/20) == 0:
                percent = int(i/int(l/20))
                if percent != last_percent:
                    logging.info("Ligation error structure: {}%".format(percent*5))
                    last_percent = percent

        def _init_gaps_and_types():
            same_count = 0
            inward_count = 0
            outward_count = 0
            gaps = []
            types = []
            for i, pair in enumerate(self):
                _log_done(i)
                if pair.is_same_fragment() and skip_self_ligations:
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
            logging.info("Pairs: %d" % l)
            logging.info("Same: %d" % same_count)
            logging.info("Inward: %d" % inward_count)
            logging.info("Outward: %d" % outward_count)
            return gaps, types

        def _sort_data(gaps, types):
            points = zip(gaps, types)
            sorted_points = sorted(points)
            gaps = [point[0] for point in sorted_points]
            types = [point[1] for point in sorted_points]
            return gaps, types

        def _guess_datapoints(data_points):
            if data_points is None:
                data_points = max(100, int(l * 0.0025))
            logging.info("Number of data points averaged per point in plot: {}".format(data_points))
            return data_points

        def _calculate_ratios(gaps, types, data_points):
            x = []
            inward_ratios = []
            outward_ratios = []
            counter = 0
            same_counter = 0
            mids = 0
            outwards = 0
            inwards = 0
            same = 0
            for typ, gap in zip(types, gaps):
                mids += gap
                if typ == type_same:
                    same += 1
                    same_counter += 1
                elif typ == type_inward:
                    inwards += 1
                else:
                    outwards += 1
                counter += 1
                if same_counter > data_points:
                    x.append(mids/counter)
                    inward_ratios.append(inwards/same)
                    outward_ratios.append(outwards/same)
                    same_counter = 0
                    counter = 0
                    mids = 0
                    outwards = 0
                    inwards = 0
                    same = 0
            return x, inward_ratios, outward_ratios

        gaps, types = _init_gaps_and_types()
        # sort data
        gaps, types = _sort_data(gaps, types)
        # best guess for number of data points
        data_points = _guess_datapoints(data_points)
        # calculate ratios
        x, inward_ratios, outward_ratios = _calculate_ratios(gaps, types, data_points)
        return x, inward_ratios, outward_ratios

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

    def _auto_dist(self, dists, ratios, threshold_ratio=0.1, threshold_std=0.1, window=3):
        """
        Function that attempts to infer sane distances for filtering inward
        and outward read pairs

        :param dists: List of distances in bp.
        :param ratios: List of ratios
        :param threshold_ratio: Threshold below which the 1+log2(ratio) must fall
                                in order to infer the corresponding distance
        :param threshold_std: Threshold below which the standard deviation of 1+log2(ratio)
                              must fall in order to infer the corresponding distance
        :param window: Window for the rolling standard deviations
        """
        def rolling_window(a, window):
            shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
            strides = a.strides + (a.strides[-1],)
            r = np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)
            return np.insert(r, 0, [r[0]]*(window-1), axis=0)

        ratios = np.array([1+np.log2(x) for x in ratios])
        stds = np.std(rolling_window(ratios, window), -1)
        ok_threshold_ratios = abs(ratios) <= threshold_ratio
        ok_threshold_stds = stds <= threshold_std
        ok_values = ok_threshold_ratios * ok_threshold_stds
        ok_indices = np.argwhere(ok_values).flatten()
        if len(ok_indices) > 0:
            return dists[ok_indices[0]]
        return None

    def filter_inward(self, minimum_distance=None, queue=False, threshold_ratio=0.1, threshold_std=0.1, window=3):
        """
        Convenience function that applies an :class:`~InwardPairsFilter`.

        :param minimum_distance: Minimum distance inward-facing read
                                 pairs must have to pass this filter
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        :param threshold_ratio: Threshold below which the 1+log2(ratio) must fall
                                in order to infer the corresponding distance
        :param threshold_std: Threshold below which the standard deviation of 1+log2(ratio)
                              must fall in order to infer the corresponding distance
        :param window: Window for the rolling standard deviations
        """
        if minimum_distance is None:
            dists, inward_ratios, _ = self.get_error_structure()
            minimum_distance = self._auto_dist(dists, inward_ratios, threshold_ratio, threshold_std, window)
        if minimum_distance:
            mask = self.add_mask_description('inward',
                                             'Mask read pairs that are inward facing and <%dbp apart' % minimum_distance)
            inward_filter = InwardPairsFilter(minimum_distance=minimum_distance, mask=mask)
            self.filter(inward_filter, queue)
        else:
            raise Exception('Could not automatically detect a sane distance threshold for filtering inward reads')

    def filter_outward(self, minimum_distance=None, queue=False, threshold_ratio=0.1, threshold_std=0.1, window=3):
        """
        Convenience function that applies an :class:`~OutwardPairsFilter`.

        :param minimum_distance: Minimum distance outward-facing read
                                 pairs must have to pass this filter
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        :param threshold_ratio: Threshold below which the 1+log2(ratio) must fall
                                in order to infer the corresponding distance
        :param threshold_std: Threshold below which the standard deviation of 1+log2(ratio)
                              must fall in order to infer the corresponding distance
        :param window: Window for the rolling standard deviations
        """
        if minimum_distance is None:
            dists, _, outward_ratios = self.get_error_structure()
            minimum_distance = self._auto_dist(dists, outward_ratios, threshold_ratio, threshold_std, window)
        if minimum_distance:
            mask = self.add_mask_description('outward',
                                             'Mask read pairs that are outward facing and <%dbp apart' % minimum_distance)
            outward_filter = OutwardPairsFilter(minimum_distance=minimum_distance, mask=mask)
            self.filter(outward_filter, queue)
        else:
            raise Exception('Could not automatically detect a sane distance threshold for filtering outward reads')
    
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
    
    def __iter__(self):
        """
        Iterate over unfiltered fragment-mapped read pairs.
        """
        return self.pairs(lazy=False)

    def pairs(self, lazy=False):
        """
        Iterate over unfiltered fragment-mapped read pairs.
        """
        this = self

        class FragmentMappedReadPairIter:
            def __init__(self):
                self.iter = iter(this._pairs)

            def __iter__(self):
                return self

            def next(self):
                return this._pair_from_row(self.iter.next(), lazy=lazy)

            def __len__(self):
                return len(this._pairs)

        return FragmentMappedReadPairIter()
    
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
    def __init__(self, fragment=None, position=None, strand=0):
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
    def __init__(self, row, parent, side='left'):
        self.row = row
        self.parent = parent
        self.side = side
        self.read_row = None
        self.fragment_row = None

    def _update_read_row(self):
        if self.read_row is None:
            ix = self.row[self.side + '_read']
            self.read_row = self.parent._reads[ix]

    def _update_fragment_row(self):
        if self.fragment_row is None:
            ix = self.row[self.side + '_fragment']
            self.fragment_row = self.parent._regions[ix]

    @property
    def position(self):
        self._update_read_row()
        return self.read_row["position"]

    @property
    def strand(self):
        self._update_read_row()
        return self.read_row["strand"]

    @property
    def fragment(self):
        self._update_fragment_row()
        return LazyGenomicRegion(self.fragment_row)


class FragmentReadPair(object):
    """
    Container for two paired :class:`~FragmentRead` objects.
    """
    def __init__(self, left_read, right_read):
        self.left = left_read
        self.right = right_read
    
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

    
class FragmentMappedReadPairFilter(MaskFilter):
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

    __metaclass__ = ABCMeta
    
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


def cmp_natural(string1, string2):
    """
    Compare two strings the same way that samtools does.
    """
    def is_digit(char):
        try:
            int(char)
            return True
        except (ValueError, TypeError):
            return False

    class CharIterator(object):
        def __init__(self, string):
            self.string = string
            self.current = 0
        
        def next(self):
            try:
                char = self.string[self.current]
                self.current += 1
                return char
            except IndexError:
                return None
    
        def current_plus(self, n):
            try:
                char = self.string[self.current+n]
                return char
            except IndexError:
                return None
    
    char_iter1 = CharIterator(string1)
    char_iter2 = CharIterator(string2)
    
    c1 = char_iter1.next()
    c2 = char_iter2.next()
    
    while c1 and c2:
        if is_digit(c1) and is_digit(c2):
            # ignore leading zeros
            while c1 == '0':
                c1 = char_iter1.next()
            while c2 == '0':
                c2 = char_iter2.next()
            
            # skip through identical digits
            while is_digit(c1) and is_digit(c2) and c1 == c2:
                c1 = char_iter1.next()
                c2 = char_iter2.next()

            if is_digit(c1) and is_digit(c2):
                # compare numbers at this point
                n = 0
                while is_digit(char_iter1.current_plus(n)) and is_digit(char_iter2.current_plus(n)):
                    n += 1
                if is_digit(char_iter1.current_plus(n)):
                    return 1
                if is_digit(char_iter2.current_plus(n)):
                    return -1
                if c1 > c2:
                    return 1
                return -1
            elif is_digit(c1):
                return 1
            elif is_digit(c2):
                return -1
            elif char_iter1.current != char_iter2.current: # TODO double-check this block!
                if char_iter1.current > char_iter2.current:
                    return 1
                return -1
        else:
            if c1 != c2:
                if c1 > c2:
                    return 1
                return -1
            c1 = char_iter1.next()
            c2 = char_iter2.next()
    
    if char_iter1.current < len(string1):
        return 1
    if char_iter1.current < len(string2):
        return -1
    return 0
