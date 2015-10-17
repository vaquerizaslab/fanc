'''
Created on Jul 13, 2015

@author: kkruse1
'''

from __future__ import division
import tables as t
import pysam
from kaic.tools.files import create_or_open_pytables_file, random_name,\
    is_sambam_file
from kaic.data.general import Maskable, MetaContainer, MaskFilter, MaskedTable,\
    FileBased
import tempfile
import os
import logging
from tables.exceptions import NoSuchNodeError
from abc import abstractmethod, ABCMeta
from bisect import bisect_right
from kaic.tools.general import bit_flags_from_int
from kaic.data.genomic import RegionsTable, GenomicRegion
from matplotlib import pyplot as plt
import subprocess

        
class Reads(FileBased, Maskable, MetaContainer):
    def __init__(self, sambam_file=None, file_name=None, group_name='reads',
                 field_sizes={'qname': 60, 'sequence': 200}):
        if (sambam_file is not None
            and file_name is None
            and not is_sambam_file(sambam_file)):
            file_name = sambam_file
            sambam_file = None
        
        FileBased.__init__(self, file_name)
        Maskable.__init__(self, self.file)
        MetaContainer.__init__(self, self.file)
        
        # try to retrieve existing table
        try:
            main_table = self.file.get_node('/' + group_name + '/main')
            try:
                self._header = main_table._v_attrs.header
            except AttributeError:
                logging.warn("No header attributes found in existing table")
                self._header = None
            try:
                self._ref = main_table._v_attrs.ref
            except AttributeError:
                logging.warn("No ref attributes found in existing table")
                self._ref = None
            
            # tags
            tags = self.file.get_node('/' + group_name + '/tags')
            
            # cigar
            cigar = self.file.get_node('/' + group_name + '/cigar')
            
            
            
        # or build table from scratch
        except NoSuchNodeError:
            
            expected_length = 50000000
            if sambam_file is not None:
                logging.info("Determining field sizes")
                field_sizes = Reads.determine_field_sizes(sambam_file)
                expected_length = int(self._sambam_size(sambam_file)*1.5)
                
            reads_defininition = {
                'ix': t.Int32Col(pos=0),
                'qname': t.StringCol(field_sizes['qname'],pos=1),
                'flag': t.Int32Col(pos=2),
                'ref': t.Int32Col(pos=3),
                'pos': t.Int64Col(pos=4),
                'mapq': t.Int32Col(pos=5),
                'rnext': t.Int32Col(pos=6),
                'pnext': t.Int32Col(pos=7),
                'tlen': t.Int32Col(pos=8),
                'seq': t.StringCol(field_sizes['sequence'],pos=9),
                'qual': t.StringCol(field_sizes['sequence'],pos=10)
            }
        
            # create data structures
            logging.info("Creating data structures...")
            
            # create reads group
            group = self.file.create_group("/", group_name, 'Read pairs group',
                                           filters=t.Filters(complib="blosc",
                                                             complevel=2, shuffle=True))
            # create main table
            main_table = MaskedTable(group, 'main', reads_defininition,
                                     expectedrows=expected_length)
            # create tags vlarrays
            tags = self.file.create_vlarray(group, 'tags', t.ObjectAtom(),
                                            expectedrows=expected_length)
            
            # create cigar vlarrays
            cigar = self.file.create_vlarray(group, 'cigar', t.VLStringAtom(),
                                             expectedrows=expected_length)    
            
        
        self._reads = main_table
        self._tags = tags
        self._cigar = cigar

        # load reads
        if sambam_file and is_sambam_file(sambam_file):
            self.load(sambam_file, ignore_duplicates=True)
    
    def _sambam_size(self, sambam):
        if type(sambam) == str:
            sambam = pysam.AlignmentFile(sambam, 'rb')  # @UndefinedVariable
        
        count = sum(1 for _ in iter(sambam))
        sambam.close()
        return count
    
    def _get_row_counter(self):
        try:
            return self._reads._v_attrs.row_counter
        except AttributeError:
            self._reads._v_attrs.row_counter = 0
            return 0
    
    def _set_row_counter(self, value):
        self._reads._v_attrs.row_counter = value
    
    def close(self):
        self.file.close()    
    
    def load(self, sambam, ignore_duplicates=True, is_sorted=False):
        logging.info("Loading SAM/BAM file")
        # get file names
        try:
            file_name = sambam.filename
        except AttributeError:
            file_name = sambam
        
        # sort files if required
        if not is_sorted:
            logging.info("Sorting...")
            tmp_file = tempfile.NamedTemporaryFile(delete=False, suffix=".bam")
            tmp_file.close()
            logging.info(file_name)
            logging.info(os.path.splitext(tmp_file.name)[0])
            #pysam.sort('-n', file_name, os.path.splitext(tmp_file.name)[0])
            try:
                subprocess.call(["samtools", "sort", "-n", file_name,
                                 os.path.splitext(tmp_file.name)[0]])
            except:
                pysam.sort('-n', file_name, os.path.splitext(tmp_file.name)[0])
            logging.info("Done. Reading sorted BAM file...")
            sambam = pysam.AlignmentFile(tmp_file.name, 'rb')
            logging.info("Done...")
        else:
            sambam = pysam.AlignmentFile(file_name, 'rb')
        
        # header
        self._reads._v_attrs.header = sambam.header
        self._header = sambam.header
        
        # references
        self._reads._v_attrs.ref = sambam.references
        self._ref = sambam.references
        
        logging.info("Reading in mapped reads...")
        last_name = ""
        for i, read in enumerate(sambam):
            if i % 10000 == 0:
                logging.info("%d" % i)
            if ignore_duplicates and read.qname == last_name:
                continue
            self.add_read(read, flush=False)
            last_name = read.qname
        self.flush()
        
        logging.info("Done.")
        
    
    def add_read(self, read, flush=True):
        reads_row = self._reads.row
        # add main read info
        reads_row['qname'] = read.qname
        reads_row['flag'] = read.flag
        reads_row['ref'] = read.reference_id
        if read.pos >= 0:
            reads_row['pos'] = read.pos+1
        else:
            reads_row['pos'] = read.pos
        reads_row['mapq'] = read.mapq
        reads_row['rnext'] = read.rnext
        reads_row['pnext'] =  read.pnext
        reads_row['tlen'] = read.tlen
        reads_row['seq'] = read.seq
        reads_row['qual'] = read.qual
        reads_row['ix'] = self._get_row_counter()
        reads_row.append()
        
        # add string info
        self._tags.append(read.tags)
        if read.cigarstring is not None:
            self._cigar.append(read.cigarstring)
        else:
            self._cigar.append('')
        self._set_row_counter(self._get_row_counter()+1)
        
        if flush:
            self.flush()
    
    def flush(self):
        self._reads.flush(update_index=True)
        self._tags.flush()
        self._cigar.flush()
    
    @staticmethod
    def determine_field_sizes(sambam, sample_size=10000):
        if type(sambam) == str:
            sambam = pysam.AlignmentFile(sambam, 'rb')  # @UndefinedVariable
            
        qname_length = 0
        seq_length = 0
        i = 0
        for r in sambam:
            i += 1
            qname_length = max(qname_length,len(r.qname))
            seq_length = max(seq_length,len(r.seq))
            if sample_size is not None and i >= sample_size:
                break
            
            if i % 100000 == 0:
                logging.info(i)
        sambam.close()
        
        return { 'qname': qname_length,
                 'sequence': seq_length }
    
    @property
    def header(self):
        return self._header
    
    def _ix2ref(self, ix):
        if self._ref is None:
            raise RuntimeError("Chromosome reference for left read not present")
        return self._ref[ix]
    
    def _row2read(self, row):
        ix = row['ix']
        tags = self._tags[ix]
        cigar = self._cigar[ix]
        ref = self._ix2ref(row['ref'])
        
        return Read(qname=row['qname'], flag=row['flag'], ref=ref,
                 pos=row['pos'], mapq=row['mapq'], cigar=cigar, rnext=row['rnext'],
                 pnext=row['pnext'], tlen=row['tlen'], seq=row['seq'], qual=row['qual'],
                 tags=tags)
        
    def __iter__(self):
        this = self
        class ReadsIter:
            def __init__(self):
                self.iter = iter(this._reads)
                  
            def __iter__(self):
                return self
              
            def next(self):
                row = self.iter.next()
                return this._row2read(row)
        return ReadsIter()
    
    def __getitem__(self, key):
        return self._row2read(self._reads[key])
    
    def __len__(self):
        return len(self._reads)
    
    def where(self, query):
        reads = []
        for row in self._reads.where(query):
            reads.append(self._row2read(row))
        return reads
    
    def filter(self, read_filter, queue=False, log_progress=False):
        read_filter.set_reads_object(self)
        if not queue:
            self._reads.filter(read_filter, _logging=log_progress)
        else:
            self._reads.queue_filter(read_filter)
    
    def filter_quality(self, cutoff=30, queue=False):
        mask = self.add_mask_description('mapq', 'Mask read pairs with a mapping quality lower than %d' % cutoff)
        quality_filter = QualityFilter(cutoff, mask)
        self.filter(quality_filter, queue)
    
    def filter_unmapped(self, queue=False):
        mask = self.add_mask_description('unmapped', 'Mask read pairs that are unmapped')
        unmapped_filter = UnmappedFilter(mask)
        self.filter(unmapped_filter, queue)
            
    def filter_non_unique(self, strict=True, queue=False):
        mask = self.add_mask_description('uniqueness', 'Mask read pairs that do not map uniquely (according to XS tag)')
        uniqueness_filter = UniquenessFilter(strict, mask)
        self.filter(uniqueness_filter, queue)
    
    def run_queued_filters(self, log_progress=False):
        self._reads.run_queued_filters(_logging=log_progress)
    
    def filtered_reads(self):
        this = self
        class MaskedReadsIter:
            def __init__(self):
                self.iter = this._reads.masked_rows()
                  
            def __iter__(self):
                return self
              
            def next(self):
                row = self.iter.next()
                read = self._row2read(row)
                
                masks = this.get_masks(row[this._mask_field])
                
                return MaskedRead(qname=read.qname, flag=read.flag, ref=read.ref,
                                  pos=read.pos, mapq=read.mapq, cigar=read.cigar, rnext=read.rnext,
                                  pnext=read.pnext, tlen=read.tlen, seq=read.seq, qual=read.qual,
                                  tags=read.tags, masks=masks)
        return MaskedReadsIter()
        
        
class ReadPairs(Maskable, MetaContainer, FileBased):
    
    def __init__(self, sambam_file1=None, sambam_file2=None, file_name=None,
                       group_name = 'reads', auto_determine_field_sizes=True,
                       field_sizes={'qname': 60, 'sequence': 200}):
        
        
        # Only one file argument and that is an h5dict file name
        if sambam_file1 is not None and sambam_file2 is None and file_name is None:
            file_name = sambam_file1
            sambam_file1 = None
        
        # create or retrieve h5dict file
        if not hasattr(self, 'file') or self.file is None:
            if file_name is None:
                file_name = random_name()
                self.file = create_or_open_pytables_file(file_name, inMemory=True)
            else:
                self.file = create_or_open_pytables_file(file_name, inMemory=False)
        else:
            if not type(self.file) == t.file.File:
                raise ValueError("Object has file attribute, but it is not a pytables File object")

        # try to retrieve existing table
        try:
            main_table = self.file.get_node('/' + group_name + '/main')
            try:
                self._header1 = main_table._v_attrs.header1
                self._header2 = main_table._v_attrs.header2
            except AttributeError:
                logging.warn("No header attributes found in existing table")
                self._header1 = None
                self._header2 = None
            
            try:
                self._ref1 = main_table._v_attrs.ref1
                self._ref2 = main_table._v_attrs.ref2
            except AttributeError:
                logging.warn("No ref attributes found in existing table")
                self._ref1 = None
                self._ref2 = None
            
            # tags
            tags_left = self.file.get_node('/' + group_name + '/tags_left')
            tags_right = self.file.get_node('/' + group_name + '/tags_right')
            
            # cigar
            cigar_left = self.file.get_node('/' + group_name + '/cigar_left')
            cigar_right = self.file.get_node('/' + group_name + '/cigar_right')
        
        # or build table from scratch
        except NoSuchNodeError:
            
            expected_length = 50000000
            if sambam_file1 is not None and sambam_file2 is not None:
                if auto_determine_field_sizes:
                    logging.info("Determining field sizes")
                    lengths1 = ReadPairs.determine_field_sizes(sambam_file1)
                    lengths2 = ReadPairs.determine_field_sizes(sambam_file2)
                    
                    field_sizes['qname'] = max(lengths1["qname"],lengths2["qname"])
                    field_sizes['sequence'] = max(lengths1["sequence"],lengths2["sequence"])
                
                expected_length = int(self._sambam_size(sambam_file1)*1.5)
                
            reads_defininition = {
                'ix': t.Int32Col(pos=0),
                'qname': t.StringCol(field_sizes['qname'],pos=1),
                'flag1': t.Int32Col(pos=2),
                'ref1': t.Int32Col(pos=3),
                'pos1': t.Int64Col(pos=4),
                'mapq1': t.Int32Col(pos=5),
                'rnext1': t.Int32Col(pos=6),
                'pnext1': t.Int32Col(pos=7),
                'tlen1': t.Int32Col(pos=8),
                'seq1': t.StringCol(field_sizes['sequence'],pos=9),
                'qual1': t.StringCol(field_sizes['sequence'],pos=10),
                'flag2': t.Int32Col(pos=11),
                'ref2': t.Int32Col(pos=12),
                'pos2': t.Int64Col(pos=13),
                'mapq2': t.Int32Col(pos=14),
                'rnext2': t.Int32Col(pos=15),
                'pnext2': t.Int32Col(pos=16),
                'tlen2': t.Int32Col(pos=17),
                'seq2': t.StringCol(field_sizes['sequence'],pos=18),
                'qual2': t.StringCol(field_sizes['sequence'],pos=19)
                
            }
            # create data structures
            logging.info("Creating data structures...")
            
            # create reads group
            group = self.file.create_group("/", group_name, 'Read pairs group',
                                           filters=t.Filters(complib="blosc",
                                                             complevel=2, shuffle=True)
                                           )
            # create main table
            main_table = MaskedTable(group, 'main', reads_defininition,
                                                expectedrows=expected_length)
            # create tags vlarrays
            tags_left = self.file.create_vlarray(group, 'tags_left', t.ObjectAtom())
            tags_right = self.file.create_vlarray(group, 'tags_right', t.ObjectAtom())
            
            # create cigar vlarrays
            cigar_left = self.file.create_vlarray(group, 'cigar_left', t.VLStringAtom())
            cigar_right = self.file.create_vlarray(group, 'cigar_right', t.VLStringAtom())
            
            
        
        # generate tables from inherited classes
        Maskable.__init__(self, self.file)
        MetaContainer.__init__(self, self.file)
        
        # make reads table maskable
        self._reads = main_table
        
        self._tags_left = tags_left
        self._tags_right = tags_right
        
        self._cigar_left = cigar_left
        self._cigar_right = cigar_right
        
        # row counter
        self._row_counter = 0

        # map reads
        if sambam_file1 is not None and sambam_file2 is not None:
            logging.info("Loading reads...")
            self.load(sambam_file1, sambam_file2, ignore_duplicates=True)
        logging.info("Done.")
        
    
    
            
    
    @staticmethod
    def determine_field_sizes(sambam, sample_size=10000):
        if type(sambam) == str:
            sambam = pysam.AlignmentFile(sambam, 'rb')  # @UndefinedVariable
            
        qname_length = 0
        seq_length = 0
        cigar_length = 0
        i = 0
        for r in sambam:
            i += 1
            qname_length = max(qname_length,len(r.qname))
            seq_length = max(seq_length,len(r.seq))
            cigar_length = max(cigar_length,len(r.cigarstring))
            if sample_size is not None and i >= sample_size:
                break
            
            if i % 100000 == 0:
                logging.info(i)
        sambam.close()
        
        return { 'qname': qname_length,
                 'sequence': seq_length,
                 'cigar': cigar_length }
        
    def _sambam_size(self, sambam):
        if type(sambam) == str:
            sambam = pysam.AlignmentFile(sambam, 'rb')  # @UndefinedVariable
        
        count = sum(1 for _ in iter(sambam))
        sambam.close()
        return count
        
    
    def load(self, sambam1, sambam2, ignore_duplicates=True, is_sorted=False):
        # get file names
        try:
            file_name1 = sambam1.filename
        except AttributeError:
            file_name1 = sambam1
            
        try:
            file_name2 = sambam2.filename
        except AttributeError:
            file_name2 = sambam2
        
        # sort files if required
        if not is_sorted:
            tmp1 = tempfile.NamedTemporaryFile(delete=False, suffix=".bam")
            tmp1.close()
            tmp2 = tempfile.NamedTemporaryFile(delete=False, suffix=".bam")
            tmp2.close()
            logging.info(tmp1.name)
            logging.info(tmp2.name)
            pysam.sort('-n', file_name1, os.path.splitext(tmp1.name)[0])  # @UndefinedVariable
            pysam.sort('-n', file_name2, os.path.splitext(tmp2.name)[0])  # @UndefinedVariable
            sambam1 = pysam.AlignmentFile(tmp1.name, 'rb')  # @UndefinedVariable
            sambam2 = pysam.AlignmentFile(tmp2.name, 'rb')  # @UndefinedVariable
        else:
            sambam1 = pysam.AlignmentFile(file_name1, 'rb')  # @UndefinedVariable
            sambam2 = pysam.AlignmentFile(file_name2, 'rb')  # @UndefinedVariable
        
        # header
        self._reads._v_attrs.header1 = sambam1.header
        self._header1 = sambam1.header
        self._reads._v_attrs.header2 = sambam2.header
        self._header2 = sambam2.header
        
        # references
        self._reads._v_attrs.ref1 = sambam1.references
        self._reads._v_attrs.ref2 = sambam2.references
        self._ref1 = sambam1.references
        self._ref2 = sambam2.references
        
        iter1 = iter(sambam1)
        iter2 = iter(sambam2)
        def get_next_read(iterator):
            try:
                r = iterator.next()
                return r
            except StopIteration:
                return None
        
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
                self._add_read_pair(r1, r2, flush=False)
                last_r1_name = r1.qname
                last_r2_name = r2.qname
                r1 = get_next_read(iter1)
                r2 = get_next_read(iter2)
                r1_count += 1
                r2_count += 1
            elif c < 0:
                self._add_left_read(r1, flush=False)
                last_r1_name = r1.qname
                r1 = get_next_read(iter1)
                r1_count += 1
            else:
                self._add_right_read(r2, flush=False)
                last_r2_name = r2.qname
                r2 = get_next_read(iter2)
                r2_count += 1
            
            if i % 100000 == 0:
                logging.info("%d reads processed" % i)
                #self._reads.flush(update_index=False)
        
        
        # add remaining unpaired reads
        while r1 is not None:
            if r1.qname == last_r1_name:
                if not ignore_duplicates:
                    raise ValueError("Duplicate left read QNAME %s" % r1.qname)
            else:
                self._add_left_read(r1, flush=False)
            last_r1_name = r1.qname
            r1 = get_next_read(iter1)
            r1_count += 1
        
        while r2 is not None:
            if r2.qname == last_r2_name:
                if not ignore_duplicates:
                    raise ValueError("Duplicate right read QNAME %s" % r2.qname)
            else:
                self._add_right_read(r2, flush=False)
            last_r2_name = r2.qname
            r2 = get_next_read(iter2)
            r2_count += 1
            
        logging.info('Counts: R1 %d R2 %d' % (r1_count,r2_count))
        
        self._reads.flush(update_index=True)
        self._tags_left.flush()
        self._tags_right.flush()
        self._cigar_left.flush()
        self._cigar_right.flush()
        
        if not is_sorted:
            os.unlink(tmp1.name)
            os.unlink(tmp2.name)
    
    
    def _add_left_read(self, r, flush=True):
        row = self._reads.row
        self._add_left_read_to_row(r, row)
        row.append()
        
        self._tags_left.append(r.tags)
        self._tags_right.append({})
        
        self._cigar_left.append(r.cigarstring)
        self._cigar_right.append('')
        
        self._row_counter += 1
        
        if flush:
            self._reads.flush(update_index=True)
            self._tags_left.flush()
            self._tags_right.flush()
            self._cigar_left.flush()
            self._cigar_right.flush()
    

    
    def _add_left_read_to_row(self, r, row):
        row['qname'] = r.qname
        row['flag1'] = r.flag
        row['ref1'] = r.reference_id
        if r.pos >= 0:
            row['pos1'] = r.pos+1
        else:
            row['pos1'] = r.pos
        row['mapq1'] = r.mapq
        row['rnext1'] = r.rnext
        row['pnext1'] =  r.pnext
        row['tlen1'] = r.tlen
        row['seq1'] = r.seq
        row['qual1'] = r.qual
        row['ix'] = self._row_counter
    
    def _add_right_read(self, r, flush=True):
        row = self._reads.row
        self._add_right_read_to_row(r, row)
        row.append()
        
        self._tags_left.append({})
        self._tags_right.append(r.tags)
        
        self._cigar_left.append('')
        self._cigar_right.append(r.cigarstring)
        
        self._row_counter += 1
        
        if flush:
            self._reads.flush(update_index=True)
            self._tags_left.flush()
            self._tags_right.flush()
            self._cigar_left.flush()
            self._cigar_right.flush()
            
    def _add_right_read_to_row(self, r, row):
        row['qname'] = r.qname
        row['flag2'] = r.flag
        row['ref2'] = r.reference_id
        if r.pos >= 0:
            row['pos2'] = r.pos+1
        else:
            row['pos2'] = r.pos
        row['mapq2'] = r.mapq
        row['rnext2'] = r.rnext
        row['pnext2'] =  r.pnext
        row['tlen2'] = r.tlen
        row['seq2'] = r.seq
        row['qual2'] = r.qual
        row['ix'] = self._row_counter        
            
    def _add_read_pair(self, r1, r2, flush=True):
        row = self._reads.row
        self._add_left_read_to_row(r1, row)
        self._add_right_read_to_row(r2, row)
        row.append()
        
        self._tags_left.append(r1.tags)
        self._tags_right.append(r2.tags)
        
        self._cigar_left.append(r1.cigarstring)
        self._cigar_right.append(r2.cigarstring)
        
        self._row_counter += 1
        
        if flush:
            self._reads.flush(update_index=True)
            self._tags_left.flush()
            self._tags_right.flush()
            self._cigar_left.flush()
            self._cigar_right.flush()
    
    def __len__(self):
        return len(self._reads)
    
    @property
    def header1(self):
        return self._header1
    
    @property
    def header2(self):
        return self._header2
    
    def ix2ref1(self, ix):
        if self._ref1 is None:
            raise RuntimeError("Chromosome reference for left read not present")
        return self._ref1[ix]
    
    def ix2ref2(self, ix):
        if self._ref2 is None:
            raise RuntimeError("Chromosome reference for right read not present")
        return self._ref2[ix]

    def filter_quality(self, cutoff=30, queue=False):
        mask = self.add_mask_description('mapq', 'Mask read pairs with a mapping quality lower than %d' % cutoff)
        quality_filter = QualityPairFilter(cutoff, mask)
        quality_filter.set_read_pairs_object(self)
        
        if not queue:
            self._reads.filter(quality_filter)
        else:
            self._reads.queue_filter(quality_filter)
            
    def filter_non_unique(self, strict=True, queue=False):
        mask = self.add_mask_description('uniqueness', 'Mask read pairs that do not map uniquely (according to XS tag)')
        uniqueness_filter = UniquenessPairFilter(strict, mask)
        uniqueness_filter.set_read_pairs_object(self)
        
        if not queue:
            self._reads.filter(uniqueness_filter)
        else:
            self._reads.queue_filter(uniqueness_filter)
            
    def filter_single(self, queue=False):
        mask = self.add_mask_description('single', 'Mask read pairs that are unpaired')
        single_filter = SinglePairFilter(mask)
        single_filter.set_read_pairs_object(self)
        
        if not queue:
            self._reads.filter(single_filter)
        else:
            self._reads.queue_filter(single_filter)
    
    def run_queued_filters(self):
        self._reads.run_queued_filters()
    
    def __iter__(self):
        this = self
        class ReadPairsIter:
            def __init__(self):
                self.iter = iter(this._reads)
                  
            def __iter__(self):
                return self
              
            def next(self):
                row = self.iter.next()
                left_read = None
                if row['pos1'] > 0:
                    left_read = ReadFromRow(row, this, is_left=True)
        
                right_read = None
                if row['pos2'] > 0:
                    right_read = ReadFromRow(row, this, is_left=False)
                
                return ReadPair(left_read=left_read, right_read=right_read)
        return ReadPairsIter()
    
    def _tuple_to_read_pair(self, t):
        main_info = t
        tags_left = self._tags_left[main_info[0]]
        tags_right = self._tags_right[main_info[0]]
        cigar_left = self._cigar_left[main_info[0]]
        cigar_right = self._cigar_right[main_info[0]]
        
        left_read = None
        if main_info[4] > 0:
            left_read = Read(qname=main_info[1], flag=main_info[2], ref=self.ix2ref1(main_info[3]),
                            pos=main_info[4], mapq=main_info[5], cigar=cigar_left, rnext=main_info[6],
                            pnext=main_info[7], tlen=main_info[8], seq=main_info[9], qual=main_info[10],
                            tags=tags_left)
        right_read = None
        if main_info[14] > 0:
            right_read = Read(qname=main_info[1], flag=main_info[11], ref=self.ix2ref2(main_info[12]),
                            pos=main_info[13], mapq=main_info[14], cigar=cigar_right, rnext=main_info[15],
                            pnext=main_info[16], tlen=main_info[17], seq=main_info[18], qual=main_info[19],
                            tags=tags_right)
        return ReadPair(left_read, right_read)
    
    def __getitem__(self, key):
        main_info = self._reads[key]
        return self._tuple_to_read_pair(main_info)
    
    def where(self, query):
        results = [x.fetch_all_fields() for x in self._reads.where(query)]
        pairs = []
        for result in results:
            pairs.append(self._tuple_to_read_pair(result))
        return pairs
    
    def filtered_reads(self):
        this = self
        class MaskedReadPairsIter:
            def __init__(self):
                self.iter = this._reads.masked_rows()
                  
            def __iter__(self):
                return self
              
            def next(self):
                row = self.iter.next()
                left_read = None
                if row['pos1'] > 0:
                    left_read = ReadFromRow(row, this, is_left=True)
        
                right_read = None
                if row['pos2'] > 0:
                    right_read = ReadFromRow(row, this, is_left=False)
                
                masks = this.get_masks(row['mask'])
                
                return MaskedReadPair(left_read=left_read, right_read=right_read, masks=masks)
        return MaskedReadPairsIter()



class Read(object):
    def __init__(self, qname="", flag=0, ref="",
                 pos=0, mapq=0, cigar="", rnext=0,
                 pnext=0, tlen=0, seq="", qual="",
                 tags={}):
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

    @property
    def strand(self):
        bit_flags = bit_flags_from_int(self.flag)
        if 4 in bit_flags:
            return -1
        return 1
    
    def __getitem__(self, key):
        try:
            value = self.__getattribute__(key)
            return value
        except:
            raise KeyError("Read does not have %s attribute" % str(key))
        
    def __repr__(self):
        return "%s, ref: %s, pos: %d" % (self.qname, self.ref, self.pos)
    
class MaskedRead(Read):
    def __init__(self, qname="", flag=0, ref="",
                 pos=0, mapq=0, cigar="", rnext=0,
                 pnext=0, tlen=0, seq="", qual="",
                 tags={}, masks=None):
        super(MaskedRead, self).__init__(qname=qname, flag=flag, ref=ref,
                                         pos=pos, mapq=mapq, cigar=cigar, rnext=rnext,
                                         pnext=pnext, tlen=tlen, seq=seq, qual=qual,
                                         tags=tags)
        self.masks = masks
    
    def __repr__(self):
        representation = super(MaskedRead, self).__repr__()
        if self.masks is not None:
            mask_names = []
            for mask in self.masks:
                mask_names.append(mask.name)
            return "%s (%s)" % (representation,", ".join(mask_names))
        return representation
        
class ReadFromRow(Read):
    def __init__(self, row, read_pairs_object, is_left=True):
        self.row = row
        self.read_pairs = read_pairs_object
        self.is_left = is_left
        
        if is_left:
            self.side = 1
        else:
            self.side = 2
    
    @property
    def qname(self): return self.row['qname']
    
    @property
    def flag(self): return self.row['flag' + str(self.side)]
    
    @property
    def ref(self):
        if self.is_left:
            return self.read_pairs.ix2ref1(self.row['ref' + str(self.side)])
        return self.read_pairs.ix2ref2(self.row['ref' + str(self.side)])
    
    @property
    def pos(self): return self.row['pos' + str(self.side)]
    
    @property
    def mapq(self): return self.row['mapq' + str(self.side)]
    
    @property
    def rnext(self): return self.row['rnext' + str(self.side)]
    
    @property
    def pnext(self): return self.row['pnext' + str(self.side)]
    
    @property
    def tlen(self): return self.row['tlen' + str(self.side)]
    
    @property
    def seq(self): return self.row['seq' + str(self.side)]
    
    @property
    def qual(self): return self.row['qual' + str(self.side)]
    
    @property
    def cigar(self):
        ix = self.row['ix']
        if self.is_left:
            return self.read_pairs._cigar_left[ix]
        return self.read_pairs._cigar_right[ix]
    
    @property
    def tags(self):
        ix = self.row['ix']
        if self.is_left:
            return self.read_pairs._tags_left[ix]
        return self.read_pairs._tags_right[ix]
    
        

class ReadPair(object):
    def __init__(self, left_read=None, right_read=None):
        self.left_read = left_read
        self.right_read = right_read
    
    def has_left_read(self):
        return self.left_read is not None
    
    def has_right_read(self):
        return self.right_read is not None
    
    @property
    def qname(self):
        if self.has_left_read():
            return self.left_read.qname
        
        if self.has_right_read():
            return self.right_read.qname
        
        return None
    
    def __repr__(self):
        left_repr = "None"
        if self.has_left_read():
            left_repr = "%s-%d" % (self.left_read.ref, self.left_read.pos)
        right_repr = "None"
        if self.has_right_read():
            right_repr = "%s-%d" % (self.right_read.ref, self.right_read.pos)
        return "%s: (%s)-(%s)" % (self.qname, left_repr, right_repr)
    
class MaskedReadPair(ReadPair):
    def __init__(self, left_read, right_read, masks=None):
        super(MaskedReadPair, self).__init__(left_read,right_read)
        self.masks = masks
        
    def __repr__(self):
        representation = super(MaskedReadPair, self).__repr__()
        if self.masks is not None:
            mask_names = []
            for mask in self.masks:
                mask_names.append(mask.name)
            return "%s (%s)" % (representation,", ".join(mask_names))
        return representation


#
# Filters Reads
#
class ReadFilter(MaskFilter):
    __metaclass__ = ABCMeta
    
    def __init__(self, mask=None):
        super(ReadFilter, self).__init__(mask)
    
    @abstractmethod
    def valid_read(self, read):
        pass
    
    def set_reads_object(self, reads_object):
        self._reads = reads_object
    
    def valid(self, row):
        read = self._reads._row2read(row)
        return self.valid_read(read)
        

class QualityFilter(ReadFilter):
    def __init__(self, cutoff=30, mask=None):
        super(QualityFilter, self).__init__(mask)
        self.cutoff = cutoff

    def valid_read(self, read):
        return read.mapq > self.cutoff

class UniquenessFilter(ReadFilter):
    def __init__(self, strict=True, mask=None):
        self.strict = strict
        super(UniquenessFilter, self).__init__(mask)
    
    def valid_read(self, read):

        for tag in read.tags:
            if tag[0] == 'XS':
                if self.strict or tag[1] == 0:
                    return False
        return True

class UnmappedFilter(ReadFilter):
    def __init__(self, mask=None):
        super(UnmappedFilter, self).__init__(mask)

    def valid_read(self, read):
        if 2 in bit_flags_from_int(read.flag):
            return False
        return True



class ReadPairFilter(MaskFilter):
    __metaclass__ = ABCMeta
    
    def __init__(self, mask=None):
        super(ReadPairFilter, self).__init__(mask)
    
    @abstractmethod
    def valid_pair(self, read_pair):
        pass
    
    def set_read_pairs_object(self, read_pairs_object):
        self._read_pairs = read_pairs_object
    
    def valid(self, row):
        
        left_read = None
        if row['pos1'] > 0:
            left_read = ReadFromRow(row, self._read_pairs, is_left=True)

        
        right_read = None
        if row['pos2'] > 0:
            right_read = ReadFromRow(row, self._read_pairs, is_left=False)
        
        pair = ReadPair(left_read=left_read, right_read=right_read)
        
        return self.valid_pair(pair)
        

class QualityPairFilter(ReadPairFilter):
    def __init__(self, cutoff=30, mask=None):
        super(QualityPairFilter, self).__init__(mask)
        self.cutoff = cutoff

    def valid_pair(self, pair):
        left_valid = True
        if pair.has_left_read() and pair.left_read.mapq < self.cutoff:
            left_valid = False
        
        right_valid = True
        if pair.has_right_read() and pair.right_read.mapq < self.cutoff:
            right_valid=False
             
        return left_valid and right_valid

class UniquenessPairFilter(ReadPairFilter):
    def __init__(self, strict=True, mask=None):
        self.strict = strict
        super(UniquenessPairFilter, self).__init__(mask)
    
    def valid_pair(self, pair):
        left_valid = True
        if pair.has_left_read():
            for tag in pair.left_read.tags:
                if tag[0] == 'XS':
                    if self.strict or tag[1] == 0:
                        left_valid = False
                    
        right_valid = True
        if pair.has_right_read() > 0:
            right_tags = pair.right_read.tags
            for tag in right_tags:
                if tag[0] == 'XS':
                    if self.strict or tag[1] == 0:
                        right_valid = False
                    
        return left_valid and right_valid


class SinglePairFilter(ReadPairFilter):
    def __init__(self, mask=None):
        super(SinglePairFilter, self).__init__(mask)
    
    def valid_pair(self, pair):
        if not pair.has_left_read() or not pair.has_right_read():
            return False
        return True





class FragmentMappedReadPairs(Maskable, MetaContainer, RegionsTable, FileBased):
    class FragmentMappedReadDescription(t.IsDescription):
        ix = t.Int32Col(pos=0)
        position = t.Int64Col(pos=2)
        strand = t.Int8Col(pos=3)
    
    class FragmentsMappedReadPairDescription(t.IsDescription):
        ix = t.Int32Col(pos=0)
        left_read = t.Int32Col(pos=1)
        left_fragment = t.Int32Col(pos=2, dflt=-1)
        right_read = t.Int32Col(pos=3)
        right_fragment = t.Int32Col(pos=4, dflt=-1)
    
    class FragmentsMappedReadSingleDescription(t.IsDescription):
        ix = t.Int32Col(pos=0)
        read = t.Int32Col(pos=1)
        fragment = t.Int32Col(pos=2, dflt=-1)
        
    def __init__(self, data=None, file_name=None,
                 group_name = 'fragment_map',
                 table_name_fragments='fragments'):
        
        # parse potential unnamed argument
        if data is not None:
            # data is file name
            if type(data) is str:
                data = os.path.expanduser(data)
                
                if file_name is None:
                    file_name = data
                    data = None
        
        if file_name is not None and isinstance(file_name, str):
            file_name = os.path.expanduser(file_name)
                
        #FileBased.__init__(self, file_name)
        RegionsTable.__init__(self, file_name=file_name, table_name_regions=table_name_fragments)
        
        # generate tables from inherited classes
        Maskable.__init__(self, self.file)
        MetaContainer.__init__(self, self.file)
        
        
        # try to retrieve existing table
        try:
            self._reads = self.file.get_node('/' + group_name + '/mapped_reads')
            self._pairs = self.file.get_node('/' + group_name + '/mapped_read_pairs')
            self._single = self.file.get_node('/' + group_name + '/mapped_read_single')
            self._read_count = len(self._reads)
            self._pair_count = len(self._pairs)
            self._single_count = len(self._single)
        # or build table from scratch
        except NoSuchNodeError:
            # create group
            group = self.file.create_group("/", group_name, 'Mapped read pairs group',
                                           filters=t.Filters(complib="blosc",
                                                             complevel=2, shuffle=True))
            # create main tables
            self._reads = t.Table(group, 'mapped_reads',
                                    FragmentMappedReadPairs.FragmentMappedReadDescription,
                                    expectedrows=10000000)
            
            self._pairs = MaskedTable(group, 'mapped_read_pairs',
                                    FragmentMappedReadPairs.FragmentsMappedReadPairDescription,
                                    expectedrows=5000000)
            
            self._single = MaskedTable(group, 'mapped_read_single',
                                    FragmentMappedReadPairs.FragmentsMappedReadSingleDescription,
                                    expectedrows=1000000)
            self._read_count = 0
            self._pair_count = 0
            self._single_count = 0
        
        try:
            self._pairs.cols.left_fragment.create_csindex()
        except ValueError:
            # Index exists, no problem!
            pass
    
    def close(self):
        self.file.close()
    
    def load(self, reads1, reads2, regions=None, ignore_duplicates=True, _in_memory_index=True):
        if regions is not None:
            logging.info("Adding regions...")
            self.add_regions(regions)
            logging.info("Done.")
        
        # generate index for fragments
        fragment_ixs = None
        fragment_ends = None
        if _in_memory_index:
            fragment_ixs = {}
            fragment_ends = {}
            for region in self.regions():
                if not fragment_ends.has_key(region.chromosome):
                    fragment_ends[region.chromosome] = []
                    fragment_ixs[region.chromosome] = []
                fragment_ixs[region.chromosome].append(region.ix)
                fragment_ends[region.chromosome].append(region.end)
        
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
                self._add_read_pair(r1, r2, flush=False, _fragment_ends=fragment_ends, _fragment_ixs=fragment_ixs)
                last_r1_name = r1.qname
                last_r2_name = r2.qname
                r1 = get_next_read(iter1)
                r2 = get_next_read(iter2)
                r1_count += 1
                r2_count += 1
            elif c < 0:
                self._add_read_single(r1, flush=False, _fragment_ends=fragment_ends, _fragment_ixs=fragment_ixs)
                last_r1_name = r1.qname
                r1 = get_next_read(iter1)
                r1_count += 1
            else:
                self._add_read_single(r2, flush=False, _fragment_ends=fragment_ends, _fragment_ixs=fragment_ixs)
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
                self._add_read_single(r1, flush=False, _fragment_ends=fragment_ends, _fragment_ixs=fragment_ixs)
            last_r1_name = r1.qname
            r1 = get_next_read(iter1)
            r1_count += 1
        
        while r2 is not None:
            if r2.qname == last_r2_name:
                if not ignore_duplicates:
                    raise ValueError("Duplicate right read QNAME %s" % r2.qname)
            else:
                self._add_read_single(r2, flush=False, _fragment_ends=fragment_ends, _fragment_ixs=fragment_ixs)
            last_r2_name = r2.qname
            r2 = get_next_read(iter2)
            r2_count += 1
            
        logging.info('Counts: R1 %d R2 %d' % (r1_count,r2_count))
        
        self._reads.flush()
        self._pairs.flush(update_index=True)
        self._single.flush(update_index=True)
        
        
    def _add_read_pair(self, read1, read2, flush=True, _fragment_ends=None, _fragment_ixs=None):
        ix1 = self._add_read(read1, flush=flush)
        ix2 = self._add_read(read2, flush=flush)
        fragment_ix1 = self._find_fragment_ix(read1.ref, read1.pos, _fragment_ends=_fragment_ends, _fragment_ixs=_fragment_ixs)
        fragment_ix2 = self._find_fragment_ix(read2.ref, read2.pos, _fragment_ends=_fragment_ends, _fragment_ixs=_fragment_ixs)
        
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
            
    
    def _add_read_single(self, read, flush=True, _fragment_ends=None, _fragment_ixs=None):
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
            for row in self._regions.where("(start <= %d) & (end >= %d) & (chromosome == '%s')" % (position, position, chromosome)):
                fragment_ix = row['ix']
        
        return fragment_ix
        
    def _add_read(self, read, flush=True):
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
    
    def _pair_from_row(self, row):
        ix1 = row['left_read']
        ix2 = row['right_read']
        fragment_ix1 = row['left_fragment']
        fragment_ix2 = row['right_fragment']
        
        read1_row = self._reads[ix1]
        fragment1_row = self._regions[fragment_ix1]
        fragment1 = GenomicRegion(fragment1_row['start'],fragment1_row['end'], fragment1_row['chromosome'])
        read1 = FragmentRead(fragment1, position = read1_row['position'], strand = read1_row['strand'])
        
        read2_row = self._reads[ix2]
        fragment2_row = self._regions[fragment_ix2]
        fragment2 = GenomicRegion(fragment2_row['start'],fragment2_row['end'], fragment2_row['chromosome'])
        read2 = FragmentRead(fragment2, position = read2_row['position'], strand = read2_row['strand'])
        
        return FragmentReadPair(left_read=read1, right_read=read2)
    
    
    def plot_error_structure(self, output=None,data_points=None,
                             skip_self_ligations=True):
        
        same_count = 0
        inward_count = 0
        outward_count = 0
        gaps = []
        types = []
        type_same = 0
        type_inward = 1
        type_outward = 2
        l = len(self)
        last_percent = -1
        for i, pair in enumerate(self):
            if i % int(l/20) == 0:
                percent = int(i/int(l/20))
                if percent != last_percent:
                    print "%d%% done" % (percent*5)
                    last_percent = percent
                    
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
        
        # sort data
        points = zip(gaps,types)
        sortedPoints = sorted(points)
        gaps = [point[0] for point in sortedPoints]
        types = [point[1] for point in sortedPoints]
        
        # best guess for number of data points
        if data_points is None:
            data_points = max(100, int(l * 0.0025))
        logging.info("Number of data points averaged per point in plot: %d" % data_points)
        
        # calculate ratios
        x = []
        inwardRatios = []
        outwardRatios = []
        counter = 0
        sameCounter = 0
        mids = 0
        outwards = 0
        inwards = 0
        same = 0
        for i in xrange(0,len(gaps)):
            mids += gaps[i]
            if types[i] == type_same:
                same += 1
                sameCounter += 1
            elif types[i] == type_inward:
                inwards += 1
            else:
                outwards += 1
            counter += 1
            
            if sameCounter > data_points:
                x.append(mids/counter)
                inwardRatios.append(inwards/same)
                outwardRatios.append(outwards/same)
                
                sameCounter = 0
                counter = 0
                mids = 0
                outwards = 0
                inwards = 0
                same = 0
                
        # plot
        if output != None:
            old_backend = plt.get_backend()
            plt.switch_backend('pdf')
            plt.ioff()
        
        fig = plt.figure()
        fig.suptitle("Error structure by distance")
        plt.plot(x,inwardRatios, 'b', label="inward/same strand")
        plt.plot(x,outwardRatios, 'r', label="outward/same strand")
        plt.xscale('log')
        plt.axhline(y=0.5,color='black',ls='dashed')
        plt.ylim(0,3)
        plt.xlabel('gap size between fragments')
        plt.ylabel('ratio of number of reads')
        plt.legend(loc='upper right')

        if output == None:
            plt.show()
        else:
            fig.savefig(output)
            plt.close(fig)
            plt.ion()
            plt.switch_backend(old_backend)
    
    def filter(self, pair_filter, queue=False):
        pair_filter.set_pairs_object(self)
        if not queue:
            self._pairs.filter(pair_filter)
        else:
            self._pairs.queue_filter(pair_filter)
    
    def run_queued_filters(self):
        self._pairs.run_queued_filters()
        
    def filter_inward(self, minimum_distance, queue=False):
        mask = self.add_mask_description('inward', 'Mask read pairs that are inward facing and <%dbp apart' % minimum_distance)
        inward_filter = InwardPairsFilter(minimum_distance=minimum_distance, mask=mask)
        self.filter(inward_filter, queue)
    
    def filter_outward(self, minimum_distance, queue=False):
        mask = self.add_mask_description('outward', 'Mask read pairs that are outward facing and <%dbp apart' % minimum_distance)
        outward_filter = OutwardPairsFilter(minimum_distance=minimum_distance, mask=mask)
        self.filter(outward_filter, queue)
    
    def filter_re_dist(self, maximum_distance, queue=False):
        mask = self.add_mask_description('re-dist', 'Mask read pairs where a read is >%dbp away from nearest RE site' % maximum_distance)
        re_filter = ReDistanceFilter(maximum_distance=maximum_distance,mask=mask)
        self.filter(re_filter, queue)
    
    def __iter__(self):
        this = self
        class FragmentMappedReadPairIter:
            def __init__(self):
                self.iter = iter(this._pairs)
                 
            def __iter__(self):
                return self
             
            def next(self):
                return this._pair_from_row(self.iter.next())
             
            def __len__(self):
                return len(this._pairs)
             
        return FragmentMappedReadPairIter()
    
    def __getitem__(self, key):
        row = self._pairs[key]
        return self._pair_from_row(row)
    
    def __len__(self):
        return len(self._pairs)
    


class FragmentRead(object):
    def __init__(self, fragment=None, position=None, strand=0):
        self.fragment = fragment
        self.position = position
        self.strand = strand
    
    def __repr__(self):
        return "%s: %d-(%d[%d])-%d" % (self.fragment.chromosome,
                                   self.fragment.start,
                                   self.position,
                                   self.strand,
                                   self.fragment.end)
        
class FragmentReadPair(object):
    def __init__(self, left_read, right_read):
        self.left = left_read
        self.right = right_read
    
    def is_same_chromosome(self):
        return self.left.fragment.chromosome == self.right.fragment.chromosome
    
    def is_inward_pair(self):
        if not self.is_same_chromosome():
            return False
        
        if self.left.strand == 1 and self.right.strand == -1:
            return True
        return False
    
    def is_outward_pair(self):
        if not self.is_same_chromosome():
            return False
        
        if self.left.strand == -1 and self.right.strand == 1:
            return True
        return False
    
    def is_same_pair(self):
        if not self.is_same_chromosome():
            return False
        
        if self.left.strand == self.right.strand:
            return True
        return False
    
    def is_same_fragment(self):
        return self.left.fragment.start == self.right.fragment.start
    
    def get_gap_size(self):
        if not self.is_same_chromosome():
            return None
        
        if self.is_same_fragment():
            return 0
        
        gap = self.right.fragment.start - self.left.fragment.end
        
        if gap == 1: # neighboring fragments
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
    __metaclass__ = ABCMeta
    
    def __init__(self, mask=None):
        super(FragmentMappedReadPairFilter, self).__init__(mask)
    
    def set_pairs_object(self, pairs):
        self.pairs = pairs
    
    @abstractmethod
    def valid_pair(self, fr_pair):
        pass
    
    def valid(self, row):
        pair = self.pairs._pair_from_row(row)
        return self.valid_pair(pair)


class InwardPairsFilter(FragmentMappedReadPairFilter):
    def __init__(self, minimum_distance=10000, mask=None):
        super(InwardPairsFilter, self).__init__(mask=mask)
        self.minimum_distance = minimum_distance
    
    def valid_pair(self, pair):        
        if not pair.is_inward_pair():
            return True
        
        if pair.get_gap_size() > self.minimum_distance:
            return True
        return False
    
class OutwardPairsFilter(FragmentMappedReadPairFilter):
    def __init__(self, minimum_distance=10000, mask=None):
        super(OutwardPairsFilter, self).__init__(mask=mask)
        self.minimum_distance = minimum_distance
    
    def valid_pair(self, pair):        
        if not pair.is_outward_pair():
            return True
        
        if pair.get_gap_size() > self.minimum_distance:
            return True
        return False
    
class ReDistanceFilter(FragmentMappedReadPairFilter):
    def __init__(self, maximum_distance=500, mask=None):
        super(ReDistanceFilter, self).__init__(mask=mask)
        self.maximum_distance = maximum_distance
    
    def valid_pair(self, pair):
        for read in [pair.left, pair.right]:
            if (read.position - read.fragment.start > self.maximum_distance
                and read.fragment.end - read.position > self.maximum_distance):
                return False
        return True


def cmp_natural(string1, string2):
    
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

        
        
        
        
