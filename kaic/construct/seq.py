'''
Created on Jul 13, 2015

@author: kkruse1
'''

import tables as t
import pysam
from kaic.tools.files import create_or_open_pytables_file, random_name
from kaic.data.general import Maskable, MetaContainer, MaskFilter, MaskedTable
import pickle
import tempfile
import os
import logging
from tables.exceptions import NoSuchNodeError



class ReadPairs(Maskable, MetaContainer):
    
    def __init__(self, sambam_file1=None, sambam_file2=None, file_name=None,
                       table_name = 'reads', auto_determine_field_sizes=True,
                       field_sizes={'qname': 60,
                                    'sequence': 200,
                                    'tags': 250,
                                    'cigar': 50}):
        
        
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
            self._reads = self.file.get_node('/' + table_name)
            try:
                self._header1 = self._reads._v_attrs.header1
                self._header2 = self._reads._v_attrs.header2
            except AttributeError:
                self._header1 = None
                self._header2 = None
            
            try:
                self._ref1 = self._reads._v_attrs.ref1
                self._ref2 = self._reads._v_attrs.ref2
            except AttributeError:
                self._ref1 = None
                self._ref2 = None
        
        # or build table from scratch
        except NoSuchNodeError:
            if sambam_file1 is not None and sambam_file2 is not None:
                if auto_determine_field_sizes:
                    logging.info("Determining field sizes")
                    lengths1 = ReadPairs.determine_field_sizes(sambam_file1)
                    lengths2 = ReadPairs.determine_field_sizes(sambam_file2)
                    
                    field_sizes['qname'] = max(lengths1["qname"],lengths2["qname"])
                    field_sizes['sequence'] = max(lengths1["sequence"],lengths2["sequence"])
                    field_sizes['cigar'] = max(lengths1["cigar"],lengths2["cigar"])
                    field_sizes['tags'] = max(lengths1["tags"],lengths2["tags"])
                
            reads_defininition = {
                'qname': t.StringCol(field_sizes['qname'],pos=0),
                'flag1': t.Int16Col(pos=1),
                'ref1': t.Int32Col(pos=2),
                'pos1': t.Int64Col(pos=3),
                'mapq1': t.Int32Col(pos=4),
                'cigar1': t.StringCol(field_sizes['cigar'],pos=5),
                'rnext1': t.Int32Col(pos=6),
                'pnext1': t.Int32Col(pos=7),
                'tlen1': t.Int32Col(pos=8),
                'seq1': t.StringCol(field_sizes['sequence'],pos=9),
                'qual1': t.StringCol(field_sizes['sequence'],pos=10),
                'tags1': t.StringCol(field_sizes['tags'],pos=11),
                'flag2': t.Int16Col(pos=12),
                'ref2': t.Int32Col(pos=13),
                'pos2': t.Int64Col(pos=14),
                'mapq2': t.Int32Col(pos=15),
                'cigar2': t.StringCol(field_sizes['cigar'],pos=16),
                'rnext2': t.Int32Col(pos=17),
                'pnext2': t.Int32Col(pos=18),
                'tlen2': t.Int32Col(pos=19),
                'seq2': t.StringCol(field_sizes['sequence'],pos=20),
                'qual2': t.StringCol(field_sizes['sequence'],pos=21),
                'tags2': t.StringCol(field_sizes['tags'],pos=22),
                'mask': t.Int16Col(pos=23),
                'mask_ix': t.Int32Col(pos=24)
            }
            # create table
            logging.info("Creating tables...")
            self._reads = self.file.create_table("/", table_name, reads_defininition)
        
        # generate tables from inherited classes
        Maskable.__init__(self, self.file)
        MetaContainer.__init__(self, self.file)
        
        # make reads table maskable
        self._reads = MaskedTable(self._reads)

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
        tags_length = 0
        i = 0
        for r in sambam:
            i += 1
            qname_length = max(qname_length,len(r.qname))
            seq_length = max(seq_length,len(r.seq))
            cigar_length = max(cigar_length,len(r.cigarstring))
            tags_length = max(tags_length,len(pickle.dumps(r.tags)))
            if sample_size is not None and i >= sample_size:
                break
        sambam.close()
        
        return { 'qname': qname_length,
                 'sequence': seq_length,
                 'tags': tags_length,
                 'cigar': cigar_length }
        
        
    
    
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
        
        # turn off auto-indexing (minimal speed gain)
        self._reads.autoindex = False
        
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
                self._reads.flush(update_index=False)
        
        
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
        
        self._reads.flush()
        self._reads.autoindex = True
        
        if not is_sorted:
            os.unlink(tmp1.name)
            os.unlink(tmp2.name)
    
    
    def _add_left_read(self, r, flush=True):
        row = self._reads.row
        self._add_left_read_to_row(r, row)
        row.append()
        
        if flush:
            self._reads.flush()
            
    def _add_left_read_to_row(self, r, row):
        row['qname'] = r.qname
        row['flag1'] = r.flag
        row['ref1'] = r.reference_id
        if r.pos >= 0:
            row['pos1'] = r.pos+1
        else:
            row['pos1'] = r.pos
        row['mapq1'] = r.mapq
        row['cigar1'] = r.cigarstring
        row['rnext1'] = r.rnext
        row['pnext1'] =  r.pnext
        row['tlen1'] = r.tlen
        row['seq1'] = r.seq
        row['qual1'] = r.qual
        row['tags1'] = pickle.dumps(r.tags)
    
    def _add_right_read(self, r, flush=True):
        row = self._reads.row
        self._add_right_read_to_row(r, row)
        row.append()
        
        if flush:
            self._reads.flush()
            
    def _add_right_read_to_row(self, r, row):
        row['qname'] = r.qname
        row['flag2'] = r.flag
        row['ref2'] = r.reference_id
        if r.pos >= 0:
            row['pos2'] = r.pos+1
        else:
            row['pos2'] = r.pos
        row['mapq2'] = r.mapq
        row['cigar2'] = r.cigarstring
        row['rnext2'] = r.rnext
        row['pnext2'] =  r.pnext
        row['tlen2'] = r.tlen
        row['seq2'] = r.seq
        row['qual2'] = r.qual
        row['tags2'] = pickle.dumps(r.tags)
        
            
    def _add_read_pair(self, r1, r2, flush=True):
        row = self._reads.row
        self._add_left_read_to_row(r1, row)
        self._add_right_read_to_row(r2, row)
        row.append()
        
        if flush:
            self._reads.flush()
    
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
        quality_filter = QualityFilter(cutoff, mask)
        
        if not queue:
            self._reads.filter(quality_filter)
        else:
            self._reads.queue_filter(quality_filter)
            
    def filter_non_unique(self, strict=True, queue=False):
        mask = self.add_mask_description('uniqueness', 'Mask read pairs that do not map uniquely (according to XS tag)')
        uniqueness_filter = UniquenessFilter(strict, mask)
        
        if not queue:
            self._reads.filter(uniqueness_filter)
        else:
            self._reads.queue_filter(uniqueness_filter)
            
    def filter_single(self, queue=False):
        mask = self.add_mask_description('single', 'Mask read pairs that are unpaired')
        print "MASK: %d" % mask.ix
        single_filter = SingleFilter(mask)
        
        if not queue:
            self._reads.filter(single_filter)
        else:
            self._reads.queue_filter(single_filter)
    
    def run_queued_filters(self):
        self._reads.run_queued_filters()
    
    def __iter__(self):
        return iter(self._reads)
                
    
    def __getitem__(self, key):
        return self._reads.__getitem__(key)
    
    
    def filtered_reads(self):
        return self._reads.masked_rows()
    
    
    

#
# Filters
#


class QualityFilter(MaskFilter):
    def __init__(self, cutoff=30, mask=None):
        super(QualityFilter, self).__init__(mask)
        self.cutoff = cutoff

    def valid(self, row):
        left_valid = True
        if row['pos1'] > 0 and row['mapq1'] < self.cutoff:
            left_valid = False
        
        right_valid = True
        if row['pos2'] > 0 and row['mapq2'] < self.cutoff:
            right_valid=False
             
        return left_valid and right_valid

class UniquenessFilter(MaskFilter):
    def __init__(self, strict=True, mask=None):
        self.strict = strict
        super(UniquenessFilter, self).__init__(mask)
    
    def valid(self, row):
        left_valid = True
        if row['pos1'] > 0:
            left_tags = pickle.loads(row['tags1'])
            for tag in left_tags:
                if tag[0] == 'XS':
                    if self.strict or tag[1] == 0:
                        left_valid = False
                    
        right_valid = True
        if row['pos2'] > 0:
            right_tags = pickle.loads(row['tags2'])
            for tag in right_tags:
                if tag[0] == 'XS':
                    if self.strict or tag[1] == 0:
                        right_valid = False
                    
        return left_valid and right_valid


class SingleFilter(MaskFilter):
    def __init__(self, mask=None):
        super(SingleFilter, self).__init__(mask)
    
    def valid(self, row):
        if row['pos1'] == 0 or row['pos2'] == 0:
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

        
        
        
        
