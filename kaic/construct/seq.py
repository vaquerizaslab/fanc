'''
Created on Jul 13, 2015

@author: kkruse1
'''

import tables as t
import pysam
from kaic.tools.files import create_or_open_pytables_file, random_name
from kaic.data.general import Maskable, MetaContainer, MaskFilter, MaskedTable
import tempfile
import os
import logging
from tables.exceptions import NoSuchNodeError
from abc import abstractmethod, ABCMeta



class ReadPairs(Maskable, MetaContainer):
    
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
        
        self._reads.flush()
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
            self._reads.flush()
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
            self._reads.flush()
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
            self._reads.flush()
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
        quality_filter = QualityFilter(cutoff, mask)
        quality_filter.set_read_pairs_object(self)
        
        if not queue:
            self._reads.filter(quality_filter)
        else:
            self._reads.queue_filter(quality_filter)
            
    def filter_non_unique(self, strict=True, queue=False):
        mask = self.add_mask_description('uniqueness', 'Mask read pairs that do not map uniquely (according to XS tag)')
        uniqueness_filter = UniquenessFilter(strict, mask)
        uniqueness_filter.set_read_pairs_object(self)
        
        if not queue:
            self._reads.filter(uniqueness_filter)
        else:
            self._reads.queue_filter(uniqueness_filter)
            
    def filter_single(self, queue=False):
        mask = self.add_mask_description('single', 'Mask read pairs that are unpaired')
        print "MASK: %d" % mask.ix
        single_filter = SingleFilter(mask)
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

    def __getitem__(self, key):
        try:
            value = self.__getattribute__(key)
            return value
        except:
            raise KeyError("Read does not have %s attribute" % str(key))
        
    def __repr__(self):
        return "%s, ref: %s, pos: %d" % (self.qname, self.ref, self.pos)
        
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
# Filters
#
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
        

class QualityFilter(ReadPairFilter):
    def __init__(self, cutoff=30, mask=None):
        super(QualityFilter, self).__init__(mask)
        self.cutoff = cutoff

    def valid_pair(self, pair):
        left_valid = True
        if pair.has_left_read() and pair.left_read.mapq < self.cutoff:
            left_valid = False
        
        right_valid = True
        if pair.has_right_read() and pair.right_read.mapq < self.cutoff:
            right_valid=False
             
        return left_valid and right_valid

class UniquenessFilter(ReadPairFilter):
    def __init__(self, strict=True, mask=None):
        self.strict = strict
        super(UniquenessFilter, self).__init__(mask)
    
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


class SingleFilter(ReadPairFilter):
    def __init__(self, mask=None):
        super(SingleFilter, self).__init__(mask)
    
    def valid_pair(self, pair):
        if not pair.has_left_read() or not pair.has_right_read():
            return False
        return True





class FragmentPairs(object):
    def __init__(self):
        pass




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

        
        
        
        
