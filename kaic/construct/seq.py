'''
Created on Jul 13, 2015

@author: kkruse1
'''

import tables as t
import pysam
from kaic.tools.files import create_or_open_pytables_file, random_name
import pickle
import logging

class ReadPairs(object):
    
    def __init__(self, sambam_file1, sambam_file2=None, file_name=None,
                       table_name_reads = 'reads', max_qname_len=60,
                       max_seq_len=100, max_tags_len=100, max_cigar_len=10):
        if sambam_file2 is None:
            file_name = sambam_file1
            self.file = create_or_open_pytables_file(file_name, inMemory=False)
            self._reads = self.file.get_node('/' + table_name_reads)
            self._header1 = pickle.loads(self._reads._v_attrs.header1)
            self._header2 = pickle.loads(self._reads._v_attrs.header2)
        else:
            
            if file_name is None:
                file_name = random_name()
                self.file = create_or_open_pytables_file(file_name, inMemory=True)
            else:
                self.file = create_or_open_pytables_file(file_name, inMemory=False)            
            
#             logging.info("Determining field sizes...")
#             lengths1 = self._determine_field_sizes(sambam_file1)
#             lengths2 = self._determine_field_sizes(sambam_file2)
#             
#             max_qname = max(lengths1["qname"],lengths2["qname"])
#             max_seq = max(lengths1["sequence"],lengths2["sequence"])
#             max_cigar = max(lengths1["cigar"],lengths2["cigar"])
#             max_tags = max(lengths1["tags"],lengths2["tags"])
            
            reads_defininition = {
                'qname': t.StringCol(max_qname_len,pos=0),
                'flag1': t.Int16Col(pos=1),
                'ref1': t.Int32Col(pos=2),
                'pos1': t.Int64Col(pos=3),
                'mapq1': t.Int32Col(pos=4),
                'cigar1': t.StringCol(max_cigar_len,pos=5),
                'rnext1': t.Int32Col(pos=6),
                'pnext1': t.Int32Col(pos=7),
                'tlen1': t.Int32Col(pos=8),
                'seq1': t.StringCol(max_seq_len,pos=9),
                'qual1': t.StringCol(max_seq_len,pos=10),
                'tags1': t.StringCol(max_tags_len,pos=11),
                'flag2': t.Int16Col(pos=12),
                'ref2': t.Int32Col(pos=13),
                'pos2': t.Int64Col(pos=14),
                'mapq2': t.Int32Col(pos=15),
                'cigar2': t.StringCol(max_cigar_len,pos=16),
                'rnext2': t.Int32Col(pos=17),
                'pnext2': t.Int32Col(pos=18),
                'tlen2': t.Int32Col(pos=19),
                'seq2': t.StringCol(max_seq_len,pos=20),
                'qual2': t.StringCol(max_seq_len,pos=21),
                'tags2': t.StringCol(max_tags_len,pos=22)
            }
            
            # create table
            logging.info("Creating reads table...")
            self._reads = self.file.create_table("/", table_name_reads, reads_defininition)
            sambam1 = pysam.AlignmentFile(sambam_file1, 'rb')  # @UndefinedVariable
            self._reads._v_attrs.header2 = pickle.dumps(sambam1.header)
            sambam2 = pysam.AlignmentFile(sambam_file2, 'rb')  # @UndefinedVariable
            self._reads._v_attrs.header2 = pickle.dumps(sambam2.header)
            
            
            
            # index node table
            logging.info("Creating index...")
            try:
                self._reads.cols.qname.create_index()
            except ValueError:
                # Index exists, no problem!
                pass
            
            # map reads
            logging.info("Loading left reads...")
            self._load_left(sambam1, check_duplicates=False)
            logging.info("Loading right reads...")
            self._load_right(sambam2, check_duplicates=True)
            logging.info("Done.")
            
        
    def _determine_field_sizes(self, sambam):
        if type(sambam) == str:
            sambam = pysam.AlignmentFile(sambam, 'rb')  # @UndefinedVariable
            
        qname_length = 0
        seq_length = 0
        cigar_length = 0
        tags_length = 0
        for r in sambam:
            qname_length = max(qname_length,len(r.qname))
            seq_length = max(seq_length,len(r.seq))
            cigar_length = max(cigar_length,len(r.cigarstring))
            tags_length = max(tags_length,len(pickle.dumps(r.tags)))
        sambam.close()
        
        return { 'qname': qname_length,
                 'sequence': seq_length,
                 'tags': tags_length,
                 'cigar': cigar_length }
    
    
    def _load_left(self, sambam, check_duplicates=True):
        if type(sambam) == str:
            sambam = pysam.AlignmentFile(sambam, 'rb')  # @UndefinedVariable
        
        i = 0
        for r in sambam:
            if i % 100000 == 0:
                logging.info("\t%d" % i)
                self._reads.flush()
            qname = r.qname
            
            if check_duplicates:
                res = [x for x in self._reads.where("qname == '%s'" % qname)]
                if len(res) >= 1:
                    raise ValueError("Duplicate read QNAME %s" % qname)
            
            row = self._reads.row
            row['qname'] = qname
            row['flag1'] = r.flag
            row['ref1'] = r.reference_id
            row['pos1'] = r.pos
            row['mapq1'] = r.mapq
            row['cigar1'] = r.cigarstring
            row['rnext1'] = r.rnext
            row['pnext1'] =  r.pnext
            row['tlen1'] = r.tlen
            row['seq1'] = r.seq
            row['qual1'] = r.qual
            row['tags1'] = pickle.dumps(r.tags)
            row.append()
            i += 1
        self._reads.flush()
        
    def load_right(self, sambam, check_duplicates=True):
        if type(sambam) == str:
            sambam = pysam.AlignmentFile(sambam, 'rb')  # @UndefinedVariable
        
        i = 0
        for r in sambam:
            if i % 100000 == 0:
                logging.info("\t%d" % i)
                self._reads.flush()
                
            qname = r.qname
            
            row = self._reads.row
            update = False
            if check_duplicates:
                res = self.find_read_row(qname)
                if len(res) > 1:
                    raise ValueError("Duplicate read QNAME %s" % qname)
                if res is not None:
                    update = True
                    row = res
             
            row['qname'] = qname
            row['flag2'] = r.flag
            row['ref2'] = r.reference_id
            row['pos2'] = r.pos
            row['mapq2'] = r.mapq
            row['cigar2'] = r.cigarstring
            row['rnext2'] = r.rnext
            row['pnext2'] =  r.pnext
            row['tlen2'] = r.tlen
            row['seq2'] = r.seq
            row['qual2'] = r.qual
            row['tags2'] = pickle.dumps(r.tags)
                
            if update:
                row.update()
            else:
                row.append()
            i += 1
        self._reads.flush()
                
                
                
                
    def find_read_row(self, name):
        res = [x for x in self._reads.where("qname == '%s'" % name)]
        if len(res) > 1:
            raise ValueError("Multiple reads with QNAME %s" % name)
        if len(res) == 1:
            return res[0]
        return None
        