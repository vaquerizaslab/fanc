'''
Created on Jul 15, 2015

@author: kkruse1
'''

import os.path
import pickle
from numpy import array_equal
from kaic.construct.seq import ReadPairs

class TestReadPairs:
    
    @classmethod
    def setup_method(self, method):
        self.dir = os.path.dirname(os.path.realpath(__file__))
        self.sam1_file = self.dir + "/test_seq/test1.sam"
        self.sam2_file = self.dir + "/test_seq/test2.sam"
#         self.pairs = ReadPairs(self.sam1_file,self.sam2_file)
    
    def test_load_from_sam_file(self):
        def compare(read, values):
            assert read.qname == values[0]
            assert read.flag == values[1]
            assert read.ref == values[2]
            assert read.pos == values[3]
            assert read.mapq == values[4]
            assert read.cigar == values[5]
            assert read.rnext == values[6]
            assert read.pnext == values[7]
            assert read.tlen == values[8]
            assert read.seq == values[9]
            assert read.qual == values[10]
            assert len(read.tags) == values[11]

        pairs = ReadPairs()
        pairs.load(self.sam1_file,self.sam2_file, is_sorted=False)
        
        # SRR038105.1    0    chrXI    128390    35    15M    *    0    0    GATATGATGGATTTG    FDFFFFFFFFFFFCF    AS:i:0    XS:i:-5    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:15    YT:Z:UU
        read_pair1 = pairs[0]
        compare(read_pair1.left_read, ['SRR038105.1',0,'chrXI',128390,35,'15M',-1,-1,0,'GATATGATGGATTTG','FDFFFFFFFFFFFCF',9])
        compare(read_pair1.right_read, ['SRR038105.1',16,'chrIX',387161,40,'19M',-1,-1,0,'GAAAAAAAAAAAGAAGGGC','..A70AA*;A;A.;A4A*A',8])
        
        
        # SRR038105.1000167    0    chrXIV    703158    42    15M    *    0    0    TACGGTATTGGTCGG    FFFFCFFFFFFFFCF    AS:i:0    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:15    YT:Z:UU
        res = pairs.where("qname == 'SRR038105.1000167'")
        read_pair1 = res[0]
        compare(read_pair1.left_read, ['SRR038105.1000167',0,'chrXIV',703158,42,'15M',-1,-1,0,'TACGGTATTGGTCGG','FFFFCFFFFFFFFCF',8])
        compare(read_pair1.right_read, ['SRR038105.1000167',0,'chrIII',130079,35,'16M',-1,-1,0,'CATTTTATATGAATTA','IIIIIIIDIIIIDIII',9])
        
        
        # SRR038105.1000320    0    chrXVI    577162    35    15M    *    0    0    TTGATAAAATAGTCC    <<@FF<FFFFAFAFA    AS:i:0    XS:i:-5    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:15    YT:Z:UU
        res = pairs.where("qname == 'SRR038105.1000320'")
        read_pair1 = res[0]
        compare(read_pair1.left_read, ['SRR038105.1000320',0,'chrXVI',577162,35,'15M',-1,-1,0,'TTGATAAAATAGTCC','<<@FF<FFFFAFAFA',9])
        # SRR038105.1000320    16    chrXVI    576849    42    15M    *    0    0    CCAACAGAGTACACT    9>7BC???B=<AC<.    AS:i:0    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:15    YT:Z:UU
        compare(read_pair1.right_read,['SRR038105.1000320',16,'chrXVI',576849,42,'15M',-1,-1,0,'CCAACAGAGTACACT','9>7BC???B=<AC<.',8])
        

        
        # check unpaired right
        # SRR038105.1000002    16    chrIV    203242    42    16M    *    0    0    ACCCATTATTTCTCGA    IIIIIFIICIFIIIII    AS:i:0    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:16    YT:Z:UU
        read_pair1 = pairs.where("qname == 'SRR038105.1000002'")[0]
        assert read_pair1.left_read is None
        compare(read_pair1.right_read, ['SRR038105.1000002',16,'chrIV',203242,42,'16M',-1,-1,0,'ACCCATTATTTCTCGA','IIIIIFIICIFIIIII',8])
        
        # check unpaired left
        # SRR038105.1000011    16    chrIV    526796    42    16M    *    0    0    GGTGAATTAGAAGATA    FFFFFFFFFFFFFFFF    AS:i:0    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:16    YT:Z:UU
        read_pair1 = pairs.where("qname == 'SRR038105.1000011'")[0]
        compare(read_pair1.left_read, ['SRR038105.1000011',16,'chrIV',526796,42,'16M',-1,-1,0,'GGTGAATTAGAAGATA','FFFFFFFFFFFFFFFF',8])
        assert read_pair1.right_read is None
        
    
    def test_ix(self):
        pairs = ReadPairs(self.sam1_file,self.sam2_file)
        i = 0
        for pair in pairs._reads:
            assert pair['ix'] == i
            i += 1
    
    def test_iter(self):
        pairs = ReadPairs(self.sam1_file,self.sam2_file)
        counter = 0
        for _ in pairs:
            counter += 1
        
        assert counter == 315
        
        pairs.filter_non_unique()
        after_counter = 0
        for _ in pairs:
            after_counter += 1
            
        assert after_counter < counter
    
    def test_select(self):
        pairs = ReadPairs(self.sam1_file,self.sam2_file)
        
        assert pairs[0].qname == 'SRR038105.1'
        pairs.filter_non_unique()
        assert pairs[0].qname == 'SRR038105.1000'
        
        
    def test_build_from_scratch(self):
        
        field_sizes = ReadPairs.determine_field_sizes(self.sam1_file, sample_size=10000)
        pairs = ReadPairs(field_sizes=field_sizes)
        pairs.load(self.sam1_file,self.sam2_file)
        
    def test_quality_filter(self):
        pairs = ReadPairs()
        pairs.load(self.sam1_file,self.sam2_file)
        
        pairs.filter_quality(30, queue=False)
        for row in pairs._reads:
            if (row['pos1'] > 0 and row['mapq1'] < 30) or (row['pos2'] > 0 and row['mapq2'] < 30):
                assert row[pairs._reads._mask_field] == 2
            else:
                assert row[pairs._reads._mask_field] == 0
        
    def test_uniqueness_filter(self):
        pairs = ReadPairs()
        pairs.load(self.sam1_file,self.sam2_file)
        
        pairs.filter_non_unique(strict=True)
        for row in pairs._reads:
            if row['pos1'] > 0:
                tags = pairs._tags_left[row['ix']]
                has_xs = False
                for tag in tags:
                    if tag[0] == 'XS':
                        has_xs = True
                        break
                if has_xs:
                    assert row[pairs._reads._mask_field] == 2
            elif row['pos2'] > 0:
                tags = pairs._tags_right[row['ix']]
                has_xs = False
                for tag in tags:
                    if tag[0] == 'XS':
                        has_xs = True
                        break
                if has_xs:
                    assert row[pairs._reads._mask_field] == 2
            else:
                assert row[pairs._reads._mask_field] == 0
                
    def test_single_filter(self):
        pairs = ReadPairs()
        pairs.load(self.sam1_file,self.sam2_file)
        
        pairs.filter_single()
        for row in pairs._reads:
            if row['pos1'] == 0  or row['pos2'] == 0:
                assert row[pairs._reads._mask_field] == 2
            else:
                assert row[pairs._reads._mask_field] == 0
                
    def test_queue_filters(self):
        pairs = ReadPairs()
        pairs.load(self.sam1_file,self.sam2_file)
        
        l = len(pairs)
        
        pairs.filter_single(queue=True)
        pairs.filter_quality(30, queue=True)
        pairs.filter_non_unique(strict=True, queue=True)
        
        assert len(pairs) == l
        
        pairs.run_queued_filters()
        
        assert len(pairs) < l
        
        

        
        