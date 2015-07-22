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
        pairs = ReadPairs()
        pairs.load(self.sam1_file,self.sam2_file, is_sorted=False)
        
#         assert len(pairs) == 324
        
        # SRR038105.1    0    chrXI    128390    35    15M    *    0    0    GATATGATGGATTTG    FDFFFFFFFFFFFCF    AS:i:0    XS:i:-5    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:15    YT:Z:UU
        read_pair1 = pairs._reads[0]
        assert read_pair1[0] == 'SRR038105.1'
        assert read_pair1[1] == 0
        assert read_pair1[2] == 12
        assert pairs.ix2ref1(read_pair1[2]) == 'chrXI'
        assert read_pair1[3] == 128390
        assert read_pair1[4] == 35
        assert read_pair1[5] == '15M'
        assert read_pair1[6] == -1
        assert read_pair1[7] == -1
        assert read_pair1[8] == 0
        assert read_pair1[9] == 'GATATGATGGATTTG'
        assert read_pair1[10] == 'FDFFFFFFFFFFFCF'
        tags1 = pickle.loads(read_pair1[11])
        assert len(tags1) == 9
        
        # SRR038105.1    16    chrIX    387161    40    19M    *    0    0    GAAAAAAAAAAAGAAGGGC    ..A70AA*;A;A.;A4A*A    AS:i:-3    XN:i:0    XM:i:1    XO:i:0    XG:i:0    NM:i:1    MD:Z:0A18    YT:Z:UU
        assert read_pair1[12] == 16
        assert read_pair1[13] == 4
        assert pairs.ix2ref2(read_pair1[13]) == 'chrIX'
        assert read_pair1[14] == 387161
        assert read_pair1[15] == 40
        assert read_pair1[16] == '19M'
        assert read_pair1[17] == -1
        assert read_pair1[18] == -1
        assert read_pair1[19] == 0
        assert read_pair1[20] == 'GAAAAAAAAAAAGAAGGGC'
        assert read_pair1[21] == '..A70AA*;A;A.;A4A*A'
        tags1 = pickle.loads(read_pair1[22])
        assert len(tags1) == 8
        
        # SRR038105.1000167    0    chrXIV    703158    42    15M    *    0    0    TACGGTATTGGTCGG    FFFFCFFFFFFFFCF    AS:i:0    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:15    YT:Z:UU
        res = [x.fetch_all_fields() for x in pairs._reads.where("qname == 'SRR038105.1000167'")]
        read_pair1 = res[0]
        assert read_pair1[0] == 'SRR038105.1000167'
        assert read_pair1[1] == 0
        assert read_pair1[2] == 13
        assert pairs.ix2ref1(read_pair1[2]) == 'chrXIV'
        assert read_pair1[3] == 703158
        assert read_pair1[4] == 42
        assert read_pair1[5] == '15M'
        assert read_pair1[6] == -1
        assert read_pair1[7] == -1
        assert read_pair1[8] == 0
        assert read_pair1[9] == 'TACGGTATTGGTCGG'
        assert read_pair1[10] == 'FFFFCFFFFFFFFCF'
        tags1 = pickle.loads(read_pair1[11])
        assert len(tags1) == 8
        
        # SRR038105.1000167    0    chrIII    130079    35    16M    *    0    0    CATTTTATATGAATTA    IIIIIIIDIIIIDIII    AS:i:0    XS:i:-6    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:16    YT:Z:UU
        assert read_pair1[12] == 0
        assert read_pair1[13] == 0
        assert pairs.ix2ref2(read_pair1[13]) == 'chrIII'
        assert read_pair1[14] == 130079
        assert read_pair1[15] == 35
        assert read_pair1[16] == '16M'
        assert read_pair1[17] == -1
        assert read_pair1[18] == -1
        assert read_pair1[19] == 0
        assert read_pair1[20] == 'CATTTTATATGAATTA'
        assert read_pair1[21] == 'IIIIIIIDIIIIDIII'
        tags1 = pickle.loads(read_pair1[22])
        assert len(tags1) == 9
        
        # SRR038105.1000320    0    chrXVI    577162    35    15M    *    0    0    TTGATAAAATAGTCC    <<@FF<FFFFAFAFA    AS:i:0    XS:i:-5    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:15    YT:Z:UU
        res = [x.fetch_all_fields() for x in pairs._reads.where("qname == 'SRR038105.1000320'")]
        read_pair1 = res[0]
        assert read_pair1[0] == 'SRR038105.1000320'
        assert read_pair1[1] == 0
        assert read_pair1[2] == 15
        assert pairs.ix2ref1(read_pair1[2]) == 'chrXVI'
        assert read_pair1[3] == 577162
        assert read_pair1[4] == 35
        assert read_pair1[5] == '15M'
        assert read_pair1[6] == -1
        assert read_pair1[7] == -1
        assert read_pair1[8] == 0
        assert read_pair1[9] == 'TTGATAAAATAGTCC'
        assert read_pair1[10] == '<<@FF<FFFFAFAFA'
        tags1 = pickle.loads(read_pair1[11])
        assert len(tags1) == 9
        
        # SRR038105.1000320    16    chrXVI    576849    42    15M    *    0    0    CCAACAGAGTACACT    9>7BC???B=<AC<.    AS:i:0    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:15    YT:Z:UU
        assert read_pair1[12] == 16
        assert read_pair1[13] == 15
        assert pairs.ix2ref2(read_pair1[13]) == 'chrXVI'
        assert read_pair1[14] == 576849
        assert read_pair1[15] == 42
        assert read_pair1[16] == '15M'
        assert read_pair1[17] == -1
        assert read_pair1[18] == -1
        assert read_pair1[19] == 0
        assert read_pair1[20] == 'CCAACAGAGTACACT'
        assert read_pair1[21] == '9>7BC???B=<AC<.'
        tags1 = pickle.loads(read_pair1[22])
        assert len(tags1) == 8
        
        # check unpaired right
        # SRR038105.1000002    16    chrIV    203242    42    16M    *    0    0    ACCCATTATTTCTCGA    IIIIIFIICIFIIIII    AS:i:0    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:16    YT:Z:UU
        read_pair1 = [x.fetch_all_fields() for x in pairs._reads.where("qname == 'SRR038105.1000002'")][0]
        assert read_pair1[0] == 'SRR038105.1000002'
        assert array_equal(tuple(read_pair1)[1:11],[0,0,0,0,'',0,0,0,'',''])
        assert array_equal(tuple(read_pair1)[12:22],[16,3,203242,42,'16M',-1,-1,0,'ACCCATTATTTCTCGA','IIIIIFIICIFIIIII'])
        
        # check unpaired left
        # SRR038105.1000011    16    chrIV    526796    42    16M    *    0    0    GGTGAATTAGAAGATA    FFFFFFFFFFFFFFFF    AS:i:0    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:16    YT:Z:UU
        read_pair1 = [x.fetch_all_fields() for x in pairs._reads.where("qname == 'SRR038105.1000011'")][0]
        assert read_pair1[0] == 'SRR038105.1000011'
        assert array_equal(tuple(read_pair1)[1:11],[16,3,526796,42,'16M',-1,-1,0,'GGTGAATTAGAAGATA','FFFFFFFFFFFFFFFF'])
        assert array_equal(tuple(read_pair1)[12:22],[0,0,0,0,'',0,0,0,'',''])
        
        
    def test_iter(self):
        pairs = ReadPairs(self.sam1_file,self.sam2_file)
        counter = 0
        for pair in pairs:
            counter += 1
        
        assert counter == 315
        
        pairs.filter_non_unique()
        after_counter = 0
        for pair in pairs:
            after_counter += 1
            
        assert after_counter < counter
    
    def test_select(self):
        pairs = ReadPairs(self.sam1_file,self.sam2_file)
        
        for pair in pairs:
            print pair[24]
        
        assert pairs[0][0] == 'SRR038105.1'
        pairs.filter_non_unique()
        assert pairs[0][0] == 'SRR038105.1000'
        
        
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
                assert row['mask'] == 2
            else:
                assert row['mask'] == 0
        
    def test_uniqueness_filter(self):
        pairs = ReadPairs()
        pairs.load(self.sam1_file,self.sam2_file)
        
        pairs.filter_non_unique(strict=True)
        for row in pairs._reads:
            if row['pos1'] > 0:
                tags = pickle.loads(row['tags1'])
                has_xs = False
                for tag in tags:
                    if tag[0] == 'XS':
                        has_xs = True
                        break
                if has_xs:
                    assert row['mask'] == 2
            elif row['pos2'] > 0:
                tags = pickle.loads(row['tags2'])
                has_xs = False
                for tag in tags:
                    if tag[0] == 'XS':
                        has_xs = True
                        break
                if has_xs:
                    assert row['mask'] == 2
            else:
                assert row['mask'] == 0
                
    def test_single_filter(self):
        pairs = ReadPairs()
        pairs.load(self.sam1_file,self.sam2_file)
        
        pairs.filter_single()
        for row in pairs._reads:
            if row['pos1'] == 0  or row['pos2'] == 0:
                assert row['mask'] == 2
            else:
                assert row['mask'] == 0
                
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
        
        

        
        