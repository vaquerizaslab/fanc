'''
Created on Jul 15, 2015

@author: kkruse1
'''

import os.path
from kaic.construct.seq import Reads, FragmentMappedReadPairs,\
    FragmentRead, InwardPairsFilter, UnmappedFilter, OutwardPairsFilter,\
    ReDistanceFilter, FragmentReadPair
from kaic.data.genomic import Genome, GenomicRegion, Chromosome


class TestReads:
    
    @classmethod
    def setup_method(self, method):
        self.dir = os.path.dirname(os.path.realpath(__file__))
        self.sam1_file = self.dir + "/test_seq/test1.sam"
        self.sam2_file = self.dir + "/test_seq/test2.sam"
        self.lambda_sam1_file = self.dir + "/test_seq/lambda_reads1.sam"
        self.bwamem_sam1_file = self.dir + "/test_seq/test_bwa1.sam"
        self.bwamem_sam2_file = self.dir + "/test_seq/test_bwa2.sam"
        
    def test_load(self):
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
            assert read.strand == values[12]
        
        reads = Reads()
        reads.load(self.sam1_file, is_sorted=False)
        
        compare(reads[0], ['SRR038105.1',0,'chrXI',128390,35,'15M',-1,-1,0,'GATATGATGGATTTG','FDFFFFFFFFFFFCF',9,1])
        
        # SRR038105.1000167    0    chrXIV    703158    42    15M    *    0    0    TACGGTATTGGTCGG    FFFFCFFFFFFFFCF    AS:i:0    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:15    YT:Z:UU
        res = reads.where("qname == 'SRR038105.1000167'")
        compare(res[0], ['SRR038105.1000167',0,'chrXIV',703158,42,'15M',-1,-1,0,'TACGGTATTGGTCGG','FFFFCFFFFFFFFCF',8,1])        
        
        # SRR038105.1000320    0    chrXVI    577162    35    15M    *    0    0    TTGATAAAATAGTCC    <<@FF<FFFFAFAFA    AS:i:0    XS:i:-5    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:15    YT:Z:UU
        res = reads.where("qname == 'SRR038105.1000320'")
        compare(res[0], ['SRR038105.1000320',0,'chrXVI',577162,35,'15M',-1,-1,0,'TTGATAAAATAGTCC','<<@FF<FFFFAFAFA',9,1])

        # check unpaired right
        # SRR038105.1000002    16    chrIV    203242    42    16M    *    0    0    ACCCATTATTTCTCGA    IIIIIFIICIFIIIII    AS:i:0    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:16    YT:Z:UU
        res = reads.where("qname == 'SRR038105.1000002'")
        assert len(res) == 0
        
        # check unpaired left
        # SRR038105.1000011    16    chrIV    526796    42    16M    *    0    0    GGTGAATTAGAAGATA    FFFFFFFFFFFFFFFF    AS:i:0    XN:i:0    XM:i:0    XO:i:0    XG:i:0    NM:i:0    MD:Z:16    YT:Z:UU
        res = reads.where("qname == 'SRR038105.1000011'")
        compare(res[0], ['SRR038105.1000011',16,'chrIV',526796,42,'16M',-1,-1,0,'GGTGAATTAGAAGATA','FFFFFFFFFFFFFFFF',8,-1])
        
    def test_ix(self):
        pairs = Reads(self.sam1_file)
        i = 0
        for pair in pairs._reads:
            assert pair['ix'] == i
            i += 1
    
    def test_strand(self):
        reads = Reads(self.sam1_file)
        for read in reads:
            if read.flag == 16:
                assert read.strand == -1
            if read.flag == 4:
                assert read.strand == 1
    
    def test_iter(self):
        pairs = Reads(self.sam1_file)
        counter = 0
        for _ in pairs:
            counter += 1
        
        assert counter == 271
        
        pairs.filter_non_unique()
        after_counter = 0
        for _ in pairs:
            after_counter += 1
            
        assert after_counter < counter
    
    def test_select(self):
        reads = Reads(self.sam1_file)
        
        assert reads[0].qname == 'SRR038105.1'
        reads.filter_non_unique()
        assert reads[0].qname == 'SRR038105.10'
        
    def test_build_from_scratch(self):
        
        qname_length, seq_length = Reads.determine_field_sizes(self.sam1_file, sample_size=10000)
        reads = Reads(qname_length=qname_length, seq_length=seq_length)
        reads.load(self.sam1_file)
        
        assert len(reads) == 271

    def test_infer_mapper(self):
        reads = Reads(self.bwamem_sam1_file)
        assert reads._mapper == 'bwa'
        reads = Reads(self.bwamem_sam1_file, mapper='bowtie2')
        assert reads._mapper == 'bowtie2'

    def test_parse_cigar(self):
        reads = Reads()
        parsed = reads.parse_cigar('22S42M1D18M4I5M')
        assert parsed ==[('S',22), ('M',42), ('D',1), ('M',18), ('I',4), ('M',5)]

    def test_parse_read_cigar(self):
        reads = Reads(self.sam1_file)
        assert reads[0].get_cigar == [('M',15)]
        assert reads[1].get_cigar == [('M',20)]

    def test_read_alen(self):
        reads = Reads(self.sam1_file)
        read = reads[0]
        read.cigar = '22S42M1D18M4I5M'
        assert read.alen == 42 + 18 + 5

    def test_read_get_tag(self):
        reads = Reads(self.sam1_file)
        read = reads[0]
        assert read.get_tag('AS') == 0
        assert read.get_tag('MD') == '15'
        assert read.get_tag('X0') == None
    
    def test_quality_filter(self):
        reads = Reads(self.sam1_file)
        
        reads.filter_quality(30, queue=False)
        for row in reads._reads.all():
            if row['mapq'] < 30:
                assert row[reads._reads._mask_field] == 2
            else:
                assert row[reads._reads._mask_field] == 0
        
    def test_uniqueness_filter(self):
        reads = Reads(self.sam1_file)
        
        reads.filter_non_unique(strict=True)
        for row in reads._reads.all():
            if row['pos'] > 0:
                tags = reads._tags[row['ix']]
                has_xs = False
                for tag in tags:
                    if tag[0] == 'XS':
                        has_xs = True
                        break
                if has_xs:
                    assert row[reads._reads._mask_field] == 2
            else:
                assert row[reads._reads._mask_field] == 0
    
    def test_unmapped_filter(self):
        reads = Reads(self.lambda_sam1_file)
        
        l = len(reads)
        
        unmapped_filter = UnmappedFilter(reads.add_mask_description("unmapped", "Filter unmapped reads"))
        reads.filter(unmapped_filter)
        
        assert len(reads) < l

    def test_bwamem_quality_filter(self):
        reads = Reads(self.bwamem_sam1_file)
        assert len(reads) == 919
        reads.filter_quality(0.90, queue=False)
        assert len(reads) == 809
        for read in reads:
            assert float(read.get_tag('AS')) / read.alen >= 0.90
        assert len(reads.cache_maker._cache['parse_cigar'].data) == 118

    def test_queue_filters(self):
        reads = Reads(self.sam1_file)
        
        l = len(reads)
        
        reads.filter_quality(30, queue=True)
        reads.filter_non_unique(strict=True, queue=True)
        
        assert len(reads) == l
        
        reads.run_queued_filters()
        
        assert len(reads) < l


class TestFragmentMappedReads:
    @classmethod
    def setup_method(self, method):
        self.dir = os.path.dirname(os.path.realpath(__file__))
        sam1_file = self.dir + "/test_seq/lambda_reads1.sam"
        sam2_file = self.dir + "/test_seq/lambda_reads2.sam"
        self.reads1 = Reads(sam1_file)
        self.reads2 = Reads(sam2_file)
        self.reads1.filter_unmapped()
        self.reads2.filter_unmapped()
        self.genome = Genome.from_folder(self.dir + "/test_seq/lambda_genome/")
        
        self.pairs = FragmentMappedReadPairs()
        self.pairs.load(self.reads1, self.reads2, regions=self.genome.get_regions(1000))
        
    
    def test_select(self):
        pair = self.pairs[0]
        assert isinstance(pair, FragmentReadPair)
        assert pair.left.position == 18401
        assert pair.right.position == 18430
        assert pair.left.strand == 1
        assert pair.right.strand == -1
        assert isinstance(pair.left.fragment, GenomicRegion)
        assert isinstance(pair.right.fragment, GenomicRegion)
        assert pair.left.fragment.start == 18001
        assert pair.left.fragment.end == 19000
        assert pair.left.fragment.chromosome == 'gi|9626243|ref|NC_001416.1|'
        assert pair.right.fragment.start == 18001
        assert pair.right.fragment.end == 19000
        assert pair.right.fragment.chromosome == 'gi|9626243|ref|NC_001416.1|'
        
        pair = self.pairs[1]
        assert isinstance(pair, FragmentReadPair)
        assert pair.left.position == 40075
        assert pair.right.position == 40211
        assert pair.left.strand == 1
        assert pair.right.strand == -1
        assert isinstance(pair.left.fragment, GenomicRegion)
        assert isinstance(pair.right.fragment, GenomicRegion)
        assert pair.left.fragment.start == 40001
        assert pair.left.fragment.end == 41000
        assert pair.left.fragment.chromosome == 'gi|9626243|ref|NC_001416.1|'
        assert pair.right.fragment.start == 40001
        assert pair.right.fragment.end == 41000
        assert pair.right.fragment.chromosome == 'gi|9626243|ref|NC_001416.1|'
        
        pair = self.pairs[-1]
        assert isinstance(pair, FragmentReadPair)
        assert pair.left.position == 5067
        assert pair.right.position == 5200
        assert pair.left.strand == 1
        assert pair.right.strand == -1
        assert isinstance(pair.left.fragment, GenomicRegion)
        assert isinstance(pair.right.fragment, GenomicRegion)
        assert pair.left.fragment.start == 5001
        assert pair.left.fragment.end == 6000
        assert pair.left.fragment.chromosome == 'gi|9626243|ref|NC_001416.1|'
        assert pair.right.fragment.start == 5001
        assert pair.right.fragment.end == 6000
        assert pair.right.fragment.chromosome == 'gi|9626243|ref|NC_001416.1|'
    
    def test_iter(self):
        for pair in self.pairs:
            assert isinstance(pair, FragmentReadPair)
            assert isinstance(pair.left, FragmentRead)
            assert isinstance(pair.right, FragmentRead)
            assert isinstance(pair.left.fragment, GenomicRegion)
            assert isinstance(pair.right.fragment, GenomicRegion)
            assert pair.left.position > 0 or pair.left.position == -1
            assert pair.right.position > 0 or pair.right.position == -1
            assert pair.left.strand == -1 or pair.left.strand == 1
            assert pair.right.strand == -1 or pair.right.strand == 1
            
            assert 0 < pair.left.fragment.start <= pair.left.fragment.end
            assert 0 < pair.right.fragment.start <= pair.right.fragment.end
            if pair.left.position > 0:
                assert pair.left.fragment.start <= pair.left.position <= pair.left.fragment.end
            if pair.right.position > 0:
                assert pair.right.fragment.start <= pair.right.position <= pair.right.fragment.end
            
    def test_len(self):
        assert len(self.pairs) == 44
            
    def test_single(self):
        assert len(self.pairs._single) == 6

    def test_filter_inward(self):
        mask = self.pairs.add_mask_description('inwards', 'Mask read pairs that inward facing and closer than 100bp')
        in_filter = InwardPairsFilter(minimum_distance=100, mask=mask)
        
        assert len(self.pairs) == 44
        self.pairs.filter(in_filter)
        
#         print "Valid pairs:"
#         for pair in self.pairs:
#             print pair[0]
#             print pair[1]
        assert len(self.pairs) == 18
        
    def test_filter_outward(self):
        mask = self.pairs.add_mask_description('outwards', 'Mask read pairs that outward facing and closer than 100bp')
        out_filter = OutwardPairsFilter(minimum_distance=100, mask=mask)
        
        assert len(self.pairs) == 44
        self.pairs.filter(out_filter)
        assert len(self.pairs) == 28
        
    def test_filter_redist(self):
        mask = self.pairs.add_mask_description('re-dist', 'Mask read pairs where one half maps more than 100bp away from both RE sites')
        re_filter = ReDistanceFilter(maximum_distance=300, mask=mask)
        
        assert len(self.pairs) == 44
        self.pairs.filter(re_filter)
        assert len(self.pairs) == 13

    def test_get_error_structure(self):
        reads1 = Reads(self.dir + "/../data/test_genomic/yeast.sample.chrI.1.sam")
        reads2 = Reads(self.dir + "/../data/test_genomic/yeast.sample.chrI.2.sam")
        chrI = Chromosome.from_fasta(self.dir + "/../data/test_genomic/chrI.fa")
        genome = Genome(chromosomes=[chrI])
        pairs = FragmentMappedReadPairs()
        pairs.load(reads1, reads2, genome.get_regions('HindIII'))
        x, i, o = pairs.get_error_structure(data_points=200, skip_self_ligations=False)
        assert len(x) == len(i) == len(o) == 3
        assert x == [494.03856041131104, 4487.5800970873788, 19399.908018867925]
        assert i == [2.616915422885572, 0.8059701492537313, 0.6417910447761194]
        assert o == [0.2537313432835821, 0.24378109452736318, 0.46766169154228854]

    def test_re_dist(self):
        read1 = FragmentRead(GenomicRegion(chromosome='chr1', start=1, end=1000), position=200, strand=-1)
        assert read1.re_distance() == 199
        read2 = FragmentRead(GenomicRegion(chromosome='chr1', start=1, end=1000), position=990, strand=-1)
        assert read2.re_distance() == 10
