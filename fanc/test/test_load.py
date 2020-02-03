import kaic
from kaic.registry import class_name_dict
import genomic_regions as gr
import pytest
import pysam
import os


class TestAuto:
    def test_auto_identification(self, tmpdir):
        for class_name in ('Hic', 'LegacyHic', 'ReadPairs',
                           'AggregateMatrix', 'ComparisonMatrix',
                           'FoldChangeMatrix', 'DifferenceMatrix',
                           'ComparisonRegions', 'FoldChangeRegions',
                           'DifferenceRegions', 'DirectionalityIndex'
                           ):
            file_name = str(tmpdir) + '/{}.h5'.format(class_name)
            cls_ = class_name_dict[class_name]
            x = cls_(file_name=file_name, mode='w')
            x.close()

            x = kaic.load(file_name, mode='r')
            assert isinstance(x, cls_)
            x.close()

    # def test_hic_based_auto_identification(self, tmpdir):
    #     with dummy.sample_hic() as hic:
    #         for class_name in ('AggregateMatrix', 'ComparisonMatrix',
    #                            'FoldChangeMatrix', 'DifferenceMatrix',
    #                            'ComparisonScores', 'FoldChangeScores',
    #                            'DifferenceScores',
    #                            'ComparisonRegions', 'FoldChangeRegions',
    #                            'DifferenceRegions', 'DirectionalityIndex'):
    #             file_name = str(tmpdir) + '/{}.h5'.format(class_name)
    #             cls_ = class_name_dict[class_name]
    #             x = cls_(hic, file_name=file_name, mode='w')
    #             x.close()
    #
    #             x = kaic.load(file_name, mode='r')
    #             assert isinstance(x, cls_)
    #             x.close()
    #         for class_name in ('FoldChangeMatrix',):
    #             file_name = str(tmpdir) + '/{}.h5'.format(class_name)
    #             cls_ = class_name_dict[class_name]
    #             x = cls_(hic, hic, file_name=file_name, mode='w')
    #             x.close()
    #
    #             x = kaic.load(file_name, mode='r')
    #             assert isinstance(x, cls_)
    #             x.close()
    #
    # def test_conversion(self, tmpdir):
    #     file_name = str(tmpdir) + '/x.hic'
    #     with dummy.sample_hic(file_name=file_name) as hic:
    #         # simulate old-style object
    #         hic.file.remove_node('/meta_information', recursive=True)
    #
    #     hic = kaic.load(file_name, mode='r')
    #     assert isinstance(hic, kaic.Hic)
    #     hic.close()
    #
    #     hic = kaic.Hic(file_name)
    #     hic.close()
    #
    #     hic = kaic.load(file_name, mode='r')
    #     hic.close()
    #     assert isinstance(hic, kaic.Hic)
    #
    # def test_old_style_index(self, tmpdir):
    #     with dummy.sample_hic() as hic:
    #         for class_name in ('ABDomains', 'ABDomainMatrix', 'ExpectedContacts', 'ObservedExpectedRatio',
    #                            'ABDomains', 'PossibleContacts', 'RegionContactAverage',
    #                            'InsulationIndex', 'DirectionalityIndex'):
    #             file_name = str(tmpdir) + '/{}.h5'.format(class_name)
    #             cls_ = class_name_dict[class_name]
    #             x = cls_(hic, file_name=file_name, mode='w')
    #             # simulate missing meta-information
    #             x.close()
    #
    #             x = kaic.load(file_name, mode='r')
    #             assert isinstance(x, cls_)
    #             x.close()
    #
    #         for class_name in ('FoldChangeMatrix',):
    #             file_name = str(tmpdir) + '/{}.h5'.format(class_name)
    #             cls_ = class_name_dict[class_name]
    #             x = cls_(hic, hic, file_name=file_name, mode='w')
    #             # simulate missing meta-information
    #             x.close()
    #
    #             x = kaic.load(file_name, mode='r')
    #             assert isinstance(x, cls_)
    #             x.close()
    #
    #     for class_name in ('Hic', 'AccessOptimisedHic', 'FragmentMappedReadPairs', 'Reads', 'GenomicTrack'):
    #         file_name = str(tmpdir) + '/{}.h5'.format(class_name)
    #         cls_ = class_name_dict[class_name]
    #         x = cls_(file_name=file_name, mode='w')
    #         # simulate missing meta-information
    #         x.file.remove_node('/meta_information', recursive=True)
    #         x.close()
    #
    #         x = kaic.load(file_name, mode='r')
    #         assert isinstance(x, cls_)
    #         x.close()

    def test_bed(self):
        this_dir = os.path.dirname(os.path.realpath(__file__))
        bed_file = this_dir + '/test_load/test.bed'

        with kaic.load(bed_file) as bed:
            assert isinstance(bed, gr.Bed)

        with pytest.raises(ValueError):
            foo_file = this_dir + '/test_load/foo.txt'
            kaic.load(foo_file)

    def test_bigwig(self):
        this_dir = os.path.dirname(os.path.realpath(__file__))
        bw_file = this_dir + '/test_load/test.bw'

        with kaic.load(bw_file) as bw:
            assert isinstance(bw, gr.BigWig)

    def test_sambam(self):
        this_dir = os.path.dirname(os.path.realpath(__file__))
        sam_file = this_dir + '/test_load/test.sam'

        with kaic.load(sam_file, mode='r') as bw:
            assert isinstance(bw, pysam.AlignmentFile)

        bam_file = this_dir + '/test_load/test.bam'

        with kaic.load(bam_file, mode='r') as bw:
            assert isinstance(bw, pysam.AlignmentFile)
