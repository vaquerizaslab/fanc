import kaic
from kaic.data.genomic import Hic
from kaic.data.registry import class_name_dict

class TestAuto:
    def test_auto_identification(self, tmpdir):
        for class_name in ('Hic', 'AccessOptimisedHic', 'FragmentMappedReadPairs', 'Reads'):
            file_name = str(tmpdir) + '/{}.h5'.format(class_name)
            cls_ = class_name_dict[class_name]
            x = cls_(file_name=file_name, mode='w')
            x.close()

            x = kaic.load(file_name, mode='r')
            assert isinstance(x, cls_)
            x.close()

    def test_hic_based_auto_identification(self, tmpdir):
        hic = kaic.sample_hic()
        for class_name in ('ABDomains', 'ABDomainMatrix', 'ExpectedContacts', 'ObservedExpectedRatio',
                           'FoldChangeMatrix', 'ABDomains', 'PossibleContacts', 'RegionContactAverage',
                           'InsulationIndex', 'DirectionalityIndex'):
            print class_name
            file_name = str(tmpdir) + '/{}.h5'.format(class_name)
            cls_ = class_name_dict[class_name]
            x = cls_(hic, file_name=file_name, mode='w')
            x.close()

            x = kaic.load(file_name, mode='r')
            assert isinstance(x, cls_)
            x.close()
