from kaic.data.general import FileGroup
from abc import abstractmethod, ABCMeta


class ArchitecturalFeature(FileGroup):

    __metaclass__ = ABCMeta

    def __init__(self, group, file_name=None, mode='a', tmpdir=None):
        FileGroup.__init__(self, group, file_name=file_name, mode=mode, tmpdir=tmpdir)

    @classmethod
    @abstractmethod
    def from_file(cls, file_name):
        pass


class HicArchitecture(object):
    def __init__(self, hic):
        self.hic = hic

    def distance_decay(self, per_chromosome=False):
        pass