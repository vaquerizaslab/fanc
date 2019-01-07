from .matrix import RegionMatrixTable
from abc import abstractmethod, ABCMeta
from future.utils import with_metaclass, string_types
from .general import MaskFilter
import numpy as np
import logging

logger = logging.getLogger(__name__)


class Hic(RegionMatrixTable):
    def __init__(self, file_name=None, mode='a', tmpdir=None,
                 partitioning_strategy='chromosome',
                 _table_name_regions='regions', _table_name_edges='edges',
                 _edge_buffer_size=1000000):
        RegionMatrixTable.__init__(self, file_name=file_name,
                                   mode=mode, tmpdir=tmpdir,
                                   partitioning_strategy=partitioning_strategy,
                                   _table_name_regions=_table_name_regions,
                                   _table_name_edges=_table_name_edges,
                                   _edge_buffer_size=_edge_buffer_size)


class HicEdgeFilter(with_metaclass(ABCMeta, MaskFilter)):
    """
    Abstract class that provides filtering functionality for the
    edges/contacts in a :class:`~Hic` object.

    Extends MaskFilter and overrides valid(self, row) to make
    :class:`~HicEdge` filtering more "natural".

    To create custom filters for the :class:`~Hic` object, extend this
    class and override the valid_edge(self, edge) method.
    valid_edge should return False for a specific :class:`~HicEdge` object
    if the object is supposed to be filtered/masked and True
    otherwise. See :class:`~DiagonalFilter` for an example.

    Pass a custom filter to the :func:`~Hic.filter` method in :class:`~Hic`
    to apply it.
    """

    def __init__(self, mask=None):
        """
        Initialize HicEdgeFilter.

        :param mask: The Mask object that should be used to mask
                     filtered :class:`~HicEdge` objects. If None the default
                     Mask will be used.
        """
        super(HicEdgeFilter, self).__init__(mask)
        self._hic = None

    @abstractmethod
    def valid_edge(self, edge):
        """
        Determine if a :class:`~HicEdge` object is valid or should
        be filtered.

        When implementing custom HicEdgeFilter this method must be
        overridden. It should return False for :class:`~HicEdge` objects that
        are to be fitered and True otherwise.

        Internally, the :class:`~Hic` object will iterate over all HicEdge
        instances to determine their validity on an individual
        basis.

        :param edge: A :class:`~HicEdge` object
        :return: True if :class:`~HicEdge` is valid, False otherwise
        """
        pass

    def set_hic_object(self, hic_object):
        """
        Set the :class:`~Hic` instance to be filtered by this
        HicEdgeFilter.

        Used internally by :class:`~Hic` instance.

        :param hic_object: :class:`~Hic` object
        """
        self._hic = hic_object

    def valid(self, row):
        """
        Map valid_edge to MaskFilter.valid(self, row).

        :param row: A pytables Table row.
        :return: The boolean value returned by valid_edge.
        """
        edge = self._hic._row_to_edge(row, lazy=True)
        return self.valid_edge(edge)


class DiagonalFilter(HicEdgeFilter):
    """
    Filter contacts in the diagonal of a :class:`~Hic` matrix.
    """
    def __init__(self, distance=0, mask=None):
        """
        Initialize filter with chosen parameters.

        :param distance: Distance from the diagonal up to which
                         contacts will be filtered
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered edges.
        """
        HicEdgeFilter.__init__(self, mask=mask)
        self.distance = distance

    def valid_edge(self, edge):
        """
        Check if an edge is on (or near) the diagonal of the :class:`~Hic` matrix.
        """
        if abs(edge.source-edge.sink) <= self.distance:
            return False
        return True


class LowCoverageFilter(HicEdgeFilter):
    """
    Filter a :class:`~HicEdge` if it connects a region that
    does not have a contact count larger than a specified
    cutoff.

    If the cutoff is not provided, it is automatically
    chosen at 10% of the mean contact count of all regions.
    """
    def __init__(self, hic_object, cutoff=None, rel_cutoff=None, mask=None):
        """
        Initialize filter with these settings.

        The cutoff can be provided in two ways:
        1. As an absolute threshold. Regions with contact count below this
        absolute threshold are filtered
        2. As a fraction relative to the median contact count of all regions.

        If both is supplied, whichever threshold is lower will be selected.

        :param hic_object: The :class:`~Hic` object that this
                           filter will be called on. Needed for
                           contact count calculation.
        :param rel_cutoff: A cutoff as a fraction (0-1) of the median contact count of all
                           regions. If cutoff and rel_cutoff are None, will be set to 10%
        :param cutoff: A cutoff in absolute contact counts (can be float) below
                       which regions are considered "low coverage"
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered edges.
        """
        HicEdgeFilter.__init__(self, mask=mask)

        self._marginals = hic_object.marginals()
        if cutoff is None and rel_cutoff is None:
            rel_cutoff = 0.1
            logger.info("Using default 10 percent relative coverage as cutoff")

        if cutoff is not None and rel_cutoff is not None:
            cutoff = min(cutoff if cutoff else float("inf"),
                         self.calculate_cutoffs(rel_cutoff)[0] if rel_cutoff else float("inf"))
        elif rel_cutoff is not None:
            cutoff = self.calculate_cutoffs(rel_cutoff)[0]
        logger.info("Final absolute cutoff threshold is {:.4}".format(float(cutoff)))

        self._regions_to_mask = set()
        for i, contacts in enumerate(self._marginals):
            if contacts < cutoff:
                self._regions_to_mask.add(i)
        logger.info("Selected a total of {} ({:.1%}) regions to be masked".format(
            len(self._regions_to_mask), len(self._regions_to_mask)/len(hic_object.regions)))

    def calculate_cutoffs(self, fraction_threshold=0.05):
        lower = np.median(self._marginals[self._marginals > 0])*fraction_threshold
        upper = np.median(self._marginals[self._marginals > 0])+lower
        return lower, upper

    def valid_edge(self, edge):
        """
        Check if an edge falls into a low-coverage region.
        """
        if edge.source in self._regions_to_mask:
            return False
        if edge.sink in self._regions_to_mask:
            return False
        return True
