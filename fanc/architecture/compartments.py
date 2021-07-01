import tables

import genomic_regions as gr
from .helpers import vector_enrichment_profile
from ..matrix import RegionMatrixTable, Edge
import numpy as np
from ..tools.general import RareUpdateProgressBar
from ..config import config
from ..regions import Genome
from future.utils import string_types
from Bio.SeqUtils import GC as calculate_gc_content
import warnings
import logging

logger = logging.getLogger(__name__)


class ABCompartmentMatrix(RegionMatrixTable):
    """
    Class representing O/E correlation matrix used to derive AB compartments.

    You can generate an :class:`~ABCompartmentMatrix` from a Hic object
    using the :func:`~ABCompartmentMatrix.from_hic` class method.

    .. code::

        hic = fanc.load("path/to/file.hic")
        ab = ABCompartmentMatrix.from_hic(hic)

    The :code:`ab` object can then be used to calculate compartmentalisation
    eigenvectors and A/B compartment assignments:

    .. code::

        ev = ab.eigenvector()
        domains = ab.domains()

    For more robust A and B calls, you can use a genome file (FASTA) to
    orient the eigenvector so that regions with higher GC content on average
    get assigned positive EV values:

    .. code::

        ev = ab.eigenvector(genome="path/to/genome.fa")
        domains = ab.domains(genome="path/to/genome.fa")

    Finally, you can calculate an AB compertment enrichment profile using

    .. code::

        profile, ev_cutoffs = ab.enrichment_profile(hic)

        # or with genome to orient the EV
        profile, ev_cutoffs = ab.enrichment_profile(hic, genome="path/to/genome.fa")

    """

    _classid = 'ABCOMPARTMENTMATRIX'

    def __init__(self, file_name=None, mode='r', tmpdir=None):
        RegionMatrixTable.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir,
                                   additional_region_fields={
                                       'name': tables.StringCol(1, pos=0),
                                       'score': tables.Float32Col(pos=1),
                                   })
        if 'per_chromosome' not in self.meta:
            self.meta['per_chromosome'] = True
        if 'oe_per_chromosome' not in self.meta:
            self.meta['oe_per_chromosome'] = True
        if 'has_ev' not in self.meta:
            self.meta['has_ev'] = False
        if 'eigenvector' not in self.meta:
            self.meta['eigenvector'] = 0
        if 'gc' not in self.meta:
            self.meta['gc'] = False

    @classmethod
    def from_hic(cls, hic, file_name=None, tmpdir=None,
                 per_chromosome=True, oe_per_chromosome=None,
                 exclude_chromosomes=None):
        """
        Generate an AB compartment matrix from a Hi-C object.

        :param hic: Hi-C object (FAN-C, Juicer, Cooler)
        :param file_name: Path to output file. If not specified, creates file in memory.
        :param tmpdir: Optional. Work in temporary directory until file is closed.
        :param per_chromosome: If ``True`` (default) calculate compartment profile
                               on a per-chromosome basis (recommended). Otherwise calculates
                               profile on the whole matrix - make sure your normalisation
                               is suitable for this (i.e. whole matrix!)
        :param oe_per_chromosome: Use the expected value vector matching each chromosome.
                                  Do not modify this unless you know what you are doing.+
        :param exclude_chromosomes: Exclude these chromosomes from compartment calculations.
        :return: :class:`~fanc.architecture.compartments.ABCompartmentMatrix` object
        """
        ab_matrix = cls(file_name=file_name, mode='w', tmpdir=tmpdir)
        ab_matrix.add_regions(hic.regions, preserve_attributes=False)
        ab_matrix.meta.per_chromosome = per_chromosome
        ab_matrix.meta.oe_per_chromosome = oe_per_chromosome

        if per_chromosome:
            if oe_per_chromosome is None:
                oe_per_chromosome = True
            chromosomes = hic.chromosomes()
            with RareUpdateProgressBar(max_value=len(chromosomes), silent=config.hide_progressbars,
                                       prefix="AB") as pb:
                for chr_ix, chromosome in enumerate(chromosomes):
                    if exclude_chromosomes is not None and chromosome in exclude_chromosomes:
                        continue
                    
                    m = hic.matrix((chromosome, chromosome), oe=True, oe_per_chromosome=oe_per_chromosome)
                    with np.errstate(divide='ignore', invalid='ignore'):
                        corr_m = np.corrcoef(m)

                    logger.debug("Chromosome {}".format(chromosome))
                    for i, row_region in enumerate(m.row_regions):
                        for j, col_region in enumerate(m.col_regions):
                            if j < i:
                                continue

                            source = row_region.ix
                            sink = col_region.ix

                            try:
                                if np.isnan(corr_m[i, j]):
                                    continue
                                ab_matrix.add_edge(Edge(source=source, sink=sink, weight=corr_m[i, j]))
                            except IndexError:
                                pass
                    pb.update(chr_ix)
        else:
            if oe_per_chromosome is None:
                oe_per_chromosome = False
            m = hic.matrix(oe=True, oe_per_chromosome=oe_per_chromosome)
            corr_m = np.corrcoef(m)
            with RareUpdateProgressBar(max_value=m.shape[0], silent=config.hide_progressbars,
                                       prefix="AB") as pb:
                for i, row_region in enumerate(m.row_regions):
                    for j in range(i, len(m.col_regions)):
                        if j < i:
                            continue
                        col_region = m.row_regions[j]
                        source = row_region.ix
                        sink = col_region.ix
                        ab_matrix.add_edge(Edge(source=source, sink=sink, weight=corr_m[i, j]))
                    pb.update(i)

        ab_matrix.flush()

        mappable = hic.mappable()
        ab_matrix.region_data('valid', mappable)

        return ab_matrix

    def eigenvector(self, sub_region=None, genome=None, eigenvector=0,
                    per_chromosome=None, oe_per_chromosome=None,
                    exclude_chromosomes=None, force=False):
        """
        Calculate the eigenvector (EV) of this AB matrix.

        :param sub_region: Optional region string to only output the EV of
                           that region.
        :param genome: A :class:`~fanc.regions.Genome` object or path to a
                       FASTA file. Used to orient EV value signs so that the "A"
                       compartment corresponds to the regions with higher
                       GC content. It is recommended to make use of this,
                       as otherwise the sign of the EV is arbitrary and
                       will not allow for between-sample comparisons.
        :param eigenvector: Index of the eigenvector to calculate. This
                            parameter is 0-based! Always try "0" first,
                            and if that EV does not seem to reflect A/B
                            compartments, try increasing that value.
        :param per_chromosome: Calculate the eigenvector on a per-chromosome
                               basis (``True`` by default). If your matrix
                               is whole-genome normalised and you know what
                               you are doing, set this to ``False`` to calculate
                               the EV on the whole matrix.
        :param oe_per_chromosome: Use the expected value vector matching each chromosome.
                                  Do not modify this unless you know what you are doing.
        :param exclude_chromosomes: List of chromosome names to exclude from the EV
                                    calculation. Can sometimes be useful if certain
                                    chromosomes do not produce reasonable compartment
                                    profiles.
        :param force: Force EV recalculation, even if the EV has already been previously
                      calculated with the same parameters and is stored in the object.
        :return: :class:`~numpy.array` of eigenvector values
        """
        if per_chromosome is None:
            per_chromosome = self.meta.per_chromosome

        if oe_per_chromosome is None:
            oe_per_chromosome = self.meta.oe_per_chromosome

        if exclude_chromosomes is None:
            exclude_chromosomes = set()
        else:
            exclude_chromosomes = set(exclude_chromosomes)
        
        logger.debug("Excluding chromosomes: {}".format(exclude_chromosomes))

        if (not force and self.meta.has_ev and
                oe_per_chromosome == self.meta.oe_per_chromosome and
                per_chromosome == self.meta.per_chromosome and
                eigenvector == self.meta.eigenvector and
                (genome is not None) == self.meta.gc):
            logger.info("Returning pre-calculated eigenvector!")
            ev = np.array([region.score for region in self.regions])
        else:
            ev = np.zeros(len(self.regions))
            if per_chromosome:
                for chromosome_sub in self.chromosomes():
                    m = self.matrix((chromosome_sub, chromosome_sub))
                    m[np.isnan(m)] = 0
                    w, v = np.linalg.eig(m)
                    # v might end up being masked
                    if hasattr(v, 'mask'):
                        v.mask = False

                    if v.shape[1] <= eigenvector:
                        warning_text = 'Cannot extract EV {} from chromosome {}, ' \
                                       'because it is too small. ' \
                                       'Setting EV values on this chromosome ' \
                                       'to NaN.'.format(eigenvector, chromosome_sub)
                        warnings.warn(warning_text)
                        logger.warning(warning_text)
                        for i, region in enumerate(m.row_regions):
                            ev[region.ix] = np.nan
                    else:
                        ab_vector = v[:, eigenvector]
                        for i, region in enumerate(m.row_regions):
                            ev[region.ix] = ab_vector[i]
            else:
                m = self.matrix()
                m[np.isnan(m)] = 0
                w, v = np.linalg.eig(m)
                # v might end up being masked
                if hasattr(v, 'mask'):
                    v.mask = False
                ab_vector = v[:, eigenvector]
                for i, region in enumerate(m.row_regions):
                    ev[region.ix] = ab_vector[i]

            if genome is not None:
                logger.info("Using GC content to orient eigenvector...")
                close_genome = False
                if isinstance(genome, string_types):
                    genome = Genome.from_string(genome, mode='r')
                    close_genome = True
                genome_chromosomes = genome.chromosomes()

                ignored_chromosomes = []
                gc_content = [np.nan] * len(self.regions)
                for chromosome_sub in self.chromosomes():
                    if chromosome_sub in exclude_chromosomes:
                        continue
                    logger.debug("{}".format(chromosome_sub))

                    if chromosome_sub not in genome_chromosomes:
                        if chromosome_sub.startswith('chr') and chromosome_sub[3:] in genome_chromosomes:
                            logger.warning("Chromosome {} not found in genome, "
                                           "but found {} - assuming these are the same.".format(chromosome_sub, 
                                                                                                chromosome_sub[3:]))
                            chromosome_sub = chromosome_sub[3:]
                        elif 'chr{}'.format(chromosome_sub) in genome_chromosomes:
                            logger.warning("Chromosome {} not found in genome, "
                                           "but found {} - assuming these are the same.".format(chromosome_sub, 
                                                                                                'chr{}'.format(chromosome_sub)))
                            chromosome_sub = 'chr{}'.format(chromosome_sub)
                        else:
                            ignored_chromosomes.append(chromosome_sub)
                            exclude_chromosomes.add(chromosome_sub)
                            continue
                    
                    chromosome_sequence = genome[chromosome_sub].sequence
                    for region in self.regions(chromosome_sub):
                        s = chromosome_sequence[region.start - 1:region.end]
                        gc_content[region.ix] = calculate_gc_content(s)
                
                if len(ignored_chromosomes) > 0:
                    warnings.warn("Several chromosomes could not be found in the provided genome: "
                                "{chromosomes}. GC eigenvector orientation has not been run for "
                                "these chromosomes! Please ensure that this is what you intended!"
                                "".format(chromosomes=", ".join(ignored_chromosomes)))
                    logger.warning("Several chromosomes could not be found in the provided genome: "
                                "{chromosomes}. GC eigenvector orientation has not been run for "
                                "these chromosomes! Please ensure that this is what you intended!"
                                "".format(chromosomes=", ".join(ignored_chromosomes)))
                
                gc_content = np.array(gc_content)

                if close_genome:
                    genome.close()

                # use gc content to orient AB domain vector per chromosome
                cb = self.chromosome_bins
                for chromosome_sub, (start, end) in cb.items():
                    ev_sub = ev[start:end]
                    gc_sub = gc_content[start:end]
                    a_ixs = np.where(ev_sub >= 0.)
                    b_ixs = np.where(ev_sub < 0.)
                    
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=RuntimeWarning)
                        gc_a = np.nanmean(gc_sub[a_ixs])
                        gc_b = np.nanmean(gc_sub[b_ixs])

                    logger.debug("Chromosome {} GC content A: {}, B: {}".format(
                        chromosome_sub, gc_a, gc_b
                    ))

                    if gc_a < gc_b:  # AB compartments are reversed!
                        logger.debug("Reversing EV for chromosome {}".format(chromosome_sub))
                        ev[start:end] = -1 * ev_sub

            try:
                self.region_data('score', ev)
                self.meta.per_chromosome = per_chromosome
                self.meta.oe_per_chromosome = oe_per_chromosome
                self.meta.has_ev = True
                self.meta.eigenvector = eigenvector
                self.meta.gc = genome is not None
                compartments = ['A' if s >= 0 else 'B' for s in ev]
                self.region_data('name', compartments)
            except OSError:
                logger.warning("Cannot save eigenvector to ABCompartmentMatrix since the file "
                               "is opened in read-only mode.")

        if sub_region is None:
            return ev

        regions = list(self.regions(sub_region))
        return ev[regions[0].ix:regions[-1].ix + 1]

    def domains(self, *args, **kwargs):
        """
        Get the AB domain regions of the compartment matrix.

        This returns a :class:`~genomic_regions.RegionWrapper` object, where you can
        iterate over the domains using

        .. code::

            for region in domains.regions:
                print(region.name)  # A or B

        :param args: Positional arguments for
                     :func:`~fanc.architecture.compartments.ABCompartmentMatrix.eigenvector`
        :param kwargs: Keyword arguments for
                      :func:`~fanc.architecture.compartments.ABCompartmentMatrix.eigenvector`
        :return: A :class:`~genomic_regions.RegionWrapper` object
        """
        ev = self.eigenvector(*args, **kwargs)

        region_subset = kwargs.get('region', None)

        current_ev_scores = []
        domains = []
        current_domain = None
        for i, region in enumerate(self.regions(region_subset)):
            if current_domain is None:
                current_domain = gr.GenomicRegion(chromosome=region.chromosome,
                                                  start=region.start, end=region.end,
                                                  name='A' if ev[i] >= 0 else 'B')
                current_ev_scores = [ev[i]]
            else:
                if current_domain.chromosome == region.chromosome and \
                        (ev[i] < 0) == (current_ev_scores[0] < 0):
                    current_domain.end = region.end
                    current_ev_scores.append(ev[i])
                else:
                    current_domain.score = np.nanmean(current_ev_scores)
                    domains.append(current_domain)
                    current_domain = gr.GenomicRegion(chromosome=region.chromosome,
                                                      start=region.start, end=region.end,
                                                      name='A' if ev[i] >= 0 else 'B')
                    current_ev_scores = [ev[i]]
        if current_domain is not None:
            current_domain.score = np.nanmean(current_ev_scores)
            domains.append(current_domain)

        return gr.RegionWrapper(domains)

    def enrichment_profile(self, hic, percentiles=(20.0, 40.0, 60.0, 80.0, 100.0),
                           only_gc=False, symmetric_at=None,
                           exclude_chromosomes=(),
                           intra_chromosomal=True, inter_chromosomal=False,
                           eigenvector=None, collapse_identical_breakpoints=False,
                           *args, **kwargs):
        """
        Generate a compartment enrichment profile for the compartment matrix.

        This returns a :class:`~numpy.ndarray` with the enrichment profile matrix,
        and a list of cutoffs used to bin regions according to the eigenvector (EV)
        values. These cutoffs are determined by the ``percentiles argument``.

        The returned objects can be used to generate a saddle plot, for example using
        :func:`~fanc.plotting.saddle_plot`

        :param hic: A Hi-C matrix
        :param percentiles: The percentiles at which to split the EV, and bin genomic
                            regions accordingly into ranges of EV values.
        :param only_gc: If True, use only the region's GC content, and not the EV,
                        to calculate the enrichment profile.
        :param symmetric_at: If set to a float, splits the genomic regions into two
                             groups with EV below and above this value. Percentiles
                             are then calculated on each group separately, and it is
                             ensured that the ``symmetric_at`` breakpoint is in the
                             centre of the enrichment profile. Note that this doubles
                             the number of bins, and that the number of regions to the
                             left and right of the breakpoint are likely not the same.
        :param exclude_chromosomes: List of chromosome names to exclude from the
                                    profile calculation.
        :param intra_chromosomal: If ``True`` (default), include intra-chromosomal contacts
                                  in the calculation
        :param inter_chromosomal: If ``True``, include inter-chromosomal contacts
                                  in the calculation. This is disabled by defaults, due to
                                  the way matrices are typically normalised (per-chromosome)
        :param eigenvector: Optional. A custom eigenvector of the same length as genomic
                            regions in the Hi-C matrix. This will skip the eigenvector
                            calculation and just use the values in this vector instead.
                            In principle, you could even use this to supply a completely
                            different type of data, such as expression values, for the
                            enrichment analysis.
        :param collapse_identical_breakpoints: (experimental) If ``True``, will merge
                                               all breakpoints with the same values
                                               (such as multiple bins with EV=0) into
                                               one. This can make the saddle plot look
                                               cleaner.
        :param args: Positional arguments for
                     :func:`~fanc.architecture.compartments.ABCompartmentMatrix.eigenvector`
        :param kwargs: Keyword arguments for
                     :func:`~fanc.architecture.compartments.ABCompartmentMatrix.eigenvector`
        :return: a :class:`~numpy.ndarray` with the enrichment profile matrix, a list of cutoffs
        """

        if eigenvector is not None:
            ev = np.array(eigenvector)
            if not len(ev) == len(hic.regions):
                raise ValueError("Eigenvector must have same number of "
                                 "entries ({}) as regions in genome ({})!"
                                 .format(len(ev), len(hic.regions)))
        else:
            logger.info("Generating profile...")
            if only_gc:
                genome = kwargs.get('genome', None)
                if genome is None:
                    raise ValueError("genome cannot be 'None' when using GC content "
                                     "for compartment analysis.")
                # calculate GC content
                if isinstance(genome, string_types):
                    logger.info("Loading genome...")
                    genome = Genome.from_string(genome)

                logger.info("Calculating GC content...")
                gc_content = [np.nan] * len(self.regions)
                for chromosome in self.chromosomes():
                    logger.info("{}".format(chromosome))
                    chromosome_sequence = genome[chromosome].sequence
                    for region in self.regions(chromosome):
                        s = chromosome_sequence[region.start - 1:region.end]
                        gc_content[region.ix] = calculate_gc_content(s)
                gc_content = np.array(gc_content)
                ev = gc_content
            else:
                ev = self.eigenvector(exclude_chromosomes=exclude_chromosomes, *args, **kwargs)

        mappable = self.mappable()
        return vector_enrichment_profile(hic, ev, mappable=mappable,
                                         per_chromosome=kwargs.get('per_chromosome', True),
                                         percentiles=percentiles, symmetric_at=symmetric_at,
                                         exclude_chromosomes=exclude_chromosomes,
                                         intra_chromosomal=intra_chromosomal,
                                         inter_chromosomal=inter_chromosomal,
                                         collapse_identical_breakpoints=collapse_identical_breakpoints)
