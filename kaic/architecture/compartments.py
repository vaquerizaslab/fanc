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
import logging

logger = logging.getLogger(__name__)


class ABCompartmentMatrix(RegionMatrixTable):

    _classid = 'ABCOMPARTMENTMATRIX'

    def __init__(self, file_name=None, tmpdir=None):
        RegionMatrixTable.__init__(self, file_name=file_name, tmpdir=tmpdir,
                                   additional_region_fields={
                                       'name': tables.StringCol(1, pos=0),
                                       'score': tables.Float32Col(pos=1),
                                   })
        self.meta['per_chromosome'] = True
        self.meta['oe_per_chromosome'] = True
        self.meta['has_ev'] = False
        self.meta['eigenvector'] = 0

    @classmethod
    def from_hic(cls, hic, file_name=None, tmpdir=None,
                 per_chromosome=True, oe_per_chromosome=None):
        ab_matrix = cls(file_name=file_name, tmpdir=tmpdir)
        ab_matrix.add_regions(hic.regions, preserve_attributes=False)
        ab_matrix.meta.per_chromosome = per_chromosome
        ab_matrix.meta.oe_per_chromosome = oe_per_chromosome

        if per_chromosome:
            if oe_per_chromosome is None:
                oe_per_chromosome = True
            chromosomes = hic.chromosomes()
            for chromosome in chromosomes:
                m = hic.matrix((chromosome, chromosome), oe=True, oe_per_chromosome=oe_per_chromosome)
                corr_m = np.corrcoef(m)
                print(corr_m)

                logger.info("Chromosome {}".format(chromosome))
                with RareUpdateProgressBar(max_value=m.shape[0], silent=config.hide_progressbars) as pb:
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
                        pb.update(i)
        else:
            if oe_per_chromosome is None:
                oe_per_chromosome = False
            m = hic.matrix(oe=True, oe_per_chromosome=oe_per_chromosome)
            corr_m = np.corrcoef(m)
            with RareUpdateProgressBar(max_value=m.shape[0], silent=config.hide_progressbars) as pb:
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

        ab_matrix.region_data('valid', list(hic.region_data('valid')))

        return ab_matrix

    def eigenvector(self, region=None, genome=None, eigenvector=0,
                    per_chromosome=None, oe_per_chromosome=None,
                    force=False):
        if per_chromosome is None:
            per_chromosome = self.meta.per_chromosome

        if oe_per_chromosome is None:
            oe_per_chromosome = self.meta.oe_per_chromosome

        if (not force and self.meta.has_ev and
                oe_per_chromosome == self.meta.oe_per_chromosome and
                per_chromosome == self.meta.per_chromosome and
                eigenvector == self.meta.eigenvector):
            ev = np.array([region.score for region in self.regions])
        else:
            ev = np.zeros(len(self.regions))
            if per_chromosome:
                for chromosome_sub in self.chromosomes():
                    m = self.matrix((chromosome_sub, chromosome_sub))
                    m[np.isnan(m)] = 0
                    w, v = np.linalg.eig(m)
                    ab_vector = v[:, eigenvector]
                    for i, region in enumerate(m.row_regions):
                        ev[region.ix] = ab_vector[i]
            else:
                m = self.matrix()
                m[np.isnan(m)] = 0
                w, v = np.linalg.eig(m)
                ab_vector = v[:, eigenvector]
                for i, region in enumerate(m.row_regions):
                    ev[region.ix] = ab_vector[i]

            if genome is not None:
                logger.info("Using GC content to orient eigenvector...")
                if isinstance(genome, string_types):
                    genome = Genome.from_string(self.genome, mode='r')
                else:
                    genome = genome

                gc_content = [np.nan] * len(self.regions)
                for chromosome_sub in self.chromosomes():
                    logger.info("{}".format(chromosome_sub))
                    chromosome_sequence = genome[chromosome_sub].sequence
                    for region in self.regions(chromosome_sub):
                        s = chromosome_sequence[region.start - 1:region.end]
                        gc_content[region.ix] = calculate_gc_content(s)
                gc_content = np.array(gc_content)

                # use gc content to orient AB domain vector per chromosome
                cb = self.chromosome_bins
                for chromosome_sub, (start, end) in cb.items():
                    ev_sub = ev[start:end]
                    gc_sub = gc_content[start:end]
                    a_ixs = np.where(ev_sub >= 0.)
                    b_ixs = np.where(ev_sub < 0.)
                    gc_a = np.nanmean(gc_sub[a_ixs])
                    gc_b = np.nanmean(gc_sub[b_ixs])

                    if gc_a < gc_b:  # AB compartments are reversed!
                        ev[start:end] = -1 * ev_sub

            try:
                self.region_data('score', ev)
                self.meta.per_chromosome = per_chromosome
                self.meta.oe_per_chromosome = oe_per_chromosome
                self.meta.has_ev = True
                self.meta.eigenvector = eigenvector
                compartments = ['A' if s >= 0 else 'B' for s in ev]
                self.region_data('name', compartments)
            except OSError:
                logger.warning("Cannot save eigenvector to ABCompartmentMatrix since the file "
                               "is opened in read-only mode.")

        if region is None:
            return ev

        regions = list(self.regions(region))
        return ev[regions[0].ix:regions[-1].ix + 1]

    def domains(self, *args, **kwargs):
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
                current_ev_scores.append(ev[i])
            else:
                if current_domain.chromosome == region.chromosome and \
                        ev[i] < 0 == current_ev_scores[0] < 0:
                    current_domain.end = region.end
                    current_ev_scores.append(ev[i])
                else:
                    current_domain.score = np.nanmean(current_ev_scores)
                    domains.append(current_domain)
                    current_domain = gr.GenomicRegion(chromosome=region.chromosome,
                                                      start=region.start, end=region.end,
                                                      name='A' if ev[i] >= 0 else 'B')
                    current_ev_scores.append(ev[i])
        if current_domain is not None:
            current_domain.score = np.nanmean(current_ev_scores)
            domains.append(current_domain)

        return gr.RegionWrapper(domains)

    def enrichment_profile(self, hic, percentiles=(20.0, 40.0, 60.0, 80.0, 100.0),
                           only_gc=False, symmetric_at=None,
                           exclude_chromosomes=(), *args, **kwargs):

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
            ev = self.eigenvector(*args, **kwargs)

        mappable = self.mappable()
        return vector_enrichment_profile(hic, ev, mappable=mappable,
                                         per_chromosome=kwargs.get('per_chromosome', True),
                                         percentiles=percentiles, symmetric_at=symmetric_at,
                                         exclude_chromosomes=exclude_chromosomes)
