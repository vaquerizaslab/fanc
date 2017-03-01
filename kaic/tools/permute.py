from collections import defaultdict
from kaic.data.genomic import GenomicRegion
from kaic.tools.files import read_chromosome_sizes
from kaic.tools.general import RareUpdateProgressBar
from future.utils import string_types
import random


def iter_randomized_regions(original_regions, iterations=1, chromosome_sizes=None, method='unconstrained',
                            preserve_attributes=False, sort=False, silent=True):
    if method == 'unconstrained':
        if chromosome_sizes is None:
            raise ValueError("Must provide chromosome_sizes dict when using unconstrained randomization method")
        for i in range(iterations):
            yield _random_regions_unconstrained(original_regions, chromosome_sizes,
                                                preserve_attributes=preserve_attributes)
    elif method == 'spacing':
        chromosome_regions = _chromosome_regions(original_regions, sort=sort)
        for i in range(iterations):
            yield _random_regions_spacing(chromosome_regions, sort=False,
                                          preserve_attributes=preserve_attributes, silent=silent)
    else:
        raise ValueError("Unknown randomization method '{}'".format(method))


def randomize_regions(original_regions, chromosome_sizes=None, method='unconstrained',
                      preserve_attributes=False, sort=False, silent=True):
    random_regions = []

    if method == 'unconstrained':
        if chromosome_sizes is None:
            raise ValueError("Must provide chromosome_sizes dict when using unconstrained randomization method")
        random_regions = _random_regions_unconstrained(original_regions, chromosome_sizes,
                                                       preserve_attributes=preserve_attributes)
    elif method == 'spacing':
        chromosome_regions = _chromosome_regions(original_regions, sort=sort)
        random_regions = _random_regions_spacing(chromosome_regions,
                                                 preserve_attributes=preserve_attributes, silent=silent)

    return random_regions


def _chromosome_regions(original_regions, sort=True):
    if not isinstance(original_regions, dict):
        chromosome_regions = defaultdict(list)
        for region in original_regions:
            chromosome_regions[region.chromosome].append(region)
    else:
        chromosome_regions = original_regions

    if sort:
        for chromosome, regions in chromosome_regions.items():
            regions.sort(key=lambda x: x.start)
    return chromosome_regions


def _random_regions_unconstrained(original_regions, chromosome_sizes, preserve_attributes=False):
    random_regions = []

    if isinstance(chromosome_sizes, string_types):
        chromosome_sizes = read_chromosome_sizes(chromosome_sizes)

    protected_attributes = {'chromosome', 'start', 'end'}
    for region in original_regions:
        if region.chromosome not in chromosome_sizes:
            continue
        max_d = chromosome_sizes[region.chromosome] - len(region)
        random_start = random.randint(1, max_d)
        random_end = random_start + len(region)
        attributes = {}
        if preserve_attributes:
            for a in region.attributes:
                if a not in protected_attributes:
                    attributes[a] = getattr(region, a)

        random_region = GenomicRegion(chromosome=region.chromosome, start=random_start, end=random_end)

        random_regions.append(random_region)
    return random_regions


def _random_regions_spacing(original_regions, sort=False, preserve_attributes=False, silent=True):
    random_regions = []
    if isinstance(original_regions, dict):
        chromosome_regions = original_regions
        if sort:
            for chromosome, regions in chromosome_regions.items():
                regions.sort(key=lambda x: x.start)
    else:
        chromosome_regions = _chromosome_regions(original_regions, sort=sort)

    for chromosome, regions in chromosome_regions.items():
        spacing_lens = []
        for i in range(len(regions) - 1):
            spacing_lens.append(regions[i + 1].start - regions[i].end)

        current_start = regions[0].start
        random.shuffle(regions)
        random.shuffle(spacing_lens)
        protected_attributes = {'chromosome', 'start', 'end'}
        with RareUpdateProgressBar(max_value=len(regions), silent=silent) as pb:
            for i in range(len(regions)):
                region_len = len(regions[i])
                attributes = {}
                if preserve_attributes:
                    for a in regions[i].attributes:
                        if a not in protected_attributes:
                            attributes[a] = getattr(regions[i], a)
                random_region = GenomicRegion(start=current_start, end=current_start + region_len,
                                              chromosome=chromosome, **attributes)

                random_regions.append(random_region)
                if i < len(spacing_lens):
                    current_start += region_len + spacing_lens[i]
                pb.update(i)
    return random_regions
