from collections import defaultdict
from kaic.data.genomic import GenomicRegion
import random


def iter_randiomized_regions(original_regions, iterations=1, chromosome_sizes=None, method='unconstrained'):
    if method == 'unconstrained':
        if chromosome_sizes is None:
            raise ValueError("Must provide chromosome_sizes list when using unconstrained randomization method")
        for i in xrange(iterations):
            yield _random_regions_unconstrained(original_regions, chromosome_sizes)
    elif method == 'spacing':
        chromosome_regions = _chromosome_regions(original_regions)
        for i in xrange(iterations):
            yield _random_regions_spacing(chromosome_regions, sort=False)
    else:
        raise ValueError("Unknown randomization method '{}'".format(method))


def randomize_regions(original_regions, chromosome_sizes=None, method='unconstrained'):
    random_regions = []

    if method == 'unconstrained':
        if chromosome_sizes is None:
            raise ValueError("Must provide chromosome_sizes list when using unconstrained randomization method")
        random_regions = _random_regions_unconstrained(original_regions, chromosome_sizes)
    elif method == 'spacing':
        chromosome_regions = _chromosome_regions(original_regions)
        random_regions = _random_regions_spacing(chromosome_regions, sort=False)

    return random_regions


def _chromosome_regions(original_regions, sort=True):
    if not isinstance(original_regions, dict):
        chromosome_regions = defaultdict(list)
        for region in original_regions:
            chromosome_regions[region.chromosome].append(region)
    else:
        chromosome_regions = original_regions

    if sort:
        for chromosome, regions in chromosome_regions.iteritems():
            regions.sort(key=lambda x: x.start)

    return chromosome_regions


def _random_regions_unconstrained(original_regions, chromosome_sizes):
    random_regions = []

    for region in original_regions:
        if region.chromosome not in chromosome_sizes:
            continue
        max_d = chromosome_sizes[region.chromosome] - len(region)
        random_start = random.randint(1, max_d)
        random_end = random_start + len(region)

        random_region = GenomicRegion(chromosome=region.chromosome, start=random_start, end=random_end)
        random_regions.append(random_region)
    return random_regions


def _random_regions_spacing(original_regions, sort=True):
    random_regions = []
    if isinstance(original_regions, dict):
        chromosome_regions = original_regions
        if sort:
            for chromosome, regions in chromosome_regions.iteritems():
                regions.sort(key=lambda x: x.start)
    else:
        chromosome_regions = _chromosome_regions(original_regions, sort=sort)

    for chromosome, regions in chromosome_regions.iteritems():
        region_lens = []
        spacing_lens = []
        for i in xrange(len(regions) - 1):
            region_lens.append(len(regions[i]))
            spacing_lens.append(regions[i + 1].start - regions[i].end)

        random.shuffle(region_lens)
        random.shuffle(spacing_lens)
        current_start = regions[0].start
        for i in xrange(len(region_lens)):
            random_region = GenomicRegion(start=current_start, end=current_start + region_lens[i],
                                          chromosome=chromosome)
            random_regions.append(random_region)
            current_start += region_lens[i] + spacing_lens[i]
    return random_regions
