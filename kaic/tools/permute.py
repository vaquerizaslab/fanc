from collections import defaultdict
from kaic.data.genomic import GenomicRegion
from kaic.tools.files import read_chromosome_sizes
from kaic.tools.general import RareUpdateProgressBar
from future.utils import string_types
from pybedtools.helpers import chromsizes
import os
import random


def iter_randomized_regions(original_regions, iterations=1, chromosome_sizes=None, method='unconstrained',
                            attribute='score',
                            preserve_attributes=False, sort=False, silent=True, _chromosome_regions=None):
    if method == 'unconstrained':
        if chromosome_sizes is None:
            raise ValueError("Must provide chromosome_sizes dict when using unconstrained randomization method")
        for i in range(iterations):
            yield _random_regions_unconstrained(original_regions, chromosome_sizes,
                                                preserve_attributes=preserve_attributes)
    elif method == 'spacing':
        if _chromosome_regions is None:
            _chromosome_regions = chromosome_regions(original_regions, sort=sort)
        for i in range(iterations):
            yield _random_regions_spacing(_chromosome_regions, sort=False,
                                          preserve_attributes=preserve_attributes, silent=silent)
    elif method == 'attribute':
        attributes = []
        for region in original_regions:
            a = getattr(region, attribute)
            attributes.append(a)
        for i in range(iterations):
            yield _random_regions_attribute(original_regions, attribute=attribute,
                                            preserve_attributes=preserve_attributes,
                                            _attributes=attributes)
    else:
        raise ValueError("Unknown randomization method '{}'".format(method))


def randomize_regions(original_regions, chromosome_sizes=None, method='unconstrained',
                      preserve_attributes=False, sort=False, silent=True, _chromosome_regions=None):
    random_regions = []

    if method == 'unconstrained':
        if chromosome_sizes is None:
            raise ValueError("Must provide chromosome_sizes dict when using unconstrained randomization method")
        random_regions = _random_regions_unconstrained(original_regions, chromosome_sizes,
                                                       preserve_attributes=preserve_attributes)
    elif method == 'spacing':
        if _chromosome_regions is None:
            _chromosome_regions = chromosome_regions(original_regions, sort=sort)
        random_regions = _random_regions_spacing(_chromosome_regions,
                                                 preserve_attributes=preserve_attributes, silent=silent)

    return random_regions


def chromosome_regions(original_regions, sort=True):
    if not isinstance(original_regions, dict):
        cr = defaultdict(list)
        for region in original_regions:
            cr[region.chromosome].append(region.copy())
    else:
        cr = original_regions

    if sort:
        for chromosome, regions in cr.items():
            regions.sort(key=lambda x: x.start)
    return cr


def _random_regions_unconstrained(original_regions, chromosome_sizes, preserve_attributes=False):
    random_regions = []

    if isinstance(chromosome_sizes, string_types):
        if os.path.isfile(os.path.expanduser(chromosome_sizes)):
            chromosome_sizes = read_chromosome_sizes(chromosome_sizes)
        else:
            genome = chromosome_sizes
            chromosome_sizes = dict()
            for chromosome, (start, end) in chromsizes(genome).items():
                chromosome_sizes[chromosome] = end

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
        cr = original_regions
        if sort:
            for chromosome, regions in cr.items():
                regions.sort(key=lambda x: x.start)
    else:
        cr = chromosome_regions(original_regions, sort=sort)

    for chromosome, regions in cr.items():
        spacing_lens = []
        for i in range(len(regions) - 1):
            l = regions[i + 1].start - regions[i].end
            spacing_lens.append(l)

        current_start = regions[0].start
        shuffled_regions = sorted(regions, key=lambda *args: random.random())
        shuffled_spacings = sorted(spacing_lens, key=lambda *args: random.random())
        protected_attributes = {'chromosome', 'start', 'end'}
        with RareUpdateProgressBar(max_value=len(regions), silent=silent) as pb:
            for i in range(len(shuffled_regions)):
                region_len = len(shuffled_regions[i])
                attributes = {}
                if isinstance(preserve_attributes, bool):
                    if preserve_attributes:
                        for a in shuffled_regions[i].attributes:
                            if a not in protected_attributes:
                                attributes[a] = getattr(shuffled_regions[i], a)
                elif preserve_attributes is not None:
                    for a in shuffled_regions[i].attributes:
                        if a not in protected_attributes and a in preserve_attributes:
                            attributes[a] = getattr(shuffled_regions[i], a)
                random_region = GenomicRegion(start=current_start, end=current_start + region_len,
                                              chromosome=chromosome, **attributes)

                random_regions.append(random_region)
                if i < len(spacing_lens):
                    current_start += region_len + shuffled_spacings[i]
                pb.update(i)
    return random_regions


def _random_regions_attribute(original_regions, attribute='score', preserve_attributes=False,
                              _attributes=None):
    random_regions = []
    if _attributes is None:
        attributes = []
        for region in original_regions:
            a = getattr(region, attribute)
            attributes.append(a)
    else:
        attributes = _attributes

    attributes = sorted(attributes, key=lambda *args: random.random())

    for i, region in enumerate(original_regions):
        if preserve_attributes:
            random_region = region.copy()
        else:
            random_region = GenomicRegion(chromosome=region.chromosome, start=region.start, end=region.end)
        random_region.set_attribute(attribute, attributes[i])

        random_regions.append(random_region)
    return random_regions
