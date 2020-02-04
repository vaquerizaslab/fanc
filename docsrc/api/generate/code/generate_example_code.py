import os
import fanc
from fanc.tools.general import mkdir


# start snippet logging
import logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')
# end snippet logging

output_folder = 'api_output'

# start snippet import map
import fanc.map as map
# end snippet import map

# start snippet mapper
mapper = map.Bowtie2Mapper('hg19_chr18_19/hg19_chr18_19',
                           threads=4, min_quality=30)
# end snippet mapper

# start snippet iterative mapping
sam_folder = mkdir(os.path.join(output_folder, 'sam'))
sam_1_file = map.iterative_mapping('SRR4271982_chr18_19_1.fastq.gzip',
                                   os.path.join(sam_folder, 'SRR4271982_1.bam'),
                                   mapper,
                                   threads=1, min_size=15, step_size=10,
                                   restriction_enzyme='HindIII')
sam_2_file = map.iterative_mapping('SRR4271982_chr18_19_2.fastq.gzip',
                                   os.path.join(sam_folder, 'SRR4271982_2.bam'),
                                   mapper,
                                   threads=1, min_size=15, step_size=10,
                                   restriction_enzyme='HindIII')
# end snippet iterative mapping

# start snippet sort sam
from fanc.tools.files import sort_natural_sam
sorted_sam_1_file = sort_natural_sam(sam_1_file)
sorted_sam_2_file = sort_natural_sam(sam_2_file)
# end snippet sort sam


# start snippet genome
from fanc.regions import genome_regions
genome_file = 'hg19_chr18_19.fa'
fragments = genome_regions(genome_file, restriction_enzyme='HindIII')
# end snippet genome

# start snippet read filters
from fanc.pairs import QualityFilter, UniquenessFilter
quality_filter = QualityFilter(30, mask='MAPQ')
uniqueness_filter = UniquenessFilter(strict=True, mask='unique')
# end snippet read filters

# start snippet import pairs
from fanc.pairs import generate_pairs_split as generate_pairs
pairs_folder = mkdir(os.path.join(output_folder, 'pairs'))
pairs = generate_pairs(sam_1_file, sam_2_file, fragments,
                       read_filters=(quality_filter, uniqueness_filter),
                       output_file=os.path.join(pairs_folder, 'example.pairs'),
                       check_sorted=True, threads=4)
# end snippet import pairs

# start snippet pairs info
chromosomes = pairs.chromosomes()
print(chromosomes)
# ['chr18', 'chr19']

pair = pairs[0]
print(pair)
# chr18: 3187827-(3191583[1])-3192106 -- chr18: 3187827-(3192073[-1])-3192106
type(pair)
# fanc.pairs.FragmentReadPair

print(pair.left)
# chr18: 3187827-(3191583[1])-3192106

print(pair.right)
# chr18: 3187827-(3192073[-1])-3192106
type(pair.right)
# fanc.pairs.FragmentRead
print(pair.right.fragment)
# chr18:3187827-3192106
type(pair.right.fragment)
# genomic_regions.regions.GenomicRegion
print(pair.right.position)
# 3192073
print(pair.right.re_distance())
# 33
print(pair.right.strand)
# -1
# end snippet pairs info

# start snippet pairs iter
# all pairs
for pair in pairs.pairs:
    print(pair)
# pair subset (only pairs with both fragments on chromosome 18)
for pair in pairs.pairs(('chr18', 'chr18')):
    print(pair)
# end snippet pairs iter

# start snippet pairs reset
pairs.reset_filters()
# end snippet pairs reset

# start snippet pairs filter example
from fanc.pairs import SelfLigationFilter
sl_filter = SelfLigationFilter(mask='self-ligation')
# alternative:
# from fanc.general import Mask
# sl_filter = SelfLigationFilter(mask=Mask(name='self-ligation', description="Filter for self-ligated fragments")
pairs.filter(sl_filter)
# end snippet pairs filter example

pairs.reset_filters()

# start snippet pairs filter queue
from fanc.pairs import SelfLigationFilter, ReDistanceFilter
sl_filter = SelfLigationFilter(mask='self-ligation')
rd_filter = ReDistanceFilter(500, mask='re-site-distance')
pairs.filter(sl_filter, queue=True)
pairs.filter(rd_filter, queue=True)
pairs.run_queued_filters()
# end snippet pairs filter queue

# start snippet pairs filter stats
statistics = pairs.filter_statistics()
# end snippet pairs filter stats

# start snippet pairs filter masked
pair = pairs[0]
print(pair)
# chr18: 899140-(899308[-1])-899476 -- chr18: 1509911-(1510021[1])-1510076
# end snippet pairs filter masked

# start snippet pairs filter exclude
for pair in pairs.pairs(excluded_filters=['self-ligation']):
    print(pair)
# end snippet pairs filter exclude

# start snippet hic convert
hic_folder = mkdir(os.path.join(output_folder, 'hic'))
hic_file = os.path.join(hic_folder, 'example.hic')
hic = pairs.to_hic(file_name=hic_file)
# end snippet hic convert

hic.close()
hic = fanc.load(hic_file)

# start snippet hic bin
binned_hic = hic.bin(1000000,
                     file_name=os.path.join(hic_folder, 'binned_1mb.hic'),
                     threads=4)
# end snippet hic bin

# start snippet hic filter
from fanc.hic import LowCoverageFilter
lc_filter = LowCoverageFilter(binned_hic, rel_cutoff=0.2)
binned_hic.filter(lc_filter)
binned_hic.run_queued_filters()
# end snippet hic filter

# start snippet hic balance
from fanc.hic import kr_balancing
kr_balancing(binned_hic, whole_matrix=False,
             restore_coverage=False)
# end snippet hic balance