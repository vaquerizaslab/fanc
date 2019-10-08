import os
import logging
import kaic
import kaic.map as map
from kaic.tools.general import mkdir
from kaic.tools.files import sort_natural_sam

# start snippet logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')
# end snippet logging

output_folder = 'api_output'


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
sorted_sam_1_file = sort_natural_sam(sam_1_file)
sorted_sam_2_file = sort_natural_sam(sam_2_file)
# end snippet sort sam