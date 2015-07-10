'''
Created on Jul 3, 2015

@author: kkruse1
'''

from kaic.data.genomic import Genome
from os import unlink
import tempfile
import subprocess
import re
from gridmap import Job, process_jobs
import logging

logging.basicConfig(level=logging.INFO)


def _do_map(tmp_input_file, bowtie_index, 
            chromosome, quality_threshold=30, 
            #bowtie_parameters='--very-sensitive --score-min "C,0,-1"'):
            bowtie_parameters='--very-sensitive'):
    bowtie_executable_path = subprocess.Popen("which bowtie2", shell=True, stdout=subprocess.PIPE).stdout.read().rstrip();
    
    tmp_output_file = tempfile.NamedTemporaryFile(delete=False)
    tmp_output_file.close()
    
    #logging.info("SAM file: %s" % tmp_output_file.name)
    
    bowtieMapCommand = '%s --no-unal %s -x %s -q -U %s -S %s' % (bowtie_executable_path,bowtie_parameters,bowtie_index,tmp_input_file,tmp_output_file.name);
    subprocess.call(bowtieMapCommand, shell=True)
    
    mappable = []
    with open(tmp_output_file.name, 'r') as f:
        for line in f:
            #logging.info(line)
            if line.startswith("@"):
                continue
            
            fields = line.split("\t")
            
            if fields[1] == '4':
                #logging.info("unmapped")
                continue
            
            try:
                if int(fields[4]) < quality_threshold:
                    #logging.info("quality")
                    continue
            except ValueError:
                continue
            
            xs = False
            for i in range(11,len(fields)):
                if fields[i].startswith('XS'):
                    #logging.info("XS")
                    xs = True
                    break
            if xs:
                #logging.info("XS")
                continue
            
            m = re.search('chr_(\w+)_pos_(\d+)_reg_(\d+)', fields[0])
            if m:
                chrm = m.group(1)
                ix = m.group(2)
                reg = m.group(3)
                if ix == fields[3] and chrm == fields[2]:
                    mappable.append([int(reg),int(ix)])
                else:
                    logging.info("Mismatch: %s-%s, %s-%s" %(chrm, fields[2], ix, fields[3]))
            else:
                raise ValueError("Cannot identify read position")
                
            
    
    unlink(tmp_input_file)
    unlink(tmp_output_file.name)
    
    return mappable, chromosome


def unique_mappability_at_regions(genome, regions, bowtie_index, 
                                  read_length, offset=1, 
                                  chunk_size=500000, max_jobs=50, 
                                  quality_threshold=30, 
                                  bowtie_parameters='--very-sensitive'):
    logging.info("Maximum number of jobs: %d" % max_jobs)
    
    if type(genome) is str:
        genome =  Genome.from_folder(genome)
    
    jobs = []
    mappable = {}
    
    def prepare(jobs, reads, chromosome):        
        tmp_input_file = tempfile.NamedTemporaryFile(dir="./", delete=False)
        for r in reads:
            tmp_input_file.write(r)
        tmp_input_file.close()
        
        # set up job
        largs = [tmp_input_file.name, bowtie_index, chromosome]
        kwargs = {'quality_threshold': quality_threshold, 'bowtie_parameters': bowtie_parameters}
        job = Job(_do_map,largs,kwlist=kwargs,queue='all.q')
        jobs.append(job)
            
        reads = []

        
    def submit_and_collect(jobs):
        # do the actual mapping
        job_outputs = process_jobs(jobs,max_processes=2)
        
        # retrieve results
        tmp_mappable = {}

                
        for (i, result) in enumerate(job_outputs): # @UnusedVariable
            m = result[0]
            chromosome = result[1]
            
            if not chromosome in tmp_mappable:
                tmp_mappable[chromosome] = {}
            
            for reg, ix in m:
                if not reg in tmp_mappable[chromosome]:
                    tmp_mappable[chromosome][reg] = []
                tmp_mappable[chromosome][reg].append(ix)
        
        for chromosome in tmp_mappable:
            for reg in tmp_mappable[chromosome]:
                logging.info("Length %s: %d" % (chromosome, len(tmp_mappable[chromosome])))
                tmp_mappable[chromosome][reg].sort()
                if len(tmp_mappable[chromosome][reg]) > 0:
                    logging.info("min: %d, max: %d" % (tmp_mappable[chromosome][0], tmp_mappable[chromosome][-1]))
                mappable[chromosome][reg] = mappable[chromosome][reg]+tmp_mappable[chromosome][reg]
        
        jobs = []
    
    
    #for chromosome in [genome["chrV"]]:
    for chromosome in genome:
        logging.info("Processing regions for chromosome %s" % chromosome.name)
        mappable[chromosome.name] = []
        
        reads = []
        region_counter = 0
        for region in regions[chromosome.name]:
            mappable[chromosome.name].append([])
            
            start = max(0,region[0])
            end = min(chromosome.length-read_length, region[1]-read_length)
            for i in range(start, end, offset):

                read = chromosome.sequence[i:i+read_length]
                r = "@chr_%s_pos_%d_reg_%d\n" % (chromosome.name, i+1, region_counter)
                r += read + '\n'
                r += '+\n'
                r += 'G' * len(read) + '\n'
                reads.append(r)
                
                if len(reads) > chunk_size:
                    prepare(jobs, reads, chromosome.name)
                    if len(jobs) == max_jobs:
                        submit_and_collect(jobs)
                        jobs = []
                    reads = []
                
            if len(reads) > 0:
                prepare(jobs, reads, chromosome.name)
                if len(jobs) == max_jobs:
                    submit_and_collect(jobs)
                    jobs = []
            region_counter += 1
                
    if len(jobs) > 0:
        submit_and_collect(jobs)
        jobs = []
    
        
    mappable_ranges = {} 
    for chrm in mappable:
        mappable_ranges[chrm] = []
        for i in range(0,len(mappable[chrm])):
            mappable_ranges[chrm].append([])
            
            if len(mappable[chrm][i]) > 0:
                current_start = mappable[chrm][i][0]
                previous = mappable[chrm][i][0]
                for ix in mappable[chrm][i]:
                    if ix-previous > offset:
                        mappable_ranges[chrm][i].append([current_start, previous])
                        current_start = ix
                    previous = ix
                mappable_ranges[chrm][i].append([current_start, previous])
            
    return mappable_ranges
    
    
    

def unique_mappability(genome, bowtie_index, 
                       read_length, offset=1, 
                       chunk_size=500000, max_jobs=50, 
                       quality_threshold=30, 
                       #bowtie_parameters='--very-sensitive --score-min "C,0,-1"'):
                       bowtie_parameters='--very-sensitive'):
    logging.info("Maximum number of jobs: %d" % max_jobs)
    
    if type(genome) is str:
        genome =  Genome.from_folder(genome)
    
    regions = {}
    
    for chromosome in genome:
        regions[chromosome.name] = [[0,chromosome.length]]
    
    mappable_regions = unique_mappability_at_regions(genome, regions, bowtie_index, read_length, offset, chunk_size, max_jobs, quality_threshold, bowtie_parameters)
    mappable = {}
    for chromosome in mappable_regions:
        mappable[chromosome] = mappable_regions[chromosome][0]
        # optimize memory usage
        del mappable_regions[chromosome]
    
    return mappable
    
def unique_mappability_at_restriction_sites(genome, bowtie_index,
                                            read_length, offset=1, 
                                            chunk_size=500000, max_jobs=50, 
                                            quality_threshold=30, 
                                            #bowtie_parameters='--very-sensitive --score-min "C,0,-1"'):
                                            bowtie_parameters='--very-sensitive'):
    pass
        