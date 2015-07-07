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


def _do_map(tmp_input_file, bowtie_index, chromosome, quality_threshold=30):
    bowtie_executable_path = subprocess.Popen("which bowtie2", shell=True, stdout=subprocess.PIPE).stdout.read().rstrip();
    
    tmp_output_file = tempfile.NamedTemporaryFile(dir="./",delete=False)
    tmp_output_file.close()
    
    logging.info("SAM file: %s" % tmp_output_file.name)
    
    bowtieMapCommand = '%s --very-sensitive --no-unal --score-min "C,0,-1" -x %s -q -U %s -S %s' % (bowtie_executable_path,bowtie_index,tmp_input_file,tmp_output_file.name);
    subprocess.call(bowtieMapCommand, shell=True)
    
    mappable = []
    with open(tmp_output_file.name, 'r') as f:
        for line in f:
            logging.info(line)
            if line.startswith("@"):
                continue
            
            fields = line.split("\t")
            
            if fields[1] == '4':
                logging.info("unmapped")
                continue
            
            try:
                if int(fields[4]) < quality_threshold:
                    logging.info("quality")
                    continue
            except ValueError:
                continue
            
            xs = False
            for i in range(11,len(fields)):
                if fields[i].startswith('XS'):
                    logging.info("XS")
                    xs = True
                    break
            if xs:
                logging.info("XS")
                continue
            
            m = re.search('chr_(\w+)_pos_(\d+)', fields[0])
            if m:
                chrm = m.group(1)
                ix = m.group(2)
                if ix == fields[3] and chrm == fields[2]:
                    mappable.append(int(ix))
                else:
                    logging.info("Mismatch: %s-%s, %s-%s" %(chrm, fields[2], ix, fields[3]))
            else:
                raise ValueError("Cannot identify read position")
                
            
    
    #unlink(tmp_input_file)
    #unlink(tmp_output_file.name)
    
    return mappable, chromosome

def unique_mappability(genome, bowtie_index, read_length, offset=1, chunk_size=500000, max_jobs=50, quality_threshold=30):
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
        kwargs = {'quality_threshold': quality_threshold}
        job = Job(_do_map,largs,kwlist=kwargs,queue='all.q')
        jobs.append(job)
            
        reads = []

        
    def submit_and_collect(jobs):
        # do the actual mapping
        job_outputs = process_jobs(jobs,max_processes=2)
        
        # retrieve results
        tmp_mappable = {}
        for chromosome in mappable:
            tmp_mappable[chromosome] = []
        
        logging.info("Collecting")
        
        for (i, result) in enumerate(job_outputs): # @UnusedVariable
            m = result[0]
            chromosome = result[1]
            for ix in m:
                tmp_mappable[chromosome].append(ix)
        
        for chromosome in tmp_mappable:
            logging.info("Length %s: %d" % (chromosome, len(tmp_mappable[chromosome])))
            tmp_mappable[chromosome].sort()
            if len(tmp_mappable[chromosome]) > 0:
                logging.info("min: %d, max: %d" % (tmp_mappable[chromosome][0], tmp_mappable[chromosome][-1]))
            mappable[chromosome] = mappable[chromosome]+tmp_mappable[chromosome]
        
        jobs = []
    
    
    for chromosome in [genome["chrV"]]:
    #for chromosome in genome:
        logging.info("Cutting chromosome %s into reads" % chromosome.name)
        mappable[chromosome.name] = []
        
        reads = []
        l = len(chromosome.sequence)
        for i in range(3545,3550,offset):
        #for i in range(0,l,offset):
            if i >= l-read_length:
                continue
            
            read = chromosome.sequence[i:i+read_length]
            r = "@chr_%s_pos_%d\n" % (chromosome.name, i+1)
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
                
    if len(jobs) > 0:
        submit_and_collect(jobs)
        jobs = []
    
    
    print mappable
    
    mappable_ranges = {} 
    for chrm in mappable:
        mappable_ranges[chrm] = []
        
        if len(mappable[chrm]) > 0:
            current_start = mappable[chrm][0]
            previous = mappable[chrm][0]
            for ix in mappable[chrm]:
                if ix-previous > offset:
                    mappable_ranges[chrm].append([current_start, previous])
                    current_start = ix
                previous = ix
            mappable_ranges[chrm].append([current_start, previous])
            
    print mappable_ranges

    return mappable_ranges
    
        
        