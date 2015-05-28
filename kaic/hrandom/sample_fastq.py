'''
Created on May 28, 2015

@author: kkruse1
'''

import random
from kaic.tools.files import get_number_of_lines

def sample_fastq(fastq_file1, fastq_file2, out_fastq1, out_fastq2, sample_size=300000000):
    # get file size
    n_lines = get_number_of_lines(fastq_file1)
    n_entries = n_lines/4
    # make sure we are not sampling more than we already have
    sample_size = min(sample_size, n_entries)
    # get sample indexes
    sample_ixs = sorted(random.sample(range(0, n_entries), sample_size))
    
    
    with open(fastq_file1, 'r') as f1:
        with open(fastq_file2, 'r') as f2:
            with open(out_fastq1, 'w') as o1:
                with open(out_fastq2, 'w') as o2:
                    line_counter = 0
                    sample_counter = 0
                    write_line = False
                    
                    f1_line = f1.readline()
                    f2_line =  f2.readline()
                    
                    while f1_line != '':
                        
                        if line_counter % 4 == 0:
                            write_line = False
                            if (sample_counter < len(sample_ixs) and 
                                line_counter/4 == sample_ixs[sample_counter]):
                                write_line = True
                                sample_counter += 1
                        
                        if write_line:
                            o1.write(f1_line)
                            o2.write(f2_line)
                            
                        
                        # update counters and lines 
                        f1_line = f1.readline()
                        f2_line =  f2.readline()
                        line_counter += 1
