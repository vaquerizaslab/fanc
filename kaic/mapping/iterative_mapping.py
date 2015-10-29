from hiclib import mapping
import os
import logging
import glob


def iterative_mapping(input_files, output_files,
                      bowtie_index, bowtie=None, 
                      min_length=25, step_size=2, 
                      threads=4, tmp_dir="/tmp", 
                      quality=30, do_join=True,
                      do_filter=False, do_clean=True,
                      bowtie_options='--very-sensitive --no-unal'):
    if not len(input_files) == len(output_files):
        raise ValueError("Input and output arguments must have the same number of elements!")
    if bowtie is None:
        bc = os.popen('which bowtie2')
        bowtie = bc.read().rstrip()
        
    for i in range(0, len(input_files)):
        logging.info("Mapping %s" % input_files[i])
        logging.info("Saving to %s " % (output_files[i]))

        mapping.iterative_mapping(
            bowtie_path=bowtie,
            bowtie_index_path=bowtie_index,
            fastq_path=input_files[i],
            out_sam_path=output_files[i],
            min_seq_len=min_length,
            len_step=step_size,
            nthreads=threads,
            # max_reads_per_chunk = 10000000,  #optional, on low-memory machines
            temp_dir=tmp_dir,
            bowtie_flags=" %s " % bowtie_options
        )
        
        logging.info("Iterative mapping done.")
    
        if do_join:
            logging.info("Joining output files into %s" % (output_files[i]))

            sam_files = glob.glob("%s.*" % output_files[i])


            if do_filter:
                logging.info("Filtering...")
                # merge and filter output sam files for:
                # - mappability
                # - mapping quality
                # - uniqueness
                os.popen('{ head -n 20000 -q %s.* | grep "^@[HD|SQ]" | sort | uniq & cat %s.* | awk \'{ duplicate=0; for(i=12; i<=NF; ++i) if($i ~ /^XS/) duplicate=1; if(duplicate == 0 && $1 !~ /^@/ && $2 != 4 && $5 > %d) print $0; }\' | sort -k1,1 -k5,5rn | awk -v last="" \'{ if(last != $1) print $0; last=$1}\'; } > %s' % (output_files[i], output_files[i], quality, output_files[i]));
            else:
                os.popen('{ head -n 20000 -q %s.* | grep "^@[HD|SQ]" | sort | uniq & cat %s.* | awk \'!/^@/\' | sort -k1,1 -k5,5rn | awk -v last="" \'{ if(last != $1) print $0; last=$1}\'; } > %s' % (output_files[i],output_files[i],output_files[i]));
            # remove partial sam files
            
            if do_clean:
                os.popen('rm -f %s.*' % (output_files[i]));

