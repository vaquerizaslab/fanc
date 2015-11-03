from hiclib import mapping
import logging
import glob
import shutil
import tempfile
import subprocess
import os


def iterative_mapping(input_files, output_files,
                      bowtie_index, bowtie=None, 
                      min_length=25, step_size=2, 
                      threads=4, tmp_dir="/tmp", 
                      do_join=True, do_clean=True,
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
            output_path = output_files[i]
            first_sam_file = sam_files.pop(0)

            shutil.copy(first_sam_file, output_path)
            if do_clean:
                os.unlink(first_sam_file)

            header_command = 'head -n 20000 -q %s %s | grep "^@[HD|SQ]" | sort | uniq > %s'
            merge_command = 'cat %s %s | awk \'!/^@/\' | sort -k1,1 -k5,5rn | ' + \
                            'awk -v last="" \'{ if(last != $1) print $0; last=$1}\' >> %s'

            # join two files and copy back to output_path
            for sam_file in sam_files:
                tmp_file = tempfile.NamedTemporaryFile(delete=False)
                tmp_file.close()
                tmp_path = tmp_file.name
                # tmp_path = "%s.tmp" % output_path

                logging.info("Merging headers")
                try:
                    logging.info(header_command % (output_path, sam_file, tmp_path))
                    subprocess.check_call(header_command % (output_path, sam_file, tmp_path), shell=True)
                except subprocess.CalledProcessError:
                    logging.error("Could not join headers of %s into file!" % sam_file)
                    os.unlink(tmp_path)
                    continue

                logging.info("Merging reads")
                try:
                    logging.info(merge_command % (output_path, sam_file, tmp_path))
                    subprocess.check_call(merge_command % (output_path, sam_file, tmp_path), shell=True)
                except subprocess.CalledProcessError:
                    logging.error("Could not join %s into file!" % sam_file)
                    os.unlink(tmp_path)
                    continue

                if do_clean:
                    os.unlink(sam_file)
                os.unlink(output_path)
                shutil.move(tmp_path, output_path)
