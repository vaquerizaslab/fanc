import argparse;
import time;
from hiclib import mapping;
import os;

def splitList(thisList):
    return thisList.split(",");



def iterative_mapping(input_files, output_files,
                      bowtie_index, bowtie=None, 
                      min_length=25, step_size=2, 
                      threads=4, tmp_dir="/tmp", 
                      quality=30, do_join=True,
                      do_filter=True, do_clean=True):
    if not len(input_files) == len(output_files):
        raise ValueError("Input and output arguments must have the same number of elements!");
    if bowtie is None:
        bc = os.popen('which bowtie2')
        bowtie = bc.read().rstrip()
        
    for i in range(0,len(input_files)):
        print("Mapping " + input_files[i]);
        print("Saving to %s " % (output_files[i]));
        time.sleep(3);
        
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
            bowtie_flags=" --very-sensitive --no-unal "
        )
        
        print("Iterative mapping done.");
    
        if do_join:
            print("Joining output files into %s" % (output_files[i]));
            if(do_filter):
                print("Filtering...");
                # merge and filter output sam files for:
                # - mappability
                # - mapping quality
                # - uniqueness
                os.popen('{ head -n 100 -q %s.* | grep "^@[HD|SQ]" | sort | uniq & cat %s.* | awk \'{ duplicate=0; for(i=12; i<=NF; ++i) if($i ~ /^XS/) duplicate=1; if(duplicate == 0 && $1 !~ /^@/ && $2 != 4 && $5 > %d) print $0; }\' | sort -k1,1 -k5,5rn | awk -v last="" \'{ if(last != $1) print $0; last=$1}\'; } > %s' % (output_files[i], output_files[i], quality, output_files[i]));
            else:
                os.popen('{ head -n 100 -q %s.* | grep "^@[HD|SQ]" | sort | uniq & cat %s.* | awk \'!/^@/\' | sort -k1,1 -k5,5rn | awk -v last="" \'{ if(last != $1) print $0; last=$1}\'; } > %s' % (output_files[i],output_files[i],output_files[i]));
            # remove partial sam files
            
            if do_clean:
                os.popen('rm -f %s.*' % (output_files[i]));

        
        
def main(args):
    if not len(args.input) == len(args.output): raise ValueError("Input and output arguments must have the same number of elements!");
    print("Performing iterative mapping on " + ', '.join(args.input));
        
    print("Using the following settings");
    for arg, val in args.__dict__.iteritems():
        print arg, " = ", val;
        
    
    
    time.sleep(5);
    
    for i in range(0,len(args.input)):
        print("Mapping " + args.input[i]);
        print("Saving to %s " % (args.output[i]));
        time.sleep(3);
        
        mapping.iterative_mapping(
            bowtie_path=args.bowtie,
            bowtie_index_path=args.index,
            fastq_path=args.input[i],
            out_sam_path=args.output[i],
            min_seq_len=args.minLength,
            len_step=args.step,
            nthreads=args.threads,
            # max_reads_per_chunk = 10000000,  #optional, on low-memory machines
            temp_dir=args.tmp,
            bowtie_flags=" --very-sensitive --no-unal "
        )
        
        print("Iterative mapping done.");
    
        if(args.join):
            print("Joining output files into %s" % (args.output[i]));
            if(args.filter):
                print("Filtering...");
                # merge and filter output sam files for:
                # - mappability
                # - mapping quality
                # - uniqueness
                os.popen('{ head -n 100 -q %s.* | grep "^@[HD|SQ]" | sort | uniq & cat %s.* | awk \'{ duplicate=0; for(i=12; i<=NF; ++i) if($i ~ /^XS/) duplicate=1; if(duplicate == 0 && $1 !~ /^@/ && $2 != 4 && $5 > %d) print $0; }\' | sort -k1,1 -k5,5rn | awk -v last="" \'{ if(last != $1) print $0; last=$1}\'; } > %s' % (args.output[i], args.output[i], args.quality, args.output[i]));
            else:
                os.popen('{ head -n 100 -q %s.* | grep "^@[HD|SQ]" | sort | uniq & cat %s.* | awk \'!/^@/\' | sort -k1,1 -k5,5rn | awk -v last="" \'{ if(last != $1) print $0; last=$1}\'; } > %s' % (args.output[i],args.output[i],args.output[i]));
            # remove partial sam files
            
            if args.clean:
                os.popen('rm -f %s.*' % (args.output[i]));

if __name__ == '__main__':
    parser = argparse.ArgumentParser();
    bc = os.popen('which bowtie2');
    bowtie_path = bc.read().rstrip();
    print bowtie_path;
    
    parser.add_argument(
        'input',
        type=splitList,
        help='''Input FASTQ files (comma-separated)'''
    );
    
    parser.add_argument(
        'index',
        help='''Bowtie 2 index (include prefix)'''
    );
    
    parser.add_argument(
        'output',
        type=splitList,
        help='''Output files (must be same number as input files)'''
    );
    
    parser.add_argument(
        '-b', '--bowtie', dest='bowtie',
        default=bowtie_path,
        help='''Bowtie 2 executable path (will check PATH variable by default)'''
    );
    
    parser.add_argument(
        '-m', '--min-length', dest='minLength',
        type=int,
        default=25,
        help='''Minimum sequence length to attempt the mapping'''
    );
    
    parser.add_argument(
        '-s', '--step', dest='step',
        type=int,
        default=2,
        help='''Step size to increase mapped sequence length'''
    );
    
    parser.add_argument(
        '-t', '--threads', dest='threads',
        type=int,
        default=8,
        help='''Number of threads'''
    );
    
    parser.add_argument(
        '-tmp', '--temp-dir', dest='tmp',
        default="/tmp",
        help='''Temporary directory'''
    );
    
    parser.add_argument(
        '-nj', '--no-join', dest='join',
        action='store_false',
        help='''Do not join partial SAM files into a single file'''
    );
    
    parser.add_argument(
        '-nf', '--no-filter', dest='filter',
        action='store_false',
        help='''Do not filter output SAM file by mappability, mapping quality, and uniqueness'''
    );
    
    parser.add_argument(
        '-nc', '--no-clean', dest='clean',
        action='store_false',
        help='''Do not delete partial SAM files (*.sam.\d)'''
    );
    
    parser.add_argument(
        '-q', '--quality', dest='quality',
        type=int,
        default=30,
        help='''Mapping quality to filter mapped reads'''
    );
    
    parser.set_defaults(join=True);
    parser.set_defaults(filter=True);
    parser.set_defaults(clean=True);
    
    
    main(parser.parse_args());
