import argparse
import logging
import os
import os.path
import textwrap
import shutil
import tempfile
import subprocess
import kaic.commands.auto


# configure logging
logger = logging.getLogger(__name__)


def kaic_parser():
    usage = '''\
        kaic <command> [options]

        Commands:
            auto                Automatically process an entire Hi-C data set
            dirs                Create default folder structure for kaic
            stats               Get statistics for kaic pipeline files

            --- Mapping
            iterative_mapping   Iteratively map a FASTQ file to a Bowtie 2 index

            --- Reads
            load_reads          Load a SAM/BAM file into a Reads object
            filter_reads        Filter a Reads object

            -- Genome
            build_genome        Convenience command to build a Genome object

            --- Pairs
            reads_to_pairs      Convert a Reads object into a Pairs object
            filter_pairs        Filter a Pairs object
            pairs_from_hicpro   Load pairs from a HiC-Pro valid pairs file
            pairs_to_homer      Write pairs in Homer compatible HiCSummary format

            --- Hic
            pairs_to_hic        Convert a pairs object into a Hic object
            filter_hic          Filter a Hic object
            merge_hic           Merge multiple Hic objects
            bin_hic             Bin a Hic object into same-size regions
            correct_hic         Correct a Hic object for biases
            hic_pca             Do a PCA on multiple Hi-C objects
            optimise            Optimise an existing Hic object for faster access
            subset_hic          Create a new Hic object by subsetting
            cis_trans           Calculate cis/trans ratio
            dump                Dump Hic file to txt file(s)
            hic_from_juicer     Convert juicer .hic file to Kai-C Hic
            hic_to_cooler       Convert Hic file into cooler format

            --- Network
            call_peaks          Call enriched peaks in a Hic object
            filter_peaks        Filter peaks called with 'call_peaks'
            merge_peaks         Merge peaks
            filter_merged_peaks Filter merged peaks
            overlap_peaks       Overlap peaks from multiple samples

            --- Plotting
            plot_ligation_err   Plot the ligation error of a Pairs object
            plot_re_dist        Plot the distance of reads to the nearest RE site
            plot_hic_corr       Plot the correlation of two Hic objects
            plot_hic_marginals  Plot marginals in a Hic object
            plot_diff           Plot the difference between two Hic matrices

            --- Architecture
            structure_tracks   Calculate structural features of a Hic object
            boundaries         Call boundaries in an Hic object
            fold_change        Create pairwise fold-change Hi-C comparison maps
            average_tracks     Calculate average Hi-C contact profiles per region
            directionality     Calculate directionality index for Hic object
            insulation         Calculate insulation index for Hic object
            ii_to_bw           Convert insulation index object to BigWig 
            ab                 Calculate AB domain matrix for a Hi-C object
            ab_domains         Assign A or B compartment to each region in matrix
            distance_decay     Calculate distance decay for a Hi-C object
            diff               Calculate difference between two vectors
            aggregate_tads     Make a TAD aggregate plot
            aggregate_loops    Make a loop aggregate plot
            ab_profile         Plot A-B compartment enrichment profiles

            --- Other
            write_config       Write default config file

        Run kaic <command> -h for help on a specific command.
        '''
    parser = argparse.ArgumentParser(
        description="kaic processing tool for Hi-C data",
        usage=textwrap.dedent(usage)
    )

    parser.add_argument(
        '--version', dest='print_version',
        action='store_true',
        help='''Print version information'''
    )
    parser.set_defaults(print_version=False)

    parser.add_argument(
        '--verbose', '-v', dest='verbosity',
        action='count',
        default=0,
        help='''Set verbosity level: Can be chained like '-vvv' to increase verbosity. Default is to show
                        errors, warnings, and info messages (same as '-vv'). '-v' shows only errors and warnings,
                        '-vvv' shows errors, warnings, info, and debug messages in addition.'''
    )

    parser.add_argument(
        '-s', '--silent', dest='silent',
        action='store_true',
        help='''if set, do not print log messages to to command line.'''
    )
    parser.set_defaults(silent=False)

    parser.add_argument(
        '-l', '--log-file', dest='log_file',
        help='''Path to file in which to save log.'''
    )

    parser.add_argument(
        '-m', '--email', dest='email_to_address',
        help='''Email address for kaic command summary.'''
    )

    parser.add_argument(
        '--smtp-server', dest='smtp_server',
        help='''SMTP server in the form smtp.server.com[:port].'''
    )

    parser.add_argument(
        '--smtp-username', dest='smtp_username',
        help='''SMTP username.'''
    )

    parser.add_argument(
        '--smtp-password', dest='smtp_password',
        help='''SMTP password.'''
    )

    parser.add_argument(
        '--smtp-sender-address', dest='email_from_address',
        help='''SMTP sender email address.'''
    )

    parser.add_argument('command', nargs='?', help='Subcommand to run')

    return parser


def dirs_parser():
    parser = argparse.ArgumentParser(
        prog="kaic dirs",
        description='Automatically process an entire Hi-C data set'
    )

    parser.add_argument(
        'root_directory',
        help='''Root directory in which to create kaic folders'''
    )
    return parser


def dirs(argv):
    parser = dirs_parser()

    args = parser.parse_args(argv[2:])

    root_dir = os.path.expanduser(args.root_directory)

    import errno

    try:
        os.makedirs(root_dir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    try:
        os.makedirs(root_dir + '/fastq')
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    try:
        os.makedirs(root_dir + '/sam')
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    try:
        os.makedirs(root_dir + '/reads/filtered')
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    try:
        os.makedirs(root_dir + '/pairs/filtered')
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    try:
        os.makedirs(root_dir + '/hic/filtered')
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    try:
        os.makedirs(root_dir + '/hic/binned')
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    try:
        os.makedirs(root_dir + '/hic/corrected')
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    try:
        os.makedirs(root_dir + '/plots/stats')
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    try:
        os.makedirs(root_dir + '/plots/matrix')
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def auto_parser():
    return kaic.commands.auto.auto_parser()


def auto(argv):
    return kaic.commands.auto.auto(argv)


def map_parser():
    parser = argparse.ArgumentParser(
        prog="kaic map",
        description='Map reads in a FASTQ file to a reference genome'
    )

    parser.add_argument(
        'input',
        nargs='+',
        help='''File name of the input FASTQ file (or gzipped FASTQ)'''
    )

    parser.add_argument(
        'index',
        help='''Bowtie 2 genome index'''
    )

    parser.add_argument(
        'output',
        help='''Output file name (or folder name if multiple input files provided)'''
    )

    parser.add_argument(
        '-m', '--min-size', dest='min_size',
        type=int,
        default=30,
        help='''Minimum length of read before extension. Default 25.'''
    )

    parser.add_argument(
        '-s', '--step-size', dest='step_size',
        type=int,
        default=10,
        help='''Number of base pairs to extend at each round of mapping. Default is 10.'''
    )

    parser.add_argument(
        '--trim-front', dest='trim_front',
        action='store_true',
        help='''Trim reads from front instead of back.'''
    )
    parser.set_defaults(trim_front=False)

    parser.add_argument(
        '-t', '--threads', dest='threads',
        type=int,
        default=1,
        help='''Number of threads used for mapping'''
    )

    parser.add_argument(
        '-q', '--quality', dest='quality',
        type=int,
        default=30,
        help='''Mapping quality cutoff for reads to be sent to another iteration'''
    )

    parser.add_argument(
        '-r', '--restriction-enzyme', dest='restriction_enzyme',
        help='''Name (case sensitive) of restriction enzyme used in Hi-C experiment.
                Will be used to split reads by predicted ligation junction before mapping.'''
    )

    parser.add_argument(
        '-k', '--max-alignments', dest='max_alignments',
        type=int,
        help='''Maximum number of alignments per read to be reported. Default: 1'''
    )

    parser.add_argument(
        '-a', '--all-alignments', dest='all_alignments',
        action='store_true',
        help='''Report all valid alignments of a read (very slow!).'''
    )
    parser.set_defaults(all_alignments=False)

    parser.add_argument(
        '-b', '--batch-size', dest='batch_size',
        type=int,
        default=100000,
        help='''Number of reads processed (mapped and merged) in one go per worker.
                The default (100000) works well for large indexes (e.g. human, mouse).
                Smaller indexes (e.g. yeast) will finish individual bowtie2 processes
                very quickly - set this number higher to spawn new processes 
                less frequently.
                '''
    )

    parser.add_argument(
        '--bowtie-parallel', dest='bowtie_parallel',
        action='store_true',
        help='''Use bowtie parallelisation rather than spawning multiple Bowtie2 processes.
                This is slower, but consumes potentially less memory.'''
    )
    parser.set_defaults(bowtie_parallel=False)

    parser.add_argument(
        '--split-fastq', dest='split_fastq',
        action='store_true',
        help='''Split FASTQ file into 10M chunks before mapping. Easier on tmp partitions.'''
    )
    parser.set_defaults(split_fastq=False)

    parser.add_argument(
        '--no-memory-map', dest='memory_map',
        action='store_false',
        help='''Do not map Bowtie2 index to memory. Only enable if you know what you are doing.'''
    )
    parser.set_defaults(memory_map=True)

    parser.add_argument(
        '--simple', dest='iterative',
        action='store_false',
        help='''Do not use iterative strategy (much faster, less sensitive).'''
    )
    parser.set_defaults(iterative=True)

    parser.add_argument(
        '--bwa', dest='bwa',
        action='store_true',
        help='''Use BWA as mapper.'''
    )
    parser.set_defaults(bwa=False)

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Copy original file to working directory (see -w option). Reduces network I/O.'''
    )
    parser.set_defaults(tmp=False)

    return parser


def map(argv):
    parser = map_parser()
    args = parser.parse_args(argv[2:])

    # check arguments
    input_files = args.input
    index_path = os.path.expanduser(args.index)
    output_folder = os.path.expanduser(args.output)

    step_size = args.step_size
    min_size = args.min_size
    trim_front = args.trim_front
    batch_size = args.batch_size
    min_quality = args.quality
    bowtie_parallel = args.bowtie_parallel
    memory_map = args.memory_map
    iterative = args.iterative
    restriction_enzyme = args.restriction_enzyme
    max_alignments = args.max_alignments
    all_alignments = args.all_alignments
    bwa = args.bwa
    tmp = args.tmp

    if bowtie_parallel:
        threads, bowtie_threads = 1, args.threads
    else:
        threads, bowtie_threads = args.threads, 1

    import kaic.mapping.map as map
    from kaic.tools.general import mkdir
    from kaic.tools.files import create_temporary_copy
    import subprocess
    import tempfile
    import shutil
    import glob

    additional_arguments = []
    if memory_map:
        additional_arguments += ['--mm']
    if all_alignments:
        additional_arguments += ['-a']
    elif max_alignments is not None:
        additional_arguments += ['-k', str(max_alignments)]

    index_dir = None
    try:
        if tmp:
            tmp = False
            index_dir = tempfile.mkdtemp()
            index_base = os.path.basename(index_path)
            if not bwa:
                for file_name in glob.glob(index_path + '*.bt2'):
                    shutil.copy(file_name, index_dir)
            else:
                for ending in ('amb', 'ann', 'bwt', 'pac', 'sa'):
                    file_name = index_path + '.{}'.format(ending)
                    shutil.copy(file_name, index_dir)

            index_path = os.path.join(index_dir, index_base)
            tmp = True

        if not bwa:
            if iterative:
                mapper = map.Bowtie2Mapper(index_path, min_quality=min_quality,
                                           additional_arguments=additional_arguments,
                                           threads=bowtie_threads)
            else:
                mapper = map.SimpleBowtie2Mapper(index_path, additional_arguments=additional_arguments,
                                                 threads=bowtie_threads)
        else:
            if iterative:
                mapper = map.BwaMapper(index_path, min_quality=min_quality,
                                       threads=bowtie_threads)
            else:
                mapper = map.SimpleBwaMapper(index_path,
                                             threads=bowtie_threads)

        for input_file in input_files:
            input_file = os.path.expanduser(input_file)
            if len(args.input) == 1 and not os.path.isdir(output_folder):
                output_file = output_folder
            else:
                output_folder = mkdir(output_folder)
                basename, extension = os.path.splitext(os.path.basename(input_file))
                if basename.endswith('.fastq'):
                    basename = basename[:-6]
                if basename.endswith('.fq'):
                    basename = basename[:-3]
                output_file = os.path.join(output_folder, basename + '.bam')

            if not args.split_fastq:
                original_output_file = output_file
                try:
                    if tmp:
                        tmp = False
                        input_file = create_temporary_copy(input_file, preserve_extension=True)
                        tmp_file = tempfile.NamedTemporaryFile(suffix=os.path.splitext(output_file)[1],
                                                               prefix='kaic_', delete=False)
                        tmp_file_name = tmp_file.name
                        output_file = tmp_file_name
                        tmp = True

                    logger.info("Starting mapping for {}".format(input_file))
                    map.iterative_mapping(input_file, output_file, mapper, threads=threads,
                                          min_size=min_size, step_size=step_size, batch_size=batch_size,
                                          trim_front=trim_front, restriction_enzyme=restriction_enzyme)
                finally:
                    if tmp:
                        os.remove(input_file)
                        shutil.copy(output_file, original_output_file)
                        os.remove(output_file)
            else:
                from kaic.tools.files import split_fastq, merge_sam, gzip_splitext

                logger.info("Splitting FASTQ files for mapping")
                if os.path.isdir(output_folder):
                    split_tmpdir = tempfile.mkdtemp(dir=output_folder)
                else:
                    split_tmpdir = tempfile.mkdtemp(dir=os.path.dirname(output_folder))
                try:
                    split_fastq_tmpdir = mkdir(os.path.join(split_tmpdir, 'fastq'))
                    split_sam_tmpdir = mkdir(os.path.join(split_tmpdir, 'sam'))
                    split_bam_files = []
                    split_fastq_results = []
                    for split_file in split_fastq(input_file, split_fastq_tmpdir):
                        basename = os.path.basename(gzip_splitext(split_file)[0])
                        split_bam_file = split_sam_tmpdir + '/{}.bam'.format(basename)
                        split_bam_files.append(split_bam_file)

                        split_command = ['kaic', 'map', split_file, index_path, split_bam_file,
                                         '-m', str(min_size), '-s', str(step_size), '-t', str(threads),
                                         '-q', str(args.quality), '-b', str(batch_size)]
                        if not iterative:
                            split_command += ['--simple']
                        if tmp:
                            split_command += ['-tmp']
                        if bowtie_parallel:
                            split_command += ['--bowtie-parallel']
                        if not memory_map:
                            split_command += ['--no-memory-map']
                        if trim_front:
                            split_command += ['--trim-front']
                        if restriction_enzyme is not None:
                            split_command += ['--restriction-enzyme', restriction_enzyme]

                        rt = subprocess.call(split_command)
                        split_fastq_results.append(rt)

                    for rt in split_fastq_results:
                        if rt != 0:
                            raise RuntimeError("Bowtie mapping had non-zero exit status")

                    logger.info("Merging BAM files into {}".format(output_file))
                    merge_sam(split_bam_files, output_file, tmp=args.copy)
                finally:
                    shutil.rmtree(split_tmpdir)
    finally:
        if tmp:
            shutil.rmtree(index_dir)


def iterative_mapping_parser():
    parser = argparse.ArgumentParser(
        prog="kaic iterative_mapping",
        description='Iteratively map a FASTQ file to a Bowtie 2 index'
    )

    parser.add_argument(
        'input',
        nargs='+',
        help='''File name of the input FASTQ file (or gzipped FASTQ)'''
    )

    parser.add_argument(
        'index',
        help='''Bowtie 2 genome index'''
    )

    parser.add_argument(
        'output',
        help='''Output file name (or folder name if multiple input files provided)'''
    )

    parser.add_argument(
        '-m', '--min-size', dest='min_size',
        type=int,
        default=25,
        help='''Minimum length of read before extension. Default 25.'''
    )

    parser.add_argument(
        '-s', '--step-size', dest='step_size',
        type=int,
        default=2,
        help='''Number of base pairs to extend at each round of mapping. Default is 2.'''
    )

    parser.add_argument(
        '-t', '--threads', dest='threads',
        type=int,
        default=1,
        help='''Number of threads used for mapping'''
    )

    parser.add_argument(
        '-q', '--quality', dest='quality',
        type=int,
        default=30,
        help='''Mapping quality cutoff for reads to be sent to another iteration'''
    )

    parser.add_argument(
        '-w', '--work-dir', dest='work_dir',
        help='''Working directory, defaults to the system temporary folder'''
    )

    parser.add_argument(
        '-r', '--restriction-enzyme', dest='restriction_enzyme',
        help='''Name of restriction enzyme used in experiment.
                                If provided, will trim reads at resulting ligation junction.'''
    )

    parser.add_argument(
        '-b', '--batch-size', dest='batch_size',
        type=int,
        default=1000000,
        help='''Number of reads processed (mapped and merged) in one go.'''
    )

    parser.add_argument(
        '--bowtie-parallel', dest='bowtie_parallel',
        action='store_true',
        help='''Use bowtie parallelisation rather than spawning multiple Bowtie2 processes.
                This is slower, but consumes a lot less memory.'''
    )
    parser.set_defaults(bowtie_parallel=False)

    parser.add_argument(
        '--split-fastq', dest='split_fastq',
        action='store_true',
        help='''Split FASTQ file into 10M chunks before mapping. Easier on tmp partitions.'''
    )
    parser.set_defaults(split_fastq=False)

    parser.add_argument(
        '--memory-map', dest='memory_map',
        action='store_true',
        help='''Map Bowtie2 index to memory (recommended if running on medium-memory systems and using many
                parallel threads. Use instead of --bowtie-parallel).'''
    )
    parser.set_defaults(memory_map=False)

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='copy',
        action='store_true',
        help='''Copy original file to working directory (see -w option). Reduces network I/O.'''
    )
    parser.set_defaults(copy=False)

    return parser


def iterative_mapping(argv):
    parser = iterative_mapping_parser()
    args = parser.parse_args(argv[2:])

    # check arguments
    index_path = os.path.expanduser(args.index)
    output_folder = os.path.expanduser(args.output)

    step_size = args.step_size
    min_size = args.min_size
    threads = args.threads
    batch_size = args.batch_size
    bowtie_parallel = args.bowtie_parallel
    memory_map = args.memory_map

    from kaic.mapping.iterative_mapping import split_iteratively_map_reads
    from kaic.tools.general import mkdir
    import subprocess

    for input_file in args.input:
        input_file = os.path.expanduser(input_file)
        if len(args.input) == 1 and not os.path.isdir(output_folder):
            output_file = output_folder
        else:
            output_folder = mkdir(output_folder)
            basename, extension = os.path.splitext(os.path.basename(input_file))
            output_file = output_folder + basename + '.bam'

        if not args.split_fastq:
            split_iteratively_map_reads(input_file, output_file, index_path, work_dir=args.work_dir,
                                        quality_cutoff=args.quality, batch_size=batch_size, threads=threads,
                                        min_size=min_size, step_size=step_size, copy=args.copy,
                                        restriction_enzyme=args.restriction_enzyme,
                                        adjust_batch_size=True,
                                        bowtie_parallel=bowtie_parallel,
                                        memory_map=memory_map)
        else:
            from kaic.tools.files import split_fastq, merge_sam, gzip_splitext

            logger.info("Splitting FASTQ files for mapping")
            if os.path.isdir(output_folder):
                split_tmpdir = tempfile.mkdtemp(dir=output_folder)
            else:
                split_tmpdir = tempfile.mkdtemp(dir=os.path.dirname(output_folder))
            try:
                split_fastq_tmpdir = mkdir(os.path.join(split_tmpdir, 'fastq'))
                split_sam_tmpdir = mkdir(os.path.join(split_tmpdir, 'sam'))
                split_bam_files = []
                split_fastq_results = []
                for split_file in split_fastq(input_file, split_fastq_tmpdir):
                    basename = os.path.basename(gzip_splitext(split_file)[0])
                    split_bam_file = split_sam_tmpdir + '/{}.bam'.format(basename)
                    split_bam_files.append(split_bam_file)

                    split_command = ['kaic', 'iterative_mapping', split_file, index_path, split_bam_file,
                                     '-m', str(min_size), '-s', str(step_size), '-t', str(threads),
                                     '-q', str(args.quality), '-b', str(batch_size)]
                    if args.restriction_enzyme is not None:
                        split_command += ['-r', args.restriction_enzyme]
                    if args.work_dir is not None:
                        split_command += ['-w', args.work_dir]
                    if args.copy:
                        split_command += ['-tmp']
                    if bowtie_parallel:
                        split_command += ['--bowtie-parallel']
                    if memory_map:
                        split_command += ['--memory-map']

                    rt = subprocess.call(split_command)
                    split_fastq_results.append(rt)

                for rt in split_fastq_results:
                    if rt != 0:
                        raise RuntimeError("Bowtie mapping had non-zero exit status")

                logger.info("Merging BAM files into {}".format(output_file))
                merge_sam(split_bam_files, output_file, tmp=args.copy)
            finally:
                shutil.rmtree(split_tmpdir)


def load_reads_parser():
    parser = argparse.ArgumentParser(
        prog="kaic load_reads",
        description='Load a SAM/BAM file into a Reads object'
    )

    parser.add_argument(
        'input',
        nargs='+',
        help='''Input SAM file'''
    )

    parser.add_argument(
        'output',
        help='''Output file'''
    )

    parser.add_argument(
        '-N', '--ignore-qname', dest='qname',
        action='store_false',
        help='''Do not store a read's qname, only a hashed version will be stored internally.'''
    )
    parser.set_defaults(qname=True)

    parser.add_argument(
        '-Q', '--ignore-qual', dest='qual',
        action='store_false',
        help='''Do not store a read's quality string.'''
    )
    parser.set_defaults(qual=True)

    parser.add_argument(
        '-S', '--ignore-seq', dest='seq',
        action='store_false',
        help='''Do not store a read's sequence string.'''
    )
    parser.set_defaults(seq=True)

    parser.add_argument(
        '-C', '--ignore-cigar', dest='cigar',
        action='store_false',
        help='''Do not store a read's cigar string. Warning: Some filters rely on this attribute.'''
    )
    parser.set_defaults(cigar=True)

    parser.add_argument(
        '-T', '--ignore-tags', dest='tags',
        action='store_false',
        help='''Do not store a read's tags. Warning: Some filters rely on this attribute.'''
    )
    parser.set_defaults(tags=True)

    parser.add_argument(
        '-D', '--ignore-default', dest='ignore_default',
        action='store_true',
        help='''Ignore qname, seq, and qual information to speed up read loading.'''
    )
    parser.set_defaults(ignore_default=False)

    parser.add_argument(
        '--split-sam', dest='split_sam',
        action='store_true',
        help='''Split SAM/BAM files into chunks of 10M alignments before loading. Useful in combination with tmp flag
                    to reduce tmp disk space usage.'''
    )
    parser.set_defaults(split_sam=False)

    parser.add_argument(
        '--append', dest='append',
        action='store_true',
        help='''Do not overwrite existing Reads object, but append to it instead.'''
    )
    parser.set_defaults(split_sam=False)

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)

    return parser


def load_reads(argv):
    parser = load_reads_parser()
    args = parser.parse_args(argv[2:])

    import kaic
    import glob

    input_paths = []
    for input_path in args.input:
        input_path = os.path.expanduser(input_path)
        if os.path.isdir(input_path):
            input_paths += glob.glob(os.path.join(input_path, '*.[sb]am'))
        else:
            input_paths.append(input_path)

    output_path = os.path.expanduser(args.output)
    original_output_path = output_path
    if args.tmp:
        tmp_file = tempfile.NamedTemporaryFile(delete=False)
        tmp_file.close()
        output_path = tmp_file.name
        logger.info("Output temporarily redirected to %s" % output_path)

    if args.ignore_default is True:
        store_qname = False
        store_seq = False
        store_qual = False
        store_cigar = True
        store_tags = True
    else:
        store_qname = args.qname
        store_seq = args.seq
        store_qual = args.qual
        store_cigar = args.cigar
        store_tags = args.tags

    mode = 'w' if not args.append else 'a'
    reads = kaic.Reads(file_name=output_path, mode=mode)

    from kaic.tools.files import split_sam

    if not args.split_sam:
        reads.load(sambam=input_paths, store_cigar=store_cigar, store_seq=store_seq, store_qname=store_qname,
                   store_qual=store_qual, store_tags=store_tags, sample_size=100000, tmp=args.tmp)
    else:
        output_folder = os.path.dirname(original_output_path)
        split_bams = []
        split_tmpdir = tempfile.mkdtemp(dir=output_folder)
        logger.info("Splitting SAM/BAM files before loading reads. Split directory: {}".format(split_tmpdir))
        try:
            for file_name in input_paths:
                split_bams += split_sam(file_name, split_tmpdir)
            reads.load(sambam=split_bams, store_cigar=store_cigar, store_seq=store_seq, store_qname=store_qname,
                       store_qual=store_qual, store_tags=store_tags, sample_size=100000, tmp=args.tmp)
        finally:
            shutil.rmtree(split_tmpdir)
    reads.close()

    if args.tmp:
        logger.info("Moving output file to destination...")
        shutil.move(output_path, original_output_path)
    logger.info("All done.")


def filter_reads_parser():
    parser = argparse.ArgumentParser(
        prog="kaic filter_reads",
        description='Filter a Reads object'
    )

    parser.add_argument(
        'input',
        help='''Input Reads file'''
    )

    parser.add_argument(
        'output',
        nargs="?",
        help='''Output Reads file. If not provided will filter existing file directly.'''
    )

    parser.add_argument(
        '-m', '--mapped', dest='mapped',
        action='store_true',
        help='''Filter unmapped reads'''
    )
    parser.set_defaults(mapped=False)

    parser.add_argument(
        '-u', '--unique', dest='unique',
        action='store_true',
        help='''Filter reads that map multiple times (with a lower score)'''
    )
    parser.set_defaults(unique=False)

    parser.add_argument(
        '-us', '--unique-strict', dest='unique_strict',
        action='store_true',
        help='''Strictly filter reads that map multiple times (XS tag)'''
    )
    parser.set_defaults(unique_strict=False)

    parser.add_argument(
        '-q', '--quality', dest='quality',
        type=int,
        help='''Cutoff for the minimum mapping quality of a read'''
    )

    parser.add_argument(
        '-c', '--contaminant', dest='contaminant',
        help='''A Reads file with contaminating reads. Will filter out reads with the same name.'''
    )

    parser.add_argument(
        '-s', '--stats', dest='stats',
        help='''Path for saving stats pdf'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)

    return parser


def filter_reads(argv):
    parser = filter_reads_parser()

    args = parser.parse_args(argv[2:])

    import kaic
    from kaic.tools.files import copy_or_expand, create_temporary_copy
    from kaic.construct.seq import ContaminantFilter

    # copy file if required
    original_input_path = os.path.expanduser(args.input)
    if args.tmp:
        logger.info("Copying Reads object to temporary file...")
        input_path = create_temporary_copy(original_input_path)
        logger.info("Temporarily working in %s" % input_path)
    else:
        input_path = copy_or_expand(args.input, args.output)

    reads = kaic.Reads(file_name=input_path, mode='a')

    if args.mapped:
        logger.info("Unmapped filter enabled")
        reads.filter_unmapped(queue=True)

    if args.unique_strict:
        logger.info("Strict multi-map filter enabled")
        reads.filter_non_unique(strict=True, queue=True)
    elif args.unique:
        logger.info("Soft multi-map filter enabled")
        reads.filter_non_unique(strict=False, queue=True)

    if args.quality:
        logger.info("Quality filter enabled (%d)" % args.quality)
        reads.filter_quality(args.quality, queue=True)

    if args.contaminant:
        contaminant_file = os.path.expanduser(args.contaminant)
        logger.info("Contaminant filter enabled %s" % contaminant_file)
        contaminant = kaic.Reads(contaminant_file)
        contaminant_filter = ContaminantFilter(contaminant,
                                               reads.add_mask_description("contaminant",
                                                                          "Filter contaminating reads"))
        reads.filter(contaminant_filter, queue=True)

    logger.info("Running filters...")
    reads.run_queued_filters(log_progress=True)
    logger.info("Done.")

    if args.stats:
        logger.info("Plotting filter statistics")
        from kaic.plotting.plot_statistics import plot_mask_statistics
        plot_mask_statistics(reads, reads._reads, output=args.stats)
        logger.info("Done.")

    reads.close()

    if args.tmp:
        output_path = os.path.expanduser(args.output)
        if os.path.isdir(output_path):
            output_path = "%s/%s" % (output_path, os.path.basename(original_input_path))
        logger.info("Moving temporary output file to destination %s..." % output_path)
        shutil.move(input_path, output_path)

    logger.info("All done.")


def build_genome_parser():
    parser = argparse.ArgumentParser(
        prog="kaic build_genome",
        description='Convenience command to build a Genome object'
    )

    parser.add_argument(
        'input',
        nargs='+',
        help=textwrap.dedent('''\
                             Can be a FASTA file,
                             a folder with FASTA files, or a
                             list of FASTA files.
                             ''')
    )

    parser.add_argument(
        'output',
        help='''Output file for Genome object'''
    )

    return parser


def build_genome(argv):
    parser = build_genome_parser()
    args = parser.parse_args(argv[2:])

    genome_string = ','.join(args.input)
    output_path = os.path.expanduser(args.output)

    import kaic

    logger.info("Building Genome...")
    genome = kaic.Genome.from_string(genome_string=genome_string, file_name=output_path)
    genome.close()
    logger.info("All done.")


def reads_to_pairs_parser():
    parser = argparse.ArgumentParser(
        prog="kaic reads_to_pairs",
        description='Convert a Reads object into a Pairs object'
    )

    parser.add_argument(
        'reads1',
        help='''First half of input reads'''
    )

    parser.add_argument(
        'reads2',
        help='''Second half of input reads'''
    )

    parser.add_argument(
        'genome',
        help=textwrap.dedent('''\
                                     Can be an HDF5 Genome object, a FASTA file,
                                     a folder with FASTA files, or a comma-separated
                                     list of FASTA files.
                                     ''')
    )

    parser.add_argument(
        'restriction_enzyme',
        help='''Restriction enzyme used in the experiment, e.g. HindIII'''
    )

    parser.add_argument(
        'output',
        help='''Output file for mapped pairs'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)

    return parser


def reads_to_pairs(argv):
    parser = reads_to_pairs_parser()
    args = parser.parse_args(argv[2:])

    import kaic
    from kaic.tools.files import create_temporary_copy

    reads1_path = os.path.expanduser(args.reads1)
    # copy file if required
    if args.tmp:
        logger.info("Creating temporary copy of first half of reads...")
        reads1_path = create_temporary_copy(reads1_path)
        logger.info("Working with temporary copy %s" % reads1_path)

    reads2_path = os.path.expanduser(args.reads2)
    # copy file if required
    if args.tmp:
        logger.info("Creating temporary copy of second half of reads...")
        reads2_path = create_temporary_copy(reads2_path)
        logger.info("Working with temporary copy %s" % reads2_path)

    genome_path = os.path.expanduser(args.genome)

    output_path = os.path.expanduser(args.output)
    original_output_path = output_path
    if args.tmp:
        tmp_file = tempfile.NamedTemporaryFile(delete=False)
        tmp_file.close()
        output_path = tmp_file.name
        logger.info("Working in temporary output file %s" % output_path)

    logger.info("Loading left side of reads...")
    reads1 = kaic.Reads(reads1_path, mode='r')
    logger.info("Loading right side of reads...")
    reads2 = kaic.Reads(reads2_path, mode='r')
    logger.info("Building genome...")
    genome = kaic.Genome.from_string(genome_path)
    logger.info("Getting fragments...")
    nodes = genome.get_regions(args.restriction_enzyme)

    logger.info("Building pairs...")
    pairs = kaic.Pairs(file_name=output_path, mode='w')
    logger.info("Mapping reads...")
    pairs.load(reads1, reads2, nodes)

    reads1.close()
    reads2.close()
    pairs.close()

    if args.tmp:
        logger.info("Removing temporary input files...")
        os.unlink(reads1_path)
        os.unlink(reads2_path)
        logger.info("Moving output file to destination %s" % original_output_path)
        shutil.move(output_path, original_output_path)

    logger.info("All done.")


def sort_sam_parser():
    parser = argparse.ArgumentParser(
        prog="kaic sort_sam",
        description="Convenience function to sort a SAM file by name "
                    "(exactly the same as 'samtools sort -n'!)"
    )

    parser.add_argument(
        'sam',
        help='''Input SAM/BAM'''
    )

    parser.add_argument(
        'output',
        nargs='?',
        help='''Output SAM/BAM. If not provided, will replace input 
                file with sorted version after sorting.'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)

    return parser


def sort_sam(argv):
    parser = sort_sam_parser()
    args = parser.parse_args(argv[2:])

    sam_file = os.path.expanduser(args.sam)
    output_file = None if args.output is None else os.path.expanduser(args.output)
    tmp = args.tmp

    from kaic.tools.files import sort_natural_sam, create_temporary_copy
    import tempfile
    import shutil

    success = False
    original_input_file = sam_file
    original_output_file = output_file
    try:
        if tmp:
            tmp = False
            sam_file = create_temporary_copy(sam_file)
            if output_file is not None:
                basename, extension = os.path.splitext(output_file)
                with tempfile.NamedTemporaryFile(delete=False, prefix='kaic_', suffix=extension) as f:
                    output_file = f.name
            tmp = True
            logger.info("Working in tmp: {}, ".format(sam_file, output_file))

        output_file = sort_natural_sam(sam_file, output_file)
        success = True
    finally:
        if tmp:
            if success:
                if original_output_file is not None:
                    shutil.copy(output_file, original_output_file)
                else:
                    shutil.copy(output_file, original_input_file)
            os.remove(sam_file)
            if os.path.normpath(sam_file) != os.path.normpath(output_file):
                os.remove(output_file)


def sam_to_pairs_parser():
    parser = argparse.ArgumentParser(
        prog="kaic sam_to_pairs",
        description='Convert a Reads object into a Pairs object'
    )

    parser.add_argument(
        'sam1',
        help='''First half of input reads'''
    )

    parser.add_argument(
        'sam2',
        help='''Second half of input reads'''
    )

    parser.add_argument(
        'genome',
        help="Can be an HDF5 Genome object, a FASTA file, "
             "a folder with FASTA files, or a "
             "comma-separated list of FASTA files."
    )

    parser.add_argument(
        'restriction_enzyme',
        help='''Restriction enzyme used in the experiment, e.g. HindIII'''
    )

    parser.add_argument(
        'output',
        help='''Output file for mapped pairs'''
    )

    parser.add_argument(
        '-m', '--mapped', dest='mapped',
        action='store_true',
        help='''Filter unmapped reads'''
    )
    parser.set_defaults(mapped=False)

    parser.add_argument(
        '-u', '--unique', dest='unique',
        action='store_true',
        help='''Filter reads that map multiple times (with a lower score)'''
    )
    parser.set_defaults(unique=False)

    parser.add_argument(
        '-us', '--unique-strict', dest='unique_strict',
        action='store_true',
        help='''Strictly filter reads that map multiple times (XS tag)'''
    )
    parser.set_defaults(unique_strict=False)

    parser.add_argument(
        '-q', '--quality', dest='quality',
        type=int,
        help='''Cutoff for the minimum mapping quality of a read'''
    )

    parser.add_argument(
        '-c', '--contaminant', dest='contaminant',
        help='''A Reads file with contaminating reads. Will filter out reads with the same name.'''
    )

    parser.add_argument(
        '-s', '--stats', dest='stats',
        help='''Path for saving stats pdf'''
    )

    parser.add_argument(
        '-S', '--no-check-sorted', dest='check_sorted',
        action='store_false',
        help='''Assume SAM files are sorted and do not check if that is actually the case'''
    )
    parser.set_defaults(check_sorted=True)

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)

    return parser


def sam_to_pairs(argv):
    parser = sam_to_pairs_parser()
    args = parser.parse_args(argv[2:])

    import kaic
    from kaic.tools.files import create_temporary_copy
    from kaic.data.general import Mask
    from kaic.construct.seq import SamBamReadPairGenerator, ReadPairs, \
        UnmappedFilter, UniquenessFilter, QualityFilter, ContaminantFilter
    import tempfile
    import shutil

    sam1_file = os.path.expanduser(args.sam1)
    sam2_file = os.path.expanduser(args.sam2)
    genome_file = os.path.expanduser(args.genome)
    restriction_enzyme = args.restriction_enzyme
    output_file = args.output
    filter_mapped = args.mapped
    filter_unique = args.unique
    filter_unique_strict = args.unique_strict
    filter_quality = args.quality
    filter_contaminant = args.contaminant
    stats_file = args.stats
    check_sorted = args.check_sorted
    tmp = args.tmp

    logger.info("Preparing filters...")
    filters = []
    if filter_mapped:
        f = UnmappedFilter(mask=Mask(ix=0, name='unmapped'))
        filters.append(f)

    if filter_unique or filter_unique_strict:
        f = UniquenessFilter(strict=filter_unique_strict, mask=Mask(ix=1, name='unique'))
        filters.append(f)

    if filter_quality is not None:
        f = QualityFilter(filter_quality, mask=Mask(ix=2, name='quality'))
        filters.append(f)

    if filter_contaminant is not None:
        f = ContaminantFilter(filter_contaminant, mask=Mask(ix=3, name='contaminant'))
        filters.append(f)

    original_output_file = output_file
    try:
        if tmp:
            tmp = False
            sam1_file = create_temporary_copy(sam1_file)
            sam2_file = create_temporary_copy(sam2_file)
            genome_file = create_temporary_copy(genome_file)
            with tempfile.NamedTemporaryFile(prefix='kaic_pairs_', suffix='.pairs') as f:
                output_file = f.name
            tmp = True

        logger.info("Getting regions")
        genome = kaic.Genome.from_string(genome_file)
        regions = genome.get_regions(restriction_enzyme)
        genome.close()

        sb = SamBamReadPairGenerator(sam1_file, sam2_file, check_sorted=check_sorted)
        for f in filters:
            sb.add_filter(f)

        pairs = ReadPairs(file_name=output_file, mode='w')

        pairs.add_regions(regions)
        pairs.add_read_pairs(sb)

        statistics = pairs.filter_statistics()
        pairs.close()
        logger.info("Done creating pairs.")

        if stats_file is not None:
            logger.info("Saving statistics...")
            import matplotlib
            matplotlib.use('agg')
            import matplotlib.pyplot as plt
            from kaic.plotting.plot_statistics import statistics_plot
            stats_file = os.path.expanduser(stats_file)
            fig, ax = plt.subplots()
            statistics_plot(statistics)
            fig.savefig(stats_file)
            plt.close(fig)
    finally:
        if tmp:
            shutil.copy(output_file, original_output_file)
            os.remove(output_file)

    logger.info("All done.")


def filter_pairs_parser():
    parser = argparse.ArgumentParser(
        prog="kaic filter_pairs",
        description='Filter a Pairs object'
    )

    parser.add_argument(
        'input',
        help='''Input FragmentMappedPairs file'''
    )

    parser.add_argument(
        'output',
        nargs="?",
        help='''Output FragmentMappedPairs file. If not provided will filter input file in place.'''
    )

    parser.add_argument(
        '-i', '--inward', dest='inward',
        type=int,
        help='''Minimum distance for inward-facing read pairs'''
    )

    parser.add_argument(
        '-o', '--outward', dest='outward',
        type=int,
        help='''Minimum distance for outward-facing read pairs'''
    )

    parser.add_argument(
        '--auto',
        action='store_true',
        help='''Auto-guess settings for inward/outward read pair filters.
                        Overrides --outward and --inward if set.'''
    )
    parser.set_defaults(auto=False)

    parser.add_argument(
        '-r', '--re-distance', dest='redist',
        type=int,
        help='''Maximum distance for a read to the nearest restriction site'''
    )

    parser.add_argument(
        '-l', '--self-ligated', dest='self_ligated',
        action='store_true',
        help='''Remove read pairs representing self-ligated fragments'''
    )
    parser.set_defaults(self_ligated=True)

    parser.add_argument(
        '-d', '--duplicate', dest='dup_thresh',
        type=int,
        help='''If specified, filter read pairs for PCR duplicates. Parameter determines
                        distance between alignment starts below which they are considered Starting
                        at same position. Sensible values are between 1 and 5.'''
    )

    parser.add_argument(
        '-s', '--stats', dest='stats',
        help='''Path for saving stats pdf'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)

    return parser


def filter_pairs(argv):
    parser = filter_pairs_parser()

    args = parser.parse_args(argv[2:])

    import kaic
    from kaic.tools.files import copy_or_expand, create_temporary_copy

    # copy file if required
    original_input_path = os.path.expanduser(args.input)
    if args.tmp:
        logger.info("Creating temporary copy of input file...")
        input_path = create_temporary_copy(original_input_path)
        logger.info("Working from copy in %s" % input_path)
    else:
        input_path = copy_or_expand(args.input, args.output)

    pairs = kaic.load(file_name=input_path, mode='a')

    if args.auto:
        logger.info("Filtering inward- and outward-facing reads using automatically"
                    "determined thresholds.")
        pairs.filter_ligation_products(queue=True)
    else:
        if args.inward:
            logger.info("Filtering inward-facing reads at %dbp" % args.inward)
            pairs.filter_inward(minimum_distance=args.inward, queue=True)

        if args.outward:
            logger.info("Filtering outward-facing reads at %dbp" % args.outward)
            pairs.filter_outward(minimum_distance=args.outward, queue=True)

    if args.redist:
        logger.info("Filtering reads with RE distance >%dbp" % args.redist)
        pairs.filter_re_dist(args.redist, queue=True)

    if args.self_ligated:
        logger.info("Filtering self-ligated read pairs")
        pairs.filter_self_ligated(queue=True)

    if args.dup_thresh:
        logger.info("Filtering PCR duplicates, threshold <=%dbp" % args.dup_thresh)
        pairs.filter_pcr_duplicates(threshold=args.dup_thresh, queue=True)

    logger.info("Running filters...")
    pairs.run_queued_filters(log_progress=True)
    logger.info("Done.")

    statistics = pairs.filter_statistics()
    pairs.close()
    logger.info("Done creating pairs.")

    if args.stats is not None:
        logger.info("Saving statistics...")
        import matplotlib
        matplotlib.use('agg')
        import matplotlib.pyplot as plt
        from kaic.plotting.plot_statistics import statistics_plot
        stats_file = os.path.expanduser(args.stats)
        fig, ax = plt.subplots()
        statistics_plot(statistics)
        fig.savefig(stats_file)
        plt.close(fig)

    if args.tmp:
        output_path = os.path.expanduser(args.output)
        if os.path.isdir(output_path):
            output_path = "%s/%s" % (output_path, os.path.basename(original_input_path))
        logger.info("Moving temporary output file to destination %s..." % output_path)
        shutil.move(input_path, output_path)

    logger.info("All done.")


def pairs_from_hicpro_parser():
    parser = argparse.ArgumentParser(
        prog="kaic pairs_from_hicpro",
        description='Import a Pairs object from HiC-Pro'
    )

    parser.add_argument(
        'input',
        help='''Valid pairs file from HiC-Pro. Format:
                name\\tchr1\\tpos1\\tstrand1\\tchr2\\tpos2\\tstrand2'''
    )

    parser.add_argument(
        'genome',
        help='''Path to genome (FASTA, folder with FASTA, hdf5 file)'''
    )

    parser.add_argument(
        'restriction_enzyme',
        help='''Name of restriction enzyme (case sensitive)'''
    )

    parser.add_argument(
        'output',
        help='''Output Pairs file.'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)
    return parser


def pairs_from_hicpro(argv):
    parser = pairs_from_hicpro_parser()

    args = parser.parse_args(argv[2:])
    input_file = os.path.expanduser(args.input)
    output_file = os.path.expanduser(args.output)
    genome_path = os.path.expanduser(args.genome)
    restriction_enzyme = args.restriction_enzyme
    tmp = args.tmp

    import kaic
    from kaic.construct.seq import TxtReadPairGenerator, ReadPairs
    from kaic.tools.files import create_temporary_copy

    original_input_file = input_file
    original_output_file = output_file
    try:
        if tmp:  # copy file if required
            tmp = False  # to prevent deleting input file should this be interrupted at this point
            logger.info("Copying input file...")
            input_file = create_temporary_copy(original_input_file)
            logger.info("New input file: {}".format(input_file))
            tmp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.pairs')
            tmp_file.close()
            output_file = tmp_file.name
            logger.info("Temporary output file: %s" % output_file)
            tmp = True

        logger.info("Getting regions")
        genome = kaic.Genome.from_string(genome_path)
        regions = genome.get_regions(restriction_enzyme)
        genome.close()

        sb = TxtReadPairGenerator(valid_pairs_file=input_file)
        pairs = ReadPairs(file_name=output_file, mode='w')

        pairs.add_regions(regions)
        pairs.add_read_pairs(sb)

        pairs.close()
        logger.info("Done creating pairs.")

    finally:
        if tmp:
            logger.info("Removing tmp files...")
            os.remove(input_file)
            shutil.copy(output_file, original_output_file)
            os.remove(output_file)

    logger.info("All done.")

def pairs_to_homer_parser():
    parser = argparse.ArgumentParser(
        prog="kaic pairs_to_homer",
        description='''Write pairs in Homer compatible "Hi-C Summary" format
                       http://homer.ucsd.edu/homer/interactions/HiCtagDirectory.html

                       Hi-C Summary Format (columns tab separated):
                       1. Read Name (can be blank)
                       2. chromosome for read 1
                       3. positions for read 1 (5' end of read, one-indexed)
                       4. strand of read 1 (+ or -)
                       5. chromosome for read 2
                       6. positions for read 2 (5' end of read, one-indexed)
                       7. strand of read 2 (+ or -)'''
    )

    parser.add_argument(
        'input',
        help='''Kaic pairs file'''
    )

    parser.add_argument(
        'output',
        help='''Path to output file. If extension is .gz will be gzipped.'''
    )

    parser.add_argument(
        '-i', '--include-filtered', dest='include_filtered',
        action='store_true',
        help='''Include filtered read pairs in output.'''
    )
    parser.set_defaults(include_filtered=False)

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)
    return parser


def pairs_to_homer(argv):
    parser = pairs_to_homer_parser()

    args = parser.parse_args(argv[2:])
    input_file = os.path.expanduser(args.input)
    output_file = os.path.expanduser(args.output)
    tmp = args.tmp
    gz = output_file.endswith(".gz")

    import kaic
    from kaic.tools.files import create_temporary_copy

    original_input_file = input_file
    original_output_file = output_file
    try:
        if tmp:  # copy file if required
            tmp = False  # to prevent deleting input file should this be interrupted at this point
            logger.info("Copying input file...")
            input_file = create_temporary_copy(original_input_file)
            logger.info("New input file: {}".format(input_file))
            tmp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.tsv.gz' if gz else '.tsv')
            tmp_file.close()
            output_file = tmp_file.name
            logger.info("Temporary output file: %s" % output_file)
            tmp = True

        pairs = kaic.load(input_file, mode="r")
        pairs.to_homer(output_file, include_filtered=args.include_filtered)
        pairs.close()
    finally:
        if tmp:
            logger.info("Removing tmp files...")
            os.remove(input_file)
            shutil.copy(output_file, original_output_file)
            os.remove(output_file)

    logger.info("All done.")

def pairs_to_hic_parser():
    parser = argparse.ArgumentParser(
        prog="kaic pairs_to_hic",
        description='Convert a pairs object into a Hic object'
    )

    parser.add_argument(
        'pairs',
        help='''Input Pairs file'''
    )

    parser.add_argument(
        'hic',
        help='''Output path for Hic file'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)
    return parser


def pairs_to_hic(argv):
    parser = pairs_to_hic_parser()
    args = parser.parse_args(argv[2:])

    import kaic
    from kaic.tools.files import create_temporary_copy

    original_pairs_path = os.path.expanduser(args.pairs)
    original_hic_path = os.path.expanduser(args.hic)
    if args.tmp:
        logger.info("Creating temporary copy of input file...")
        pairs_path = create_temporary_copy(original_pairs_path)
        logger.info("Working from temporary input file %s" % pairs_path)
        tmp_file = tempfile.NamedTemporaryFile(delete=False)
        tmp_file.close()
        hic_path = tmp_file.name
        logger.info("Working in temporary output file %s" % hic_path)
    else:
        pairs_path = original_pairs_path
        hic_path = original_hic_path

    logger.info("Loading pairs...")
    pairs = kaic.load(pairs_path, mode='r')
    logger.info("Done.")

    hic = pairs.to_hic(file_name=hic_path)
    logger.info("Done.")

    pairs.close()
    hic.close()

    if args.tmp:
        logger.info("Removing temporary input file...")
        os.unlink(pairs_path)
        logger.info("Moving output file to destination %s" % original_hic_path)
        shutil.move(hic_path, original_hic_path)

    logger.info("All done.")


def merge_hic_parser():
    parser = argparse.ArgumentParser(
        prog="kaic merge_hic",
        description='Merge multiple Hic objects'
    )

    parser.add_argument(
        'hic',
        nargs='+',
        help='''Input Hic files'''
    )

    parser.add_argument(
        'output',
        help='''Output binned Hic object'''
    )

    parser.add_argument(
        '-O', '--no-optimise', dest='optimise',
        action='store_false',
        help='''Produce a Hi-C object optimised for fast access times. May impact compatibility.'''
    )
    parser.set_defaults(optimise=True)

    parser.add_argument(
        '--intra', dest='intra',
        action='store_true',
        help='''Only merge intra-chromosomal contacts and set inter-chromosomal to 0.'''
    )
    parser.set_defaults(intra=False)

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)

    return parser


def merge_hic(argv):
    parser = merge_hic_parser()
    args = parser.parse_args(argv[2:])
    import tempfile
    tmpdir = tempfile.gettempdir() if args.tmp else None

    import kaic

    output_path = os.path.expanduser(args.output)
    paths = [os.path.expanduser(path) for path in args.hic]
    hics = [kaic.load_hic(path, mode='r', tmpdir=tmpdir) for path in paths]

    # try fast, basic loading first:
    try:
        if args.optimise:
            merged = kaic.AccessOptimisedHic.from_hic(hics, file_name=output_path, tmpdir=tmpdir,
                                                      only_intrachromosomal=args.intra)
        else:
            merged = kaic.Hic.from_hic(hics, file_name=output_path, tmpdir=tmpdir,
                                       only_intrachromosomal=args.intra)

        merged.close()
    except ValueError:
        logging.warning("The regions in your Hi-C objects do not appear to be identical. This will slow down"
                        "merging significantly.")
        first_hic = hics.pop(0)

        if args.optimise:
            merged = kaic.AccessOptimisedHic(data=first_hic, file_name=output_path, mode='w', tmpdir=tmpdir)
        else:
            merged = kaic.Hic(data=first_hic, file_name=output_path, mode='w', tmpdir=tmpdir)

        merged.merge(hics)

        merged.close()

    for hic in hics:
        hic.close()

    logger.info("All done")


def filter_hic_parser():
    parser = argparse.ArgumentParser(
        prog="kaic filter_hic",
        description='Filter a Hic object'
    )

    parser.add_argument(
        'input',
        help='''Input Hic file'''
    )

    parser.add_argument(
        'output',
        nargs="?",
        help='''Output Hic file. If not provided will filter input file in place.'''
    )

    parser.add_argument(
        '-l', '--low-coverage', dest='low',
        type=float,
        help='''Filter bins with "low coverage" (lower than specified absolute contact threshold)'''
    )

    parser.add_argument(
        '-rl', '--relative-low-coverage', dest='rel_low',
        type=float,
        help='''Filter bins using a relative "low coverage" threshold
                    (lower than the specified fraction of the median contact count)'''
    )

    parser.add_argument(
        '-ld', '--low-coverage-default', dest='low_default',
        action='store_true',
        help='''Filter bins with "low coverage" (under 10%% of median coverage for all non-zero bins)'''
    )
    parser.set_defaults(low_default=False)

    parser.add_argument(
        '-d', '--diagonal', dest='diagonal',
        type=int,
        help='''Filter bins along the diagonal up to this specified distance'''
    )

    parser.add_argument(
        '-s', '--stats', dest='stats',
        help='''Path for saving stats pdf'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)

    return parser


def filter_hic(argv):
    parser = filter_hic_parser()
    args = parser.parse_args(argv[2:])
    import kaic
    from kaic.tools.files import copy_or_expand, create_temporary_copy

    original_input_path = os.path.expanduser(args.input)
    if args.tmp:
        logger.info("Creating temporary copy of input file...")
        input_path = create_temporary_copy(original_input_path)
        logger.info("Working from copy in %s" % input_path)
    else:
        input_path = copy_or_expand(args.input, args.output)

    with kaic.load_hic(file_name=input_path, mode='a') as hic:
        if args.low_default:
            if args.low or args.rel_low:
                logger.info("Already specified an cutoff with -l or -rl, skipping -ld...")
            else:
                logger.info("Filtering low-coverage bins at 10%%")
                hic.filter_low_coverage_regions(rel_cutoff=0.1, cutoff=None, queue=True)

        if (args.low is not None or args.rel_low is not None) and args.low_default is False:
            logger.info("Filtering low-coverage bins using absolute cutoff {:.4}, "
                        "relative cutoff {:.1%}".format(float(args.low) if args.low else 0.,
                                                        float(args.rel_low) if args.rel_low else 0.))
            hic.filter_low_coverage_regions(cutoff=args.low, rel_cutoff=args.rel_low, queue=True)

            if args.diagonal is not None:
                logger.info("Filtering diagonal at distance %d" % args.diagonal)
                hic.filter_diagonal(distance=args.diagonal, queue=True)

            logger.info("Running filters...")
            hic.run_queued_filters(log_progress=True)
            logger.info("Done.")

            if args.stats:
                logger.info("Plotting filter statistics")
                from kaic.plotting.plot_statistics import plot_mask_statistics
                plot_mask_statistics(hic, hic._edges, output=args.stats)
                logger.info("Done.")

    if args.tmp:
        output_path = os.path.expanduser(args.output)
        if os.path.isdir(output_path):
            output_path = "%s/%s" % (output_path, os.path.basename(original_input_path))
        logger.info("Moving temporary output file to destination %s..." % output_path)
        shutil.move(input_path, output_path)

    logger.info("All done.")


def bin_hic_parser():
    parser = argparse.ArgumentParser(
        prog="kaic bin_hic",
        description='Bin a Hic object into same-size regions'
    )

    parser.add_argument(
        'hic',
        help='''Input Hic file'''
    )

    parser.add_argument(
        'output',
        help='''Output binned Hic object'''
    )

    parser.add_argument(
        'bin_size',
        type=int,
        help='''Bin size in base pairs'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)
    return parser


def bin_hic(argv):
    parser = bin_hic_parser()

    args = parser.parse_args(argv[2:])

    import kaic
    from kaic.tools.files import create_temporary_copy

    original_output_path = os.path.expanduser(args.output)
    if args.tmp:
        input_path = create_temporary_copy(args.hic)
        logger.info("Working in temporary directory...")
        tmp_file = tempfile.NamedTemporaryFile(delete=False)
        tmp_file.close()
        output_path = tmp_file.name
        logger.info("Temporary output file: %s" % output_path)
    else:
        input_path = os.path.expanduser(args.hic)
        output_path = original_output_path

    hic = kaic.load_hic(file_name=input_path, mode='r')

    logger.info("Binning at %dbp" % args.bin_size)
    binned = hic.bin(args.bin_size, file_name=output_path)

    hic.close()
    binned.close()

    if args.tmp:
        logger.info("Moving temporary output file to destination %s" % original_output_path)
        os.unlink(input_path)
        shutil.move(output_path, original_output_path)

    logger.info("All done.")


def correct_hic_parser():
    parser = argparse.ArgumentParser(
        prog="kaic correct_hic",
        description='Correct a Hic object for biases'
    )

    parser.add_argument(
        'input',
        help='''Input Hic file'''
    )

    parser.add_argument(
        'output',
        help='''Output Hic file.'''
    )

    parser.add_argument(
        '-i', '--ice', dest='ice',
        action='store_true',
        help='''Use ICE iterative correction instead of Knight matrix balancing'''
    )
    parser.set_defaults(ice=False)

    parser.add_argument(
        '-c', '--chromosome', dest='chromosome',
        action='store_true',
        help='''Correct intra-chromosomal data individually, ignore inter-chromosomal data'''
    )
    parser.set_defaults(chromosome=False)

    parser.add_argument(
        '-O', '--no-optimise', dest='optimise',
        action='store_false',
        help='''Produce a Hi-C object optimised for fast access times. May impact compatibility.'''
    )
    parser.set_defaults(optimise=True)

    parser.add_argument(
        '-r', '--restore-coverage', dest='restore_coverage',
        action='store_true',
        help='''Restore coverage to the original total number of reads. 
                    Otherwise matrix entries will be contact probabilities.
                    Only available for KR matrix balancing.'''
    )
    parser.set_defaults(restore_coverage=False)

    parser.add_argument(
        '--only-inter', dest='only_inter',
        action='store_true',
        help='''Only correct inter-chromosomal contacts. Sets intra-chromosomal contacts to 0.
                Always uses whole-matrix balancing, therefore incompatible with -c.'''
    )
    parser.set_defaults(only_inter=False)

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)
    return parser


def correct_hic(argv):
    parser = correct_hic_parser()

    args = parser.parse_args(argv[2:])

    if args.only_inter and args.chromosome:
        raise RuntimeError("--only-inter incompatible with -c! Aborting...")

    import kaic
    from kaic.tools.files import create_temporary_copy

    # copy file if required
    original_input_path = os.path.expanduser(args.input)
    original_output_path = os.path.expanduser(args.output)
    if args.tmp:
        logger.info("Copying data to temporary file...")
        input_path = create_temporary_copy(original_input_path)
        logger.info("Working from temporary file %s" % input_path)
        tmp_file = tempfile.NamedTemporaryFile(delete=False)
        tmp_file.close()
        output_path = tmp_file.name
        logger.info("Temporary output file: %s" % output_path)
    else:
        input_path = os.path.expanduser(args.input)
        output_path = os.path.expanduser(args.output)

    if args.ice:
        import kaic.correcting.ice_matrix_balancing as ice

        if args.restore_coverage:
            raise ValueError("Restoring coverage (-r) only available for KR balancing!")

        hic = kaic.load_hic(file_name=input_path, mode='r')
        ice.correct(hic, whole_matrix=not args.chromosome, copy=True, optimise=args.optimise,
                    file_name=output_path, intra_chromosomal=not args.only_inter)
        hic.close()
        if args.tmp:
            logger.info("Moving temporary output file to destination %s" % original_output_path)
            os.unlink(input_path)
            shutil.move(output_path, original_output_path)
    else:
        import kaic.correcting.knight_matrix_balancing as knight

        hic = kaic.load_hic(file_name=input_path, mode='r')
        hic_new = knight.correct(hic, whole_matrix=not args.chromosome,
                                 copy=True, file_name=output_path, optimise=args.optimise,
                                 restore_coverage=args.restore_coverage,
                                 intra_chromosomal=not args.only_inter)
        hic.close()
        hic_new.close()
        if args.tmp:
            logger.info("Moving temporary output file to destination %s" % original_output_path)
            os.unlink(input_path)
            shutil.move(output_path, original_output_path)

    logger.info("All done.")


def hic_from_juicer_parser():
    parser = argparse.ArgumentParser(
        prog="kaic hic_from_juicer",
        description='Import a Hi-C object from juicer (Aiden lab)'
    )

    parser.add_argument(
        'input',
        help='''Input .hic file, juicer format'''
    )

    parser.add_argument(
        'jar',
        help='''path to juicer_tools jar file'''
    )

    parser.add_argument(
        'genome',
        help='''Path to genome (FASTA, folder with FASTA, hdf5 file)'''
    )

    parser.add_argument(
        'resolution',
        type=int,
        help='''Resolution in base pairs'''
    )

    parser.add_argument(
        'output',
        help='''Output Hic file.'''
    )

    parser.add_argument(
        '-c', '--chromosomes', dest='chromosomes',
        nargs='+',
        help='''List of chromosomes to extract. Extracts all chromosomes in genome by default.'''
    )
    parser.set_defaults(ice=False)

    parser.add_argument(
        '--no-inter-chromosomal', dest='inter_chromosomal',
        action='store_false',
        help='''Do not extract inter-chromosomal matrices'''
    )
    parser.set_defaults(inter_chromosomal=True)

    parser.add_argument(
        '--juicer-norm', dest='juicer_norm',
        default='NONE',
        help='''Juicer normalisation method. Default: NONE, 
                see juicer documentation for alternatives.'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)
    return parser


def hic_from_juicer(argv):
    parser = hic_from_juicer_parser()

    args = parser.parse_args(argv[2:])
    input_file = os.path.expanduser(args.input)
    output_file = os.path.expanduser(args.output)
    genome_path = os.path.expanduser(args.genome)
    jar_path = args.jar
    resolution = args.resolution
    chromosomes = args.chromosomes
    inter_chromosomal = args.inter_chromosomal
    juicer_norm = args.juicer_norm
    tmp = args.tmp

    import kaic
    from kaic.tools.files import create_temporary_copy

    original_input_file = input_file
    original_output_file = output_file
    try:
        if tmp:  # copy file if required
            tmp = False  # to prevent deleting input file should this be interrupted at this point
            logger.info("Copying input file...")
            input_file = create_temporary_copy(original_input_file)
            logger.info("New input file: {}".format(input_file))
            tmp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.hic')
            tmp_file.close()
            output_file = tmp_file.name
            logger.info("Temporary output file: %s" % output_file)
            tmp = True

        hic = kaic.AccessOptimisedHic.from_juicer(input_file, jar_path, genome_path, resolution,
                                                  norm=juicer_norm, output_file=output_file,
                                                  inter_chromosomal=inter_chromosomal,
                                                  chromosomes=chromosomes)
        hic.close()
    finally:
        if tmp:
            logger.info("Removing tmp files...")
            os.remove(input_file)
            shutil.copy(output_file, original_output_file)
            os.remove(output_file)

    logger.info("All done.")


def to_cooler_parser():
    parser = argparse.ArgumentParser(
        prog="kaic hic_to_cooler",
        description="""Convert a Hic file into cooler format.
                       See https://github.com/mirnylab/cooler for details.
                       If input Hi-C matrix is uncorrected, the uncorrected matrix is stored.
                       If it is corrected, the uncorrected matrix is stored and the bias vector.
                       Cooler always calculates corrected matrix on-the-fly from the uncorrected
                       matrix and the bias vector."""
    )

    parser.add_argument(
        'input',
        help='''Input .hic file, kaic format.'''
    )

    parser.add_argument(
        'output',
        help='''Output cooler file.'''
    )
    return parser


def hic_to_cooler(argv):
    parser = to_cooler_parser()

    args = parser.parse_args(argv[2:])
    input_file = os.path.expanduser(args.input)
    output_file = os.path.expanduser(args.output)

    import kaic
    try:
        import cooler
    except ImportError:
        logger.error("Cannot import cooler. Install cooler 'pip install cooler'.")
        raise

    hic = kaic.load_hic(input_file, mode='r')
    hic.to_cooler(output_file)
    logger.info("All done.")

def dump_parser():
    parser = argparse.ArgumentParser(
        prog="kaic dump",
        description='Dump Hic file to txt file(s).'
    )

    parser.add_argument(
        'hic',
        help='''Hic file'''
    )

    parser.add_argument(
        'matrix',
        nargs='?',
        help='''Output file for matrix entries. If not provided, will write to stdout.'''
    )

    parser.add_argument(
        'regions',
        nargs='?',
        help='''Output file for Hic regions. If not provided, will write regions into matrix file.'''
    )

    parser.add_argument(
        '-s', '--subset', dest='subset',
        help='''Only output this matrix subset. Format:
                <chr>[:<start>-<end>][--<chr>[:<start><end>]],
                e.g.:
                "chr1--chr1" to extract only the chromosome 1 submatrix;
                "chr2:3400000-4200000" to extract contacts of this region 
                on chromosome 2 to all other regions in the genome;'''
    )

    parser.add_argument(
        '-S', '--no-sparse', dest='sparse',
        action='store_false',
        help='''Store full, square matrix instead of sparse format.'''
    )
    parser.set_defaults(sparse=True)
    return parser


def dump(argv):
    parser = dump_parser()
    args = parser.parse_args(argv[2:])
    hic_file = os.path.expanduser(args.hic)
    output_matrix = None if args.matrix is None else os.path.expanduser(args.matrix)
    output_regions = None if args.regions is None else os.path.expanduser(args.regions)
    subset_string = args.subset
    sparse = args.sparse

    import kaic

    col_subset_region = None
    row_subset_region = None
    if subset_string is not None:
        col_subset_string = None
        if '--' in subset_string:
            row_subset_string, col_subset_string = subset_string.split('--')
        else:
            row_subset_string = subset_string
        row_subset_region = kaic.GenomicRegion.from_string(row_subset_string)
        if col_subset_string is not None:
            col_subset_region = kaic.GenomicRegion.from_string(col_subset_string)

    hic = kaic.load(hic_file, mode='r')

    row_regions = [region for region in hic.regions(row_subset_region)]
    col_regions = [region for region in hic.regions(col_subset_region)]
    row_regions_dict = {region.ix: (region, i) for i, region in enumerate(row_regions)}
    col_regions_dict = {region.ix: (region, i) for i, region in enumerate(col_regions)}

    if not sparse:
        if output_matrix is None or output_regions is None:
            raise ValueError("Cannot write matrix to stdout, must provide "
                             "both matrix and regions file for output")
        m = hic.as_matrix(key=(row_subset_region, col_subset_region))
        import numpy as np
        np.savetxt(output_matrix, m)
    else:
        if output_matrix is None:
            for edge in hic.edges(key=(row_subset_region, col_subset_region), lazy=True):
                source, i = row_regions_dict[edge.source]
                sink, j = col_regions_dict[edge.sink]
                print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                    source.chromosome, source.start, source.end,
                    sink.chromosome, sink.start, sink.end,
                    edge.weight
                ))
        else:
            with open(output_matrix, 'w') as o:
                if output_regions is None:
                    for edge in hic.edges(key=(row_subset_region, col_subset_region), lazy=True):
                        source, i = row_regions_dict[edge.source]
                        sink, j = col_regions_dict[edge.sink]
                        o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            source.chromosome, source.start, source.end,
                            sink.chromosome, sink.start, sink.end,
                            edge.weight
                        ))
                else:
                    for edge in hic.edges(key=(row_subset_region, col_subset_region), lazy=True):
                        source, i = row_regions_dict[edge.source]
                        sink, j = col_regions_dict[edge.sink]
                        o.write("{}\t{}\t{}\n".format(
                            i, j, edge.weight
                        ))

    # write regions to file
    if output_regions is not None:
        if row_subset_region == col_subset_region:
            with open(output_regions, 'w') as o:
                for region in row_regions:
                    o.write("{}\t{}\t{}\n".format(region.chromosome, region.start, region.end))
        else:
            basepath, extension = os.path.splitext(output_regions)
            with open(basepath + '_row' + extension, 'w') as o:
                for region in row_regions:
                    o.write("{}\t{}\t{}\n".format(region.chromosome, region.start, region.end))
            with open(basepath + '_col' + extension, 'w') as o:
                for region in col_regions:
                    o.write("{}\t{}\t{}\n".format(region.chromosome, region.start, region.end))

    logger.info("All done.")


def hic_pca_parser():
    parser = argparse.ArgumentParser(
        prog="kaic hic_pca",
        description='Do a PCA on multiple Hi-C objects'
    )

    parser.add_argument(
        'input',
        nargs='+',
        help='''Input Hic files'''
    )

    parser.add_argument(
        'output_folder',
        help='''Output folder for PCA results.'''
    )

    parser.add_argument(
        '-s', '--sample-sizes', dest='sample_sizes',
        nargs='+',
        type=int,
        default=[50000],
        help='''Sample sizes for contacts to do the PCA on.'''
    )

    parser.add_argument(
        '-i', '--intra', dest='intra',
        action='store_true',
        help='''Only do PCA on intra-chromosomal contacts'''
    )
    parser.set_defaults(intra=False)

    parser.add_argument(
        '-d', '--divide', dest='divide',
        action='store_true',
        help='''Divide PCAs into individual chromosomes'''
    )
    parser.set_defaults(divide=False)

    parser.add_argument(
        '-e', '--expected-filter', dest='expected_filter',
        type=float,
        help='''Cutoff for expected/observed ratio of a contact to be considered for PCA. Default: no filter.'''
    )

    parser.add_argument(
        '-b', '--background-filter', dest='background_filter',
        type=float,
        help='''Cutoff for ratio of average inter-chromosomal to
                    observed contact to be considered for PCA. Default: no filter.'''
    )

    parser.add_argument(
        '-w', '--window-filter', dest='window_filter',
        nargs=2,
        type=int,
        help='''Min and max values in base pairs defining a window of
                    contact distances that are retained for analysis.'''
    )

    parser.add_argument(
        '-n', '--names', dest='names',
        nargs='+',
        help='''Sample names for plot labelling.'''
    )

    parser.add_argument(
        '-p', '--pair-selection', dest='pair_selection',
        default='variance',
        help='''Mechanism to select pairs from Hi-C matrix. Default: variance.
                    Possible values are:
                    variance: Selects pairs with the largest variance across samples first.
                    fc: Select pairs with the largest fold-change across samples first.
                    passthrough: Selects pairs without preference.
                 '''
    )

    parser.add_argument(
        '-c', '--colors', dest='colors',
        nargs='+',
        help='''Colors for plotting.'''
    )

    parser.add_argument(
        '-m', '--markers', dest='markers',
        nargs='+',
        help='''Markers for plotting. Follows Matplotlib marker
                    definitions: http://matplotlib.org/api/markers_api.html'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)
    return parser


def hic_pca(argv):
    parser = hic_pca_parser()

    args = parser.parse_args(argv[2:])

    import errno
    import matplotlib
    matplotlib.use('pdf')
    import kaic
    from kaic.plotting.plot_statistics import pca_plot
    from kaic.architecture.pca import do_pca, HicCollectionWeightMeanVariance, PassthroughPairSelection, \
        LargestVariancePairSelection, LargestFoldChangePairSelection
    from kaic.architecture.hic_architecture import ExpectedObservedCollectionFilter,\
        BackgroundLigationCollectionFilter, MinMaxDistanceCollectionFilter
    from kaic.tools.files import create_temporary_copy
    import shutil

    hics = []
    output_folder = os.path.expanduser(args.output_folder)

    try:
        os.makedirs(output_folder)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    if args.tmp:
        tmp_file = tempfile.NamedTemporaryFile(delete=False)
        tmp_file.close()
        output_path = tmp_file.name
    else:
        output_path = output_folder + '/hics.coll'

    sample_names = args.names if args.names is not None else []
    try:
        if len(args.input) == 1:
            shutil.copy(args.input[0], output_path)
            coll = HicCollectionWeightMeanVariance(output_path)
        else:

            for file_name in args.input:
                if args.names is None:
                    sample_names.append(os.path.splitext(os.path.basename(file_name))[0])

                if args.tmp:
                    input_path = create_temporary_copy(file_name)
                else:
                    input_path = os.path.expanduser(file_name)
                hics.append(kaic.load_hic(input_path, mode='r'))

            coll = HicCollectionWeightMeanVariance(hics, file_name=output_path, only_intra_chromosomal=args.intra,
                                                   scale_libraries=True, mode='w')
            coll.calculate()

        # filtering
        if args.expected_filter is not None:
            eof = ExpectedObservedCollectionFilter(coll, fold_change=args.expected_filter)
            coll.filter(eof, queue=True)
        if args.background_filter is not None:
            bgf = BackgroundLigationCollectionFilter(coll, all_contacts=True, fold_change=args.background_filter)
            coll.filter(bgf, queue=True)
        if args.window_filter is not None:
            mmdf = MinMaxDistanceCollectionFilter(coll, min_distance=args.window_filter[0],
                                                  max_distance=args.window_filter[1])
            coll.filter(mmdf, queue=True)
        coll.run_queued_filters()

        if args.pair_selection == 'variance':
            pair_selector = LargestVariancePairSelection()
        elif args.pair_selection == 'fc':
            pair_selector = LargestFoldChangePairSelection()
        elif args.pair_selection == 'passthrough':
            pair_selector = PassthroughPairSelection()
        else:
            raise ValueError("Pair selection mechanism {} is not valid".format(args.pair_selection))

        if args.divide:
            regions = coll.chromosomes()
        else:
            regions = [None]

        for sample_size in args.sample_sizes:
            for region in regions:
                logger.info("Sample size: %d" % sample_size)
                pca_info, pca_res = do_pca(coll, pair_selection=pair_selector, sample_size=sample_size,
                                           regions=region)

                with open(output_folder + "/explained_variance_{}_{}.txt".format(region, sample_size), 'w') as var:
                    for i, variance in enumerate(pca_info.explained_variance_ratio_):
                        var.write(str(variance))
                        if i == len(pca_info.explained_variance_ratio_)-1:
                            var.write("\n")
                        else:
                            var.write("\t")

                with open(output_folder + "/pca_results_{}_{}.txt".format(region, sample_size), 'w') as res:
                    for i, row in enumerate(pca_res):
                        for j, value in enumerate(row):
                            res.write(str(value))
                            if j == len(row)-1:
                                res.write("\n")
                            else:
                                res.write("\t")

                fig, ax = pca_plot(pca_res, pca_info=pca_info, colors=args.colors, names=sample_names,
                                   markers=args.markers)
                ax.set_title("PCA sample size %d" % sample_size)
                fig.savefig(output_folder + "/pca_plot_{}_{}.pdf".format(region, sample_size))
                fig.clf()
    finally:
        if args.tmp:
            for hic in hics:
                file_name = hic.file.filename
                hic.close()
                os.unlink(file_name)
            shutil.move(output_path, output_folder + '/hics.coll')


def call_peaks_parser():
    parser = argparse.ArgumentParser(
        prog="kaic call_peaks",
        description='Call enriched peaks in a Hic object'
    )

    parser.add_argument(
        'input',
        help='''Input Hic file'''
    )

    parser.add_argument(
        'output',
        help='''Output HDF5 file'''
    )

    parser.add_argument(
        '-c', '--chromosomes', dest='chromosomes',
        nargs='+',
        help='''Chromosomes to be investigated.'''
    )

    parser.add_argument(
        '-p', '--peak-size', dest='peak_size',
        type=int,
        help='''Size of the expected peak in pixels. If not set, will be estimated to correspond to ~ 25kb.'''
    )

    parser.add_argument(
        '-w', '--width', dest='width',
        type=int,
        help='''Width of the investigated area surrounding a peak in pixels. If not set, will be estimated at p+3'''
    )

    parser.add_argument(
        '-m', '--min-dist', dest='min_dist',
        type=int,
        default=3,
        help='''Minimum distance in pixels for two loci to be considered as peaks. Default: 3'''
    )

    parser.add_argument(
        '-t', '--threads', dest='threads',
        type=int,
        default=4,
        help='''Number of threads for parallel processing. Default: 4'''
    )

    parser.add_argument(
        '-s', '--slice-size', dest='slice_size',
        type=int,
        default=200,
        help='''Width of submatrix examined per process. Default: 200'''
    )

    parser.add_argument(
        '--minimum_mappability', dest='minimum_mappability',
        type=float,
        default=0.7,
        help='''Minimum mappable fraction of a pixel neighborhood to consider pixel as peak'''
    )

    parser.add_argument(
        '-i', '--inter-chromosomal', dest='inter',
        action='store_true',
        help='''If set, also find peaks in inter-chromosomal data.'''
    )
    parser.set_defaults(inter=False)

    parser.add_argument(
        '--sge', dest='sge',
        action='store_true',
        help='''Run on SGE cluster'''
    )
    parser.set_defaults(sge=False)

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)
    return parser


def call_peaks(argv):
    parser = call_peaks_parser()

    args = parser.parse_args(argv[2:])

    sge = args.sge

    import kaic
    import kaic.data.network as kn
    from kaic.tools.files import create_temporary_copy

    # copy file if required
    original_input_path = os.path.expanduser(args.input)
    original_output_path = os.path.expanduser(args.output)
    if args.tmp:
        logger.info("Copying data to temporary file...")
        input_path = create_temporary_copy(original_input_path)
        logger.info("Working from temporary file %s" % input_path)
        tmp_file = tempfile.NamedTemporaryFile(delete=False)
        tmp_file.close()
        output_path = tmp_file.name
        logger.info("Temporary output file: %s" % output_path)
    else:
        input_path = os.path.expanduser(args.input)
        output_path = os.path.expanduser(args.output)

    pk = kn.RaoPeakCaller(p=args.peak_size, w_init=args.width, min_locus_dist=args.min_dist,
                          n_processes=args.threads, slice_size=args.slice_size,
                          process_inter=args.inter, cluster=sge,
                          min_mappable_fraction=args.minimum_mappability)

    hic = kaic.load_hic(input_path, mode='r')

    if args.chromosomes is None:
        chromosome_pairs = None
    else:
        chromosome_pairs = []
        for i in range(len(args.chromosomes)):
            chromosome1 = args.chromosomes[i]
            for j in range(i, len(args.chromosomes)):
                chromosome2 = args.chromosomes[j]
                chromosome_pairs.append((chromosome1, chromosome2))

    peaks = pk.call_peaks(hic, chromosome_pairs=chromosome_pairs, file_name=output_path)

    logger.info("Found %d potential peaks" % len(peaks))
    peaks.close()

    if args.tmp:
        os.unlink(input_path)
        logger.info("Moving temporary output file to destination %s" % original_output_path)
        shutil.move(output_path, original_output_path)


def filter_peaks_parser():
    parser = argparse.ArgumentParser(
        prog="kaic filter_peaks",
        description='Filter peaks called with call_peaks. By default, this is a pass-through filter.'
                    'Enable different filters using the below parameters.'
    )

    parser.add_argument(
        'input',
        help='''Input Peaks file'''
    )

    parser.add_argument(
        'output',
        nargs='?',
        help='''Output filtered Peaks file'''
    )

    parser.add_argument(
        '-f', '--fdr', dest='fdr_cutoff',
        type=float,
        help='''Global FDR cutoff for all neighborhoods. Value between 0 and 1.'''
    )

    parser.add_argument(
        '-fd', '--fdr-donut', dest='fdr_donut_cutoff',
        type=float,
        help='''Donut neighborhood FDR cutoff. Value between 0 and 1.'''
    )

    parser.add_argument(
        '-fh', '--fdr-horizontal', dest='fdr_horizontal_cutoff',
        type=float,
        help='''Horizontal neighborhood FDR cutoff. Value between 0 and 1.'''
    )

    parser.add_argument(
        '-fv', '--fdr-vertical', dest='fdr_vertical_cutoff',
        type=float,
        help='''Vertical neighborhood FDR cutoff. Value between 0 and 1.'''
    )

    parser.add_argument(
        '-fl', '--fdr-lower-left', dest='fdr_lower_left_cutoff',
        type=float,
        help='''Lower-left neighborhood FDR cutoff. Value between 0 and 1.'''
    )

    parser.add_argument(
        '-e', '--enrichment', dest='enrichment_cutoff',
        type=float,
        help='''Global enrichment cutoff. Value between 0 and infinity,
                    e.g. 2.0 means two-fold enrichment over every contact neighborhood.'''
    )

    parser.add_argument(
        '-ed', '--enrichment-donut', dest='enrichment_donut_cutoff',
        type=float,
        help='''Donut enrichment cutoff. Value between 0 and infinity.'''
    )

    parser.add_argument(
        '-eh', '--enrichment-horizontal', dest='enrichment_horizontal_cutoff',
        type=float,
        help='''Horizontal enrichment cutoff. Value between 0 and infinity.'''
    )

    parser.add_argument(
        '-ev', '--enrichment-vertical', dest='enrichment_vertical_cutoff',
        type=float,
        help='''Vertical enrichment cutoff. Value between 0 and infinity.'''
    )

    parser.add_argument(
        '-el', '--enrichment-lower_left', dest='enrichment_lower_left_cutoff',
        type=float,
        help='''Lower left enrichment cutoff. Value between 0 and infinity.'''
    )

    parser.add_argument(
        '-m', '--mappability', dest='mappability_cutoff',
        type=float,
        help='''Global mappability cutoff for all neighborhoods. Value between 0 and 1.'''
    )

    parser.add_argument(
        '-md', '--mappability-donut', dest='mappability_donut_cutoff',
        type=float,
        help='''Donut neighborhood mappability cutoff. Value between 0 and 1.'''
    )

    parser.add_argument(
        '-mh', '--mappability-horizontal', dest='mappability_horizontal_cutoff',
        type=float,
        help='''Horizontal neighborhood mappability cutoff. Value between 0 and 1.'''
    )

    parser.add_argument(
        '-mv', '--mappability-vertical', dest='mappability_vertical_cutoff',
        type=float,
        help='''Vertical neighborhood mappability cutoff. Value between 0 and 1.'''
    )

    parser.add_argument(
        '-ml', '--mappability-lower-left', dest='mappability_lower_left_cutoff',
        type=float,
        help='''Lower-left neighborhood mappability cutoff. Value between 0 and 1.'''
    )

    parser.add_argument(
        '-o', '--observed', dest='observed',
        type=int,
        help='''Minimum observed value (integer, uncorrected). Default: 1'''
    )

    parser.add_argument(
        '-r', '--rao', dest='rao',
        action='store_true',
        help='''Filter peaks as Rao et al. (2014) does. It only retains peaks that

                        1. are at least 2-fold enriched over either the donut or lower-left neighborhood
                        2. are at least 1.5-fold enriched over the horizontal and vertical neighborhoods
                        3. are at least 1.75-fold enriched over both the donut and lower-left neighborhood
                        4. have an FDR <= 0.1 in every neighborhood
            '''
    )
    parser.set_defaults(rao=False)

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)
    return parser


def filter_peaks(argv):
    parser = filter_peaks_parser()

    args = parser.parse_args(argv[2:])

    rao = args.rao
    fdr_cutoff = args.fdr_cutoff
    fdr_cutoff_ll = args.fdr_lower_left_cutoff
    fdr_cutoff_h = args.fdr_horizontal_cutoff
    fdr_cutoff_v = args.fdr_vertical_cutoff
    fdr_cutoff_d = args.fdr_donut_cutoff
    mappability_cutoff = args.mappability_cutoff
    mappability_cutoff_ll = args.mappability_lower_left_cutoff
    mappability_cutoff_h = args.mappability_horizontal_cutoff
    mappability_cutoff_v = args.mappability_vertical_cutoff
    mappability_cutoff_d = args.mappability_donut_cutoff
    enrichment_cutoff = args.enrichment_cutoff
    enrichment_cutoff_ll = args.enrichment_lower_left_cutoff
    enrichment_cutoff_h = args.enrichment_horizontal_cutoff
    enrichment_cutoff_v = args.enrichment_vertical_cutoff
    enrichment_cutoff_d = args.enrichment_donut_cutoff
    observed_cutoff = args.observed
    tmp = args.tmp

    import kaic.data.network as kn
    from kaic.tools.files import create_temporary_copy, copy_or_expand

    # copy file if required
    original_input_path = os.path.expanduser(args.input)
    if tmp:
        logger.info("Copying data to temporary file...")
        input_path = create_temporary_copy(original_input_path)
    else:
        input_path = copy_or_expand(args.input, args.output)

    peaks = kn.RaoPeakInfo(input_path, mode='a')

    if rao:
        logger.info("Rao filter")
        peaks.filter_rao(queue=True)

    if fdr_cutoff is not None:
        logger.info("Global FDR filter at {}".format(fdr_cutoff))
        peaks.filter_fdr(fdr_cutoff, queue=True)
    elif (fdr_cutoff_d is not None or fdr_cutoff_h is not None
          or fdr_cutoff_v is not None or fdr_cutoff_ll is not None):
        logger.info("Local FDR filter")
        fdr_mask = peaks.add_mask_description('fdr', 'FDR cutoff filter')
        fdr_filter = kn.FdrPeakFilter(fdr_ll_cutoff=fdr_cutoff_ll,
                                      fdr_d_cutoff=fdr_cutoff_d,
                                      fdr_h_cutoff=fdr_cutoff_h,
                                      fdr_v_cutoff=fdr_cutoff_v,
                                      mask=fdr_mask)
        peaks.filter(fdr_filter, queue=True)

    if mappability_cutoff is not None:
        logger.info("Global mappability filter at {}".format(mappability_cutoff))
        peaks.filter_mappability(mappability_cutoff, queue=True)
    elif (mappability_cutoff_d is not None or mappability_cutoff_h is not None
          or mappability_cutoff_v is not None or mappability_cutoff_ll is not None):
        logger.info("Local mappability filter")
        mappability_mask = peaks.add_mask_description('mappability', 'mappability cutoff filter')
        mappability_filter = kn.MappabilityPeakFilter(mappability_ll_cutoff=mappability_cutoff_ll,
                                                      mappability_d_cutoff=mappability_cutoff_d,
                                                      mappability_h_cutoff=mappability_cutoff_h,
                                                      mappability_v_cutoff=mappability_cutoff_v,
                                                      mask=mappability_mask)
        peaks.filter(mappability_filter, queue=True)

    if enrichment_cutoff is not None:
        logger.info("Global enrichment filter at {}".format(enrichment_cutoff))
        peaks.filter_enrichment(enrichment_cutoff, queue=True)
    elif (enrichment_cutoff_d is not None or enrichment_cutoff_h is not None
          or enrichment_cutoff_v is not None or enrichment_cutoff_ll is not None):
        logger.info("Local enrichment filter")
        enrichment_mask = peaks.add_mask_description('enrichment', 'enrichment cutoff filter')
        enrichment_filter = kn.EnrichmentPeakFilter(enrichment_ll_cutoff=enrichment_cutoff_ll,
                                                    enrichment_d_cutoff=enrichment_cutoff_d,
                                                    enrichment_h_cutoff=enrichment_cutoff_h,
                                                    enrichment_v_cutoff=enrichment_cutoff_v,
                                                    mask=enrichment_mask)
        peaks.filter(enrichment_filter, queue=True)

    if observed_cutoff is not None:
        peaks.filter_observed(observed_cutoff, queue=True)

    peaks.run_queued_filters()
    peaks.close()

    if args.tmp:
        output_path = os.path.expanduser(args.output)
        if os.path.isdir(output_path):
            output_path = "%s/%s" % (output_path, os.path.basename(original_input_path))
        logger.info("Moving temporary output file to destination %s..." % output_path)
        shutil.move(input_path, output_path)

    logger.info("All done.")


def merge_peaks_parser():
    parser = argparse.ArgumentParser(
        prog="kaic merge_peaks",
        description='Filter peaks called with call_peaks'
    )

    parser.add_argument(
        'input',
        help='''Input Peaks file'''
    )

    parser.add_argument(
        'output',
        help='''Output merged Peaks file'''
    )

    parser.add_argument(
        '-d', '--distance', dest='distance',
        type=int,
        default=20000,
        help='''Maximum distance in base pairs at which to merge two peaks. Default 20000bp'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)
    return parser


def merge_peaks(argv):
    parser = merge_peaks_parser()

    args = parser.parse_args(argv[2:])

    import kaic.data.network as kn
    from kaic.tools.files import create_temporary_copy

    # copy file if required
    original_input_path = os.path.expanduser(args.input)
    original_output_path = os.path.expanduser(args.output)
    if args.tmp:
        logger.info("Copying data to temporary file...")
        input_path = create_temporary_copy(original_input_path)
        logger.info("Working from temporary file %s" % input_path)
        tmp_file = tempfile.NamedTemporaryFile(delete=False)
        tmp_file.close()
        output_path = tmp_file.name
        logger.info("Temporary output file: %s" % output_path)
    else:
        input_path = os.path.expanduser(args.input)
        output_path = os.path.expanduser(args.output)

    peaks = kn.RaoPeakInfo(input_path, mode='r')
    peaks.merged_peaks(output_path, euclidian_distance=args.distance)

    if args.tmp:
        os.unlink(input_path)
        logger.info("Moving temporary output file to destination %s" % original_output_path)
        shutil.move(output_path, original_output_path)


def filter_merged_peaks_parser():
    parser = argparse.ArgumentParser(
        prog="kaic filter_merged_peaks",
        description='Filter merged peaks'
    )

    parser.add_argument(
        'input',
        help='''Input merged Peaks file'''
    )

    parser.add_argument(
        'output',
        nargs='?',
        help='''Output filtered merged Peaks file'''
    )

    parser.add_argument(
        '-r', '--rao', dest='rao',
        action='store_true',
        help='''Filter peaks as Rao et al. (2014) does.
                It removes peaks that are singlets and have a q-value sum >.02.
            '''
    )
    parser.set_defaults(rao=False)

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)
    return parser


def filter_merged_peaks(argv):
    parser = filter_merged_peaks_parser()

    args = parser.parse_args(argv[2:])

    import kaic.data.network as kn
    from kaic.tools.files import create_temporary_copy, copy_or_expand

    # copy file if required
    original_input_path = os.path.expanduser(args.input)
    if args.tmp:
        logger.info("Copying data to temporary file...")
        input_path = create_temporary_copy(original_input_path)
    else:
        input_path = copy_or_expand(args.input, args.output)

    merged_peaks = kn.PeakInfo(input_path, mode='a')

    if args.rao:
        logger.info("Running Rao filter")
        merged_peaks.filter_rao()

    if args.tmp:
        output_path = os.path.expanduser(args.output)
        if os.path.isdir(output_path):
            output_path = "%s/%s" % (output_path, os.path.basename(original_input_path))
        logger.info("Moving temporary output file to destination %s..." % output_path)
        shutil.move(input_path, output_path)

    logger.info("All done.")


def overlap_peaks_parser():
    parser = argparse.ArgumentParser(
        prog="kaic overlap_peaks",
        description='Overlap peaks from multiple samples'
    )

    parser.add_argument(
        'input',
        nargs='+',
        help='''Input Peak files. Two or more.'''
    )

    parser.add_argument(
        'output',
        help='''Output directory. Overlapped peaks and stats are written there.'''
    )

    parser.add_argument(
        '-d', '--distance', dest='distance',
        type=int,
        help='''Maximum distance between peaks for merging them. Default=3x bin size'''
    )
    parser.set_defaults(distance=None)

    parser.add_argument(
        '-n', '--names', dest='names',
        nargs='*',
        help='''Names for input Peak samples. Default: Use file names'''
    )
    parser.set_defaults(names=None)

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)
    return parser


def overlap_peaks(argv):
    parser = overlap_peaks_parser()

    args = parser.parse_args(argv[2:])

    import kaic.data.network as kn
    from kaic.tools.files import create_temporary_copy

    original_input_paths = [os.path.expanduser(i) for i in args.input]
    if not args.names:
        names = [os.path.splitext(os.path.basename(i))[0] for i in original_input_paths]
    else:
        names = args.names

    if len(original_input_paths) <= 2:
        raise ValueError("Need 2 or more inputs.")
    if len(names) != len(original_input_paths):
        raise ValueError("Number of inputs and names is different.")
    if len(set(names)) < len(names):
        raise ValueError("Names not unique.")

    # copy file if required
    if args.tmp:
        logger.info("Copying data to temporary file...")
        input_paths = [create_temporary_copy(i) for i in original_input_paths]
    else:
        input_paths = original_input_paths

    peaks = [kaic.load(i, mode="r") for i in input_paths]

    if not args.distance:
        distance = peaks[0].bin_size*3
    else:
        distance = args.distance

    stats, merged = kn.overlap_peaks({n:p for n, p in zip(names, peaks)}, max_distance=distance)

    if args.tmp:
        for i in input_paths:
            os.remove(i)

    output_path = os.path.expanduser(args.output)
    logger.info("Writing files...")
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    for name_set, peaks in merged.items():
        peaks.file.copy_file(os.path.join(output_path, "_".join(n for n in names if n in name_set) + ".peaks"), overwrite=True)
    stats.to_csv(os.path.join(output_path, "stats_" + "_".join(names) + ".tsv"), sep="\t", index=False)

    logger.info("All done.")


def plot_ligation_err_parser():
    parser = argparse.ArgumentParser(
        prog="kaic plot_ligation_err",
        description='Plot the ligation structure biases of a Pairs object'
    )

    parser.add_argument(
        'input',
        help='''Input Pairs file'''
    )

    parser.add_argument(
        'output',
        nargs='?',
        help='''Output pdf'''
    )

    parser.add_argument(
        '-p', '--points', dest='points',
        type=int,
        help='''Data points that make up one increment of the x axis. More=smoother=less detail.'''
    )

    return parser


def plot_ligation_err(argv):
    parser = plot_ligation_err_parser()
    args = parser.parse_args(argv[2:])

    import kaic
    from kaic.plotting.plot_statistics import hic_ligation_structure_biases_plot

    input_path = os.path.expanduser(args.input)
    output_path = None
    if args.output:
        output_path = os.path.expanduser(args.output)

    pairs = kaic.load(file_name=input_path, mode='r')
    hic_ligation_structure_biases_plot(pairs, output=output_path, sampling=args.points)
    pairs.close()

    logger.info("All done.")


def plot_re_dist_parser():
    parser = argparse.ArgumentParser(
        prog="kaic plot_re_dist",
        description='Plot the restriction site distance of reads in a Pairs object'
    )

    parser.add_argument(
        'input',
        help='''Input Pairs file'''
    )

    parser.add_argument(
        'output',
        nargs='?',
        help='''Output pdf'''
    )

    parser.add_argument(
        '-l', '--limit', dest='limit',
        type=int,
        default=10000,
        help='''Limit the plot to the first LIMIT read pairs for the sake of speed. Default 10000'''
    )

    parser.add_argument(
        '-m', '--max-dist', dest='max_dist',
        type=int,
        help='''Maximum RE site distance to include in the plot. Default: no max'''
    )

    return parser


def plot_re_dist(argv):
    parser = plot_re_dist_parser()
    args = parser.parse_args(argv[2:])

    import kaic
    from kaic.plotting.plot_statistics import pairs_re_distance_plot

    input_path = os.path.expanduser(args.input)
    output_path = None
    if args.output:
        output_path = os.path.expanduser(args.output)

    pairs = kaic.load(file_name=input_path, mode='r')
    pairs_re_distance_plot(pairs, output=output_path, limit=args.limit, max_distance=args.max_dist)
    pairs.close()

    logger.info("All done.")


def plot_hic_corr_parser():
    parser = argparse.ArgumentParser(
        prog="kaic plot_hic_corr",
        description='Plot the correlation of two Hic objects'
    )

    parser.add_argument(
        'hic1',
        help='''First Hi-C file'''
    )

    parser.add_argument(
        'hic2',
        help='''Second Hi-C file'''
    )

    parser.add_argument(
        'output',
        nargs="?",
        help='''Output PDF file'''
    )

    parser.add_argument(
        '-c', '--colormap', dest='colormap',
        help='''Matplotlib colormap'''
    )
    return parser


def plot_hic_corr(argv):
    parser = plot_hic_corr_parser()

    args = parser.parse_args(argv[2:])

    import kaic
    from kaic.config import config
    from kaic.plotting.plot_genomic_data import hic_correlation_plot

    colormap = config.colormap_hic if args.colormap is None else args.colormap

    hic1_path = os.path.expanduser(args.hic1)
    hic2_path = os.path.expanduser(args.hic2)

    hic1 = kaic.load_hic(hic1_path, mode='r')
    hic2 = kaic.load_hic(hic2_path, mode='r')

    output_path = None
    if args.output:
        output_path = os.path.expanduser(args.output)

    hic_correlation_plot(hic1, hic2, output=output_path, colormap=colormap, size=15)

    hic1.close()
    hic2.close()
    logger.info("All done.")


def plot_hic_marginals_parser():
    parser = argparse.ArgumentParser(
        prog="kaic plot_hic_marginals",
        description='Plot Hic matrix marginals'
    )

    parser.add_argument(
        'input',
        help='''Input Hi-C file'''
    )

    parser.add_argument(
        'output',
        nargs="?",
        help='''Output PDF file'''
    )

    parser.add_argument(
        '-l', '--lower', dest='lower',
        type=float,
        help='''Plot lower coverage bound at this level'''
    )

    parser.add_argument(
        '-u', '--upper', dest='upper',
        type=float,
        help='''Plot lower coverage bound at this level'''
    )
    return parser


def plot_hic_marginals(argv):
    parser = plot_hic_marginals_parser()
    args = parser.parse_args(argv[2:])

    import kaic
    from kaic.plotting.plot_genomic_data import hic_marginals_plot

    input_path = os.path.expanduser(args.input)

    hic = kaic.load_hic(input_path, mode='r')

    output_path = None
    if args.output:
        output_path = os.path.expanduser(args.output)

    hic_marginals_plot(hic, output=output_path, lower=args.lower, upper=args.upper)
    hic.close()
    logger.info("All done.")


def structure_tracks_parser():
    parser = argparse.ArgumentParser(
        prog="kaic structure_tracks",
        description='Calculate genomic tracks about structural features of the Hi-C map'
    )
    parser.add_argument(
        'hic',
        help='Input Hic file'
    )

    parser.add_argument(
        'output',
        help='Output path for genomic track'
    )

    parser.add_argument(
        'window_size',
        type=int,
        nargs='+',
        help='Window sizes (in base pairs) used for directionality index,'
             'insulation index and relative insulation index.'
    )

    parser.add_argument(
        '-oe', '--observed_expected',
        action='store_true',
        help='Use observed over expected heatmap for calculations.'
    )

    parser.add_argument(
        '--no_imputation',
        action="store_true",
        help='Do not use imputation to guess value of unmappable bins. '
             'Turning off imputation may lead to artifacts '
             'near unmappable bins. The mask threshold should '
             'be set to a very low value (.1 or so) in this case.'
    )

    parser.add_argument(
        '-n', '--normalise', dest='normalise',
        action='store_true',
        help='''Normalise index values'''
    )
    parser.set_defaults(normalise=False)

    parser.add_argument(
        '-nw', '--normalisation_window',
        type=int, default=300,
        help='Window size for calculating long-range mean for normalization of insulation_index,'
             ' relative_insulation_index, contact_band.'
             'Default 300 bins.'
    )

    parser.add_argument(
        '-o', '--offset',
        type=int, default=0,
        help='Offset of insulation index window from the diagonal in base pairs.'
    )

    parser.add_argument(
        '-w', '--smoothing_window',
        type=int, default=15,
        help='Window size for smoothing derivatives in Savitzky Golay filter (in bins).'
             'Default 15. Must be an odd number.'
    )

    parser.add_argument(
        '-p', '--poly_order',
        type=int, default=2,
        help='Order of polynomial used for smoothing derivatives. Default 2.'
    )

    parser.add_argument(
        '-d', '--derivative',
        type=int,
        nargs='+',
        help='Optionally save derivatives of the specified order (>1).'
    )

    parser.add_argument(
        '--delta',
        type=int,
        nargs='+',
        help='Save delta transformation of metrics according to Crane et al. 2015. '
             'Specify window size in bins. Sensible values for 5kb Drosophila 5-10.'
    )

    parser.add_argument(
        '-ii', '--insulation_index',
        action='store_true',
        help='Calculate insulation index for the given distances (in bins).'
    )

    parser.add_argument(
        '-di', '--directionality_index',
        action='store_true',
        help='Calculate the directionality index for the given distances (in bp)'
    )

    parser.add_argument(
        '-r', '--relative',
        action='store_true',
        help='Calculate the relative insulation indices for the given distances (in bins)'
    )
    return parser


def structure_tracks(argv):
    parser = structure_tracks_parser()

    args = parser.parse_args(argv[2:])

    import kaic
    import numpy as np
    from kaic.architecture.hic_architecture import InsulationIndex, DirectionalityIndex, ObservedExpectedRatio
    from kaic.architecture.genome_architecture import GenomicTrack
    from kaic.tools.matrix import delta_window
    from scipy.signal import savgol_filter

    input_path = os.path.expanduser(args.hic)
    output_path = os.path.expanduser(args.output)

    logger.info("Fetching Hi-C matrix")
    hic_matrix = kaic.load(input_path)

    if args.observed_expected:
        logger.info("Calculating observed/expected ratio")
        hic_matrix = ObservedExpectedRatio(hic_matrix)

    try:
        os.remove(output_path)
    except OSError:
        pass

    # prepare genomic track object
    gt = GenomicTrack(output_path, regions=hic_matrix.regions)

    # calculate insulation index
    if args.insulation_index:
        with InsulationIndex(hic_matrix, relative=args.relative, offset=args.offset,
                             normalise=args.normalise, window_sizes=args.window_size,
                             _normalisation_window=args.normalisation_window) as ii:
            for window_size in args.window_size:
                insulation_index = ii.insulation_index(window_size)
                gt.add_data("insulation_index_{}".format(window_size), insulation_index)

    # calculate directioality index
    if args.directionality_index:
        with DirectionalityIndex(hic_matrix, window_sizes=args.window_size) as di:
            for window_size in args.window_size:
                directionality_index = di.directionality_index(window_size)
                gt.add_data("directionality_index_{}".format(window_size), directionality_index)

    # calculate derivatives, if requested
    if args.derivative or args.delta:
        for k, v in gt.tracks.items():
            if args.derivative:
                for i in args.derivative:
                    if "matrix" in k:
                        deriv_matrix = np.vstack([savgol_filter(x, window_length=args.smoothing_window,
                                                                polyorder=args.poly_order, deriv=i) for x in v.T]).T
                        gt.add_data("{}_d{}".format(k, i), deriv_matrix)
                    else:
                        d = savgol_filter(v, window_length=args.smoothing_window,
                                          polyorder=args.poly_order, deriv=i)
                        gt.add_data("{}_d{}".format(k, i), d)
            if args.delta:
                for i in args.delta:
                    if "matrix" in k:
                        delta_matrix = np.vstack([delta_window(x, i) for x in v.T]).T
                        gt.add_data("{}_delta{}".format(k, i), delta_matrix)
                    else:
                        gt.add_data("{}_delta{}".format(k, i), delta_window(v, i))
    logger.info("All done.")


def boundaries_parser():
    parser = argparse.ArgumentParser(
        prog="kaic boundaries",
        description='Determine structural boundaries'
    )
    parser.add_argument(
        'architecture',
        help='Input InsulationIndex file'
    )
    parser.add_argument(
        'output',
        help="Output folder for boundary BED files (default or when using '-r' option) or "
             "path for boundary BED file (when using -w option)."
    )
    parser.add_argument(
        '-r', '--range', dest='range',
        type=int,
        nargs=2,
        help='Range of insulation index window sizes (<low> <high>) to calculate boundaries on.'
    )
    parser.add_argument(
        '-w', '--window', dest='window',
        type=int,
        help='Insulation index window size to calculate boundaries on'
    )
    parser.add_argument(
        '-d', '--delta', dest='delta',
        type=int, default=3,
        help='Window size for calculating the delta vector (in bins). Calculation takes into '
             'account d bins upstream and d bins downstream for a total '
             'window size of 2*d + 1 bins. Default 3.'
    )
    parser.add_argument(
        '-s', '--min-score', dest='min_score',
        type=float,
        help='Report only peaks where the two surrounding extrema of the delta vector have '
             'at least this difference in height. Default: no threshold.'
    )
    parser.add_argument(
        '-x', '--sub-bin-precision', dest='sub_bin_precision',
        action='store_true',
        help='Report boundary positions with sub-bin precision. This works because the minimum '
             'the insulation score can be determined with sub-bin precision. Default: False'
    )
    parser.add_argument(
        '-p', '--prefix', dest='prefix',
        default='boundaries',
        help='''Output file prefix. Not necessary when using 'w' modus. Default: boundaries'''
    )
    parser.add_argument(
        '-l', '--log', dest='log',
        action='store_true',
        help='''log-transform index values before boundary calling.'''
    )
    parser.set_defaults(log=False)
    return parser


def boundaries(argv):
    parser = boundaries_parser()

    args = parser.parse_args(argv[2:])

    import kaic
    from kaic.tools.general import mkdir

    input_file = os.path.expanduser(args.architecture)
    output_path = os.path.expanduser(args.output)

    array = kaic.load(input_file, mode='r')

    single = False
    window_sizes = []
    if args.range is not None:
        for window_size in array.window_sizes:
            if args.range[0] <= window_size <= args.range[1]:
                window_sizes.append(window_size)
    elif args.window is not None:
        if args.window in array.window_sizes:
            window_sizes.append(args.window)
        single = True
    else:
        window_sizes = array.window_sizes

    if len(window_sizes) == 0:
        raise ValueError("No valid window size specified!")

    def _to_bed(bs, file_name):
        with open(file_name, 'w') as bed:
            for b in bs:
                bed.write("{}\t{}\t{}\t.\t{}\n".format(b.chromosome, b.start, b.end, b.score))

    if not single:
        mkdir(output_path)
        for window_size in window_sizes:
            logger.info("Processing window size: {}".format(window_size))
            boundaries = array.boundaries(window_size, min_score=args.min_score,
                                          delta_window=args.delta, log=args.log,
                                          sub_bin_precision=args.sub_bin_precision)
            _to_bed(boundaries, output_path + "/{}_{}.bed".format(args.prefix, window_size))
    else:
        boundaries = array.boundaries(window_sizes[0], min_score=args.min_score,
                                      delta_window=args.delta, log=args.log,
                                      sub_bin_precision=args.sub_bin_precision)
        _to_bed(boundaries, output_path)

    logger.info("All done.")


def fold_change_parser():
    parser = argparse.ArgumentParser(
        prog="kaic fold_change",
        description='Create pairwise fold-change Hi-C comparison maps'
    )
    parser.add_argument(
        'input',
        nargs=2,
        help='Input Hic files'
    )
    parser.add_argument(
        'output',
        help='Output FoldChangeMatrix file.'
    )

    parser.add_argument(
        '-S', '--no-scale', dest='scale',
        action='store_false',
        help='''Do not scale input matrices'''
    )
    parser.set_defaults(scale=True)

    parser.add_argument(
        '-l', '--log2', dest='log',
        action='store_true',
        help='''Log2-convert fold-change values'''
    )
    parser.set_defaults(log=False)

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)
    return parser


def fold_change(argv):
    parser = fold_change_parser()

    args = parser.parse_args(argv[2:])

    import kaic
    from kaic.architecture.hic_architecture import FoldChangeMatrix
    import os.path

    tmpdir = None
    if args.tmp:
        import tempfile
        tmpdir = tempfile.gettempdir()

    hic1 = kaic.load_hic(os.path.expanduser(args.input[0]), mode='r')
    hic2 = kaic.load_hic(os.path.expanduser(args.input[1]), mode='r')

    output_file = os.path.expanduser(args.output)
    with FoldChangeMatrix(hic1, hic2, file_name=output_file, tmpdir=tmpdir, mode='w',
                          scale_matrices=args.scale, log2=args.log) as fcm:
        fcm.calculate()


def average_tracks_parser():
    parser = argparse.ArgumentParser(
        prog="kaic average_tracks",
        description='Calculate average Hi-C contact profiles per region'
    )
    parser.add_argument(
        'input',
        help='Input matrix (Hi-C, fold-change map, ...)'
    )
    parser.add_argument(
        'output',
        help='Output RegionContactAverage file.'
    )

    parser.add_argument(
        '-w', '--window-sizes', dest='window_sizes',
        nargs='+',
        type=int,
        default=[200000, 400000, 600000, 1000000],
        help='''Window sizes in base pairs to calculate region average in.
                    The total window size is composed of the left window plus the right window, i.e. 2x this value.'''
    )

    parser.add_argument(
        '-o', '--offset', dest='offset',
        type=int,
        default=0,
        help='''Window offset in base pairs from the diagonal.'''
    )

    parser.add_argument(
        '-p', '--padding', dest='padding',
        type=int,
        default=1,
        help='''Padding (in number of regions) to calculate average on larger regions.
                    Acts similarly to curve smooting'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)

    parser.add_argument(
        '-i', '--impute', dest='impute',
        action='store_true',
        help='''Impute missing values in matrix'''
    )
    parser.set_defaults(impute=False)
    return parser


def average_tracks(argv):
    parser = average_tracks_parser()

    args = parser.parse_args(argv[2:])

    import kaic
    from kaic.architecture.hic_architecture import RegionContactAverage
    import os.path
    import tempfile

    input_file = os.path.expanduser(args.input)
    output_file = os.path.expanduser(args.output)
    tmpdir = None
    if args.tmp:
        tmpdir = tempfile.gettempdir()

    matrix = kaic.load(input_file, mode='r')

    with RegionContactAverage(matrix, file_name=output_file, tmpdir=tmpdir, window_sizes=args.window_sizes,
                              offset=args.offset, padding=args.padding, impute_missing=args.impute) as rca:
        rca.calculate()


def directionality_parser():
    parser = argparse.ArgumentParser(
        prog="kaic directionality",
        description='Calculate directionality index for Hic object'
    )
    parser.add_argument(
        'input',
        help='Input matrix (Hi-C, fold-change map, ...)'
    )
    parser.add_argument(
        'output',
        help='Output DirectionalityIndex file.'
    )

    parser.add_argument(
        '-w', '--window-sizes', dest='window_sizes',
        nargs='+',
        type=int,
        default=[200000, 400000, 600000, 1000000],
        help='''Window sizes in base pairs to calculate directionality index on.
                    The total window size is composed of the left window plus the right window, i.e. 2x this value.'''
    )

    parser.add_argument(
        '-r', '--region', dest='region',
        help='''Region selector (<chr>:<start>-<end>) to only calculate DI for this region.'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)

    parser.add_argument(
        '-i', '--impute', dest='impute',
        action='store_true',
        help='''Impute missing values in matrix'''
    )
    parser.set_defaults(impute=False)
    return parser


def directionality(argv):
    parser = directionality_parser()

    args = parser.parse_args(argv[2:])

    import kaic
    from kaic.architecture.hic_architecture import DirectionalityIndex
    import os.path
    import tempfile

    input_file = os.path.expanduser(args.input)
    output_file = os.path.expanduser(args.output)
    tmpdir = None
    if args.tmp:
        tmpdir = tempfile.gettempdir()

    matrix = kaic.load(input_file, mode='r')

    with DirectionalityIndex(matrix, file_name=output_file, tmpdir=tmpdir, window_sizes=args.window_sizes) as di:
        di.calculate()


def insulation_parser():
    parser = argparse.ArgumentParser(
        prog="kaic insulation",
        description='Calculate insulation index for Hic object'
    )
    parser.add_argument(
        'input',
        help='Input matrix (Hi-C, fold-change map, ...)'
    )
    parser.add_argument(
        'output',
        help='Output InsulationIndex file.'
    )

    parser.add_argument(
        '-w', '--window-sizes', dest='window_sizes',
        nargs='+',
        type=int,
        default=[200000, 400000, 600000, 1000000],
        help='''Window sizes in base pairs to calculate insulation index on.
                    The total window size is composed of the left window plus the right window, i.e. 2x this value.'''
    )

    parser.add_argument(
        '-r', '--region', dest='region',
        help='''Region selector (<chr>:<start>-<end>) to only calculate II for this region.'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)

    parser.add_argument(
        '-i', '--impute', dest='impute',
        action='store_true',
        help='''Impute missing values in matrix'''
    )
    parser.set_defaults(impute=False)

    parser.add_argument(
        '-o', '--offset', dest='offset',
        type=int,
        default=0,
        help='''Window offset in base pairs from the diagonal.'''
    )

    parser.add_argument(
        '-l', '--log', dest='log',
        action='store_true',
        help='''Log2-transform insulation index after normalisation (roughly centers values around 0).'''
    )
    parser.set_defaults(log=False)

    parser.add_argument(
        '-n', '--normalise', dest='normalise',
        action='store_true',
        help='''Normalise index to insulation average (default is per-chromosome - to normalise to
                    smaller regions, use -nw).'''
    )
    parser.set_defaults(normalise=False)

    parser.add_argument(
        '-nw', '--normalisation-window', dest='normalisation_window',
        type=int,
        help='''Size of the normalisation window (moving average) in bins. Default: whole chromosome.'''
    )

    parser.add_argument(
        '-s', '--subtract-mean', dest='subtract',
        action='store_true',
        help='''Subtract mean instead of dividing by it when '--normalise' is enabled.
                You probably don't want this, unless you are working with 
                log-transformed matrices (e.g. fold-change matrices)'''
    )
    parser.set_defaults(subtract=False)
    return parser


def insulation(argv):
    parser = insulation_parser()

    args = parser.parse_args(argv[2:])

    import kaic
    from kaic.architecture.hic_architecture import InsulationIndex
    import os.path
    import tempfile

    input_file = os.path.expanduser(args.input)
    output_file = os.path.expanduser(args.output)
    tmpdir = None
    if args.tmp:
        tmpdir = tempfile.gettempdir()

    matrix = kaic.load(input_file, mode='r')

    with InsulationIndex(matrix, file_name=output_file, tmpdir=tmpdir, window_sizes=args.window_sizes,
                         impute_missing=args.impute, normalise=args.normalise, offset=args.offset,
                         mode='w', subtract_mean=args.subtract, log=args.log,
                         _normalisation_window=args.normalisation_window) as ii:
        ii.calculate()


def ii_to_bw_parser():
    parser = argparse.ArgumentParser(
        prog="kaic ii_to_bw",
        description='Convert insulation index object to BigWig'
    )
    parser.add_argument(
        'input',
        help='InsulationIndex file'
    )
    parser.add_argument(
        'output',
        help='Folder for BigWig output or BigWig filename (if -w option in use and only a single size'
             'was specified).'
    )

    parser.add_argument(
        '-w', '--window-sizes', dest='window_sizes',
        nargs='+',
        type=int,
        help='''Window sizes in base pairs to convert to BigWig.'''
    )

    parser.add_argument(
        '-r', '--region', dest='region',
        help='''Region selector (<chr>:<start>-<end>) to only write values for these regions.'''
    )

    parser.add_argument(
        '-p', '--prefix', dest='prefix',
        help='''Prefix for files when using multiple window sizes.'''
    )

    return parser


def ii_to_bw(argv):
    parser = ii_to_bw_parser()

    args = parser.parse_args(argv[2:])

    import os
    import kaic

    input_file = os.path.expanduser(args.input)
    output = os.path.expanduser(args.output)
    window_sizes = args.window_sizes
    region = args.region
    prefix = args.prefix

    if prefix is None:
        prefix = os.path.basename(os.path.splitext(input_file)[0])

    ii = kaic.load(input_file, mode='r')

    if window_sizes is None:
        window_sizes = ii.window_sizes

    for window_size in window_sizes:
        logger.info("Window size {}".format(window_size))
        output_file = os.path.join(output, prefix + '_{}.bw'.format(window_size)) if len(window_sizes) > 1 else output
        ii.to_bigwig(output_file, subset=region, score_field='ii_{}'.format(window_size))
    ii.close()


def ab_parser():
    parser = argparse.ArgumentParser(
        prog="kaic ab",
        description='Calculate AB compartment matrix'
    )
    parser.add_argument(
        'input',
        help='Input matrix (Hi-C, fold-change map, ...)'
    )
    parser.add_argument(
        'output',
        help='Output AB matrix file.'
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)

    return parser


def ab(argv):
    parser = ab_parser()

    args = parser.parse_args(argv[2:])

    import kaic
    from kaic.architecture.hic_architecture import ABDomainMatrix
    import os.path
    import tempfile

    input_file = os.path.expanduser(args.input)
    output_file = os.path.expanduser(args.output)
    tmpdir = None
    if args.tmp:
        tmpdir = tempfile.gettempdir()

    matrix = kaic.load(input_file, mode='r')

    with ABDomainMatrix(matrix, file_name=output_file, mode='w', tmpdir=tmpdir) as ab:
        ab.calculate()


def ab_domains_parser():
    parser = argparse.ArgumentParser(
        prog="kaic ab_domains",
        description='Assign A or B compartment to each region in Hi-C object'
    )
    parser.add_argument(
        'input',
        help='Input matrix (Hi-C or ABDomains)'
    )
    parser.add_argument(
        'output',
        help='Output BED.'
    )

    parser.add_argument(
        '-g', '--genome', dest='genome',
        help='''Genome file (FASTA, folder with FASTA, comma-separated list of FASTA, Genome object)
                used to change sign of eigenvector based on GC content.'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)

    return parser


def ab_domains(argv):
    parser = ab_domains_parser()

    args = parser.parse_args(argv[2:])
    import os
    input_file = os.path.expanduser(args.input)
    output_file = os.path.expanduser(args.output)
    genome_file = args.genome
    tmp = args.tmp

    import kaic
    from kaic.tools.files import write_bed

    matrix = kaic.load(input_file, mode='r')
    if isinstance(matrix, kaic.Hic):
        matrix = kaic.ABDomainMatrix(matrix, tmpdir=tmp)
        matrix.calculate()
    elif not isinstance(matrix, kaic.ABDomainMatrix):
        raise ValueError("Supplied input must be either a Hic object or an ABDomainMatrix")

    with kaic.ABDomains(matrix, genome=genome_file, tmpdir=tmp) as abd:
        regions = abd.ab_regions()
        write_bed(output_file, regions)
    matrix.close()


def distance_decay_parser():
    parser = argparse.ArgumentParser(
        prog="kaic distance_decay",
        description='Calculate Hi-C distance decay (expected values)'
    )
    parser.add_argument(
        'input',
        help='Input matrix (Hi-C, fold-change map, ...)'
    )
    parser.add_argument(
        'output',
        help='Output distance decay file.'
    )

    parser.add_argument(
        '-r', '--regions', dest='regions',
        help='''Region subset for expected value calculation.'''
    )

    parser.add_argument(
        '-s', '--smooth', dest='smooth',
        action='store_true',
        help='''Smoothe expected value curve'''
    )
    parser.set_defaults(smooth=False)

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)

    return parser


def distance_decay(argv):
    parser = distance_decay_parser()

    args = parser.parse_args(argv[2:])

    import kaic
    from kaic.architecture.hic_architecture import ExpectedContacts
    import os.path
    import tempfile

    input_file = os.path.expanduser(args.input)
    output_file = os.path.expanduser(args.output)
    tmpdir = None
    if args.tmp:
        tmpdir = tempfile.gettempdir()

    matrix = kaic.load(input_file, mode='r')

    with ExpectedContacts(matrix, file_name=output_file, smooth=args.smooth, mode='w',
                          regions=args.regions, tmpdir=tmpdir) as ex:
        ex.calculate()


def optimise_parser():
    parser = argparse.ArgumentParser(
        prog="kaic optimise",
        description='Optimise a Hic object for faster access'
    )
    parser.add_argument(
        'input',
        help='Input Hic file'
    )
    parser.add_argument(
        'output',
        help='Output AccessOptimisedHic file.'
    )
    return parser


def optimise(argv):
    parser = optimise_parser()
    args = parser.parse_args(argv[2:])

    import kaic
    import os.path
    old_hic = kaic.load_hic(os.path.expanduser(args.input), mode='r')
    new_hic = kaic.AccessOptimisedHic(old_hic, file_name=os.path.expanduser(args.output))
    new_hic.close()
    old_hic.close()


def subset_hic_parser():
    parser = argparse.ArgumentParser(
        prog="kaic subset",
        description='Create a new Hic object by subsetting'
    )
    parser.add_argument(
        'input',
        help='Input Hic file'
    )
    parser.add_argument(
        'output',
        help='Output Hic file.'
    )

    parser.add_argument(
        'regions',
        nargs='+'
    )
    return parser


def subset_hic(argv):
    parser = subset_hic_parser()
    args = parser.parse_args(argv[2:])

    import os.path
    import kaic

    input_file = os.path.expanduser(args.input)
    output_file = os.path.expanduser(args.output)

    old_hic = kaic.load_hic(input_file, mode='r')
    new_hic = kaic.AccessOptimisedHic(file_name=output_file, mode='w')

    bias_vector = old_hic.bias_vector()

    new_bias_vector = []
    ix_converter = {}
    ix = 0
    for region_string in args.regions:
        for region in old_hic.subset(region_string):
            ix_converter[region.ix] = ix
            ix += 1

            new_hic.add_region(region, flush=False)
            new_bias_vector.append(bias_vector[region.ix])
    new_hic.flush()

    for region_string in args.regions:
        for edge in old_hic.edge_subset(key=(region_string, region_string), lazy=True):
            source = ix_converter[edge.source]
            sink = ix_converter[edge.sink]
            new_hic.add_edge([source, sink, edge.weight], flush=False)
    new_hic.flush()

    new_hic.bias_vector(new_bias_vector)


def diff_parser():
    parser = argparse.ArgumentParser(
        prog="kaic diff",
        description='Calculate difference between two vectors (v1-v2)'
    )

    parser.add_argument(
        'vector1',
        help='First vector (/array, e.g. InsulationIndex)'
    )

    parser.add_argument(
        'vector2',
        help='Second vector (/array, e.g. InsulationIndex)'
    )

    parser.add_argument(
        'output',
        help='Output VectorDifference file.'
    )

    parser.add_argument(
        '-a', '--absolute', dest='absolute',
        action='store_true',
        help='''Output absolute difference'''
    )
    parser.set_defaults(absolute=False)
    return parser


def diff(argv):
    parser = diff_parser()

    args = parser.parse_args(argv[2:])

    import os.path
    import kaic

    v1 = kaic.load(args.vector1, mode='r')
    v2 = kaic.load(args.vector2, mode='r')

    output_file = os.path.expanduser(args.output)

    with kaic.VectorDifference(v1, v2, absolute=args.absolute, file_name=output_file, mode='w') as d:
        d.calculate()


def aggregate_tads_parser():
    parser = argparse.ArgumentParser(
        prog="kaic aggregate_tads",
        description='Make a TAD aggregate plot'
    )

    parser.add_argument(
        'hic',
        help='Hic file'
    )

    parser.add_argument(
        'tads',
        help='File with TAD regions (BED, GFF, Tabix, ...)'
    )

    parser.add_argument(
        'output',
        help='Output image file (extension determines file format)'
    )

    parser.add_argument(
        '-p', '--pixels', dest='pixels',
        type=int,
        default=90,
        help='''Width of the output image in pixels. Default: 90'''
    )

    parser.add_argument(
        '-i', '--interpolation', dest='interpolation',
        default='nearest',
        help='''Type of interpolation performed to shrink/expand matrices. Default: nearest'''
    )

    parser.add_argument(
        '-r', '--relative', dest='relative',
        type=float,
        default=1.0,
        help='''Extension (e) of each region as fraction of TAD length (l). 
                    Final region in the image will be: 
                    <start of TAD - e*l> to <end of TAD + e*l>.
                    Default: 1.0 (results in 3 times TAD size image). 
                    Additive with '-a' parameter!'''
    )

    parser.add_argument(
        '-a', '--absolute', dest='absolute',
        type=int,
        default=0,
        help='''Extension (e) of each region in base pairs. 
                    Final region in the image will be: 
                    <start of TAD - e> to <end of TAD + e>.
                    Default: 0 (no extension). 
                    Additive with '-r' parameter!'''
    )

    parser.add_argument(
        '-n', '--norm', dest='norm',
        action='store_true',
        help='''Normalize matrix to expected values'''
    )
    parser.set_defaults(norm=False)

    parser.add_argument(
        '-L', '--no-log', dest='log',
        action='store_false',
        help='''Do not log2-transform normalized matrices.
                Only used in conjunction with '--norm'.'''
    )
    parser.set_defaults(log=True)

    parser.add_argument(
        '-s', '--rescale', dest='rescale',
        action='store_true',
        help='''Do not rescale normalized contact matrices using an a=-0.25 power law.
                Only used in conjunction with '--norm'.'''
    )
    parser.set_defaults(rescale=False)

    parser.add_argument(
        '-C', '--no-cache', dest='cache',
        action='store_false',
        help='''Do not cache chromosome matrices (slower, but saves a lot of memory)'''
    )
    parser.set_defaults(cache=True)

    parser.add_argument(
        '-m', '--save-matrix', dest='matrix_file',
        help='''Path to save aggregate matrix (numpy txt format)'''
    )

    parser.add_argument(
        '-c', '--colormap', dest='colormap',
        help='''Matplotlib colormap to use for matrix'''
    )

    parser.add_argument(
        '--vmin', dest='vmin',
        type=float,
        help='''Minimum saturation value in image'''
    )

    parser.add_argument(
        '--vmax', dest='vmax',
        type=float,
        help='''Maximum saturation value in image'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)

    parser.add_argument(
        '--preset-flyamer', dest='preset_flyamer',
        action='store_true',
        help='''Use presets to create matrices as in Flyamer et al. (2017).
                Overrides the following parameters:
                --norm -L -s -r 1.0 -a 0'''
    )
    parser.set_defaults(preset_flyamer=False)

    return parser


def aggregate_tads(argv):
    parser = aggregate_tads_parser()

    args = parser.parse_args(argv[2:])
    import os

    hic_file = os.path.expanduser(args.hic)
    tads_file = os.path.expanduser(args.tads)
    output_file = os.path.expanduser(args.output)
    pixels = args.pixels
    interpolation = args.interpolation
    relative = args.relative
    absolute = args.absolute
    norm = args.norm
    log = args.log
    rescale = args.rescale
    cache = args.cache
    matrix_file = None if args.matrix_file is None else os.path.expanduser(args.matrix_file)
    cmap = args.colormap
    vmin = args.vmin
    vmax = args.vmax
    tmp = args.tmp
    preset_flyamer = args.preset_flyamer

    if preset_flyamer:
        logging.info("Using Flyamer et al. (2017) preset")
        norm = True
        rescale = True
        relative = 1.0
        absolute = 0
        log = False

    import kaic
    from kaic.config import config
    from kaic.architecture.hic_architecture import aggregate_tads
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import kaic.plotting

    if cmap is None:
        cmap = 'bwr' if norm and log else config.colormap_hic

    with kaic.load(hic_file, mode='r', tmpdir=tmp) as hic:
        tads = kaic.load(tads_file)
        m = aggregate_tads(hic, tads.regions, pixels=pixels, interpolation=interpolation,
                           relative_extension=relative, absolute_extension=absolute,
                           rescale=rescale, cache=cache, norm=norm, log=log)

    if matrix_file is not None:
        import numpy as np
        np.savetxt(matrix_file, m)

    fig, ax = plt.subplots()
    im = ax.imshow(m, cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest')
    plt.colorbar(im)
    ax.set_axis_off()
    fig.savefig(output_file)
    plt.close(fig)


def aggregate_loops_parser():
    parser = argparse.ArgumentParser(
        prog="kaic aggregate_loop",
        description='Make a loop aggregate plot'
    )

    parser.add_argument(
        'hic',
        help='Hic file'
    )

    parser.add_argument(
        'loops',
        help='File with loop anchor regions (BEDPE)'
    )

    parser.add_argument(
        'output',
        help='Output image file (extension determines file format)'
    )

    parser.add_argument(
        '-p', '--pixels', dest='pixels',
        type=int,
        default=16,
        help='''Width of the output image in pixels. 
                Equal to the amount of bins in Hi-C matrix around loop
                Default: 16'''
    )

    parser.add_argument(
        '-N', '--no-norm', dest='norm',
        action='store_false',
        help='''Do not normalize matrix to expected values'''
    )
    parser.set_defaults(norm=True)

    parser.add_argument(
        '-L', '--no-log', dest='log',
        action='store_false',
        help='''Do not log2-transform normalized matrices.
                Only used in conjunction with '--norm'.'''
    )
    parser.set_defaults(log=True)

    parser.add_argument(
        '-C', '--no-cache', dest='cache',
        action='store_false',
        help='''Do not cache chromosome matrices (slower, but saves a lot of memory)'''
    )
    parser.set_defaults(cache=True)

    parser.add_argument(
        '-m', '--save-matrix', dest='matrix_file',
        help='''Path to save aggregate matrix (numpy txt format)'''
    )

    parser.add_argument(
        '-c', '--colormap', dest='colormap',
        help='''Matplotlib colormap to use for matrix'''
    )

    parser.add_argument(
        '--vmin', dest='vmin',
        type=float,
        help='''Minimum saturation value in image'''
    )

    parser.add_argument(
        '--vmax', dest='vmax',
        type=float,
        help='''Maximum saturation value in image'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)

    return parser


def aggregate_loops(argv):
    parser = aggregate_loops_parser()

    args = parser.parse_args(argv[2:])
    import os

    hic_file = os.path.expanduser(args.hic)
    loops_file = os.path.expanduser(args.loops)
    output_file = os.path.expanduser(args.output)
    pixels = args.pixels
    norm = args.norm
    log = args.log
    cache = args.cache
    matrix_file = None if args.matrix_file is None else os.path.expanduser(args.matrix_file)
    cmap = args.colormap
    vmin = args.vmin
    vmax = args.vmax
    tmp = args.tmp

    import kaic
    from kaic.config import config
    from kaic.architecture.hic_architecture import aggregate_loops
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import kaic.plotting
    from kaic.tools.general import human_format

    if cmap is None:
        cmap = 'bwr' if norm and log else config.colormap_hic

    with kaic.load(hic_file, mode='r', tmpdir=tmp) as hic:
        b = hic.bin_size
        loops = kaic.Bedpe(loops_file)
        m = aggregate_loops(hic, loops, pixels=pixels, cache=cache, norm=norm, log=log)

    if matrix_file is not None:
        import numpy as np
        np.savetxt(matrix_file, m)

    left = int(pixels / 2)
    right = left if pixels % 2 == 1 else left - 1

    fig, ax = plt.subplots()
    im = ax.imshow(m, cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest')
    plt.colorbar(im)
    ax.set_xticks([0, left, pixels - 1])
    ax.set_xticklabels(['-{}b'.format(human_format(left * b)), '', '+{}b'.format(human_format(right * b))])
    ax.set_yticks([0, left, pixels - 1])
    ax.set_yticklabels(['-{}b'.format(human_format(left * b)), '', '+{}b'.format(human_format(right * b))])
    fig.savefig(output_file)
    plt.close(fig)


def ab_profile_parser():
    parser = argparse.ArgumentParser(
        prog="kaic ab_profile",
        description='Make an A-B compartment interaction profile'
    )

    parser.add_argument(
        'hic',
        help='Hic file'
    )

    parser.add_argument(
        'genome',
        help="Can be an HDF5 Genome object, a FASTA file, "
             "a folder with FASTA files, or a "
             "comma-separated list of FASTA files."
    )

    parser.add_argument(
        'output',
        help='Output image file (extension determines file format)'
    )

    parser.add_argument(
        '-p', '--percentiles', dest='percentiles',
        type=float,
        nargs='+',
        default=(20.0, 40.0, 60.0, 80.0, 100.0),
        help='''Percentiles to use for calculation'''
    )

    parser.add_argument(
        '-c', '--colormap', dest='colormap',
        default='RdBu_r',
        help='''Matplotlib colormap to use for matrix'''
    )

    parser.add_argument(
        '--vmin', dest='vmin',
        type=float,
        default=-1,
        help='''Minimum saturation value in image'''
    )

    parser.add_argument(
        '--vmax', dest='vmax',
        type=float,
        default=1,
        help='''Maximum saturation value in image'''
    )

    parser.add_argument(
        '-C', '--no-chromosome', dest='per_chromosome',
        action='store_false',
        help='''Do not restrict calculation to intra-chromosomal regions'''
    )
    parser.set_defaults(per_chromosome=True)

    parser.add_argument(
        '-G', '--only-gc', dest='only_gc',
        action='store_true',
        help='''Only use GC content for domain calculation, 
                not the correlation matrix eigenvector.'''
    )
    parser.set_defaults(only_gc=False)

    parser.add_argument(
        '-m', '--save-matrix', dest='matrix_file',
        help='''Path to save aggregate matrix (numpy txt format)'''
    )

    parser.add_argument(
        '-x', '--exclude', dest='exclude',
        nargs='+',
        help='''Chromosome names to exclude from analysis'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Work in temporary directory'''
    )
    parser.set_defaults(tmp=False)

    return parser


def ab_profile(argv):
    parser = ab_profile_parser()

    args = parser.parse_args(argv[2:])
    import os

    hic_file = os.path.expanduser(args.hic)
    genome_file = os.path.expanduser(args.genome)
    output_file = os.path.expanduser(args.output)
    percentiles = args.percentiles
    per_chromosome = args.per_chromosome
    cmap = args.colormap
    matrix_file = None if args.matrix_file is None else os.path.expanduser(args.matrix_file)
    vmin = args.vmin
    vmax = args.vmax
    only_gc = args.only_gc
    exclude = args.exclude if args.exclude is not None else []
    tmp = args.tmp

    import kaic
    from kaic.architecture.hic_architecture import ab_enrichment_profile
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt

    with kaic.load(hic_file, mode='r', tmpdir=tmp) as hic:
        with kaic.Genome.from_string(genome_file, tmpdir=tmp, mode='r') as genome:
            m = ab_enrichment_profile(hic, genome, percentiles=percentiles,
                                      per_chromosome=per_chromosome, only_gc=only_gc,
                                      exclude_chromosomes=exclude)

    if matrix_file is not None:
        import numpy as np
        np.savetxt(matrix_file, m)

    fig, ax = plt.subplots()
    im = ax.imshow(m, cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest')
    cb = plt.colorbar(im)
    cb.set_label('log enrichment')
    ax.set_xticks([0, 4])
    ax.set_xticklabels(['active', 'inactive'])
    ax.set_yticks([0, 4])
    ax.set_yticklabels(['active', 'inactive'])
    fig.savefig(output_file)
    plt.close(fig)


def stats_parser():
    parser = argparse.ArgumentParser(
        prog="kaic stats",
        description='Get statistics on number of reads used at each step of a pipeline.'
    )

    parser.add_argument(
        'output',
        help="Output file (.txt) to store statistics."
    )

    parser.add_argument(
        '-f', '--fastq', dest='fastq',
        nargs='+',
        help='''List of FASTQ files or folders containing FASTQ files.'''
    )

    parser.add_argument(
        '-r', '--reads', dest='reads',
        nargs='+',
        help='''List of Reads files or folders containing Reads files ('.reads ending').'''
    )

    parser.add_argument(
        '-p', '--pairs', dest='pairs',
        nargs='+',
        help='''List of Pairs files or folders containing Pairs files ('.pairs ending').'''
    )

    parser.add_argument(
        '-c', '--hic', dest='hic',
        nargs='+',
        help='''List of Hic files or folders containing Hic files ('.hic ending').'''
    )
    return parser


def stats(argv):
    parser = stats_parser()

    args = parser.parse_args(argv[2:])

    from collections import defaultdict

    output_file = os.path.expanduser(args.output)
    with open(output_file, 'w') as o:
        o.write("type\tfile\tproperty\tcount\n")

    def get_files(paths, endings=()):
        files = []
        for p in paths:
            ppath = os.path.expanduser(p)
            if os.path.isdir(ppath):
                for path in os.listdir(ppath):
                    full_path = ppath + '/{}'.format(path)
                    if os.path.isfile(full_path):
                        for ending in endings:
                            if path.endswith(ending):
                                files.append(full_path)
                                continue
            elif os.path.isfile(ppath):
                files.append(ppath)
        return files

    def stats(maskable, masked_table):
        import tables as t
        statistics = maskable.mask_statistics(masked_table)

        # calculate total
        if isinstance(masked_table, t.Group):
            total = 0
            for table in masked_table:
                total += table._original_len()
        else:
            total = masked_table._original_len()
        return statistics, total

    # 1. Count number of reads
    if args.fastq is not None:
        logger.info("Processing FASTQ files.")
        import gzip

        fastq_files = get_files(args.fastq, ('.fq', '.fastq', '.fq.gz', '.fastq.gz', 'fq.gzip', 'fastq.gzip'))

        total_count = 0
        for fastq_file in fastq_files:
            logger.info("{}".format(fastq_file))
            if fastq_file.endswith('gz') or fastq_file.endswith('gzip'):
                read = gzip.open
            else:
                read = open

            with read(fastq_file, 'r') as f:
                line_count = sum(1 for line in f)
            total_count += line_count/4

            with open(output_file, 'a') as o:
                o.write("fastq\t{}\tcount\t{}\n".format(fastq_file, line_count/4))

        with open(output_file, 'a') as o:
            o.write("fastq\ttotal\tcount\t{}\n".format(total_count))

    # 2. Reads statistics
    if args.reads is not None:
        logger.info("Processing Reads files.")

        from kaic.construct.seq import Reads
        reads_files = get_files(args.reads, ('.reads',))

        reads_summary = defaultdict(int)
        for reads_file in reads_files:
            logger.info("{}".format(reads_file))
            reads = Reads(reads_file, mode='r')
            statistics, total = stats(reads, reads._reads)

            with open(output_file, 'a') as o:
                for key in sorted(statistics.keys()):
                    o.write("reads\t{}\t{}\t{}\n".format(reads_file, key, statistics[key]))
                    reads_summary[key] += statistics[key]

            with open(output_file, 'a') as o:
                o.write("reads\t{}\ttotal\t{}\n".format(reads_file, total))
                reads_summary['total'] += total
            reads_summary['filtered'] += total - statistics['unmasked']
            reads_summary['remaining'] += statistics['unmasked']

        with open(output_file, 'a') as o:
            for key in sorted(reads_summary.keys()):
                if key != 'filtered' and key != 'remaining':
                    o.write("reads\ttotal\t{}\t{}\n".format(key, reads_summary[key]))
            o.write("reads\ttotal\tfiltered\t{}\n".format(reads_summary['filtered']))
            o.write("reads\ttotal\tremaining\t{}\n".format(reads_summary['remaining']))

    # 3. Pairs statistics
    if args.pairs is not None:
        logger.info("Processing Pairs files.")
        import kaic
        pairs_files = get_files(args.pairs, ('.pairs',))
        pairs_summary = defaultdict(int)
        for pairs_file in pairs_files:
            logger.info("{}".format(pairs_file))
            pairs = kaic.Pairs(pairs_file, mode='r')
            statistics, total = stats(pairs, pairs._pairs)

            with open(output_file, 'a') as o:
                for key in sorted(statistics.keys()):
                    o.write("pairs\t{}\t{}\t{}\n".format(pairs_file, key, statistics[key]))
                    pairs_summary[key] += statistics[key]

            with open(output_file, 'a') as o:
                o.write("pairs\t{}\ttotal\t{}\n".format(pairs_file, total))
                pairs_summary['total'] += total

            pairs_summary['filtered'] += total - statistics['unmasked']
            pairs_summary['remaining'] += statistics['unmasked']

        with open(output_file, 'a') as o:
            for key in sorted(pairs_summary.keys()):
                if key != 'filtered' and key != 'remaining':
                    o.write("pairs\ttotal\t{}\t{}\n".format(key, pairs_summary[key]))
            o.write("pairs\ttotal\tfiltered\t{}\n".format(pairs_summary['filtered']))
            o.write("pairs\ttotal\tremaining\t{}\n".format(pairs_summary['remaining']))

    # 3. Hic statistics
    if args.hic is not None:
        logger.info("Processing Hic files.")
        from kaic.data.genomic import load_hic
        hic_files = get_files(args.hic, ('.hic',))

        hic_summary = defaultdict(int)
        for hic_file in hic_files:
            logger.info("{}".format(hic_file))
            hic = load_hic(hic_file, mode='r')
            statistics, total = stats(hic, hic._edges)

            with open(output_file, 'a') as o:
                for key in sorted(statistics.keys()):
                    o.write("hic\t{}\t{}\t{}\n".format(hic_file, key, statistics[key]))
                    hic_summary[key] += statistics[key]

            with open(output_file, 'a') as o:
                o.write("hic\t{}\ttotal\t{}\n".format(hic_file, total))
                hic_summary['total'] = total

            hic_summary['filtered'] += total - statistics['unmasked']
            hic_summary['remaining'] += statistics['unmasked']

        with open(output_file, 'a') as o:
            for key in sorted(hic_summary.keys()):
                if key != 'filtered' and key != 'remaining':
                    o.write("hic\ttotal\t{}\t{}\n".format(key, hic_summary[key]))
            o.write("hic\ttotal\tfiltered\t{}\n".format(hic_summary['filtered']))
            o.write("hic\ttotal\tremaining\t{}\n".format(hic_summary['remaining']))


def write_config_parser():
    parser = argparse.ArgumentParser(
        prog="kaic write_config",
        description='Write default config file to specified location.'
    )

    parser.add_argument(
        'config_file',
        help="Output file for default configuration."
    )

    parser.add_argument(
        '-f', '--force', dest='force',
        action='store_true',
        help='''Force overwrite of existing config file.'''
    )
    parser.set_defaults(force=False)
    return parser


def write_config(argv):
    parser = write_config_parser()

    args = parser.parse_args(argv[2:])

    from kaic.config import write_default_config
    write_default_config(os.path.expanduser(args.config_file), overwrite=args.force)


def cis_trans_parser():
    parser = argparse.ArgumentParser(
        prog="kaic cis_trans",
        description='Calculate cis/trans ratio of this Hi-C object.'
    )

    parser.add_argument(
        'hic',
        nargs='+',
        help="Hic object(s) for cis/trans calculation."
    )

    parser.add_argument(
        '-o', '--output', dest='output',
        help="Output file."
    )

    parser.add_argument(
        '-n', '--norm', dest='normalise',
        action='store_true',
        help='''Normalise ratio to the prior ratio of possible cis / trans contacts.'''
    )
    parser.set_defaults(normalise=False)
    return parser


def cis_trans(argv):
    parser = cis_trans_parser()

    args = parser.parse_args(argv[2:])
    hic_files = [os.path.expanduser(f) for f in args.hic]
    output_file = os.path.expanduser(args.output) if args.output is not None else None
    normalise = args.normalise

    import kaic
    from kaic.architecture.hic_architecture import cis_trans_ratio

    if output_file:
        with open(output_file, 'w') as o:
            o.write("file\tcis\ttrans\tratio\tfactor\n")

    for hic_file in hic_files:
        hic = kaic.load(hic_file, mode='r')

        r, cis, trans, f = cis_trans_ratio(hic, normalise)

        if output_file:
            with open(output_file, 'a') as o:
                o.write("{}\t{}\t{}\t{:.3f}\t{:.3f}\n".format(hic_file, cis, trans, r, f))
        print("{}".format(hic_file))
        print("\tcis: {}".format(cis))
        print("\ttrans: {}".format(trans))
        print("\tratio: {:.3f}".format(r))
        print("\tfactor: {:.3f}".format(f))
        hic.close()
