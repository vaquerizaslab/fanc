import argparse
import logging
import subprocess

# configure logging
logger = logging.getLogger(__name__)


def auto_parser():
    parser = argparse.ArgumentParser(
        prog="kaic auto",
        description='Automatically process an entire Hi-C data set'
    )

    parser.add_argument(
        'input',
        nargs='+',
        help='''Input files. kaic will try to guess the file by its extension.'''
    )

    parser.add_argument(
        'output_folder',
        help='''Folder where output files and sub-folders will be generated'''
    )

    parser.add_argument(
        '-g', '--genome', dest='genome',
        help='''Can be an HDF5 Genome object, a FASTA file,
                a folder with FASTA files, or a comma-separated
                list of FASTA files.'''
    )

    parser.add_argument(
        '-r', '--restriction-enzyme', dest='restriction_enzyme',
        help='''Restriction enzyme used for digestion (e.g. HindIII, case-sensitive)'''
    )

    parser.add_argument(
        '-i', '--genome-index', dest='genome_index',
        help='''Bowtie 2 genome index. Only required when passing FASTQ files as input'''
    )

    parser.add_argument(
        '-n', '--basename', dest='basename',
        help='''Basename for output files. If not provided, will be guessed based on input file names'''
    )

    parser.add_argument(
        '-s', '--step-size', dest='step_size',
        type=int,
        default=3,
        help='''Step size for iterative mapping. Default: 3'''
    )

    parser.add_argument(
        '-b', '--bin-sizes', dest='bin_sizes',
        type=int,
        nargs='+',
        default=[5000000, 2000000, 1000000, 500000, 250000, 100000, 50000,
                 25000, 10000, 5000],
        help='''Bin sizes for Hi-C matrix generation. Default: 5000000, 2000000, 1000000,
                500000, 250000, 100000, 50000, 25000, 10000, 5000'''
    )

    parser.add_argument(
        '-t', '--threads', dest='threads',
        type=int,
        default=1,
        help='''Maximum number of threads to use for the analysis.'''
    )

    parser.add_argument(
        '-O', '--no-optimise', dest='optimise',
        action='store_false',
        help='''Produce a Hi-C object optimised for fast access times. May impact compatibility.'''
    )
    parser.set_defaults(optimise=True)

    parser.add_argument(
        '--bowtie_parallel', dest='bowtie_parallel',
        action='store_true',
        help='''Use Bowtie2 parallelisation instead of spawning multiple Bowtie 2 processes.
                Slower, but less of a memory overhead.'''
    )
    parser.set_defaults(bowtie_parallel=False)

    parser.add_argument(
        '--split-fastq', dest='split_fastq',
        action='store_true',
        help='''Split fastq files into chunks of 10M reads. These will be merged again on the Reads level.
                Splitting and merging bypasses the -tmp flag. This option reduces disk usage in tmp, in case
                the system has a small tmp partition.'''
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
        '--ice', dest='ice',
        action='store_true',
        help='''Use ICE iterative matrix balancing rather than Knight-Ruiz.'''
    )
    parser.set_defaults(ice=False)

    parser.add_argument(
        '--le-inward-cutoff', dest='inward_cutoff',
        type=int,
        help='''Ligation error inward cutoff (if not set, will be determined automatically).'''
    )

    parser.add_argument(
        '--le-outward-cutoff', dest='outward_cutoff',
        type=int,
        help='''Ligation error outward cutoff (if not set, will be determined automatically).'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        help='''Copy original files to temporary directory. Reduces network I/O.'''
    )
    parser.set_defaults(copy=False)

    parser.add_argument(
        '--reads-intermediate', dest='reads_intermediate',
        action='store_true',
        help='''Use '.reads' file intermediates. Will use more time and disk space,
                    but could provide more control over filtering and downstream processing.'''
    )
    parser.set_defaults(reads_intermediate=False)

    parser.add_argument(
        '--no-iterative', dest='iterative',
        action='store_false',
        help='''Do not map reads iteratively, use simple mapping.'''
    )
    parser.set_defaults(iterative=True)

    parser.add_argument(
        '--no-sam-sort', dest='sam_sort',
        action='store_false',
        help='''Do not sort SAM/BAM files.'''
    )
    parser.set_defaults(sam_sort=True)

    parser.add_argument(
        '--restore-coverage', dest='restore_coverage',
        action='store_true',
        help='''Restore coverage to the original total number of reads. 
                    Otherwise matrix entries will be contact probabilities.
                    Only available for KR matrix balancing.'''
    )
    parser.set_defaults(restore_coverage=False)

    parser.add_argument(
        '--split-ligation-junction', dest='split_ligation_junction',
        action='store_true',
        help='''Split reads at predicted ligation junction before mapping.'''
    )
    parser.set_defaults(split_ligation_junction=False)

    parser.add_argument(
        '--bwa', dest='bwa',
        action='store_true',
        help='''Use BWA as mapper.'''
    )
    parser.set_defaults(bwa=False)

    return parser


def reads_worker(file_names, reads_file, args):
    load_reads_command = ['kaic', 'load_reads', '-D']
    if args.tmp:
        load_reads_command.append('-tmp')

    if args.split_fastq:
        load_reads_command.append('--split-sam')

    for file_name in file_names:
        load_reads_command.append(file_name)
    load_reads_command.append(reads_file)

    return subprocess.call(load_reads_command)


def filtered_reads_worker(reads_file, filtered_reads_file, filtered_reads_stats_file, args):
    filter_reads_command = ['kaic', 'filter_reads', '-m', '-q', '30']
    if not args.bwa:
        filter_reads_command.append('-us')
    if args.tmp:
        filter_reads_command.append('-tmp')
    filter_reads_command.append('-s')

    return subprocess.call(filter_reads_command + [filtered_reads_stats_file, reads_file, filtered_reads_file])


def pairs_worker(pairs_file, filtered_reads_file1, filtered_reads_file2, genome, restriction_enzyme, args):
    logger.info("Creating Pairs object...")
    pairs_command = ['kaic', 'reads_to_pairs', filtered_reads_file1, filtered_reads_file2,
                     genome, restriction_enzyme, pairs_file]
    if args.tmp:
        pairs_command.append('-tmp')
    return subprocess.call(pairs_command)


def sam_to_pairs_worker(sam1_file, sam2_file, genome_file, restriction_enzyme, pairs_file, args):
    logger.info("Creating Pairs object...")
    pairs_command = ['kaic', 'sam_to_pairs', sam1_file, sam2_file, genome_file,
                     restriction_enzyme, pairs_file,
                     '-m', '-q', '30']
    if not args.bwa:
        pairs_command.append('-us')

    if args.tmp:
        pairs_command.append('-tmp')

    if args.sam_sort:
        pairs_command.append('-S')

    return subprocess.call(pairs_command)


def sam_sort_worker(sam_file, output_file, args):
    logger.info("Sorting SAM file...")
    sort_command = ['kaic', 'sort_sam', sam_file, output_file]
    if args.tmp:
        sort_command.append('-tmp')

    return subprocess.call(sort_command)


def pairs_ligation_error_worker(pairs_file, ligation_error_file):
    ligation_error_command = ['kaic', 'plot_ligation_err']
    return subprocess.call(ligation_error_command + [pairs_file, ligation_error_file])


def pairs_re_dist_worker(pairs_file, re_dist_file):
    re_dist_command = ['kaic', 'plot_re_dist']
    return subprocess.call(re_dist_command + [pairs_file, re_dist_file])


def filtered_pairs_worker(pairs_file, filtered_pairs_file, filtered_pairs_stats_file, args):
    filter_pairs_command = ['kaic', 'filter_pairs', '-r', '5000', '-l', '-d', '2']

    if args.inward_cutoff is not None:
        filter_pairs_command += ['-i', str(args.inward_cutoff)]
        if args.outward_cutoff is None:
            filter_pairs_command += ['-o', '0']

    if args.outward_cutoff is not None:
        filter_pairs_command += ['-o', str(args.outward_cutoff)]
        if args.inward_cutoff is None:
            filter_pairs_command += ['-i', '0']

    if args.inward_cutoff is None and args.outward_cutoff is None:
        filter_pairs_command += ['--auto']

    if args.tmp:
        filter_pairs_command.append('-tmp')
    filter_pairs_command.append('-s')

    p1 = subprocess.call(filter_pairs_command + [filtered_pairs_stats_file, pairs_file,
                                                 filtered_pairs_file])
    if p1 != 0:
        logger.error("Filtering failed for some reason, trying again with fixed thresholds...")
        filter_pairs_command = ['kaic', 'filter_pairs', '-i', '10000',
                                '-o', '10000', '-r', '5000', '-l', '-d', '2']
        if args.tmp:
            filter_pairs_command.append('-tmp')
        filter_pairs_command.append('-s')
        p1 = subprocess.call(filter_pairs_command + [filtered_pairs_stats_file, pairs_file,
                                                     filtered_pairs_file])
    return p1


def hic_worker(pairs_file, hic_file, args):
    hic_command = ['kaic', 'pairs_to_hic']
    if args.tmp:
        hic_command.append('-tmp')

    return subprocess.call(hic_command + [pairs_file, hic_file])


def batch_hic_worker(hic_file, bin_size, binned_hic_file, filtered_hic_file, filtered_hic_stats_file,
                     corrected_hic_file, chromosome_corrected_hic_file, args):
    """
    batch worker to:
    * bin
    * filter
    * correct
    """
    import kaic
    logger.info("Binning Hic {} at {}...".format(hic_file, bin_size))
    bin_hic_command = ['kaic', 'bin_hic']
    if args.tmp:
        bin_hic_command.append('-tmp')

    ret1 = subprocess.call(bin_hic_command + [hic_file, binned_hic_file, str(bin_size)])
    if ret1 != 0:
        return ret1

    logger.info("Filtering Hic {}...".format(binned_hic_file))
    filter_hic_command = ['kaic', 'filter_hic']
    if args.tmp:
        filter_hic_command.append('-tmp')
    filter_hic_command.append('-rl')
    filter_hic_command.append('0.1')
    filter_hic_command.append('-s')

    ret2 = subprocess.call(filter_hic_command + [filtered_hic_stats_file, binned_hic_file, filtered_hic_file])
    if ret2 != 0:
        return ret2

    logger.info("Correcting Hic {}...".format(filtered_hic_file))
    correct_hic_command = ['kaic', 'correct_hic']
    if args.tmp:
        correct_hic_command.append('-tmp')
    if not args.optimise:
        correct_hic_command.append('-O')
    if args.ice:
        correct_hic_command.append('-i')
    if args.restore_coverage:
        correct_hic_command.append('-r')

    ret3 = subprocess.call(correct_hic_command + ['-c', filtered_hic_file, chromosome_corrected_hic_file])
    if ret3 != 0:
        return ret3

    hic = kaic.load(filtered_hic_file, mode='r')
    n_regions = len(hic.regions)
    hic.close()

    if n_regions <= 50000:
        ret4 = subprocess.call(correct_hic_command + [filtered_hic_file, corrected_hic_file])
    else:
        ret4 = 0

    if ret4 != 0:
        return ret4

    return 0


def auto(argv):
    import os
    from kaic.tools.general import which

    parser = auto_parser()
    args = parser.parse_args(argv[2:])

    bin_sizes = args.bin_sizes
    split_ligation_junction = args.split_ligation_junction
    restriction_enzyme = args.restriction_enzyme
    bwa = args.bwa

    def is_fastq_file(file_name):
        base, extension = os.path.splitext(file_name)
        if extension in ['.gz', '.gzip']:
            base, extension = os.path.splitext(base)

        if extension in ['.fq', '.fastq']:
            return True
        return False

    def is_sam_or_bam_file(file_name):
        _, extension = os.path.splitext(file_name)
        return extension in ['.sam', '.bam']

    def is_reads_file(file_name):
        _, extension = os.path.splitext(file_name)
        return extension in ['.reads']

    def is_pairs_file(file_name):
        _, extension = os.path.splitext(file_name)
        return extension in ['.pairs']

    def is_hic_file(file_name):
        _, extension = os.path.splitext(file_name)
        return extension in ['.hic']

    def file_type(file_name):
        if is_fastq_file(file_name):
            return 'fastq'
        if is_sam_or_bam_file(file_name):
            return 'sam'
        if is_reads_file(file_name):
            return 'reads'
        if is_pairs_file(file_name):
            return 'pairs'
        if is_hic_file(file_name):
            return 'hic'
        return None

    def file_basename(file_name):
        basename = os.path.basename(os.path.splitext(file_name)[0])
        if file_name.endswith('.gz') or file_name.endswith('.gzip'):
            basename = os.path.splitext(basename)[0]
        return basename

    file_names = [os.path.expanduser(file_name) for file_name in args.input]
    file_types = [file_type(file_name) for file_name in file_names]
    file_basenames = [file_basename(file_name) for file_name in file_names]

    for i in range(len(file_types)):
        if file_types[i] not in ('fastq', 'sam', 'reads', 'pairs', 'hic'):
            import kaic
            try:
                ds = kaic.load(file_names[i], mode='r')
                if isinstance(ds, kaic.Hic) or isinstance(ds, kaic.AccessOptimisedHic):
                    file_types[i] = 'hic'
                elif isinstance(ds, kaic.Pairs) or isinstance(ds, kaic.FragmentMappedReadPairs):
                    file_types[i] = 'pairs'
                elif isinstance(ds, kaic.Reads):
                    file_types[i] = 'reads'
                else:
                    raise ValueError("Could not detect file type using kaic load.")
            except ValueError:
                raise ValueError("Not a valid input file type: {}".format(file_type))

    if args.basename is None:
        if len(file_basenames) == 1:
            basename = file_basenames[0]
        else:
            basename = []
            for pair in zip(*file_basenames):
                if pair[0] == pair[1]:
                    basename.append(pair[0])
                else:
                    break
            if len(basename) == 0:
                basename = file_basenames[0]
            else:
                if basename[-1] in ['.', '_']:
                    basename = "".join(basename[:-1])
                else:
                    basename = "".join(basename)
    else:
        basename = args.basename

    output_folder = os.path.expanduser(args.output_folder)
    if not output_folder[-1] == '/':
        output_folder += '/'

    threads = args.threads

    # 0. Do some sanity checks on required flags
    if 'fastq' in file_types:
        if args.genome_index is None:
            print("Error: Must provide genome index (-i) when mapping FASTQ files!")
            quit(1)
        else:
            check_path = os.path.expanduser(args.genome_index)
            if check_path.endswith('.'):
                check_path = check_path[:-1]

            if not bwa:
                for i in range(1, 5):
                    if not os.path.exists(check_path + '.{}.bt2'.format(i)):
                        raise ValueError("Bowtie2 index incomplete, check index files for completeness.")
                for i in range(1, 3):
                    if not os.path.exists(check_path + '.rev.{}.bt2'.format(i)):
                        raise ValueError("Bowtie2 index incomplete, check index files for completeness.")
            else:
                for ending in ('amb', 'ann', 'bwt', 'pac', 'sa'):
                    if not os.path.exists(check_path + '.{}'.format(ending)):
                        raise ValueError("BWA index incomplete, check index files for completeness.")

        if not which('bowtie2'):
            raise ValueError("bowtie2 must be in PATH for iterative mapping!")

    if 'fastq' in file_types or 'sam' in file_types or 'reads' in file_types:
        if args.genome is None:
            print("Error: Must provide genome (-g) to process read pair files!")
            quit(1)

        if args.restriction_enzyme is None:
            print("Error: Must provide restriction enzyme (-r) to process read pair files!")
            quit(1)
        else:
            from Bio import Restriction
            try:
                getattr(Restriction, args.restriction_enzyme)
            except AttributeError:
                raise ValueError("restriction_enzyme string is not recognized: %s" % args.restriction_enzyme)

    logger.info("Output folder: %s" % output_folder)
    logger.info("Input files: %s" % str(file_names))
    logger.info("Input file types: %s" % str(file_types))

    if args.basename:
        logger.info("Final basename: %s" % basename)
    else:
        logger.info("Final basename: %s (you can change this with the -n option!)" % basename)

    import subprocess
    # from multiprocessing.pool import ThreadPool as Pool
    from multiprocessing import Pool

    # 1. create default folders in root directory
    logger.info("Creating output folders...")
    rc = subprocess.call(['kaic', 'dirs', output_folder])
    if rc != 0:
        print("Creating folders failed for some reason, aborting...")
        quit(rc)

    # 2. If input files are (gzipped) FASTQ, map them iteratively first
    def mapping_worker(file_name, index, bam_file, mapping_threads=1):
        iterative_mapping_command = ['kaic', 'map',
                                     '-m', '25', '-s', str(args.step_size),
                                     '-t', str(mapping_threads)]

        if not bwa:
            iterative_mapping_command += ['-q', '30']
        else:
            iterative_mapping_command.append('--bwa')

        if args.tmp:
            iterative_mapping_command.append('-tmp')

        if args.bowtie_parallel:
            iterative_mapping_command.append('--bowtie-parallel')
        if args.split_fastq:
            iterative_mapping_command.append('--split-fastq')
        if not args.memory_map:
            iterative_mapping_command.append('--no-memory-map')
        if not args.iterative:
            iterative_mapping_command.append('--simple')
        if split_ligation_junction:
            iterative_mapping_command.append('--restriction-enzyme')
            iterative_mapping_command.append(restriction_enzyme)

        return subprocess.call(iterative_mapping_command + [file_name, index, bam_file])

    fastq_files = []
    for i in range(len(file_names)):
        if file_types[i] != 'fastq':
            continue
        fastq_files.append(i)

    if len(fastq_files) > 0:
        mapping_processes = threads

        index = os.path.expanduser(args.genome_index)
        if index.endswith('.'):
            index = index[:-1]
        logger.info("Iteratively mapping FASTQ files...")

        bam_files = []

        fastq_results = []
        for i, ix in enumerate(fastq_files):
            bam_file = output_folder + 'sam/' + file_basenames[ix] + '.bam'
            bam_files.append(bam_file)
            fastq_results.append(mapping_worker(file_names[ix], index, bam_file, mapping_processes))

        for rt in fastq_results:
            if rt != 0:
                raise RuntimeError("Bowtie mapping had non-zero exit status")

        for ix, i in enumerate(fastq_files):
            file_names[i] = bam_files[ix]
            file_types[i] = 'sam'

    if args.reads_intermediate:
        # 3. SAM/BAM to Reads object conversion
        sam_files = []
        for i in range(len(file_names)):
            if file_types[i] != 'sam':
                continue
            sam_files.append(i)

        if len(sam_files) > 0:
            tp = Pool(threads if not args.split_fastq else 1)

            reads_files = []
            reads_results = []
            for ix in sam_files:
                reads_file = output_folder + 'reads/' + file_basenames[ix] + '.reads'
                reads_files.append(reads_file)
                rt = tp.apply_async(reads_worker, ([file_names[ix]], reads_file, args))
                reads_results.append(rt)
            tp.close()
            tp.join()

            for rt in reads_results:
                if rt.get() != 0:
                    raise RuntimeError("Read loading from SAM/BAM had non-zero exit status")

            for ix, i in enumerate(sam_files):
                file_names[i] = reads_files[ix]
                file_types[i] = 'reads'

        # 4. Filter reads
        reads_files = []
        for i in range(len(file_names)):
            if file_types[i] != 'reads':
                continue
            reads_files.append(i)

        if len(reads_files) > 0:
            tp = Pool(threads)

            filtered_reads_files = []
            filter_reads_results = []
            for ix in reads_files:
                filtered_reads_file = output_folder + 'reads/filtered/' + file_basenames[ix] + '_filtered.reads'
                filtered_reads_stats_file = output_folder + 'plots/stats/' + file_basenames[ix] + '.reads.stats.pdf'
                filtered_reads_files.append(filtered_reads_file)

                rt = tp.apply_async(filtered_reads_worker,
                                    (file_names[ix], filtered_reads_file, filtered_reads_stats_file, args))
                filter_reads_results.append(rt)
            tp.close()
            tp.join()

            for rt in filter_reads_results:
                if rt.get() != 0:
                    raise RuntimeError("Read filtering had non-zero exit status")

            for ix, i in enumerate(reads_files):
                file_names[i] = filtered_reads_files[ix]

        # 5. Reads to Pairs
        reads_file_pairs = []
        i = 0
        while i < len(file_names):
            if file_types[i] == 'reads':
                if not file_types[i + 1] == 'reads':
                    raise RuntimeError("Cannot create read pairs, because %s is missing a partner file" % file_names[i])
                reads_file_pairs.append((i, i + 1))
                i += 1
            i += 1

        # get reads file pair basenames
        pair_basenames = [basename + '_' + str(i) for i in range(len(reads_file_pairs))]

        if len(reads_file_pairs) > 0:
            tp = Pool(threads)
            genome = args.genome
            restriction_enzyme = args.restriction_enzyme

            pairs_files = []
            pairs_results = []
            for i, j in reads_file_pairs:
                if len(reads_file_pairs) > 1:
                    pairs_file = output_folder + 'pairs/' + pair_basenames[len(pairs_files)] + '.pairs'
                else:
                    pairs_file = output_folder + 'pairs/' + basename + '.pairs'
                rt = tp.apply_async(pairs_worker,
                                    (pairs_file, file_names[i], file_names[j], genome, restriction_enzyme, args))
                pairs_results.append(rt)
                pairs_files.append(pairs_file)
            tp.close()
            tp.join()

            for rt in pairs_results:
                if rt.get() != 0:
                    raise RuntimeError("Pairs loading from reads had non-zero exit status")

            for ix, read_pair in enumerate(reversed(reads_file_pairs)):
                file_names[read_pair[0]] = pairs_files[ix]
                del file_names[read_pair[1]]
                file_types[read_pair[0]] = 'pairs'
                del file_types[read_pair[1]]
    else:
        if args.sam_sort:
            # sort SAM files
            sam_files = []
            for i in range(len(file_names)):
                if file_types[i] != 'sam':
                    continue
                sam_files.append(i)

            if len(sam_files) > 0:
                tp = Pool(threads)

                sorted_sam_files = []
                sort_results = []
                for ix in sam_files:
                    sam_path, sam_extension = os.path.splitext(file_names[ix])
                    sam_basename = os.path.basename(sam_path)
                    sorted_sam_file = os.path.join(output_folder, 'sam', sam_basename + '_sort' + sam_extension)
                    sorted_sam_files.append(sorted_sam_file)
                    rt = tp.apply_async(sam_sort_worker, (file_names[ix], sorted_sam_file, args))
                    sort_results.append(rt)
                tp.close()
                tp.join()

                for rt in sort_results:
                    if rt.get() != 0:
                        raise RuntimeError("Read loading from SAM/BAM had non-zero exit status")

                for ix, i in enumerate(sam_files):
                    file_names[i] = sorted_sam_files[ix]

        # load pairs directly from SAM
        sam_file_pairs = []
        i = 0
        while i < len(file_names):
            if file_types[i] == 'sam':
                if not file_types[i + 1] == 'sam':
                    raise RuntimeError("Cannot create SAM pairs, because %s is missing a partner file" % file_names[i])
                sam_file_pairs.append((i, i + 1))
                i += 1
            i += 1

        pair_basenames = [basename + '_' + str(i) for i in range(len(sam_file_pairs))]

        tp = Pool(threads)
        genome = args.genome
        restriction_enzyme = args.restriction_enzyme

        pairs_files = []
        pairs_results = []
        for i, j in sam_file_pairs:
            if len(sam_file_pairs) > 1:
                pairs_file = output_folder + 'pairs/' + pair_basenames[len(pairs_files)] + '.pairs'
            else:
                pairs_file = output_folder + 'pairs/' + basename + '.pairs'
            rt = tp.apply_async(sam_to_pairs_worker,
                                (file_names[i], file_names[j], genome, restriction_enzyme, pairs_file, args))
            pairs_results.append(rt)
            pairs_files.append(pairs_file)
        tp.close()
        tp.join()

        for rt in pairs_results:
            if rt.get() != 0:
                raise RuntimeError("Pairs loading from reads had non-zero exit status")

        for ix, sam_pair in enumerate(reversed(sam_file_pairs)):
            file_names[sam_pair[0]] = pairs_files[ix]
            del file_names[sam_pair[1]]
            file_types[sam_pair[0]] = 'pairs'
            del file_types[sam_pair[1]]

    # 7. Pairs stats and filtering
    pairs_files = []
    for i in range(len(file_names)):
        if file_types[i] != 'pairs':
            continue
        pairs_files.append(i)

    if len(pairs_files) > 0:
        tp = Pool(threads)

        filtered_pairs_files = []
        filter_pairs_results = []
        for ix in pairs_files:
            pair_basename = os.path.basename(os.path.splitext(file_names[ix])[0])
            filtered_pairs_file = output_folder + 'pairs/filtered/' + pair_basename + '_filtered.pairs'
            filtered_pairs_stats_file = output_folder + 'plots/stats/' + pair_basename + '.pairs.stats.pdf'
            ligation_error_file = output_folder + 'plots/stats/' + pair_basename + '.pairs.ligation_error.pdf'
            re_dist_file = output_folder + 'plots/stats/' + pair_basename + '.pairs.re_dist.pdf'

            tp.apply_async(pairs_ligation_error_worker, (file_names[ix], ligation_error_file))
            tp.apply_async(pairs_re_dist_worker, (file_names[ix], re_dist_file))
            rt = tp.apply_async(filtered_pairs_worker,
                                (file_names[ix], filtered_pairs_file, filtered_pairs_stats_file, args))
            filter_pairs_results.append(rt)

            filtered_pairs_files.append(filtered_pairs_file)
        tp.close()
        tp.join()

        for rt in filter_pairs_results:
            if rt.get() != 0:
                raise RuntimeError("Pair filtering had non-zero exit status")

        for ix, i in enumerate(pairs_files):
            file_names[i] = filtered_pairs_files[ix]

    # 8. Pairs to Hic
    pairs_files = []
    for i in range(len(file_names)):
        if file_types[i] != 'pairs':
            continue
        pairs_files.append(i)

    if len(pairs_files) > 0:
        tp = Pool(threads)

        hic_files = []
        hic_results = []
        for ix in pairs_files:
            hic_basename = os.path.basename(os.path.splitext(file_names[ix])[0])
            if hic_basename.endswith('_filtered'):
                hic_basename = hic_basename[:-9]
            hic_file = output_folder + 'hic/' + hic_basename + '.hic'

            rt = tp.apply_async(hic_worker, (file_names[ix], hic_file, args))
            hic_results.append(rt)

            hic_files.append(hic_file)
        tp.close()
        tp.join()

        for rt in hic_results:
            if rt.get() != 0:
                raise RuntimeError("Hi-C conversion had non-zero exit status")

        for ix, i in enumerate(pairs_files):
            file_names[i] = hic_files[ix]
            file_types[i] = 'hic'

    # 9. Merge Hic
    hic_files = []
    for i in range(len(file_names)):
        if file_types[i] != 'hic':
            continue
        hic_files.append(i)

    if len(hic_files) > 1:
        output_hic = output_folder + 'hic/' + basename + '.hic'
        logger.info("Merging Hi-C files...")
        merge_hic_command = ['kaic', 'merge_hic']
        if args.tmp:
            merge_hic_command.append('-tmp')

        if not args.optimise:
            merge_hic_command.append('-O')

        hics = [file_names[i] for i in hic_files]
        rt = subprocess.call(merge_hic_command + hics + [output_hic])

        if rt != 0:
            raise RuntimeError("Hi-C merge had non-zero exit status")

        file_names[hic_files[0]] = output_hic
        hic_files.pop(0)
        for ix, i in enumerate(reversed(hic_files)):
            del file_names[i]
            del file_types[i]

    from kaic.tools.general import human_format

    hic_files = []
    for i in range(len(file_names)):
        if file_types[i] != 'hic':
            continue
        hic_files.append(i)

    if len(hic_files) > 0:
        tp = Pool(threads)

        hic_results = []
        for ix in hic_files:
            binned_hic_file_base = output_folder + 'hic/binned/' + basename + '_'

            for bin_size in bin_sizes:
                bin_size_str = human_format(bin_size, 0).lower() + 'b'
                binned_hic_file = binned_hic_file_base + bin_size_str + '.hic'
                hic_basename = os.path.basename(os.path.splitext(binned_hic_file)[0])
                filtered_hic_file = output_folder + 'hic/filtered/' + hic_basename + '_filtered.hic'
                filtered_hic_stats_file = output_folder + 'plots/stats/' + hic_basename + '_filtered.stats.pdf'
                chromosome_corrected_hic_file = output_folder + 'hic/corrected/' + hic_basename + '_corrected_pc.hic'
                corrected_hic_file = output_folder + 'hic/corrected/' + hic_basename + '_corrected.hic'

                rt = tp.apply_async(batch_hic_worker, (file_names[ix],  # hic_file
                                                       bin_size,
                                                       binned_hic_file,
                                                       filtered_hic_file,
                                                       filtered_hic_stats_file,
                                                       corrected_hic_file,
                                                       chromosome_corrected_hic_file,
                                                       args))
                hic_results.append(rt)
        tp.close()
        tp.join()

    for rt in hic_results:
        if rt.get() != 0:
            raise RuntimeError("Hi-C binning/filtering/correcting had non-zero exit status")

    return 0
