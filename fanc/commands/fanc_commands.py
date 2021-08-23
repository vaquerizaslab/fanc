import argparse
import logging
import os
import os.path
import textwrap
import shutil
import tempfile
import warnings

# configure logging
logger = logging.getLogger(__name__)


def fanc_parser():
    usage = "fanc <command> [options]\n\n"

    command_descriptions = dict()
    for name, function in globals().items():
        if name.endswith("_parser") and name != 'fanc_parser':
            parser = function()
            short_name = name[:-7].replace('_', '-')
            command_descriptions[short_name] = parser.description.split(".")[0]

    max_len = max([len(name) for name in command_descriptions.keys()]) + 4

    usage += "-- Matrix generation --\n"
    for name in ['auto', 'map', 'pairs', 'hic']:
        padding = ' ' * (max_len - len(name))
        usage += "{}{}{}\n".format(name, padding, command_descriptions.pop(name))

    usage += "\n-- Matrix analysis --\n"
    for name in ['cis-trans', 'expected', 'pca', 'compartments',
                 'insulation', 'directionality', 'boundaries',
                 'compare', 'loops', 'aggregate']:
        padding = ' ' * (max_len - len(name))
        usage += "{}{}{}\n".format(name, padding, command_descriptions.pop(name))

    usage += "\n-- Other helpers --\n"
    for name in command_descriptions.keys():
        padding = ' ' * (max_len - len(name))
        usage += "{}{}{}\n".format(name, padding, command_descriptions.get(name))

    parser = argparse.ArgumentParser(
        description="fanc processing tool for Hi-C data",
        usage=textwrap.dedent(usage)
    )

    parser.add_argument(
        '-V', '--version', dest='print_version',
        action='store_true',
        help='Print version information'
    )
    parser.set_defaults(print_version=False)

    parser.add_argument(
        '--verbose', '-v', dest='verbosity',
        action='count',
        default=0,
        help='Set verbosity level: Can be chained like '
             '"-vvv" to increase verbosity. Default is to show '
             'errors, warnings, and info messages (same as "-vv"). '
             '"-v" shows only errors and warnings,  "-vvv" shows '
             'errors, warnings, info, and debug messages.'
    )

    parser.add_argument(
        '-s', '--silent', dest='silent',
        action='store_true',
        help='Do not print log messages to command line.'
    )
    parser.set_defaults(silent=False)

    parser.add_argument(
        '-l', '--log-file', dest='log_file',
        help='Path to file in which to save log.'
    )

    parser.add_argument(
        '-m', '--email', dest='email_to_address',
        help='Email address for fanc command summary.'
    )

    parser.add_argument(
        '--smtp-server', dest='smtp_server',
        help='SMTP server in the form smtp.server.com[:port].'
    )

    parser.add_argument(
        '--smtp-username', dest='smtp_username',
        help='SMTP username.'
    )

    parser.add_argument(
        '--smtp-password', dest='smtp_password',
        help='SMTP password.'
    )

    parser.add_argument(
        '--smtp-sender-address', dest='email_from_address',
        help='SMTP sender email address.'
    )

    parser.add_argument('command', nargs='?', help='Subcommand to run')

    return parser


def auto_parser():
    import fanc.commands.auto
    return fanc.commands.auto.auto_parser()


def auto(argv, **kwargs):
    import fanc.commands.auto
    return fanc.commands.auto.auto(argv, **kwargs)


def map_parser():
    parser = argparse.ArgumentParser(
        prog="fanc map",
        description='Map reads in a FASTQ file to a reference genome.'
    )

    parser.add_argument(
        'input',
        nargs='+',
        help='File name of the input FASTQ file (or gzipped FASTQ)'
    )

    parser.add_argument(
        'index',
        help='Bowtie 2 or BWA genome index base. Index type will be '
             'determined automatically.'
    )

    parser.add_argument(
        'output',
        help='Output file or folder. '
             'When providing multiple input files, this must be the '
             'path to an output folder.'
    )

    parser.add_argument(
        '-m', '--min-size', dest='min_size',
        type=int,
        default=25,
        help='Minimum length of read before extension. '
             'Default %(default)d.'
    )

    parser.add_argument(
        '-s', '--step-size', dest='step_size',
        type=int,
        default=10,
        help='Number of base pairs to extend at each round of mapping. '
             'Default is %(default)d.'
    )

    parser.add_argument(
        '--trim-front', dest='trim_front',
        action='store_true',
        default=False,
        help='Trim reads from front instead of back.'
    )

    parser.add_argument(
        '-t', '--threads', dest='threads',
        type=int,
        default=1,
        help='Number of threads used for mapping. Default: %(default)d'
    )

    parser.add_argument(
        '-q', '--quality', dest='quality',
        type=int,
        help='Mapping quality cutoff. '
             'Alignments with a quality score lower than this will be '
             'sent to another mapping iteration. '
             'Default: 3 (BWA), 30 (Bowtie2)'
    )

    parser.add_argument(
        '-r', '--restriction-enzyme', dest='restriction_enzyme',
        help='Name (case sensitive) of restriction enzyme used in Hi-C experiment. '
             'Will be used to split reads by predicted ligation junction before mapping. '
             'You can omit this if you do not want to split your reads by ligation junction. '
             'Restriction names can be any supported by Biopython, which obtains data '
             'from REBASE (http://rebase.neb.com/rebase/rebase.html). '
             'For restriction enzyme cocktails, separate enzyme names with ","'
    )

    parser.add_argument(
        '-k', '--max-alignments', dest='max_alignments',
        type=int,
        help='Maximum number of alignments per read to be reported.'
    )

    parser.add_argument(
        '-a', '--all-alignments', dest='all_alignments',
        action='store_true',
        default=False,
        help='Report all valid alignments of a read '
             'Warning: very slow!.'
    )

    parser.add_argument(
        '-b', '--batch-size', dest='batch_size',
        type=int,
        default=100000,
        help='Number of reads processed (mapped and merged) in one go per worker. '
             'The default %(default)d works well for large indexes (e.g. human, mouse). '
             'Smaller indexes (e.g. yeast) will finish individual bowtie2 processes '
             'very quickly - set this number higher to spawn new processes '
             'less frequently. '
    )

    parser.add_argument(
        '--fanc-parallel', dest='mapper_parallel',
        action='store_false',
        default=True,
        help='Use FAN-C parallelisation, which launches multiple mapper jobs. '
             'This may be faster in some cases than relying '
             'on the internal paralellisation of the mapper, '
             'but has potentially high disk I/O and memory usage.'
    )

    parser.add_argument(
        '--split-fastq', dest='split_fastq',
        action='store_true',
        default=False,
        help='Split FASTQ file into 10M chunks before mapping. '
             'Easier on tmp partitions.'
    )

    parser.add_argument(
        '--memory-map', dest='memory_map',
        action='store_true',
        default=False,
        help='Map Bowtie2 index to memory. ' 
             'Only enable if your system has enough memory '
             'to hold the entire Bowtie2 index.'
    )

    parser.add_argument(
        '--no-iterative', dest='iterative',
        action='store_false',
        default=True,
        help='Do not use iterative mapping strategy. '
             '(much faster, less sensitive).'
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        default=False,
        help='Copy original file to temporary directory.'
             'Reduces network I/O.'
    )

    return parser


def map(argv, **kwargs):
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
    mapper_parallel = args.mapper_parallel
    memory_map = args.memory_map
    iterative = args.iterative
    restriction_enzyme = args.restriction_enzyme
    max_alignments = args.max_alignments
    all_alignments = args.all_alignments
    tmp = args.tmp

    if mapper_parallel:
        threads, mapper_threads = 1, args.threads
    else:
        threads, mapper_threads = args.threads, 1

    if restriction_enzyme is not None:
        restriction_enzyme = restriction_enzyme.split(",")

    import fanc.map as map
    from fanc.tools.general import mkdir
    from genomic_regions.files import create_temporary_copy
    from fanc.tools.files import random_name
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

    if index_path.endswith('.'):
        index_path = index_path[:-1]

    mapper_type = 'bwa'
    for ending in ('amb', 'ann', 'bwt', 'pac', 'sa'):
        file_name = index_path + '.{}'.format(ending)
        if not os.path.isfile(file_name):
            mapper_type = None
            break

    if mapper_type is None:
        mapper_type = 'bowtie2'
        for i in range(1, 5):
            if not os.path.exists(index_path + '.{}.bt2'.format(i)):
                mapper_type = None
                break
        for i in range(1, 3):
            if not os.path.exists(index_path + '.rev.{}.bt2'.format(i)):
                mapper_type = None
                break

    if mapper_type is None:
        raise RuntimeError("Cannot detect mapper type from index (supported are Bowtie2 and BWA)")
    elif min_quality is None:
        if mapper_type == 'bwa':
            min_quality = 3
        elif mapper_type == 'bowtie2':
            min_quality = 30

    index_dir = None
    try:
        if tmp:
            tmp = False

            index_dir = tempfile.mkdtemp()
            index_base = os.path.basename(index_path)
            if mapper_type == 'bowtie2':
                for file_name in glob.glob(index_path + '*.bt2'):
                    shutil.copy(file_name, index_dir)
            elif mapper_type == 'bwa':
                index_base = random_name()
                for ending in ('amb', 'ann', 'bwt', 'pac', 'sa'):
                    file_name = index_path + '.{}'.format(ending)
                    shutil.copy(file_name, os.path.join(index_dir, '{}.{}'.format(index_base, ending)))

            index_path = os.path.join(index_dir, index_base)
            logger.debug('Index path: {}'.format(index_path))
            tmp = True

        mapper = None
        if mapper_type == 'bowtie2':
            if iterative:
                mapper = map.Bowtie2Mapper(index_path, min_quality=min_quality,
                                           additional_arguments=additional_arguments,
                                           threads=mapper_threads)
            else:
                mapper = map.SimpleBowtie2Mapper(index_path, additional_arguments=additional_arguments,
                                                 threads=mapper_threads)
        elif mapper_type == 'bwa':
            if iterative:
                mapper = map.BwaMapper(index_path, min_quality=min_quality,
                                       threads=mapper_threads, memory_map=memory_map)
            else:
                mapper = map.SimpleBwaMapper(index_path,
                                             threads=mapper_threads, memory_map=memory_map)

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
                                                               prefix='fanc_', delete=False)
                        tmp_file_name = tmp_file.name
                        output_file = tmp_file_name
                        tmp = True

                    logger.info("Starting mapping for {}".format(input_file))
                    map.iterative_mapping(input_file, output_file, mapper, threads=threads,
                                          min_size=min_size, step_size=step_size, batch_size=batch_size,
                                          trim_front=trim_front, restriction_enzyme=restriction_enzyme)
                    logger.debug("Mapping complete")
                finally:
                    if tmp:
                        logger.debug("Removing tmp files")
                        os.remove(input_file)
                        shutil.copy(output_file, original_output_file)
                        os.remove(output_file)
                    logger.debug("Closing mapper")
                    mapper.close()
            else:
                from fanc.tools.files import split_fastq, merge_sam, gzip_splitext

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

                        split_command = ['fanc', 'map', split_file, index_path, split_bam_file,
                                         '-m', str(min_size), '-s', str(step_size), '-t', str(threads),
                                         '-q', str(args.quality), '-b', str(batch_size)]
                        if not iterative:
                            split_command += ['--no-iterative']
                        if tmp:
                            split_command += ['-tmp']
                        if not mapper_parallel:
                            split_command += ['--fanc-parallel']
                        if memory_map:
                            split_command += ['--memory-map']
                        if trim_front:
                            split_command += ['--trim-front']
                        if restriction_enzyme is not None:
                            split_command += ['--restriction-enzyme', ",".join(restriction_enzyme)]

                        rt = subprocess.call(split_command)
                        split_fastq_results.append(rt)

                    for rt in split_fastq_results:
                        if rt != 0:
                            raise RuntimeError("Mapping had non-zero exit status")

                    logger.info("Merging BAM files into {}".format(output_file))
                    merge_sam(split_bam_files, output_file, tmp=tmp)
                finally:
                    shutil.rmtree(split_tmpdir, ignore_errors=True)
    finally:
        if tmp:
            logger.debug("Removing index dir")
            shutil.rmtree(index_dir, ignore_errors=True)


def fragments_parser():
    parser = argparse.ArgumentParser(
        prog="fanc fragments",
        description='In-silico genome digestion'
    )

    parser.add_argument(
        'input',
        help="Path to genome file (FASTA, folder with FASTA, hdf5 file), "
             "which will be used in conjunction with the type of restriction enzyme to "
             "calculate fragments directly."
    )

    parser.add_argument(
        're_or_bin_size',
        help="Restriction enzyme name or bin size to divide genome into fragments. "
             "Restriction names can be any supported by Biopython, which obtains data "
             "from REBASE (http://rebase.neb.com/rebase/rebase.html). "
             "Use commas to separate multiple restriction enzymes, e.g. 'HindIII,MboI'"
    )

    parser.add_argument(
        'output',
        help='Output file with restriction fragments in BED format.'
    )

    parser.add_argument(
        '-c', '--chromosomes', dest='chromosomes',
        help='Comma-separated list of chromosomes to include in fragments BED file. '
             'Other chromosomes will be excluded. The order of chromosomes will '
             'be as stated in the list.'
    )

    return parser


def fragments(argv, **kwargs):
    parser = fragments_parser()
    args = parser.parse_args(argv[2:])

    import os

    genome_file = os.path.expanduser(args.input)
    re_or_bin_size = args.re_or_bin_size.split(",")
    output_file = os.path.expanduser(args.output)
    chromosomes = args.chromosomes

    from fanc.regions import genome_regions
    regions = genome_regions(genome_file, re_or_bin_size)

    from fanc.tools.files import write_bed
    import warnings
    from collections import defaultdict

    if chromosomes is not None:
        chromosomes = chromosomes.split(",")
        chromosome_regions = defaultdict(list)
        for region in regions:
            chromosome_regions[region.chromosome].append(region)

        regions = []
        for chromosome in chromosomes:
            if chromosome not in chromosome_regions:
                raise ValueError("Chromosome {} is not found in genome. Please check your chromosome list!")
            regions += chromosome_regions[chromosome]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        write_bed(output_file, regions)


def sort_sam_parser():
    parser = argparse.ArgumentParser(
        prog="fanc sort_sam",
        description="Convenience function to sort a SAM file by name. "
                    "Exactly the same as 'samtools sort -n', but potentially"
                    "faster if sambamba is available."
    )

    parser.add_argument(
        'sam',
        help='Input SAM/BAM'
    )

    parser.add_argument(
        'output',
        nargs='?',
        help='Output SAM/BAM. If not provided, will replace input ' 
             'file with sorted version after sorting.'
    )

    parser.add_argument(
        '-t', '--threads', dest='threads',
        default=1,
        type=int,
        help='Number of sorting threads (only when sambamba is available). '
             'Default: %(default)d'
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        default=False,
        help='Work in temporary directory'
    )

    return parser


def sort_sam(argv, **kwargs):
    parser = sort_sam_parser()
    args = parser.parse_args(argv[2:])

    sam_file = os.path.expanduser(args.sam)
    output_file = None if args.output is None else os.path.expanduser(args.output)
    threads = args.threads
    tmp = args.tmp

    from genomic_regions.files import create_temporary_copy
    from fanc.tools.files import sort_natural_sam
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
                with tempfile.NamedTemporaryFile(delete=False, prefix='fanc_', suffix=extension) as f:
                    output_file = f.name
            tmp = True
            logger.info("Working in tmp: {}, ".format(sam_file, output_file))

        output_file = sort_natural_sam(sam_file, output_file, threads=threads)
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


def pairs_parser():
    parser = argparse.ArgumentParser(
        prog="fanc pairs",
        description='Process and filter read pairs'
    )

    parser.add_argument(
        'input',
        nargs='+',
        help='IMPORTANT: The last positional argument will be '
             'the output file, unless only a single Pairs object '
             'is provided. In that case, filtering and correcting ' 
             'will be done in place. '
             'Possible inputs are: two SAM/BAM files (paired-end reads, ' 
             'sorted by read name using "samtools sort -n" or equivalent) and an output file; ' 
             'a HiC-Pro pairs file (format: ' 
             'name<tab>chr1<tab>pos1<tab>strand1<tab>chr2<tab>pos2<tab>strand2) ' 
             'and an output file; a pairs file in 4D Nucleome format '
             '(https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md) '
             'and an output file, ' 
             'or an existing fanc Pairs object. ' 
             'In case of SAM/BAM, HiC-Pro, or 4D Nucleome you must also provide the ' 
             '--genome argument, and if --genome is not a file with '
             'restriction fragments (or Hi-C bins), you must also provide the ' 
             '--restriction-enzyme argument.'
    )

    parser.add_argument(
        '-g', '--genome', dest='genome',
        help='Path to region-based file (BED, GFF, ...) containing the non-overlapping '
             'regions to be used for Hi-C binning. Typically restriction-enzyme fragments. '
             'Alternatively: Path to genome file (FASTA, folder with FASTA, HDF5 file), '
             'which will be used in conjunction with the type of restriction enzyme to ' 
             'calculate fragments directly.'
    )

    parser.add_argument(
        '-r', '--restriction_enzyme',
        help='Name of the restriction enzyme used in the '
             'experiment, e.g. HindIII, or MboI. Case-sensitive, '
             'only necessary when --genome is provided as FASTA. '
             'Restriction names can be any supported by Biopython, which obtains data '
             'from REBASE (http://rebase.neb.com/rebase/rebase.html). '
             'Separate multiple restriction enzymes with ","'
    )

    parser.add_argument(
        '-m', '--filter-unmappable', dest='unmappable',
        action='store_true',
        default=False,
        help='Filter read pairs where one or both halves are unmappable. ' 
             'Only applies to SAM/BAM input!'
    )

    parser.add_argument(
        '-u', '--filter-multimapping', dest='unique',
        action='store_true',
        default=False,
        help='Filter reads that map multiple times. If the other '
              'mapping locations have a lower score than the best one, ' 
              'the best read is kept. '
              'Only applies to SAM/BAM input!'
    )

    parser.add_argument(
        '-us', '--filter-multimapping-strict', dest='unique_strict',
        action='store_true',
        default=False,
        help='Strictly filter reads that map multiple times. ' 
             'Only applies to SAM/BAM input!'
    )

    parser.add_argument(
        '-q', '--filter-quality', dest='quality',
        type=float,
        help='Cutoff for the minimum mapping quality of a read. '
             'For numbers larger than 1, will filter on MAPQ. '
             'If a number between 0 and 1 is provided, will filter '
             'on the AS tag instead of mapping quality (only BWA). '
             'The quality cutoff is then interpreted as the '
             'fraction of bases that have to be matched for any '
             'given read. Only applies to SAM/BAM input! '
             'Default: no mapping quality filter.'
    )

    parser.add_argument(
        '-c', '--filter-contaminant', dest='contaminant',
        help='Filter contaminating reads from other organism. '
             'Path to mapped SAM/BAM file. ' 
             'Will filter out reads with the same name. ' 
             'Only applies to SAM/BAM input! '
             'Default: no contaminant filter'
    )

    parser.add_argument(
        '-i', '--filter-inward', dest='inward',
        type=int,
        help='Minimum distance for inward-facing read pairs. '
             'Default: no inward ligation error filter'
    )

    parser.add_argument(
        '-o', '--filter-outward', dest='outward',
        type=int,
        help='Minimum distance for outward-facing read pairs. '
             'Default: no outward ligation error filter'
    )

    parser.add_argument(
        '--filter-ligation-auto', dest='filter_le_auto',
        action='store_true',
        default=False,
        help='Auto-guess settings for inward/outward read pair filters. '
             'Overrides --filter-outward and --filter-inward if set. This '
             'is highly experimental and known to overshoot in some cases. '
             'It is generally recommended to specify cutoffs manually.'
    )

    parser.add_argument(
        '-d', '--filter-re-distance', dest='redist',
        type=int,
        help='Maximum distance for a read to the nearest restriction site. '
             'Default: no RE distance filter'
    )

    parser.add_argument(
        '-l', '--filter-self-ligations', dest='self_ligated',
        action='store_true',
        default=False,
        help='Remove read pairs representing self-ligated fragments.'
             'Default: no self-ligation filter.'
    )

    parser.add_argument(
        '-p', '--filter-pcr-duplicates', dest='dup_thresh',
        type=int,
        help='If specified, filter read pairs for PCR duplicates. Parameter determines '
             'distance between alignment starts below which they are considered starting '
             'at same position. Sensible values are between 1 and 5. Default: '
             'no PCR duplicates filter'
    )

    parser.add_argument(
        '-s', '--statistics', dest='stats',
        help='Path for saving filter statistics'
    )

    parser.add_argument(
        '--reset-filters', dest='reset_filters',
        action='store_true',
        default=False,
        help='Remove all filters from the ReadPairs object.'
    )

    parser.add_argument(
        '--statistics-plot', dest='stats_plot',
        help='Path for saving filter statistics plot (PDF)'
    )

    parser.add_argument(
        '--re-dist-plot', dest='re_dist_plot',
        help='Plot the distribution of restriction site distances '
             'of all read pairs (sum left and right read).'
    )

    parser.add_argument(
        '--ligation-error-plot', dest='ligation_error_plot',
        help='Plot the relative orientation of read pairs mapped '
             'to the reference genome as a fraction of reads oriented '
             'in the same direction. Allows the identification of '
             'ligation errors as a function of genomic distance.'
    )

    parser.add_argument(
        '-t', '--threads', dest='threads',
        type=int,
        default=1,
        help='Number of threads to use for extracting fragment information. '
             'Default: %(default)d'
    )

    parser.add_argument(
        '-b', '--batch-size', dest='batch_size',
        type=int,
        default=1000000,
        help='Batch size for read pairs to be submitted to individual processes. '
             'Default: %(default)d'
    )

    parser.add_argument(
        '-S', '--no-check-sorted', dest='check_sorted',
        action='store_false',
        default=True,
        help='Assume SAM files are sorted and do not '
             'check if that is actually the case'
    )

    parser.add_argument(
        '-f', '--force-overwrite', dest='force_overwrite',
        action='store_true',
        default=False,
        help='If the specified output file exists, it will be ' 
             'overwritten without warning.'
    )

    parser.add_argument(
        '--bwa', dest='bwa',
        action='store_true',
        default=False,
        help='Use filters appropriate for BWA and not Bowtie2. '
             'This will typically be identified automatically ' 
             'from the SAM/BAM header. Set this flag if you are ' 
             'having problems during filtering (typically 0 reads ' 
             'pass the filtering threshold). '
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        default=False,
        help='Work in temporary directory'
    )

    return parser


def pairs(argv, **kwargs):
    parser = pairs_parser()
    args = parser.parse_args(argv[2:])

    import os

    input_files = [os.path.expanduser(file_name) for file_name in args.input]
    genome_file = os.path.expanduser(args.genome) if args.genome is not None else None
    restriction_enzyme = args.restriction_enzyme
    filter_unmappable = args.unmappable
    filter_unique = args.unique
    filter_unique_strict = args.unique_strict
    filter_quality = args.quality
    filter_contaminant = args.contaminant
    filter_inward = args.inward
    filter_outward = args.outward
    filter_le_auto = args.filter_le_auto
    filter_re_distance = args.redist
    filter_self_ligations = args.self_ligated
    filter_pcr_duplicates = args.dup_thresh
    statistics_file = os.path.expanduser(args.stats) if args.stats is not None else None
    statistics_plot_file = os.path.expanduser(args.stats_plot) if args.stats_plot is not None else None
    re_dist_plot_file = os.path.expanduser(args.re_dist_plot) if args.re_dist_plot is not None else None
    ligation_error_plot_file = os.path.expanduser(args.ligation_error_plot) \
        if args.ligation_error_plot is not None else None
    threads = args.threads
    batch_size = args.batch_size
    check_sam_sorted = args.check_sorted
    force_overwrite = args.force_overwrite
    reset_filters = args.reset_filters
    tmp = args.tmp

    from genomic_regions.files import create_temporary_copy, create_temporary_output

    if restriction_enzyme is not None:
        restriction_enzyme = restriction_enzyme.split(",")

    tmp_input_files = []
    original_pairs_file = None
    pairs_file = None
    pairs = None
    try:
        regions = None
        if 2 <= len(input_files) <= 3:
            if not force_overwrite and os.path.exists(input_files[-1]):
                parser.error("Output file {} exists! Use -f to force "
                             "overwriting it!".format(input_files[-1]))

            if genome_file is None:
                parser.error("Must provide genome file (-g) when loading reads or pairs!")

            logger.info("Getting genome regions (fragments or bins)")
            from fanc.regions import genome_regions
            try:
                if tmp:
                    tmp = False
                    genome_file = create_temporary_copy(os.path.expanduser(genome_file))
                    tmp_input_files.append(genome_file)
                    tmp = True

                regions = genome_regions(genome_file, restriction_enzyme=restriction_enzyme)
            except ValueError:
                if restriction_enzyme is None:
                    parser.error("Must provide --restriction-enzyme when --genome is "
                                 "not a list of fragments or genomic bins!")
                else:
                    raise

        import fanc

        if len(input_files) == 3:
            logger.info("Three arguments detected, assuming SAM/BAM input.")
            sam1_file = os.path.expanduser(input_files[0])
            sam2_file = os.path.expanduser(input_files[1])
            pairs_file = os.path.expanduser(input_files[2])

            if tmp:
                tmp = False
                sam1_file = create_temporary_copy(sam1_file)
                sam2_file = create_temporary_copy(sam2_file)
                tmp_input_files += [sam1_file, sam2_file]
                original_pairs_file = pairs_file
                pairs_file = create_temporary_output(pairs_file)
                tmp = True

            from fanc.tools.general import get_sam_mapper
            bwa = get_sam_mapper(sam1_file) == 'bwa' or args.bwa
            logger.info("Using filters appropriate for {}.".format('BWA' if bwa else 'Bowtie2'))

            from fanc.pairs import BwaMemQualityFilter, BwaMemUniquenessFilter, \
                UniquenessFilter, ContaminantFilter, QualityFilter, UnmappedFilter
            from fanc.general import Mask

            read_filters = []
            if filter_unmappable:
                f = UnmappedFilter(mask=Mask(ix=0, name='unmappable'))
                read_filters.append(f)

            if filter_unique or filter_unique_strict:
                if bwa:
                    f = BwaMemUniquenessFilter(strict=filter_unique_strict,
                                               mask=Mask(ix=1, name='multi-mapping'))
                else:
                    f = UniquenessFilter(strict=filter_unique_strict,
                                         mask=Mask(ix=1, name='multi-mapping'))
                read_filters.append(f)

            if filter_quality is not None:
                if 0 < filter_quality < 1:
                    f = BwaMemQualityFilter(filter_quality, mask=Mask(ix=2, name='alignment score'))
                else:
                    f = QualityFilter(int(filter_quality), mask=Mask(ix=2, name='MAPQ'))
                read_filters.append(f)

            if filter_contaminant is not None:
                f = ContaminantFilter(filter_contaminant, mask=Mask(ix=3, name='contaminant'))
                read_filters.append(f)

            from fanc.pairs import generate_pairs_split as generate_pairs
            pairs = generate_pairs(sam1_file, sam2_file, regions,
                                   restriction_enzyme=restriction_enzyme,
                                   read_filters=read_filters, output_file=pairs_file,
                                   check_sorted=check_sam_sorted, threads=threads,
                                   batch_size=batch_size)
            pairs.close()
        elif len(input_files) == 2:
            logger.info("Two arguments detected, assuming HiC-Pro or 4D Nucleome input.")
            from fanc.pairs import HicProPairGenerator, FourDNucleomePairGenerator, ReadPairs
            from fanc.regions import genome_regions
            import genomic_regions as gr

            input_file = os.path.expanduser(input_files[0])
            pairs_file = os.path.expanduser(input_files[1])

            if tmp:
                tmp = False
                from genomic_regions.files import create_temporary_copy, create_temporary_output
                input_file = create_temporary_copy(input_file)
                tmp_input_files = [input_file]
                original_pairs_file = pairs_file
                pairs_file = create_temporary_output(pairs_file)
                tmp = True

            try:
                logger.debug("Trying 4D nucleome format...")
                sb = FourDNucleomePairGenerator(input_file)
            except ValueError:
                logger.debug("Trying HiC-Pro format...")
                sb = HicProPairGenerator(input_file)

            pairs = ReadPairs(file_name=pairs_file, mode='w')

            if isinstance(regions, gr.RegionBased):
                pairs.add_regions(regions.regions, preserve_attributes=False)
            else:
                pairs.add_regions(regions, preserve_attributes=False)
            pairs.add_read_pairs(sb, threads=threads, batch_size=batch_size)
            pairs.close()
        elif len(input_files) == 1:
            logger.info("One argument received, assuming existing Pairs object.")
            pairs_file = os.path.expanduser(input_files[0])
            if tmp:
                tmp = False
                from genomic_regions.files import create_temporary_copy
                original_pairs_file = pairs_file
                pairs_file = create_temporary_copy(pairs_file)
                tmp = True
        else:
            parser.error("Number of input arguments ({}) cannot be parsed. "
                         "See help for details.".format(len(input_files)))

        if reset_filters:
            logger.info("Resetting all filters")
            pairs = fanc.load(pairs_file, mode='a')
            pairs.reset_filters()
            pairs.close()

        if (filter_le_auto or filter_inward or filter_outward or filter_re_distance or
                filter_self_ligations or filter_pcr_duplicates):
            pairs = fanc.load(pairs_file, mode='a')

            logger.debug("Preparing read pair filters")
            if filter_le_auto:
                logger.info("Filtering inward- and outward-facing reads using automatically"
                            "determined thresholds.")
                pairs.filter_ligation_products(queue=True)
            else:
                if filter_inward:
                    logger.info("Filtering inward-facing reads at %dbp" % filter_inward)
                    pairs.filter_inward(minimum_distance=filter_inward, queue=True)

                if filter_outward:
                    logger.info("Filtering outward-facing reads at %dbp" % filter_outward)
                    pairs.filter_outward(minimum_distance=filter_outward, queue=True)

            if filter_re_distance:
                logger.info("Filtering reads with RE distance > %dbp" % filter_re_distance)
                pairs.filter_re_dist(filter_re_distance, queue=True)

            if filter_self_ligations:
                logger.info("Filtering self-ligated read pairs")
                pairs.filter_self_ligated(queue=True)

            if filter_pcr_duplicates:
                logger.info("Filtering PCR duplicates, threshold <= %dbp" % filter_pcr_duplicates)
                pairs.filter_pcr_duplicates(threshold=filter_pcr_duplicates, queue=True)

            logger.info("Running filters...")
            pairs.run_queued_filters(log_progress=True)
            logger.info("Done.")

            pairs.close()

        if (statistics_file is not None or statistics_plot_file is not None or
                re_dist_plot_file is not None or ligation_error_plot_file is not None):
            import matplotlib
            matplotlib.use('agg')
            import matplotlib.pyplot as plt

            alignment_filter_keys = ['unmappable', 'multi-mapping', 'alignment score', 'MAPQ', 'contaminant']

            pairs = fanc.load(pairs_file, mode='r')
            if statistics_file is not None or statistics_plot_file is not None:
                statistics = pairs.filter_statistics()

                if statistics_file is not None:
                    with open(statistics_file, 'w') as o:
                        for name, value in statistics.items():
                            o.write("{}\t{}\n".format(name, value))

                if statistics_plot_file is not None:
                    logger.info("Saving statistics...")
                    from fanc.plotting.statistics import summary_statistics_plot
                    statistics_plot_file = os.path.expanduser(statistics_plot_file)
                    fig, axes = plt.subplots(1, 2)
                    summary_statistics_plot(statistics, exclude=alignment_filter_keys, ax=axes[1])
                    summary_statistics_plot(statistics, include=alignment_filter_keys, exclude=['total', 'valid'],
                                            ax=axes[0])
                    axes[0].set_ylabel("Number of filtered reads (R1+R2)")
                    axes[1].set_title("Read pair filter statistics")
                    axes[0].set_title("Alignment filter statistics")
                    fig.savefig(statistics_plot_file)
                    plt.close(fig)

            if re_dist_plot_file is not None:
                from fanc.plotting.statistics import restriction_site_distance_plot
                fig, ax = plt.subplots()
                restriction_site_distance_plot(pairs, ax=ax)
                fig.savefig(re_dist_plot_file)
                plt.close(fig)

            if ligation_error_plot_file is not None:
                from fanc.plotting.statistics import ligation_bias_plot
                fig, ax = plt.subplots()
                ligation_bias_plot(pairs, ax=ax)
                fig.savefig(ligation_error_plot_file)
                plt.close(fig)

            pairs.close()
    finally:
        if tmp:
            for tmp_input_file in tmp_input_files:
                os.remove(tmp_input_file)

            shutil.copy(pairs_file, original_pairs_file)
            os.remove(pairs_file)

    logger.info("All done.")


def hic_parser():
    parser = argparse.ArgumentParser(
        prog="fanc hic",
        description='Process, filter, and correct Hic files'
    )

    parser.add_argument(
        'input',
        nargs='+',
        help='IMPORTANT: The last positional argument will be '
             'the output file, unless only a single Hic object '
             'is provided. In that case, binning, filtering and '
             'correcting will be done in place. '
             'Input files. If these are FAN-C Pairs objects '
             '(see "fanc pairs"), they will be turned into '
             'Hic objects. Hic objects (also the ones converted '
             'from Pairs) will first be merged and the merged '
             'object will be binned, filtered and corrected as '
             'specified in the remaining parameters.'
    )

    parser.add_argument(
        '-b', '--bin-size', dest='bin_size',
        help='Bin size in base pairs. You can use human-readable formats,'
             'such as 10k, or 1mb. If omitted, the command will '
             'end after the merging step.'
    )

    parser.add_argument(
        '-l', '--filter-low-coverage', dest='filter_low_coverage',
        type=float,
        help='Filter bins with low coverage (lower than '
             'specified absolute number of contacts)'
    )

    parser.add_argument(
        '-r', '--filter-low-coverage-relative', dest='filter_low_coverage_relative',
        type=float,
        help='Filter bins using a relative low coverage threshold '
             '(lower than the specified fraction of the median contact count)'
    )

    parser.add_argument(
        '-a', '--low-coverage-auto', dest='filter_low_coverage_auto',
        action='store_true',
        default=False,
        help='Filter bins with "low coverage" (under ' 
             '10%% of median coverage for all non-zero bins)'
    )

    parser.add_argument(
        '-d', '--diagonal', dest='filter_diagonal',
        type=int,
        help='Filter bins along the diagonal up to this specified distance. '
             'Use 0 for only filtering the diagonal.'
    )

    parser.add_argument(
        '--marginals-plot', dest='marginals_plot',
        help='Plot Hi-C marginals to determine low coverage thresholds.'
    )

    parser.add_argument(
        '--reset-filters', dest='reset_filters',
        action='store_true',
        default=False,
        help='Remove all filters from the Hic object.'
    )

    parser.add_argument(
        '--downsample', dest='downsample',
        help="Downsample a binned Hi-C object before filtering and correcting. "
             "Sample size or reference Hi-C object. If sample size is < 1,"
             "will be interpreted as a fraction of valid pairs."
    )

    parser.add_argument(
        '--subset', dest='subset',
        help='Comma-separated list of regions that will be used in the output '
             'Hic object. All contacts between these regions '
             'will be in the output object. For example, '
             '"chr1,chr3" will result in a Hic object with '
             'all regions in chromosomes 1 and 3, plus all '
             'contacts within chromosome 1, all contacts within '
             'chromosome 3, and all contacts between chromosome '
             '1 and 3. "chr1" will only contain regions and contacts '
             'within chromosome 1.'
    )

    parser.add_argument(
        '-i', '--ice-correct', dest='ice',
        action='store_true',
        default=False,
        help='DEPRECATED. Use ICE iterative correction on the binned Hic matrix'
    )

    parser.add_argument(
        '-k', '--kr-correct', dest='kr',
        action='store_true',
        default=False,
        help='DEPRECATED. Use Knight-Ruiz matrix balancing to correct '
             'the binned Hic matrix'
    )

    parser.add_argument(
        '-n', '--normalise', dest='normalise',
        action='store_true',
        default=False,
        help='Normalise Hi-C matrix according to --norm-method'
    )

    parser.add_argument(
        '-m', '--norm-method', dest='norm_method',
        default='kr',
        help='Normalisation method used for -n. Options are: '
             'KR (default) = Knight-Ruiz matrix balancing '
             '(Fast, accurate, but memory-intensive normalisation); '
             'ICE = ICE matrix balancing (less accurate, but more memory-efficient); '
             'VC = vanilla coverage (a single round of ICE balancing);'
             'VC-SQRT = vanilla coverage square root (reduces overcorrection compared to VC)'
    )

    parser.add_argument(
        '-w', '--whole-matrix', dest='whole_matrix',
        action='store_true',
        default=False,
        help='Correct the whole matrix at once, rather than individual chromosomes.'
    )

    parser.add_argument(
        '-c', '--restore-coverage', dest='restore_coverage',
        action='store_true',
        default=False,
        help='Restore coverage to the original total number of reads. '
             'Otherwise matrix entries will be contact probabilities.'
    )

    parser.add_argument(
        '--only-inter', dest='only_inter',
        action='store_true',
        default=False,
        help="Calculate bias vector only on inter-chromosomal contacts. "
             "Ignores all intra-chromosomal contacts. "
             "Always uses whole-matrix balancing, i.e. implicitly sets -w"
    )

    parser.add_argument(
        '-s', '--statistics', dest='stats',
        help='Path for saving filter statistics'
    )

    parser.add_argument(
        '--statistics-plot', dest='stats_plot',
        help='Path for saving filter statistics plot (PDF)'
    )

    parser.add_argument(
        '--chromosomes', dest='chromosomes',
        nargs='+',
        help='Limit output Hic object to these chromosomes. '
             'Only available in conjunction with "-b" option.'
    )

    parser.add_argument(
        '-f', '--force-overwrite', dest='force_overwrite',
        action='store_true',
        default=False,
        help='If the specified output file exists, it will be ' 
             'overwritten without warning.'
    )

    parser.add_argument(
        '-t', '--threads', dest='threads',
        type=int,
        default=1,
        help="Number of threads (currently used for binning only)"
    )

    parser.add_argument(
        '--deepcopy', dest='deepcopy',
        action='store_true',
        default=False,
        help='Deep copy Hi-C file. Copies a Hi-C file to FAN-C format '
             'by duplicating individual bins, pixels, and bias information. '
             'Can be used to upgrade an existing FAN-C file with an older '
             'version or to convert Cooler or Juicer files to FAN-C format. '
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        default=False,
        help='Work in temporary directory'
    )

    return parser


def hic(argv, **kwargs):
    parser = hic_parser()
    args = parser.parse_args(argv[2:])

    import os
    from fanc.tools.general import str_to_int

    input_files = [os.path.expanduser(file_name) for file_name in args.input]
    bin_size = str_to_int(args.bin_size) if args.bin_size is not None else None
    filter_low_coverage = args.filter_low_coverage
    filter_low_coverage_relative = args.filter_low_coverage_relative
    filter_low_coverage_auto = args.filter_low_coverage_auto
    filter_diagonal = args.filter_diagonal
    downsample = args.downsample
    subset = args.subset

    do_norm = args.normalise
    norm_method = args.norm_method
    ice = args.ice
    kr = args.kr

    whole_matrix = args.whole_matrix
    restore_coverage = args.restore_coverage
    only_interchromosomal = args.only_inter
    statistics_file = os.path.expanduser(args.stats) if args.stats is not None else None
    statistics_plot_file = os.path.expanduser(args.stats_plot) if args.stats_plot is not None else None
    force_overwrite = args.force_overwrite
    reset_filters = args.reset_filters
    marginals_plot_file = args.marginals_plot
    limit_chromosomes = args.chromosomes
    threads = args.threads
    deepcopy = args.deepcopy
    tmp = args.tmp

    if kr or ice:
        warnings.warn("-k and -i have been deprecated in favor of -n and --norm-method. "
                      "-k and -i will be removed in a future version. "
                      "Please change your scripts accordingly!")

    if kr and ice:
        parser.error("The arguments --ice and --kr are mutually exclusive")

    if (kr or ice) and do_norm:
        parser.error("-n and -k (or -i) are incompatible. "
                     "Please use -n in combination with "
                     "--norm-method!")

    if kr:
        do_norm = True
        norm_method = 'kr'

    if ice:
        do_norm = True
        norm_method = 'ice'

    coverage_args = 0
    if filter_low_coverage_auto:
        coverage_args += 1
    if filter_low_coverage_relative is not None:
        coverage_args += 1
    if filter_low_coverage is not None:
        coverage_args += 1

    if coverage_args > 1:
        parser.error("The arguments -l, -r, and -a are mutually exclusive")

    if only_interchromosomal:
        whole_matrix = True

    import tempfile
    import fanc
    from genomic_regions.files import create_temporary_copy, create_temporary_output

    original_output_file = None
    output_file = None
    if len(input_files) > 1:
        output_file = input_files.pop()
        original_output_file = output_file

        if not force_overwrite and os.path.exists(output_file):
            parser.error("Output file {} exists! Use -f to force "
                         "overwriting it!".format(output_file))

        if tmp:
            tmp = False
            output_file = create_temporary_output(output_file)
            tmp = True

    hic_files = []
    pairs_files = []
    tmp_input_files = []
    try:
        for input_file in input_files:
            original_input_file = input_file
            if tmp:
                tmp = False
                if output_file is None:
                    original_output_file = input_file

                input_file = create_temporary_copy(input_file)
                tmp_input_files.append(input_file)
                tmp = True

            o = None
            try:
                o = fanc.load(input_file)
                if isinstance(o, fanc.Hic) or isinstance(o, fanc.hic.LegacyHic):
                    hic_files.append(input_file)
                elif isinstance(o, fanc.ReadPairs):
                    pairs_files.append(input_file)
                else:
                    from fanc.compatibility.cooler import CoolerHic
                    from fanc.compatibility.juicer import JuicerHic
                    if isinstance(o, JuicerHic) or isinstance(o, CoolerHic):
                        if deepcopy:
                            hic_files.append(input_file)
                        else:
                            parser.error("Can only work on Juicer or Cooler files when --deepcopy is enabled!")
                    else:
                        parser.error("File ({}) type {} not supported."
                                     "Provide Pairs or Hic files!".format(original_input_file, type(o)))
            finally:
                if o is not None and hasattr(o, 'close'):
                    o.close()

        if len(pairs_files) > 0:
            logger.info("Converting Pairs files")

            for pairs_file in pairs_files:
                with fanc.load(pairs_file) as pairs:
                    tmp_output_file = tempfile.NamedTemporaryFile(suffix='.hic', delete=False)
                    tmp_input_files.append(tmp_output_file.name)
                    logger.info("Converting {} to Hic".format(pairs_file))
                    hic = pairs.to_hic(file_name=tmp_output_file.name)
                    hic.close()
                    hic_files.append(tmp_output_file.name)

        if len(hic_files) > 1:
            if deepcopy:
                parser.error("Deep copy can only be performed on a single input Hi-C file!")
            logger.info("Merging {} Hic files into {}".format(len(hic_files), original_output_file))
            hics = []
            try:
                for hic_file in hic_files:
                    hics.append(fanc.load(hic_file))

                tmp_output_file = tempfile.NamedTemporaryFile(suffix='.hic', delete=False)
                merged_hic_file = tmp_output_file.name
                tmp_input_files.append(merged_hic_file)
                merged_hic = fanc.Hic.merge(hics, file_name=merged_hic_file, mode='w')
                merged_hic.close()
            finally:
                for hic in hics:
                    hic.close()
        else:
            if deepcopy:
                tmp_output_file = tempfile.NamedTemporaryFile(suffix='.hic', delete=False)
                merged_hic_file = tmp_output_file.name
                tmp_input_files.append(merged_hic_file)
                hic = fanc.load(hic_files[0])
                copy_hic = hic.deepcopy(target_class=fanc.Hic, file_name=merged_hic_file, mode='w')
                copy_hic.close()
            else:
                merged_hic_file = hic_files[0]

        if bin_size is not None:
            merged_hic = fanc.load(merged_hic_file)

            logger.info("Binning Hic file ({})".format(bin_size))
            if downsample is None and subset is None:
                output_binned_file = output_file
            else:
                f = tempfile.NamedTemporaryFile(delete=False, suffix='.hic')
                output_binned_file = f.name
                tmp_input_files.append(output_binned_file)
            binned_hic = merged_hic.bin(bin_size, file_name=output_binned_file,
                                        threads=threads, chromosomes=limit_chromosomes)
        else:
            if downsample is None and subset is None and output_file is not None:
                shutil.copy(merged_hic_file, output_file)
                merged_hic_file = output_file
            binned_hic = fanc.load(merged_hic_file, mode='a')

        if downsample is not None:
            downsampled_hic = binned_hic.downsample(downsample, file_name=output_file)
            binned_hic = downsampled_hic

        if subset is not None:
            subset_regions = subset.split(",")
            subset_hic = binned_hic.subset(*subset_regions, file_name=output_file)
            binned_hic = subset_hic

        if reset_filters:
            logger.info("Resetting all filters")
            binned_hic.reset_filters()

        filters = []
        if filter_low_coverage_auto:
            from fanc.hic import LowCoverageFilter
            logger.info("Filtering low-coverage bins at 10%%")
            mask = binned_hic.add_mask_description('low_coverage',
                                                   'Mask low coverage regions in the Hic matrix '
                                                   '(relative cutoff {:.1%}'.format(0.1))

            low_coverage_auto_filter = LowCoverageFilter(binned_hic, rel_cutoff=0.1,
                                                         cutoff=None, mask=mask)
            filters.append(low_coverage_auto_filter)

        if filter_low_coverage is not None or filter_low_coverage_relative is not None:
            from fanc.hic import LowCoverageFilter
            logger.info("Filtering low-coverage bins using absolute cutoff {:.4}, "
                        "relative cutoff {:.1%}".format(float(filter_low_coverage)
                                                        if filter_low_coverage else 0.,
                                                        float(filter_low_coverage_relative)
                                                        if filter_low_coverage_relative else 0.))
            mask = binned_hic.add_mask_description('low_coverage',
                                                   'Mask low coverage regions in the Hic matrix '
                                                   '(absolute cutoff {:.4}, '
                                                   'relative cutoff {:.1%}'.format(
                                                       float(filter_low_coverage)
                                                       if filter_low_coverage else 0.,
                                                       float(filter_low_coverage_relative)
                                                       if filter_low_coverage_relative else 0.)
                                                   )

            low_coverage_filter = LowCoverageFilter(binned_hic, rel_cutoff=filter_low_coverage_relative,
                                                    cutoff=filter_low_coverage, mask=mask)
            filters.append(low_coverage_filter)

        if filter_diagonal is not None:
            from fanc.hic import DiagonalFilter
            logger.info("Filtering diagonal at distance {}".format(filter_diagonal))
            mask = binned_hic.add_mask_description('diagonal',
                                                   'Mask the diagonal of the Hic matrix '
                                                   '(up to distance {})'.format(filter_diagonal))
            diagonal_filter = DiagonalFilter(binned_hic, distance=filter_diagonal, mask=mask)
            filters.append(diagonal_filter)

        if len(filters) > 0:
            logger.info("Running filters...")
            for f in filters:
                binned_hic.filter(f, queue=True)
            binned_hic.run_queued_filters(log_progress=True)
            logger.info("Done.")

        if statistics_file is not None or statistics_plot_file is not None or marginals_plot_file is not None:
            import matplotlib
            matplotlib.use('agg')
            import matplotlib.pyplot as plt
            if statistics_file is not None or statistics_plot_file is not None:
                statistics = binned_hic.filter_statistics()

                if statistics_file is not None:
                    with open(statistics_file, 'w') as o:
                        for name, value in statistics.items():
                            o.write("{}\t{}\n".format(name, value))

                if statistics_plot_file is not None:
                    logger.info("Saving statistics...")
                    from fanc.plotting.statistics import summary_statistics_plot
                    statistics_plot_file = os.path.expanduser(statistics_plot_file)
                    fig, ax = plt.subplots()
                    summary_statistics_plot(statistics, ax=ax)
                    ax.set_ylabel("Non-zero pixels / Positive contacts between region pairs")
                    fig.savefig(statistics_plot_file)
                    plt.close(fig)

            if marginals_plot_file is not None:
                from fanc.plotting.statistics import marginals_plot
                chromosomes = binned_hic.chromosomes()
                cols = min(len(chromosomes), 4)
                rows, remainder = divmod(len(chromosomes), cols)
                if remainder > 0:
                    rows += 1
                # fig, axes = plt.subplots(rows, cols)
                fig = plt.figure(figsize=(cols*2, rows*2))
                for i, chromosome in enumerate(chromosomes):
                    row, col = divmod(i, cols)
                    ax = plt.subplot2grid((rows, cols), (row, col))
                    marginals_plot(binned_hic, chromosome, lower=filter_low_coverage,
                                   rel_cutoff=filter_low_coverage_relative, ax=ax)
                    ax.set_title(chromosome)
                fig.savefig(marginals_plot_file)
                plt.close(fig)

        if do_norm:
            logger.info("Normalising binned Hic file")

            binned_hic.normalise(norm_method, whole_matrix=whole_matrix,
                                 intra_chromosomal=not only_interchromosomal,
                                 restore_coverage=restore_coverage)

        binned_hic.close()
    finally:
        if tmp and original_output_file is not None:
            if output_file is not None:
                shutil.copy(output_file, original_output_file)
                os.remove(output_file)
            elif len(tmp_input_files) == 1:
                shutil.copy(tmp_input_files[0], original_output_file)

        for tmp_input_file in tmp_input_files:
            os.remove(tmp_input_file)


def from_juicer_parser():
    parser = argparse.ArgumentParser(
        prog="fanc from_juicer",
        description='Import a Hi-C object from juicer (Aiden lab)'
    )

    parser.add_argument(
        'input',
        help='''Input .hic file, juicer format'''
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
        '--juicer-tools-jar', dest='juicer_tools_jar_path',
        help='Path to juicer jar. You can also specify this in fanc.conf'
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        default=False,
        help='''Work in temporary directory'''
    )
    return parser


def from_juicer(argv, **kwargs):
    parser = from_juicer_parser()

    args = parser.parse_args(argv[2:])
    input_file = os.path.expanduser(args.input)
    output_file = os.path.expanduser(args.output)
    genome_path = os.path.expanduser(args.genome)
    juicer_tools_jar_path = args.juicer_tools_jar_path
    resolution = args.resolution
    chromosomes = args.chromosomes
    inter_chromosomal = args.inter_chromosomal
    juicer_norm = args.juicer_norm
    tmp = args.tmp

    from genomic_regions.files import create_temporary_copy

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

        from fanc.compatibility.juicer import convert_juicer_to_hic
        hic = convert_juicer_to_hic(input_file, genome_path, resolution,
                                    juicer_tools_jar_path=juicer_tools_jar_path,
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


def from_txt_parser():
    parser = argparse.ArgumentParser(
        prog="fanc from-txt",
        description='Import a Hi-C object from a sparse matrix txt format'
    )

    parser.add_argument(
        'contacts',
        help='Contacts file in sparse matrix format, i.e. each row '
             'should contain <bin1><tab><bin2><tab><weight>.'
    )

    parser.add_argument(
        'regions',
        help='Path to file with genomic regions, for example in BED format: '
             '<chromosome><tab><start><tab><end>. The BED can optionally contain '
             'the bin index, as corresponding to the index used in the contacts file.'
    )

    parser.add_argument(
        'output',
        help='''Output Hic file.'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        default=False,
        help='''Work in temporary directory'''
    )
    return parser


def from_txt(argv, **kwargs):
    parser = from_txt_parser()

    args = parser.parse_args(argv[2:])
    contacts_file = os.path.expanduser(args.contacts)
    regions_file = os.path.expanduser(args.regions)
    output_file = os.path.expanduser(args.output)
    tmp = args.tmp

    from genomic_regions.files import create_temporary_copy

    original_contacts_file = contacts_file
    original_regions_file = regions_file
    original_output_file = output_file
    try:
        if tmp:  # copy file if required
            tmp = False  # to prevent deleting input file should this be interrupted at this point
            logger.info("Copying input files...")
            contacts_file = create_temporary_copy(original_contacts_file)
            regions_file = create_temporary_copy(original_regions_file)

            tmp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.hic')
            tmp_file.close()
            output_file = tmp_file.name
            logger.info("Temporary output file: %s" % output_file)
            tmp = True

        from fanc.compatibility.txt import load_regions, load_contacts
        import fanc
        regions, ix_converter = load_regions(regions_file)

        hic = fanc.Hic(file_name=output_file, mode='w')
        hic.add_regions(regions)

        for source, sink, weight in load_contacts(contacts_file, ix_converter=ix_converter):
            hic.add_edge_simple(source, sink, weight=weight)
        hic.flush()
        hic.close()
    finally:
        if tmp:
            logger.info("Removing tmp files...")
            os.remove(contacts_file)
            os.remove(regions_file)
            shutil.copy(output_file, original_output_file)
            os.remove(output_file)

    logger.info("All done.")


def from_cooler_parser():
    parser = argparse.ArgumentParser(
        prog="fanc from-cooler",
        description="Convert a Cooler (.cool or .mcool) file to FAN-C format."
    )

    parser.add_argument(
        'input',
        help="Input .cool or .mcool file. You must specify the "
             "resolution you want to export to FAN-C using the format "
             "/path/to/hic.mcool@<resolution>"
    )

    parser.add_argument(
        'output',
        help='Output FAN-C file.'
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        default=False,
        help='Work in temporary directory'
    )

    return parser


def from_cooler(argv, **kwargs):
    parser = from_cooler_parser()

    args = parser.parse_args(argv[2:])
    input_file = os.path.expanduser(args.input)
    output_file = os.path.expanduser(args.output)
    tmp = args.tmp

    import fanc
    from fanc.tools.general import str_to_int

    try:
        import cooler
    except ImportError:
        parser.error("Cannot import cooler. Install cooler with 'pip install cooler'.")

    tmp_files = []
    original_output_file = output_file
    try:
        if tmp:  # copy file if required
            tmp = False  # to prevent deleting input file should this be interrupted at this point
            tmp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.hic')
            tmp_file.close()
            output_file = tmp_file.name
            logger.info("Temporary output file: %s" % output_file)
            tmp_files = [output_file]
            tmp = True

        with fanc.load(input_file, mode='r', tmpdir=tmp) as cool:
            cool.deepcopy(fanc.Hic, file_name=output_file, mode='w')
    finally:
        if tmp:
            shutil.copy(output_file, original_output_file)
        for tmp_file in tmp_files:
            os.remove(tmp_file)
    logger.info("All done.")


def to_cooler_parser():
    parser = argparse.ArgumentParser(
        prog="fanc hic_to_cooler",
        description="""Convert a Hic file into cooler format.
                       See https://github.com/mirnylab/cooler for details.
                       If input Hi-C matrix is uncorrected, the uncorrected matrix is stored.
                       If it is corrected, the uncorrected matrix is stored and the bias vector.
                       Cooler always calculates corrected matrix on-the-fly from the uncorrected
                       matrix and the bias vector."""
    )

    parser.add_argument(
        'input',
        help='''Input .hic file, fanc format.'''
    )

    parser.add_argument(
        'output',
        help='''Output cooler file.'''
    )

    parser.add_argument(
        '-u', '--uncorrected', dest='norm',
        action='store_false',
        default=True,
        help='Output uncorrected matrix.'
    )

    parser.add_argument(
        '-t', '--threads', dest='threads',
        type=int,
        default=1,
        help='Number of threads used for balancing.'
    )

    parser.add_argument(
        '-M', '--no-multi', dest='multi',
        action='store_false',
        default=True,
        help='Do not produce a multi-resolution file. This is fast, '
             'as it does not "coarsen" the matrix at multiple resolutions, '
             'but the resulting file will be incompatible with HiGlass!'
    )

    parser.add_argument(
        '-r', '--resolutions', dest='resolutions',
        nargs='+',
        default=None,
        help='Resolutions in bp at which to "coarsen" the cooler matrix. '
             'Default resolutions are calculated as '
             'base-resolution * 2 ** z, where z is an increasing integer '
             'zoom level.'
    )

    parser.add_argument(
        '-S', '--no-natural-sort', dest='natural_sort',
        action='store_false',
        default=True,
        help='Do not sort regions by their natural chromosome order. '
             'When using this option, chromosomes will appear in the Cooler '
             'file in the order they are listed in the FAN-C file.'
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        default=False,
        help='Work in temporary directory'
    )

    return parser


def to_cooler(argv, **kwargs):
    parser = to_cooler_parser()

    args = parser.parse_args(argv[2:])
    input_file = os.path.expanduser(args.input)
    output_file = os.path.expanduser(args.output)
    norm = args.norm
    multi = args.multi
    threads = args.threads
    natural_sort = args.natural_sort
    tmp = args.tmp

    resolutions = args.resolutions

    import fanc
    from fanc.tools.general import str_to_int

    try:
        import cooler
    except ImportError:
        parser.error("Cannot import cooler. Install cooler with 'pip install cooler'.")
    from fanc.compatibility.cooler import to_cooler

    if resolutions is not None:
        resolutions = [str_to_int(r) for r in resolutions]

    tmp_files = []
    original_output_file = output_file
    try:
        if tmp:  # copy file if required
            tmp = False  # to prevent deleting input file should this be interrupted at this point
            tmp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.hic')
            tmp_file.close()
            output_file = tmp_file.name
            logger.info("Temporary output file: %s" % output_file)
            tmp_files = [output_file]
            tmp = True

        with fanc.load(input_file, mode='r', tmpdir=tmp) as hic:
            to_cooler(hic, output_file, balance=norm, multires=multi, resolutions=resolutions,
                      threads=threads, natural_order=natural_sort)
    finally:
        if tmp:
            shutil.copy(output_file, original_output_file)
        for tmp_file in tmp_files:
            os.remove(tmp_file)
    logger.info("All done.")


def to_juicer_parser():
    parser = argparse.ArgumentParser(
        prog="fanc hic_to_juicer",
        description="Convert a ReadPairs file to Juicer .hic format"
    )

    parser.add_argument(
        'input',
        nargs='+',
        help='Input .pairs file(s), FAN-C format.'
    )

    parser.add_argument(
        'output',
        help='Output Juicer file.'
    )

    parser.add_argument(
        '--juicer-tools-jar', dest='juicer_tools_jar_path',
        help='Path to juicer jar. You can also specify this in fanc.conf'
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        default=False,
        help='Work in temporary directory'
    )

    parser.add_argument(
        '-r', '--resolutions', dest='resolutions',
        nargs='+',
        help='Resolutions in bp at which to "zoom" the juicer matrix.'
    )

    return parser


def to_juicer(argv, **kwargs):
    parser = to_juicer_parser()

    args = parser.parse_args(argv[2:])
    input_files = [os.path.expanduser(f) for f in args.input]
    output_file = os.path.expanduser(args.output)
    juicer_tools_jar_path = args.juicer_tools_jar_path
    resolutions = args.resolutions
    tmp = args.tmp

    import fanc
    from fanc.tools.general import str_to_int
    from fanc.compatibility.juicer import to_juicer
    if resolutions is not None:
        resolutions = [str_to_int(r) for r in resolutions]

    pairs = [fanc.load(f, tmpdir=tmp) for f in input_files]
    try:
        to_juicer(pairs, output_file, resolutions=resolutions,
                  juicer_tools_jar_path=juicer_tools_jar_path,
                  tmp=tmp)
    finally:
        for p in pairs:
            p.close()
    logger.info("All done.")


def dump_parser():
    parser = argparse.ArgumentParser(
        prog="fanc dump",
        description='Dump Hic file to txt file(s).'
    )

    parser.add_argument(
        'hic',
        help='''Hic file'''
    )

    parser.add_argument(
        'matrix',
        nargs='?',
        help='Output file for matrix entries. If not provided, will write to stdout.'
    )

    parser.add_argument(
        'regions',
        nargs='?',
        help='Output file for Hic regions. If not provided, will write regions into matrix file.'
    )

    parser.add_argument(
        '-s', '--subset', dest='subset',
        help='Only output this matrix subset. Format: '
             '<chr>[:<start>-<end>][--<chr>[:<start><end>]], '
             'e.g.: '
             '"chr1--chr1" to extract only the chromosome 1 submatrix; '
             '"chr2:3400000-4200000" to extract contacts of this region '
             'on chromosome 2 to all other regions in the genome;'
    )

    parser.add_argument(
        '-S', '--no-sparse', dest='sparse',
        action='store_false',
        default=True,
        help='Store full, square matrix instead of sparse format.'
    )

    parser.add_argument(
        '--only-intra', dest='only_intra',
        default=False,
        action='store_true',
        help='Only dump intra-chromosomal data. ' 
             'Dumps everything by default.'
    )

    parser.add_argument(
        '-e', '--observed-expected', dest='oe',
        action='store_true',
        default=False,
        help='O/E transform matrix values.'
    )

    parser.add_argument(
        '-l', '--log2', dest='log2',
        action='store_true',
        default=False,
        help='Log2-transform matrix values. '
             'Useful for O/E matrices (-e option)'
    )

    parser.add_argument(
        '-u', '--uncorrected', dest='norm',
        action='store_false',
        default=True,
        help='Output uncorrected (not normalised) matrix values).'
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        default=False,
        help='Work in temporary directory'
    )

    return parser


def dump(argv, **kwargs):
    parser = dump_parser()
    args = parser.parse_args(argv[2:])
    hic_file = os.path.expanduser(args.hic)
    output_matrix = None if args.matrix is None else os.path.expanduser(args.matrix)
    output_regions = None if args.regions is None else os.path.expanduser(args.regions)
    subset_string = args.subset
    sparse = args.sparse
    oe = args.oe
    log2 = args.log2
    only_intra = args.only_intra
    norm = args.norm
    tmp = args.tmp

    import fanc
    import sys
    import numpy as np
    import signal
    # prevent BrokenPipeError message
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    col_subset_region = None
    row_subset_region = None
    if subset_string is not None:
        col_subset_string = None
        if '--' in subset_string:
            row_subset_string, col_subset_string = subset_string.split('--')
        else:
            row_subset_string = subset_string
        row_subset_region = fanc.GenomicRegion.from_string(row_subset_string)
        if col_subset_string is not None:
            col_subset_region = fanc.GenomicRegion.from_string(col_subset_string)

    logger.info("Extracting the following matrix region: {} vs {}".format(row_subset_region, col_subset_region))

    original_output_matrix = output_matrix
    try:
        if tmp:
            tmp = False
            tmp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.mat')
            tmp_file.close()
            output_matrix = tmp_file.name
            logger.info("Temporary output file: {}".format(output_matrix))
            tmp = True

        with fanc.load(hic_file, mode='r', tmpdir=tmp) as hic:
            ix_to_chromosome = dict()
            for i, region in enumerate(hic.regions):
                ix_to_chromosome[region.ix] = region.chromosome

            row_regions = [region for region in hic.regions(row_subset_region)]
            col_regions = [region for region in hic.regions(col_subset_region)]
            row_regions_dict = {region.ix: (region, i) for i, region in enumerate(row_regions)}
            col_regions_dict = {region.ix: (region, i) for i, region in enumerate(col_regions)}

            if oe:
                _, expected_intra, expected_inter = hic.expected_values(norm=norm)
            else:
                expected_intra, expected_inter = {}, 1

            if log2:
                transform = np.log2
            else:
                transform = float

            if not sparse:
                if output_matrix is None or output_regions is None:
                    raise ValueError("Cannot write matrix to stdout, must provide "
                                     "both matrix and regions file for output")
                m = hic.matrix(key=(row_subset_region, col_subset_region), oe=oe, log=log2, norm=norm)
                np.savetxt(output_matrix, m)
            else:
                if output_matrix is None:
                    o = sys.stdout
                else:
                    o = open(output_matrix, 'w')

                try:
                    for edge in hic.edges(key=(row_subset_region, col_subset_region), lazy=True, norm=norm):
                        source, i = row_regions_dict[edge.source]
                        sink, j = col_regions_dict[edge.sink]
                        weight = getattr(edge, hic._default_score_field)
                        source_ix, sink_ix = source.ix, sink.ix
                        source_chromosome, sink_chromosome = ix_to_chromosome[source_ix], ix_to_chromosome[sink_ix]

                        if source_chromosome == sink_chromosome:
                            if oe:
                                weight /= expected_intra[source_chromosome][abs(sink_ix - source_ix)]
                        else:
                            if only_intra:
                                continue
                            if oe:
                                weight /= expected_inter

                        if output_regions is None:
                            print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                                source.chromosome, source.start, source.end,
                                sink.chromosome, sink.start, sink.end,
                                transform(weight)
                            ), file=o)
                        else:
                            print("{}\t{}\t{}".format(
                                i, j, transform(weight)
                            ), file=o)
                except BrokenPipeError:
                    pass
                finally:
                    if output_matrix is not None:
                        o.close()
    finally:
        if tmp and original_output_matrix is not None:
            shutil.copy(output_matrix, original_output_matrix)
            os.remove(output_matrix)

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


def pca_parser():
    parser = argparse.ArgumentParser(
        prog="fanc pca",
        description='Do a PCA on multiple Hi-C objects'
    )

    parser.add_argument(
        'input',
        nargs='+',
        help='Input Hic files'
    )

    parser.add_argument(
        'output',
        help='Output file with PCA results.'
    )

    parser.add_argument(
        '-p', '--plot', dest='plot',
        help='Output plot. Path to PDF file where '
             'the PCA plot will be saved.'
    )

    parser.add_argument(
        '-s', '--sample-size', dest='sample_size',
        type=int,
        default=50000,
        help='Sample size for contacts to do the PCA on.'
             'Default: 50000'
    )

    parser.add_argument(
        '--inter-chromosomal', dest='inter',
        action='store_true',
        default=False,
        help='Also include inter-chromosomal contacts in PCA. '
             'By default, only intra-schromosomal contacts are '
             'considered.'
    )

    parser.add_argument(
        '-r', '--region', dest='region',
        help='Region to do PCA on. You could put a specific '
             'chromosome here, for example. By default, the whole '
             'genome is considered. Comma-separate multiple regions.'
    )

    parser.add_argument(
        '-e', '--expected-filter', dest='expected_filter',
        type=float,
        help='Minimum fold-enrichment over expected value. Contacts '
             'with a strength lower than <b>*E(d), where d is the '
             'distance between two loci and E is the corresponding expected '
             'contact strength, are filtered out before PCA. Default: no filter.'
    )

    parser.add_argument(
        '-b', '--background-filter', dest='background_filter',
        type=float,
        help='Minimum fold-enrichment over average inter-chromosomal '
             'contacts. Default: no filter.'
    )

    parser.add_argument(
        '--min-distance', dest='min_distance',
        help='Minimum distance of matrix bins in base pairs. '
             'You can use abbreviated formats such as 1mb, 10k, etc.'
    )

    parser.add_argument(
        '--max-distance', dest='max_distance',
        help='Maximum distance of matrix bins in base pairs. '
             'You can use abbreviated formats such as 1mb, 10k, etc.'
    )

    parser.add_argument(
        '-n', '--names', dest='names',
        nargs='+',
        help='''Sample names for plot labelling.'''
    )

    parser.add_argument(
        '--strategy', dest='strategy',
        default='variance',
        help='Mechanism to select pairs from Hi-C matrix. '
             'Default: variance. Possible values are: variance '
             '(select contacts with the largest variance in '
             'strength across samples first), fold-change '
             '(select pairs with the largest fold-change '
             'across samples first), and passthrough '
             '(no preference on pairs). '
    )

    parser.add_argument(
        '-c', '--colors', dest='colors',
        nargs='+',
        help='Colors for plotting.'
    )

    parser.add_argument(
        '-m', '--markers', dest='markers',
        nargs='+',
        help='Markers for plotting. Follows Matplotlib marker '
             'definitions: http://matplotlib.org/api/markers_api.html'
    )

    parser.add_argument(
        '-v', '--eigenvectors', dest='eigenvectors',
        nargs=2,
        default=[1, 2],
        help='Which eigenvectors to plot. Default: 1 2'
    )

    parser.add_argument(
        '-Z', '--no-zeros', dest='no_zeros',
        action='store_true',
        default=False,
        help='''Ignore pixels with no contacts in any sample.'''
    )

    parser.add_argument(
        '-S', '--no-scaling', dest='scaling',
        action='store_false',
        default=True,
        help='Do not scale input matrices to the '
             'same number of valid pairs. Use this only '
             'if you are sure matrices are directly '
             'comparable.'
    )

    parser.add_argument(
        '-f', '--force-overwrite', dest='force_overwrite',
        action='store_true',
        default=False,
        help='''If the specified output file exists, it will be 
                        overwritten without warning.'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        default=False,
        help='''Work in temporary directory'''
    )
    return parser


def pca(argv, **kwargs):
    parser = pca_parser()

    args = parser.parse_args(argv[2:])

    import os

    input_files = [os.path.expanduser(file_name) for file_name in args.input]
    output_file = os.path.expanduser(args.output)
    plot_file = os.path.expanduser(args.plot) if args.plot is not None else None
    argument_sample_names = args.names
    sample_size = args.sample_size
    inter_chromosomal = args.inter
    sub_region = args.region.split(",") if args.region is not None else None
    expected_filter = args.expected_filter
    background_filter = args.background_filter
    min_distance = args.min_distance
    max_distance = args.max_distance
    ignore_zeros = args.no_zeros
    strategy = args.strategy
    colors = args.colors
    markers = args.markers
    eigenvectors = args.eigenvectors
    scale = args.scaling
    force = args.force_overwrite
    tmp = args.tmp

    if argument_sample_names is None:
        sample_names = ['matrix_{}'.format(i) for i in range(len(input_files))]
    else:
        sample_names = argument_sample_names

    if not force and os.path.exists(output_file):
        parser.error("Output file {} already exists, use -f to overwrite!".format(output_file))

    if not force and plot_file is not None and os.path.exists(plot_file):
        parser.error("Plot file {} already exists, use -f to overwrite!".format(plot_file))

    if len(sample_names) != len(input_files):
        parser.error("Number of samples ({}) does not equal number "
                     "of sample names ({})".format(len(input_files), len(sample_names)))

    if colors is not None and len(colors) != len(input_files):
        parser.error("Number of samples ({}) does not equal number "
                     "of colors ({})".format(len(input_files), len(colors)))

    if markers is not None and len(markers) != len(input_files):
        parser.error("Number of samples ({}) does not equal number "
                     "of markers ({})".format(len(input_files), len(markers)))

    eigenvectors = [eigenvectors[0] - 1, eigenvectors[1] - 1]

    from fanc.tools.general import str_to_int
    min_distance = str_to_int(min_distance) if min_distance is not None else None
    max_distance = str_to_int(max_distance) if max_distance is not None else None

    import fanc
    matrices = []
    try:
        if len(input_files) > 1:
            matrices = [fanc.load(file_name, tmpdir=tmp) for file_name in input_files]

            from fanc.architecture.comparisons import hic_pca

            pca_info, pca_res = hic_pca(*matrices, sample_size=sample_size, region=sub_region,
                                        strategy=strategy, ignore_zeros=ignore_zeros,
                                        oe_enrichment=expected_filter,
                                        background_ligation=background_filter,
                                        min_distance=min_distance, max_distance=max_distance,
                                        min_libraries_above_background=1,
                                        inter_chromosomal=inter_chromosomal,
                                        scale=scale)
            variance = pca_info.explained_variance_ratio_

            with open(output_file, 'w') as o:
                o.write("sample\t{}\n".format("\t".join(["ev_{}".format(i+1) for i in range(len(sample_names))])))
                for i, v in enumerate(pca_res):
                    o.write("{}\t".format(sample_names[i]) + "\t".join([str(x) for x in v]) + "\n")
                o.write("variance\t" + "\t".join([str(v) for v in variance]) + "\n")

            pca_res_file = output_file
        else:
            pca_res_file = input_files[0]

        if plot_file is not None:
            import pandas
            import numpy as np
            pca_info = pandas.read_csv(pca_res_file, sep="\t")
            pca_res = pca_info.iloc[:-1, 1:]
            variance = pca_info.iloc[-1, 1:]

            if argument_sample_names is None:
                sample_names = pca_info.iloc[:-1, 0]

            import matplotlib
            matplotlib.use("agg")
            from fanc.plotting.statistics import pca_plot
            fig, ax = pca_plot(np.array(pca_res), variance=variance, names=sample_names,
                               markers=args.markers, colors=args.colors,
                               eigenvectors=eigenvectors)
            ax.set_title("PCA sample size {}".format(sample_size))
            fig.savefig(plot_file)

    finally:
        for matrix in matrices:
            matrix.close()


def loops_parser():
    parser = argparse.ArgumentParser(
        prog="fanc loops",
        description='Call loops in a Hic object using FAN-C '
                    'implementation of HICCUPS. See. Rao, Huntley et '
                    'al. (2014), Cell, for details.'
    )

    parser.add_argument(
        'input',
        help='Input Hic file'
    )

    parser.add_argument(
        'output',
        nargs='?',
        help='Output file. If input file is already '
             'a FAN-C compatible loops object, filtering '
             'can also be done in place.'
    )

    parser.add_argument(
        '-c', '--chromosomes', dest='chromosomes',
        nargs='+',
        help='Chromosomes to be investigated.'
    )

    parser.add_argument(
        '-p', '--peak-size', dest='peak_size',
        type=int,
        help='Size of the expected peak in pixels. '
             'If not set, will be estimated to '
             'correspond to ~ 25kb.'
    )

    parser.add_argument(
        '-w', '--neighborhood-width', dest='width',
        type=int,
        help='Width of the investigated area surrounding '
             'pixels. If not set, will be estimated at p+3'
    )

    parser.add_argument(
        '-t', '--threads', dest='threads',
        type=int,
        default=1,
        help='Number of threads for parallel processing. ' 
             'Default: 1 - it is advised to set this as '
             'high as possible, since loop calling is '
             'very computationally expensive!'
    )

    parser.add_argument(
        '-d', '--min-distance', dest='min_dist',
        type=int,
        help='Minimum distance in bins for two loci '
             'to be considered as loops. Default: peak size. '
             'Set this value higher to exclude loops '
             'close to the diagonal.'
    )

    parser.add_argument(
        '-m', '--mappability', dest='mappability_global_cutoff',
        type=float,
        help='Global mappability filter for all neighborhoods.'
             'Minimum mappable fraction of a pixel ' 
             'neighborhood to consider pixel as loop. '
             'Can be overridden by filters for local neighborhoods.'
    )

    parser.add_argument(
        '--mappability-donut', dest='mappability_donut_cutoff',
        type=float,
        help='Mappability filter for donut neighborhood. '
             'Value between 0 and 1.'
    )

    parser.add_argument(
        '--mappability-horizontal', dest='mappability_horizontal_cutoff',
        type=float,
        help='Mappability filter for horizontal neighborhood. '
             'Value between 0 and 1.'
    )

    parser.add_argument(
        '--mappability-vertical', dest='mappability_vertical_cutoff',
        type=float,
        help='Mappability filter for vertical neighborhood. '
             'Value between 0 and 1.'
    )

    parser.add_argument(
        '--mappability-lower-left', dest='mappability_lower_left_cutoff',
        type=float,
        help='Mappability filter for lower-left neighborhood. '
             'Value between 0 and 1.'
    )

    parser.add_argument(
        '-q', '--fdr', dest='fdr_global_cutoff',
        type=float,
        help='Global FDR filter all neighborhoods. '
             'Individual neighborhood filters can '
             'override this global setting. '
             'Value between 0 and 1.'
    )

    parser.add_argument(
        '--fdr-donut', dest='fdr_donut_cutoff',
        type=float,
        help='FDR filter for donut neighborhood. '
             'Value between 0 and 1.'
    )

    parser.add_argument(
        '--fdr-horizontal', dest='fdr_horizontal_cutoff',
        type=float,
        help='FDR filter for horizontal neighborhood. '
             'Value between 0 and 1.'
    )

    parser.add_argument(
        '--fdr-vertical', dest='fdr_vertical_cutoff',
        type=float,
        help='FDR filter for vertical neighborhood. '
             'Value between 0 and 1.'
    )

    parser.add_argument(
        '--fdr-lower-left', dest='fdr_lower_left_cutoff',
        type=float,
        help='FDR filter for lower-left neighborhood. '
             'Value between 0 and 1.'
    )

    parser.add_argument(
        '-e', '--enrichment', dest='enrichment_global_cutoff',
        type=float,
        help='Global observed/expected filter all neighborhoods. '
             'Individual neighborhood filters can '
             'override this global setting. '
    )

    parser.add_argument(
        '--enrichment-donut', dest='enrichment_donut_cutoff',
        type=float,
        help='Observed/expected enrichment filter for donut neighborhood. '
    )

    parser.add_argument(
        '--enrichment-horizontal', dest='enrichment_horizontal_cutoff',
        type=float,
        help='Observed/expected enrichment filter for horizontal neighborhood. '
    )

    parser.add_argument(
        '--enrichment-vertical', dest='enrichment_vertical_cutoff',
        type=float,
        help='Observed/expected enrichment filter for vertical neighborhood. '
    )

    parser.add_argument(
        '--enrichment-lower-left', dest='enrichment_lower_left_cutoff',
        type=float,
        help='Observed/expected enrichment filter for lower-left neighborhood. '
    )

    parser.add_argument(
        '-o', '--observed', dest='observed',
        type=int,
        help='Minimum observed value '
             '(integer, uncorrected). Default: 1'
    )

    parser.add_argument(
        '--rh-filter', dest='rh_filter',
        action='store_true',
        default=False,
        help='Filter peaks as in Rao, Huntley et al. (2014), Cell. '
             'It only retains peaks that are at least 2-fold enriched '
             'over either the donut or lower-left neighborhood, '
             'at least 1.5-fold enriched over the horizontal and '
             'vertical neighborhoods, at least 1.75-fold enriched '
             'over both the donut and lower-left neighborhood, and '
             'have an FDR <= 0.1 in every neighborhood '
    )

    parser.add_argument(
        '--sge', dest='sge',
        action='store_true',
        default=False,
        help='Run on SGE cluster. This option is highly '
             'recommended if you are running "fanc loops" on '
             'a Sun/Oracle Grid Engine Cluster. The "-t" option '
             'specifies the number of SGE slots if this flag '
             'is set. The main process will then submit jobs '
             'to the grid using "gridmap" and collect the '
             'results. (https://gridmap.readthedocs.io/) '
    )

    parser.add_argument(
        '--batch-size', dest='batch_size',
        type=int,
        default=200,
        help='Width of submatrix examined per '
             'process. Default: 200'
    )

    parser.add_argument(
        '-j', '--merge-pixels', dest='merge',
        action='store_true',
        default=False,
        help='Merge individual pixels into peaks after filtering.'
    )

    parser.add_argument(
        '--merge-distance', dest='merge_distance',
        type=int,
        default=20000,
        help='Maximum distance in base pairs at '
             'which to merge two pixels. '
             'Default 20000'
    )

    parser.add_argument(
        '-s', '--remove-singlets', dest='remove_singlets',
        action='store_true',
        default=False,
        help='Remove isolated pixels after merging step.'
    )

    parser.add_argument(
        '--fdr-sum', dest='fdr_sum',
        type=float,
        help='FDR sum filter for merged peaks. '
             'Merged peaks where the sum of donut FDR values of '
             'all constituent pixels is larger than the specified '
             'cutoff are filtered.'
    )

    parser.add_argument(
        '-b', '--bedpe', dest='bedpe',
        help='BEDPE output file. When set, merged loops '
             'will be written to this file after all '
             'filtering steps have completed.'
    )

    parser.add_argument(
        '-f', '--force-overwrite', dest='force_overwrite',
        action='store_true',
        default=False,
        help='''If the specified output file exists, it will be 
                        overwritten without warning.'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        default=False,
        help='Work in temporary directory'
    )

    return parser


def loops(argv, **kwargs):
    parser = loops_parser()

    args = parser.parse_args(argv[2:])

    import os

    input_file = os.path.expanduser(args.input)
    output_file = os.path.expanduser(args.output) if args.output is not None else None
    chromosomes = args.chromosomes
    peak_size = args.peak_size
    width = args.width
    threads = args.threads
    min_dist = args.min_dist

    mappability_cutoff_global = args.mappability_global_cutoff
    mappability_cutoff_donut = args.mappability_donut_cutoff
    mappability_cutoff_horizontal = args.mappability_horizontal_cutoff
    mappability_cutoff_vertical = args.mappability_vertical_cutoff
    mappability_cutoff_lower_left = args.mappability_lower_left_cutoff

    enrichment_cutoff_global = args.enrichment_global_cutoff
    enrichment_cutoff_donut = args.enrichment_donut_cutoff
    enrichment_cutoff_horizontal = args.enrichment_horizontal_cutoff
    enrichment_cutoff_vertical = args.enrichment_vertical_cutoff
    enrichment_cutoff_lower_left = args.enrichment_lower_left_cutoff

    fdr_cutoff_global = args.fdr_global_cutoff
    fdr_cutoff_donut = args.fdr_donut_cutoff
    fdr_cutoff_horizontal = args.fdr_horizontal_cutoff
    fdr_cutoff_vertical = args.fdr_vertical_cutoff
    fdr_cutoff_lower_left = args.fdr_lower_left_cutoff

    min_observed = args.observed

    rh_filter = args.rh_filter
    sge = args.sge
    batch_size = args.batch_size
    merge = args.merge
    merge_distance = args.merge_distance
    remove_singlets = args.remove_singlets
    fdr_sum = args.fdr_sum
    bedpe_file = os.path.expanduser(args.bedpe) if args.bedpe is not None else None

    force_overwrite = args.force_overwrite
    tmp = args.tmp

    if not force_overwrite and output_file is not None and os.path.exists(output_file):
        parser.error("Output file {} exists! Use -f to force "
                     "overwriting it!".format(output_file))

    import fanc
    import fanc.peaks
    import shutil
    from genomic_regions.files import create_temporary_output

    tmp_input_files = []

    original_output_file = output_file
    if tmp:
        tmp = False
        if output_file is not None:
            output_file = create_temporary_output(output_file)
        tmp = True

    matrix = None
    try:
        matrix = fanc.load(input_file, mode='r', tmpdir=tmp)
        is_rh_peaks = isinstance(matrix, fanc.peaks.RaoPeakInfo)
        is_merged_peaks = isinstance(matrix, fanc.peaks.PeakInfo)
        is_matrix = isinstance(matrix, fanc.matrix.RegionMatrixContainer)

        if not is_matrix:
            parser.error("Input type {} not supported".format(type(matrix)))

        if (merge and not is_merged_peaks) or (is_matrix and not is_rh_peaks and not is_merged_peaks):
            if output_file is None:
                parser.error("Must provide output file when calling or merging peaks!")

        # perform pixel loop probability estimate
        ran_peak_calling = False
        if is_matrix and not is_rh_peaks and not is_merged_peaks:
            pk = fanc.peaks.RaoPeakCaller(p=peak_size, w_init=width, min_locus_dist=peak_size,
                                          n_processes=threads, slice_size=batch_size,
                                          cluster=sge, min_mappable_fraction=0.0)
            if chromosomes is not None:
                chromosome_pairs = [(chromosome, chromosome) for chromosome in chromosomes]
            else:
                chromosome_pairs = None

            if not merge:
                o = output_file
            else:
                o = create_temporary_output('test.peaks')
                tmp_input_files.append(o)

            peaks = pk.call_peaks(matrix, chromosome_pairs=chromosome_pairs, file_name=o)
            matrix.close()
            matrix = peaks
            is_rh_peaks = True
            ran_peak_calling = True

        # filter pixels based on loop probability
        if is_rh_peaks and not is_merged_peaks:
            if not ran_peak_calling and not merge and output_file is not None:
                matrix.close()
                import shutil
                shutil.copy(matrix.file.filename, output_file)
                matrix = fanc.load(output_file, mode='a')

            filters = []
            if (fdr_cutoff_global is not None or fdr_cutoff_donut is not None or
                    fdr_cutoff_horizontal is not None or fdr_cutoff_lower_left is not None or
                    fdr_cutoff_vertical is not None):
                fdr_mask = matrix.add_mask_description('fdr', 'FDR cutoff filter')
                fdr_filter = fanc.peaks.FdrPeakFilter(fdr_cutoff=fdr_cutoff_global,
                                                      fdr_ll_cutoff=fdr_cutoff_lower_left,
                                                      fdr_d_cutoff=fdr_cutoff_donut,
                                                      fdr_h_cutoff=fdr_cutoff_horizontal,
                                                      fdr_v_cutoff=fdr_cutoff_vertical,
                                                      mask=fdr_mask)
                filters.append(fdr_filter)
            
            if (mappability_cutoff_global is not None or
                    mappability_cutoff_donut is not None or
                    mappability_cutoff_horizontal is not None or
                    mappability_cutoff_lower_left is not None or
                    mappability_cutoff_vertical is not None):
                mappability_mask = matrix.add_mask_description('mappability', 'Mappability filter')
                mappability_filter = fanc.peaks.MappabilityPeakFilter(
                    mappability_cutoff=mappability_cutoff_global,
                    mappability_ll_cutoff=mappability_cutoff_lower_left,
                    mappability_d_cutoff=mappability_cutoff_donut,
                    mappability_h_cutoff=mappability_cutoff_horizontal,
                    mappability_v_cutoff=mappability_cutoff_vertical,
                    mask=mappability_mask)
                filters.append(mappability_filter)
                
            if (enrichment_cutoff_global is not None or
                    enrichment_cutoff_donut is not None or
                    enrichment_cutoff_horizontal is not None or
                    enrichment_cutoff_lower_left is not None or
                    enrichment_cutoff_vertical is not None):
                enrichment_mask = matrix.add_mask_description('enrichment', 'O/E filter')
                enrichment_filter = fanc.peaks.EnrichmentPeakFilter(
                    enrichment_cutoff=enrichment_cutoff_global,
                    enrichment_ll_cutoff=enrichment_cutoff_lower_left,
                    enrichment_d_cutoff=enrichment_cutoff_donut,
                    enrichment_h_cutoff=enrichment_cutoff_horizontal,
                    enrichment_v_cutoff=enrichment_cutoff_vertical,
                    mask=enrichment_mask)
                filters.append(enrichment_filter)

            if min_dist is not None:
                if peak_size is None or min_dist > peak_size:
                    distance_mask = matrix.add_mask_description('distance', 'Min distance filter')
                    distance_filter = fanc.peaks.DistancePeakFilter(min_dist,
                                                                    mask=distance_mask)
                    filters.append(distance_filter)

            if min_observed is not None:
                observed_mask = matrix.add_mask_description('observed', 'Min observed filter')
                observed_filter = fanc.peaks.ObservedPeakFilter(min_observed,
                                                                mask=observed_mask)
                filters.append(observed_filter)

            if rh_filter:
                rh_mask = matrix.add_mask_description('rao_huntley', "Rao Huntley filter")
                rh = fanc.peaks.RaoPeakFilter(mask=rh_mask)
                filters.append(rh)

            if len(filters) > 0:
                for f in filters:
                    matrix.filter(f, queue=True)
                matrix.run_queued_filters()

        if merge and not is_merged_peaks:
            merged_peaks = matrix.merged_peaks(output_file, euclidian_distance=merge_distance)
            matrix.close()
            matrix = merged_peaks
            is_merged_peaks = True

        if is_merged_peaks:
            filters = []
            if remove_singlets:
                mask = matrix.add_mask_description('rao', 'Mask singlet peaks'
                                                          'with a q-value sum > .02')
                rao_filter = fanc.peaks.RaoMergedPeakFilter(mask=mask)
                filters.append(rao_filter)
            if fdr_sum is not None:
                mask = matrix.add_mask_description('fdr_sum', 'Mask merged peaks '
                                                              'with a q-value sum > cutoff')
                fdr_sum_filter = fanc.peaks.FdrSumFilter(cutoff=fdr_sum, mask=mask)
                filters.append(fdr_sum_filter)

            if len(filters) > 0:
                for f in filters:
                    matrix.filter(f, queue=True)
                matrix.run_queued_filters()

        if is_merged_peaks and bedpe_file is not None:
            logger.info("Exporting to BEDPE")
            matrix.to_bedpe(bedpe_file)

    finally:
        if matrix is not None:
            matrix.close()

        for tmp_input_file in tmp_input_files:
            os.remove(tmp_input_file)

        if tmp and output_file is not None and original_output_file is not None:
            shutil.copy(output_file, original_output_file)
            os.remove(output_file)


def overlap_peaks_parser():
    parser = argparse.ArgumentParser(
        prog="fanc overlap-peaks",
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


def overlap_peaks(argv, **kwargs):
    parser = overlap_peaks_parser()

    args = parser.parse_args(argv[2:])

    import fanc.peaks as kn
    from genomic_regions.files import create_temporary_copy

    original_input_paths = [os.path.expanduser(i) for i in args.input]
    if not args.names:
        names = [os.path.splitext(os.path.basename(i))[0] for i in original_input_paths]
    else:
        names = args.names

    if len(original_input_paths) < 2:
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
    import fanc

    peaks = [fanc.load(i, mode="r") for i in input_paths]

    if not args.distance:
        distance = peaks[0].bin_size*3
    else:
        distance = args.distance

    stats, merged = kn.overlap_peaks({n: p for n, p in zip(names, peaks)}, max_distance=distance)

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


def boundaries_parser():
    parser = argparse.ArgumentParser(
        prog="fanc boundaries",
        description='Determine domain boundaries'
    )

    parser.add_argument(
        'input',
        help='Input InsulationScores or region-based file with scores ' 
             '(BED, BigWig, ...)'
    )

    parser.add_argument(
        'output',
        help="Path for boundary BED file. When specifying "
             "multiple window sizes or if input file "
             "has multiple scores in it, this forms "
             "the output file prefix and will be appended by "
             "'<window size>.bed'"
    )

    parser.add_argument(
        '-w', '--window-sizes', dest='window',
        help='Insulation index window size to calculate boundaries on. '
             'Separate multiple window sizes with comma, e.g. '
             '1mb,500kb,100kb'
    )

    parser.add_argument(
        '-d', '--delta', dest='delta',
        type=int, default=3,
        help='Window size for calculating the delta '
             'vector (in bins). Calculation takes into '
             'account d bins upstream and d bins '
             'downstream for a total window size of '
             '2*d + 1 bins. Default 3.'
    )

    parser.add_argument(
        '-s', '--min-score', dest='min_score',
        type=float,
        help='Report only peaks where the two '
             'surrounding extrema of the delta '
             'vector have at least this difference '
             'in height. Default: no threshold.'
    )

    parser.add_argument(
        '-x', '--sub-bin-precision', dest='sub_bin_precision',
        action='store_true',
        help='Report boundary positions with sub-bin precision. '
             'This works because the minimum or the '
             'the insulation score track can be determined '
             'with sub-bin precision. Default: False'
    )

    parser.add_argument(
        '-l', '--log', dest='log',
        action='store_true',
        default=False,
        help='log-transform index values before '
             'boundary calling.'
    )

    parser.add_argument(
        '-m', '--maxima', dest='maxima',
        action='store_true',
        default=False,
        help='Call maxima of the insulation score '
             'instead of minima.'
    )

    return parser


def boundaries(argv, **kwargs):
    parser = boundaries_parser()

    args = parser.parse_args(argv[2:])

    import os

    input_file = os.path.expanduser(args.input)
    output_path = os.path.expanduser(args.output)
    windows = args.window
    delta = args.delta
    min_score = args.min_score
    sub_bin_precision = args.sub_bin_precision
    maxima = args.maxima
    log = args.log

    import fanc
    from fanc.tools.general import human_format, str_to_int
    import genomic_regions as gr
    from fanc.architecture.domains import InsulationScores, InsulationScore, Boundaries
    from fanc.architecture.comparisons import ComparisonScores

    insulation = None
    try:
        insulation = fanc.load(input_file, mode='r')
        if isinstance(insulation, InsulationScores) or isinstance(insulation, ComparisonScores):
            if windows is None:
                windows = insulation.window_sizes
            else:
                windows = [str_to_int(window) for window in windows.split(",")]

            if len(windows) > 1:
                for window in windows:
                    regions = insulation.score_regions(window)
                    try:
                        b = Boundaries.from_insulation_score(regions, min_score=min_score,
                                                            delta_window=delta,
                                                            sub_bin_precision=sub_bin_precision,
                                                            log=log, call_maxima=maxima)
                        regions.close()
                        b.to_bed('{}_{}b.bed'.format(output_path, human_format(window, lowercase=True)))
                    finally:
                        b.close()
            else:
                regions = insulation.score_regions(windows[0])
                try:
                    b = Boundaries.from_insulation_score(regions, min_score=min_score,
                                                        delta_window=delta,
                                                        sub_bin_precision=sub_bin_precision,
                                                        log=log, call_maxima=maxima)
                    b.to_bed(output_path)
                finally:
                    b.close()

        elif isinstance(insulation, InsulationScore) or isinstance(insulation, gr.RegionBased):
            regions = insulation.regions
            try:
                b = Boundaries.from_insulation_score(regions, min_score=min_score,
                                                    delta_window=delta,
                                                    sub_bin_precision=sub_bin_precision,
                                                    log=log, call_maxima=maxima)
                b.to_bed(output_path)
            finally:
                b.close()
        else:
            parser.error("Cannot recognise input file format!")
    finally:
        if insulation is not None:
            insulation.close()

    logger.info("All done.")


def compare_parser():
    parser = argparse.ArgumentParser(
        prog="fanc compare",
        description='Create pairwise comparisons of Hi-C comparison maps'
    )
    parser.add_argument(
        'input',
        nargs=2,
        help='Input matrix (e.g. Hic) files.'
    )
    parser.add_argument(
        'output',
        help='Output ComparisonMatrix file.'
    )

    parser.add_argument(
        '-c', '--comparison', dest='comparison',
        default='fold-change',
        help='Type of comparison. Default: fold-change, '
             'other options are: difference'
    )

    parser.add_argument(
        '-o', '--output-format', dest='output_format',
        help='Output format for region-based comparisons. '
             'Only relevant when using BED, GFF, or another '
             'region-based format as input.'
    )

    parser.add_argument(
        '-S', '--no-scale', dest='scale',
        action='store_false',
        default=True,
        help='Do not scale input matrices to the ' 
             'same number of valid pairs'
    )

    parser.add_argument(
        '-l', '--log', dest='log',
        action='store_true',
        default=False,
        help='''Log2-convert comparison values (AFTER the comparison)'''
    )

    parser.add_argument(
        '--log-matrix', dest='log_matrix',
        action='store_true',
        default=False,
        help='''Log2-convert matrices (BEFORE the comparison)'''
    )

    parser.add_argument(
        '-Z', '--ignore-zero', dest='ignore_zero',
        action='store_true',
        default=False,
        help='''Do not consider pixels where one matrix entry is zero'''
    )

    parser.add_argument(
        '-I', '--ignore-infinite', dest='ignore_infinite',
        action='store_true',
        default=False,
        help='Do not consider pixels where the comparison yields "inf"'
    )

    parser.add_argument(
        '-e', '--observed-expected', dest='oe',
        action='store_true',
        default=False,
        help='O/E transform matrix values before comparison. '
             'Only has an effect on matrix comparisons.'
    )

    parser.add_argument(
        '-u', '--uncorrected', dest='norm',
        action='store_false',
        default=True,
        help='Compare uncorrected matrices. Only has an '
             'effect on matrix comparisons.'
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        default=False,
        help='Work in temporary directory'
    )
    return parser


def compare(argv, **kwargs):
    parser = compare_parser()

    args = parser.parse_args(argv[2:])

    import os
    input_file1 = os.path.expanduser(args.input[0])
    input_file2 = os.path.expanduser(args.input[1])
    output_file = os.path.expanduser(args.output)
    scale = args.scale
    comparison = args.comparison
    filter_zero = args.ignore_zero
    filter_infinite = args.ignore_infinite
    output_format = args.output_format
    norm = args.norm
    oe = args.oe
    log = args.log
    log_matrix = args.log_matrix
    tmp = args.tmp

    import fanc
    from fanc.architecture.comparisons import FoldChangeMatrix, DifferenceMatrix, NonzeroFilter, \
        FoldChangeScores, DifferenceScores, FoldChangeRegions, DifferenceRegions
    from fanc.matrix import RegionMatrixContainer
    from fanc.architecture.domains import RegionScoreParameterTable
    from genomic_regions import RegionBased

    matrix1 = fanc.load(input_file1, mode='r')
    matrix2 = fanc.load(input_file2, mode='r')

    if isinstance(matrix1, RegionMatrixContainer) and isinstance(matrix2, RegionMatrixContainer):
        if filter_zero:
            logger.info("Enabling zero filter")

        ComparisonMatrix = None
        if comparison == 'fold-change' or comparison == 'fc':
            ComparisonMatrix = FoldChangeMatrix
        elif comparison == 'difference' or comparison == 'diff':
            ComparisonMatrix = DifferenceMatrix
        else:
            parser.error("Comparison type -c {} not recognised!".format(comparison))

        cmp = ComparisonMatrix.from_matrices(matrix1, matrix2, file_name=output_file,
                                             tmpdir=tmp, mode='w',
                                             scale=scale,
                                             log_cmp=log,
                                             ignore_infinite=filter_infinite,
                                             ignore_zeros=filter_zero,
                                             oe=oe, norm=norm, log=log_matrix)
        cmp.close()
    elif isinstance(matrix1, RegionScoreParameterTable) and isinstance(matrix2, RegionScoreParameterTable):
        ComparisonScores = None
        if comparison == 'fold-change' or comparison == 'fc':
            ComparisonScores = FoldChangeScores
        elif comparison == 'difference' or comparison == 'diff':
            ComparisonScores = DifferenceScores
        else:
            parser.error("Comparison type -c {} not recognised!".format(comparison))

        cmp = ComparisonScores.from_scores(matrix1, matrix2, file_name=output_file,
                                           tmpdir=tmp, mode='w', log=log)
        cmp.close()
    elif isinstance(matrix1, RegionBased) and isinstance(matrix2, RegionBased):
        ComparisonRegions = None
        if comparison == 'fold-change' or comparison == 'fc':
            ComparisonRegions = FoldChangeRegions
        elif comparison == 'difference' or comparison == 'diff':
            ComparisonRegions = DifferenceRegions
        else:
            parser.error("Comparison type -c {} not recognised!".format(comparison))

        if output_format is None:
            comparison_output = output_file
        else:
            output_format = output_format.lower()
            comparison_output = None
        cmp = ComparisonRegions.from_regions(matrix1, matrix2, file_name=comparison_output,
                                             tmpdir=tmp, mode='w', log=log)

        if output_format == 'bed':
            from fanc.tools.files import write_bed
            write_bed(output_file, cmp.regions)
        elif output_format == 'bigwig' or output_format == 'bw':
            from fanc.tools.files import write_bigwig
            write_bigwig(output_file, cmp.regions)
        elif output_format == 'gff':
            from fanc.tools.files import write_gff
            write_gff(output_file, cmp.regions)
        cmp.close()


def directionality_parser():
    parser = argparse.ArgumentParser(
        prog="fanc directionality",
        description='Calculate directionality index for Hic object'
    )
    parser.add_argument(
        'input',
        help='Input matrix (Hi-C, fold-change map, ...)'
    )
    parser.add_argument(
        'output',
        nargs='?',
        help='Output file. Format will be determined by '
             '"-o". By default, this is a FAN-C DirectionalityIndexes '
             'object, for maximum compatibility with other analyses. '
             'If you choose a text-based output format (BED, GFF, '
             'BigWig), this parameter will be the file prefix, and '
             'the window size will be appended.'
             'If not specified and output format is one of '
             'bed, gff, or bigwig, the input file name forms the '
             'output file prefix.'
    )

    parser.add_argument(
        '-o', '--output-format', dest='output_format',
        default='directionality_index',
        help='Format of the output file. '
             'By default, this is a FAN-C DirectionalityIndex '
             'object, for maximum compatibility with other '
             'analyses. Other options are "bed", "bigwig", '
             'and "gff"'
    )

    parser.add_argument(
        '-w', '--window-sizes', dest='window_sizes',
        nargs='+',
        help='Window sizes in base pairs. You can also use '
             'abbreviated number format (i.e. 1.5M, 250kb, etc). '
             'If not specified, will choose the window sizes based '
             'on the matrix resolution r when calculating scores. '
             'Specifically: r*3, r*5, r*7, r*10, and r*15'
    )

    parser.add_argument(
        '-r', '--region', dest='region',
        help='Region selector (<chr>:<start>-<end>) '
             'to only calculate directionality index '
             'for this region.'
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        default=False,
        help='Work in temporary directory'
    )

    return parser


def directionality(argv, **kwargs):
    parser = directionality_parser()

    args = parser.parse_args(argv[2:])

    from fanc.architecture.domains import DirectionalityIndexes
    import os

    input_file = os.path.expanduser(args.input)
    output_file = os.path.expanduser(args.output) if args.output is not None else None
    output_format = args.output_format
    sub_region = args.region
    window_sizes = args.window_sizes
    tmp = args.tmp

    if output_format.lower() not in ['directionality_index', 'bed', 'gff', 'bigwig', 'bw']:
        parser.error("Output format must be one of directionality_index, bed, gff, or bigwig!")

    _domain_scores(parser, input_file, output_file, output_format,
                   'directionality_index', DirectionalityIndexes,
                   sub_region=sub_region, tmp=tmp,
                   tmpdir=tmp, window_sizes=window_sizes)


def insulation_parser():
    parser = argparse.ArgumentParser(
        prog="fanc insulation",
        description='Calculate insulation scores for Hic object'
    )
    parser.add_argument(
        'input',
        help='Input matrix (Hi-C, fold-change map, ...)'
    )
    parser.add_argument(
        'output',
        nargs='?',
        help='Output file. Format will be determined by '
             '"-o". By default, this is a FAN-C InsulationScores '
             'object, for maximum compatibility with other analyses. '
             'If you choose a text-based output format (BED, GFF, '
             'BigWig), this parameter will be the file prefix, and '
             'the window size will be appended.'
             'If not specified and output format is one of '
             'bed, gff, or bigwig, the input file name forms the '
             'output file prefix.'
    )

    parser.add_argument(
        '-o', '--output-format', dest='output_format',
        default='insulation_score',
        help='Format of the output file. '
             'By default, this is a FAN-C InsulationScore '
             'object, for maximum compatibility with other '
             'analyses. Other options are "bed", "bigwig", '
             'and "gff"'
    )

    parser.add_argument(
        '-w', '--window-sizes', dest='window_sizes',
        nargs='+',
        help='Window sizes in base pairs. You can also use '
             'abbreviated number format (i.e. 1.5M, 250kb, etc). '
             'If not specified, will choose the window sizes based '
             'on the matrix resolution r when calculating scores. '
             'Specifically: r*3, r*5, r*7, r*10, and r*15'
    )

    parser.add_argument(
        '-r', '--region', dest='region',
        help='Region selector (<chr>:<start>-<end>) '
             'to only calculate II for this region.'
    )

    parser.add_argument(
        '-i', '--impute', dest='impute',
        action='store_true',
        default=False,
        help='Impute missing values in matrix. If set, '
             'missing matrix values (where an entire Hi-C bin '
             'has 0 contacts) will be replaced by the '
             'expected value at the given distance.'
    )

    parser.add_argument(
        '--offset', dest='offset',
        type=int,
        default=0,
        help='Window offset in base pairs from the diagonal.'
    )

    parser.add_argument(
        '-L', '--no-log', dest='log',
        action='store_false',
        default=True,
        help='Do not log2-transform insulation index after '
             'normalisation. Log-transformation roughly centers values '
             'around 0, but if you need this to be exactly '
             'centered, use the "--geom-mean" option.'
    )

    parser.add_argument(
        '-N', '--no-norm', dest='normalise',
        action='store_false',
        default=True,
        help='Do not normalise index to insulation average '
             'Default is whole chromosome normalisation - to normalise '
             'to smaller regions, use --normalisation-window.'
    )

    parser.add_argument(
        '--normalisation-window', dest='normalisation_window',
        type=int,
        help='Size of the normalisation window (moving average) '
             'in bins. Default: whole chromosome.'
    )

    parser.add_argument(
        '-s', '--subtract-mean', dest='subtract',
        action='store_true',
        default=False,
        help='Subtract mean instead of dividing by it when "-n" is enabled. '
             'You probably do not want this, unless you are working with '
             'log-transformed matrices (e.g. fold-change matrices)'
    )

    parser.add_argument(
        '-g', '--geom-mean', dest='geom_mean',
        action='store_true',
        default=False,
        help='Use geometric mean for normalisation (rather than arithmetic mean). '
             'Useful in conjunction with --log to center the distribution at 0. '
             'This is very important when comparing insulation scores, for example '
             'using the "fanc compare" command!'
    )

    parser.add_argument(
        '--trim-mean', dest='trim_mean',
        type=float,
        default=0.0,
        help='Use a trimmed mean for insulation index normalisation with '
             'this cutoff (fraction of scores)'
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        default=False,
        help='Work in temporary directory'
    )

    return parser


def insulation(argv, **kwargs):
    parser = insulation_parser()

    args = parser.parse_args(argv[2:])

    from fanc.architecture.domains import InsulationScores
    import os

    input_file = os.path.expanduser(args.input)
    output_file = os.path.expanduser(args.output) if args.output is not None else None
    output_format = args.output_format
    sub_region = args.region
    tmp = args.tmp
    window_sizes = args.window_sizes
    impute = args.impute
    normalise = args.normalise
    offset = args.offset
    subtract_mean = args.subtract
    trim_mean = args.trim_mean
    geom_mean = args.geom_mean
    log = args.log
    normalisation_window = args.normalisation_window

    if output_format.lower() not in ['insulation_score', 'bed', 'gff', 'bigwig', 'bw']:
        parser.error("Output format must be one of insulation_score, bed, gff, or bigwig!")

    _domain_scores(parser, input_file, output_file, output_format,
                   'insulation_score', InsulationScores,
                   sub_region=sub_region, tmp=tmp,
                   tmpdir=tmp, window_sizes=window_sizes,
                   impute_missing=impute, normalise=normalise, window_offset=offset,
                   subtract_mean=subtract_mean, log=log,
                   trim_mean_proportion=trim_mean, geometric_mean=geom_mean,
                   normalisation_window=normalisation_window)


def _domain_scores(parser, input_file, output_file, output_format, default_output_format,
                   score_class, sub_region=None, tmp=False,
                   **kwargs):
    import fanc
    from fanc.tools.general import str_to_int, human_format
    output_format = output_format.lower()
    window_sizes = kwargs.get('window_sizes', None)

    if not output_format == default_output_format:
        if output_file is None:
            output_prefix = input_file
        else:
            output_prefix = output_file
        output_file = None
        sub_output_file = None
    else:
        if sub_region is not None:
            sub_output_file = output_file
            output_file = None
            output_prefix = None
        else:
            output_prefix = None
            sub_output_file = None

    scores = None
    try:
        matrix = fanc.load(input_file, mode='r')

        if not isinstance(matrix, score_class):
            if output_file is None and output_format == default_output_format:
                parser.error("Output file cannot be empty when "
                             "choosing default output format!")

            if window_sizes is not None:
                window_sizes = [str_to_int(window_size) for window_size in window_sizes]
            else:
                bin_size = matrix.bin_size
                window_sizes = [
                    bin_size * 3, bin_size * 5, bin_size * 7, bin_size * 10, bin_size * 15
                ]
            kwargs['window_sizes'] = window_sizes
            logger.info("Chosen window sizes: {}".format(" ".join([str(w) for w in window_sizes])))
            scores = score_class.from_hic(matrix, file_name=output_file,
                                          **kwargs)
        else:
            scores = matrix

            if window_sizes is None:
                window_sizes = scores.window_sizes
            else:
                window_sizes = [str_to_int(window_size) for window_size in window_sizes]

        if output_format in ['bed', 'gff', 'bigwig', 'bw']:
            for window_size in window_sizes:
                window_size_human = human_format(window_size).lower() + 'b'

                output_file = output_prefix + '_{}.{}'.format(window_size_human, output_format.lower())

                if output_format in ['bw', 'bigwig']:
                    scores.to_bigwig(output_file, window_size, subset=sub_region)

                elif output_format in ['bed']:
                    scores.to_bed(output_file, window_size, subset=sub_region)

                elif output_format in ['gff']:
                    scores.to_gff(output_file, window_size, subset=sub_region)
        elif sub_region is not None:
            if output_file is None:
                parser.error("Output file cannot be empty when "
                             "choosing default output format!")
            sub_scores = score_class(sub_output_file, mode='w', tmpdir=tmp)
            sub_scores.add_regions(scores.regions(sub_region), preserve_attributes=True)
            sub_scores.close()
        elif output_file is None:
            logger.debug("No action specified, printing insulation score info")
            print("Window sizes available in object:")
            print("{}".format(" ".join([human_format(window_size, lowercase=True) + 'b'
                                        for window_size in scores.window_sizes])))

    finally:
        if scores is not None:
            scores.close()


def compartments_parser():
    parser = argparse.ArgumentParser(
        prog="fanc compartments",
        description='Calculate AB compartment matrix'
    )
    parser.add_argument(
        'matrix',
        help='Input matrix (Hi-C, fold-change map, ...) or '
             'existing AB compartment matrix.'
    )
    parser.add_argument(
        'ab_compartments',
        nargs='?',
        help='AB compartment matrix file.'
    )

    parser.add_argument(
        '-d', '--domains', dest='domains',
        help='Write AB domains to this file. AB domains '
             'are output in BED format, and include the '
             'domains type (A/B) in the name field, and the '
             'eigenvector values (averaged across all bins '
             'in the domain) in the score field'
    )

    parser.add_argument(
        '-v', '--eigenvector', dest='eigenvector',
        help='Write eigenvector values to this file.'
             'Output format is BED, containing of each '
             'matrix bin. The score field contains '
             'the eigenvector value of the bin.'
    )

    parser.add_argument(
        '-e', '--enrichment-profile', dest='enrichment_file',
        help='Plot AB enrichment profile to this file.'
    )

    parser.add_argument(
        '-m', '--enrichment-matrix', dest='matrix_file',
        help='Path to save enrichment profile matrix '
             '(numpy txt format)'
    )

    parser.add_argument(
        '-g', '--genome', dest='genome',
        help='Genome file. Used to "orient" the '
             'eigenvector values (change sign) using '
             'the average GC content of domains. '
             'Possible input files are '
             'FASTA, folder with FASTA, '
             'comma-separated list of FASTA) '
             'used to change sign of eigenvector '
             'based on GC content.'
    )

    parser.add_argument(
        '-w', '--whole-genome', dest='whole_genome',
        action='store_true',
        default=False,
        help='Calculate AB compartments on the whole genome '
             'matrix, instead of individual chromosomes. '
             'Only enable if you are sure your use-case requires it. '
             'This is likely to introduce artefacts when working '
             'with matrices that have been normalised per-chromosome.'
    )

    parser.add_argument(
        '-r', '--region', dest='region',
        help='Only outputs domains / eigenvector values in '
             'this region. Only works with the -d and -e '
             'arguments. Compartmentalisation is always '
             'calculated on the whole genome.'
    )

    parser.add_argument(
        '-i', '--eigenvector-index', dest='eigenvector_index',
        type=int,
        default=1,
        help='Eigenvector index. By default, the first '
             'eigenvector is output for eigenvector and domain '
             'analysis. Sometimes, it is useful to choose a '
             'higher eigenvector. E.g. for second eigenvector,'
             'specify "-i 2".'
    )

    parser.add_argument(
        '-p', '--enrichment-percentiles', dest='percentiles',
        type=float,
        nargs='+',
        default=(20.0, 40.0, 60.0, 80.0, 100.0),
        help='Percentiles to use for calculation of the ' 
             'enrichment profile. By default uses '
             '20, 40, 60, 80, 100. The 0 percentile is '
             'included by default.'
    )

    parser.add_argument(
        '-c', '--enrichment-colormap', dest='colormap',
        default='RdBu_r',
        help='Matplotlib colormap to use for plotting ' 
             'the enrichment profile.'
    )

    parser.add_argument(
        '-s', '--enrichment-symmetric-at', dest='symmetric_at',
        type=int,
        help='Make enrichment profile plot symmetric around this value '
             '(e.g. use 0 to ensure that 0 is in the center of the plot).'
    )

    parser.add_argument(
        '--enrichment-min', dest='vmin',
        type=float,
        default=-1,
        help='Minimum saturation value in enrichment profile. '
             'Default -1'
    )

    parser.add_argument(
        '--enrichment-max', dest='vmax',
        type=float,
        default=1,
        help='Maximum saturation value in enrichment profile. '
             'Default: 1'
    )

    parser.add_argument(
        '-G', '--only-gc', dest='only_gc',
        action='store_true',
        help='Only use GC content for enrichment profile calculation, '
             'not the correlation matrix eigenvector.'
    )
    parser.set_defaults(only_gc=False)

    parser.add_argument(
        '-x', '--enrichment-exclude', dest='exclude',
        nargs='+',
        help='Chromosome names to exclude from '
             'AB compartment and '
             'enrichment profile calculation'
    )

    parser.add_argument(
        '--compartment-strength', dest='compartment_strength_file',
        help='File for saving the compartment strength as defined by '
             'Flyamer, Gassler, and Imakaev et. al (2017) '
             '(https://pubmed.ncbi.nlm.nih.gov/28355183).'
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        default=False,
        help='Work in temporary directory'
    )

    parser.add_argument(
        '-f', '--force', dest='force',
        action='store_true',
        default=False,
        help='Force overwriting of output files.'
    )

    parser.add_argument(
        '--recalculate', dest='recalculate',
        action='store_true',
        default=False,
        help='Force recalculation of eigenvector even if a '
             'vector with the same parameters has previously '
             'been calculated.'
    )

    return parser


def compartments(argv, **kwargs):
    parser = compartments_parser()

    args = parser.parse_args(argv[2:])

    import os

    input_file = os.path.expanduser(args.matrix)
    output_file = os.path.expanduser(args.ab_compartments) if args.ab_compartments is not None else None
    domains_file = os.path.expanduser(args.domains) if args.domains is not None else None
    genome_file = os.path.expanduser(args.genome) if args.genome is not None else None
    eigenvector_file = os.path.expanduser(args.eigenvector) if args.eigenvector is not None else None
    compartment_strength_file = os.path.expanduser(args.compartment_strength_file) if args.compartment_strength_file is not None else None
    ev_index = args.eigenvector_index - 1
    region_subset = args.region
    whole_genome = args.whole_genome

    matrix_file = os.path.expanduser(args.matrix_file) if args.matrix_file is not None else None
    enrichment_file = os.path.expanduser(args.enrichment_file) if args.enrichment_file is not None else None

    percentiles = args.percentiles
    symmetric_at = args.symmetric_at
    cmap = args.colormap
    vmin = args.vmin
    vmax = args.vmax
    only_gc = args.only_gc
    exclude = args.exclude if args.exclude is not None else []
    recalculate = args.recalculate

    force = args.force
    tmp = args.tmp

    if not force and domains_file is not None and os.path.exists(domains_file):
        parser.error("Output file {} already exists. Use -f to override!".format(domains_file))
    if not force and matrix_file is not None and os.path.exists(matrix_file):
        parser.error("Output file {} already exists. Use -f to override!".format(matrix_file))
    if not force and enrichment_file is not None and os.path.exists(enrichment_file):
        parser.error("Output file {} already exists. Use -f to override!".format(enrichment_file))
    if eigenvector_file is not None and not force and os.path.exists(eigenvector_file):
            if matrix_file is not None or enrichment_file is not None:
                logger.warning("Found existing eigenvector at {}. "
                               "Will use the values in the file for "
                               "enrichment profile calculations. "
                               "This may not be what you want!".format(eigenvector_file))
            else:
                parser.error("Output file {} already exists. Use -f to override!".format(eigenvector_file))

    if matrix_file is not None or enrichment_file is not None:
        if output_file is None:
            parser.error("Need Input matrix for AB compartment profile!")

    import fanc
    import warnings
    from fanc.tools.files import write_bed
    from fanc.matrix import RegionMatrixContainer
    from fanc.architecture.compartments import ABCompartmentMatrix

    matrix = None
    ab_matrix = None
    try:
        matrix = fanc.load(input_file, tmpdir=tmp)
        # input is Hic or other matrix
        if not isinstance(matrix, ABCompartmentMatrix):
            if not isinstance(matrix, RegionMatrixContainer):
                parser.error("Input must be FAN-C matrix (e.g. Hic)")

            if output_file is None:
                parser.error("Must provide an output file if calculating or loading "
                             "AB compartment matrix. You can only omit this if your "
                             "first argument is already an AB compartment matrix!")

            if os.path.exists(output_file) and not recalculate:
                ab_matrix = fanc.load(output_file)
                if not isinstance(ab_matrix, ABCompartmentMatrix):
                    parser.error("Found existing file {}, but it is not an AB compartment matrix."
                                 "Use -f to overwrite it.")
                logger.warning("Found existing AB compartment matrix at {}. Will not recalculate - "
                               "use --recalculate to overwrite the existing file!".format(output_file))
            else:
                logger.info("Computing AB compartment matrix")
                ab_matrix = ABCompartmentMatrix.from_hic(matrix,
                                                         file_name=output_file, tmpdir=tmp,
                                                         per_chromosome=not whole_genome,
                                                         exclude_chromosomes=exclude)
        else:
            matrix.close()
            matrix = None
            ab_matrix = fanc.load(input_file, mode='a', tmpdir=tmp)

        # calculate eigenvector
        ev = None
        if eigenvector_file is not None:
            if os.path.exists(eigenvector_file) and not force:
                ev_regions = fanc.load(eigenvector_file)
                ev = [r.score for r in ev_regions.regions]
            else:
                ev = ab_matrix.eigenvector(sub_region=region_subset, genome=genome_file,
                                           eigenvector=ev_index, force=recalculate,
                                           exclude_chromosomes=exclude)
                regions = []
                for i, region in enumerate(ab_matrix.regions(region_subset)):
                    r = region.copy()
                    r.score = ev[i]
                    r.name = 'A' if ev[i] >= 0 else 'B'
                    regions.append(r)
                write_bed(eigenvector_file, regions)

        # Calculate domains
        if domains_file is not None:
            domains = ab_matrix.domains(sub_region=region_subset, genome=genome_file,
                                        eigenvector=ev_index, force=recalculate)
            write_bed(domains_file, domains)

        if compartment_strength_file is not None:
            import numpy as np
            m, cutoffs = ab_matrix.enrichment_profile(matrix,
                                                      percentiles=[0, 20, 40, 60, 80, 100],
                                                      per_chromosome=not whole_genome,
                                                      only_gc=only_gc,
                                                      symmetric_at=symmetric_at,
                                                      exclude_chromosomes=exclude,
                                                      eigenvector=ev,
                                                      genome=genome_file)
            aa = 2**m[0, 0]
            bb = 2**m[4, 4]
            ab = 2**m[0, 4]
            s = np.log((aa * bb) / ab ** 2)
            with open(compartment_strength_file, 'w') as o:
                o.write("AB-strength\t{}\n".format(s))
                o.write("AA\t{}\n".format(aa))
                o.write("BB\t{}\n".format(bb))
                o.write("AB\t{}\n".format(ab))

        # Calculate enrichment profile
        if matrix_file is not None or enrichment_file is not None:
            m, cutoffs = ab_matrix.enrichment_profile(matrix,
                                                      percentiles=percentiles,
                                                      per_chromosome=not whole_genome,
                                                      only_gc=only_gc,
                                                      symmetric_at=symmetric_at,
                                                      exclude_chromosomes=exclude,
                                                      eigenvector=ev,
                                                      genome=genome_file)

            if matrix_file is not None:
                import numpy as np
                np.savetxt(matrix_file, m)

            if enrichment_file is not None:
                import matplotlib
                matplotlib.use('agg')
                import matplotlib.pyplot as plt
                import matplotlib.gridspec as grd
                from fanc.plotting.statistics import saddle_plot

                fig = plt.figure(figsize=(5, 5), dpi=300)
                fig, axes = saddle_plot(m, cutoffs, colormap=cmap, vmin=vmin, vmax=vmax, only_gc=only_gc,
                                        fig=fig)
                fig.savefig(enrichment_file)
                plt.close(fig)
    finally:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if matrix is not None:
                matrix.close()
            if ab_matrix is not None:
                ab_matrix.close()


def expected_parser():
    parser = argparse.ArgumentParser(
        prog="fanc expected",
        description='Calculate Hi-C expected values '
                    '(distance decay)'
    )
    parser.add_argument(
        'input',
        nargs='+',
        help='Input matrix (Hi-C, fold-change map, ...)'
    )
    parser.add_argument(
        'output',
        help='Output expected contacts (tsv).'
    )

    parser.add_argument(
        '-p', '--plot', dest='plot_file',
        help='Output file for distance decay plot (pdf).'
    )

    parser.add_argument(
        '-l', '--labels', dest='labels',
        nargs='+',
        help='Labels for input objects.'
    )

    parser.add_argument(
        '-c', '--chromosome', dest='chromosome',
        help='Specific chromosome to calculate '
             'expected values for.'
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        default=False,
        help='Work in temporary directory'
    )

    parser.add_argument(
        '--recalculate', dest='recalculate',
        action='store_true',
        default=False,
        help='Recalculate expected values regardless of whether '
             'they are already stored in the matrix object.'
    )

    parser.add_argument(
        '-N', '--no-norm', dest='norm',
        action='store_false',
        default=True,
        help='Calculate expected values on unnormalised data.'
    )

    return parser


def expected(argv, **kwargs):
    parser = expected_parser()

    args = parser.parse_args(argv[2:])

    import fanc
    import os.path

    input_files = [os.path.expanduser(file_name) for file_name in args.input]
    output_file = os.path.expanduser(args.output)
    chromosome = args.chromosome
    plot_file = os.path.expanduser(args.plot_file) if args.plot_file is not None else None
    labels = args.labels
    recalculate = args.recalculate
    norm = args.norm
    tmp = args.tmp

    if labels is None:
        labels = ['Matrix_{}'.format(i) for i in range(len(input_files))]

    if len(labels) != len(input_files):
        parser.error("Must provide the same number of labels "
                     "({}) as input files ({})!".format(len(labels), len(input_files)))

    expected_values = dict()
    distances = dict()
    equal_distances = True
    for label, input_file in zip(labels, input_files):
        with fanc.load(input_file, mode='a', tmpdir=tmp) as matrix:
            intra_expected, intra_expected_chromosome, inter_expected = matrix.expected_values(
                force=recalculate, norm=norm)

            logger.info("Inter-chromosomal expected value: {}".format(inter_expected))

            if chromosome is not None:
                expected = intra_expected_chromosome[chromosome]
            else:
                expected = intra_expected

            bin_size = matrix.bin_size
            distance = [i * bin_size for i in range(len(expected))]

            expected_values[label] = expected
            distances[label] = distance

            equal_d = True
            for i in range(len(distance)):
                try:
                    if distance[i] != distances[labels[0]][i]:
                        equal_d = False
                        break
                except KeyError:
                    equal_d = False
                    break

            if not equal_d:
                equal_distances = False
                logger.warning("Distances of matrices are not equal. "
                               "({} vs {}). This may be due to "
                               "different binning of matrices.".format(labels[0], label))

    with open(output_file, 'w') as o:
        if equal_distances:
            header = "distance"
            for label in labels:
                header += "\t{}".format(label)
            o.write(header + "\n")

            for i in range(len(distances[labels[0]])):
                line = "{}".format(distances[labels[0]][i])
                for label in labels:
                    line += "\t{}".format(expected_values[label][i])
                o.write(line + "\n")
        else:
            logger.warning("Cannot write output file due to unequal distances / bin sizes!")

    if plot_file is not None:
        import matplotlib
        matplotlib.use("agg")
        import matplotlib.pyplot as plt
        from fanc.plotting.statistics import distance_decay_plot

        hics = [fanc.load(file_name) for file_name in input_files]

        fig, ax = plt.subplots(figsize=(4, 3), dpi=300)
        distance_decay_plot(*hics, chromosome=chromosome, ax=ax, labels=labels)
        fig.savefig(plot_file)
        plt.close(fig)

        for hic in hics:
            hic.close()


def subset_parser():
    parser = argparse.ArgumentParser(
        prog="fanc subset",
        description='Create a new Hic object by subsetting.'
    )
    parser.add_argument(
        'input',
        help='Input Hic file.'
    )
    parser.add_argument(
        'output',
        help='Output Hic file.'
    )

    parser.add_argument(
        'regions',
        nargs='+',
        help='List of regions that will be used in the output '
             'Hic object. All contacts between these regions '
             'will be in the output object. For example, '
             '"chr1 chr3" will result in a Hic object with '
             'all regions in chromosomes 1 and 3, plus all '
             'contacts within chromosome 1, all contacts within '
             'chromosome 3, and all contacts between chromosome '
             '1 and 3. "chr1" will only contain regions and contacts'
             'within chromosome 1.'
    )
    return parser


def subset(argv, **kwargs):
    parser = subset_parser()
    args = parser.parse_args(argv[2:])

    import os.path
    import fanc
    from fanc.compatibility.juicer import JuicerHic
    from fanc.compatibility.cooler import CoolerHic

    input_file = os.path.expanduser(args.input)
    output_file = os.path.expanduser(args.output)
    regions = args.regions

    old_hic = fanc.load(input_file, mode='r')
    if isinstance(old_hic, JuicerHic) or isinstance(old_hic, CoolerHic):
        parser.error("Cannot use subset on Juicer or Cooler files. Use FAN-C format if you want to subset.")
    new_hic = old_hic.subset(*regions, file_name=output_file)
    new_hic.close()


def aggregate_parser():
    parser = argparse.ArgumentParser(
        prog="fanc aggregate",
        description='Make aggregate plots with FAN-C'
    )

    parser.add_argument(
        'input',
        help='FAN-C matrix file (e.g. Hic)'
    )

    parser.add_argument(
        'regions',
        nargs='?',
        help='File with regions '
             '(BED, GFF, Tabix, ...) or '
             'region pairs (BEDPE)'
    )

    parser.add_argument(
        'output',
        nargs='?',
        help='Output AggregateMatrix file for further processing.'
             'See -p and -m option for aggregate plot and matrix, '
             'respectively.'
    )

    parser.add_argument(
        '-m', '--save-matrix', dest='matrix_file',
        help='Path to save aggregate matrix (numpy txt format)'
    )

    parser.add_argument(
        '-p', '--save-plot', dest='plot_file',
        help='Path to save aggregate plot (PDF)'
    )

    parser.add_argument(
        '--tads', dest='tads_preset',
        default=False,
        action='store_true',
        help='Use presets for aggregate TADs: '
             '--relative 1.0 --expected '
             '--log --vmin -1 --vmax 1'
    )

    parser.add_argument(
        '--tads-imakaev', dest='tads_imakaev_preset',
        default=False,
        action='store_true',
        help='Use presets for aggregate TADs: '
             '--relative 1.0 --expected'
             '--rescale'
    )

    parser.add_argument(
        '--loops', dest='loops_preset',
        default=False,
        action='store_true',
        help='Use presets for aggregate loops: '
             '--pixels 16 -l '
    )

    parser.add_argument(
        '--loop-strength', dest='loop_strength_file',
        help='Calculate loop strengths and save to file. Only works when providing BEDPE file, '
             'and Hi-C matrix.'
    )

    parser.add_argument(
        '--tad-strength', dest='tad_strength_file',
        help='Calculate tad strengths and save to file. Only works with --tads preset'
    )

    parser.add_argument(
        '-w', '--window', dest='window',
        help='Width of the region window used for aggregation. '
             'If set, will only use the center position from '
             'the input regions and extract a submatrix of width -w '
             'around this region.'
    )

    parser.add_argument(
        '--pixels', dest='pixels',
        type=int,
        help='Width of the output image in pixels. '
             'Default: 90'
    )

    parser.add_argument(
        '-v', '--region-viewpoint', dest='region_viewpoint',
        default='center',
        help='Viewpoint relative to region when using -w. '
             'By default, this measures the window from the '
             'region center. You can change this to other locations '
             'within each region using this parameter. Possible values:'
             'start, end, five_prime, three_prime, center'
    )

    parser.add_argument(
        '-b', '--boundary-mode', dest='boundary_mode',
        default='reflect',
        help='Points outside the boundaries of the input '
             'are filled according to the given mode. Options are'
             'constant, edge, symmetrix, reflect, and warp.'
             'Default: reflect.'
    )

    parser.add_argument(
        '-i', '--interpolation', dest='interpolation',
        type=int,
        default=0,
        help='Type of interpolation to use for resizing. '
             '0: Nearest-neighbor (default), 1: Bi-linear, '
             '2: Bi-quadratic, 3: Bi-cubic, 4: Bi-quartic, '
             '5: Bi-quintic'
    )

    parser.add_argument(
        '-r', '--relative', dest='relative',
        type=float,
        help='Relative extension of each region as fraction '
             'of region length (l). Final region in the '
             'image will be: <start of region - e*l> to '
             '<end of region + e*l>. Default: 1.0 (results '
             'in 3 times region size image). Additive '
             'with "-a" parameter!'
    )

    parser.add_argument(
        '-a', '--absolute', dest='absolute',
        type=int,
        help='Extension (e) of each region in base pairs. '
             'Final region in the image will be:  <start '
             'of TAD - e> to <end of TAD + e>. Default: 0 '
             '(no extension). Additive with "-r" parameter!'
    )

    parser.add_argument(
        '-e', '--expected-norm', dest='oe',
        action='store_true',
        default=False,
        help='Normalize matrix to expected values'
    )

    parser.add_argument(
        '-l', '--log', dest='log',
        action='store_true',
        default=False,
        help='log2-transform normalized matrices. '
             'Only used in conjunction with "-e".'
    )

    parser.add_argument(
        '--rescale', dest='rescale',
        action='store_true',
        default=False,
        help='Rescale normalized contact matrices '
             'using an a=-0.25 power law. '
             'Only used in conjunction with "-e".'
    )

    parser.add_argument(
        '--colormap', dest='colormap',
        help='Matplotlib colormap to use for matrix'
    )

    parser.add_argument(
        '--vmin', dest='vmin',
        type=float,
        help='Minimum saturation value in image'
    )

    parser.add_argument(
        '--vmax', dest='vmax',
        type=float,
        help='Maximum saturation value in image'
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        default=False,
        help='Work in temporary directory'
    )

    parser.add_argument(
        '-C', '--no-cache', dest='cache',
        action='store_false',
        default=True,
        help='Do not cache chromosome matrices. '
             'Slower, but saves a lot of memory. '
             'Use this if you are having trouble '
             'with memory usage.'
    )

    parser.add_argument(
        '--keep-submatrices', dest='keep_submatrices',
        action='store_true',
        default=False,
        help='Save all the individual matrices that make '
             'up the aggregate matrix in the output object. '
             'Useful for debugging and downstream processing. '
             'Potentially uses a lot of memory and/or disk space.'
    )

    parser.add_argument(
        '-s', '--orient-by-strand', dest='orient_strand',
        action='store_true',
        default=False,
        help='Flip submatrix if region is on the negative strand.'
    )
    
    parser.add_argument(
        '--lower-triangular-plot', dest='lower_triangular',
        action='store_true',
        default=False,
        help='Invert the aggregate plot so it corresponds to '
             'the lower triangle of the Hi-C matrix.'
    )

    parser.add_argument(
        '--labels', dest='labels',
        help='Labels for the left, center, and right edge of the '
             'matrix (comma-separated).'
    )

    parser.add_argument(
        '--label-locations', dest='label_locations',
        default="0,0.5,1",
        help='Relative location of ticks on bottom and left of aggregate plot '
             '(comma-separated). '
             'Ranges from 0 (left/bottom) to 1.0 (right/top). '
             'Default: 0,0.5,1.0'
    )

    return parser


def aggregate(argv, **kwargs):
    parser = aggregate_parser()

    args = parser.parse_args(argv[2:])
    import os

    input_file = os.path.expanduser(args.input)
    regions_file = os.path.expanduser(args.regions) if args.regions is not None else None
    output_file = os.path.expanduser(args.output) if args.output is not None else None
    matrix_file = os.path.expanduser(args.matrix_file) if args.matrix_file is not None else None
    plot_file = os.path.expanduser(args.plot_file) if args.plot_file is not None else None
    loop_strength_file = os.path.expanduser(args.loop_strength_file) if args.loop_strength_file is not None else None
    tad_strength_file = os.path.expanduser(args.tad_strength_file) if args.tad_strength_file is not None else None
    tads_preset = args.tads_preset
    tads_imakaev_preset = args.tads_imakaev_preset
    loops_preset = args.loops_preset
    window = args.window
    pixels = args.pixels
    interpolation = args.interpolation
    boundary_mode = args.boundary_mode
    relative = args.relative
    absolute = args.absolute
    oe = args.oe
    log = args.log
    rescale = args.rescale
    colormap = args.colormap
    vmin = args.vmin
    vmax = args.vmax
    cache = args.cache
    keep_submatrices = args.keep_submatrices
    region_viewpoint = args.region_viewpoint
    orient_strand = args.orient_strand
    labels = args.labels if args.labels is None else args.labels.split(",")
    label_locations = [float(loc) for loc in args.label_locations.split(",")]
    lower_triangular = args.lower_triangular
    tmp = args.tmp

    presets = sum([tads_preset, tads_imakaev_preset, loops_preset])
    if presets > 1:
        parser.error("--tads, --tads-imakaev, and --loops are mutually exclusive!")

    if labels is not None:
        if len(labels) != len(label_locations):
            parser.error("Number of labels ({}) must be the same as the number of ticks ({})".format(
                len(labels), len(label_locations)
            ))

    if tads_preset:
        if relative is None:
            relative = 1.0
        oe = True
        log = True
        rescale = False
        if vmin is None:
            vmin = -1
        if vmax is None:
            vmax = 1

    if tads_imakaev_preset:
        relative = 1.0
        oe = True
        log = False
        rescale = True

    if loops_preset:
        oe = True
        log = True
        rescale = False
        if vmin is None:
            vmin = -1
        if vmax is None:
            vmax = 1
        if pixels is None:
            pixels = 16

    if colormap is None:
        if oe and not rescale:
            if log:
                colormap = 'RdBu_r'
            else:
                colormap = 'Reds'
        else:
            colormap = 'germany'

    if pixels is None:
        pixels = 90

    import fanc
    import numpy as np
    import genomic_regions as gr
    import warnings
    from fanc.architecture.aggregate import AggregateMatrix, loop_strength, tad_strength
    from fanc.tools.general import human_format, str_to_int

    if window is not None:
        window = str_to_int(window)

    aggregate_matrix = None
    regions = None
    try:
        with fanc.load(input_file, mode='r') as matrix:
            if not isinstance(matrix, AggregateMatrix):
                regions = fanc.load(regions_file)

                b = matrix.bin_size

                if isinstance(regions, gr.Bedpe):
                    logger.info("Detected BEDPE. Running pairwise region extraction")

                    region_pairs = []
                    for region in regions.regions:
                        a1 = gr.GenomicRegion(chromosome=region.chromosome1,
                                              start=region.start1, end=region.end1)
                        a2 = gr.GenomicRegion(chromosome=region.chromosome2,
                                              start=region.start2, end=region.end2)
                        region_pairs.append((a1, a2))

                    aggregate_matrix = AggregateMatrix.from_center_pairs(matrix, region_pairs,
                                                                         window=window, pixels=pixels,
                                                                         keep_components=keep_submatrices,
                                                                         file_name=output_file, tmpdir=tmp,
                                                                         oe=oe, log=log,
                                                                         orient_strand=orient_strand,
                                                                         cache=cache,
                                                                         region_viewpoint=region_viewpoint)

                    if labels is None:
                        left = int(pixels / 2)
                        right = left if pixels % 2 == 1 else left - 1
                        labels = ['-{}b'.format(human_format(left * b)), '',
                                  '+{}b'.format(human_format(right * b))]

                    if loop_strength_file is not None:
                        loop_strengths = loop_strength(matrix, region_pairs)
                        with open(loop_strength_file, 'w') as o:
                            for s, (r1, r2) in zip(loop_strengths, region_pairs):
                                o.write("{}\t{}\t{}\t{}\t{}\t{}\t.\t{}\n".format(
                                    r1.chromosome, r1.start, r1.end,
                                    r2.chromosome, r2.start, r2.end,
                                    s
                                ))

                elif isinstance(regions, gr.RegionBased):
                    if window is not None:
                        logger.info("Creating aggregate matrix with fixed window size of {}".format(window))
                        aggregate_matrix = AggregateMatrix.from_center(matrix, regions.regions,
                                                                       window=window, rescale=rescale,
                                                                       keep_components=keep_submatrices,
                                                                       file_name=output_file, tmpdir=tmp,
                                                                       oe=oe, log=log,
                                                                       orient_strand=orient_strand,
                                                                       cache=cache,
                                                                       region_viewpoint=region_viewpoint)

                        if labels is None:
                            wh = int(window / 2)
                            labels = ['-{}b'.format(human_format(-wh)), '',
                                      '+{}b'.format(human_format(wh))]

                        pixels = int(window/b)
                        if pixels % 2 == 0:
                            pixels += 1
                    else:
                        logger.info("Creating aggregate matrix from differently-sized regions")
                        aggregate_matrix = AggregateMatrix.from_regions(matrix, regions.regions,
                                                                        pixels=pixels, rescale=rescale,
                                                                        interpolation=interpolation,
                                                                        boundary_mode=boundary_mode,
                                                                        absolute_extension=absolute,
                                                                        relative_extension=relative,
                                                                        keep_components=keep_submatrices,
                                                                        file_name=output_file,
                                                                        tmpdir=tmp,
                                                                        oe=oe, log=log,
                                                                        orient_strand=orient_strand,
                                                                        cache=cache)

                        if tad_strength_file is not None:
                            logger.info("Calculating TAD strength")
                            tad_strengths = tad_strength(matrix, regions.regions)
                            with open(tad_strength_file, 'w') as o:
                                for s, r in zip(tad_strengths, regions.regions):
                                    o.write("{}\t{}\t{}\t.\t{}\n".format(
                                        r.chromosome, r.start, r.end,
                                        s
                                    ))
            else:
                aggregate_matrix = matrix

            if matrix_file is not None:
                np.savetxt(matrix_file, aggregate_matrix.matrix())

            if plot_file is not None:
                import matplotlib
                matplotlib.use('agg')
                import matplotlib.pyplot as plt
                from fanc.plotting.statistics import aggregate_plot

                fig, ax = plt.subplots()
                aggregate_plot(aggregate_matrix, labels=labels, vmin=vmin, vmax=vmax,
                               oe=oe, log=log, colormap=colormap, ax=ax,
                               relative_label_locations=label_locations,
                               lower_triangular=lower_triangular)
                fig.savefig(plot_file)
                plt.close(fig)
    finally:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if aggregate_matrix is not None:
                aggregate_matrix.close()


def stats_parser():
    parser = argparse.ArgumentParser(
        prog="fanc stats",
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


def stats(argv, **kwargs):
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
                line_count = sum(1 for _ in f)
            total_count += line_count/4

            with open(output_file, 'a') as o:
                o.write("fastq\t{}\tcount\t{}\n".format(fastq_file, line_count/4))

        with open(output_file, 'a') as o:
            o.write("fastq\ttotal\tcount\t{}\n".format(total_count))

    # 2. Pairs statistics
    if args.pairs is not None:
        logger.info("Processing Pairs files.")
        import fanc
        pairs_files = get_files(args.pairs, ('.pairs',))
        pairs_summary = defaultdict(int)
        for pairs_file in pairs_files:
            logger.info("{}".format(pairs_file))
            pairs = fanc.load(pairs_file, mode='r')
            statistics, total = stats(pairs, pairs._pairs)

            with open(output_file, 'a') as o:
                for key in sorted(statistics.keys()):
                    o.write("pairs\t{}\t{}\t{}\n".format(pairs_file, key, statistics[key]))
                    pairs_summary[key] += statistics[key]

            with open(output_file, 'a') as o:
                o.write("pairs\t{}\ttotal\t{}\n".format(pairs_file, total))
                pairs_summary['total'] += total

            pairs_summary['filtered'] += total - statistics['valid']
            pairs_summary['remaining'] += statistics['valid']

        with open(output_file, 'a') as o:
            for key in sorted(pairs_summary.keys()):
                if key != 'filtered' and key != 'remaining':
                    o.write("pairs\ttotal\t{}\t{}\n".format(key, pairs_summary[key]))
            o.write("pairs\ttotal\tfiltered\t{}\n".format(pairs_summary['filtered']))
            o.write("pairs\ttotal\tremaining\t{}\n".format(pairs_summary['remaining']))

    # 3. Hic statistics
    if args.hic is not None:
        logger.info("Processing Hic files.")
        import fanc
        hic_files = get_files(args.hic, ('.hic',))

        hic_summary = defaultdict(int)
        for hic_file in hic_files:
            logger.info("{}".format(hic_file))
            hic = fanc.load(hic_file, mode='r')
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
        prog="fanc write_config",
        description='Write default config file to specified location.'
    )

    parser.add_argument(
        'config_file',
        nargs='?',
        help="Output file for default configuration."
    )

    parser.add_argument(
        '-f', '--force', dest='force',
        action='store_true',
        help='''Force overwrite of existing config file.'''
    )
    parser.set_defaults(force=False)
    return parser


def write_config(argv, **kwargs):
    parser = write_config_parser()

    args = parser.parse_args(argv[2:])
    file_name = args.config_file

    from fanc.config import write_default_config
    write_default_config(file_name, overwrite=args.force)


def cis_trans_parser():
    parser = argparse.ArgumentParser(
        prog="fanc cis_trans",
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


def cis_trans(argv, **kwargs):
    parser = cis_trans_parser()

    args = parser.parse_args(argv[2:])
    hic_files = [os.path.expanduser(f) for f in args.hic]
    output_file = os.path.expanduser(args.output) if args.output is not None else None
    normalise = args.normalise

    import fanc
    from fanc.architecture.stats import cis_trans_ratio

    if output_file:
        with open(output_file, 'w') as o:
            o.write("file\tcis\ttrans\tratio\tfactor\n")

    for hic_file in hic_files:
        hic = fanc.load(hic_file, mode='r')

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


def downsample_parser():
    parser = argparse.ArgumentParser(
        prog="fanc downsample",
        description='Downsample contacts from a Hic object.'
    )

    parser.add_argument(
        'hic',
        help="Hic object to be downsampled."
    )

    parser.add_argument(
        'n',
        help="Sample size or reference Hi-C object. If sample size is < 1,"
             "will be interpreted as a fraction of valid pairs."
    )

    parser.add_argument(
        'output',
        help="Downsampled Hic output."
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        default=False,
        help='Work in temporary directory'
    )

    return parser


def downsample(argv, **kwargs):
    parser = downsample_parser()

    print("*** fanc downsample is deprecated. Please use fanc hic --downsample instead! ***")

    args = parser.parse_args(argv[2:])

    import os
    import shutil
    import fanc
    from genomic_regions.files import create_temporary_copy, create_temporary_output

    hic_file = args.hic
    tmp = args.tmp
    n = args.n
    output_file = args.output

    original_output_file = None
    tmp_files = []
    try:

        if os.path.exists(os.path.expanduser(n)):
            if tmp:
                tmp = False
                n = create_temporary_copy(n)
                tmp_files.append(n)
                tmp = True
            n = fanc.load(n)
        else:
            n = float(n)

        if tmp:
            tmp = False
            hic_file = create_temporary_copy(hic_file)
            tmp_files.append(hic_file)
            original_output_file = output_file
            output_file = create_temporary_output(original_output_file)
            tmp_files.append(output_file)
            tmp = True

        with fanc.load(hic_file) as hic:
            output_hic = hic.downsample(n, file_name=output_file)
            output_hic.close()
    finally:
        if original_output_file is not None:
            shutil.copy(output_file, original_output_file)

        for tmp_file in tmp_files:
            os.remove(tmp_file)


def upgrade_parser():
    parser = argparse.ArgumentParser(
        prog="fanc upgrade",
        description='Upgrade objects from old FAN-C versions.'
    )

    parser.add_argument(
        'hic',
        help="Hic object to be upgraded."
    )

    parser.add_argument(
        'output',
        nargs='?',
        help="Output file. If omitted, will try to perform upgrade in place."
    )

    parser.add_argument(
        '-f', '--force', dest='force',
        action='store_true',
        default=False,
        help='''Force upgrade even if object can be loaded.'''
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        default=False,
        help='Work in temporary directory'
    )

    return parser


def upgrade(argv, **kwargs):
    parser = upgrade_parser()

    args = parser.parse_args(argv[2:])

    input_file = os.path.expanduser(args.hic)
    output_file = os.path.expanduser(args.output) if args.output is not None else None
    force = args.force
    tmp = args.tmp

    import fanc
    import numpy as np
    from fanc.hic import LegacyHic
    import tables
    from fanc.tools.general import RareUpdateProgressBar

    original_output_file = None
    tmp_files = []
    try:
        if tmp:
            tmp = False
            from genomic_regions.files import create_temporary_copy, create_temporary_output
            input_file = create_temporary_copy(input_file)
            tmp_files.append(input_file)
            original_output_file = output_file
            output_file = create_temporary_output(output_file)
            tmp = True

        if not force:
            try:
                f = fanc.load(input_file)
                f.close()

                parser.error("Can already load input file. "
                             "There does not seem to be a need for an upgrade? "
                             "You can force an upgrade with -f")
            except (ValueError, TypeError):
                pass

        try:
            old_fanc = fanc.load(input_file, mode='a')

            if isinstance(old_fanc, LegacyHic):
                import fanc

                if output_file is None:
                    raise ValueError("Must provide an output file to upgrade to latest Hi-C version!")

                new_hic = fanc.Hic(file_name=output_file, mode='w')
                new_hic.add_regions(old_fanc.regions)

                bv = [row['bias']**2 for row in old_fanc.file.get_node('/', 'node_annot').iterrows()]

                with RareUpdateProgressBar(prefix="Upgrade", max_value=len(old_fanc.edges)) as pb:
                    for i, edge in enumerate(old_fanc.edges(lazy=True)):
                        source = edge.source
                        sink = edge.sink
                        weight = int(np.round(edge.weight / bv[source] / bv[sink]))
                        new_hic.add_edge_simple(source, sink, weight=weight)
                        pb.update(i)
                new_hic.flush(update_mappability=False)
                new_hic.bias_vector(np.sqrt(bv))

                old_fanc.close()
                new_hic.close()
                return
            elif isinstance(old_fanc, fanc.ReadPairs) or isinstance(old_fanc, fanc.Hic):
                if old_fanc._chromosomes_info is None:
                    if output_file is not None:
                        raise ValueError("Can only do upgrade in place for chromosome table")
                    logger.info("Updating chromosome info table")
                    old_fanc._update_chromosomes_info()
                    old_fanc.close()
                return
        except (ValueError, TypeError):
            pass

        file_based = fanc.FileBased(input_file, mode='r')
        class_id = file_based.meta._classid
        logger.info("Detected class ID '{}'".format(class_id))

        target_class = None
        if class_id == 'ACCESSOPTIMISEDHIC' or class_id == 'HIC':
            target_class = fanc.Hic
        elif class_id == 'RAOPEAKINFO':
            target_class = fanc.RaoPeakInfo
        else:
            parser.error("No suitable upgrade method for {} - "
                         "please consult the developer!".format(class_id))

        bias = [row['bias'] for row in file_based.file.get_node('/', 'node_annot').iterrows()]

        regions = []
        nodes_table = file_based.file.get_node('/', 'nodes')
        try:
            region_fields = nodes_table.coldescrs
        except tables.NoSuchNodeError:
            nodes_table = nodes_table.regions
            region_fields = nodes_table.coldescrs

        for i, row in enumerate(nodes_table.iterrows()):
            kwargs = {name: row[name] for name in region_fields.keys()}
            kwargs['bias'] = bias[i]
            r = fanc.GenomicRegion(**kwargs)
            regions.append(r)

        edges_table = file_based.file.get_node('/', 'edges')
        if isinstance(edges_table, tables.Table):
            edge_fields = edges_table.coldescrs
        else:
            edge_fields = edges_table.chrpair_0_0.coldescrs

        upgraded_hic = target_class(output_file, mode='w',
                                    additional_region_fields=region_fields,
                                    additional_edge_fields=edge_fields)
        upgraded_hic.add_regions(regions, preserve_attributes=False)

        if isinstance(edges_table, tables.Table):
            for row in edges_table.iterrows():
                kwargs = {name: row[name] for name in edge_fields.keys()}
                source = kwargs['source']
                sink = kwargs['sink']
                weight = kwargs['weight']
                # uncorrect matrix
                uncorrected = weight / bias[source] / bias[sink]
                kwargs['weight'] = uncorrected
                edge = fanc.Edge(**kwargs)
                upgraded_hic.add_edge(edge)
        else:
            for table in edges_table:
                for row in table.iterrows():
                    kwargs = {name: row[name] for name in edge_fields.keys()}
                    source = kwargs['source']
                    sink = kwargs['sink']
                    weight = kwargs['weight']
                    # uncorrect matrix
                    uncorrected = weight / bias[source] / bias[sink]
                    kwargs['weight'] = uncorrected
                    edge = fanc.Edge(**kwargs)
                    upgraded_hic.add_edge(edge)
        upgraded_hic.flush()
    finally:
        if tmp:
            if output_file is not None and original_output_file is not None:
                shutil.copy(output_file, original_output_file)
                os.remove(output_file)

            for file_name in tmp_files:
                os.remove(file_name)
