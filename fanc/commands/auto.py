import os
import argparse
import logging
import subprocess
import shlex
import uuid
import threading
import multiprocessing as mp
import tempfile
from future.utils import string_types
from fanc.config import config
import re
import warnings


# configure logging
logger = logging.getLogger(__name__)


class CommandTask(object):
    def __init__(self, command):
        if isinstance(command, string_types):
            self.command = shlex.split(command, posix=False)
        else:
            self.command = command
        self.id = uuid.uuid4()


class QueuedTask(object):
    def __init__(self, task, wait_for=None, threads=1):
        self.task = task
        self.id = task.id
        self.command = task.command

        if isinstance(wait_for, CommandTask):
            wait_for = [wait_for]
        elif wait_for is None:
            wait_for = []

        self.wait_for = wait_for
        self.threads = threads


class TaskRunner(object):
    def __init__(self):
        self._tasks = []
        self._task_ixs = dict()

    def add_task(self, task, wait_for=None, threads=1, *args, **kwargs):
        self._tasks.append(QueuedTask(task, wait_for=wait_for, threads=threads))
        self._task_ixs[task.id] = len(self._task_ixs)

    def run(self, *args, **kwargs):
        raise NotImplementedError("Subclasses of TaskWorker must implement 'run'")


def _run_task(task, process_id, completed_queue):
    res = subprocess.call(task.command)
    completed_queue.put((task, process_id, res))


class ParallelTaskRunner(TaskRunner):
    def __init__(self, threads=1, test=False):
        TaskRunner.__init__(self)
        self._threads = threads
        self._active_threads_lock = threading.Lock()
        self._active_threads = 0
        self._test = test

    def _change_active_threads(self, threads):
        with self._active_threads_lock:
            self._active_threads += threads

    def run(self, *args, **kwargs):
        completed = set()
        running = set()
        completed_queue = mp.Queue()
        processes = dict()

        remaining_tasks = len(self._tasks)
        while remaining_tasks > 0:

            # find tasks that can be executed
            for task in self._tasks:
                if task.id in completed:
                    continue

                if task.id in running:
                    continue

                if task.threads > self._threads:
                    raise RuntimeError("Task requests more threads ({}) than "
                                       "available in parallel task runner ({})!".format(task.threads,
                                                                                        self._threads))

                available_threads = self._threads - self._active_threads
                if task.threads > available_threads:
                    continue

                is_executable = True
                for wait_task in task.wait_for:
                    if not wait_task.id in completed:
                        is_executable = False
                        break

                if not is_executable:
                    continue

                p_id = uuid.uuid4()
                if not self._test:
                    logger.debug("Running task {}".format(task.id))
                    p = mp.Process(target=_run_task, args=(task, p_id, completed_queue))
                    p.start()
                    processes[p_id] = p
                else:
                    print("{} (threads: {}, depends on: {}): {}".format(
                        self._task_ixs[task.id],
                        task.threads,
                        ", ".join([str(self._task_ixs[t.id])
                                   for t in task.wait_for]),
                        " ".join(task.command))
                    )

                    completed_queue.put((task, p_id, 0))
                running.add(task.id)
                self._change_active_threads(task.threads)

            # wait for any task to finish
            task, p_id, res = completed_queue.get(block=True)
            completed.add(task.id)
            running.remove(task.id)
            self._change_active_threads(-task.threads)
            if not self._test:
                processes[p_id].join()
                del processes[p_id]

            if res != 0:
                for p in processes.values():
                    p.terminate()
                    p.join()
                raise RuntimeError("Task {} had non-zero exit status. "
                                   "Cancelling execution of all tasks.".format(task.id))
            else:
                remaining_tasks -= 1


class SgeTaskRunner(TaskRunner):
    def __init__(self, task_prefix=None,
                 log_dir=None, trap_sigusr=True,
                 startup_commands_file=None,
                 cleanup_commands_file=None,
                 buffer_python_output=False):
        TaskRunner.__init__(self)
        if task_prefix is None:
            from fanc.tools.files import random_name
            self._task_prefix = 'fanc_' + random_name(6) + '_'
        else:
            self._task_prefix = task_prefix

        self._log_dir = log_dir
        if log_dir is None:
            self._log_dir = config.sge_log_dir
        else:
            self._log_dir = os.path.expanduser(log_dir)
        self._trap_sigusr = trap_sigusr
        self._startup_commands_file = startup_commands_file
        self._cleanup_commands_file = cleanup_commands_file
        self._buffer_python_output = buffer_python_output

    def _submit_task(self, task, kill=True):
        task_ix = self._task_ixs[task.id]
        job_id = self._task_prefix + '{}'.format(task_ix)
        with tempfile.NamedTemporaryFile('w', prefix='fanc_auto_', suffix='.sh') as tmp_file:
            if self._trap_sigusr:
                tmp_file.write('function notify_handler() {\n  '
                               '(>&2 echo "Received termination notice")\n}\n')
                tmp_file.write("trap notify_handler SIGUSR1\n")
                tmp_file.write("trap notify_handler SIGUSR2\n\n")
            if not self._buffer_python_output:
                tmp_file.write("export PYTHONUNBUFFERED=1\n\n")

            if self._startup_commands_file is not None:
                with open(self._startup_commands_file) as f:
                    startup_commands = f.read()
                tmp_file.write(startup_commands)
                tmp_file.write("\n")

            tmp_file.write(" ".join(task.command) + "\n")
            tmp_file.write("OUT=$?\n")
            tmp_file.write("if [ $OUT -ne 0 ]; then\n")
            tmp_file.write('    (>&2 echo "Job $JOB_ID / $JOB_NAME had non-zero exit status")\n')
            if kill is not None:
                if self._log_dir is not None:
                    tmp_file.write(
                        '    echo "Failed job ID: {}" >> {}\n'.format(
                            job_id, os.path.join(self._log_dir, self._task_prefix + "FAILED")
                        ))
                tmp_file.write('    {} "{}*"\n'.format(config.sge_qdel_path, self._task_prefix))
            tmp_file.write('    res="failed"\n')
            tmp_file.write('else\n')
            tmp_file.write('    res="succeeded"\n')
            tmp_file.write("fi\n\n")

            if self._cleanup_commands_file is not None:
                with open(self._cleanup_commands_file) as f:
                    cleanup_commands = f.read()
                tmp_file.write(cleanup_commands)
                tmp_file.write("\n")

            tmp_file.flush()

            command = [config.sge_qsub_path,
                       '-N', job_id, '-cwd',
                       '-pe', config.sge_parallel_environment,
                       str(task.threads)] + config.sge_qsub_options.split()

            if config.sge_default_queue is not None:
                command += ['-q', config.sge_default_queue]

            if config.sge_shell is None:
                shell = '/bin/bash'
            else:
                shell = config.sge_shell
            command += ['-S', shell]

            if task.wait_for is not None and len(task.wait_for) > 0:
                hold_ids = ",".join([self._task_prefix + '{}'.format(self._task_ixs[t.id])
                                     for t in task.wait_for])
                command += ['-hold_jid', hold_ids]

            if self._log_dir is not None:
                command += ['-o', os.path.join(self._log_dir, job_id + '_o')]
                command += ['-e', os.path.join(self._log_dir, job_id + '_e')]
            else:
                command += ['-o', '/dev/null']
                command += ['-e', '/dev/null']

            command += [tmp_file.name]

            logger.info("Submitting {}".format(" ".join(command)))
            subprocess.call(command)

    def run(self, *args, **kwargs):
        for task in self._tasks:
            self._submit_task(task)


class SlurmTaskRunner(TaskRunner):
    def __init__(self, task_prefix=None,
                 log_dir=None, trap_sigusr=True,
                 startup_commands_file=None,
                 cleanup_commands_file=None,
                 buffer_python_output=False):
        TaskRunner.__init__(self)
        if task_prefix is None:
            from fanc.tools.files import random_name
            self._task_prefix = 'fanc_' + random_name(6) + '_'
        else:
            self._task_prefix = task_prefix

        self._log_dir = log_dir
        if log_dir is None:
            self._log_dir = config.slurm_log_dir
        else:
            self._log_dir = os.path.expanduser(log_dir)
        self._trap_sigusr = trap_sigusr

        self._slurm_job_ids = dict()
        self._startup_commands_file = startup_commands_file
        self._cleanup_commands_file = cleanup_commands_file
        self._buffer_python_output = buffer_python_output

    def _submit_batch_command(self, task, kill=True):
        job_name = self._task_prefix
        with tempfile.NamedTemporaryFile('w', prefix='fanc_auto_', suffix='.sh') as tmp_file:
            tmp_file.write("#!{}\n".format(config.slurm_shell if config.slurm_shell is not None else "/bin/bash"))
            tmp_file.write("#SBATCH --nodes=1\n")
            tmp_file.write("#SBATCH --ntasks=1\n")
            tmp_file.write("#SBATCH --cpus-per-task={}\n".format(task.threads))
            tmp_file.write("#SBATCH --job-name={}\n".format(job_name))
            if self._log_dir is not None:
                task_ix = self._task_ixs[task.id]
                log_prefix = self._task_prefix + str(task_ix)
                tmp_file.write("#SBATCH --output={}\n".format(os.path.join(self._log_dir, log_prefix + '_o')))
                tmp_file.write("#SBATCH --error={}\n".format(os.path.join(self._log_dir, log_prefix + '_e')))
            else:
                tmp_file.write("#SBATCH --output={}\n".format('/dev/null'))
                tmp_file.write("#SBATCH --error={}\n".format('/dev/null'))

            if task.wait_for is not None and len(task.wait_for) > 0:
                wait_job_ids = []
                for t in task.wait_for:
                    try:
                        wait_job_ids.append(self._slurm_job_ids[t.id])
                    except KeyError:
                        raise RuntimeError("Task {} has not been submitted to Slurm yet!".format(t.id))

                dependencies = "afterok:{}".format(":".join(wait_job_ids))
                tmp_file.write("#SBATCH --dependency={}\n".format(dependencies))

            tmp_file.write("#SBATCH --export=ALL\n")

            if not self._buffer_python_output:
                tmp_file.write("export PYTHONUNBUFFERED=1\n\n")

            if self._startup_commands_file is not None:
                with open(self._startup_commands_file) as f:
                    startup_commands = f.read()
                tmp_file.write(startup_commands)
                tmp_file.write("\n")

            if self._trap_sigusr:
                tmp_file.write('function notify_handler() {\n  '
                               '(>&2 echo "Received termination notice")\n}\n')
                tmp_file.write("trap notify_handler SIGUSR1\n")
                tmp_file.write("trap notify_handler SIGUSR2\n\n")

            tmp_file.write(" ".join(task.command) + "\n")
            tmp_file.write("OUT=$?\n")
            tmp_file.write("if [ $OUT -ne 0 ]; then\n")
            tmp_file.write('    (>&2 echo "Job $JOB_ID / $JOB_NAME had non-zero exit status")\n')
            if kill is not None:
                if self._log_dir is not None:
                    tmp_file.write(
                        '    echo "Failed job ID: {}" >> {}\n'.format(
                            job_name, os.path.join(self._log_dir,
                                                   self._task_prefix + "FAILED")
                        ))
                tmp_file.write('    {} --signal=USR1 -n "{}"\n'.format(config.slurm_scancel_path,
                                                                       self._task_prefix))
            tmp_file.write('    res="failed"\n')
            tmp_file.write('else\n')
            tmp_file.write('    res="succeeded"\n')
            tmp_file.write("fi\n\n")

            if self._cleanup_commands_file is not None:
                with open(self._cleanup_commands_file) as f:
                    cleanup_commands = f.read()
                tmp_file.write(cleanup_commands)
                tmp_file.write("\n")

            tmp_file.flush()

            with open(tmp_file.name, 'r') as f:
                for line in f:
                    line = line.rstrip()
                    print(line)

            command = [config.slurm_sbatch_path] + config.slurm_sbatch_options.split()

            command += [tmp_file.name]

            slurm_job_id = subprocess.check_output(command)
            slurm_job_id = slurm_job_id.rstrip()
            if isinstance(slurm_job_id, bytes):
                slurm_job_id = slurm_job_id.decode('utf-8')
            m = re.search(r'(\d+)', slurm_job_id)
            slurm_job_id = m.group(1)

            return slurm_job_id

    def run(self, *args, **kwargs):
        completed = set()

        remaining_tasks = len(self._tasks)
        while remaining_tasks > 0:

            # find tasks that can be executed
            for task in self._tasks:
                if task.id in completed:
                    continue

                is_executable = True
                for wait_task in task.wait_for:
                    if wait_task.id not in completed:
                        is_executable = False
                        break

                if not is_executable:
                    continue

                # all dependencies submitted
                self._slurm_job_ids[task.id] = self._submit_batch_command(task)
                completed.add(task.id)
                remaining_tasks -= 1


def auto_parser():
    parser = argparse.ArgumentParser(
        prog="fanc auto",
        description='Automatically process an entire Hi-C data set.'
    )

    parser.add_argument(
        'input',
        nargs='+',
        help="Input files. "
             "fanc will try to guess the file "
             "type by its extension."
    )

    parser.add_argument(
        'output_folder',
        help="Output folder. "
             "All output files and folders " 
             "will be generated under this directory."
    )

    parser.add_argument(
        '-g', '--genome', dest='genome',
        help="Genome for the Hi-C object." 
             "Path to region-based file (BED, GFF, ...) containing "
             "the non-overlapping regions to be used for Hi-C "
             "object construction. Typically restriction-enzyme fragments. "
             "Alternatively: Path to genome file (FASTA, folder with "
             "FASTA, hdf5 file), which will be used in conjunction "
             "with the type of restriction enzyme (-r) to calculate "
             "fragments directly."
    )

    parser.add_argument(
        '-r', '--restriction-enzyme', dest='restriction_enzyme',
        help="Restriction enzyme name. "
             "Used for in silico digestion "
             "of genomic sequences and splitting of reads at Hi-C "
             "ligation junctions. (e.g. HindIII, case-sensitive). "
             "Separate multiple enzymes with ','. "
             "Restriction names can be any supported by Biopython, which obtains data "
             "from REBASE (http://rebase.neb.com/rebase/rebase.html). "
    )

    parser.add_argument(
        '-i', '--genome-index', dest='genome_index',
        help="Bowtie 2 or BWA genome index. "
             "Only required when passing FASTQ "
             "files as input."
    )

    parser.add_argument(
        '-n', '--basename', dest='basename',
        help="Basename for output files. " 
             "If not provided, will be guessed based " 
             "on input file names."
    )

    parser.add_argument(
        '-s', '--step-size', dest='step_size',
        type=int,
        default=3,
        help="Step size for iterative mapping. " 
             "Default: %(default)d"
    )

    parser.add_argument(
        '-b', '--bin-sizes', dest='bin_sizes',
        nargs='+',
        default=['5mb', '2mb', '1mb', '500kb', '250kb', '100kb', '50kb',
                 '25kb', '10kb', '5kb'],
        help="Bin sizes for Hi-C matrix generation. " 
             "Default: 5mb, 2mb, 1mb, "
             "500kb, 250kb, 100kb, 50kb, 25kb, " 
             "10kb, 5kb."
    )

    parser.add_argument(
        '-t', '--threads', dest='threads',
        type=int,
        default=1,
        help="Maximum number of threads. "
             "The number provided here will not be exceeded "
             "by analysis steps running alone or in parallel. "
             "Default: %(default)d"
    )

    parser.add_argument(
        '--max-restriction-site-distance', dest='max_restriction_site_distance',
        type=int,
        default=10000,
        help="Insert / ligation fragment sizes are inferred from the "
             "sum of distances of both reads to the nearest restriction sites. "
             "Fragment larger than this value are filtered from the Pairs object. "
             "The default value of %(default)d only removes extremely large fragments, "
             "and is thus considered conservative."
             "Default: %(default)d"
    )

    parser.add_argument(
        '--fanc-parallel', dest='mapper_parallel',
        action='store_false',
        default=True,
        help='Use FAN-C parallelisation, which launches multiple mapper jobs. '
             'This may be faster in some cases than relying '
             'on the internal parallelisation of the mapper, '
             'but has potentially high disk I/O and memory usage.'
    )

    parser.add_argument(
        '--split-fastq', dest='split_fastq',
        action='store_true',
        default=False,
        help="Split fastq files into chunks of 10M reads. " 
             "Reads will be merged again on the SAM level. "
             "Splitting and merging bypasses the -tmp flag. " 
             "This option reduces disk usage in tmp, in case "
             "the system has a small tmp partition. "
    )

    parser.add_argument(
        '--memory-map', dest='memory_map',
        action='store_true',
        default=False,
        help="Map Bowtie2 index to memory. Recommended " 
             "if running on medium-memory systems and using many "
             "parallel threads)."
    )

    parser.add_argument(
        '--ice', dest='ice',
        action='store_true',
        default=False,
        help="DEPRECATED. Correct Hi-C matrices using ICE instead of Knight-Ruiz "
             "matrix balancing. Slower, but much more memory-friendly."
    )

    parser.add_argument(
        '--norm-method', dest='norm_method',
        default='kr',
        help='Normalisation method. Options are: '
             'KR (default) = Knight-Ruiz matrix balancing '
             '(Fast, accurate, but memory-intensive); '
             'ICE = ICE matrix balancing (more CPU-intensive, but also more memory-efficient); '
             'VC = vanilla coverage (a single round of ICE balancing); '
             'VC-SQRT = vanilla coverage square root (reduces overcorrection compared to VC)'
    )

    parser.add_argument(
        '-q', '--quality-cutoff', dest='quality_cutoff',
        type=float,
        help='Cutoff for the minimum mapping quality of a read. '
             'For numbers larger than 1, will filter on MAPQ. '
             'If a number between 0 and 1 is provided, will filter '
             'on the AS tag instead of mapping quality (only BWA). '
             'The quality cutoff is then interpreted as the '
             'fraction of bases that have to be matched for any '
             'given read. Only applies to SAM/BAM input!. '
             'Default is not to filter on mapping quality.'
    )

    parser.add_argument(
        '--iterative-quality-cutoff', dest='iterative_quality_cutoff',
        type=int,
        help='MAPQ cutoff for mapped reads. Only applies when iterative '
             'mapping is enabled: if a mapped read has MAPQ below this cutoff,'
             'it will be sent to another iteration in an attempt to find a '
             'higher quality alignment. Default is 3 for BWA and 30 for Bowtie2.'
    )

    parser.add_argument(
        '--le-inward-cutoff', dest='inward_cutoff',
        type=int,
        help="Ligation error inward cutoff. "
             "Default: no ligation error filtering."
    )

    parser.add_argument(
        '--le-outward-cutoff', dest='outward_cutoff',
        type=int,
        help="Ligation error outward cutoff. "
             "Default: no ligation error filtering."
    )

    parser.add_argument(
        '--auto-le-cutoff', dest='auto_le_cutoff',
        action='store_true',
        default=False,
        help="Automatically determine ligation error cutoffs. "
             "Use with caution, this setting has a tendency to "
             "choose large cutoffs, removing many pairs close to " 
             "the diagonal."
    )

    parser.add_argument(
        '-tmp', '--work-in-tmp', dest='tmp',
        action='store_true',
        default=False,
        help="Work in temporary directory. Copies input files and " 
             "generates output files in a temporary directory. "
             "Files will be moved to their intended destination once " 
             "an analysis step finishes. "
             "Reduces network I/O if using remote file systems."
    )

    parser.add_argument(
        '--iterative', dest='iterative',
        action='store_true',
        default=False,
        help="Map reads iteratively. Can improve mappability, "
             "especially with low-quality reads. Reads are initially " 
             "trimmed to 25bp and mapped to the reference genome. "
             "If no unique mapping location is found, the read is " 
             "extended by 3bp and the process is repeated until the " 
             "full length of the read is reached or a unique mapping " 
             "location is found."
    )

    parser.add_argument(
        '--no-sam-sort', dest='sam_sort',
        action='store_false',
        default=True,
        help="Do not sort SAM/BAM files. Sorted files are required " 
             "for the pair generating step. Only omit this if you are " 
             "supplying presorted (by read name) SAM/BAM files."
    )

    parser.add_argument(
        '--restore-coverage', dest='restore_coverage',
        action='store_true',
        default=False,
        help="Restore coverage to the original total number of reads. " 
             "Otherwise matrix entries will be contact probabilities. "
             "Only available for KR matrix balancing."
    )

    parser.add_argument(
        '--split-ligation-junction', dest='split_ligation_junction',
        action='store_true',
        default=False,
        help="Split reads at predicted ligation junction before mapping. "
             "Requires the -r argument."
    )

    parser.add_argument(
        '--no-filter-pairs', dest='filter_pairs',
        action='store_false',
        default=True,
        help="Do not filter read pairs. By default, the following "
             "filters are applied: self-ligations, PCR duplicates,"
             "restriction distance (>10kb)"
    )

    parser.add_argument(
        '--no-hic', dest='process_hic',
        action='store_false',
        default=True,
        help="Do not process pairs into Hi-C maps "
             "(stop after read pairing step)."
    )

    parser.add_argument(
        '--run-with', dest='run_with',
        default='parallel',
        help="Choose how to run the commands in fanc auto. Options: "
             "'parallel' (default): Run fanc commands on local machine, "
             "use multiprocessing parallelisation. "
             "'sge': Submit fanc commands to a Sun/Oracle Grid Engine cluster. "
             "'slurm': Submit fanc commands to a Slurm cluster. "
             "'test': Do not run fanc commands but print all commands " 
             "and their dependencies to stdout for review."
    )

    parser.add_argument(
        '--job-prefix', dest='job_prefix',
        help="Job Prefix for SGE and Slurm. "
             "Works with '--run-with sge' and --run-with slurm. "
             "Default: 'fanc_<6 random letters>_'"
    )

    parser.add_argument(
        '--grid-startup-commands', dest='grid_startup_commands',
        help="Path to a file with BASH commands that are executed "
             "before every FAN-C command that is run on a grid engine / cluster. "
             "This could, for example, include environment-specific settings, "
             "such as activation of a Python virtualenv."
    )

    parser.add_argument(
        '--grid-cleanup-commands', dest='grid_cleanup_commands',
        help="Path to a file with BASH commands that are executed "
             "after every FAN-C command that is run on a grid engine / cluster. "
             "Use this to clean the file system or environment set up with "
             "--grid-startup-commands."
    )

    parser.add_argument(
        '-f', '--force-overwrite', dest='force_overwrite',
        action='store_true',
        default=False,
        help="Force overwriting of existing files. Otherwise you will " 
             "be prompted before files are overwritten."
    )

    return parser


def file_type(file_name):
    # try to see if this is a valid pairs file
    from fanc.pairs import HicProPairGenerator, FourDNucleomePairGenerator
    try:
        _ = FourDNucleomePairGenerator(file_name)
        return 'pairs_txt'
    except ValueError:
        pass

    try:
        _ = HicProPairGenerator(file_name)
        return 'pairs_txt'
    except ValueError:
        pass

    base, extension = os.path.splitext(file_name)
    if extension in ['.sam', '.bam']:
        return 'sam'
    if extension in ['.pairs']:
        return 'pairs'
    if extension in ['.hic']:
        return 'hic'
    if extension in ['.gz', '.gzip']:
        base, extension = os.path.splitext(base)
    if extension in ['.fq', '.fastq']:
        return 'fastq'
    return None


def file_basename(file_name):
    basename = os.path.basename(os.path.splitext(file_name)[0])
    if file_name.endswith('.gz') or file_name.endswith('.gzip'):
        basename = os.path.splitext(basename)[0]
    return basename


def auto(argv, **kwargs):
    from fanc.tools.general import which

    parser = auto_parser()
    args = parser.parse_args(argv[2:])

    verbosity = kwargs.get("verbosity", 2)
    if verbosity > 0:
        verbosity_flag = '-' + 'v' * verbosity
        fanc_base_command = ['fanc', verbosity_flag]
    else:
        fanc_base_command = ['fanc']

    log_file = kwargs.get("log_file", None)
    if log_file is not None:
        fanc_base_command += ['-l', log_file]

    bin_sizes = args.bin_sizes
    split_ligation_junction = args.split_ligation_junction
    restriction_enzyme = args.restriction_enzyme
    threads = args.threads
    genome = args.genome
    genome_index = args.genome_index
    basename = args.basename
    quality_cutoff = args.quality_cutoff
    iterative_quality_cutoff = args.iterative_quality_cutoff
    tmp = args.tmp
    mapper_parallel = args.mapper_parallel
    split_fastq = args.split_fastq
    memory_map = args.memory_map
    iterative = args.iterative
    step_size = args.step_size
    sam_sort = args.sam_sort
    filter_pairs = args.filter_pairs
    inward_cutoff = args.inward_cutoff
    outward_cutoff = args.outward_cutoff
    auto_le_cutoff = args.auto_le_cutoff
    process_hic = args.process_hic
    ice = args.ice
    norm_method = args.norm_method
    restore_coverage = args.restore_coverage
    run_with = args.run_with
    job_prefix = args.job_prefix
    max_restriction_site_distance = args.max_restriction_site_distance
    grid_startup_commands = os.path.expanduser(args.grid_startup_commands) \
        if args.grid_startup_commands is not None else None
    grid_cleanup_commands = os.path.expanduser(args.grid_cleanup_commands) \
        if args.grid_cleanup_commands is not None else None
    force_overwrite = args.force_overwrite
    output_folder = os.path.expanduser(args.output_folder)

    file_names = [os.path.expanduser(file_name) for file_name in args.input]
    file_types = [file_type(file_name) for file_name in file_names]
    file_basenames = [file_basename(file_name) for file_name in file_names]

    if ice:
        warnings.warn("The --ice option is deprecated. Please use '--norm-method ice' instead!")
        norm_method = 'ice'

    for file_name in file_names:
        if not os.path.exists(file_name):
            parser.error("File '{}' does not exist!".format(file_name))

    runner = None
    if run_with == 'parallel':
        runner = ParallelTaskRunner(threads)
    elif run_with == 'sge':
        from fanc.config import config
        if which(config.sge_qsub_path) is None:
            parser.error("Using SGE not possible: "
                         "Cannot find 'qsub' at path '{}'. You can change "
                         "this path using fanc config files and the "
                         "'sge_qsub_path' parameter".format(config.sge_qsub_path))
        from fanc.tools.files import mkdir
        sge_log_dir = mkdir(output_folder, 'sge_logs')
        runner = SgeTaskRunner(log_dir=sge_log_dir, task_prefix=job_prefix,
                               startup_commands_file=grid_startup_commands,
                               cleanup_commands_file=grid_cleanup_commands)
    elif run_with == 'slurm':
        from fanc.config import config
        if which(config.slurm_sbatch_path) is None:
            parser.error("Using Slurm not possible: "
                         "Cannot find 'sbatch' at path '{}'. You can change "
                         "this path using fanc config files and the "
                         "'slurm_sbatch_path' parameter".format(config.slurm_sbatch_path))
        from fanc.tools.files import mkdir
        slurm_log_dir = mkdir(output_folder, 'slurm_logs')
        runner = SlurmTaskRunner(log_dir=slurm_log_dir, task_prefix=job_prefix,
                               startup_commands_file=grid_startup_commands,
                               cleanup_commands_file=grid_cleanup_commands)
    elif run_with == 'test':
        runner = ParallelTaskRunner(threads, test=True)
    else:
        parser.error("Runner '{}' is not valid. See --run-with "
                     "parameter for options".format(run_with))

    for i in range(len(file_types)):
        if file_types[i] not in ('fastq', 'sam', 'pairs_txt', 'pairs', 'hic'):
            import fanc
            try:
                ds = fanc.load(file_names[i], mode='r')
                if isinstance(ds, fanc.Hic):
                    file_types[i] = 'hic'
                elif isinstance(ds, fanc.ReadPairs):
                    file_types[i] = 'pairs'
                else:
                    raise ValueError("Could not detect file type using fanc load.")
            except ValueError:
                parser.error("Not a valid input file type: {}".format(file_types[i]))

    if basename is None:
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

    if not output_folder[-1] == '/':
        output_folder += '/'

    # 0. Do some sanity checks on required flags
    is_bwa = False
    is_bowtie2 = False
    if 'fastq' in file_types:
        if args.genome_index is None:
            parser.error("Must provide genome index (-i) when mapping FASTQ files!")
        else:
            check_path = os.path.expanduser(genome_index)
            if check_path.endswith('.'):
                check_path = check_path[:-1]

            is_bowtie2 = True
            for i in range(1, 5):
                if not os.path.exists(check_path + '.{}.bt2'.format(i)):
                    is_bowtie2 = False
            for i in range(1, 3):
                if not os.path.exists(check_path + '.rev.{}.bt2'.format(i)):
                    is_bowtie2 = False

            is_bwa = True
            for ending in ('amb', 'ann', 'bwt', 'pac', 'sa'):
                if not os.path.exists(check_path + '.{}'.format(ending)):
                    is_bwa = False

            if not is_bowtie2 and not is_bwa:
                parser.error("Cannot detect Bowtie2 or BWA index.")

            if is_bowtie2 and not which('bowtie2'):
                parser.error("bowtie2 must be in PATH for mapping!")

            if is_bwa and not which('bwa'):
                parser.error("bwa must be in PATH for mapping!")

    if 'fastq' in file_types or 'sam' in file_types:
        if genome is None:
            parser.error("Must provide genome (-g) to process read pair files!")

        if restriction_enzyme is None:
            from fanc.regions import genome_regions
            try:
                genome_regions(genome)
            except ValueError:
                parser.error("Must provide restriction enzyme (-r) to process read pair files!")
        else:
            res = restriction_enzyme.split(",")
            from Bio import Restriction
            for r in res:
                try:
                    getattr(Restriction, r)
                except AttributeError:
                    parser.error("Restriction enzyme string '{}' is not recognized".format(restriction_enzyme))

    logger.info("Output folder: {}".format(output_folder))
    logger.info("Input files: {}".format(", ".join(file_names)))
    logger.info("Input file types: {}".format(", ".join(file_types)))

    logger.info("Final basename: %s (you can change this with the -n option!)" % basename)

    # 1. create default folders in root directory
    if run_with != 'test':
        logger.info("Creating output folders...")
        from ..tools.general import mkdir
        mkdir(output_folder, 'fastq')
        mkdir(output_folder, 'sam')
        mkdir(output_folder, 'pairs/')
        mkdir(output_folder, 'hic/binned')
        mkdir(output_folder, 'plots/stats')

    # 2. If input files are (gzipped) FASTQ, map them iteratively first
    fastq_files = []
    for i in range(len(file_names)):
        if file_types[i] != 'fastq':
            continue
        fastq_files.append(i)

    sam_created = [False] * len(file_names)
    mapping_tasks = []
    if len(fastq_files) > 0:
        if genome_index.endswith('.'):
            genome_index = genome_index[:-1]

        bam_files = []
        for i, ix in enumerate(fastq_files):
            bam_file = output_folder + 'sam/' + file_basenames[ix] + '.bam'
            if not force_overwrite and os.path.exists(bam_file):
                parser.error("File exists ({}), use -f to force overwriting it.".format(bam_file))
            bam_files.append(bam_file)

            mapping_command = fanc_base_command + ['map', '-m', '25',
                                                   '-s', str(step_size),
                                                   '-t', str(threads)]

            if iterative_quality_cutoff is not None:
                mapping_command += ['-q', str(iterative_quality_cutoff)]
            if tmp:
                mapping_command.append('-tmp')
            if not mapper_parallel:
                mapping_command.append('--fanc-parallel')
            if split_fastq:
                mapping_command.append('--split-fastq')
            if memory_map:
                mapping_command.append('--memory-map')
            if not iterative:
                mapping_command.append('--no-iterative')
            if split_ligation_junction:
                mapping_command.append('--restriction-enzyme')
                mapping_command.append(restriction_enzyme)

            mapping_command += [file_names[ix], genome_index, bam_file]

            mapping_task = CommandTask(mapping_command)
            runner.add_task(mapping_task, threads=threads)
            mapping_tasks.append(mapping_task)

        for ix, i in enumerate(fastq_files):
            file_names[i] = bam_files[ix]
            file_types[i] = 'sam'
            sam_created[i] = True

    if sam_sort:
        sort_threads = min(4, threads)

        sam_sort_tasks = []
        # sort SAM files
        sam_files = []
        in_place = []
        for i in range(len(file_names)):
            if file_types[i] != 'sam':
                continue
            sam_files.append(i)
            in_place.append(sam_created[i])

        if len(sam_files) > 0:
            sorted_sam_files = []
            for i, ix in enumerate(sam_files):
                sort_command = fanc_base_command + ['sort_sam', '-t', str(sort_threads),
                                                    file_names[ix]]
                if in_place[i]:
                    sorted_sam_files.append(file_names[ix])
                else:
                    sam_path, sam_extension = os.path.splitext(file_names[ix])
                    sam_basename = os.path.basename(sam_path)
                    sorted_sam_file = os.path.join(output_folder, 'sam', sam_basename + '_sort' + sam_extension)
                    if not force_overwrite and os.path.exists(sorted_sam_file):
                        parser.error("File exists ({}), use -f to force overwriting it.".format(sorted_sam_file))

                    sorted_sam_files.append(sorted_sam_file)
                    sort_command.append(sorted_sam_file)

                if tmp:
                    sort_command.append('-tmp')

                sam_sort_task = CommandTask(sort_command)
                runner.add_task(sam_sort_task, wait_for=mapping_tasks, threads=1)
                sam_sort_tasks.append(sam_sort_task)

            for ix, i in enumerate(sam_files):
                file_names[i] = sorted_sam_files[ix]
    else:
        sam_sort_tasks = mapping_tasks

    total_pairs = 0
    pairs_txt_tasks = []
    # sort SAM files
    pairs_txt_files = []
    for i in range(len(file_names)):
        if file_types[i] != 'pairs_txt':
            continue
        pairs_txt_files.append(i)

    if len(pairs_txt_files) > 0:
        load_threads = max(int(threads / len(pairs_txt_files)), 1)

        pairs_files = []
        for ix in pairs_txt_files:
            pairs_txt_file = file_names[ix]
            pairs_file = os.path.join(output_folder, 'pairs', '{}_{}.pairs'.format(basename, total_pairs))
            total_pairs += 1
            if not force_overwrite and os.path.exists(pairs_file):
                parser.error("File exists ({}), use -f to force overwriting it.".format(pairs_file))

            pairs_files.append(pairs_file)

            pairs_command = fanc_base_command + ['pairs', '-f',
                                                 # loading
                                                 '-g', genome,
                                                 '-t', str(load_threads)]

            if restriction_enzyme is not None:
                pairs_command += ['-r', restriction_enzyme]
            if is_bwa:
                pairs_command.append('--bwa')
            if tmp:
                pairs_command.append('-tmp')
            if sam_sort:
                pairs_command.append('-S')

            pairs_command += [pairs_txt_file, pairs_file]

            pairs_task = CommandTask(pairs_command)
            runner.add_task(pairs_task, wait_for=[], threads=load_threads)
            pairs_txt_tasks.append(pairs_task)

            pairs_files.append(pairs_file)

        for ix, i in enumerate(pairs_txt_files):
            file_names[i] = pairs_files[ix]
            file_types[i] = 'pairs'

    # load pairs directly from SAM
    sam_file_pairs = []
    i = 0
    while i < len(file_names):
        if file_types[i] == 'sam':
            if not file_types[i + 1] == 'sam':
                parser.error("Cannot create SAM pairs, because {} "
                             "is missing a partner file".format(file_names[i]))
            sam_file_pairs.append((i, i + 1))
            i += 1
        i += 1

    if len(sam_file_pairs) > 0:
        sam_to_pairs_tasks = pairs_txt_tasks
        load_threads = max(int(threads/len(sam_file_pairs)), 1)

        pairs_files = []
        for i, j in sam_file_pairs:
            if len(sam_file_pairs) > 1 or total_pairs > 0:
                pairs_file = os.path.join(output_folder, 'pairs', '{}_{}.pairs'.format(basename, total_pairs))
                total_pairs += 1
            else:
                pairs_file = output_folder + 'pairs/' + basename + '.pairs'

            if not force_overwrite and os.path.exists(pairs_file):
                parser.error("File exists ({}), use -f to force overwriting it.".format(pairs_file))

            pairs_command = fanc_base_command + ['pairs', '-f',
                                                 # loading
                                                 '-g', genome,
                                                 '-t', str(load_threads),
                                                 # filtering
                                                 '-us']
            if restriction_enzyme is not None:
                pairs_command += ['-r', restriction_enzyme]
            if quality_cutoff is not None:
                pairs_command += ['-q', str(quality_cutoff)]
            if is_bwa:
                pairs_command.append('--bwa')
            if tmp:
                pairs_command.append('-tmp')
            if sam_sort:
                pairs_command.append('-S')

            pairs_command += [file_names[i], file_names[j], pairs_file]

            pairs_task = CommandTask(pairs_command)
            runner.add_task(pairs_task, wait_for=sam_sort_tasks, threads=load_threads)
            sam_to_pairs_tasks.append(pairs_task)

            pairs_files.append(pairs_file)

        for ix, sam_pair in enumerate(reversed(sam_file_pairs)):
            file_names[sam_pair[0]] = pairs_files[ix]
            del file_names[sam_pair[1]]
            file_types[sam_pair[0]] = 'pairs'
            del file_types[sam_pair[1]]
    else:
        sam_to_pairs_tasks = pairs_txt_tasks + sam_sort_tasks

    # 7. Pairs stats and filtering
    pairs_files = []
    for i in range(len(file_names)):
        if file_types[i] != 'pairs':
            continue
        pairs_files.append(i)

    if len(pairs_files) > 0 and filter_pairs:
        pairs_tasks = []
        for ix in pairs_files:
            pair_basename = os.path.basename(os.path.splitext(file_names[ix])[0])

            pairs_stats_file = output_folder + 'plots/stats/' + pair_basename + '.pairs.stats.pdf'
            ligation_error_file = output_folder + 'plots/stats/' + pair_basename + '.pairs.ligation_error.pdf'
            re_dist_file = output_folder + 'plots/stats/' + pair_basename + '.pairs.re_dist.pdf'

            pairs_command = fanc_base_command + ['pairs',
                                                 # filtering
                                                 '-d', str(max_restriction_site_distance),
                                                 '-l', '-p', '2']

            if tmp:
                pairs_command.append('-tmp')

            if inward_cutoff is not None:
                pairs_command += ['-i', str(inward_cutoff)]
                if outward_cutoff is None:
                    pairs_command += ['-o', '0']
            if outward_cutoff is not None:
                pairs_command += ['-o', str(outward_cutoff)]
                if inward_cutoff is None:
                    pairs_command += ['-i', '0']
            if inward_cutoff is None and outward_cutoff is None and auto_le_cutoff:
                pairs_command += ['--filter-ligation-auto']

            pairs_command += ['--statistics-plot', pairs_stats_file, file_names[ix]]
            pairs_task = CommandTask(pairs_command)
            runner.add_task(pairs_task, wait_for=sam_to_pairs_tasks, threads=1)
            pairs_tasks.append(pairs_task)

            ligation_error_command = fanc_base_command + ['pairs', '--ligation-error-plot',
                                                          ligation_error_file, file_names[ix]]
            ligation_error_task = CommandTask(ligation_error_command)
            runner.add_task(ligation_error_task, wait_for=pairs_task, threads=1)

            re_dist_command = fanc_base_command + ['pairs', '--re-dist-plot',
                                                   re_dist_file, file_names[ix]]
            re_dist_task = CommandTask(re_dist_command)
            runner.add_task(re_dist_task, wait_for=pairs_task, threads=1)
    else:
        pairs_tasks = sam_to_pairs_tasks

    # 8. Pairs to Hic
    pairs_files = []
    for i in range(len(file_names)):
        if file_types[i] != 'pairs':
            continue
        pairs_files.append(i)

    if len(pairs_files) > 0 and process_hic:
        hic_tasks = []
        hic_files = []
        for ix in pairs_files:
            hic_basename = os.path.basename(os.path.splitext(file_names[ix])[0])
            if hic_basename.endswith('_filtered'):
                hic_basename = hic_basename[:-9]
            hic_file = output_folder + 'hic/' + hic_basename + '.hic'

            if not force_overwrite and os.path.exists(hic_file):
                parser.error("File exists ({}), use -f to force overwriting it.".format(hic_file))

            hic_command = fanc_base_command + ['hic', '-f']
            if tmp:
                hic_command.append('-tmp')
            hic_command += [file_names[ix], hic_file]
            hic_task = CommandTask(hic_command)
            runner.add_task(hic_task, wait_for=pairs_tasks, threads=1)
            hic_tasks.append(hic_task)
            hic_files.append(hic_file)

        for ix, i in enumerate(pairs_files):
            file_names[i] = hic_files[ix]
            file_types[i] = 'hic'
    else:
        hic_tasks = pairs_tasks

    # 9. Merge Hic
    hic_files = []
    for i in range(len(file_names)):
        if file_types[i] != 'hic':
            continue
        hic_files.append(i)

    if len(hic_files) > 1:
        merge_hic_tasks = []
        output_hic = output_folder + 'hic/' + basename + '.hic'

        if not force_overwrite and os.path.exists(output_hic):
            parser.error("File exists ({}), use -f to force overwriting it.".format(output_hic))

        merge_hic_command = fanc_base_command + ['hic', '-f']
        if tmp:
            merge_hic_command.append('-tmp')

        hics = [file_names[i] for i in hic_files]
        merge_hic_command += hics + [output_hic]
        merge_hic_task = CommandTask(merge_hic_command)
        runner.add_task(merge_hic_task, wait_for=hic_tasks, threads=1)
        merge_hic_tasks.append(merge_hic_task)

        file_names[hic_files[0]] = output_hic
        hic_files.pop(0)
        for ix, i in enumerate(reversed(hic_files)):
            del file_names[i]
            del file_types[i]
    else:
        merge_hic_tasks = hic_tasks

    from fanc.tools.general import human_format, str_to_int

    hic_files = []
    for i in range(len(file_names)):
        if file_types[i] != 'hic':
            continue
        hic_files.append(i)

    if len(hic_files) > 0:
        for ix in hic_files:
            hic_file = file_names[ix]
            binned_hic_file_base = output_folder + 'hic/binned/' + basename + '_'
            bin_threads = min(4, threads)

            for bin_size in bin_sizes:
                bin_size = str_to_int(str(bin_size))
                bin_size_str = human_format(bin_size, 0).lower() + 'b'

                binned_hic_file = binned_hic_file_base + bin_size_str + '.hic'

                if not force_overwrite and os.path.exists(binned_hic_file):
                    parser.error("File exists ({}), use -f to force overwriting it.".format(binned_hic_file))

                hic_basename = os.path.basename(os.path.splitext(binned_hic_file)[0])
                hic_stats_file = output_folder + 'plots/stats/' + \
                                 hic_basename + '.stats.pdf'

                hic_command = fanc_base_command + ['hic', '-f', '-b', str(bin_size),
                                                   '-r', '0.1', '-t', str(bin_threads),
                                                   '--statistics-plot', hic_stats_file,
                                                   '-n', '--norm-method', norm_method]
                if tmp:
                    hic_command.append('-tmp')
                if restore_coverage:
                    hic_command.append('-c')

                hic_command += [hic_file, binned_hic_file]

                hic_task = CommandTask(hic_command)
                runner.add_task(hic_task, wait_for=merge_hic_tasks, threads=bin_threads)

    runner.run()
    return 0
