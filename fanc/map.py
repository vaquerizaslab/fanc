"""
This module performs iterative mapping of reads in FASTQ files to a reference genome.

The main function is :func:`~iterative_mapping`, which requires an input FASTQ file,
an output SAM file, and a suitable mapper. Mapper options are :class:`~Bowtie2Mapper`
and :class:`~BwaMapper` for the time being, but by subclassing :class:`~Mapper` a user
can easily write their own mapper implementations that can fully leverage the iterative
mapping capabilities of FAN-C. Take a look at the code of :class:`~Bowtie2Mapper` for
an example.

Example usage:

.. code::

    import fanc
    mapper = fanc.BwaMapper("bwa-index/hg19_chr18_19.fa", min_quality=3)
    fanc.iterative_mapping("SRR4271982_chr18_19_1.fastq.gzip", "SRR4271982_chr18_19_1.bam",
                           mapper, threads=4, restriction_enzyme="HindIII")

"""

import os
import io
import subprocess
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from queue import Empty, Full
import traceback
import gzip
import re
import threading
import pysam
import tempfile
import shutil
import uuid
from future.utils import string_types
from .config import config
from .tools.files import which
from .tools.general import ligation_site_pattern, split_at_ligation_junction, WorkerMonitor
import logging
logger = logging.getLogger(__name__)

__all__ = ['Bowtie2Mapper', 'BwaMapper',
           'SimpleBowtie2Mapper', 'SimpleBwaMapper',
           'iterative_mapping']

_read_name_re = re.compile("^@(.+?)\s(.*)$")
_read_name_nospace_re = re.compile("^@(.+)$")


def read_name(line):
    """
    Extract read name from FASTQ line.

    :param line: str of the FASTQ line, starting with "@"
    :return: tuple of read name and other info. info can be None
    """
    matches = _read_name_re.match(line)
    if matches is None:
        matches = _read_name_nospace_re.match(line)
        name = matches.group(1)
        info = None
    else:
        name = matches.group(1)
        info = matches.group(2)
    return name, info


class Monitor(WorkerMonitor):
    """
    Monitor class keeping an eye on mapping workers.

    Thread-safe by using locks on each attribute.

    .. attribute:: value

        Counter used to track worker progress.
    """
    def __init__(self, value=0, manager=None):
        if manager is None:
            manager = mp.Manager()
        WorkerMonitor.__init__(self, value=value, manager=manager)

        self.resubmitting_lock = manager.Lock()
        self.submitting_lock = manager.Lock()

        with self.resubmitting_lock:
            self.resubmitting = False

        with self.submitting_lock:
            self.submitting = False

    def set_resubmitting(self, value):
        """
        Set the resubmission status of a mapping worker.

        :param value: True if worker status changed to resubmitting,
                      False otherwise
        """
        with self.resubmitting_lock:
            self.resubmitting = value

    def is_resubmitting(self):
        """
        Query resubission status of a worker.

        :return: True if worker is in the process of resubmitting,
                 False otherwise
        """
        with self.resubmitting_lock:
            return self.resubmitting

    def set_submitting(self, value):
        """
        Set the submission status of a mapping worker.

        :param value: True if worker status changed to submitting,
                      False otherwise
        """
        with self.submitting_lock:
            self.submitting = value

    def is_submitting(self):
        """
        Query subission status of a worker.

        :return: True if worker is in the process of submitting,
                 False otherwise
        """
        with self.submitting_lock:
            return self.submitting


class Mapper(object):
    """
    Mapper base class for mapping sequencing reads to a reference genome.

    This base class handles the resubmission checks for iterative mapping,
    spawns mapping processes, and compiles a FASTQ that can be used as input
    for external mapping software.

    To actually do any work, this class must be inherited, and the _map and
    _resubmit methods must be overridden.

    :func:`~Mapper._map` must handle the actual mapping given an input file

    :func:`~Mapper._resubmit` determines whether a read should be resubmitted
    (after expanding it)

    See :class:`~Bowtie2Mapper` or :class:`~BwaMapper` for examples.

    .. attribute:: resubmit_unmappable

        Set to True (default) if unmappable reads should be resubmitted

    .. attribute:: attempt_resubmit

        Try to resubmit a read if it failed the :func:`~Mapper._resubmit`
        test.
    """

    def __init__(self, resubmit_unmappable=True, attempt_resubmit=True):
        self.resubmit_unmappable = resubmit_unmappable
        self.attempt_resubmit = attempt_resubmit

    def _map(self, input_file, output_file, *args, **kwargs):
        """
        Map reads in the given FASTQ file to a reference genome.

        Must be overridden by classes implementing :class:`~Mapper`.
        Arguments are provided at runtime by the iterative mapping
        process.

        :param input_file: Path to FASTQ file
        :param output_file: Path to SAM file
        :param args: Optional positional arguments
        :param kwargs: Optional keyword arguments
        :return: Must return 0 if mapping successful,
                 any other integer otherwise
        """
        raise NotImplementedError("Must implement _map method!")

    def map(self, input_file, output_folder=None):
        """
        Map reads in the given FASTQ file using :func:`~Mapper._map` implementation.

        Will internally map the FASTQ reads to a SAM file in a
        temporary folder (use output_folder to choose a specific folder),
        split the SAM output into (i) valid alignments according to
        :func:`~Mapper._resubmit` and (ii) invalid alignments that get
        caught by the resubmission filter. A FASTQ file will be constructed
        from the invalid alignments, extending the reads by a given step size,
        which can then be used to iteratively repeat the mapping process until
        a valid alignment is found or the full length of the read has been
        restored.

        :param input_file: Path to FASTQ file
        :param output_folder: (optional) path to temporary folder for SAM output
        :return: tuple, path to valid SAM alignments, path to resubmission FASTQ
        """
        if output_folder is None:
            output_folder = tempfile.mkdtemp()

        logger.debug('Output folder for SAM process: {}'.format(output_folder))

        with tempfile.NamedTemporaryFile(prefix='output', suffix='.sam', dir=output_folder,
                                         delete=False) as tmp:
            sam_output_file = tmp.name

            ret = self._map(input_file, sam_output_file)

            if isinstance(ret, tuple) and len(ret) == 2:
                ret, error_message = ret
            else:
                error_message = ''

        if ret != 0:
            raise RuntimeError('Mapping had non-zero exit status {}. {}'.format(ret, error_message))

        logger.debug('Done mapping')

        if not self.resubmit_unmappable and self.attempt_resubmit:
            sam_valid_file, resubmission_file = sam_output_file, None
        else:
            sam_valid_file, resubmit = self._valid_and_resubmissions(sam_output_file, output_folder)
            os.remove(sam_output_file)
            resubmission_file = self._resubmission_fastq(input_file, resubmit, output_folder)

        logger.debug('Mapper done.')
        return sam_valid_file, resubmission_file

    def _resubmit(self, sam_fields):
        """
        Determine if an alignment should be resubmitted.

        :param sam_fields: The individual fields in a SAM line (split by tab)
        :return: True if read should be extended and aligned to the reference again,
                 False if the read passes the validity criteria.
        """
        raise NotImplementedError("Mapper must implement '_resubmit'")

    def resubmit(self, sam_fields):
        """
        Determine if an alignment should be resubmitted.

        Filters unmappable reads by default. Additional criteria can be
        implemented using the :func:`~Mapper._resubmit` method.

        :param sam_fields: The individual fields in a SAM line (split by tab)
        :return: True if read should be extended and aligned to the reference again,
                 False if the read passes the validity criteria.
        """
        if self.resubmit_unmappable and int(sam_fields[1]) & 4:
            return True
        return self._resubmit(sam_fields)

    def _valid_and_resubmissions(self, input_sam_file, output_folder):
        """
        Determine valid alignments and resubmissions from SAM file.

        :param input_sam_file: Path to SAM file
        :param output_folder: Path to temporary output folder
        :return: tuple with path to valid SAM file, dictionary with
                 read names that should be resubmitted
        """
        logger.debug('Getting mapped reads and resubmissions')
        resubmit = dict()
        with io.open(input_sam_file) as f:
            with tempfile.NamedTemporaryFile(mode='w', prefix='valid', suffix='.sam',
                                             delete=False, dir=output_folder) as o:
                output_sam_file = o.name
                for i, line in enumerate(f):
                    line = line.rstrip()
                    if line == '':
                        continue
                    if line.startswith('@'):
                        o.write(line + '\n')
                        continue

                    sam_fields = line.split("\t")

                    if self._resubmit(sam_fields):
                        if not sam_fields[0] in resubmit:
                            resubmit[sam_fields[0]] = True
                    else:
                        resubmit[sam_fields[0]] = False
                        o.write(line + '\n')

        return output_sam_file, resubmit

    def _resubmission_fastq(self, input_fastq, resubmit, output_folder):
        """
        Create the FASTQ file for resubmissions.

        :param input_fastq: Path to original FASTQ file
        :param resubmit: Dictionary with read names that should be expanded
                         and then resubmitted
        :param output_folder: temporary folder
        :return: Path to resubmission FASTQ
        """
        logger.debug('Creating FASTQ file for resubmission')

        resubmission_counter = 0
        total_counter = 0
        with io.open(input_fastq) as f:
            with tempfile.NamedTemporaryFile(mode='w', prefix='resubmission_', suffix='.fastq',
                                             delete=False, dir=output_folder) as o:
                output_fastq = o.name
                resubmit_current = False
                for i, line in enumerate(f):
                    line = line.rstrip()
                    if line == '':
                        continue

                    if i % 4 == 0:
                        name, _ = read_name(line)
                        if name not in resubmit:
                            if self.resubmit_unmappable:
                                resubmit_current = True
                        elif resubmit[name]:
                            resubmit_current = True
                        else:
                            resubmit_current = False

                    if resubmit_current:
                        resubmission_counter += 1
                        o.write(line + '\n')
                    total_counter += 1

        if resubmission_counter == 0:
            os.remove(output_fastq)
            return None

        logger.debug("Resubmitting {}/{}".format(resubmission_counter/4, total_counter/4))
        return output_fastq

    def close(self):
        """
        Final operations after mapping completes.
        """
        pass


class Bowtie2Mapper(Mapper):
    """
    Bowtie2 Mapper for aligning reads against a reference genome.

    Implements :class:`~Mapper` by calling the command line "bowtie2"
    program.

    .. attribute:: bowtie2_index

        Path to the bowtie2 index of the reference genome of choice.

    .. attribute:: min_quality

        Minimum MAPQ of an alignment so that it won't be resubmitted in
        iterative mapping.

    .. attribute:: additional_arguments

        Arguments passed to the "bowtie2" command in addition to
        -x, -U, --no-unal, --threads, and -S.

    .. attribute:: threads

        Number of threads for this mapping process.

    """
    def __init__(self, bowtie2_index, min_quality=30, additional_arguments=(),
                 threads=1, _bowtie2_path=config.bowtie2_path, **kwargs):
        Mapper.__init__(self, **kwargs)
        self.index = os.path.expanduser(bowtie2_index)
        if self.index.endswith('.'):
            self.index = self.index[:-1]
        self.args = [a for a in additional_arguments]
        self._path = _bowtie2_path
        if self._path is None:
            self._path = "bowtie2"
        if which(self._path) is None:
            raise ValueError("Cannot find {}".format(self._path))
        self.min_quality = min_quality
        self.threads = threads
        self.attempt_resubmit = (self.min_quality is not None and self.min_quality > 0)

    def _map(self, input_file, output_file, *args, **kwargs):
        bowtie2_command = [self._path, '-x', self.index, '-U', input_file, '--no-unal',
                           '--threads', str(self.threads), '-S', output_file] + self.args
        logger.debug('Bowtie2 command: {}'.format(bowtie2_command))

        proc = subprocess.Popen(bowtie2_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, universal_newlines=True)
        proc.wait()

        if proc.returncode != 0:
            error_message = proc.stderr.read()
        else:
            error_message = ''

        return proc.returncode, error_message

    def _resubmit(self, sam_fields):
        if int(sam_fields[4]) < self.min_quality:
            return True
        return False


class SimpleBowtie2Mapper(Bowtie2Mapper):
    """
    Bowtie2 Mapper for aligning reads against a reference genome without resubmission.

    Implements :class:`~Mapper` by calling the command line "bowtie2"
    program. Does not resubmit reads under any circumstance.

    .. attribute:: bowtie2_index

        Path to the bowtie2 index of the reference genome of choice.

    .. attribute:: additional_arguments

        Arguments passed to the "bowtie2" command in addition to
        -x, -U, --no-unal, --threads, and -S.

    .. attribute:: threads

        Number of threads for this mapping process.
    """
    def __init__(self, bowtie2_index, additional_arguments=(),
                 threads=1, _bowtie2_path='bowtie2'):
        Bowtie2Mapper.__init__(self, bowtie2_index, min_quality=0,
                               additional_arguments=additional_arguments,
                               threads=threads,
                               _bowtie2_path=_bowtie2_path,
                               resubmit_unmappable=False)

    def _resubmit(self, sam_fields):
        return False


class BwaMapper(Mapper):
    """
    BWA Mapper for aligning reads against a reference genome.

    Implements :class:`~Mapper` by calling the command line "bwa"
    program.

    .. attribute:: bwa_index

        Path to the BWA index of the reference genome of choice.

    .. attribute:: min_quality

        Minimum MAPQ of an alignment so that it won't be resubmitted in
        iterative mapping.

    .. attribute:: additional_arguments

        Arguments passed to the "bowtie2" command in addition to
        -t and -o.

    .. attribute:: threads

        Number of threads for this mapping process.

    .. attribute:: algorithm

        BWA algorithm to use for mapping. Uses "mem" by default.
        See http://bio-bwa.sourceforge.net/bwa.shtml for other options.

    .. attribute: memory_map

        Set to true if you want to commit the entire BWA index to memory.
        This sometimes leads to clashes between different BWA threads.

    """
    def __init__(self, bwa_index, min_quality=0, additional_arguments=(),
                 threads=1, algorithm='mem', memory_map=False, _bwa_path=config.bwa_path):
        Mapper.__init__(self)
        self.index = os.path.expanduser(bwa_index)
        if self.index.endswith('.'):
            self.index = self.index[:-1]
        self.args = [a for a in additional_arguments]
        self._path = _bwa_path
        if self._path is None:
            self._path = 'bwa'
        if which(self._path) is None:
            raise ValueError("Cannot find {}".format(self._path))
        self.algorithm = algorithm
        self.min_quality = min_quality
        self.threads = threads
        self.attempt_resubmit = (self.min_quality is not None and self.min_quality > 0)
        self.memory_map = memory_map

    def _map(self, input_file, output_file, *args, **kwargs):
        if self.memory_map:
            with open(os.devnull, 'w') as f_null:
                ret = subprocess.call(['bwa', 'shm', self.index], stderr=subprocess.STDOUT, stdout=f_null)
                if ret == 0:
                    logger.debug("Memory mapped BWA index")

        bwa_command = [self._path, self.algorithm, '-t', str(self.threads), '-o', output_file] + \
                      self.args + \
                      [self.index, input_file]
        logger.debug('BWA command: {}'.format(bwa_command))

        proc = subprocess.Popen(bwa_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, universal_newlines=True)
        proc.wait()

        if proc.returncode != 0:
            error_message = proc.stderr.read()
        else:
            error_message = ''

        return proc.returncode, error_message

    def _resubmit(self, sam_fields):
        if int(sam_fields[4]) < self.min_quality:
            return True
        return False

    def close(self):
        if self.memory_map:
            logger.debug("Deleting shared memory index")
            try:
                with open(os.devnull, 'w') as f_null:
                    ret = subprocess.call(['pgrep', '-x', 'bwa'], stderr=subprocess.STDOUT, stdout=f_null)
                    if ret == 1:
                        ret2 = subprocess.call(['bwa', 'shm', '-d'], stderr=subprocess.STDOUT, stdout=f_null)
                        if ret2 == 0:
                            logger.debug("Deleted memory mapped BWA index")

            except FileNotFoundError:
                pass


class SimpleBwaMapper(BwaMapper):
    """
    BWA Mapper for aligning reads against a reference genome without resubmission.

    Implements :class:`~Mapper` by calling the command line "bwa"
    program. Does not resubmit reads under any circumstance, i.e.
    does not perform iterative mapping.

    .. attribute:: bwa_index

        Path to the BWA index of the reference genome of choice.

    .. attribute:: min_quality

        Minimum MAPQ of an alignment so that it won't be resubmitted in
        iterative mapping.

    .. attribute:: additional_arguments

        Arguments passed to the "bowtie2" command in addition to
        -t and -o.

    .. attribute:: threads

        Number of threads for this mapping process.

    .. attribute:: algorithm

        BWA algorithm to use for mapping. Uses "mem" by default.
        See http://bio-bwa.sourceforge.net/bwa.shtml for other options.

    .. attribute: memory_map

        Set to true if you want to commit the entire BWA index to memory.
        This sometimes leads to clashes between different BWA threads.

    """
    def __init__(self, bwa_index, additional_arguments=(),
                 threads=1, memory_map=False, _bwa_path='bwa'):
        BwaMapper.__init__(self, bwa_index, min_quality=0,
                           additional_arguments=additional_arguments,
                           threads=threads, memory_map=memory_map,
                           _bwa_path=_bwa_path)
        self.resubmit_unmappable = False
        self.attempt_resubmit = False

    def _resubmit(self, sam_fields):
        return False


def _trim_read(input_seq, step_size=5, min_size=25, front=False):
    """
    Trim an input sequence by step_size.

    :param input_seq: input DNA sequence
    :param step_size: truncation step size
    :param min_size: minimum size of input DNA before error is raised
    :param front: If True, trim from front instead of back.
    :return: tuple: read name, trimmed sequence, '+', and trimmed quality
                    string
    """
    name, seq, plus, qual = input_seq
    if len(seq) <= min_size:
        raise ValueError("Already reached minimum size, cannot truncate read further")

    if len(seq) - step_size < min_size:
        final_length = min_size
    else:
        final_length = len(seq) - step_size

    if front:
        sl = len(seq)
        seq = seq[sl - final_length:]
        qual = qual[sl - final_length:]
    else:
        seq = seq[:final_length]
        qual = qual[:final_length]
    return name, seq, plus, qual


def _iterative_mapping_worker(mapper, input_queue, output_folder, output_queue,
                              resubmission_queue, monitor, exception_queue):
    """
    Worker performing the iterative mapping of FASTQ files.

    :param mapper: :class:`~Mapper` implementation
    :param input_queue: queue with input FASTQ files
    :param output_folder: Output folder for intermediate results
    :param output_queue: queue for valid output alignments
    :param resubmission_queue: queue for resubmissions (iterative mapping)
    :param monitor: worker monitor
    :param exception_queue: queue for exceptions
    """
    try:
        # generate worker's uuid
        worker_uuid = uuid.uuid4()

        while True:
            # set worker state to idle
            monitor.set_worker_idle(worker_uuid)
            logger.debug('Worker {} idle'.format(worker_uuid))

            # wait for input
            input_file = input_queue.get(True)
            monitor.set_worker_busy(worker_uuid)
            logger.debug('Mapper {} busy, got input file'.format(worker_uuid))

            # mapping file
            sam_file, unmapped_file = mapper.map(input_file, output_folder=output_folder)
            logger.debug('{} done mapping'.format(worker_uuid))

            # clean up
            os.remove(input_file)
            logger.debug('{} waiting to put SAM in output queue'.format(worker_uuid))
            output_queue.put(sam_file)

            # send resubmissions back to writing thread
            logger.debug('{} waiting to put FASTQ in resubmission queue'.format(worker_uuid))
            resubmission_queue.put(unmapped_file)
            logger.debug('{} finished mapping round.'.format(worker_uuid))
    except Exception:
        import sys
        exception_queue.put("".join(traceback.format_exception(*sys.exc_info())))


def _fastq_to_queue(fastq_file, output_folder, batch_size, input_queue, monitor,
                    exception_queue=None, worker_pool=None, restriction_enzyme=None):
    """
    Submit FASTQ batches to input queue.

    Also does ligation site splitting, if requested.

    :param fastq_file: Path to input FASTQ file
    :param output_folder: Path to output folder (tmp)
    :param batch_size: Number of reads submitted to a monitor
    :param input_queue: queue for batch FASTQ files
    :param monitor: worker monitor
    :param exception_queue: queue for exceptions
    :param worker_pool: Pool of mapping workers "_iterative_mapping_worker"
    :param restriction_enzyme: Name of restriction enzyme of restriction pattern
                               (e.g. A^AGCT_T)
    """
    monitor.set_submitting(True)

    if restriction_enzyme is not None:
        if isinstance(restriction_enzyme, string_types):
            restriction_enzyme = [restriction_enzyme]

        ligation_patterns = [ligation_site_pattern(r) for r in restriction_enzyme]
    else:
        ligation_patterns = None

    tmp_output_file = None
    line_counter = 0
    submission_counter = 0
    try:
        if fastq_file.endswith('.gz') or fastq_file.endswith('.gzip'):
            open_file = lambda x: io.BufferedReader(gzip.open(x, 'r'),
                                                    buffer_size=4*io.DEFAULT_BUFFER_SIZE)
        else:
            open_file = io.open

        read_counter = 0
        with open_file(fastq_file) as f:
            tmp_output_file = tempfile.NamedTemporaryFile(suffix='.fastq', delete=False,
                                                          dir=output_folder, mode='w+b')

            current_fastq = []
            for i, line in enumerate(f):
                line = line.encode('utf-8') if not isinstance(line, bytes) else line
                current_fastq.append(line)
                line_counter += 1

                # line = line.decode() if isinstance(line, bytes) else line
                if line_counter % 4 == 0:
                    if not ligation_patterns:
                        for fastq_line in current_fastq:
                            tmp_output_file.write(fastq_line)
                        read_counter += 1
                    else:
                        name, info = read_name(current_fastq[0].rstrip().decode())
                        seq = current_fastq[1].rstrip()
                        qual = current_fastq[3].rstrip()
                        seqs = split_at_ligation_junction(seq.decode(), ligation_patterns)

                        qualities = []
                        current_pos = 0
                        for s in seqs:
                            qualities.append(qual[current_pos:current_pos+len(s)])
                            current_pos += len(s)

                        for j in range(len(seqs)):
                            if len(seqs) > 1:
                                new_name = '@' + name + '__{}'.format(j)
                            else:
                                new_name = '@' + name

                            if info is not None:
                                new_name += ' ' + info
                            tmp_output_file.write(new_name.encode('utf-8') + b'\n')
                            tmp_output_file.write(seqs[j].encode('utf-8') + b'\n')
                            tmp_output_file.write(b'+\n')
                            tmp_output_file.write(qualities[j] + b'\n')
                            read_counter += 1
                    current_fastq = []

                    if line_counter/4 >= batch_size:
                        tmp_output_file.close()
                        line_counter = 0

                        if exception_queue is not None:
                            while True:
                                try:
                                    input_queue.put(tmp_output_file.name, True, 5)
                                    monitor.increment()
                                    submission_counter += 1
                                    break
                                except Full:
                                    pass
                                try:
                                    exc = exception_queue.get(block=False)
                                except Empty:
                                    pass
                                else:
                                    worker_pool.terminate()
                                    raise Exception(exc)
                        else:
                            input_queue.put(tmp_output_file.name, True, 5)
                            monitor.increment()
                            submission_counter += 1
                        tmp_output_file = tempfile.NamedTemporaryFile(suffix='.fastq', delete=False,
                                                                      dir=output_folder, mode='w+b')
        tmp_output_file.close()
        if line_counter > 0:
            input_queue.put(tmp_output_file.name, True)
            monitor.increment()
            submission_counter += 1
    except Exception:
        import sys
        stacktrace = "".join(traceback.format_exception(*sys.exc_info()))
        logger.error(stacktrace)
        exception_queue.put(stacktrace)
    finally:
        if tmp_output_file is not None and not tmp_output_file.closed:
            tmp_output_file.close()
    logger.debug("Submitted {} FASTQ chunks".format(submission_counter))
    monitor.set_submitting(False)


def _resubmissions_to_queue(resubmission_queue, output_folder, batch_size,
                            input_queue, monitor, step_size=5, min_size=25,
                            trim_front=False,
                            exception_queue=None, worker_pool=None):
    """
    Collect resubmission files, process the FASTQ, and resubmit to input queue.

    :param resubmission_queue: queue for resubmission files
    :param output_folder: output folder (tmp)
    :param batch_size: Number of reads submitted to a monitor
    :param input_queue: queue for input (batch) FASTQ files
    :param monitor: worker monitor
    :param step_size: Step size for read trimming
    :param min_size: Minimum read size
    :param trim_front: If True, trims reads from front instead of back
    :param exception_queue: queue for exceptions
    :param worker_pool: Pool of mapping workers "_iterative_mapping_worker"
    :return:
    """
    collected_resubmits = 0
    tmp_output_file = None
    read_counter = 0
    try:
        tmp_output_file = tempfile.NamedTemporaryFile(suffix='.fastq', delete=False,
                                                      dir=output_folder, mode='w+b')

        while monitor.is_submitting() or collected_resubmits < monitor.value():
            logger.debug("Status: {}/{}".format(collected_resubmits, monitor.value()))
            monitor.set_resubmitting(False)
            resubmission_file = resubmission_queue.get(block=True)
            monitor.set_resubmitting(True)
            logger.debug('Got resubmission file {}'.format(resubmission_file))

            collected_resubmits += 1
            if resubmission_file is None:
                continue

            with io.open(resubmission_file) as f:

                current_fastq = []
                for i, line in enumerate(f):
                    current_fastq.append(line.rstrip())

                    if i % 4 == 3:  # FASTQ record is complete
                        try:
                            new_fastq = _trim_read(current_fastq, step_size=step_size,
                                                   min_size=min_size, front=trim_front)
                        except ValueError:
                            continue

                        if len(new_fastq[1]) != len(current_fastq[1]):
                            for fastq_line in new_fastq:
                                new_line = fastq_line + '\n'
                                new_line = new_line.encode('utf-8')
                                tmp_output_file.write(new_line)
                            read_counter += 1
                        current_fastq = []

                        if read_counter >= batch_size:
                            logger.debug("Resubmitting because batch is full ({})".format(read_counter))
                            tmp_output_file.close()
                            read_counter = 0

                            if exception_queue is not None:
                                while True:
                                    try:
                                        input_queue.put(tmp_output_file.name, True, 5)
                                        monitor.increment()
                                        break
                                    except Full:
                                        pass
                                    try:
                                        exc = exception_queue.get(block=False)
                                    except Empty:
                                        pass
                                    else:
                                        worker_pool.terminate()
                                        raise Exception(exc)
                            else:
                                input_queue.put(tmp_output_file.name, True)
                                monitor.increment()
                            tmp_output_file = tempfile.NamedTemporaryFile(suffix='.fastq', delete=False,
                                                                          dir=output_folder, mode='w+b')
            os.remove(resubmission_file)

            if (read_counter > 0 and monitor.workers_idle()
                    and not monitor.is_submitting() and resubmission_queue.empty()):
                logger.debug("Resubmitting prematurely ({}) because workers are waiting for input".format(read_counter))
                tmp_output_file.close()
                input_queue.put(tmp_output_file.name, True)
                monitor.increment()
                tmp_output_file = tempfile.NamedTemporaryFile(suffix='.fastq', delete=False,
                                                              dir=output_folder, mode='w+b')
    except Exception:
        import sys
        stacktrace = "".join(traceback.format_exception(*sys.exc_info()))
        logger.error(stacktrace)
        exception_queue.put(stacktrace)
    finally:
        monitor.set_resubmitting(False)
        if tmp_output_file is not None and not tmp_output_file.closed:
            tmp_output_file.close()
        os.remove(tmp_output_file.name)


def iterative_mapping(fastq_file, sam_file, mapper, tmp_folder=None, threads=1, min_size=25, step_size=5,
                      batch_size=200000, trim_front=False, restriction_enzyme=None):
    """
    Iteratively map sequencing reads using the provided mapper.

    Will attempt to align a read using mapper. If unsuccessful, will
    truncate the read by step_size and attempt to align again. This is
    repeated until a successful alignment is found or the read gets
    truncated below min_size.

    :param fastq_file: An input FASTQ file path with reds to align
    :param sam_file: An output file path for sequencing results.
                     If it ends with '.bam' will compress output
                     in bam format.
    :param mapper: An instance of :class:`Mapper`, e.g. :class:`Bowtie2Mapper`.
                   Override :class:`Mapper` for creating your own custom mappers.
    :param tmp_folder: A temporary folder for outputting subsets of FASTQ files
    :param threads: Number of mapper threads to use in parallel.
    :param min_size: Minimum length of read for which an alignment is attempted.
    :param step_size: Number of base pairs by which to truncate read.
    :param batch_size: Maximum number of reads processed in one batch
    :param trim_front: Trim bases from front of read instead of back
    :param restriction_enzyme: If provided, will calculate the expected ligation
                               junction between reads and split reads accordingly.
                               Both ends will be attempted to map. Can be the name
                               of a restriction enzyme or a restriction pattern
                               (e.g. A^AGCT_T)
    """
    if tmp_folder is None:
        tmp_folder = tempfile.mkdtemp()
    else:
        tmp_folder = tempfile.mkdtemp(dir=tmp_folder)
    logger.debug('Tmp folder: {}'.format(tmp_folder))

    try:
        mp_manager = mp.Manager()
        input_queue = mp_manager.Queue(maxsize=threads)
        resubmission_queue = mp_manager.Queue()
        output_queue = mp_manager.Queue()
        exception_queue = mp_manager.Queue()
        monitor = Monitor()

        monitor.set_submitting(True)

        with mp.get_context("spawn").Pool(threads, _iterative_mapping_worker,
                              (mapper, input_queue, tmp_folder, output_queue,
                               resubmission_queue, monitor, exception_queue)) as worker_pool:

            t_resub = threading.Thread(target=_resubmissions_to_queue, args=(resubmission_queue, tmp_folder,
                                                                             batch_size, input_queue, monitor,
                                                                             step_size, min_size, trim_front,
                                                                             exception_queue, worker_pool))
            t_resub.daemon = True
            t_resub.start()

            t_sub = threading.Thread(target=_fastq_to_queue, args=(fastq_file, tmp_folder,
                                                                   batch_size, input_queue,
                                                                   monitor, exception_queue,
                                                                   worker_pool, restriction_enzyme))
            t_sub.daemon = True
            t_sub.start()

            # check BAM status
            if not sam_file.endswith('bam'):
                intermediate_sam_file = sam_file
                convert_to_bam = False
                logger.info("Starting to output alignments to SAM file {}".format(sam_file))
            else:
                intermediate_sam_file = os.path.join(tmp_folder, 'intermediate.sam')
                convert_to_bam = True
                logger.info("Starting to output alignments to intermediate SAM file {}".format(intermediate_sam_file))

            ligation_name_pattern = re.compile(r'(.+)__(\d)+$')
            sam_counter = 0
            with open(intermediate_sam_file, 'w') as o:
                while (sam_counter < monitor.value() or monitor.is_resubmitting()
                       or monitor.is_submitting() or not monitor.workers_idle()):
                    try:
                        exc = exception_queue.get(block=False)
                    except Empty:
                        pass
                    else:
                        worker_pool.terminate()
                        raise Exception(exc)
    
                    try:
                        partial_sam_file = output_queue.get(block=True, timeout=10)
                        logger.debug('Processing output file {}'.format(partial_sam_file))

                        with open(partial_sam_file, 'r') as f:
                            for line in f:
                                if line.startswith('@'):
                                    if sam_counter == 0:
                                        o.write(line)
                                    continue

                                line = line.rstrip()
                                if line == '':
                                    continue

                                fields = line.split("\t")
                                m = ligation_name_pattern.match(fields[0])
                                if m is not None:
                                    fields[0] = m.group(1)
                                    fields[-1] += '\tZL:i:{}'.format(m.group(2))

                                o.write("\t".join(fields) + '\n')

                        os.remove(partial_sam_file)

                        sam_counter += 1
                        logger.debug('Got {}/{} SAM files'.format(sam_counter, monitor.value()))
                    except Empty:
                        pass

            t_sub.join()
            t_resub.join()

        if convert_to_bam:
            logger.info("Converting intermediate SAM file to BAM ({})".format(sam_file))
            try:
                pysam.view('-b', '-o', sam_file, intermediate_sam_file, catch_stdout=False)
                logger.info("Success, removing intermediate.")
                os.remove(intermediate_sam_file)
            except pysam.SamtoolsError:
                logger.error("Could not convert to BAM, but your "
                             "output is still in {}".format(intermediate_sam_file))
    finally:
        logger.debug(tmp_folder)
        shutil.rmtree(tmp_folder, ignore_errors=True)

    return sam_file
