import os
import io
import subprocess
import multiprocessing as mp
from queue import Empty, Full
import traceback
import gzip
import re
import threading
import tempfile
import shutil
import uuid
from kaic.tools.general import which
import logging
logger = logging.getLogger(__name__)


class Monitor(object):
    def __init__(self, value=0):
        self.counter_lock = threading.Lock()
        self.resubmitting_lock = threading.Lock()
        self.submitting_lock = threading.Lock()
        self.worker_lock = threading.Lock()

        with self.counter_lock:
            self.val = value

        with self.resubmitting_lock:
            self.resubmitting = False

        with self.submitting_lock:
            self.submitting = False

        with self.worker_lock:
            self.worker_states = dict()

    def increment(self):
        with self.counter_lock:
            self.val += 1

    def value(self):
        with self.counter_lock:
            return self.val

    def set_resubmitting(self, value):
        with self.resubmitting_lock:
            self.resubmitting = value

    def is_resubmitting(self):
        with self.resubmitting_lock:
            return self.resubmitting

    def set_submitting(self, value):
        with self.submitting_lock:
            self.submitting = value

    def is_submitting(self):
        with self.submitting_lock:
            return self.submitting

    def set_worker_busy(self, worker_uuid):
        with self.worker_lock:
            self.worker_states[worker_uuid] = True

    def set_worker_idle(self, worker_uuid):
        with self.worker_lock:
            self.worker_states[worker_uuid] = False

    def get_worker_state(self, worker_uuid):
        with self.worker_lock:
            return self.worker_states[worker_uuid]

    def workers_idle(self):
        with self.worker_lock:
            for busy in self.worker_states.values():
                if busy:
                    return False
            return True


class Mapper(object):
    def __init__(self):
        self.resubmit_unmappable = True

    def map(self, input_file, output_folder=None):
        raise NotImplementedError("Mapper must implement 'map'")

    def _resubmit(self, sam_fields):
        raise NotImplementedError("Mapper must implement 'resubmit'")

    def resubmit(self, sam_fields):
        if self.resubmit_unmappable and int(sam_fields[1]) & 4:
            return True
        return self._resubmit(sam_fields)


class Bowtie2Mapper(Mapper):
    def __init__(self, bowtie2_index, min_quality=30, additional_arguments=(),
                 threads=1, _bowtie2_path='bowtie2'):
        Mapper.__init__(self)
        self.index = os.path.expanduser(bowtie2_index)
        if self.index.endswith('.'):
            self.index = self.index[:-1]
        self.args = [a for a in additional_arguments]
        self._path = _bowtie2_path
        if which(self._path) is None:
            raise ValueError("Cannot find {}".format(self._path))
        self.min_quality = min_quality
        self.threads = threads

    def map(self, input_file, output_folder=None):
        if output_folder is None:
            output_folder = tempfile.mkdtemp()

        logger.debug('Output folder for SAM process: {}'.format(output_folder))

        with tempfile.NamedTemporaryFile(prefix='output', suffix='.sam', dir=output_folder,
                                         delete=False) as tmp:
            sam_output_file = tmp.name

        bowtie2_command = [self._path, '-x', self.index, '-U', input_file, '--no-unal',
                           '--threads', str(self.threads), '-S', sam_output_file] + self.args
        logger.debug('Bowtie2 command: {}'.format(bowtie2_command))

        f_null = open(os.devnull, 'w')
        ret = subprocess.call(bowtie2_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=f_null, universal_newlines=True)
        f_null.close()

        if ret != 0:
            raise RuntimeError('Bowtie2 had non-zero exit status {}'.format(ret))

        logger.debug('Done mapping')

        if not self.resubmit_unmappable and (self.min_quality is None or self.min_quality == 0):
            sam_valid_file, resubmission_file = sam_output_file, None
        else:
            sam_valid_file, resubmit = self._valid_and_resubmissions(sam_output_file, output_folder)
            os.remove(sam_output_file)
            resubmission_file = self._resubmission_fastq(input_file, resubmit, output_folder)

        logger.debug('Mapper done.')
        return sam_valid_file, resubmission_file

    def _valid_and_resubmissions(self, input_sam_file, output_folder):
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
                        resubmit[sam_fields[0]] = True
                    else:
                        resubmit[sam_fields[0]] = False
                        o.write(line + '\n')

        return output_sam_file, resubmit

    def _resubmission_fastq(self, input_fastq, resubmit, output_folder):
        logger.debug('Creating FASTQ file for resubmission')

        name_re = re.compile("^@(.+?)\s.*$")

        resubmission_counter = 0
        total_counter = 0
        with io.open(input_fastq) as f:
            with tempfile.NamedTemporaryFile(mode='w', prefix='resubmission', suffix='.fastq',
                                             delete=False, dir=output_folder) as o:
                output_fastq = o.name
                resubmit_current = False
                for i, line in enumerate(f):
                    line = line.rstrip()
                    if line == '':
                        continue

                    if i % 4 == 0:
                        matches = name_re.match(line)
                        name = matches.group(1)
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

    def _resubmit(self, sam_fields):
        if int(sam_fields[4]) < self.min_quality:
            return True
        return False


class SimpleBowtie2Mapper(Bowtie2Mapper):
    def __init__(self, bowtie2_index, additional_arguments=(),
                 threads=1, _bowtie2_path='bowtie2'):
        Bowtie2Mapper.__init__(self, bowtie2_index, min_quality=0,
                               additional_arguments=additional_arguments,
                               threads=threads,
                               _bowtie2_path=_bowtie2_path)
        self.resubmit_unmappable = False

    def _resubmit(self, sam_fields):
        return False


def _trim_read(input_seq, step_size=5, min_size=25):
    name, seq, plus, qual = input_seq
    if len(seq) == min_size:
        raise ValueError("Already reached minimum size, cannot truncate read further")
    if len(seq) - step_size < min_size:
        final_length = min_size
    else:
        final_length = len(seq) - step_size
    seq = seq[:final_length]
    qual = qual[:final_length]
    return name, seq, plus, qual


def _iterative_mapping_worker(mapper, input_queue, output_folder, output_queue,
                              resubmission_queue, monitor, exception_queue):
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
            logger.debug('{} waiting to put SAM in output queue')
            output_queue.put(sam_file)

            # send resubmissions back to writing thread
            logger.debug('{} waiting to put FASTQ in resubmission queue')
            resubmission_queue.put(unmapped_file)
            logger.debug('{} finished mapping round.')
    except Exception:
        import sys
        exception_queue.put("".join(traceback.format_exception(*sys.exc_info())))


def _fastq_to_queue(fastq_file, output_folder, batch_size, input_queue, monitor,
                    exception_queue=None, worker_pool=None):
    monitor.set_submitting(True)
    tmp_output_file = None
    line_counter = 0
    submission_counter = 0
    try:
        if fastq_file.endswith('.gz') or fastq_file.endswith('.gzip'):
            open_file = lambda x: io.BufferedReader(gzip.open(x, 'r'), buffer_size=4*io.DEFAULT_BUFFER_SIZE)
        else:
            open_file = io.open

        read_counter = 0
        with open_file(fastq_file) as f:
            tmp_output_file = tempfile.NamedTemporaryFile(suffix='.fastq', delete=False,
                                                          dir=output_folder)

            for i, line in enumerate(f):
                tmp_output_file.write(line)
                line_counter += 1
                read_counter += 1

                # line = line.decode() if isinstance(line, bytes) else line
                if line_counter % 4 == 0 and line_counter/4 >= batch_size:
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
                                                                  dir=output_folder)
        tmp_output_file.close()
        if line_counter > 0:
            input_queue.put(tmp_output_file.name, True)
            monitor.increment()
            submission_counter += 1
    except Exception:
        import sys
        exception_queue.put("".join(traceback.format_exception(*sys.exc_info())))
    finally:
        if tmp_output_file is not None and not tmp_output_file.closed:
            tmp_output_file.close()
    logger.debug("Submitted {} FASTQ chunks".format(submission_counter))
    monitor.set_submitting(False)


def _resubmissions_to_queue(resubmission_queue, output_folder, batch_size,
                            input_queue, monitor, step_size=5, min_size=25,
                            exception_queue=None, worker_pool=None):
    collected_resubmits = 0
    tmp_output_file = None
    read_counter = 0
    try:
        tmp_output_file = tempfile.NamedTemporaryFile(suffix='.fastq', delete=False,
                                                      dir=output_folder, mode='w')

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
                                                   min_size=min_size)
                        except ValueError:
                            continue

                        if len(new_fastq[1]) != len(current_fastq[1]):
                            for fastq_line in new_fastq:
                                tmp_output_file.write(fastq_line + '\n')
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
                                                                          dir=output_folder, mode='w')
            os.remove(resubmission_file)

            if read_counter > 0 and monitor.workers_idle() and not monitor.is_submitting():
                logger.debug("Resubmitting prematurely ({}) because workers are waiting for input".format(read_counter))
                tmp_output_file.close()
                input_queue.put(tmp_output_file.name, True)
                monitor.increment()
                tmp_output_file = tempfile.NamedTemporaryFile(suffix='.fastq', delete=False,
                                                              dir=output_folder, mode='w')
    except Exception:
        import sys
        exception_queue.put("".join(traceback.format_exception(*sys.exc_info())))
    finally:
        monitor.set_resubmitting(False)
        if tmp_output_file is not None and not tmp_output_file.closed:
            tmp_output_file.close()
        os.remove(tmp_output_file.name)


def iterative_mapping(fastq_file, sam_file, mapper, tmp_folder=None, threads=1, min_size=25, step_size=5,
                      batch_size=200000):
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
    :return:
    """
    if tmp_folder is None:
        tmp_folder = tempfile.mkdtemp()
    else:
        tmp_folder = tempfile.mkdtemp(dir=tmp_folder)
    logger.debug('Tmp folder: {}'.format(tmp_folder))

    try:
        input_queue = mp.Queue(maxsize=threads)
        resubmission_queue = mp.Queue()
        output_queue = mp.Queue()
        exception_queue = mp.Queue()
        monitor = Monitor()

        monitor.set_submitting(True)

        worker_pool = mp.Pool(threads, _iterative_mapping_worker,
                              (mapper, input_queue, tmp_folder, output_queue,
                               resubmission_queue, monitor, exception_queue))

        t_resub = threading.Thread(target=_resubmissions_to_queue, args=(resubmission_queue, tmp_folder,
                                                                         batch_size, input_queue, monitor,
                                                                         step_size, min_size,
                                                                         exception_queue, worker_pool))
        t_resub.daemon = True
        t_resub.start()

        t_sub = threading.Thread(target=_fastq_to_queue, args=(fastq_file, tmp_folder,
                                                               batch_size, input_queue,
                                                               monitor, exception_queue,
                                                               worker_pool))
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
                            if sam_counter > 0 and line.startswith('@'):
                                continue
                            if line.rstrip() == '':
                                continue
                            o.write(line)

                    os.remove(partial_sam_file)

                    sam_counter += 1
                    logger.debug('Got {}/{} SAM files'.format(sam_counter, monitor.value()))
                except Empty:
                    pass

        t_sub.join()
        t_resub.join()

        if convert_to_bam:
            logger.info("Converting intermediate SAM file to BAM ({})".format(sam_file))
            bam_command = ['samtools', 'view', '-b', '-o', sam_file, intermediate_sam_file]
            res = subprocess.call(bam_command)
            if res == 0:
                logger.info("Success, removing intermediate.")
                os.remove(intermediate_sam_file)
            else:
                logger.error("Could not convert to BAM, but your output is still in {}".format(intermediate_sam_file))
    finally:
        logger.debug(tmp_folder)
        shutil.rmtree(tmp_folder)
