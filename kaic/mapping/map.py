import os
import subprocess
import multiprocessing as mp
from queue import Empty, Full
import msgpack
import traceback
import gzip
import threading
from kaic.tools.general import which
from kaic.tools.general import RareUpdateProgressBar
import logging
logger = logging.getLogger(__name__)


class Mapper(object):
    def __init__(self):
        self.resubmit_unmappable = True

    def map(self, fastq_strings):
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

    def map(self, fastq_strings):
        bowtie2_command = [self._path, '-x', self.index, '-U', '-', '--no-unal',
                           '--threads', str(self.threads)] + self.args
        f_null = open(os.devnull, 'w')
        bp = subprocess.Popen(bowtie2_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=f_null, universal_newlines=True)
        output = bp.communicate(input=''.join(fastq_strings))[0]
        f_null.close()
        bp.kill()
        return output

    def _resubmit(self, sam_fields):
        if int(sam_fields[4]) < self.min_quality:
            return True
        return False


class SimpleBowtie2Mapper(Bowtie2Mapper):
    def __init__(self, bowtie2_index, min_quality=30, additional_arguments=(),
                 threads=1, _bowtie2_path='bowtie2'):
        Bowtie2Mapper.__init__(self, bowtie2_index, min_quality=min_quality,
                               additional_arguments=additional_arguments,
                               threads=threads,
                               _bowtie2_path=_bowtie2_path)
        self.resubmit_unmappable = False

    def _resubmit(self, sam_fields):
        return False


def _trim_read(input_seq, step_size=5, min_size=25):
    name, seq, plus, qual, blank = input_seq.split("\n")
    if len(seq) == min_size:
        raise ValueError("Already reached minimum size, cannot truncate read further")
    if len(seq) - step_size < min_size:
        final_length = min_size
    else:
        final_length = len(seq) - step_size
    seq = seq[:final_length]
    qual = qual[:final_length]
    return len(seq), "\n".join([name, seq, plus, qual, blank])


def _iterative_mapping_worker(mapper, input_queue, resubmission_queue, output_queue, header_queue,
                              exception_queue, batch_size=100, min_size=25, step_size=5):

    try:
        lines = dict()
        while True:
            try:
                # first process resubmissions
                input_string = resubmission_queue.get(False)
            except Empty:
                input_string = input_queue.get(True)

            if input_string is not None:
                s = msgpack.loads(input_string)
                input_string = s.decode() if isinstance(s, bytes) else s
                input_name = input_string.split("\n")[0].split()[0][1:]
                lines[input_name] = input_string

            if len(lines) == batch_size or (len(lines) > 0 and input_string is None):
                output = mapper.map(lines.values())
                sam_lines = output.split("\n")

                remaining = set(lines.keys())
                seen = set()
                header = ""

                line_iter = iter(sam_lines)
                try:
                    line = next(line_iter)
                    # process header
                    while line.startswith("@"):
                        header += "{}\n".format(line)
                        line = next(line_iter)
                    # send header to main thread
                    if header_queue.empty():
                        header_queue.put(header)
                    # process alignment lines
                    while True:
                        if line != '':
                            fields = line.split("\t")
                            name = fields[0]
                            seen.add(name)

                            if not mapper.resubmit(fields):
                                remaining.remove(name)
                                output_queue.put(msgpack.dumps(line))
                        line = next(line_iter)
                except StopIteration:
                    pass

                # process remaining (unaligned or resubmitted) reads
                for name in remaining:
                    # if the read was not in SAM output it was unmappable
                    if name not in seen and not mapper.resubmit_unmappable:
                        output_queue.put(None)
                        continue

                    # else check if we reached minimum length
                    try:
                        new_len, new_seq = _trim_read(lines[name], step_size, min_size)
                    except ValueError:
                        output_queue.put(None)
                        continue

                    resubmission_queue.put(msgpack.dumps(new_seq))

                lines = dict()
            if input_string is None:
                input_queue.put(None)
                if resubmission_queue.empty() and input_queue.empty():
                    logger.info('Worker thread is done, no more reads to process.')
                    break
    except Exception:
        import sys
        exception_queue.put("".join(traceback.format_exception(*sys.exc_info())))


def _fastq_to_queue(fastq_file, input_queue, counter_queue,
                    exception_queue=None, worker_pool=None):
    if fastq_file.endswith('.gz') or fastq_file.endswith('.gzip'):
        open_file = gzip.open
    else:
        open_file = open

    current_fastq = None
    read_counter = 0
    with open_file(fastq_file, 'r') as f:
        for i, line in enumerate(f):
            line = line.decode() if isinstance(line, bytes) else line
            if i % 4 == 0:
                if current_fastq is not None:
                    if exception_queue is not None:
                        while True:
                            try:
                                input_queue.put(msgpack.dumps(current_fastq), True, 5)
                                read_counter += 1
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
                        input_queue.put(msgpack.dumps(current_fastq), True, 5)
                        read_counter += 1
                current_fastq = ''
            current_fastq += line

        input_queue.put(msgpack.dumps(current_fastq))
        input_queue.put(None)
        read_counter += 1
    counter_queue.put(read_counter)


def iterative_mapping(fastq_file, sam_file, mapper, threads=1, min_size=25, step_size=5,
                      batch_size=20000, header_timeout=1800, max_queue_size=20000):
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
    :param threads: Number of mapper threads to use in parallel.
    :param min_size: Minimum length of read for which an alignment is attempted.
    :param step_size: Number of base pairs by which to truncate read.
    :param batch_size:
    :param header_timeout:
    :param max_queue_size:
    :return:
    """
    input_queue = mp.Queue(maxsize=max_queue_size)
    resubmission_queue = mp.Queue()
    output_queue = mp.Queue()
    header_queue = mp.Queue()
    exception_queue = mp.Queue()
    counter_queue = mp.Queue()

    worker_pool = mp.Pool(threads, _iterative_mapping_worker,
                          (mapper, input_queue, resubmission_queue, output_queue,
                           header_queue, exception_queue, batch_size, min_size,
                           step_size))

    t = threading.Thread(target=_fastq_to_queue, args=(fastq_file, input_queue,
                                                       counter_queue, exception_queue,
                                                       worker_pool))
    t.daemon = True
    t.start()

    logger.debug("Waiting for first mapping process to return to grab header...")
    try:
        header = header_queue.get(True, header_timeout)
    except Empty:
        raise RuntimeError("Waited too long ({} min) for first mapping process to end".format(header_timeout/60/60))

    def _open_output(file_name, bam=True):
        if not bam:
            return open(file_name, 'w')
        process = subprocess.Popen(['samtools', 'view', '-b', '-o', file_name, '-'],
                                   stdin=subprocess.PIPE,
                                   universal_newlines=True)
        return process.stdin

    def write_next_result_from_queue():
        try:
            exc = exception_queue.get(block=False)
        except Empty:
            pass
        else:
            worker_pool.terminate()
            raise Exception(exc)
        try:
            result = output_queue.get(block=True, timeout=10)
            if result is not None:
                line = msgpack.loads(result)
                line = line.decode() if isinstance(line, bytes) else line
                o.write("{}\n".format(line))
        except Empty:
            pass

    logger.info("Starting to output alignments to file {}".format(sam_file))
    with _open_output(sam_file, sam_file.endswith('.bam')) as o:
        o.write(header)
        alignment_counter = 0

        while counter_queue.empty():
            write_next_result_from_queue()
            alignment_counter += 1
        else:
            read_counter = counter_queue.get(True)

        logger.info("Reading FASTQ file {} complete".format(fastq_file))
        with RareUpdateProgressBar(max_value=read_counter) as pb:
            while alignment_counter < read_counter:
                write_next_result_from_queue()
                alignment_counter += 1
                pb.update(alignment_counter)
    t.join()
