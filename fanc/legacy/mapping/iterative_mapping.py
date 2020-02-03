import os.path
import tempfile
import shutil
import subprocess
import multiprocessing as mp
import re
from Bio import Restriction
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from collections import defaultdict
import glob
from queue import Empty
from future.utils import string_types
import pysam
from fanc.tools.files import fastq_reader
import logging
logger = logging.getLogger(__name__)


def ligation_site_pattern(restriction_enzyme):
    if isinstance(restriction_enzyme, string_types):
        restriction_enzyme = getattr(Restriction, restriction_enzyme)

    cut_pattern = restriction_enzyme.elucidate()

    left_side = []
    right_side = []
    for character in cut_pattern:
        if not (len(left_side) > 0 and left_side[-1] == '_'):
            if not character == '^':
                left_side.append(character)

        if character == '^' or len(right_side) > 0:
            if not character == '_':
                right_side.append(character)

    left_side_re = "".join(left_side[:-1])
    right_side_re = "".join(right_side[1:])

    left_side_re = left_side_re.replace('N', '[ACGT]')
    right_side_re = right_side_re.replace('N', '[ACGT]')

    pattern = "(^.+?" + left_side_re + ')' + right_side_re
    return re.compile(pattern)


class SequenceMapper(object):
    """
    Abstract class for mapping sequencing reads and
    comparing alignment qualities.
    """
    PERFECT_ALIGNMENT = 0
    IMPROVABLE_ALIGNMENT = 1
    BAD_ALIGNMENT = 2

    def __init__(self):
        pass

    def alignment_quality(self, alignment):
        """
        Get the quality of the given alignment to the
        reference genome.

        :param alignment: List of fields in the alignment
                          of a SAM file. See `SAM format
                          <https://samtools.github.io/hts-specs/SAMv1.pdf/>`_
        :return: :attr:`~SequenceMapper.PERFECT_ALIGNMENT` if the alignment
                 satisfies an internal set of quality standards;
                 :attr:`~SequenceMapper.IMPROVABLE_ALIGNMENT` if the alignment
                 does not satisfy an internal set of quality standards, but
                 a different alignment strategy could improve this result
                 :attr:`~SequenceMapper.BAD_ALIGNMENT` if the alignment
                 does not satisfy an internal set of quality standards
        """
        pass

    def get_better_alignment(self, alignment1, alignment2):
        """
        Return the better of the two alignments given an internal
        set of quality standards.

        :param alignment1: List of fields in the alignment
                           of a SAM file. See `SAM format
                           <https://samtools.github.io/hts-specs/SAMv1.pdf/>`_
        :param alignment2: List of fields in the alignment
                           of a SAM file. See `SAM format
                           <https://samtools.github.io/hts-specs/SAMv1.pdf/>`_
        :return: The better of the two alignments
        """
        pass

    def map(self, fastq):
        """
        Map a fastq file to a reference genome.

        :param fastq: FASTQ file
        :return: Tuple (header, alignments), where header is a
                 list of header lines and alignments is a list
                 of lists representing SAM alignments. Fields
                 in each SAM list correspond to fields in the
                 alignment of a SAM file. See `SAM format
                 <https://samtools.github.io/hts-specs/SAMv1.pdf/>`_
        """
        pass


class Bowtie2Mapper(SequenceMapper):
    """
    :class:`~SequenceMapper` for Bowtie2 aligner.
    """
    def __init__(self, index, executable="bowtie2", options=('--very-sensitive', '--no-unal'),
                 threads=1, quality_cutoff=0, memory_map=False, verbose=False):
        """
        Initialize Bowtie2 aligner with reference genome index
        and options.

        :param index: Location of the reference genome Bowtie2 index
                      in the form /path/to/index/prefix
        :param executable: Location of the bowtie2 executable. Leave
                           default to search in PATH.
        :param options: List of options passed the bowtie2
        :param threads: Number of threads bowtie2 should use.
        :param quality_cutoff: Cutoff to distinguish improvable from
                               perfect alignments
        """
        SequenceMapper.__init__(self)
        self.index = index
        self.threads = threads
        self.executable = executable
        self.options = options
        self.quality_cutoff = quality_cutoff
        self.memory_map = memory_map
        self.verbose = verbose

    def threads(self, threads):
        """
        Set the number of threads used by bowtie2
        """
        self.threads = threads

    def options(self, options):
        """
        Set a list of options for the bowtie2 executable
        """
        self.options = options

    def quality_cutoff(self, quality_cutoff):
        """
        Set the mapping quality cutoff
        """
        self.quality_cutoff = quality_cutoff

    def verbose(self, verbose):
        """
        Set the mapping quality cutoff
        """
        self.verbose = verbose

    def alignment_quality(self, alignment):
        """
        Determine the alignment quality of a SAM entry.

        See :func:`~SequenceMapper.alignment_quality` for details.

        :param alignment: List of fields in the alignment
                          of a SAM file. See `SAM format
                          <https://samtools.github.io/hts-specs/SAMv1.pdf/>`_
        :return: See :func:`~SequenceMapper.alignment_quality`
        """
        for i in range(11, len(alignment)):
            if alignment[i].startswith(b'XS:'):
                return SequenceMapper.BAD_ALIGNMENT

        try:
            if int(alignment[4]) < self.quality_cutoff:
                return SequenceMapper.IMPROVABLE_ALIGNMENT
        except (IndexError, ValueError):
            return SequenceMapper.IMPROVABLE_ALIGNMENT

        return SequenceMapper.PERFECT_ALIGNMENT

    def get_better_alignment(self, sam_read1, sam_read2):
        """
        Return the alignment with the larger mapq score.
        """
        try:
            mapq1 = int(sam_read1[4])
        except (IndexError, ValueError):
            return sam_read2

        try:
            mapq2 = int(sam_read2[4])
        except (IndexError, ValueError):
            return sam_read1

        if mapq1 >= mapq2:
            return sam_read1
        return sam_read2

    def map(self, fastq, output=None):
        """
        Map fastq to the bowtie2 index.

        If output is not None, output will be written directly to file.
        Otherwise output will be stored in a (header, alignments) tuple.

        See :func:`~SequenceMapper.map` for details
        """
        mapping_command = [self.executable, '-p', str(self.threads),
                           '-x', self.index, '-U', fastq]
        if output is not None:
            mapping_command += ['-S', output]
        mapping_command += self.options

        if self.memory_map:
            mapping_command += ['--mm']

        logger.debug("Mapping command: %s" % " ".join(mapping_command))

        if not self.verbose:
            stderr = open(os.devnull, 'w')
        else:
            stderr = subprocess.PIPE
        if output is not None:
            subprocess.call(mapping_command, stdout=subprocess.PIPE, stderr=stderr)
            if not self.verbose:
                stderr.close()
            return None, output

        mapping_process = subprocess.Popen(mapping_command, stdout=subprocess.PIPE, stderr=stderr)

        header = []
        alignments = defaultdict(list)
        while True:
            line = mapping_process.stdout.readline()
            if line != b'':
                if line.startswith(b'@'):
                    header.append(line)
                else:
                    fields = line.split(b'\t')
                    alignments[fields[0]].append(fields)
            else:
                break
        if not self.verbose:
            stderr.close()

        logger.debug("Aligned reads: %d" % len(alignments))

        return header, alignments


def _iteratively_map_reads(file_name, mapper=None, steps=None, min_read_length=None, step_size=2,
                           work_dir=None, output_file=None, write_header=False):
    """
    Iteratively map reads in a FASTQ or gzipped FASTQ file.

    :param file_name: Location of the FASTQ or gzipped FASTQ file.
    :param mapper: A :class:`~SequenceMapper` instance.
    :param steps: An iterable of read lengths to map. Overrides
                  min_read_length and step_size.
    :param min_read_length: Minimum read length to start iterative
                            mapping.
    :param step_size: Base pairs by which to extend unmapped reads
                      at each iteration.
    :param work_dir: Working directory. If None assumes the same
                     directory as FASTQ file.
    :param output_file: If provided, and write_header is True,
                        creates or overwrites a SAM file at this
                        location. If write_header is False, appends
                        alignments to this file.
    :param write_header: Determines if SAM header should be written
                         to file. Ignored if output_file is None.
    :return: Tuple (header, alignments), where header is a
             list of header lines and alignments is a list
             of lists representing SAM alignments. Fields
             in each SAM list correspond to fields in the
             alignment of a SAM file. See `SAM format
             <https://samtools.github.io/hts-specs/SAMv1.pdf/>`_
    """
    if work_dir is None:
        work_dir = os.path.split(file_name)[0]

    reader = fastq_reader(file_name)

    if steps is None:
        logger.debug("Finding maximum read length...")

        max_len = 0
        with reader(file_name, 'r') as fastq:
            for title, seq, qual in FastqGeneralIterator(fastq):
                max_len = max(len(seq), max_len)

        logger.debug("Maximum read length: %d" % max_len)

        if min_read_length is None:
            min_read_length = max_len

        steps = list(range(min_read_length, max_len+1, step_size))
        if len(steps) == 0 or steps[-1] != max_len:
            steps.append(max_len)

        ixs = [0]
        current = 1
        for i in range(len(steps)-1):
            if i % 2 == 0:
                ixs.append(-1*current)
            else:
                ixs.append(current)
                current += 1
        steps = [steps[ix] for ix in ixs]
    else:
        if len(steps) <= 1:
            step_size = 0
        else:
            step_size = abs(steps[0]-steps[1])
            for i in range(2, len(steps)):
                step_size = min(step_size, abs(steps[0]-steps[i]))
            min_read_length = min(steps)

    trimmed_file = work_dir + '/trimmed.fastq'
    perfect_alignments = {}
    improvable_alignments = {}
    for i, size in enumerate(steps):
        fastq_counter = 0
        with reader(file_name, 'rt') as fastq:
            with open(trimmed_file, 'w') as trimmed:
                for title, seq, qual in FastqGeneralIterator(fastq):
                    name = title.split(" ")[0]
                    if name not in perfect_alignments or name in improvable_alignments:
                        if len(seq)+step_size >= size:
                            trimmed.write("@%s\n%s\n+\n%s\n" % (title, seq[:size], qual[:size]))
                            fastq_counter += 1

        logger.debug("Sending %d reads to the next iteration (length %d)" % (fastq_counter, size))

        header, alignments_trimmed = mapper.map(trimmed_file)

        # merge alignments
        try:
            # noinspection PyCompatibility
            alignments_trimmed_items = alignments_trimmed.iteritems()
        except AttributeError:
            alignments_trimmed_items = alignments_trimmed.items()

        for name, fields_array in alignments_trimmed_items:
            worst_quality = SequenceMapper.PERFECT_ALIGNMENT
            for fields in fields_array:
                worst_quality = max(worst_quality, mapper.alignment_quality(fields))

            if worst_quality == SequenceMapper.PERFECT_ALIGNMENT:
                perfect_alignments[name] = fields_array
                if name in improvable_alignments:
                    del improvable_alignments[name]
            elif worst_quality == SequenceMapper.IMPROVABLE_ALIGNMENT:
                if name in improvable_alignments:
                    if len(improvable_alignments[name]) < len(fields_array):
                        improvable_alignments[name] = fields_array
                    elif len(improvable_alignments[name]) == len(fields_array):
                        new_fields_array = []
                        for j in range(len(fields_array)):
                            fields = fields_array[j]
                            existing_fields = improvable_alignments[name][j]
                            new_fields_array.append(mapper.get_better_alignment(fields, existing_fields))
                else:
                    improvable_alignments[name] = fields_array

        logger.debug("Resubmitting %d improvable alignments" % len(improvable_alignments))

    # clean
    os.unlink(trimmed_file)

    # merge alignments into one
    perfect_alignments.update(improvable_alignments)
    if output_file is not None:
        mode = 'wb' if write_header else 'ab'
        # flush results to file
        with open(output_file, mode) as o:
            if write_header:
                for header_line in header:
                    o.write(header_line)

            try:
                # noinspection PyCompatibility
                perfect_alignments_items = perfect_alignments.iteritems()
            except AttributeError:
                perfect_alignments_items = perfect_alignments.items()
            for _, fields_array in perfect_alignments_items:
                for fields in fields_array:
                    alignment_line = b'\t'.join(fields)
                    o.write(alignment_line)
        return output_file

    return header, perfect_alignments


def split_iteratively_map_reads(input_file, output_file, index_path, work_dir=None, quality_cutoff=30,
                                batch_size=1000000, threads=1, min_size=25, step_size=2, copy=False,
                                restriction_enzyme=None, adjust_batch_size=False, mapper=None,
                                bowtie_parallel=True, memory_map=False):

    check_path = os.path.expanduser(index_path)
    if check_path.endswith('.'):
        check_path = check_path[:-1]
    for i in range(1, 5):
        if not os.path.exists(check_path + '.{}.bt2'.format(i)):
            raise ValueError("Cannot find bowtie2 path!")
    for i in range(1, 3):
        if not os.path.exists(check_path + '.rev.{}.bt2'.format(i)):
            raise ValueError("Bowtie2 index incomplete, check index files for completeness.")

    bowtie_threads, worker_threads = (threads, 1) if bowtie_parallel else (1, threads)

    if work_dir is not None:
        work_dir = tempfile.mkdtemp(dir=os.path.expanduser(work_dir))
    else:
        work_dir = tempfile.mkdtemp()

    header, alignments = None, None
    working_output_file = None
    working_input_file = None
    try:
        logger.info("Working directory: %s" % work_dir)

        if index_path.endswith('.'):
            index_path = index_path[:-1]

        if copy:
            working_input_file = work_dir + '/' + os.path.basename(input_file)
            shutil.copyfile(input_file, working_input_file)
            working_output_file = work_dir + '/' + os.path.basename(output_file)
            os.makedirs(work_dir + '/index')
            index_base = os.path.basename(index_path)

            for file_name in glob.glob(index_path + '*.bt2'):
                shutil.copy(file_name, work_dir + '/index')
            index_path = work_dir + '/index/' + index_base
        else:
            working_input_file = input_file
            working_output_file = output_file

        reader = fastq_reader(working_input_file)
        working_file = gzip.open(work_dir + '/full_reads_0.fastq.gz', 'wb')
        working_files = [working_file.name]

        if mapper is None:
            mapper = Bowtie2Mapper(index=index_path, quality_cutoff=quality_cutoff, threads=bowtie_threads,
                                   memory_map=memory_map)

        logger.info("Splitting files...")
        re_pattern = None
        if restriction_enzyme is not None:
            re_pattern = ligation_site_pattern(restriction_enzyme)

        if adjust_batch_size and threads > 1:
            logger.info("Counting lines to adjust batch size...")
            with reader(working_input_file, 'rb') as fastq:
                n_lines = sum(1 for _ in fastq)/4

            if n_lines/threads < batch_size*threads:
                batch_size = int(n_lines/threads)+threads
                logger.info("Adjusted batch size to: %d" % batch_size)

        def _mapping_process_with_queue(input_queue, output_queue):
            while True:
                logger.debug("Waiting for input...")
                p_number, file_name, mapper, min_size, max_length, step_size, work_dir = input_queue.get(True)
                steps = list(range(min_size, max_length+1, step_size))
                if len(steps) == 0 or steps[-1] != max_length:
                    steps.append(max_length)

                ixs = [0]
                current = 1
                for i in range(len(steps)-1):
                    if i % 2 == 0:
                        ixs.append(-1*current)
                    else:
                        ixs.append(current)
                        current += 1
                steps = [steps[ix] for ix in ixs]

                logger.debug("Mapping %s" % file_name)
                partial_output_file = work_dir + '/mapped_reads_' + str(p_number) + '.sam'
                process_work_dir = work_dir + "/mapping_%d/" % p_number
                os.makedirs(process_work_dir)

                _iteratively_map_reads(file_name, mapper, steps, None, None, process_work_dir,
                                       partial_output_file, True)

                logger.debug("Done mapping %s" % file_name)

                os.unlink(file_name)
                shutil.rmtree(process_work_dir)

                output_queue.put(partial_output_file)

        input_queue = mp.Queue()
        output_queue = mp.Queue()
        worker_pool = mp.Pool(worker_threads, _mapping_process_with_queue, (input_queue, output_queue))

        max_length = 0
        trimmed_count = 0
        batch_count = 0
        batch_reads_count = 0
        output_count = 0
        o = None

        try:
            with reader(working_input_file, 'rt') as fastq:
                for title, seq, qual in FastqGeneralIterator(fastq):
                    # check if ligation junction is in this read
                    if re_pattern is not None:
                        m = re_pattern.search(seq)
                        if m is not None:
                            seq = m.group(1)
                            qual = qual[:len(seq)]
                            trimmed_count += 1
                    max_length = max(max_length, len(seq))
                    line = "@{}\n{}\n+\n{}\n".format(title, seq, qual)
                    working_file.write(line.encode())
                    batch_reads_count += 1

                    if batch_reads_count > batch_size:
                        # reset
                        batch_reads_count = 0
                        working_file.close()

                        # prepare process
                        input_queue.put((batch_count, working_file.name, mapper, min_size,
                                         max_length, step_size, work_dir))

                        # write output if any
                        logger.debug("Merging output files...")
                        while True:
                            logger.debug("Trying to collect results")
                            try:
                                partial_output_file = output_queue.get(False)
                                logger.debug("Processing %s..." % partial_output_file)
                                with pysam.AlignmentFile(partial_output_file) as p:
                                    if o is None:
                                        if os.path.splitext(output_file)[1] == '.bam':
                                            o = pysam.AlignmentFile(working_output_file, 'wb', template=p)
                                        else:
                                            o = pysam.AlignmentFile(working_output_file, 'wh', template=p)
                                    for alignment in p:
                                        o.write(alignment)
                                output_count += 1
                                os.unlink(partial_output_file)
                            except Empty:
                                logger.debug("No results found")
                                break

                        max_length = 0
                        batch_count += 1
                        working_file = gzip.open(work_dir + '/full_reads_' + str(batch_count) + '.fastq.gz', 'wb')
                        working_files.append(working_file.name)

                working_file.close()

                # prepare last process
                if batch_reads_count > 0:
                    # prepare process
                    input_queue.put((batch_count, working_file.name, mapper, min_size,
                                     max_length, step_size, work_dir))

            if copy:
                os.unlink(working_input_file)

            logger.info("Trimmed %d reads at ligation junction" % trimmed_count)

            # merge files
            while output_count < batch_count + 1:
                partial_output_file = output_queue.get(True)
                logger.info("Processing %s..." % partial_output_file)
                with pysam.AlignmentFile(partial_output_file, 'r') as p:
                    if o is None:
                        if os.path.splitext(output_file)[1] == '.bam':
                            o = pysam.AlignmentFile(working_output_file, 'wb', template=p)
                        else:
                            o = pysam.AlignmentFile(working_output_file, 'wh', template=p)

                    for alignment in p:
                        o.write(alignment)
                output_count += 1
                logger.info("%d/%d" % (output_count, batch_count + 1))
                os.unlink(partial_output_file)
        finally:
            try:
                o.close()
            except AttributeError:
                pass

        logger.info("All done.")
    finally:
        # clean up

        # move output file into place
        if copy:
            try:
                shutil.move(working_output_file, output_file)
            except (OSError, IOError):
                logger.error("Cannot move working output file to required location.")
        # remove temporary directory
        shutil.rmtree(work_dir)

    return header, alignments