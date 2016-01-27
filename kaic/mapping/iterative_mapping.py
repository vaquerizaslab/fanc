#!/usr/bin/env python

import argparse
import os.path
import tempfile
import shutil
import logging
import subprocess
import multiprocessing as mp
from functools import partial
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from collections import defaultdict
logging.basicConfig(level=logging.DEBUG)


def _get_fastq_reader(file_name):
    """
    Return appropriate 'open' method by filename extension.

    :param file_name: Filename of the FASTQ or gzipped FASTQ
                      file.
    :return: gzip.open for '.gz' and '.gzip' files, 'os.open'
             otherwise.
    """
    input_extension = os.path.splitext(file_name)[1]
    if input_extension == ".gz" or input_extension == ".gzip":
        return gzip.open
    return open


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
                 threads=1, quality_cutoff=0, verbose=False):
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
        for i in xrange(11, len(alignment)):
            if alignment[i].startswith("XS:"):
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

        logging.info("Mapping command: %s" % " ".join(mapping_command))

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
            if line != '':
                if line.startswith("@"):
                    header.append(line)
                else:
                    fields = line.split("\t")
                    alignments[fields[0]].append(fields)
            else:
                break
        if not self.verbose:
            stderr.close()

        logging.info("Aligned reads: %d" % len(alignments))

        return header, alignments


def iteratively_map_reads(file_name, mapper=None, min_read_length=None, step_size=2,
                          work_dir=None, output_file=None, write_header=False):
    """
    Iteratively map reads in a FASTQ or gzipped FASTQ file.

    :param file_name: Location of the FASTQ or gzipped FASTQ file.
    :param mapper: A :class:`~SequenceMapper` instance.
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

    reader = _get_fastq_reader(file_name)
    max_len = 0
    with reader(file_name, 'r') as fastq:
        for title, seq, qual in FastqGeneralIterator(fastq):
            max_len = len(seq)
            break

    logging.info("Maximum read length: %d" % max_len)

    if min_read_length is None:
        min_read_length = max_len

    steps = list(xrange(min_read_length, max_len+1, step_size))
    if steps[-1] != max_len:
        steps.append(max_len)

    perfect_alignments = {}
    improvable_alignments = {}
    for i, size in enumerate(steps):
        trimmed_file = work_dir + '/trimmed.fastq'

        fastq_counter = 0
        with reader(file_name, 'r') as fastq:
            with open(trimmed_file, 'w') as trimmed:
                for title, seq, qual in FastqGeneralIterator(fastq):
                    name = title.split(" ")[0]
                    if name not in perfect_alignments or name in improvable_alignments:
                        trimmed.write("@%s\n%s\n+\n%s\n" % (title, seq[:size], qual[:size]))
                        fastq_counter += 1

        logging.info("Sending %d reads to the next iteration (length %d)" % (fastq_counter, size))

        header, alignments_trimmed = mapper.map(trimmed_file)

        # merge alignments
        for name, fields_array in alignments_trimmed.iteritems():
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
                        for j in xrange(len(fields_array)):
                            fields = fields_array[j]
                            existing_fields = improvable_alignments[name][j]
                            new_fields_array.append(mapper.get_better_alignment(fields, existing_fields))
                else:
                    improvable_alignments[name] = fields_array

        logging.info("Resubmitting %d improvable alignments" % len(improvable_alignments))

    # merge alignments into one
    perfect_alignments.update(improvable_alignments)
    if output_file is not None:
        mode = 'w' if write_header else 'a'
        # flush results to file
        with open(output_file, mode) as o:
            if write_header:
                for header_line in header:
                    o.write(header_line)
            for _, fields_array in perfect_alignments.iteritems():
                for fields in fields_array:
                    alignment_line = "\t".join(fields)
                    o.write(alignment_line)
        return output_file

    return header, perfect_alignments


def split_iteratively_map_reads(input_file, output_file, index_path, work_dir=None, quality_cutoff=30,
                                batch_size=250000, threads=1, min_size=25, step_size=2, copy=False):
    if work_dir is not None:
        work_dir = tempfile.mkdtemp(dir=os.path.expanduser(work_dir))
    else:
        work_dir = tempfile.mkdtemp()

    header, alignments = None, None
    working_output_file = None
    working_input_file = None
    try:
        logging.info("Working directory: %s" % work_dir)

        if copy:
            working_input_file = work_dir + '/' + os.path.basename(input_file)
            shutil.copyfile(input_file, working_input_file)
            working_output_file = work_dir + '/' + os.path.basename(output_file)
        else:
            working_input_file = input_file
            working_output_file = output_file

        reader = _get_fastq_reader(working_input_file)
        working_file = gzip.open(work_dir + '/full_reads_0.fastq.gz', 'w')
        working_files = [working_file.name]

        mapper = Bowtie2Mapper(index=index_path, quality_cutoff=quality_cutoff, threads=1)

        logging.info("Splitting files...")
        batch_count = 0
        batch_reads_count = 0
        with reader(working_input_file, 'r') as fastq:
            for title, seq, qual in FastqGeneralIterator(fastq):
                if batch_reads_count <= batch_size:
                    line = "@%s\n%s\n+\n%s\n" % (title, seq, qual)
                    working_file.write(line)
                    batch_reads_count += 1
                else:
                    # reset
                    batch_reads_count = 0
                    working_file.close()
                    batch_count += 1
                    working_file = gzip.open(work_dir + '/full_reads_' + str(batch_count) + '.fastq.gz', 'w')
                    working_files.append(working_file.name)

            working_file.close()

        logging.info("Starting to map...")
        output_files = []
        processes = []
        for i, working_file in enumerate(working_files):
            partial_output_file = work_dir + '/mapped_reads_' + str(i) + '.sam'
            output_files.append(partial_output_file)
            process_work_dir = work_dir + "/mapping_%d/" % i
            os.makedirs(process_work_dir)
            processes.append(mp.Process(target=iteratively_map_reads, args=(working_file, mapper,
                                                                            min_size, step_size, process_work_dir,
                                                                            partial_output_file, True)))
        current_processes = []
        for i, p in enumerate(processes):
            current_processes.append(p)
            if len(current_processes) > threads or i == len(processes)-1:
                for cp in current_processes:
                    cp.start()
                for cp in current_processes:
                    cp.join()

        print output_files
        print working_files

        # merge files
        with open(working_output_file, 'w') as o:
            for i, partial_output_file in enumerate(output_files):
                with open(partial_output_file, 'r') as p:
                    for line in p:
                        if line.startswith("@") and i > 0:
                            continue
                        o.write(line)

        if os.path.splitext(output_file)[1] == '.bam':
            logging.info("Converting to BAM...")
            success = False
            o = open(work_dir + '/output.bam', 'w')
            try:
                bam_command = ['samtools', 'view', '-bS', working_output_file]
                logging.info("BAM conversion command: %s" % " ".join(bam_command))
                exit_code = subprocess.call(bam_command, stdout=o)
                success = True if exit_code == 0 else False
            except OSError:
                logging.error("Cannot convert to BAM, samtools not in PATH!")
            finally:
                o.close()

            if success:
                if not copy:
                    shutil.move(o.name, output_file)
                working_output_file = o.name
            else:
                logging.info("BAM conversion failed.")

        logging.info("All done.")
    finally:
        # clean up

        # move output file into place
        if copy:
            try:
                shutil.move(working_output_file, output_file)
            except (OSError, IOError):
                logging.error("Cannot move working output file to required location.")
        # remove temporary directory
        shutil.rmtree(work_dir)

    return header, alignments