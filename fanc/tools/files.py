import gzip
import logging
import multiprocessing
import os.path
import pyBigWig
import random
import shutil
import string
import subprocess
import sys
import tempfile
import threading
from collections import defaultdict

import numpy as np
import pysam
import tables as t
from Bio import SeqIO
from future.utils import string_types

from .general import mkdir, which
from .sambam import natural_cmp


# configure logging
logger = logging.getLogger(__name__)


def get_number_of_lines(file_name):
    with open(file_name, 'r') as f:
        return sum(1 for _ in f)


def tmp_file_name(tmpdir, prefix='tmp_fanc', extension='h5'):
    with tempfile.NamedTemporaryFile('w', prefix=prefix, dir=tmpdir, suffix=extension) as f:
        name = f.name
    return name


def random_name(length=6):
    return ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(length))


def create_or_open_pytables_file(file_name=None, mode='a'):
    in_memory = False
    if file_name is None:
        file_name = random_name()
        in_memory = True
    
    # check if already is pytables File
    if isinstance(file_name, t.file.File):
        return file_name
    
    # check if is existing
    if mode == 'a' and os.path.isfile(file_name):
        try:
            f = t.open_file(file_name, "r", chunk_cache_size=270536704, chunk_cache_nelmts=2084)
            f.close()
        except t.HDF5ExtError:
            raise ImportError("File exists and is not an HDF5 dict")
    
    if in_memory:
        return t.open_file(file_name, mode, driver="H5FD_CORE", driver_core_backing_store=0)
    else:
        # 256Mb cache
        return t.open_file(file_name, mode, chunk_cache_size=270536704, chunk_cache_nelmts=2084)


def is_sambam_file(file_name):
    file_name = os.path.expanduser(file_name)
    if not os.path.isfile(file_name):
        return False
    
    try:
        sb = pysam.AlignmentFile(file_name, 'rb')
        sb.close()
    except (ValueError, IOError):
        return False
    return True


def is_fasta_file(file_name):
    file_name = os.path.expanduser(file_name)
    if not os.path.isfile(file_name):
        return False

    try:
        is_fasta = True
        open_ = gzip.open if file_name.endswith('.gz') or file_name.endswith('.gzip') else open

        with open_(file_name, 'rt') as f:
            fastas = SeqIO.parse(f, 'fasta')

            try:
                next(fastas)
            except StopIteration:
                is_fasta = False
    except Exception:
        return False
        
    return is_fasta 
        

def copy_or_expand(input_file, output_file=None):
    # copy file if required
    input_path = os.path.expanduser(input_file)
    if output_file:
        output_path = os.path.expanduser(output_file)
        shutil.copy(input_path, output_path)
        if os.path.isdir(output_path):
            input_path = "%s/%s" % (output_path, os.path.basename(input_path))
        else:
            input_path = output_path
    return input_path


def read_chromosome_sizes(file_name):
    chrom_sizes = {}
    with open(os.path.expanduser(file_name), 'r') as chrom_sizes_file:
        for line in chrom_sizes_file:
            line = line.rstrip()
            if line != '':
                chromosome, chromosome_length = line.split("\t")
                chrom_sizes[chromosome] = int(chromosome_length)
    return chrom_sizes


def fastq_reader(file_name):
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


def gzip_splitext(file_name):
    basepath, extension = os.path.splitext(file_name)
    if extension == '.gz' or extension == '.gzip':
        basepath, extension2 = os.path.splitext(basepath)
        extension = extension2 + extension
    return basepath, extension


def _split_fastq_worker(fastq_file, output_folder, chunk_size=10000000, unzip=False, q=None):
    fastq_file = os.path.expanduser(fastq_file)
    output_folder = mkdir(output_folder)

    basepath, extension = gzip_splitext(fastq_file)
    basename = os.path.basename(basepath)

    fastq_open = fastq_reader(fastq_file)
    fastq_write = fastq_open if not unzip else open

    if unzip:
        extension = '.fastq'

    with fastq_open(fastq_file, 'r') as fastq:
        current_chunk_number = 0
        current_chunk_size = 0
        current_file_name = os.path.join(output_folder, '{}_{}{}'.format(current_chunk_number,
                                                                         basename, extension))
        current_split_file = fastq_write(current_file_name, 'w')
        for i, line in enumerate(fastq):
            if i % 4 == 0:
                if current_chunk_size >= chunk_size:
                    current_split_file.close()
                    q.put(current_file_name)
                    current_chunk_size = 0
                    current_chunk_number += 1
                    current_file_name = os.path.join(output_folder, '{}_{}{}'.format(current_chunk_number,
                                                                                     basename, extension))
                    current_split_file = fastq_write(current_file_name, 'w')
                current_chunk_size += 1
            current_split_file.write(line)
    current_split_file.close()
    q.put(current_file_name)
    q.put(None)


def _split_fastq_pair_worker(fastq_pair, output_folder, chunk_size=10000000, unzip=False, q=None):
    fastq_file_1 = os.path.expanduser(fastq_pair[0])
    fastq_file_2 = os.path.expanduser(fastq_pair[1])
    output_folder = mkdir(output_folder)

    basepath_1, extension_1 = gzip_splitext(fastq_file_1)
    basename_1 = os.path.basename(basepath_1)

    basepath_2, extension_2 = gzip_splitext(fastq_file_2)
    basename_2 = os.path.basename(basepath_2)

    fastq_open_1 = fastq_reader(fastq_file_1)
    fastq_open_2 = fastq_reader(fastq_file_2)

    fastq_write_1 = fastq_open_1 if not unzip else open
    fastq_write_2 = fastq_open_2 if not unzip else open

    if unzip:
        extension_1 = '.fastq'
        extension_2 = '.fastq'

    with fastq_open_1(fastq_file_1, 'r') as fastq_1:
        with fastq_open_2(fastq_file_2, 'r') as fastq_2:
            current_chunk_number = 0
            current_chunk_size = 0
            current_file_name_1 = os.path.join(output_folder, '{}_{}{}'.format(current_chunk_number,
                                                                               basename_1, extension_1))
            current_file_name_2 = os.path.join(output_folder, '{}_{}{}'.format(current_chunk_number,
                                                                               basename_2, extension_2))
            current_split_file_1 = fastq_write_1(current_file_name_1, 'w')
            current_split_file_2 = fastq_write_2(current_file_name_2, 'w')
            fastq2_iter = iter(fastq_2)
            for i, line_1 in enumerate(fastq_1):
                line_2 = next(fastq2_iter)
                if i % 4 == 0:
                    if current_chunk_size >= chunk_size:
                        current_split_file_1.close()
                        current_split_file_2.close()
                        q.put((current_file_name_1, current_file_name_2))
                        current_chunk_size = 0
                        current_chunk_number += 1
                        current_file_name_1 = os.path.join(output_folder, '{}_{}{}'.format(current_chunk_number,
                                                                                           basename_1, extension_1))
                        current_file_name_2 = os.path.join(output_folder, '{}_{}{}'.format(current_chunk_number,
                                                                                           basename_2, extension_2))
                        current_split_file_1 = fastq_write_1(current_file_name_1, 'w')
                        current_split_file_2 = fastq_write_2(current_file_name_2, 'w')
                    current_chunk_size += 1
                current_split_file_1.write(line_1)
                current_split_file_2.write(line_2)
    current_split_file_1.close()
    current_split_file_2.close()
    q.put((current_file_name_1, current_file_name_2))
    q.put(None)


def split_fastq(fastq_file, output_folder, chunk_size=10000000, parallel=True, unzip=False):
    if isinstance(fastq_file, string_types):
        _worker = _split_fastq_worker
    else:
        _worker = _split_fastq_pair_worker

    q = multiprocessing.Queue()
    t = None
    if parallel:
        t = threading.Thread(target=_worker, args=(
            fastq_file, output_folder
        ), kwargs={
            'chunk_size': chunk_size,
            'q': q,
            'unzip': unzip
        })
        t.daemon = True
        t.start()
    else:
        _worker(fastq_file, output_folder, chunk_size=chunk_size, q=q, unzip=unzip)

    while True:
        file_name = q.get()
        if file_name is None:
            break
        yield file_name

    q.close()
    q.join_thread()

    if parallel:
        t.join()


def split_sam(sam_file, output_folder, chunk_size=5000000):
    sam_file = os.path.expanduser(sam_file)
    output_folder = mkdir(output_folder)

    basepath, extension = gzip_splitext(sam_file)
    basename = os.path.basename(basepath)

    split_files = []
    with pysam.AlignmentFile(sam_file) as sambam:
        mode = 'wb' if os.path.splitext(sam_file)[1] == '.bam' else 'wh'

        current_chunk_number = 0
        current_chunk_size = 0
        current_file_name = os.path.join(output_folder, '{}_{}{}'.format(current_chunk_number,
                                                                         basename, extension))
        split_files.append(current_file_name)
        current_split_file = pysam.AlignmentFile(current_file_name, mode, template=sambam)
        for i, alignment in enumerate(sambam):
            if current_chunk_size >= chunk_size:
                current_split_file.close()
                current_chunk_size = 0
                current_chunk_number += 1
                current_file_name = os.path.join(output_folder, '{}_{}{}'.format(current_chunk_number,
                                                                                 basename, extension))
                split_files.append(current_file_name)
                current_split_file = pysam.AlignmentFile(current_file_name, mode, template=sambam)
            current_chunk_size += 1
            current_split_file.write(alignment)
    current_split_file.close()

    return split_files


def merge_sam(input_sams, output_sam, tmp=None):
    output_sam = os.path.expanduser(output_sam)
    first_sam = os.path.expanduser(input_sams[0])

    if tmp is not None:
        if isinstance(tmp, string_types):
            tmp = os.path.expanduser(tmp)
        elif tmp:
            tmp = tempfile.mkdtemp()

        if tmp:
            logger.info("Working from tmp dir {}".format(tmp))

    try:
        with pysam.AlignmentFile(first_sam) as fs:
            mode = 'wb' if os.path.splitext(output_sam)[1] == '.bam' else 'wh'
            with pysam.AlignmentFile(output_sam, mode, template=fs) as o:
                for input_sam in input_sams:
                    input_sam = os.path.expanduser(input_sam)
                    if tmp:
                        shutil.copy(input_sam, tmp)
                        input_sam = os.path.join(tmp, os.path.basename(input_sam))
                        logger.info("Copied input SAM to {}".format(input_sam))
                    with pysam.AlignmentFile(input_sam) as r:
                        for alignment in r:
                            o.write(alignment)
                    if tmp:
                        os.remove(input_sam)
    finally:
        if tmp:
            shutil.rmtree(tmp)
    return output_sam


def reads_with_same_qname(iterator, last_read=None):
    if last_read is not None:
        read = last_read
    else:
        read = next(iterator)
    qname = read.qname.encode('utf-8')

    reads = []
    try:
        next_read = read
        while natural_cmp(next_read.qname.encode('utf-8'), qname) == 0:
            if not next_read.is_unmapped:
                reads.append(next_read)
            next_read = next(iterator)
    except StopIteration:
        next_read = None
        if len(reads) == 0:
            raise

    return qname, reads, next_read


def split_sam_pairs(sam_file_1, sam_file_2, output_prefix,
                    chunk_size=10000000, check_sorted=True):
    """
    Form mate pairs and write them into separate chunks of predefined size.

    :param sam_file_1: Path to SAM/BAM file or :class:`~pysam.AlignmentFile`
    :param sam_file_2: Path to SAM/BAM file or :class:`~pysam.AlignmentFile`
    :param output_prefix: prefix str that will form the output files of the form
                          <output_prefix>_<n>.bam
    :param check_sorted: If True, will raise an Exception if SAM/BAM files are not
                         sorted by read name
    :param chunk_size: Number of
    :return: iterator over tuples with path to chunk output file, valid read pairs,
             unmappable read pairs
    """
    logger.debug("Chunk size: {}".format(chunk_size))
    if isinstance(sam_file_1, pysam.AlignmentFile):
        sam1 = sam_file_1
    else:
        sam1 = pysam.AlignmentFile(sam_file_1)

    if isinstance(sam_file_2, pysam.AlignmentFile):
        sam2 = sam_file_2
    else:
        sam2 = pysam.AlignmentFile(sam_file_2)

    sam1_iter = iter(sam1)
    sam2_iter = iter(sam2)

    output_counter = 0
    output_base = output_prefix + "_{}.sam"
    output_file_name = output_base.format(output_counter)
    output_file = open(output_file_name, 'w')
    output_file.write(str(sam1.header))

    get_reads = reads_with_same_qname
    unmappable = 0
    total = 0
    try:
        qname1, reads1, next_read1 = get_reads(sam1_iter)
        qname2, reads2, next_read2 = get_reads(sam2_iter)

        while True:
            check1 = False
            check2 = False
            previous_qname1 = qname1
            previous_qname2 = qname2

            cmp = natural_cmp(qname1, qname2)
            if cmp == 0:  # read name identical
                read_counter = 0
                if len(reads1) + len(reads2) > 1:
                    total += 1
                    for read in reads1 + reads2:
                        output_file.write(read.to_string() + '\n')
                        read_counter += 1

                if read_counter == 0:
                    unmappable += 1

                if total >= chunk_size:
                    output_file.close()
                    yield output_file_name, total, unmappable
                    output_counter += 1
                    output_file_name = output_base.format(output_counter)
                    output_file = open(output_file_name, 'w')
                    output_file.write(str(sam1.header))
                    total = 0
                    unmappable = 0

                qname1, reads1, next_read1 = get_reads(sam1_iter, last_read=next_read1)
                qname2, reads2, next_read2 = get_reads(sam2_iter, last_read=next_read2)
            elif cmp < 0:  # first pointer behind
                qname1, reads1, next_read1 = get_reads(sam1_iter, last_read=next_read1)
                check1 = True
                unmappable += 1
            else:  # second pointer behind
                qname2, reads2, next_read2 = get_reads(sam2_iter, last_read=next_read2)
                check2 = True
                unmappable += 1

            # check that the files are sorted
            if check_sorted:
                if check1 and natural_cmp(previous_qname1, qname1) > 0:
                    raise ValueError("First SAM file is not sorted by "
                                     "read name (samtools sort -n)! Read names:"
                                     "{} and {}".format(previous_qname1, qname1))
                if check2 and natural_cmp(previous_qname2, qname2) > 0:
                    raise ValueError("Second SAM file is not sorted by "
                                     "read name (samtools sort -n)! Read names:"
                                     "{} and {}".format(previous_qname2, qname2))
    except StopIteration:
        pass
    finally:
        output_file.close()
        yield output_file_name, total, unmappable


def write_bed(file_name, regions, mode='w', **kwargs):
    if file_name == '-':
        bed_file = sys.stdout
        must_close = False
    elif hasattr(file_name, 'write'):
        must_close = False
        bed_file = file_name
    else:
        bed_file = open(file_name, mode)
        must_close = True

    try:
        for region in regions:
            bed_file.write(region.as_bed_line(**kwargs) + '\n')
    finally:
        if must_close:
            bed_file.close()
        else:
            bed_file.flush()

    return file_name


def write_gff(file_name, regions, mode='w', **kwargs):
    if file_name == '-':
        gff_file = sys.stdout
        must_close = False
    elif hasattr(file_name, 'write'):
        must_close = False
        gff_file = file_name
    else:
        gff_file = open(file_name, mode)
        must_close = True

    try:
        for region in regions:
            gff_file.write(region.as_gff_line(**kwargs) + '\n')
    finally:
        if must_close:
            gff_file.close()
        else:
            gff_file.flush()

    return file_name


def write_bigwig(file_name, regions, mode='w', score_field='score'):
    logger.debug("Writing output...")
    bw = pyBigWig.open(file_name, mode)

    chromosomes = []
    chromosome_lengths = defaultdict(int)
    interval_chromosomes = []
    interval_starts = []
    interval_ends = []
    interval_values = []
    for region in regions:
        if not isinstance(region.chromosome, str):
            chromosome = region.chromosome.decode() if isinstance(region.chromosome, bytes) \
                else region.chromosome.encode('ascii', 'ignore')
        else:
            chromosome = region.chromosome

        if chromosome not in chromosome_lengths:
            chromosomes.append(chromosome)
        chromosome_lengths[chromosome] = region.end

        interval_chromosomes.append(chromosome)
        interval_starts.append(region.start - 1)
        interval_ends.append(region.end)
        try:
            score = float(getattr(region, score_field))
        except AttributeError:
            score = np.nan
        interval_values.append(score)

    header = []
    for chromosome in chromosomes:
        chromosome_length = chromosome_lengths[chromosome]
        header.append((chromosome, chromosome_length))
    bw.addHeader(header)

    bw.addEntries(interval_chromosomes, interval_starts, ends=interval_ends, values=interval_values)

    bw.close()
    return file_name


def sort_natural_sam(sam_file, output_file=None, sambamba=True, threads=1, _sambamba_path='sambamba'):
    if which(_sambamba_path) is None and sambamba:
        logger.info('Cannot find {} on this machine, falling back to samtools sort. '
                    'This is not a problem, but if you want to speed up your SAM/BAM '
                    'file sorting, install sambamba and ensure it is in your PATH!'.format(_sambamba_path))
        sambamba = False

    basename, extension = os.path.splitext(sam_file)

    replace_input = False
    if output_file is None:
        with tempfile.NamedTemporaryFile(delete=False, prefix='fanc_', suffix=extension) as f:
            output_file = f.name
        replace_input = True

    if sambamba:
        sambamba_command = [_sambamba_path, 'sort', '-N', '-t', str(threads), '-o', output_file, sam_file]
        ret = subprocess.call(sambamba_command)
        if ret != 0:
            sambamba = False
            logger.warning("{} failed, falling back to pysam/samtools".format(_sambamba_path))
    if not sambamba:
        pysam.sort('-n', '-o', output_file, sam_file)

    if replace_input:
        logger.info("Replacing input SAM/BAM file {} with sorted version...".format(sam_file))
        shutil.copy(output_file, sam_file)
        os.remove(output_file)
        output_file = sam_file

    return output_file
