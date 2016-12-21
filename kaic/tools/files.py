import tables as t
import os.path
import string
import random
import pysam
import gzip
import binascii
from Bio import SeqIO
import tempfile
import shutil
from kaic.tools.general import mkdir
from future.utils import string_types
import logging

# configure logging
logger = logging.getLogger(__name__)


def create_temporary_copy(src_file_name, preserve_extension=False):
    """
    Copies the source file into a temporary file.
    Returns a _TemporaryFileWrapper, whose destructor deletes the temp file
    (i.e. the temp file is deleted when the object goes out of scope).
    """
    src_file_name = os.path.expanduser(src_file_name)
    tf_suffix = ''
    if preserve_extension:
        _, tf_suffix = os.path.splitext(src_file_name)
    tf = tempfile.NamedTemporaryFile(suffix=tf_suffix, delete=False)
    shutil.copy2(src_file_name, tf.name)
    return tf.name


def get_number_of_lines(file_name):
    with open(file_name, 'r') as f:
        return sum(1 for _ in f)


def tmp_file_name(tmpdir, prefix='tmp_kaic', extension='h5'):
    name = os.path.join(tmpdir, "{}_{}.{}".format(prefix, binascii.b2a_hex(os.urandom(15)), extension))
    while os.path.exists(name):
        name = os.path.join(tmpdir, "{}_{}.{}".format(prefix, binascii.b2a_hex(os.urandom(15)), extension))
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
        
    is_fasta = True
    with open(file_name, 'r') as f:
        fastas = SeqIO.parse(f, 'fasta')
        
        try:
            next(fastas)
        except StopIteration:
            is_fasta = False
        
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


def split_fastq(fastq_file, output_folder, chunk_size=10000000):
    fastq_file = os.path.expanduser(fastq_file)
    output_folder = mkdir(output_folder)

    basepath, extension = gzip_splitext(fastq_file)
    basename = os.path.basename(basepath)

    fastq_open = fastq_reader(fastq_file)

    split_files = []
    with fastq_open(fastq_file, 'r') as fastq:
        current_chunk_number = 0
        current_chunk_size = 0
        current_file_name = os.path.join(output_folder, '{}_{}{}'.format(current_chunk_number,
                                                                         basename, extension))
        split_files.append(current_file_name)
        current_split_file = fastq_open(current_file_name, 'w')
        for i, line in enumerate(fastq):
            if i % 4 == 0:
                if current_chunk_size >= chunk_size:
                    current_split_file.close()
                    current_chunk_size = 0
                    current_chunk_number += 1
                    current_file_name = os.path.join(output_folder, '{}_{}{}'.format(current_chunk_number,
                                                                                     basename, extension))
                    split_files.append(current_file_name)
                    current_split_file = fastq_open(current_file_name, 'w')
                current_chunk_size += 1
            current_split_file.write(line)
    current_split_file.close()

    return split_files


def split_sam(sam_file, output_folder, chunk_size=1000):
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
                        input_sam = shutil.copy(input_sam, tmp)
                    with pysam.AlignmentFile(input_sam) as r:
                        for alignment in r:
                            o.write(alignment)
                    if tmp:
                        os.remove(input_sam)
    finally:
        if tmp:
            shutil.rmtree(tmp)
    return output_sam
