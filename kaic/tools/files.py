import tables as t
import os.path
import string
import random
import pysam
import binascii
from Bio import SeqIO
import tempfile
import shutil


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
    return ''.join(random.SystemRandom().choice(string.uppercase + string.digits) for _ in xrange(length))


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
            fastas.next()
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
