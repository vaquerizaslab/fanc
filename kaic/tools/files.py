'''
Created on May 20, 2015

@author: kkruse1
'''

import tables as t
import os.path
from xml.etree.ElementTree import iterparse, ParseError
import string
import random
import h5py
import pysam
from Bio import SeqIO
import tempfile
import shutil


def without_extension(file_name):
    os.path.splitext(file_name)[0]


def get_extension(file_name):
    os.path.splitext(file_name)[1][1:]


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


def make_dir(dir_name, fail_if_exists=False, make_subdirs=True):
    if make_subdirs:
        f = os.makedirs
    else:
        f = os.mkdir
        
    try: 
        f(dir_name)
    except OSError:
        if not fail_if_exists and not os.path.isdir(dir_name):
            raise


def get_number_of_lines(file_name):
    with open(file_name,'r') as f:
        n = sum(1 for line in f)  # @UnusedVariable
    return n


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


def is_bed_file(file_name):
    file_name = os.path.expanduser(file_name)
    if not os.path.isfile(file_name):
        return False
    
    def is_bed_line(line):
        l = line.rstrip().split("\t")
        if len(l) > 2:
            try:
                int(l[1])
                int(l[2])
            except ValueError:
                return False
        else:
            return False
        return True
        
    with open(file_name, 'r') as f:
        if not is_bed_line(f.readline()):
            return is_bed_line(f.readline())
        else:
            return True


def is_bedpe_file(file_name):
    file_name = os.path.expanduser(file_name)
    if not os.path.isfile(file_name):
        return False
    
    def is_bedpe_line(line):
        l = line.rstrip().split("\t")
        if len(l) > 5:
            try:
                int(l[1])
                int(l[2])
                int(l[4])
                int(l[5])
            except ValueError:
                return False
        else:
            return False
        return True
        
    with open(file_name, 'r') as f:
        if not is_bedpe_line(f.readline()):
            return is_bedpe_line(f.readline())
        else:
            return True
        
        
def is_hic_xml_file(file_name):
    file_name = os.path.expanduser(file_name)
    if not os.path.isfile(file_name):
        return False
    
    try:
        for event, elem in iterparse(file_name):  # @UnusedVariable
            if elem.tag == 'hic':
                return True
            elem.clear()
    except ParseError:
        return False
    
    return False


def is_hdf5_file(file_name):
    file_name = os.path.expanduser(file_name)
    if not os.path.isfile(file_name):
        return False
    
    try:
        f = h5py.File(file_name,'r')
        f.close()
    except IOError:
        return False
    return True


def is_sambam_file(file_name):
    file_name = os.path.expanduser(file_name)
    if not os.path.isfile(file_name):
        return False
    
    try:
        sb = pysam.AlignmentFile(file_name)
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
        fastas = SeqIO.parse(f,'fasta')
        
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


def file_type(file_name):
    file_name = os.path.expanduser(file_name)
    if not os.path.isfile(file_name):
        raise IOError("File {} not found.".format(file_name))

    # might be HDF5
    try:
        f = t.open_file(file_name, mode='r')
        f.close()
        is_hdf5 = True
    except t.HDF5ExtError:
        is_hdf5 = False

    # if is_hdf5:
    #     with t.open_file(file_name, mode='r') as f:
    #         # Hi-C file
    #
    #
    #
    # n = f.get_node('/' + _edge_table_name)
    # if isinstance(n, MaskedTable):
    #     hic_class = Hic
    # elif isinstance(n, t.group.Group):
    #     hic_class = AccessOptimisedHic
    # else:
    #     raise ValueError("%s is not a valid Hi-C object file" % file_name)
    #
    # f.close()
    # return hic_class(file_name=file_name, mode=mode, tmpdir=tmpdir)