'''
Created on May 20, 2015

@author: kkruse1
'''

import tables as t
import os.path

def create_or_open_pytables_file(file_name, inMemory=False, mode='a'):
    
    mem = 0 if inMemory else 1
    
    # check if is existing
    if os.path.isfile(file_name):
        try:
            f = t.open_file(file_name, "r", driver="H5FD_CORE",driver_core_backing_store=mem)
            f.close()
        except t.HDF5ExtError:
            raise ImportError("File exists and is not an HDF5 dict")
        
    return t.open_file(file_name, mode, driver="H5FD_CORE",driver_core_backing_store=mem)


def is_bed_file(file_name):
    if not file_name:
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
    if not file_name:
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