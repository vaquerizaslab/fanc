'''
Created on May 20, 2015

@author: kkruse1
'''

import tables as t
import pandas as p
import numpy as np
from kaic.tools.files import create_or_open_pytables_file, is_hic_xml_file
from kaic.tools.files import is_bed_file
from kaic.tools.files import is_bedpe_file
import string
import random
from Bio import SeqIO, Restriction, Seq  # @UnusedImport
from kaic.data.general import Table, TableRow, TableArray, TableObject,\
    MetaContainer, Maskable, MaskedTable, FileBased
import os.path

import logging
from kaic.tools.general import ranges
logging.basicConfig(level=logging.INFO)
from xml.etree import ElementTree as et




class GenomicFeature(object):
    def __init__(self, data, names):
        self.data = {}
        self.id2ix = {}
        self.ix2id = {}
        
        if type(data) is dict:
            self.data = data
            j = 0
            for k in data:
                self.id2ix[k] = j
                self.ix2id[j] = k
                
        elif type(data) is list:
            self.data = {}
            for i in xrange(0,len(data)):
                self.data[names[i]] = data[i]
                self.id2ix[names[i]] = i
                self.ix2id[i] = names[i]
                


class BedImproved(Table):
    
    
    
    @staticmethod
    def col_type(name,pos=None):
        col_type = {
            'chrom': (str, t.StringCol(16,pos=pos)), # @UndefinedVariable
            'start': (int, t.Int64Col(pos=pos)), # @UndefinedVariable
            'end': (int, t.Int64Col(pos=pos)), # @UndefinedVariable
            'name': (str, t.StringCol(255,pos=pos)), # @UndefinedVariable
            'score': (float, t.Float32Col(pos=pos)), # @UndefinedVariable
            'strand': (str, t.StringCol(2,pos=pos)), # @UndefinedVariable
            'thickStart': (int, t.Int64Col(pos=pos)), # @UndefinedVariable
            'thickEnd': (int, t.Int64Col(pos=pos)), # @UndefinedVariable
            'itemRgb': (str, t.StringCol(12,pos=pos)), # @UndefinedVariable
            'blockCount': (int, t.Int64Col(pos=pos)), # @UndefinedVariable
            'blockSizes': (str, t.StringCol(255,pos=pos)), # @UndefinedVariable
            'blockStarts': (str, t.StringCol(255,pos=pos)) # @UndefinedVariable
        }
        if name in col_type:
            return col_type[name]
        return (-1, str, t.StringCol(255)) # @UndefinedVariable
    
    @classmethod
    def from_bed_file(cls, file_name, has_header=True, sep="\t", name=None):
        if not is_bed_file:
            raise ImportError("File does not appear to be a BED file")
        
        all_fields = ['chrom','start','end','name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
        
        
        
        with open(file_name, 'r') as f:
            
            # process first line, update table
            line = f.readline()
            fields = line.rstrip().split(sep)
            header = []
            col_types = []
            headerTypes = []
            if has_header:
                header = fields[:]
                line = f.readline()
                fields = line.rstrip().split(sep)
            else:
                header = all_fields[0:len(fields)]
                for i in (len(header), len(fields)):
                    header.append("feature_%d" % i)
                
            for i in xrange(0,len(header)):
                ptype, ttype = BedImproved.col_type(header[i],i+1)
                col_types.append(ttype)
                headerTypes.append(ptype)
                
            
            data = []
            while line != '':
                d = {}
                for i in xrange(0,len(fields)):
                    d[header[i]] = headerTypes[i](fields[i])
                data.append(d)
                    
                line = f.readline()
                fields = line.rstrip().split(sep)
            
            bed = cls()
            
            super(BedImproved, bed).__init__(colnames=header, col_types=col_types, data=data, name=name)
            
            return bed
        
    
    def as_data_frame(self, chrom=None, start=None, end=None):
        query = """"""
        if chrom:
            if query != '':
                query += " & "
            query += "(chrom == '%s')" % chrom
            
        if start is not None:
            if query != '':
                query += " & "
            if type(start) is list:
                query += "(start >= %d) & (start <= %d)" % (start[0], start[1])
            else: 
                query += "(start >= %d)" % (start)
        
        if end:
            if query != '':
                query += " & "
            if type(end) is list:
                query += "(end >= %d) & (end <= %d)" % (end[0], end[1])
            else:
                query += "(end <= %d)" % (end)
        
        
        # get field names
        desc = self._table.description._v_colobjects.copy()
        labels = ['chrom', 'start', 'end']
        if 'name' in desc:
            labels.append('name')
        if 'score' in desc:
            labels.append('score')
        if 'strand' in desc:
            labels.append('strand')
        if 'thickStart' in desc:
            labels.append('thickStart')
        if 'thickEnd' in desc:
            labels.append('thickEnd')
        if 'itemRgb' in desc:
            labels.append('itemRgb')
        if 'blockCount' in desc:
            labels.append('blockCount')
        if 'blockSizes' in desc:
            labels.append('blockSizes')
        if 'blockStarts' in desc:
            labels.append('blockStarts')
        for label in desc:
            if label not in labels and label is not self._rowname_field:
                labels.append(label)
                
        if query != '':
            contacts = [[x[y] for y in labels] for x in self._table.where(query)]
        else:
            contacts = [[x[y] for y in labels] for x in self._table]
            
        df = p.DataFrame(contacts, columns=labels)
        return df
    
    
    
class Bed(object):
    '''
    Bed object for genomic features
    '''


    def __init__(self, file_name=None, name=None):

        inMemory = False
        h5file_name = file_name
        isFlatFile = False
        if file_name == None:
            inMemory=True
        elif is_bed_file(file_name):
            isFlatFile = True
            inMemory=True
        
        if inMemory:
            rs = ''.join(random.SystemRandom().choice(string.uppercase + string.digits) for _ in xrange(6))  # @UndefinedVariable
            h5file_name = rs
            
            
        self.file = create_or_open_pytables_file(h5file_name, inMemory=inMemory)
            
        if not 'bed' in self.file.root:
            columns = {
                'chrom': t.StringCol(16), # @UndefinedVariable
                'start': t.Int64Col(), # @UndefinedVariable
                'end': t.Int64Col() # @UndefinedVariable
            }
            self.table = self.file.create_table("/", 'bed', columns)
        else:
            self.table = self.file.root.bed
        
        self.name = name if name else file_name
        
        if isFlatFile:
            self.load_bed_file(file_name)
    
    def close(self):
        self.file.close()
    
    def __del__(self):
        try:
            print "Closing hdf5 file"
            self.close()
        except AttributeError:
            print "Nothing to close"
    
    
    def load_bed_file(self,in_file,has_header=True):
        
        if not is_bed_file:
            raise ImportError("File does not appear to be a BED file")
        
        with open(in_file, 'r') as f:
            
            # process first line, update table
            line = f.readline()
            fields = line.rstrip().split("\t")
            
            desc = self.table.description._v_colobjects.copy()
            header = []
            headerTypes = []
            if has_header:
                for i in xrange(0,len(fields)):
                    #if fields[i] in desc:
                    #    raise ValueError("Duplicate column name! " + fields[i])
                    if fields[i] == 'name' or fields[i] == 'blockSizes' or fields[i] == 'blockStarts':
                        desc[fields[i]] = t.StringCol(255) # @UndefinedVariable
                        headerTypes.append(str)
                    elif fields[i] == 'strand':
                        desc[fields[i]] = t.StringCol(1) # @UndefinedVariable
                        headerTypes.append(str)
                    elif fields[i] == 'itemRgb':
                        desc[fields[i]] = t.StringCol(12) # @UndefinedVariable
                        headerTypes.append(str)
                    elif fields[i] == 'chrom':
                        desc[fields[i]] = t.StringCol(16) # @UndefinedVariable
                        headerTypes.append(str)
                    elif fields[i] == 'score':
                        desc[fields[i]] = t.Float32Col() # @UndefinedVariable
                        headerTypes.append(float)
                    elif (fields[i] == 'thickStart' or fields[i] == 'thickEnd' or 
                        fields[i] == 'blockCount' or 
                        fields[i] == 'start' or fields[i] == 'end'):
                        desc[fields[i]] = t.Int64Col() # @UndefinedVariable
                        headerTypes.append(int)
                    else:
                        desc[fields[i]] = t.StringCol(255) # @UndefinedVariable
                        headerTypes.append(str)
                    header.append(fields[i])
                line = f.readline()
                fields = line.rstrip().split("\t")
            else:
                
                for i in xrange(0,len(fields)):
                    if i == 0:
                        header.append('chrom')
                        headerTypes.append(str)
                        desc['chrom'] = t.StringCol(16) # @UndefinedVariable
                    elif i == 1:
                        header.append('start')
                        headerTypes.append(str)
                        desc['start'] = t.Int64Col() # @UndefinedVariable
                    elif i == 2:
                        header.append('end')
                        headerTypes.append(str)
                        desc['end'] = t.Int64Col(255) # @UndefinedVariable
                    elif i == 3:
                        header.append('name')
                        headerTypes.append(str)
                        desc['name'] = t.StringCol(255) # @UndefinedVariable
                    elif i == 4:
                        header.append('score')
                        headerTypes.append(float)
                        desc['score'] = t.Float32Col() # @UndefinedVariable
                    elif i == 5:
                        header.append('strand')
                        headerTypes.append(str)
                        desc['strand'] = t.StringCol(1) # @UndefinedVariable
                    elif i == 6:
                        header.append('thickStart')
                        headerTypes.append(int)
                        desc['thickStart'] = t.Int64Col() # @UndefinedVariable
                    elif i == 7:
                        header.append('thickEnd')
                        headerTypes.append(int)
                        desc['thickEnd'] = t.Int64Col() # @UndefinedVariable
                    elif i == 8:
                        header.append('itemRgb')
                        headerTypes.append(str)
                        desc['itemRgb'] = t.StringCol(12) # @UndefinedVariable
                    elif i == 9:
                        header.append('blockCount')
                        headerTypes.append(int)
                        desc['blockCount'] = t.Int64Col() # @UndefinedVariable
                    elif i == 10:
                        header.append('blockSizes')
                        headerTypes.append(str)
                        desc['blockSizes'] = t.StringCol(255) # @UndefinedVariable
                    elif i == 11:
                        header.append('blockStarts')
                        headerTypes.append(str)
                        desc['blockStarts'] = t.StringCol(255) # @UndefinedVariable
                    else:
                        header.append('feature_' + str(i))
                        headerTypes.append(str)
                        desc['feature_' + str(i)] = t.StringCol(255) # @UndefinedVariable
            
            table2 = self.file.create_table(self.file.root, 'table2', desc, "bed", t.Filters(1))
 
            # Copy the user attributes
            self.table.attrs._f_copy(table2)
             
            # Fill the rows of new table with default values
            for i in xrange(self.table.nrows):
                table2.row.append()
            # Flush the rows to disk
            table2.flush()
            
            # Copy the columns of source table to destination
            for col in self.table.description._v_colobjects:
                if (len(getattr(self.table.cols, col)[:]) > 0 and
                    len(getattr(table2.cols, col)[:]) > 0):
                    getattr(table2.cols, col)[:] = getattr(self.table.cols, col)[:]
            
            # fill with new data
            entry = table2.row
            while line != '':
                
                #if len(fields) == len(headerTypes):
                for i in xrange(0,len(fields)):
                    value = headerTypes[i](fields[i])
                    entry[header[i]] = value
                entry.append()
                
                line = f.readline()
                fields = line.rstrip().split("\t")
            table2.flush()
            
            # Remove the original table
            self.table.remove()
             
            # Move table2 to table
            table2.move('/','bed')
            self.table = table2
    
    
    
    def as_data_frame(self, chrom=None, start=None, end=None, as_gene=False):
        query = """"""
        if chrom:
            if query != '':
                query += " & "
            query += "(chrom == '%s')" % chrom
            
        if start is not None:
            if query != '':
                query += " & "
            if type(start) is list:
                query += "(start >= %d) & (start <= %d)" % (start[0], start[1])
            else: 
                query += "(start >= %d)" % (start)
        
        if end:
            if query != '':
                query += " & "
            if type(end) is list:
                query += "(end >= %d) & (end <= %d)" % (end[0], end[1])
            else:
                query += "(end <= %d)" % (end)
        
        
        # get field names
        desc = self.table.description._v_colobjects.copy()
        labels = ['chrom', 'start', 'end']
        if 'gene' in desc:
            labels.append('gene')
        if 'score' in desc:
            labels.append('score')
        if 'strand' in desc:
            labels.append('strand')
        if 'type' in desc:
            labels.append('type')
        if 'thickStart' in desc:
            labels.append('thickStart')
        if 'thickEnd' in desc:
            labels.append('thickEnd')
        if 'itemRgb' in desc:
            labels.append('itemRgb')
        if 'blockCount' in desc:
            labels.append('blockCount')
        if 'blockSizes' in desc:
            labels.append('blockSizes')
        if 'blockStarts' in desc:
            labels.append('blockStarts')
        for label in desc:
            if label not in labels:
                labels.append(label)
                
        if query != '':
            contacts = [[x[y] for y in labels] for x in self.table.where(query)]
        else:
            contacts = [[x[y] for y in labels] for x in self.table]
            
        df = p.DataFrame(contacts, columns=labels)
        return df
    
    
    
    
    
    
class Bedpe(object):
    '''
    Bedpe object for genomic features
    '''


    def __init__(self, file_name=None, name=None):
        
        inMemory = False
        h5file_name = file_name
        isFlatFile = False
        if file_name == None:
            inMemory=True
        elif is_bed_file(file_name):
            isFlatFile = True
            inMemory=True
        
        if inMemory:
            rs = ''.join(random.SystemRandom().choice(string.uppercase + string.digits) for _ in xrange(6))  # @UndefinedVariable
            h5file_name = rs
            
            
        self.file = create_or_open_pytables_file(h5file_name, inMemory=inMemory)
        
        if not 'bedpe' in self.file.root:
            columns = {
                'chrom1': t.StringCol(16), # @UndefinedVariable
                'start1': t.Int64Col(), # @UndefinedVariable
                'end1': t.Int64Col(), # @UndefinedVariable
                'chrom2': t.StringCol(16), # @UndefinedVariable
                'start2': t.Int64Col(), # @UndefinedVariable
                'end2': t.Int64Col() # @UndefinedVariable
            }
            self.table = self.file.create_table("/", 'bedpe', columns)
        else:
            self.table = self.file.root.bedpe
        
        self.name = name if name else file_name
        
        if isFlatFile:
            self.load_bedpe_file(file_name)
    
    def close(self):
        self.file.close()
    
    def __del__(self):
        try:
            print "Closing hdf5 file"
            self.close()
        except AttributeError:
            print "Nothing to close"
    
    
    def load_bedpe_file(self,in_file,has_header=True):
        
        if not is_bedpe_file:
            raise ImportError("File does not appear to be a BED file")
        
        with open(in_file, 'r') as f:
            
            # process first line, update table
            line = f.readline()
            fields = line.rstrip().split("\t")
            
            desc = self.table.description._v_colObjects.copy()
            header = []
            headerTypes = []
            if has_header:
                for i in xrange(0,len(fields)):
                    #if fields[i] in desc:
                    #    raise ValueError("Duplicate column name! " + fields[i])
                    if fields[i] == 'name':
                        desc[fields[i]] = t.StringCol(255) # @UndefinedVariable
                        headerTypes.append(str)
                    elif fields[i] == 'strand1' or fields[i] == 'strand2':
                        desc[fields[i]] = t.StringCol(1) # @UndefinedVariable
                        headerTypes.append(str)
                    elif fields[i] == 'chrom1' or fields[i] == 'chrom2':
                        desc[fields[i]] = t.StringCol(16) # @UndefinedVariable
                        headerTypes.append(str)
                    elif fields[i] == 'score':
                        desc[fields[i]] = t.Float32Col() # @UndefinedVariable
                        headerTypes.append(float)
                    elif (fields[i] == 'start1' or fields[i] == 'start2' or 
                        fields[i] == 'end1' or fields[i] == 'end2'):
                        desc[fields[i]] = t.Int64Col() # @UndefinedVariable
                        headerTypes.append(int)
                    else:
                        desc[fields[i]] = t.StringCol(255) # @UndefinedVariable
                        headerTypes.append(str)
                    header.append(fields[i])
                line = f.readline()
                fields = line.rstrip().split("\t")
            else:
                
                for i in xrange(0,len(fields)):
                    if i == 0:
                        header.append('chrom1')
                        headerTypes.append(str)
                        desc['chrom1'] = t.StringCol(16) # @UndefinedVariable
                    elif i == 1:
                        header.append('start1')
                        headerTypes.append(str)
                        desc['start1'] = t.Int64Col() # @UndefinedVariable
                    elif i == 2:
                        header.append('end1')
                        headerTypes.append(str)
                        desc['end1'] = t.Int64Col() # @UndefinedVariable
                    elif i == 3:
                        header.append('chrom2')
                        headerTypes.append(str)
                        desc['chrom2'] = t.StringCol(16) # @UndefinedVariable
                    elif i == 4:
                        header.append('start2')
                        headerTypes.append(int)
                        desc['start2'] = t.Int64Col() # @UndefinedVariable
                    elif i == 5:
                        header.append('end2')
                        headerTypes.append(int)
                        desc['end2'] = t.Int64Col() # @UndefinedVariable
                    elif i == 6:
                        header.append('name')
                        headerTypes.append(str)
                        desc['name'] = t.StringCol(255) # @UndefinedVariable
                    elif i == 7:
                        header.append('score')
                        headerTypes.append(float)
                        desc['score'] = t.Float32Col() # @UndefinedVariable
                    elif i == 8:
                        header.append('strand1')
                        headerTypes.append(str)
                        desc['strand1'] = t.StringCol(1) # @UndefinedVariable
                    elif i == 9:
                        header.append('strand2')
                        headerTypes.append(int)
                        desc['strand2'] = t.StringCol(1) # @UndefinedVariable
                    elif i == 10:
                        header.append('blockSizes')
                        headerTypes.append(str)
                        desc['blockSizes'] = t.StringCol(255) # @UndefinedVariable
                    else:
                        header.append('feature_' + str(i))
                        headerTypes.append(str)
                        desc['feature_' + str(i)] = t.StringCol(255) # @UndefinedVariable
            
            if not 'chrom1' in header:
                raise ImportError("File must contain chrom1 field!")
            if not 'chrom2' in header:
                raise ImportError("File must contain chrom2 field!")
            if not 'start1' in header:
                raise ImportError("File must contain start1 field!")
            if not 'start2' in header:
                raise ImportError("File must contain start2 field!")
            if not 'end1' in header:
                raise ImportError("File must contain end1 field!")
            if not 'end2' in header:
                raise ImportError("File must contain end2 field!")
            
            table2 = self.file.createTable(self.file.root, 'table2', desc, "bedpe", t.Filters(1))
 
            # Copy the user attributes
            self.table.attrs._f_copy(table2)
             
            # Fill the rows of new table with default values
            for i in xrange(self.table.nrows):
                table2.row.append()
            # Flush the rows to disk
            table2.flush()
            
            # Copy the columns of source table to destination
            for col in self.table.description._v_colObjects:
                if (len(getattr(self.table.cols, col)[:]) > 0 and
                    len(getattr(table2.cols, col)[:]) > 0):
                    getattr(table2.cols, col)[:] = getattr(self.table.cols, col)[:]
             
            # fill with new data
            entry = table2.row
            while line != '':
                
                if len(fields) == len(headerTypes):
                    for i in xrange(0,len(fields)):
                        value = headerTypes[i](fields[i])
                        entry[header[i]] = value
                    entry.append()
                
                line = f.readline()
                fields = line.rstrip().split("\t")
            table2.flush()
            
            # Remove the original table
            self.table.remove()
             
            # Move table2 to table
            table2.move('/','bedpe')
            self.table = table2
    
    def as_data_frame(self, chrom=None, lower_bound=None, upper_bound=None):
        query = """"""
        if chrom:
            if query != '':
                query += " & "
            query += "(chrom1 == '%s')" % chrom
            
        if lower_bound:
            if query != '':
                query += " & "
            query += "(end1 >= %d) & (end2 >= %d)" % (lower_bound,lower_bound)
        
        if upper_bound:
            if query != '':
                query += " & "
            query += "(start1 <= %d) & (start2 <= %d)" % (upper_bound,upper_bound)
        
        
        # get field names
        desc = self.table.description._v_colObjects.copy()
        labels = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']
        if 'name' in desc:
            labels.append('name')
        if 'score' in desc:
            labels.append('score')
        if 'strand1' in desc:
            labels.append('strand1')
        if 'strand2' in desc:
            labels.append('strand2')
        if 'blockSizes' in desc:
            labels.append('blockSizes')
        for label in desc:
            if label not in labels:
                labels.append(label)
        
        print labels
        
        print "Running query"
        if query != '':
            contacts = [[x[y] for y in labels] for x in self.table.where(query)]
        else:
            contacts = [[x[y] for y in labels] for x in self.table]
            
        df = p.DataFrame(contacts, columns=labels)
        return df
        


class Hic(Bedpe):
    def __init__(self, file_name=None, name=None):
        Bedpe.__init__(self, file_name=file_name, name=name)
        

        
    def as_data_frame(self, resolution, chrom=None, lower_bound=None, upper_bound=None):
        
        query = """"""
        if chrom:
            if query != '':
                query += " & "
            query += "(chrom1 == '%s')" % chrom
            
        if lower_bound:
            if query != '':
                query += " & "
            query += "(end1 >= %d) & (end2 >= %d)" % (lower_bound,lower_bound)
        
        if upper_bound:
            if query != '':
                query += " & "
            query += "(start1 <= %d) & (start2 <= %d)" % (upper_bound,upper_bound)
        
        print "Running query"
        if query != '':
            contacts = [[x['start1'],x['start2'],x['score']] for x in self.table.where(query)]
        else:
            contacts = [[x['start1'],x['start2'],x['score']] for x in self.table]
        
        
        if lower_bound is None:
            lower_bound = 0
        if upper_bound is None:
            upper_bound = max(max(contacts)[0], max(contacts[1])) + 1
            
            
        # TODO
        # take into account chromosome sizes
        # pull resolution from Hi-C object
        min_ix = int(lower_bound/resolution)*resolution
        max_ix = int(upper_bound/resolution)*resolution+resolution
#         print "Calculating lowest bound"
#         min_ix = min(min(contacts)[0:2])
#         print min_ix
#         print "Calculating upper bound"
#         max_ix = max(max(contacts)[0:2])
#         print max_ix
        
        
        labels = range(min_ix,max_ix+resolution,resolution)
        ix_l = int(min_ix/resolution)
        
        print "Assigning to matrix"
        M = np.zeros((len(labels),len(labels)))
        for c in contacts:
            i = int(c[0]/resolution)-ix_l
            j = int(c[1]/resolution)-ix_l
            try:
                M[i,j] = c[2]
                M[j,i] = c[2]
            except IndexError:
                raise IndexError("%d - %d (%d, %d)" % (c[0], c[1], i, j))
                
        
        print "Creating data frame"
        df = p.DataFrame(M, index=labels, columns=labels)
        
        return df
        
        
    def as_matrix(self, resolution, chrom=None, lower_bound=None, upper_bound=None):
        query = """"""
        if chrom:
            if query != '':
                query += " & "
            query += "(chrom1 == %s)" % chrom
            
        if lower_bound:
            if query != '':
                query += " & "
            query += "(end1 >= %d) & (end2 >= %d)" % (lower_bound,lower_bound)
        
        if upper_bound:
            if query != '':
                query += " & "
            query += "(start1 <= %d) & (start2 <= %d)" % (upper_bound,upper_bound)
        
        print "Running query"
        if query != '':
            contacts = [[x['start1'],x['start2'],x['score']] for x in self.table.where(query)]
        else:
            contacts = [[x['start1'],x['start2'],x['score']] for x in self.table]
        
        if lower_bound is None:
            lower_bound = 0
        if upper_bound is None:
            upper_bound = max(max(contacts)[0], max(contacts[1])) + 1
        
        min_ix = int(lower_bound/resolution)*resolution
        max_ix = int(upper_bound/resolution)*resolution+resolution
        
        labels = range(min_ix,max_ix+resolution,resolution)
        ix_l = int(min_ix/resolution)
        
        print "Assigning to matrix"
        M = np.zeros((len(labels),len(labels)))
        for c in contacts:
            i = int(c[0]/resolution)-ix_l
            j = int(c[1]/resolution)-ix_l
            try:
                M[i,j] = c[2]
                M[j,i] = c[2]
            except IndexError:
                raise IndexError("%d - %d (%d, %d)" % (c[0], c[1], i, j))
        
        return M
    
    def directionality(self, resolution, window_size=2000000):
        M = self.as_matrix(resolution)
        bin_window_size = int(window_size/resolution)
        if window_size%resolution > 0:
            bin_window_size += 1
        
        n_bins = M.shape[0]
        dis = np.zeros(n_bins,dtype='float64')
        for i in xrange(0,n_bins):
            max_window_size = min(bin_window_size, n_bins-i, i)
            start = i-max_window_size
            end = i+max_window_size
            
            A = np.sum(M[i][start:i])
            B = np.sum(M[i][i+1:end])
            E = (A+B)/2
            
            if B == A:
                dis[i] = 0
            else:
                dis[i] = ((B-A)/abs(B-A)) * ((((A-E)**2)/E) + (((B-E)**2)/E))
        
        return dis
    
    


class Chromosome(object):
    def __init__(self, name=None, length=None, sequence=None):
        self.name = name
        self.length = length
        self.sequence = sequence
        if length is None and sequence is not None:
            self.length = len(sequence)
        if sequence is None and length is not None:
            self.length = length
            
    def __repr__(self):
        return "Name: %s\nLength: %d\nSequence: %s" % (self.name if self.name else '',
                                                       self.length if self.length else -1,
                                                       self.sequence[:20] + "..." if self.sequence else '')
        
    def __len__(self):
        return self.length
    
    def __getitem__(self, key):
        if key == 'name':
            return self.name
        if key == 'length':
            return self.length
        if key == 'sequence':
            return self.sequence
        
    
    @classmethod
    def from_fasta(cls, file_name, name=None, include_sequence=True):
        if type(file_name) is file:
            fastas = SeqIO.parse(file_name,'fasta')
        else:
            fastas = SeqIO.parse(open(file_name,'r'),'fasta')
            
        try:
            fasta = fastas.next()
        except StopIteration:
            raise ValueError("File %s does not appear to be a FASTA file" % file_name)
        
        if include_sequence:
            return cls(name if name else fasta.id, length=len(fasta), sequence=str(fasta.seq))
        else:
            return cls(name if name else fasta.id, length=len(fasta))
            
        
        
    def get_restriction_sites(self, restriction_enzyme):
        logging.info("Calculating RE sites")
        try:
            re = eval('Restriction.%s' % restriction_enzyme)
        except SyntaxError:
            raise ValueError("restriction_enzyme must be a string")
        except AttributeError:
            raise ValueError("restriction_enzyme string is not recognized: %s" % restriction_enzyme)
        
        return re.search(Seq.Seq(self.sequence))
    
    
        
    
    
        
class Genome(Table):
    def __init__(self, file_name=None, chromosomes=None):        
        self.file = create_or_open_pytables_file(file_name)
            
        columns = ["ix", "name", "length"]
        column_types = [t.Int32Col(pos=0), t.StringCol(50, pos=1), t.Int32Col(pos=2)]  # @UndefinedVariable
        Table.__init__(self, colnames=columns, col_types=column_types)
        
        try:
            self._sequences = self.file.get_node('/genome_sequences')
        except t.NoSuchNodeError:
            self._sequences = self.file.create_vlarray("/", 'genome_sequences', t.VLStringAtom())
        
        if chromosomes is not None:
            for chromosome in chromosomes:
                self.add_chromosome(chromosome)
            
    @classmethod
    def from_folder(cls, folder_name, file_name=None, exclude=None, include_sequence=True):
        chromosomes = []
        folder_name = os.path.expanduser(folder_name)
        for f in os.listdir(folder_name):
            try:
                chromosome = Chromosome.from_fasta(folder_name + "/" + f, include_sequence=include_sequence)
                logging.info("Adding chromosome %s" % chromosome.name)
                if exclude is None:
                    chromosomes.append(chromosome)
                elif chromosome.name not in exclude:
                    chromosomes.append(chromosome)
            except (ValueError, IOError):
                pass
        
        return cls(chromosomes=chromosomes, file_name=file_name)
        
    
    def __getitem__(self, key):
        res = Table.__getitem__(self, key)
        
        if isinstance(res, TableRow):
            return Chromosome(name=res.name, length=res.length, sequence=self._sequences[res.ix])
        elif isinstance(res, TableArray):
            l = []
            for row in res:
                l.append(Chromosome(name=row["name"], length=row["length"], sequence=self._sequences[row["ix"]]))
            return l
        return res
    
    def __iter__(self):
        this = self
        class Iter:
            def __init__(self):
                self.current = 0
                
            def __iter__(self):
                self.current = 0
                return self
            
            def next(self):
                if self.current >= len(this):
                    raise StopIteration
                self.current += 1
                return this[self.current-1]
        return Iter()
    
    def __del__(self):
        self.file.close()
        super(Genome, self).__del__()
    
    
    def add_chromosome(self, chromosome):
        i = len(self)-1
        
        n = str(i)
        if chromosome.name is not None:
            n = chromosome.name
        i += 1
        
        l = 0
        if chromosome.length is not None:
            l = chromosome.length
        
        s = ''
        if chromosome.sequence is not None:
            s = chromosome.sequence
            if l == 0:
                l = len(s)
        
        self.append([i,n,l], rownames=[n])
        self._sequences.append(s)
        self._sequences.flush()

    
    def get_regions(self, split, file_name=None):
        regions = GenomicRegions(file_name=file_name)
        for chromosome in self:
            split_locations = []
            if isinstance(split,str):
                split_locations = chromosome.get_restriction_sites(split)
            elif isinstance(split,int):
                for i in xrange(split,len(chromosome)-1, split):
                    split_locations.append(i)
            else:
                for i in split:
                    split_locations.append(i)
            
            for i in xrange(0,len(split_locations)):
                if i == 0:
                    region = GenomicRegion(start=1, end=split_locations[i], chromosome=chromosome.name)
                else:
                    region = GenomicRegion(start=split_locations[i-1]+1, end=split_locations[i], chromosome=chromosome.name)
                
                regions.add_region(region, flush=False)
                
            # add last node
            if len(split_locations) > 0:
                region = GenomicRegion(start=split_locations[len(split_locations)-1]+1, end=chromosome.length, chromosome=chromosome.name)
            else:
                region = GenomicRegion(start=1, end=chromosome.length, chromosome=chromosome.name)
            regions.add_region(region, flush=False)
            regions._flush()
                
        return regions

class GenomicRegion(TableObject):
    def __init__(self, start, end, chromosome=None, strand=None, ix=None):
        self.start = start
        self.end = end
        self.strand = strand
        self.chromosome = chromosome
        self.ix = ix
    
    @classmethod
    def from_row(cls, row):
        strand = row['strand']
        if strand == 0:
            strand = None
        return cls(start=row["start"], end=row["end"],
                   strand=strand, chromosome=row["chromosome"])
    
    @classmethod
    def from_string(cls, region_string):
        chromosome= None
        start = None
        end = None
        strand = None
        
        # strip whitespace
        no_space_region_string = "".join(region_string.split())
        fields = no_space_region_string.split(':')
        
        if len(fields) > 3:
            raise ValueError("Genomic range string must be of the form <chromosome>[:<start>-<end>:[<strand>]]")
        
        # there is chromosome information
        if len(fields) > 0:
            chromosome = fields[0]
        
        # there is range information
        if len(fields) > 1 and fields[1] != '':
            start_end_bp = fields[1].split('-')
            if len(start_end_bp) > 0:
                try:
                    start = int(start_end_bp[0])
                except ValueError:
                    raise ValueError("Start of genomic range must be integer")
            
            if len(start_end_bp) > 1:
                try:
                    end = int(start_end_bp[1])
                except ValueError:
                    raise ValueError("End of genomic range must be integer") 
        
        # there is strand information
        if len(fields) > 2:
            if fields[2] == '+' or fields[2] == '+1' or fields[2] == '1':
                strand = 1
            elif fields[2] == '-' or fields[2] == '-1':
                strand = -1
            else:
                raise ValueError("Strand only can be one of '+', '-', '+1', '-1', and '1'")
        return cls(start=start, end=end, chromosome=chromosome, strand=strand)
    
    def to_string(self):
        region_string = ''
        if self.chromosome is not None:
            region_string += '%s' % self.chromosome
            
            if self.start is not None:
                region_string += ':%d' % self.start
                
                if self.end is not None:
                    region_string += '-%d' % self.end
                
                if self.strand is not None:
                    if self.strand == 1:
                        region_string += ':+'
                    else:
                        region_string += ':-'
        return region_string
    
    def __repr__(self):
        return self.to_string()
        
class GenomicRegions(Table):
    class GenomicRegionDescription(t.IsDescription):
        chromosome = t.StringCol(50, pos=0)
        start = t.Int64Col(pos=1)
        end = t.Int64Col(pos=2)
        strand = t.Int8Col(pos=3)
        
    def __init__(self, file_name=None, regions=None):
        if not isinstance(file_name, str) and not isinstance(file_name, t.file.File):
            regions = file_name
            file_name = None
        self.file = create_or_open_pytables_file(file_name)
        
        # create table
        columns = ["chromosome", "start", "end", "strand"]
        column_types = [t.StringCol(50, pos=0), t.Int64Col(pos=1),
                        t.Int64Col(pos=2), t.Int8Col(pos=3)]
        Table.__init__(self, colnames=columns, col_types=column_types, return_type=GenomicRegion)
        
        # load data if provided
        if regions is not None:
            for region in regions:
                self.add_region(region, flush=False)
        self._flush()
                
    def add_region(self, region, flush=True):
        
        # try access by attribute first
        try:
            chromosome = region.chromosome
            start = region.start
            end = region.end
            strand = region.strand
        # if that fails try access by item
        except AttributeError:
            chromosome = region['chromosome']
            start = region['start']
            end = region['end']
            strand = region['strand']
        
        if strand is None:
            strand = 0
        
        self._append_row_dict({
            'chromosome': chromosome,
            'start': start,
            'end': end,
            'strand': strand
        }, flush=flush)
    

class RegionsTable(FileBased):
    class RegionDescription(t.IsDescription):
        ix = t.Int32Col(pos=0)
        chromosome = t.StringCol(50,pos=1)
        start = t.Int64Col(pos=2)
        end = t.Int64Col(pos=3)
    
    def __init__(self, data=None, file_name=None,
                       table_name_regions='regions'):
        
        # parse potential unnamed argument
        if data is not None:
            # data is file name
            if type(data) is str:                
                if not os.path.isfile(data) and file_name is None:
                    file_name = data
                    data = None
        
        if file_name is not None:
            file_name = os.path.expanduser(file_name)
        
        FileBased.__init__(self, file_name)
        
        # check if this is an existing Hi-C file
        if table_name_regions in self.file.root:
            self._regions = self.file.get_node('/', table_name_regions)
            self._max_region_ix = max(row['ix'] for row in self._regions.iterrows())
        else:
            self._regions = t.Table(self.file.root, table_name_regions,
                                    RegionsTable.RegionDescription, expectedrows=10000)
            self._max_region_ix = -1
        
        if data is not None:
            self.add_regions(data)
        
    def add_region(self, region, flush=True):
        ix = -1
        
        if isinstance(region, GenomicRegion):
            if hasattr(region, 'ix') and region.ix is not None:
                ix = region.ix
            chromosome = region.chromosome
            start = region.start
            end = region.end
        elif type(region) is dict:
            if 'ix' in region:
                ix = region['ix']
            chromosome = region['chromosome']
            start = region['start']
            end = region['end']
        else:
            try:
                offset = 0
                if len(region) == 4:
                    ix = region[0]
                    offset += 1
                chromosome = region[offset]
                start = region[offset + 1]
                end = region[offset + 2]
            except TypeError:
                raise ValueError("Node parameter has to be HicNode, dict, or list")
        
        if ix == -1:
            ix = self._max_region_ix + 1
        
        # actually append
        row = self._regions.row
        row['ix'] = ix
        row['chromosome'] = chromosome
        row['start'] = start
        row['end'] = end
        row.append()
        
        if ix > self._max_region_ix:
            self._max_region_ix = ix
            
        if flush:
            self._regions.flush()
        
        return ix
    
    def add_regions(self, regions):
        for region in regions:
            self.add_region(region, flush=False)
        self._regions.flush()
    
    def _get_region_ix(self, region):
        condition = "(start == %d) & (end == %d) & (chromosome == '%s')"
        condition = condition % (region.start, region.end, region.chromosome)
        for res in self._regions.where(condition):
            return res["ix"]
        return None
            
        
    def regions(self):
        this = self
        class RegionIter:
            def __init__(self):
                self.iter = iter(this._regions)
                
            def __iter__(self):
                return self
            
            def next(self):
                return HicNode.from_row(self.iter.next())
            
            def __len__(self):
                return len(this._regions)
            
        return RegionIter()

class HicNode(GenomicRegion, TableObject):
    def __init__(self, chromosome=None, start=None, end=None, ix=None):
        self.ix = ix
        super(HicNode, self).__init__(chromosome=chromosome, start=start, end=end, ix=ix)
    
    def __repr__(self):
        if self.ix is None:
            return "%s, %d-%d" % (self.chromosome, self.start, self.end)
        else:
            return "%d: %s, %d-%d" % (self.ix, self.chromosome, self.start, self.end)
    
    @classmethod
    def from_string(cls, region_string):
        node = super(HicNode, cls).from_string(region_string)
        node.ix = None
        return node
    
    @classmethod
    def from_row(cls, row):
        return cls(chromosome=row['chromosome'], start=row['start'], end=row['end'], ix=row['ix'])
    
        
class HicEdge(TableObject):
    def __init__(self, source, sink, weight=1):
        self.source = source
        self.sink = sink
        self.weight = weight
    
    def __repr__(self):
        return "%d--%d (%.2f)" % (self.source, self.sink, self.weight)
    
    @classmethod
    def from_row(cls, row):
        return cls(source=row['source'], sink=row['sink'], weight=row['weight'])


class HicBasic(Maskable, MetaContainer, RegionsTable, FileBased):

    class HicEdgeDescription(t.IsDescription):
        source = t.Int32Col(pos=0)  
        sink = t.Int32Col(pos=1)  
        weight = t.Float64Col(pos=2)  
    
    def __init__(self, data=None, file_name=None,
                       table_name_nodes='nodes',
                       table_name_edges='edges'):
        
        # private variables
        self._max_node_ix = -1
        
        # parse potential unnamed argument
        if data is not None:
            # data is file name
            if type(data) is str:
                data = os.path.expanduser(data)
                
                if (not os.path.isfile(data) or not is_hic_xml_file(data)) and file_name is None:
                    file_name = data
                    data = None
        
        if file_name is not None:
            file_name = os.path.expanduser(file_name)
        
        FileBased.__init__(self, file_name)
        RegionsTable.__init__(self, file_name=file_name, table_name_regions=table_name_nodes)

        if table_name_edges in self.file.root:
            self._edges = self.file.get_node('/', table_name_edges)
        else:
            self._edges = MaskedTable(self.file.root, table_name_edges,
                                      HicBasic.HicEdgeDescription, expectedrows=500000)
        
        self._edges.flush()
        
        # generate tables from inherited classes
        Maskable.__init__(self, self.file)
        MetaContainer.__init__(self, self.file)
        
        # index node table
        try:
            self._regions.cols.ix.create_csindex()
        except ValueError:
            # Index exists, no problem!
            pass
        try:
            self._regions.cols.start.create_csindex()
        except ValueError:
            # Index exists, no problem!
            pass
        try:
            self._regions.cols.end.create_csindex()
        except ValueError:
            # Index exists, no problem!
            pass
        # index edge table
        try:
            self._edges.cols.source.create_csindex()
        except ValueError:
            # Index exists, no problem!
            pass
        try:
            self._edges.cols.sink.create_csindex()
        except ValueError:
            # Index exists, no problem!
            pass
        
        # add data
        if data is not None:
            if type(data) is str:
                if is_hic_xml_file(data):
                    xml = HicXmlFile(data)
                    for node in xml.nodes():
                        self.add_node(node, flush=False)
                    self.flush()
                    
                    for edge in xml.edges():
                        self.add_edge(edge, flush=False)
                    self.flush()
                else:
                    raise ValueError("File is not in Hi-C XML format")
                    
            # data is existing HicBasic object
            elif isinstance(data, HicBasic):
                for node in data.nodes():
                    self.add_node(node, flush=False)
                self.flush()
                for edge in data.edges():
                    self.add_edge(edge,flush=False)
                self.flush()
            else:
                raise ValueError("Input data type not recognized")
        
    
    
    def __del__(self):
        self.close()
        
    @classmethod
    def from_hiclib(cls, hl, file_name=None):
        hic = cls(file_name=file_name)
        
        # nodes
        chrms = {hl.genome.chrmStartsBinCont[i] : hl.genome.chrmLabels[i] for i in xrange(0,len(hl.genome.chrmLabels))}
        chromosome = ''
        for i in xrange(0,len(hl.genome.posBinCont)):
            start = hl.genome.posBinCont[i]+1
            if i in chrms:
                chromosome = chrms[i]
            
            if i < len(hl.genome.posBinCont)-1:
                end = hl.genome.posBinCont[i+1]
            else:
                ix = hl.genome.label2idx[chromosome]
                end = hl.genome.chrmLens[ix]
            
            hic.add_node([chromosome, start, end], flush=False)
        hic.flush(flush_edges=False)
        
        # edges
        for chr1, chr2 in hl.data:
            data = hl.data[(chr1, chr2)].getData()
            chr1StartBin = hl.genome.chrmStartsBinCont[chr1]
            chr2StartBin = hl.genome.chrmStartsBinCont[chr2]
            
            for i in xrange(0,data.shape[0]):
                iNode = i+chr1StartBin
                start = i
                if chr1 != chr2:
                    start = 0
                for j in xrange(start,data.shape[1]):
                    jNode = j+chr2StartBin
                    
                    if data[i,j] != 0:
                        hic.add_edge([iNode, jNode, data[i,j]], flush=False)
        hic.flush(flush_nodes=False)
        
        return hic
            
    def add_node(self, node, flush=True):
        return self.add_region(node, flush)        
    
    def add_edge(self, edge, check_nodes_exist=True, flush=True):
        weight = None
        
        if isinstance(edge, HicEdge):
            source = edge.source
            sink = edge.sink
            weight = edge.weight
        elif type(edge) is dict:
            source = edge['source']
            sink = edge['sink']
            if 'weight' in edge:
                weight = edge['weight']
        else:
            try:
                source = edge[0]
                sink = edge[1]
                if len(edge) > 2:
                    weight = edge[2]
            except TypeError:
                raise ValueError("Edge parameter has to be HicEdge, dict, or list")
        
        if weight is None:
            weight = 1.
        if source > sink:
            tmp = source
            source = sink
            sink = tmp
        
        if check_nodes_exist:
            if source >= len(self._regions) or sink >= len(self._regions):
                raise ValueError("Node index exceeds number of nodes in object")
        
        if weight != 0:
            row = self._edges.row
            row['source'] = source
            row['sink'] = sink
            row['weight'] = weight
            row.append()
        
        if flush:
            self.flush()
    
    def add_nodes(self, nodes):
        self.add_regions(nodes)
    
    
    def add_edges(self, edges):
        for edge in edges:
            self.add_edge(edge, flush=False)
        self.flush(flush_nodes=False)
    
    def merge(self, hic):
        ix_conversion = {}
        
        # merge genomic regions
        for region in hic.regions():
            print region.ix
            ix = self._get_region_ix(region)
            if ix is None:
                ix = self.add_region([region.chromosome, region.start, region.end], flush=False)
            ix_conversion[region.ix] = ix
        self._regions.flush()
                
        # merge edges
        for edge in hic.edges():
            source = ix_conversion[edge.source]
            sink = ix_conversion[edge.sink]
            if source > sink:
                tmp = source
                source = sink
                sink = tmp
            self._update_edge_weight(source, sink, edge.weight, add=True, flush=False)
        self._edges.flush()

            
    def flush(self, flush_nodes=True, flush_edges=True):
        if flush_nodes:
            self._regions.flush()
            # re-indexing not necessary when 'autoindex' is True on table
            if not self._regions.autoindex:
                # reindex node table
                self._regions.flush_rows_to_index()
        if flush_edges:
            self._edges.flush()
            if not self._edges.autoindex:
                # reindex edge table
                self._edges.flush_rows_to_index()
                
    
    
    def __getitem__(self, key):
        """
        Get a chunk of the Hi-C matrix.
        
        Possible key types are:
            Region types
                HicNode: Only the ix of this node will be used for
                    identification
                GenomicRegion: self-explanatory
                str: key is assumed to describe a genomic region
                    of the form: <chromosome>[:<start>-<end>:[<strand>]],
                    e.g.: 'chr1:1000-54232:+'
            Node types
                int: node index
                slice: node range
            List types
                list: This key type allows for a combination of all
                    of the above key types - the corresponding matrix
                    will be concatenated
            
        If the key is a 2-tuple, each entry will be treated as the 
        row and column key, respectively,
        e.g.: 'chr1:0-1000, chr4:2300-3000' will extract the Hi-C
        map of the relevant regions between chromosomes 1 and 4.
            
        """
        
        nodes_ix_row, nodes_ix_col = self._get_nodes_from_key(key, as_index=True)
        
        m = self._get_matrix(nodes_ix_row, nodes_ix_col)
        
        # select the correct output format
        # empty result: matrix
        if m.shape[0] == 0 and m.shape[1] == 0:
            return m
        # both selectors are lists: matrix
        if isinstance(nodes_ix_row, list) and isinstance(nodes_ix_col, list):
            return m
        # row selector is list: vector
        if isinstance(nodes_ix_row, list):
            return m[:,0]
        # column selector is list: vector
        if isinstance(nodes_ix_col, list):
            return m[0,:]
        # both must be indexes
        return m[0,0]
    
    def _get_nodes_from_key(self, key, as_index=False):
        nodes_ix_col = None
        if isinstance(key, tuple):
            nodes_ix_row = self._getitem_nodes(key[0], as_index=as_index)
            nodes_ix_col = self._getitem_nodes(key[1], as_index=as_index)
        else:
            nodes_ix_row = self._getitem_nodes(key, as_index=as_index)
            nodes_ix_col = []
            for row in self._regions:
                if as_index:
                    nodes_ix_col.append(row['ix'])
                else:
                    nodes_ix_col.append(HicNode.from_row(row))
        
        return nodes_ix_row, nodes_ix_col
    
    def _get_matrix(self, nodes_ix_row=None, nodes_ix_col=None):
        # calculate number of rows
        if nodes_ix_row is None:
            n_rows = len(self._regions)
        else:
            if not isinstance(nodes_ix_row, list):
                nodes_ix_row = [nodes_ix_row]
            n_rows = len(nodes_ix_row)
        
        # calculate number of columns
        if nodes_ix_col is None:
            n_cols = len(self._regions)
        else:
            if not isinstance(nodes_ix_col, list):
                nodes_ix_col = [nodes_ix_col]
            n_cols = len(nodes_ix_col)

        # create empty matrix
        m = np.zeros((n_rows, n_cols))
        
        # get row range generator
        row_ranges = ranges(nodes_ix_row)
        
        
        # fill matrix with weights
        row_offset = 0
        for row_range in row_ranges:
            n_rows_sub = row_range[1] - row_range[0] + 1
            col_offset = 0
            col_ranges = ranges(nodes_ix_col)
            for col_range in col_ranges:
                n_cols_sub = col_range[1] - col_range[0] + 1
                
                condition = "((source >= %d) & (source <= %d)) & ((sink >= %d) & (sink <= %d))"
                condition += "| ((source >= %d) & (source <= %d)) & ((sink >= %d) & (sink <= %d))"
                condition = condition % (row_range[0], row_range[1], col_range[0], col_range[1],
                                         col_range[0], col_range[1], row_range[0], row_range[1])
                
                for edge_row in self._edges.where(condition):
                    source = edge_row['source']
                    sink = edge_row['sink']
                    weight = edge_row['weight']
                    ir = source - row_range[0]
                    jr = sink - col_range[0]
                    
                    if (row_range[0] <= source <= row_range[1]
                        and col_range[0] <= sink <= col_range[1]):
                        ir = source - row_range[0]
                        jr = sink - col_range[0]
                        m[ir + row_offset,jr + col_offset] = weight
                        #m[jr + col_offset,ir + row_offset] = weight
                    if (row_range[0] <= sink <= row_range[1]
                        and col_range[0] <= source <= col_range[1]):
                        ir = sink - row_range[0]
                        jr = source - col_range[0]
                        m[ir + row_offset,jr + col_offset] = weight
                
                col_offset += n_cols_sub
            row_offset += n_rows_sub

        return m
    
    def _getitem_nodes(self, key, as_index=False):
        # chr1:1234:56789
        if isinstance(key, str):
            key = GenomicRegion.from_string(key)
        
        # HicNode('chr1', 1234, 56789, ix=0)
        if isinstance(key, HicNode):
            if as_index:
                return key.ix
            else:
                return key
        
        # GenomicRegion('chr1', 1234, 56789) 
        if isinstance(key, GenomicRegion):
            chromosome = key.chromosome
            start = key.start
            end = key.end
            
            # check defaults
            if chromosome is None:
                raise ValueError("Genomic region must provide chromosome name")
            if start is None:
                start = 0
            if end is None:
                end = max(row['end'] for row in self._regions.where("(chromosome == '%s')" % chromosome))
            
            condition = "(chromosome == '%s') & (end >= %d) & (start <= %d)" % (chromosome, start, end)
            if as_index:
                region_nodes = [row['ix'] for row in self._regions.where(condition)]
            else:
                region_nodes = [HicNode.from_row(row) for row in self._regions.where(condition)]
            
            return region_nodes
        
        # 1:453
        if isinstance(key, slice):
            if as_index:
                return [row['ix'] for row in self._regions.iterrows(key.start, key.stop, key.step)]
            else:
                return [HicNode.from_row(row) for row in self._regions.iterrows(key.start, key.stop, key.step)]
        
        # 432
        if isinstance(key, int):
            row = self._regions[key]
            if as_index:
                return row['ix']
            else:
                return HicNode.from_row(row)
        
        # [item1, item2, item3]
        all_nodes_ix = []
        for item in key:
            nodes_ix = self._getitem_nodes(item, as_index=as_index)
            if isinstance(nodes_ix, list):
                all_nodes_ix += nodes_ix
            else:
                all_nodes_ix.append(nodes_ix)
        return all_nodes_ix
    
    def as_data_frame(self, key):
        nodes_ix_row, nodes_ix_col = self._get_nodes_from_key(key, as_index=True)
        nodes_row, nodes_col = self._get_nodes_from_key(key, as_index=False)
        m = self._get_matrix(nodes_ix_row, nodes_ix_col)
        labels_row = []
        for node in nodes_row:
            labels_row.append(node.start)
        labels_col = []
        for node in nodes_col:
            labels_col.append(node.start)
        df = p.DataFrame(m, index=labels_row, columns=labels_col)
        
        return df
    
    def __setitem__(self, key, item):
        """
        Set a chunk of the Hi-C matrix.
        
        Possible key types are:
            Region types
                HicNode: Only the ix of this node will be used for
                    identification
                GenomicRegion: self-explanatory
                str: key is assumed to describe a genomic region
                    of the form: <chromosome>[:<start>-<end>:[<strand>]],
                    e.g.: 'chr1:1000-54232:+'
            Node types
                int: node index
                slice: node range
            List types
                list: This key type allows for a combination of all
                    of the above key types - the corresponding matrix
                    will be concatenated
            
        If the key is a 2-tuple, each entry will be treated as the 
        row and column key, respectively,
        e.g.: 'chr1:0-1000, chr4:2300-3000' will set the Hi-C
        map of the relevant regions between chromosomes 1 and 4.
            
        """
        
        nodes_ix_row, nodes_ix_col = self._get_nodes_from_key(key, as_index=True)
                
        # select the correct output format
        # both selectors are lists: matrix
        if isinstance(nodes_ix_row, list) and isinstance(nodes_ix_col, list):
            n_rows = len(nodes_ix_row)
            n_cols = len(nodes_ix_col)
            # check that we have a matrix with the correct dimensions
            if (not isinstance(item, np.ndarray) or 
                not np.array_equal(item.shape, [n_rows,n_cols])):
                raise ValueError("Item is not a numpy array with shape (%d,%d)!" % (n_rows,n_cols))
            
            for i in xrange(0, n_rows):
                for j in xrange(0,n_cols):
                    source = nodes_ix_row[i]
                    sink = nodes_ix_col[j]
                    weight = item[i,j]
                    self._update_edge_weight(source, sink, weight, flush=False)

        # row selector is list: vector
        elif isinstance(nodes_ix_row, list):
            n_rows = len(nodes_ix_row)
            if (not isinstance(item, np.ndarray) or 
                not np.array_equal(item.shape, [n_rows])):
                raise ValueError("Item is not a numpy vector of length %d!" % (n_rows))
            
            for i, sink in enumerate(nodes_ix_row):
                source = nodes_ix_col
                weight = item[i]
                self._update_edge_weight(source, sink, weight, flush=False)
        
        # column selector is list: vector
        elif isinstance(nodes_ix_col, list):
            n_cols = len(nodes_ix_col)
            if (not isinstance(item, np.ndarray) or 
                not np.array_equal(item.shape, [n_cols])):
                raise ValueError("Item is not a numpy vector of length %d!" % (n_cols))
            
            for i, source in enumerate(nodes_ix_col):
                sink = nodes_ix_row
                weight = item[i]
                self._update_edge_weight(source, sink, weight, flush=False)

        # both must be indexes
        else:
            weight = item
            self._update_edge_weight(nodes_ix_row, nodes_ix_col, weight, flush=False)
        
        self._edges.flush()
        self._remove_zero_edges()
        
        #return m[0,0]
    
    def _update_edge_weight(self, source, sink, weight, add=False, flush=True):
        if source > sink:
            tmp = source
            source = sink
            sink = tmp
        
        value_set = False
        for row in self._edges.where("(source == %d) & (sink == %d)" % (source, sink)):
            original = 0
            if add:
                original = row['weight']
            row['weight'] = weight + original
            row.update()
            value_set = True
            if flush:
                self._edges.flush()
        if not value_set:
            self.add_edge(HicEdge(source=source,sink=sink,weight=weight), flush=flush)
    
    def _remove_zero_edges(self, flush=True):
        zero_edge_ix = []
        ix = 0
        for row in self._edges.iterrows():
            if row['weight'] == 0:
                zero_edge_ix.append(ix)
            ix += 1
        
        for ix in reversed(zero_edge_ix):
            self._edges.remove_row(ix)        
        
        if flush:
            self._edges.flush()
    
    def autoindex(self, index=None):
        if index is not None:
            self._regions.autoindex = bool(index)
            self._edges.autoindex = bool(index)
            return index
        return self._regions.autoindex
            
        
    def save(self, file_name, table_name_nodes='nodes', table_name_edges='edges',
                              table_name_meta='meta', table_name_mask='mask'):
        self.file.copy_file(file_name)
        self.file.close()
        self.file = create_or_open_pytables_file(file_name)
        self._regions = self.file.get_node('/' + table_name_nodes)
        self._edges = self.file.get_node('/' + table_name_edges)
        self._meta = self.file.get_node('/' + table_name_meta)
        self._mask = self.file.get_node('/' + table_name_mask)
    
    
    def get_node(self, key):
        found_nodes = self.get_nodes(key)
        if isinstance(found_nodes, list): 
            raise IndexError("More than one node found matching %s" % str(key))
        return None
        
    
    def get_nodes(self, key):
        return self._getitem_nodes(key)

    def get_edge(self, ix):
        row = self._edges[ix]
        return HicEdge.from_row(row)
    
    def nodes(self):
        return self.regions()
    
    def edges(self):
        hic = self
        class EdgeIter:
            def __init__(self):
                self.iter = iter(hic._edges)
                
            def __iter__(self):
                return self
            
            def next(self):
                return HicEdge.from_row(self.iter.next())
            
            def __len__(self):
                return len(hic._edges)
        return EdgeIter()
    
    
        
        
class HicXmlFile(object):
    def __init__(self, file_name):
        self.file_name = file_name
    
    
    def nodes(self):
        file_name = self.file_name
        
        class XmlNodeIter:
            def __init__(self, file_name):
                self.iter = et.iterparse(file_name)
                
            def __iter__(self):
                return self
            
            def next(self):
                event, elem = self.iter.next()  # @UnusedVariable
                while elem.tag != "node":
                    elem.clear()
                    event, elem = self.iter.next()  # @UnusedVariable
            
                a = elem.attrib
                ix = None
                if 'ix' in a:
                    ix = int(a['ix'])
                    
                chromosome = None
                if 'chromosome' in a:
                    chromosome = a['chromosome']
                
                if 'start' not in a:
                    raise ValueError("start must be a node attribute")
                start = int(a['start'])
                
                if 'end' not in a:
                    raise ValueError("end must be a node attribute")
                end = int(a['end'])
                
                elem.clear()
                return HicNode(ix=ix, chromosome=chromosome, start=start, end=end)
            
        return XmlNodeIter(file_name)
    
    
    def edges(self):
        file_name = self.file_name
        
        class XmlEdgeIter:
            def __init__(self, file_name):
                self.iter = et.iterparse(file_name)
                
            def __iter__(self):
                return self
            
            def next(self):
                event, elem = self.iter.next()  # @UnusedVariable
                while elem.tag != "edge":
                    elem.clear()
                    event, elem = self.iter.next()  # @UnusedVariable
            
                a = elem.attrib
                    
                weight = 1.
                if 'weight' in a:
                    weight = float(a['weight'])
                
                if 'source' not in a:
                    raise ValueError("source must be an edge attribute")
                source = int(a['source'])
                
                if 'sink' not in a:
                    raise ValueError("sink must be an edge attribute")
                sink = int(a['sink'])
                
                elem.clear()
                return HicEdge(source=source, sink=sink, weight=weight)
            
        return XmlEdgeIter(file_name)
        
