'''
Created on May 20, 2015

@author: kkruse1
'''

import tables as t
import pandas as p
import numpy as np
from kaic.tools.files import create_or_open_pytables_file
from kaic.tools.files import is_bed_file
from kaic.tools.files import is_bedpe_file
import string
import random
from Bio import SeqIO, Restriction, Seq
from kaic.data.general import Table

import logging
logging.basicConfig(level=logging.INFO)

class Genome(object):
    def __init__(self, file_name=None):
        logging.info("Creating new genome object")
        
        # chromosomes
        #   length
        #   name
        #   (sequence)


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
            for i in range(0,len(data)):
                self.data[names[i]] = data[i]
                self.id2ix[names[i]] = i
                self.ix2id[i] = names[i]
                


class BedImproved(Table):
    
    @staticmethod
    def col_type(name,pos=None):
        col_type = {
            'chrom': (0, str, t.StringCol(16,pos=pos)), # @UndefinedVariable
            'start': (1, int, t.Int64Col(pos=pos)), # @UndefinedVariable
            'end': (2, int, t.Int64Col(pos=pos)), # @UndefinedVariable
            'name': (3, str, t.StringCol(255,pos=pos)), # @UndefinedVariable
            'score': (4, float, t.Float32Col(pos=pos)), # @UndefinedVariable
            'strand': (5, str, t.StringCol(2,pos=pos)), # @UndefinedVariable
            'thickStart': (6, int, t.Int64Col(pos=pos)), # @UndefinedVariable
            'thickEnd': (7, int, t.Int64Col(pos=pos)), # @UndefinedVariable
            'itemRgb': (8, str, t.StringCol(12,pos=pos)), # @UndefinedVariable
            'blockCount': (9, int, t.Int64Col(pos=pos)), # @UndefinedVariable
            'blockSizes': (10, str, t.StringCol(255,pos=pos)), # @UndefinedVariable
            'blockStarts': (11, str, t.StringCol(255,pos=pos)) # @UndefinedVariable
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
                    header.append("feature_\d" % i)
                
            for i in range(0,len(header)):
                pos, ptype, ttype = BedImproved.col_type(header[i],i+1)
                col_types.append(ttype)
                headerTypes.append(ptype)
                
            
            data = []
            while line != '':
                d = {}
                for i in range(0,len(fields)):
                    d[header[i]] = headerTypes[i](fields[i])
                data.append(d)
                    
                line = f.readline()
                fields = line.rstrip().split(sep)
            
            bed = cls()
            
            super(BedImproved, bed).__init__(col_names=header, col_types=col_types, data=data, name=name)
            
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
        desc = self._table.description._v_colObjects.copy()
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
            
            desc = self.table.description._v_colObjects.copy()
            header = []
            headerTypes = []
            if has_header:
                for i in range(0,len(fields)):
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
                
                for i in range(0,len(fields)):
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
            
            table2 = self.file.createTable(self.file.root, 'table2', desc, "bed", t.Filters(1))
 
            # Copy the user attributes
            self.table.attrs._f_copy(table2)
             
            # Fill the rows of new table with default values
            for i in xrange(self.table.nrows):
                table2.row.append()
            # Flush the rows to disk
            table2.flush()
            
            # Copy the columns of source table to destination
            for col in self.table.description._v_colObjects:
                getattr(table2.cols, col)[:] = getattr(self.table.cols, col)[:]
             
            # fill with new data
            entry = table2.row
            while line != '':
                
                if len(fields) == len(headerTypes):
                    for i in range(0,len(fields)):
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
        desc = self.table.description._v_colObjects.copy()
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
                for i in range(0,len(fields)):
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
                
                for i in range(0,len(fields)):
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
            print self.table.cols
            print table2.cols
            for col in self.table.description._v_colObjects:
                print len(self.table.cols)
                print getattr(self.table.cols, col)[:]
                print len(getattr(table2.cols, col)[:])
                print getattr(table2.cols, col)[:]
                getattr(table2.cols, col)[:] = getattr(self.table.cols, col)[:]
             
            # fill with new data
            entry = table2.row
            while line != '':
                
                if len(fields) == len(headerTypes):
                    for i in range(0,len(fields)):
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
        
        print "Calculating lowest bound"
        min_ix = min(min(contacts)[0:2])
        max_ix = max(max(contacts)[0:2])
        print "Calculating upper bound"
        
        labels = range(min_ix,max_ix+resolution,resolution)
        ix_l = int(min_ix/resolution)
        
        print "Assigning to matrix"
        M = np.zeros((len(labels),len(labels)))
        for c in contacts:
            i = int(c[0]/resolution)-ix_l
            j = int(c[1]/resolution)-ix_l
            M[i,j] = c[2]
            M[j,i] = c[2]
        
        return M
    
    def directionality(self, resolution, window_size=2000000):
        M = self.as_matrix(resolution)
        bin_window_size = int(window_size/resolution)
        if window_size%resolution > 0:
            bin_window_size += 1
        
        n_bins = M.shape[0]
        dis = np.zeros(n_bins,dtype='float64')
        for i in range(0,n_bins):
            start = max(0, i-bin_window_size)
            end = min(n_bins, i+window_size, i+i-start)
            
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
    
    @classmethod
    def from_fasta(cls, file_name, name=None, include_sequence=True):
        fastas = SeqIO.parse(open(file_name,'r'),'fasta')
        fasta = fastas.next()
        
        if include_sequence:
            return cls(name if name else fasta.id, length=len(fasta), sequence=fasta.seq.tostring())
        else:
            return cls(name if name else fasta.id, length=len(fasta))
            
        
        
    def get_restriction_sites(self, restriction_enzyme, include_chromosome_end=True, include_chromosome_start=False):
        logging.info("Calculating RE sites")
        re = eval('Restriction.%s' % restriction_enzyme)
        return re.search(Seq.Seq(self.sequence))
    
    def __repr__(self):
        return "Name: %s\nLength: %d\nSequence: %s" % (self.name if self.name else '',
                                                       self.length if self.length else -1,
                                                       self.sequence[:20] + "..." if self.sequence else '')
         
        