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

class Bed(object):
    '''
    Bed object for genomic features
    '''


    def __init__(self, file_name=None, name=None):

        inMemory = False
        h5file_name = file_name
        isFlatFile = False
        if is_bed_file(file_name):
            isFlatFile = True
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
                    if fields[i] == 'strand':
                        desc[fields[i]] = t.StringCol(1) # @UndefinedVariable
                        headerTypes.append(str)
                    if fields[i] == 'itemRgb':
                        desc[fields[i]] = t.StringCol(12) # @UndefinedVariable
                        headerTypes.append(str)
                    if fields[i] == 'chrom':
                        desc[fields[i]] = t.StringCol(16) # @UndefinedVariable
                        headerTypes.append(str)
                    if fields[i] == 'score':
                        desc[fields[i]] = t.Float32Col() # @UndefinedVariable
                        headerTypes.append(float)
                    if (fields[i] == 'thickStart' or fields[i] == 'thickEnd' or 
                        fields[i] == 'blockCount' or fields[i] == 'start' or
                        fields[i] == 'chromStart' or fields[i] == 'end' or
                        fields[i] == 'chromEnd'):
                        desc[fields[i]] = t.Int64Col() # @UndefinedVariable
                        headerTypes.append(int)
                    header.append(fields[i])
                line = f.readline()
                fields = line.rstrip().split("\t")
            else:
                
                for i in range(0,len(fields)):
                    if i == 0:
                        header.append('chrom')
                        headerTypes.append(str)
                        desc['chrom'] = t.StringCol(16) # @UndefinedVariable
                    if i == 1:
                        header.append('start')
                        headerTypes.append(str)
                        desc['start'] = t.Int64Col() # @UndefinedVariable
                    if i == 2:
                        header.append('end')
                        headerTypes.append(str)
                        desc['end'] = t.Int64Col(255) # @UndefinedVariable
                    if i == 3:
                        header.append('name')
                        headerTypes.append(str)
                        desc['name'] = t.StringCol(255) # @UndefinedVariable
                    if i == 4:
                        header.append('score')
                        headerTypes.append(float)
                        desc['score'] = t.Float32Col() # @UndefinedVariable
                    if i == 5:
                        header.append('strand')
                        headerTypes.append(str)
                        desc['strand'] = t.StringCol(1) # @UndefinedVariable
                    if i == 6:
                        header.append('thickStart')
                        headerTypes.append(int)
                        desc['thickStart'] = t.Int64Col() # @UndefinedVariable
                    if i == 7:
                        header.append('thickEnd')
                        headerTypes.append(int)
                        desc['thickEnd'] = t.Int64Col() # @UndefinedVariable
                    if i == 8:
                        header.append('itemRgb')
                        headerTypes.append(str)
                        desc['itemRgb'] = t.StringCol(12) # @UndefinedVariable
                    if i == 9:
                        header.append('blockCount')
                        headerTypes.append(int)
                        desc['blockCount'] = t.Int64Col() # @UndefinedVariable
                    if i == 10:
                        header.append('blockSizes')
                        headerTypes.append(str)
                        desc['blockSizes'] = t.StringCol(255) # @UndefinedVariable
                    if i == 11:
                        header.append('blockStarts')
                        headerTypes.append(str)
                        desc['blockStarts'] = t.StringCol(255) # @UndefinedVariable
            
            print headerTypes
            print header
            print desc
            
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
    
    
    
    
    
    
    
    
    
    
class Bedpe(object):
    '''
    Bedpe object for genomic features
    '''


    def __init__(self, file_name=None, name=None):
        
        inMemory = False
        h5file_name = file_name
        isFlatFile = False
        if is_bedpe_file(file_name):
            isFlatFile = True
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
                        fields[i] == 'chromStart1' or fields[i] == 'chromStart2' or
                        fields[i] == 'end1' or fields[i] == 'end2' or
                        fields[i] == 'chromEnd1' or fields[i] == 'chromEnd2'):
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
                        desc['chrom'] = t.StringCol(16) # @UndefinedVariable
                    elif i == 1:
                        header.append('start1')
                        headerTypes.append(str)
                        desc['start'] = t.Int64Col() # @UndefinedVariable
                    elif i == 2:
                        header.append('end1')
                        headerTypes.append(str)
                        desc['end'] = t.Int64Col() # @UndefinedVariable
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
            
            print headerTypes
            print header
            print desc
            
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
    



class Hic(Bedpe):
    def __init__(self, file_name=None, name=None):
        Bedpe.__init__(self, file_name=file_name, name=name)
        
    def as_data_frame(self):
        print "test"
        
        
    def as_structured_array(self, resolution, chrom=None, lower_bound=None, upper_bound=None):
        
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
        
        if query != '':
            contacts = [[x['start1'],x['start2'],x['score']] for x in self.table.where(query)]
        else:
            contacts = [[x['start1'],x['start2'],x['score']] for x in self.table]
        
        
        min_ix = min(min(contacts)[0:2])
        max_ix = max(max(contacts)[0:2])
        
        labels = range(min_ix,max_ix+resolution,resolution)
        
        df = p.DataFrame(np.zeros((len(labels),len(labels))), index=labels, columns=labels)
        
        for c in contacts:
            df[c[0],c[1]] = c[2]
            df[c[1],c[0]] = c[2]
        
        return df
        
        
        
        
        