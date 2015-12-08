"""
Module for working with genomic data.

This module provides classes and functions to work with objects in the context
of the genome.

:class:`~Chromosome`, :class:`~Genome`, and :class:`~GenomicRegion` simplify
working with reference sequence data by providing easy access and many convenience
functions.

Examples:

.. code:: python

    # assemble custom genome
    chr1 = Chromosome.from_fasta("/path/to/chr1_fasta_file")
    chr2 = Chromosome.from_fasta("/path/to/chr2_fasta_file")
    genome = Genome(chromosomes=[chr1,chr2])

    # extract genomic regions binned every 10000bp
    regions = genome.get_regions(10000)

.. code:: python

    # assemble genome from folder with FASTA files
    genome = Genome.from_folder("/path/to/fasta_folder/")

:class:`~Hic` is the central class for working with Hi-C data. It provides
matrix-like selectors and convenient access to specific genomic regions. In the
kaic pipeline, a Hic object is assembled at the fragment level from
:class:`~kaic.construct.seq.FragmentMappedReadPairs`. From there, it can be
binned to equi-distant genomic regions.

.. code:: python

    # use previously existing FragmentMappedReadPairs object 'pairs'
    hic = Hic(file_name="/path/to/save_file")
    hic.load_read_fragment_pairs(pairs)

    # bin Hi-C object
    binned = hic.bin(10000)

    # ... further processing


Alternatively, a Hic object can be assembled from scratch using genomic
regions and edges (contacts) between them.

Example:

.. code:: python

    hic = Hic()
    genome = Genome.from_folder("/path/to/fasta_folder")
    hic.add_regions(genome.get_regions(10000))

    hic.add_edges(list_of_edges)

"""

from __future__ import division
import tables as t
import pandas as p
import numpy as np
from kaic.tools.files import create_or_open_pytables_file, is_hic_xml_file,\
    is_fasta_file, is_hdf5_file
from kaic.tools.files import is_bed_file, is_bedpe_file
from Bio import SeqIO, Restriction, Seq
from kaic.data.general import Table, TableRow, TableArray, TableObject,\
    MetaContainer, Maskable, MaskedTable, FileBased, MaskFilter
from abc import abstractmethod, ABCMeta
import os.path
import logging
from kaic.tools.general import ranges, distribute_integer
from itertools import izip as zip
from xml.etree import ElementTree as et
logging.basicConfig(level=logging.INFO)


def _edge_overlap_split_rao(original_edge, overlap_map):
    """
    Resolve the distribution of contacts when binning using
    Rao et al. 2014 approach.
    """
    original_source = original_edge[0]
    original_sink = original_edge[1]
    original_weight = original_edge[2]

    new_source_nodes = overlap_map[original_source]
    new_sink_nodes = overlap_map[original_sink]
    
    if len(new_source_nodes) == 0:
        return []
    elif len(new_source_nodes) == 1:
        new_source_nodes = [new_source_nodes[0][0]]
    else:
        new_source_nodes = [new_source_nodes[0][0], new_source_nodes[-1][0]]
    
    if len(new_sink_nodes) == 0:
        return []
    elif len(new_sink_nodes) == 1:
        new_sink_nodes = [new_sink_nodes[0][0]]
    else:
        new_sink_nodes = [new_sink_nodes[0][0], new_sink_nodes[-1][0]]
    
    edges = {}
    for new_source in new_source_nodes:
        for new_sink in new_sink_nodes:
            if new_source <= new_sink:
                edges[(new_source, new_sink)] = 0
            else:
                edges[(new_sink, new_source)] = 0
    
    weights = distribute_integer(original_weight, len(edges))
    edges_list = []
    for i, key_pair in enumerate(edges):
        edges_list.append([key_pair[0], key_pair[1], weights[i]])

    return edges_list


class Bed(Table):
    """
    Data type representing a BED file.

    This class is an extension of :class:`~kaic.data.general.Table`, and
    can be used to load and represent BED-formatted data.
    """

    def __init__(self, data=None, colnames=None, col_types=None):
        """
        Initialize Bed object (as Table).
        """
        super(Bed, self).__init__(data=data, colnames=colnames, col_types=col_types)
    
    @staticmethod
    def col_type(name, pos=None):
        """
        Determine the column type (PyTables) by its name.

        :param name: The name of the column
        :return: PyTables column type
        """
        col_type = {
            'chrom': (str, t.StringCol(16, pos=pos)),
            'start': (int, t.Int64Col(pos=pos)),
            'end': (int, t.Int64Col(pos=pos)),
            'name': (str, t.StringCol(255, pos=pos)),
            'score': (float, t.Float32Col(pos=pos)),
            'strand': (str, t.StringCol(2, pos=pos)),
            'thickStart': (int, t.Int64Col(pos=pos)),
            'thickEnd': (int, t.Int64Col(pos=pos)),
            'itemRgb': (str, t.StringCol(12, pos=pos)),
            'blockCount': (int, t.Int64Col(pos=pos)),
            'blockSizes': (str, t.StringCol(255, pos=pos)),
            'blockStarts': (str, t.StringCol(255, pos=pos))
        }
        if name in col_type:
            return col_type[name]
        return str, t.StringCol(255)
    
    @classmethod
    def from_bed_file(cls, file_name, has_header=True, sep="\t"):
        if not is_bed_file:
            raise ImportError("File does not appear to be a BED file")
        
        all_fields = ['chrom', 'start', 'end', 'name', 'score', 'strand',
                      'thickStart', 'thickEnd', 'itemRgb', 'blockCount',
                      'blockSizes', 'blockStarts']
        
        with open(file_name, 'r') as f:
            
            # process first line, update table
            line = f.readline()
            fields = line.rstrip().split(sep)
            col_types = []
            header_types = []
            if has_header:
                header = fields[:]
                line = f.readline()
                fields = line.rstrip().split(sep)
            else:
                header = all_fields[0:len(fields)]
                for i in (len(header), len(fields)):
                    header.append("feature_%d" % i)
                
            for i, name in enumerate(header):
                ptype, ttype = Bed.col_type(name, i+1)
                col_types.append(ttype)
                header_types.append(ptype)

            data = []
            while line != '':
                d = {}
                for i, field in enumerate(fields):
                    # ignore dot by default
                    if field == '.':
                        d[header[i]] = header_types[i]()
                    else:
                        d[header[i]] = header_types[i](field)
                data.append(d)
                    
                line = f.readline()
                fields = line.rstrip().split(sep)

            bed = cls(colnames=header, col_types=col_types)
            bed.append(data)
            
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
                query += "(start >= %d)" % start
        
        if end:
            if query != '':
                query += " & "
            if type(end) is list:
                query += "(end >= %d) & (end <= %d)" % (end[0], end[1])
            else:
                query += "(end <= %d)" % end

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

    def _row_to_bed_element(self, row):
        kwargs = dict()
        for key in self.colnames:
            if key == 'chrom' or key == 'start' or key == 'end':
                continue
            kwargs[key] = row[key]
        return BedElement(row['chrom'], row['start'], row['end'], **kwargs)

    def __getitem__(self, item):
        row = super(Bed, self).__getitem__(item)
        return self._row_to_bed_element(row)

    def __iter__(self):
        this = self

        class BedIter:
            def __init__(self):
                self.iter = iter(this._table)

            def __iter__(self):
                return self

            def next(self):
                row = self.iter.next()
                return this._row_to_bed_element(row)
        return BedIter()


class Bedpe(object):
    """
    Bedpe object for genomic features
    """

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
            self.file = create_or_open_pytables_file()
        else:
            self.file = create_or_open_pytables_file(h5file_name)

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
                        desc['chrom1'] = t.StringCol(16)
                    elif i == 1:
                        header.append('start1')
                        headerTypes.append(str)
                        desc['start1'] = t.Int64Col()
                    elif i == 2:
                        header.append('end1')
                        headerTypes.append(str)
                        desc['end1'] = t.Int64Col()
                    elif i == 3:
                        header.append('chrom2')
                        headerTypes.append(str)
                        desc['chrom2'] = t.StringCol(16)
                    elif i == 4:
                        header.append('start2')
                        headerTypes.append(int)
                        desc['start2'] = t.Int64Col()
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
        

        print "Running query"
        if query != '':
            contacts = [[x[y] for y in labels] for x in self.table.where(query)]
        else:
            contacts = [[x[y] for y in labels] for x in self.table]
            
        df = p.DataFrame(contacts, columns=labels)
        return df


class Chromosome(object):
    """
    Chromosome data type.

    .. attribute:: name

        Name of the chromosome

    .. attribute:: length

        Length of the chromosome in base-pairs

    .. attribute:: sequence

        Base-pair sequence of DNA in the chromosome
    """
    def __init__(self, name=None, length=None, sequence=None):
        """
        Initialize chromosome

        :param name: Name of the chromosome
        :param length: Length of the chromosome in base-pairs
        :param sequence: Base-pair sequence of DNA in the chromosome
        """
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
        """
        Get length of the chromosome.
        """
        return self.length
    
    def __getitem__(self, key):
        """
        Get object attributes by name
        """
        if key == 'name':
            return self.name
        if key == 'length':
            return self.length
        if key == 'sequence':
            return self.sequence
    
    @classmethod
    def from_fasta(cls, file_name, name=None, include_sequence=True):
        """
        Create a :class:`~Chromosome` from a FASTA file.

        This class method will load a FASTA file and convert it into
        a :class:`~Chromosome` object. If the FASTA file contains multiple
        sequences, only the first one will be read.

        :param file_name: Path to the FASTA file
        :param name: Chromosome name. If None (default), will be read
                     from the FASTA file header
        :param include_sequence: If True (default), stores the chromosome
                                 sequence in memory. Else, the sequence
                                 attribute will be set to None.
        :return: :class:`~Chromosome` if there is only a single FASTA
                 sequence in the file, list(:class:`~Chromosome`) if
                 there are multiple sequences.
        """
        if type(file_name) is file:
            fastas = SeqIO.parse(file_name, 'fasta')
        else:
            fastas = SeqIO.parse(open(file_name, 'r'), 'fasta')

        chromosomes = []
        for fasta in fastas:
            if include_sequence:
                chromosome = cls(name if name else fasta.id, length=len(fasta), sequence=str(fasta.seq))
            else:
                chromosome = cls(name if name else fasta.id, length=len(fasta))
            chromosomes.append(chromosome)

        if len(chromosomes) == 0:
            raise ValueError("File %s does not appear to be a FASTA file" % file_name)
        if len(chromosomes) == 1:
            return chromosomes[0]
        return chromosomes

    def get_restriction_sites(self, restriction_enzyme):
        """
        Find the restriction sites of a provided enzyme in this chromosome.

        Internally uses biopython to find RE sites.

        :param restriction_enzyme: The name of the restriction enzyme
                                   (e.g. HindIII)
        :return: List of RE sites in base-pairs (1-based)
        """
        logging.info("Calculating RE sites")
        try:
            re = eval('Restriction.%s' % restriction_enzyme)
        except SyntaxError:
            raise ValueError("restriction_enzyme must be a string")
        except AttributeError:
            raise ValueError("restriction_enzyme string is not recognized: %s" % restriction_enzyme)
        
        return re.search(Seq.Seq(self.sequence))


class Genome(Table):
    """
    Class representing a collection of chromosomes.

    Extends the :class:`~kaic.data.general.Table` class and provides
    all the expected functionality. Provides some convenience batch
    methods that call :class:`~Chromosome` methods for every
    chromosome in this object.

    This object can be saved to file.
    """
    def __init__(self, file_name=None, chromosomes=None):
        """
        Build :class:`~Genome` from a list of chromosomes or load
        previously saved object.

        :param file_name: Path of the file to load or to save to.
        :param chromosomes: List of chromosomes to load into this
                            object.
        """
        self.file = create_or_open_pytables_file(file_name)
            
        columns = ["ix", "name", "length"]
        column_types = [t.Int32Col(pos=0), t.StringCol(50, pos=1), t.Int32Col(pos=2)]  # @UndefinedVariable
        Table.__init__(self, colnames=columns, col_types=column_types)
        
        try:
            self._sequences = self.file.get_node('/genome_sequences')
        except t.NoSuchNodeError:
            self._sequences = self.file.create_vlarray("/", 'genome_sequences', t.VLStringAtom())
        
        if chromosomes is not None:
            if isinstance(chromosomes, Chromosome):
                self.add_chromosome(chromosomes)
            else:
                for chromosome in chromosomes:
                    self.add_chromosome(chromosome)
            
    @classmethod
    def from_folder(cls, folder_name, file_name=None, exclude=None, include_sequence=True):
        """
        Load every FASTA file from a folder as a chromosome.

        :param folder_name: Path to the folder to load
        :param file_name: File to save Genome object to
        :param exclude: List or set of chromosome names that
                        should NOT be loaded
        :param include_sequence: If True, will save the
                                 chromosome sequences in the
                                 Genome object
        """
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

    @classmethod
    def from_string(cls, genome_string, file_name=None):
        """
        Convenience function to load a :class:`~Genome` from a string.

        :param genome_string: Path to FASTA file, path to folder with
                              FASTA files, comma-separated list of
                              paths to FASTA files, path to HDF5 file
        :return: A :class:`~Genome` object
        """
        # case 1: FASTA file = Chromosome
        if is_fasta_file(genome_string):
            chromosomes = Chromosome.from_fasta(genome_string)
            genome = cls(chromosomes=chromosomes, file_name=file_name)
        # case 2: Folder with FASTA files
        elif os.path.isdir(genome_string):
            genome = cls.from_folder(genome_string, file_name=file_name)
        # case 3: path to HDF5 file
        elif is_hdf5_file(genome_string):
            genome = cls(genome_string)
        # case 4: List of FASTA files
        else:
            chromosome_files = genome_string.split(',')
            chromosomes = []
            for chromosome_file in chromosome_files:
                chromosome = Chromosome.from_fasta(os.path.expanduser(chromosome_file))
                chromosomes.append(chromosome)
            genome = cls(chromosomes=chromosomes, file_name=file_name)

        return genome

    def __getitem__(self, key):
        """
        Get Genome table subsets.

        If the result is one or more rows, they will be converted to
        :class:`~Chromosome` objects, if the result is a column, it
        will be returned without conversion.
        """
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
        """
        Get iterator over :class:`~Chromosome` objects.
        """
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
        """
        Add a :class:`~Chromosome` to this object.

        Will choose suitable defaults for missing attributes.

        :param chromosome: :class:`~Chromosome` object or similar
                           object (e.g. dict) with the same fields
        """
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
        """
        Extract genomic regions from genome.

        Provides two options:

        - Splits chromosomes at restriction sites if the split
          parameter is the name of a restriction enzyme.

        - Splits chromosomes at equi-distant points if split
          is an integer

        :param split: Name of a restriction enzyme or positive
                      integer
        :param file_name: Name of a file if the result of this
                          method should be saved to file
        :return: :class:`~GenomicRegions`
        """

        regions = GenomicRegions(file_name=file_name)
        for chromosome in self:
            split_locations = []
            if isinstance(split, str):
                split_locations = chromosome.get_restriction_sites(split)
            elif isinstance(split, int):
                for i in xrange(split, len(chromosome)-1, split):
                    split_locations.append(i)
            else:
                for i in split:
                    split_locations.append(i)
            
            for i in xrange(0, len(split_locations)):
                if i == 0:
                    region = GenomicRegion(start=1, end=split_locations[i], chromosome=chromosome.name)
                else:
                    region = GenomicRegion(start=split_locations[i-1]+1,
                                           end=split_locations[i], chromosome=chromosome.name)
                
                regions.add_region(region, flush=False)
                
            # add last node
            if len(split_locations) > 0:
                region = GenomicRegion(start=split_locations[len(split_locations)-1]+1,
                                       end=chromosome.length, chromosome=chromosome.name)
            else:
                region = GenomicRegion(start=1, end=chromosome.length, chromosome=chromosome.name)
            regions.add_region(region, flush=False)
            regions._flush()
                
        return regions


class GenomicRegion(TableObject):
    """
    Class representing a genomic region.

    .. attribute:: chromosome

        Name of the chromosome this region is located on

    .. attribute:: start

        Start position of the region in base pairs

    .. attribute:: end

        End position of the region in base pairs

    .. attribute:: strand

        Strand this region is on (+1, -1)

    .. attribute:: ix

        Index of the region in the context of all genomic
        regions.

    """

    def __init__(self, start, end, chromosome=None, strand=None, ix=None):
        """
        Initialize this object.

        :param start: Start position of the region in base pairs
        :param end: End position of the region in base pairs
        :param chromosome: Name of the chromosome this region is located on
        :param strand: Strand this region is on (+1, -1)
        :param ix: Index of the region in the context of all genomic
                   regions.
        """
        self.start = start
        self.end = end
        self.strand = strand
        self.chromosome = chromosome
        self.ix = ix

    @classmethod
    def from_row(cls, row):
        """
        Create a :class:`~GenomicRegion` from a PyTables row.
        """
        strand = row['strand']
        if strand == 0:
            strand = None
        return cls(start=row["start"], end=row["end"],
                   strand=strand, chromosome=row["chromosome"])

    @classmethod
    def from_string(cls, region_string):
        """
        Convert a string into a :class:`~GenomicRegion`.

        This is a very useful convenience function to quickly
        define a :class:`~GenomicRegion` object from a descriptor
        string.

        :param region_string: A string of the form
                              <chromosome>[:<start>-<end>[:<strand>]]
                              (with square brackets indicating optional
                              parts of the string). If any optional
                              part of the string is omitted, intuitive
                              defaults will be chosen.
        :return: :class:`~GenomicRegion`
        """
        chromosome = None
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
        """
        Convert this :class:`~GenomicRegion` to its string representation.

        :return: str
        """
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

    def overlaps(self, region):
        """
        Check if this region overlaps with the specified region.

        :param region: :class:`~GenomicRegion` object or string
        """
        if isinstance(region, str):
            region = GenomicRegion.from_string(region)

        if region.chromosome != self.chromosome:
            return False

        if region.start <= self.end and region.end >= self.start:
            return True
        return False

    def contains(self, region):
        """
        Check if the specified region is completely contained in this region.

        :param region: :class:`~GenomicRegion` object or string
        """
        if isinstance(region, str):
            region = GenomicRegion.from_string(region)

        if region.chromosome != self.chromosome:
            return False

        if region.start >= self.start and region.end <= self.end:
            return True
        return False


class BedElement(GenomicRegion):
    def __init__(self, chromosome, start, end, **kwargs):
        super(BedElement, self).__init__(start, end, chromosome)
        for key, value in kwargs.iteritems():
            setattr(self, key, value)


class LazyGenomicRegion(GenomicRegion):
    def __init__(self, row, ix=None):
        self.row = row
        self.static_ix = ix

    @property
    def chromosome(self):
        return self.row["chromosome"]

    @property
    def start(self):
        return self.row["start"]

    @property
    def end(self):
        return self.row["end"]

    @property
    def strand(self):
        return self.row["strand"]

    @property
    def ix(self):
        if self.static_ix is None:
            return self.row["ix"]
        return self.static_ix


class GenomicRegions(Table):
    """
    A collection of :class:`~GenomicRegion` objects.
    """
    class GenomicRegionDescription(t.IsDescription):
        """
        Description for PyTables Table representing a :class:`~GenomicRegion`.
        """
        chromosome = t.StringCol(50, pos=0)
        start = t.Int64Col(pos=1)
        end = t.Int64Col(pos=2)
        strand = t.Int8Col(pos=3)
        
    def __init__(self, file_name=None, regions=None):
        """
        :param file_name: Path to a file this object will be saved to.
        :param regions: A list of :class:`~GenomicRegion` objects.
        """
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
        """
        Add a :class:`~GenomicRegion` to this object.

        :param region: A :class:`~GenomicRegion` object or any object
                       with the same attributes (at least chromosome,
                       start, end, strand)
        """

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
    """
    PyTables Table wrapper for storing genomic regions.

    This class is inherited by objects working with lists of genomic
    regions, such as equi-distant bins along chromosomes in a genome
    (:class:`~Hic`) or restriction fragments of genomic DNA
    (:class:`~kaic.construct.seq.FragmentMappedReadPairs`)
    """

    class RegionDescription(t.IsDescription):
        """
        Description of a genomic region for PyTables Table
        """
        ix = t.Int32Col(pos=0)
        chromosome = t.StringCol(50, pos=1)
        start = t.Int64Col(pos=2)
        end = t.Int64Col(pos=3)
    
    def __init__(self, data=None, file_name=None,
                 mode='a',
                 _table_name_regions='regions'):
        """
        Initialize region table.

        :param data: List of regions to load in object. Can also
                     be used to load saved regions from file by
                     providing a path to an HDF5 file and setting
                     the file_name parameter to None.
        :param file_name: Path to a save file.
        :param _table_name_regions: (Internal) name of the HDF5
                                    node that stores data for this
                                    object
        """
        
        # parse potential unnamed argument
        if data is not None:
            # data is file name
            if type(data) is str or isinstance(data, t.file.File):                
                if file_name is None:
                    file_name = data
                    data = None
        
        if file_name is not None and isinstance(file_name, str):
            file_name = os.path.expanduser(file_name)
        
        FileBased.__init__(self, file_name, mode=mode)
        
        # check if this is an existing Hi-C file
        if _table_name_regions in self.file.root:
            self._regions = self.file.get_node('/', _table_name_regions)
            if len(self._regions) > 0:
                self._max_region_ix = max(row['ix'] for row in self._regions.iterrows())
            else:
                self._max_region_ix = -1
        else:
            self._regions = t.Table(self.file.root, _table_name_regions,
                                    RegionsTable.RegionDescription)
            self._max_region_ix = -1

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

        self._ix_to_chromosome = dict()
        self._chromosome_to_ix = dict()

        if data is not None:
            self.add_regions(data)
        else:
            self._update_references()

    def add_region(self, region, flush=True):
        """
        Add a genomic region to this object.

        This method offers some flexibility in the types of objects
        that can be loaded. See below for details.

        :param region: Can be a :class:`~GenomicRegion`, a dict with
                       at least the fields 'chromosome', 'start', and
                       'end', optionally 'ix', or a list of length 3
                       (chromosome, start, end) or 4 (ix, chromosome,
                       start, end).
        :param flush: If True, data will be written to file and made
                      available immediately. For bulk inserts use
                      :func:`~RegionsTable.add_regions`.
        """
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
            self._update_references()
        
        return ix

    def _update_references(self):
        chromosomes = self.chromosomes()
        for i, chromosome in enumerate(chromosomes):
            self._ix_to_chromosome[i] = chromosome
            self._chromosome_to_ix[chromosome] = i

    def add_regions(self, regions):
        """
        Bulk insert multiple genomic regions.

        :param regions: List (or any iterator) with objects that
                        describe a genomic region. See
                        :class:`~RegionsTable.add_region` for options.
        """
        for region in regions:
            self.add_region(region, flush=False)
        self._regions.flush()
        self._update_references()

    def _get_region_ix(self, region):
        """
        Get index from other region properties (chromosome, start, end)
        """
        condition = "(start == %d) & (end == %d) & (chromosome == '%s')"
        condition = condition % (region.start, region.end, region.chromosome)
        for res in self._regions.where(condition):
            return res["ix"]
        return None

    @staticmethod
    def _row_to_region(row, lazy=False):
        if lazy:
            return LazyGenomicRegion(row)
        return GenomicRegion(chromosome=row["chromosome"], start=row["start"],
                             end=row["end"], ix=row["ix"])

    def regions(self):
        """
        Iterate over genomic regions in this object.

        Will return a :class:`~HicNode` object in every iteration.
        Can also be used to get the number of regions by calling
        len() on the object returned by this method.

        :return: RegionIter
        """
        this = self

        class RegionIter:
            def __init__(self):
                self.iter = iter(this._regions)
                
            def __iter__(self):
                return self
            
            def next(self):
                return RegionsTable._row_to_region(self.iter.next())
            
            def __len__(self):
                return len(this._regions)
            
        return RegionIter()

    def chromosomes(self):
        """
        Get a list of chromosome names.

        :return:
        """
        chromosomes_set = set()
        chromosomes = []
        for region in self.regions():
            if region.chromosome not in chromosomes_set:
                chromosomes_set.add(region.chromosome)
                chromosomes.append(region.chromosome)
        return chromosomes


class HicNode(GenomicRegion, TableObject):
    """
    Class representing a node in a :class:`~Hic` object.

    Backed by a :class:`~GenomicRegion`, this class additionally
    provides methods to access the node index in the context of
    the :class:`~Hic` object.

    .. attribute:: chromosome

        Name of the chromosome this region is located on

    .. attribute:: start

        Start position of the region in base pairs

    .. attribute:: end

        End position of the region in base pairs

    .. attribute:: strand

        Strand this region is on (+1, -1)

    .. attribute:: ix

        Index of the region in the context of all genomic
        regions.
    """
    def __init__(self, chromosome=None, start=None, end=None, ix=None):
        self.ix = ix
        super(HicNode, self).__init__(chromosome=chromosome, start=start, end=end, ix=ix)
    
    def __repr__(self):
        if self.ix is None:
            return "%s, %d-%d" % (self.chromosome, self.start, self.end)
        else:
            return "%d: %s, %d-%d" % (self.ix, self.chromosome, self.start, self.end)


class LazyHicNode(LazyGenomicRegion, HicNode):
    def __init__(self, row, ix=None):
        LazyGenomicRegion.__init__(self, row=row, ix=ix)


class HicEdge(TableObject):
    """
    A contact / an Edge between two genomic regions.

    .. attribute:: source

        The index of the "source" genomic region. By convention,
        source <= sink.

    .. attribute:: sink

        The index of the "sink" genomic region.

    .. attribute:: weight

        The weight or contact strength of the edge. Can, for
        example, be the number of reads mapping to a contact.
    """
    def __init__(self, source, sink, weight=1):
        """
        :param source: The index of the "source" genomic region
                       or :class:`~HicNode` object.
        :param sink: The index of the "sink" genomic region
                     or :class:`~HicNode` object.
        :param weight: The weight or contact strength of the edge.
        """
        self._source = source
        self._sink = sink
        self.weight = weight

    @property
    def source(self):
        if isinstance(self._source, GenomicRegion):
            return self._source.ix
        return self._source

    @property
    def sink(self):
        if isinstance(self._sink, GenomicRegion):
            return self._sink.ix
        return self._sink

    @property
    def source_node(self):
        if isinstance(self._source, GenomicRegion):
            return self._source
        raise RuntimeError("Source not not provided during object initialization!")

    @property
    def sink_node(self):
        if isinstance(self._sink, GenomicRegion):
            return self._sink
        raise RuntimeError("Sink not not provided during object initialization!")

    def __repr__(self):
        return "%d--%d (%.2f)" % (self.source, self.sink, self.weight)
    
    @classmethod
    def from_row(cls, row):
        return cls(source=row['source'], sink=row['sink'], weight=row['weight'])


class LazyHicEdge(HicEdge):
    def __init__(self, row, nodes_table):
        self._row = row
        self._nodes_table = nodes_table
        self._source_node = None
        self._sink_node = None

    @property
    def weight(self):
        return self._row['weight']

    @property
    def source(self):
        return self._row['source']

    @property
    def sink(self):
        return self._row['sink']

    @property
    def source_node(self):
        if self._source_node is None:
            source_row = self._nodes_table[self.source]
            return LazyHicNode(source_row)
        return self._source_node

    @property
    def sink_node(self):
        if self._sink_node is None:
            sink_row = self._nodes_table[self.sink]
            return LazyHicNode(sink_row)
        return self._sink_node


class Hic(Maskable, MetaContainer, RegionsTable, FileBased):
    """
    Class for working with Hi-C data.

    Generally, a Hi-C object has two components:

    - Nodes or regions: (Non-overlapping) genomic regions
      obtained by splitting the genome into distinct pieces.
      See also :class:`~GenomicRegion` and :class:`~RegionsTable`

    - Edges or contacts: Pairs of genomic regions with optionally
      associated weight or contact strength. See also
      :class:`~HicEdge`

    This is a memory-efficient implementation of a Hi-C data
    container. Internally, this is achieved by saving entries
    of the Hi-C matrix in sparse notation, i.e. in a list of
    non-zero contacts.

    Its bracket-notation access behaves like a numpy
    array and handles data retrieval and assignment in matrix-
    fashion, e.g. hic[1:3] would return rows 1 and 2 of
    the Hi-C matrix (0-based index). However, the bracket
    notation can also handle :class:`~GenomicRegion` descriptior
    strings, i.e. hic['chr1','chr5'] will extract the inter-
    chromosomal matrix between chromosomes 1 and 5 only.

    Examples:

    .. code:: python

        hic = Hic(file_name="/path/to/save/file")

        # load genomic regions
        genome = Genome.from_folder("/path/to/fasta/folder")
        regions = genome.get_regions("HindIII")
        hic.add_regions(regions)

        # load edges
        edges = []
        edges.append(HicEdge(source=10, sink=23, weight=3)
        edges.append(HicEdge(source=8, sink=9, weight=57)
        # ...
        hic.add_edges(edges)
    """

    class HicEdgeDescription(t.IsDescription):
        source = t.Int32Col(pos=0)  
        sink = t.Int32Col(pos=1)  
        weight = t.Float64Col(pos=2)  
    
    def __init__(self, data=None, file_name=None,
                 mode='a',
                 _table_name_nodes='nodes',
                 _table_name_edges='edges'):

        """
        Initialize a :class:`~Hic` object.

        :param data: Can be the path to an XML file denoting a Hic object,
                     another Hic object, a :class:`~FragmentMappedReadPairs`
                     object, or a path to a save file. In the latter case,
                     this parameter may replace file_name, but only if
                     file_name is None.
        :param file_name: Path to a save file
        :param _table_name_nodes: (Internal) name of the HDF5 node for regions
        :param _table_name_edges: (Internal) name of the HDF5 node for edges
        """
        
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
        
        FileBased.__init__(self, file_name, mode=mode)
        RegionsTable.__init__(self, file_name=self.file, _table_name_regions=_table_name_nodes)

        if _table_name_edges in self.file.root:
            self._edges = self.file.get_node('/', _table_name_edges)
        else:
            self._edges = MaskedTable(self.file.root, _table_name_edges,
                                      Hic.HicEdgeDescription)
        
        self._edges.flush()
        
        # generate tables from inherited classes
        Maskable.__init__(self, self.file)
        MetaContainer.__init__(self, self.file)


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
                    
            # data is existing Hic object
            elif isinstance(data, Hic):
                self.load_from_hic(data)
            else:
                try:
                    self.load_read_fragment_pairs(data)
                except AttributeError:
                    raise ValueError("Input data type not recognized")
    
    def __del__(self):
        self.close()
    
    def load_read_fragment_pairs(self, pairs, _max_buffer_size=5000000):
        """
        Load data from :class:`~kaic.construct.seq.FragmentMappedReadPairs`.

        This method automatically sums up reads mapping to the same
        fragment pairs and creates exactly one edge per fragment pair.

        :param pairs: A :class:`~kaic.construct.seq.FragmentMappedReadPairs`
                      object.
        :param _max_buffer_size: Number of edges kept in buffer before
                                 writing to Table.
        """
        # add regions
        if len(self._regions) != 0:
            raise RuntimeError("When importing from read pairs you MUST start from an empty data set!")
        self.add_regions(pairs.regions())

        edge_buffer = {}
        for pair in pairs._pairs:
            source = pair["left_fragment"]
            sink = pair["right_fragment"]
            if source > sink:
                tmp = source
                source = sink
                sink = tmp
            key = (source, sink)
            if key not in edge_buffer:
                edge_buffer[key] = 0
            edge_buffer[key] += 1

            if len(edge_buffer) > _max_buffer_size:
                logging.info("Flushing buffer")
                self._flush_edge_buffer(edge_buffer, replace=False, update_index=False)
                edge_buffer = {}
        logging.info("Final flush")
        self._flush_edge_buffer(edge_buffer, replace=False)

    def load_from_hic(self, hic, _edge_buffer_size=5000000,
                      _edges_by_overlap_method=_edge_overlap_split_rao):
        """
        Load data from another :class:`~Hic` object.

        :param hic: Another :class:`~Hic` object
        :param _edge_buffer_size: Number of edges in memory before writing
                                  to file
        :param _edges_by_overlap_method: A function that maps reads from
                                         one genomic region to others using
                                         a supplied overlap map. By default
                                         it uses the Rao et al. (2014) method.
                                         See :func:`~_edge_overlap_split_rao`
        """
        # if we do not have any nodes in this Hi-C object...
        if len(self.regions()) == 0:
            logging.info("Copying Hi-C")
            # ...simply import everything
            for region in hic.regions():
                self.add_region(region, flush=False)
            for edge in hic.edges():
                self.add_edge(edge, check_nodes_exist=False, flush=False)
            self.flush()
        # if already have nodes in this HiC object...
        else:
            logging.info("Binning Hi-C contacts")
            # create region "overlap map"
            overlap_map = _get_overlap_map(hic.regions(), self.regions())

            edge_buffer = {}
            for old_edge in hic._edges:
                old_source = old_edge['source']
                old_sink = old_edge['sink']
                old_weight = old_edge['weight']
                new_edges = _edges_by_overlap_method([old_source, old_sink, old_weight], overlap_map)

                for new_edge in new_edges:
                    key_pair = (new_edge[0], new_edge[1])
                    if key_pair not in edge_buffer:
                        edge_buffer[key_pair] = 0
                    edge_buffer[key_pair] += new_edge[2]

                if len(edge_buffer) > _edge_buffer_size:
                    self._flush_edge_buffer(edge_buffer, replace=False, update_index=False)
                    edge_buffer = {}
            self._flush_edge_buffer(edge_buffer)

    def bin(self, bin_size, file_name=None):
        """
        Map edges in this object to equi-distant bins.

        :param bin_size: Bin size in base pairs
        :param file_name: File name of the new, binned Hic object
        :return: :class:`~Hic` object
        """
        # find chromosome lengths
        chromosomes = self.chromosomes()
        chromosome_sizes = {chromosome: 0 for chromosome in chromosomes}
        for region in self.regions():
            if chromosome_sizes[region.chromosome] < region.end:
                chromosome_sizes[region.chromosome] = region.end

        chromosome_list = []
        for chromosome in chromosomes:
            chromosome_list.append(Chromosome(name=chromosome,length=chromosome_sizes[chromosome]))

        genome = Genome(chromosomes=chromosome_list)
        hic = Hic(file_name=file_name)
        if len(hic.regions()) > 0:
            # you are loading an existing Hic object from file - this is probably not
            # what you want to do.
            raise RuntimeError("The Hic object that you are trying to bin into already exists "
                               "and has more than 0 regions!")
        hic.add_regions(genome.get_regions(bin_size))

        hic.load_from_hic(self)

        return hic

    def bin_size(self):
        node = self.get_node(0)
        return node.end - node.start + 1

    @classmethod
    def from_hiclib(cls, hl, file_name=None):
        """
        Create :class:`~Hic` object from hiclib object.

        :param hl: hiclib object
        :param file_name: Path to save file
        :return: :class:`~Hic`
        """
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
        """
        Add a :class:`~HicNode` or :class:`~GenomicRegion`.

        :param node: :class:`~HicNode` or :class:`~GenomicRegion`,
                     see :func:`~RegionsTable.add_region` for details
        :param flush: Write data to file immediately after import.
        """
        return self.add_region(node, flush)        
    
    def add_edge(self, edge, check_nodes_exist=True, flush=True):
        """
        Add an edge to this object.

        :param edge: :class:`~HicEdge`, dict with at least the
                     attributes source and sink, optionally weight,
                     or a list of length 2 (source, sink) or 3
                     (source, sink, weight).
        :param check_nodes_exist: Make sure that there are nodes
                                  that match source and sink indexes
        :param flush: Write data to file immediately after import
        """
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
            n_regions = len(self._regions)
            if source >= n_regions or sink >= n_regions:
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
        """
        Bulk-add nodes from a list.

        :param nodes: List (or iterator) of nodes. See
                      :func:`~Hic.add_node` for details.
        """
        self.add_regions(nodes)
    
    def add_edges(self, edges):
        """
        Bulk-add edges from a list.

        :param edges: List (or iterator) of edges. See
                      :func:`~Hic.add_edge` for details
        """
        for edge in edges:
            self.add_edge(edge, flush=False)
        self.flush(flush_nodes=False)

    def merge(self, hic, _edge_buffer_size=5000000):
        """
        Merge this object with another :class:`~Hic` object.

        First merges genomic regions, then merges edges.
        It is strongly advised that the genomic regions in
        both objects are the same, although this method will attempt to
        "translate" regions from one object to the other if
        this is not the case.

        :param hic: :class:`~Hic` object to be merged into this one
        """

        ix_conversion = {}

        # check if regions are identical (saves a lot of time)
        logging.info("Checking if regions are identical")
        identical = True
        region_counter = 0
        for self_region, hic_region in zip(self.regions(), hic.regions()):
            if self_region.chromosome != hic_region.chromosome:
                identical = False
                break
            if self_region.start != hic_region.start:
                identical = False
                break
            if self_region.end != hic_region.end:
                identical = False
                break
            ix_conversion[region_counter] = region_counter
            region_counter += 1

        if region_counter < len(hic.regions()):
            identical = False

        if not identical:
            ix_conversion = {}
            # merge genomic regions
            self.log_info("Merging genomic regions...")
            for region in hic.regions():
                ix = self._get_region_ix(region)
                if ix is None:
                    ix = self.add_region([region.chromosome, region.start, region.end], flush=False)
                ix_conversion[region.ix] = ix
            self._regions.flush()

        # merge edges
        self.log_info("Merging contacts...")
        edge_buffer = {}
        l = len(hic._edges)
        last_percent = 0.0
        for i, merge_row in enumerate(hic._edges):
            merge_source = ix_conversion[merge_row["source"]]
            merge_sink = ix_conversion[merge_row["sink"]]
            merge_weight = merge_row["weight"]

            if merge_source > merge_sink:
                tmp = merge_source
                merge_source = merge_sink
                merge_sink = tmp

            edge_buffer[(merge_source, merge_sink)] = merge_weight

            if i/l > last_percent:
                logging.info("%d%%" % int(round(last_percent*100)))
                last_percent += 0.05

            if len(edge_buffer) > _edge_buffer_size:
                logging.info("Flushing buffer...")
                self._flush_edge_buffer(edge_buffer, replace=False, update_index=False)
                edge_buffer = {}

        # final flush
        self.log_info("Final flush")
        self._flush_edge_buffer(edge_buffer, replace=False)

    def _flush_edge_buffer(self, e_buffer, replace=False, update_index=True):
        # update current rows
        for row in self._edges:
            key = (row["source"], row["sink"])

            if key in e_buffer:
                if replace:
                    row["weight"] = e_buffer[key]
                else:
                    row["weight"] += e_buffer[key]
                row.update()
                del e_buffer[key]
        self._edges.flush()

        # flush remaining buffer
        row = self._edges.row
        for source, sink in e_buffer:
            weight = e_buffer[(source, sink)]
            if weight == 0:
                continue
            row["source"] = source
            row["sink"] = sink
            row["weight"] = weight
            row.append()
        self._edges.flush()
        self._remove_zero_edges(update_index=update_index)

    def flush(self, flush_nodes=True, flush_edges=True, update_index=True):
        """
        Write data to file and flush buffers.

        :param flush_nodes: Flush nodes tables
        :param flush_edges: Flush edges table
        :param update_index: Update mask indices in edges table
        """
        if flush_nodes:
            self._regions.flush()
            # re-indexing not necessary when 'autoindex' is True on table
            if not self._regions.autoindex:
                # reindex node table
                self._regions.flush_rows_to_index()
        if flush_edges:
            self._edges.flush(update_index=update_index)
            if not self._edges.autoindex:
                # reindex edge table
                self._edges.flush_rows_to_index()

    def __getitem__(self, key):
        """
        Get a chunk of the Hi-C matrix.
        
        Possible key types are:

        Region types

        - HicNode: Only the ix of this node will be used for
          identification
        - GenomicRegion: self-explanatory
        - str: key is assumed to describe a genomic region
          of the form: <chromosome>[:<start>-<end>:[<strand>]],
          e.g.: 'chr1:1000-54232:+'

        Node types

        - int: node index
        - slice: node range

        List types

        - list: This key type allows for a combination of all
          of the above key types - the corresponding matrix
          will be concatenated

            
        If the key is a 2-tuple, each entry will be treated as the 
        row and column key, respectively,
        e.g.: 'chr1:0-1000, chr4:2300-3000' will extract the Hi-C
        map of the relevant regions between chromosomes 1 and 4.

        :return: :class:`HicMatrix`
        """
        
        nodes_row, nodes_col = self._get_nodes_from_key(key, as_index=False)

        nodes_ix_row = None
        if nodes_row is not None:
            if isinstance(nodes_row, list):
                nodes_ix_row = [node.ix for node in nodes_row]
            else:
                nodes_ix_row = nodes_row.ix

        nodes_ix_col = None
        if nodes_col is not None:
            if isinstance(nodes_col, list):
                nodes_ix_col = [node.ix for node in nodes_col]
            else:
                nodes_ix_col = nodes_col.ix
        
        m = self._get_matrix(nodes_ix_row, nodes_ix_col)

        # select the correct output format
        # empty result: matrix
        if m.shape[0] == 0 and m.shape[1] == 0:
            return HicMatrix(m, col_regions=[], row_regions=[])
        # both selectors are lists: matrix
        if isinstance(nodes_ix_row, list) and isinstance(nodes_ix_col, list):
            return HicMatrix(m, col_regions=nodes_col, row_regions=nodes_row)
            #return m
        # row selector is list: vector
        if isinstance(nodes_ix_row, list):
            return HicMatrix(m[:, 0], col_regions=[nodes_ix_col], row_regions=nodes_row)
        # column selector is list: vector
        if isinstance(nodes_ix_col, list):
            return HicMatrix(m[0, :], col_regions=nodes_col, row_regions=[nodes_row])
        # both must be indexes
        return m[0, 0]
    
    def _get_nodes_from_key(self, key, as_index=False):
        if isinstance(key, tuple):
            nodes_ix_row = self._getitem_nodes(key[0], as_index=as_index)
            nodes_ix_col = self._getitem_nodes(key[1], as_index=as_index)
        else:
            nodes_ix_row = self._getitem_nodes(key, as_index=as_index)
            nodes_ix_col = []
            for region in self.regions():
                if as_index:
                    nodes_ix_col.append(region.ix)
                else:
                    nodes_ix_col.append(region)
        
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
                        m[ir + row_offset, jr + col_offset] = weight
                    if (row_range[0] <= sink <= row_range[1]
                        and col_range[0] <= source <= col_range[1]):
                        ir = sink - row_range[0]
                        jr = source - col_range[0]
                        m[ir + row_offset, jr + col_offset] = weight
                
                col_offset += n_cols_sub
            row_offset += n_rows_sub

        return m
    
    def _getitem_nodes(self, key, as_index=False):
        # 'chr1:1234:56789'
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
                region_nodes = [RegionsTable._row_to_region(row) for row in self._regions.where(condition)]
            
            return region_nodes
        
        # 1:453
        if isinstance(key, slice):
            if as_index:
                return [row['ix'] for row in self._regions.iterrows(key.start, key.stop, key.step)]
            else:
                return [RegionsTable._row_to_region(row) for row in self._regions.iterrows(key.start, key.stop, key.step)]
        
        # 432
        if isinstance(key, int):
            row = self._regions[key]
            if as_index:
                return row['ix']
            else:
                return RegionsTable._row_to_region(row)
        
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
        """
        Get a pandas data frame by key.

        For key types see :func:`~Hic.__getitem__`.

        :param key: For key types see :func:`~Hic.__getitem__`.
        :return: Pandas data frame, row and column labels are
                 corresponding node start positions
        """
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

        - HicNode: Only the ix of this node will be used for
          identification
        - GenomicRegion: self-explanatory
        - str: key is assumed to describe a genomic region
          of the form: <chromosome>[:<start>-<end>:[<strand>]],
          e.g.: 'chr1:1000-54232:+'

        Node types

        - int: node index
        - slice: node range

        List types

        - list: This key type allows for a combination of all
          of the above key types - the corresponding matrix
          will be concatenated

        If the key is a 2-tuple, each entry will be treated as the 
        row and column key, respectively,
        e.g.: 'chr1:0-1000, chr4:2300-3000' will set the Hi-C
        map of the relevant regions between chromosomes 1 and 4.
            
        """
        
        nodes_ix_row, nodes_ix_col = self._get_nodes_from_key(key, as_index=True)
        self._set_matrix(item, nodes_ix_row, nodes_ix_col)

    def _set_matrix(self, item, nodes_ix_row=None, nodes_ix_col=None):
        replacement_edges = {}

        # create new edges with updated weights
        # select the correct format:
        def swap(old_source, old_sink):
            if old_source > old_sink:
                return old_sink, old_source
            return old_source, old_sink

        # both selectors are lists: matrix
        if isinstance(nodes_ix_row, list) and isinstance(nodes_ix_col, list):
            n_rows = len(nodes_ix_row)
            n_cols = len(nodes_ix_col)
            # check that we have a matrix with the correct dimensions
            if not isinstance(item, np.ndarray) or not np.array_equal(item.shape, [n_rows, n_cols]):
                raise ValueError("Item is not a numpy array with shape (%d,%d)!" % (n_rows, n_cols))
            
            for i in xrange(0, n_rows):
                for j in xrange(0, n_cols):
                    source = nodes_ix_row[i]
                    sink = nodes_ix_col[j]
                    source, sink = swap(source, sink)
                    weight = item[i, j]
                    key = (source, sink)
                    if key not in replacement_edges:
                        replacement_edges[key] = weight

        # row selector is list: vector
        elif isinstance(nodes_ix_row, list):
            n_rows = len(nodes_ix_row)
            if not isinstance(item, np.ndarray) or not np.array_equal(item.shape, [n_rows]):
                raise ValueError("Item is not a numpy vector of length %d!" % n_rows)
            
            for i, my_sink in enumerate(nodes_ix_row):
                source = nodes_ix_col
                source, sink = swap(source, my_sink)
                weight = item[i]
                key = (source, sink)
                if key not in replacement_edges:
                    replacement_edges[key] = weight
        
        # column selector is list: vector
        elif isinstance(nodes_ix_col, list):
            n_cols = len(nodes_ix_col)
            if not isinstance(item, np.ndarray) or not np.array_equal(item.shape, [n_cols]):
                raise ValueError("Item is not a numpy vector of length %d!" % n_cols)
            
            for i, my_source in enumerate(nodes_ix_col):
                sink = nodes_ix_row
                source, sink = swap(my_source, sink)
                weight = item[i]
                key = (source, sink)
                if key not in replacement_edges:
                    replacement_edges[key] = weight

        # both must be indexes
        else:
            weight = item
            source, sink = swap(nodes_ix_row, nodes_ix_col)
            key = (source, sink)
            if key not in replacement_edges:
                replacement_edges[key] = weight

        self._flush_edge_buffer(replacement_edges, replace=True)

    def _set_matrix_old(self, item, nodes_ix_row=None, nodes_ix_col=None):
        # calculate number of rows
        if (nodes_ix_row is not None
            and not isinstance(nodes_ix_row, list)):
            range_nodes_ix_row = [nodes_ix_row]
        else:
            range_nodes_ix_row = nodes_ix_row

        # calculate number of columns
        if (nodes_ix_col is not None
            and not isinstance(nodes_ix_col, list)):
            range_nodes_ix_col = [nodes_ix_col]
        else:
            range_nodes_ix_col = nodes_ix_col

        # get row range generator
        row_ranges = ranges(range_nodes_ix_row)

        # set every edge that is to be replaced to 0
        row_offset = 0
        for row_range in row_ranges:
            n_rows_sub = row_range[1] - row_range[0] + 1

            col_ranges = ranges(range_nodes_ix_col)
            col_offset = 0
            for col_range in col_ranges:
                n_cols_sub = col_range[1] - col_range[0] + 1

                condition = "((source >= %d) & (source <= %d)) & ((sink >= %d) & (sink <= %d))"
                condition += "| ((source >= %d) & (source <= %d)) & ((sink >= %d) & (sink <= %d))"
                condition = condition % (row_range[0], row_range[1], col_range[0], col_range[1],
                                         col_range[0], col_range[1], row_range[0], row_range[1])

                # actually set weight to zero
                for edge_row in self._edges.where(condition):
                    edge_row['weight'] = 0
                    edge_row.update()

                col_offset += n_cols_sub
            row_offset += n_rows_sub

        self.flush()
        self._remove_zero_edges()

        # create new edges with updated weights
        # select the correct format
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
                    self.add_edge([source, sink, weight], flush=False)

        # row selector is list: vector
        elif isinstance(nodes_ix_row, list):
            n_rows = len(nodes_ix_row)
            if (not isinstance(item, np.ndarray) or
                not np.array_equal(item.shape, [n_rows])):
                raise ValueError("Item is not a numpy vector of length %d!" % (n_rows))

            for i, sink in enumerate(nodes_ix_row):
                source = nodes_ix_col
                weight = item[i]
                self.add_edge([source, sink, weight], flush=False)

        # column selector is list: vector
        elif isinstance(nodes_ix_col, list):
            n_cols = len(nodes_ix_col)
            if (not isinstance(item, np.ndarray) or
                not np.array_equal(item.shape, [n_cols])):
                raise ValueError("Item is not a numpy vector of length %d!" % (n_cols))

            for i, source in enumerate(nodes_ix_col):
                sink = nodes_ix_row
                weight = item[i]
                self.add_edge([source, sink, weight], flush=False)

        # both must be indexes
        else:
            weight = item
            self.add_edge([nodes_ix_row, nodes_ix_col, weight], flush=False)

        self.flush()
    
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
                self.flush()
        if not value_set:
            self.add_edge(HicEdge(source=source, sink=sink, weight=weight), flush=flush)
    
    def _remove_zero_edges(self, flush=True, update_index=True):
        zero_edge_ix = []
        ix = 0
        for row in self._edges.iterrows():
            if row['weight'] == 0:
                zero_edge_ix.append(ix)
            ix += 1
        
        for ix in reversed(zero_edge_ix):
            self._edges.remove_row(ix)        
        
        if flush:
            self.flush(update_index=update_index)
    
    def autoindex(self, index=None):
        """
        Switch on/off autoindexing.
        """
        if index is not None:
            self._regions.autoindex = bool(index)
            self._edges.autoindex = bool(index)
            return index
        return self._regions.autoindex

    def save(self, file_name, _table_name_nodes='nodes', _table_name_edges='edges',
             _table_name_meta='meta', _table_name_meta_values='meta', _table_name_mask='mask'):
        """
        Copy content of this object to a new file.

        :param file_name: Path to new save file
        """
        self.file.copy_file(file_name)
        self.file.close()
        self.file = create_or_open_pytables_file(file_name)
        self._regions = self.file.get_node('/' + _table_name_nodes)
        self._edges = self.file.get_node('/' + _table_name_edges)
        self._meta = self.file.get_node('/' + _table_name_meta)
        self._meta_values = self.file.get_node('/' + _table_name_meta_values)
        self._mask = self.file.get_node('/' + _table_name_mask)

    def _row_to_node(self, row, lazy=False):
        if lazy:
            return LazyHicNode(row)
        return HicNode(chromosome=row["chromosome"], start=row["start"],
                       end=row["end"], ix=row["ix"])

    def get_node(self, key):
        """
        Get a single node by key.

        :param key: For possible key types see :func:`~Hic.__getitem__`
        :return: A :class:`~HicNode` matching key
        """
        found_nodes = self.get_nodes(key)
        if isinstance(found_nodes, list):
            if len(found_nodes) > 1:
                raise IndexError("More than one node found matching %s" % str(key))
            if len(found_nodes) == 1:
                return found_nodes[0]
        return found_nodes

    def get_nodes(self, key):
        """
        Get multiple nodes by key.

        :param key: For possible key types see :func:`~Hic.__getitem__`
        :return: A list of :class:`~HicNode` objects matching key
        """
        return self._getitem_nodes(key)

    def _row_to_edge(self, row, lazy=False):
        if not lazy:
            source = row["source"]
            sink = row["sink"]
            weight = row["weight"]
            source_node_row = self._regions[source]
            source_node = self._row_to_node(source_node_row)
            sink_node_row = self._regions[sink]
            sink_node = self._row_to_node(sink_node_row)
            return HicEdge(source_node, sink_node, weight)
        else:
            return LazyHicEdge(row, self._regions)

    def get_edge(self, ix, lazy=False):
        """
        Get an edge from this object's edge list.

        :param ix: integer
        :return:
        """
        return self._row_to_edge(self._edges[ix], lazy=lazy)
    
    def nodes(self):
        """
        Iterator over this object's nodes/regions.

        See :func:`~RegionsTable.regions` for details.
        :return: Iterator over :class:`~GenomicRegions`
        """
        return self.regions()
    
    def edges(self, lazy=False):
        """
        Iterate over :class:`~HicEdge` objects.

        :return: Iterator over :class:`~HicEdge`
        """
        hic = self

        class EdgeIter:
            def __init__(self):
                self.iter = iter(hic._edges)
                
            def __iter__(self):
                return self
            
            def next(self):
                return hic._row_to_edge(self.iter.next(), lazy=lazy)

            def __len__(self):
                return len(hic._edges)
        return EdgeIter()

    def filter(self, edge_filter, queue=False, log_progress=False):
        """
        Filter edges in this object by using a
        :class:`~HicEdgeFilter`.

        :param edge_filter: Class implementing :class:`~HicEdgeFilter`.
                            Must override valid_edge method, ideally sets mask parameter
                            during initialization.
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        :param log_progress: If true, process iterating through all edges
                             will be continuously reported.
        """
        edge_filter.set_hic_object(self)
        if not queue:
            self._edges.filter(edge_filter, _logging=log_progress)
        else:
            self._edges.queue_filter(edge_filter)

    def run_queued_filters(self, log_progress=False):
        """
        Run queued filters.

        :param log_progress: If true, process iterating through all edges
                             will be continuously reported.
        """
        self._edges.run_queued_filters(_logging=log_progress)

    def filter_diagonal(self, distance=0, queue=False):
        """
        Convenience function that applies a :class:`~DiagonalFilter`.

        :param distance: Distance from the diagonal up to which matrix entries
                         will be filtered. The default, 0, filters only the
                         diagonal itself.
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        """
        mask = self.add_mask_description('diagonal',
                                         'Mask the diagonal of the Hic matrix (up to distance %d)' % distance)
        diagonal_filter = DiagonalFilter(distance=distance, mask=mask)
        self.filter(diagonal_filter, queue)

    def filter_low_coverage_regions(self, cutoff=None, queue=False):
        """
        Convenience function that applies a :class:`~LowCoverageFilter`.

        :param cutoff: Cutoff (contact count, float) below which a region
                       is considered to have low coverage. If not set
                       explicitly, defaults to 5% of the mean region coverage.
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        """
        if cutoff is not None:
            mask = self.add_mask_description('low_coverage',
                                             'Mask low coverage regions in the Hic matrix (cutoff %.4f)' % cutoff)
        else:
            mask = self.add_mask_description('low_coverage',
                                             'Mask low coverage regions in the Hic matrix (10%)')
        low_coverage_filter = LowCoverageFilter(self, cutoff=cutoff, mask=mask)
        self.filter(low_coverage_filter, queue)
    
    def bias_vector(self, vector=None):
        """
        Get the bias vector of this Hic matrix.

        Only works if previously corrected.
        """
        if vector is not None:
            self._edges._v_attrs.bias_vector = vector
        return self._edges._v_attrs.bias_vector

    def marginals(self):
        """
        Get the marginals vector of this Hic matrix.
        """
        # prepare marginals dict
        marginals = np.zeros(len(self.regions()), float)

        for edge in self.edges(lazy=True):
            marginals[edge.source] += edge.weight
            if edge.source != edge.sink:
                marginals[edge.sink] += edge.weight

        return marginals

    def _get_boundary_distances(self):
        n_bins = len(self.regions())
        # find distances to chromosome boundaries in bins
        boundary_dist = np.zeros(n_bins, dtype=int)
        last_chromosome = None
        last_chromosome_index = 0
        for i, node in enumerate(self.nodes()):
            chromosome = node.chromosome
            if last_chromosome is not None and chromosome != last_chromosome:
                chromosome_length = i-last_chromosome_index
                for j in xrange(chromosome_length):
                    boundary_dist[last_chromosome_index+j] = min(j, i-last_chromosome_index-1-j)
                last_chromosome_index = i
            last_chromosome = chromosome
        chromosome_length = n_bins-last_chromosome_index
        for j in xrange(chromosome_length):
            boundary_dist[last_chromosome_index+j] = min(j, n_bins-last_chromosome_index-1-j)

        return boundary_dist

    def directionality_index(self, window_size=2000000):
        """
        Calculate the directionality index according to Dixon 2012 for every bin.

        :param window_size: size of the sliding window in base pairs
        :return: numpy array same length as bins in the object
        """
        bin_size = self.bin_size()
        bin_window_size = int(window_size/bin_size)
        if window_size % bin_size > 0:
            bin_window_size += 1

        n_bins = len(self.regions())
        boundary_dist = self._get_boundary_distances()

        left_sums = np.zeros(n_bins)
        right_sums = np.zeros(n_bins)
        directionality_index = np.zeros(n_bins)
        for edge in self.edges(lazy=True):
            source = edge.source
            sink = edge.sink
            weight = edge.weight
            if source == sink:
                continue
            if sink - source <= bin_window_size:
                if boundary_dist[sink] >= sink-source:
                    left_sums[sink] += weight
                if boundary_dist[source] >= sink-source:
                    right_sums[source] += weight

        for i in xrange(n_bins):
            A = left_sums[i]
            B = right_sums[i]
            E = (A+B)/2
            if E != 0 and B-A != 0:
                directionality_index[i] = ((B-A)/abs(B-A)) * ((((A-E)**2)/E) + (((B-E)**2)/E))

        return directionality_index


class HicEdgeFilter(MaskFilter):
    """
    Abstract class that provides filtering functionality for the
    edges/contacts in a :class:`~Hic` object.

    Extends MaskFilter and overrides valid(self, row) to make
    :class:`~HicEdge` filtering more "natural".

    To create custom filters for the :class:`~Hic` object, extend this
    class and override the valid_edge(self, edge) method.
    valid_edge should return False for a specific :class:`~HicEdge` object
    if the object is supposed to be filtered/masked and True
    otherwise. See :class:`~DiagonalFilter` for an example.

    Pass a custom filter to the :func:`~Hic.filter` method in :class:`~Hic`
    to apply it.
    """

    __metaclass__ = ABCMeta

    def __init__(self, mask=None):
        """
        Initialize HicEdgeFilter.

        :param mask: The Mask object that should be used to mask
                     filtered :class:`~HicEdge` objects. If None the default
                     Mask will be used.
        """
        super(HicEdgeFilter, self).__init__(mask)
        self._hic = None

    @abstractmethod
    def valid_edge(self, edge):
        """
        Determine if a :class:`~HicEdge` object is valid or should
        be filtered.

        When implementing custom HicEdgeFilter this method must be
        overridden. It should return False for :class:`~HicEdge` objects that
        are to be fitered and True otherwise.

        Internally, the :class:`~Hic` object will iterate over all HicEdge
        instances to determine their validity on an individual
        basis.

        :param edge: A :class:`~HicEdge` object
        :return: True if :class:`~HicEdge` is valid, False otherwise
        """
        pass

    def set_hic_object(self, hic_object):
        """
        Set the :class:`~Hic` instance to be filtered by this
        HicEdgeFilter.

        Used internally by :class:`~Hic` instance.

        :param hic_object: :class:`~Hic` object
        """
        self._hic = hic_object

    def valid(self, row):
        """
        Map valid_edge to MaskFilter.valid(self, row).

        :param row: A pytables Table row.
        :return: The boolean value returned by valid_edge.
        """
        edge = self._hic._row_to_edge(row, lazy=True)
        return self.valid_edge(edge)


class DiagonalFilter(HicEdgeFilter):
    """
    Filter contacts in the diagonal of a :class:`~Hic` matrix.
    """
    def __init__(self, distance=0, mask=None):
        """
        Initialize filter with chosen parameters.

        :param distance: Distance from the diagonal up to which
                         contacts will be filtered
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered edges.
        """
        HicEdgeFilter.__init__(self, mask=mask)
        self.distance = distance

    def valid_edge(self, edge):
        """
        Check if an edge is on (or near) the diagonal of the :class:`~Hic` matrix.
        """
        if abs(edge.source-edge.sink) <= self.distance:
            return False
        return True


class LowCoverageFilter(HicEdgeFilter):
    """
    Filter a :class:`~HicEdge` if it connects a region that
    does not have a contact count larger than a specified
    cutoff.

    If the cutoff is not provided, it is automatically
    chosen at 5% of the mean contact count of all regions.
    """
    def __init__(self, hic_object, cutoff=None, mask=None):
        """
        Initialize filter with these settings.

        :param hic_object: The :class:`~Hic` object that this
                           filter will be called on. Needed for
                           contact count calculation.
        :param cutoff: A cutoff in contacts (can be float) below
                       which regions are considered "low coverage"
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered edges.
        """
        HicEdgeFilter.__init__(self, mask=mask)

        self._marginals = hic_object.marginals()
        if cutoff is None:
            cutoff, _ = self.calculate_cutoffs(0.1)

        self._regions_to_mask = set()
        for i, contacts in enumerate(self._marginals):
            if contacts < cutoff:
                self._regions_to_mask.add(i)

    def calculate_cutoffs(self, fraction_threshold=0.05):
        lower = np.median(self._marginals[self._marginals > 0])*fraction_threshold
        upper = np.median(self._marginals[self._marginals > 0])+lower
        return lower, upper

    def valid_edge(self, edge):
        """
        Check if an edge falls into a low-coverage region.
        """
        if edge.source in self._regions_to_mask:
            return False
        if edge.sink in self._regions_to_mask:
            return False
        return True


class HicMatrix(np.ndarray):
    def __new__(cls, input_matrix, col_regions=None, row_regions=None):
        obj = np.asarray(input_matrix).view(cls)
        obj.col_regions = col_regions
        obj.row_regions = row_regions
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return

        self.row_regions = getattr(obj, 'row_regions', None)
        self.col_regions = getattr(obj, 'col_regions', None)

    def __getitem__(self, index):
        self._getitem = True

        # convert string types into region indexes
        if isinstance(index, tuple):
            row_key = self._convert_key(index[0], self.row_regions)
            col_key = self._convert_key(index[1], self.col_regions)
            index = (row_key, col_key)
        else:
            row_key = self._convert_key(index, self.row_regions)
            col_key = slice(0, len(self.col_regions), 1)
            index = row_key

        try:
            out = np.ndarray.__getitem__(self, index)
        finally:
            self._getitem = False

        if not isinstance(out, np.ndarray):
            return out

        # get regions
        try:
            row_regions = self.row_regions[row_key]
        except TypeError:
            row_regions = None
            #logging.warn("Key type %s cannot yet be handeled by HicMatrix." % str(row_key) +
            #             "Falling back on setting row regions to None")

        try:
            col_regions = self.col_regions[col_key]
        except TypeError:
            col_regions = None
            #logging.warn("Key type %s cannot yet be handeled by HicMatrix." % str(col_key) +
            #             "Falling back on setting col regions to None")

        out.col_regions = col_regions
        out.row_regions = row_regions

        return out

    def __getslice__(self, start, stop):
        return self.__getitem__(slice(start, stop))

    def _convert_key(self, key, regions):
        if isinstance(key, str):
            key = GenomicRegion.from_string(key)
        if isinstance(key, GenomicRegion):
            key_start = max(0, key.start)
            key_end = key.end
            start = None
            stop = None
            for i, region in enumerate(regions):
                if region.chromosome == key.chromosome:
                    if (key_end is None or region.start <= key_end) and region.end >= key_start:
                        if start is None:
                            start = i
                        stop = i
            return slice(start, stop, 1)
        return key


class HicXmlFile(object):
    def __init__(self, file_name):
        self.file_name = file_name

    def nodes(self):
        file_name = self.file_name
        
        class XmlNodeIter:
            def __init__(self):
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
            
        return XmlNodeIter()
    
    def edges(self):
        file_name = self.file_name
        
        class XmlEdgeIter:
            def __init__(self):
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
            
        return XmlEdgeIter()


def genome_from_string(genome_string):
    """
    Convenience function to load a :class:`~Genome` from a string.

    :param genome_string: Path to FASTA file, path to folder with
                          FASTA files, comma-separated list of
                          paths to FASTA files, path to HDF5 file
    :return: A :class:`~Genome` object
    """
    genome = None
    # case 1: FASTA file = Chromosome
    if is_fasta_file(genome_string):
        chromosome = Chromosome.from_fasta(genome_string)
        genome = Genome(chromosomes=[chromosome])
    # case 2: Folder with FASTA files
    elif os.path.isdir(genome_string):
        genome = Genome.from_folder(genome_string)
    # case 3: path to HDF5 file
    elif is_hdf5_file(genome_string):
        genome = Genome(genome_string)
    # case 4: List of FASTA files
    else:
        chromosome_files = genome_string.split(',')
        chromosomes = []
        for chromosome_file in chromosome_files:
            chromosome = Chromosome.from_fasta(os.path.expanduser(chromosome_file))
            chromosomes.append(chromosome)
        genome = Genome(chromosomes=chromosomes)
    
    return genome


def _get_overlap_map(old_regions, new_regions):
    # 1. organize regions in self by chromosome
    new_region_map = {}
    for i, new_region in enumerate(new_regions):
        if not new_region.chromosome in new_region_map:
            new_region_map[new_region.chromosome] = []
        new_region_map[new_region.chromosome].append([new_region.start,new_region.end,i])
        
    # 2. iterate over regions in hic to find overlap
    def _get_overlap(new_region, old_region):
        new_region_length = new_region[1] - new_region[0] + 1
        overlap = min(old_region[1], new_region[1]) - max(old_region[0],new_region[0]) + 1
        return max(0,overlap/new_region_length)
        
    old_to_new = {}
    current_chromosome = ''
    current_ix = 0
    for i, old_region in enumerate(old_regions):
        old_to_new[i] = []
        if current_chromosome != old_region.chromosome:
            current_ix = 0
            current_chromosome = old_region.chromosome
        
        found_overlap = True
        while found_overlap:
            found_overlap = False
            if current_ix < len(new_region_map[current_chromosome]):
                new_region = new_region_map[current_chromosome][current_ix]
                overlap = _get_overlap(new_region, [old_region.start,old_region.end,i])
                if overlap > 0:
                    old_to_new[i].append([new_region[2], overlap])
                    current_ix += 1
                    found_overlap = True
                elif old_region.start > new_region[1]:
                    current_ix += 1
                    found_overlap = True
            
        current_ix -= 1
    
    return old_to_new
