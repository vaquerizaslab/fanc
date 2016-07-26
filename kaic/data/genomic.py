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

from __future__ import division, print_function
import tables as t
import pandas as p
import numpy as np
import pybedtools
from kaic.tools.files import create_or_open_pytables_file, is_hic_xml_file,\
    is_fasta_file, is_hdf5_file
from kaic.tools.files import is_bed_file, is_bedpe_file
from Bio import SeqIO, Restriction, Seq
from kaic.data.general import Table, TableRow, TableArray, TableObject,\
    MetaContainer, Maskable, MaskedTable, FileBased, MaskFilter, FileGroup
from abc import abstractmethod, ABCMeta
import os.path
import logging
from kaic.tools.general import ranges, distribute_integer, create_col_index
from itertools import izip as zip
from xml.etree import ElementTree as et
import pickle
from collections import defaultdict
import copy
from kaic.tools.general import RareUpdateProgressBar
from kaic.tools.general import range_overlap
from bisect import bisect_right
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


class Bed(pybedtools.BedTool):
    """
    Data type representing a BED file.

    Only exists to support 'with' statements
    """

    def __init__(self, *args, **kwargs):
        pybedtools.BedTool.__init__(self, *args, **kwargs)
    
    def __exit__(self, exec_type, exec_val, exec_tb):
        pass

    def __enter__(self):
        return self

    @property
    def regions(self):
        class RegionIter(object):
            def __init__(self, bed):
                self.bed = bed

            def __iter__(self):
                for region in self.bed:
                    score = float(region.score) if region.score != "." else None
                    gr = GenomicRegion(chromosome=region.chrom, start=region.start, end=region.end,
                                       strand=region.strand, score=score, fields=region.fields)
                    yield gr

            def __call__(self):
                return iter(self)

            def __len__(self):
                return len(self.bed)

        return RegionIter(self)


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

    def close(self):
        self.file.close()

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

        regions = RegionsTable(file_name=file_name)
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
            regions.add_region(region, flush=True)

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

    def __init__(self, start, end=None, chromosome=None, strand=None, ix=None, **kwargs):
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
        if end is None:
            end = start
        self.end = end
        if strand == "+":
            strand = 1
        elif strand == "-":
            strand = -1
        elif strand == "0" or strand == ".":
            strand = None
        self.strand = strand
        self.chromosome = chromosome
        self.ix = ix

        for name, value in kwargs.iteritems():
            setattr(self, name, value)

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

                if not end >= start:
                    raise ValueError("The end coordinate must be bigger than the start.")

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

        if region.start <= self.end or self.end is None or region.start is None:
            if region.end >= self.start or self.start is None or region.end is None:
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

    def _equals(self, region):
        if region.chromosome != self.chromosome:
            return False
        if region.start != self.start:
            return False
        if region.end != self.end:
            return False
        return True

    def __eq__(self, other):
        return self._equals(other)

    def __ne__(self, other):
        return not self._equals(other)

    def __len__(self):
        return self.end - self.start


class BedElement(GenomicRegion):
    def __init__(self, chromosome, start, end, **kwargs):
        super(BedElement, self).__init__(start, end, chromosome)
        for key, value in kwargs.iteritems():
            setattr(self, key, value)


class LazyGenomicRegion(GenomicRegion):
    def __init__(self, row, ix=None, auto_update=True):
        self.reserved = {'_row', 'static_ix', 'strand', 'auto_update'}
        self._row = row
        self.static_ix = ix
        self.auto_update = auto_update

    def __getattr__(self, item):
        if item == 'reserved' or item in self.reserved:
            return object.__getattribute__(self, item)
        try:
            return self._row[item]
        except KeyError:
            raise AttributeError

    def __setattr__(self, key, value):
        if key == 'reserved' or key in self.reserved:
            super(LazyGenomicRegion, self).__setattr__(key, value)
        else:
            self._row[key] = value
            if self.auto_update:
                self.update()

    def update(self):
        self._row.update()

    @property
    def strand(self):
        try:
            return self._row["strand"]
        except KeyError:
            return None

    @property
    def ix(self):
        if self.static_ix is None:
            return self._row["ix"]
        return self.static_ix


class GenomicRegions(object):

    def __init__(self, regions=None):
        self._regions = []
        self._max_region_ix = -1

        if regions is not None:
            for region in regions:
                self.add_region(region)

    def add_region(self, region):
        """
        Add a genomic region to this object.

        This method offers some flexibility in the types of objects
        that can be loaded. See below for details.

        :param region: Can be a :class:`~GenomicRegion`, a dict with
                       at least the fields 'chromosome', 'start', and
                       'end', optionally 'ix', or a list of length 3
                       (chromosome, start, end) or 4 (ix, chromosome,
                       start, end).
        """
        ix = -1

        if isinstance(region, GenomicRegion):
            return self._add_region(copy.copy(region))
        elif isinstance(region, str):
            return self._add_region(GenomicRegion.from_string(region))
        elif type(region) is dict:
            return self._add_region(GenomicRegion(**copy.copy(region)))
        else:
            try:
                offset = 0
                if len(region) == 4:
                    ix = region[0]
                    offset += 1
                chromosome = region[offset]
                start = region[offset + 1]
                end = region[offset + 2]
                strand = 1
            except TypeError:
                raise ValueError("Node parameter has to be GenomicRegion, dict, or list")

        new_region = GenomicRegion(chromosome=chromosome, start=start, end=end, strand=strand, ix=ix)
        return self._add_region(new_region)

    def _add_region(self, region):
        region.ix = self._max_region_ix + 1

        self._regions.append(region)

        if region.ix > self._max_region_ix:
            self._max_region_ix = region.ix

        return self._len()

    def _len(self):
        return len(self._regions)

    def __len__(self):
        return self._len()

    def _get_regions(self, key):
        return self._regions[key]

    @property
    def regions(self):
        """
        Iterate over genomic regions in this object.

        Will return a :class:`~GenomicRegion` object in every iteration.
        Can also be used to get the number of regions by calling
        len() on the object returned by this method.

        :return: Iterator over requested :class:`~GenomicRegion` objects
        """

        this = self

        class RegionIter:
            def __init__(self):
                self.regions = this._regions
                self.iter = iter(self.regions)

            def __iter__(self):
                return self

            def next(self):
                return self.iter.next()

            def __len__(self):
                return this._len()

            def __getitem__(self, key):
                return this._get_regions(key)

            def __call__(self):
                return this.regions

        return RegionIter()

    def __iter__(self):
        return self.regions

    def __getitem__(self, item):
        return self.regions[item]

    def region_bins(self, region):
        """
        Takes a genomic region and returns a slice of the bin
        indices that are covered by the region.

        :param region: String or class:`~GenomicRegion`
                       object for which covered bins will
                       be returned.
        :return: slice
        """
        if isinstance(region, basestring):
            region = GenomicRegion.from_string(region)
        start_ix = None
        end_ix = None
        for r in self.regions:
            if not (r.chromosome == region.chromosome and r.start <= region.end and r.end >= region.start):
                continue
            if start_ix is None:
                start_ix = r.ix
                end_ix = r.ix + 1
                continue
            end_ix = r.ix + 1
        return slice(start_ix, end_ix)

    def intersect(self, region):
        """
        Takes a class:`~GenomicRegion` and returns all regions that
        overlap with the supplied region.

        :param region: String or class:`~GenomicRegion`
                       object for which covered bins will
                       be returned.
        """
        return self.regions[self.region_bins(region)]

    def chromosomes(self):
        """
        Get a list of chromosome names.
        """
        chromosomes_set = set()
        chromosomes = []
        for region in self.regions():
            if region.chromosome not in chromosomes_set:
                chromosomes_set.add(region.chromosome)
                chromosomes.append(region.chromosome)
        return chromosomes

    @property
    def chromosome_lens(self):
        """
        Returns a dictionary of chromosomes and their length
        in bp.
        """
        chr_lens = {}
        for r in self.regions:
            if chr_lens.get(r.chromosome) is None:
                chr_lens[r.chromosome] = r.end
                continue
            if r.end > chr_lens[r.chromosome]:
                chr_lens[r.chromosome] = r.end
        return chr_lens

    @property
    def chromosome_bins(self):
        """
        Returns a dictionary of chromosomes and the start
        and end index of the bins they cover.

        Returned list is xrange-compatible, i.e. chromosome
        bins [0,5] cover chromosomes 1, 2, 3, and 4, not 5.
        """
        chr_bins = {}
        for r in self.regions:
            if chr_bins.get(r.chromosome) is None:
                chr_bins[r.chromosome] = [r.ix, r.ix + 1]
                continue
            chr_bins[r.chromosome][1] = r.ix + 1
        return chr_bins

    def range(self, range_region):
        regions = []

        for region in self.regions:
            if not range_region.chromosome == region.chromosome:
                if len(regions) == 0:
                    continue
                break

            if ((range_region.end is None or region.start <= range_region.end) and
                    (range_region.start is None and region.end >= range_region.start)):
                regions.append(region)

        return regions

    def to_bed(self, file):
        """
        Export regions as BED file
        """
        with open(file, 'w') as f:
            for i, r in enumerate(self.regions):
                print(r.chromosome, r.start - 1, r.end, i, sep="\t", file=f)

    @property
    def regions_dict(self):
        regions_dict = dict()
        for r in self.regions:
            regions_dict[r.ix] = r
        return regions_dict

    @property
    def bin_size(self):
        node = self.regions[0]
        return node.end - node.start + 1

    def distance_to_bins(self, distance):
        bin_size = self.bin_size
        bin_distance = int(distance/bin_size)
        if distance % bin_size > 0:
            bin_distance += 1
        return bin_distance

    def bins_to_distance(self, bins):
        return self.bin_size*bins


class RegionsTable(GenomicRegions, FileGroup):
    """
    PyTables Table wrapper for storing genomic regions.

    This class is inherited by objects working with lists of genomic
    regions, such as equi-distant bins along chromosomes in a genome
    (:class:`~Hic`) or restriction fragments of genomic DNA
    (:class:`~kaic.construct.seq.FragmentMappedReadPairs`)
    """

    _classid = 'REGIONSTABLE'

    class RegionDescription(t.IsDescription):
        """
        Description of a genomic region for PyTables Table
        """
        ix = t.Int32Col(pos=0)
        chromosome = t.StringCol(50, pos=1)
        start = t.Int64Col(pos=2)
        end = t.Int64Col(pos=3)
        strand = t.Int8Col(pos=4)

    def __init__(self, regions=None, file_name=None, mode='a',
                 additional_fields=None, _table_name_regions='regions',
                 tmpdir=None):
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
        if regions is not None:
            # data is file name
            if type(regions) is str or isinstance(regions, t.file.File):
                if file_name is None:
                    file_name = regions
                    regions = None

        try:
            FileGroup.__init__(self, _table_name_regions, file_name, mode=mode, tmpdir=tmpdir)
        except TypeError:
            logging.warn("RegionsTable is now a FileGroup-based object and "
                         "this object will no longer be compatible in the future")

        # check if this is an existing regions file
        try:
            group = self.file.get_node('/', _table_name_regions)

            if isinstance(group, t.table.Table):
                self._regions = group
            else:
                self._regions = self._group.regions

            if len(self._regions) > 0:
                self._max_region_ix = max(row['ix'] for row in self._regions.iterrows())
            else:
                self._max_region_ix = -1
        except t.NoSuchNodeError:
            basic_fields = RegionsTable.RegionDescription().columns.copy()
            if additional_fields is not None:
                if not isinstance(additional_fields, dict) and issubclass(additional_fields, t.IsDescription):
                    # IsDescription subclass case
                    additional_fields = additional_fields.columns

                current = len(basic_fields)
                for key, value in sorted(additional_fields.iteritems(), key=lambda x: x[1]._v_pos):
                    if key not in basic_fields:
                        if value._v_pos is not None:
                            value._v_pos = current
                            current += 1
                        basic_fields[key] = value
            self._regions = t.Table(self._group, 'regions', basic_fields)
            self._max_region_ix = -1

        # index regions table
        create_col_index(self._regions.cols.ix)
        create_col_index(self._regions.cols.start)
        create_col_index(self._regions.cols.end)

        self._ix_to_chromosome = dict()
        self._chromosome_to_ix = dict()

        if regions is not None:
            self.add_regions(regions)
        else:
            self._update_references()

    def flush(self):
        self._regions.flush()

    def add_region(self, region, flush=True):
        # super-method, calls below '_add_region'
        ix = GenomicRegions.add_region(self, region)
        if flush:
            self.flush()
            self._update_references()
        return ix

    def _add_region(self, region):
        ix = self._max_region_ix + 1
        
        # actually append
        row = self._regions.row
        row['ix'] = ix
        row['chromosome'] = region.chromosome
        row['start'] = region.start
        row['end'] = region.end
        if hasattr(region, 'strand') and region.strand is not None:
            row['strand'] = region.strand

        for name in self._regions.colnames[5:]:
            if hasattr(region, name):
                row[name] = getattr(region, name)

        row.append()
        
        if ix > self._max_region_ix:
            self._max_region_ix = ix

        return ix

    def _update_references(self):
        chromosomes = []
        for region in self.regions(lazy=True):
            if len(chromosomes) == 0 or chromosomes[-1] != region.chromosome:
                chromosomes.append(region.chromosome)

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
        try:
            l = len(regions)
            _log = True
        except TypeError:
            l = None
            _log = False

        pb = RareUpdateProgressBar(max_value=l)
        if _log:
            pb.start()

        for i, region in enumerate(regions):
            self.add_region(region, flush=False)
            if _log:
                pb.update(i)
        if _log:
            pb.finish()

        self.flush()
        self._update_references()

    def data(self, key, value=None):
        """
        Retrieve or add vector-data to this object. If there is exsting data in this
        object with the same name, it will be replaced

        :param key: Name of the data column
        :param value: vector with region-based data (one entry per region)
        """
        if key not in self._regions.colnames:
            raise KeyError("%s is unknown region attribute" % key)

        if value is not None:
            for i, row in enumerate(self._regions):
                row[key] = value[i]
                row.update()
            self._regions.flush()

        return (row[key] for row in self._regions)

    def _get_region_ix(self, region):
        """
        Get index from other region properties (chromosome, start, end)
        """
        condition = "(start == %d) & (end == %d) & (chromosome == '%s')"
        condition = condition % (region.start, region.end, region.chromosome)
        for res in self._regions.where(condition):
            return res["ix"]
        return None

    def _row_to_region(self, row, lazy=False, auto_update=True):
        if lazy:
            return LazyGenomicRegion(row, auto_update=auto_update)

        kwargs = {}
        for name in self._regions.colnames:
            if name not in RegionsTable.RegionDescription().columns.keys():
                kwargs[name] = row[name]
        return GenomicRegion(chromosome=row["chromosome"], start=row["start"],
                             end=row["end"], ix=row["ix"], **kwargs)

    @property
    def regions(self):
        """
        Iterate over genomic regions in this object.

        Will return a :class:`~Node` object in every iteration.
        Can also be used to get the number of regions by calling
        len() on the object returned by this method.

        :param lazy: If True, will only retrieve properties in
                     a lazy fashion, i.e. on request

        :return: RegionIter
        """
        this = self

        class RegionIter:
            def __init__(self):
                self.iter = iter(this._regions)
                self.lazy = False
                self.auto_update = True

            def __iter__(self):
                return self
            
            def next(self):
                return this._row_to_region(self.iter.next(), lazy=self.lazy)
            
            def __len__(self):
                return len(this._regions)

            def __call__(self, lazy=False, auto_update=True):
                self.lazy = lazy
                self.auto_update = auto_update
                return iter(self)

            def __getitem__(self, item):
                res = this._regions[item]

                if isinstance(res, np.ndarray):
                    regions = []
                    for region in res:
                        regions.append(this._row_to_region(region, lazy=self.lazy,
                                                           auto_update=self.auto_update))
                    return regions
                else:
                    return this._row_to_region(res, lazy=self.lazy,
                                               auto_update=self.auto_update)

        return RegionIter()

    def subset(self, region, lazy=False, auto_update=True):
        """
        Iterate over a subset of regions given the specified key.

        :param region: A :class:`~kaic.data.genomic.GenomicRegion` object,
                       or a list of the former.
        :param lazy: Load region attributes on demand only.
        :param auto_update: Auto update regions upon modification
        :return: Iterator over the specified subset of regions
        """
        if isinstance(region, slice):
            for row in self._regions.where("(ix >= {}) & (ix < {})".format(region.start, region.stop)):
                sub_region = self._row_to_region(row, lazy=lazy, auto_update=auto_update)
                yield sub_region
        elif isinstance(region, int):
            sub_region = self._row_to_region(self._regions[region], lazy=lazy, auto_update=auto_update)
            yield sub_region
        elif isinstance(region, list) and len(region) > 0 and isinstance(region[0], int):
            for ix in region:
                sub_region = self._row_to_region(self._regions[ix], lazy=lazy, auto_update=auto_update)
                yield sub_region
        else:
            if isinstance(region, str):
                region = GenomicRegion.from_string(region)

            if isinstance(region, GenomicRegion):
                regions = [region]
            else:
                regions = region

            for r in regions:
                if isinstance(r, str):
                    r = GenomicRegion.from_string(r)

                query = '('
                if r.chromosome is not None:
                    query += "(chromosome == '%s') & " % r.chromosome
                if r.end is not None:
                    query += "(start <= %d) & " % r.end
                if r.start is not None:
                    query += "(end >= %d) & " % r.start
                if query.endswith(' & '):
                    query = query[:-3]
                query += ')'

                if len(query) == 2:
                    for region in self.regions(lazy=lazy, auto_update=auto_update):
                        yield region
                else:
                    for row in self._regions.where(query):
                        sub_region = self._row_to_region(row, lazy=lazy, auto_update=auto_update)
                        yield sub_region


class Node(GenomicRegion, TableObject):
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
        super(Node, self).__init__(chromosome=chromosome, start=start, end=end, ix=ix)
    
    def __repr__(self):
        if self.ix is None:
            return "%s, %d-%d" % (self.chromosome, self.start, self.end)
        else:
            return "%d: %s, %d-%d" % (self.ix, self.chromosome, self.start, self.end)


class LazyNode(LazyGenomicRegion, Node):
    def __init__(self, row, ix=None):
        LazyGenomicRegion.__init__(self, row=row, ix=ix)


class Edge(TableObject):
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
    def __init__(self, source, sink, **kwargs):
        """
        :param source: The index of the "source" genomic region
                       or :class:`~Node` object.
        :param sink: The index of the "sink" genomic region
                     or :class:`~Node` object.
        :param data: The weight or of the edge or a dictionary with
                     other fields
        """
        self._source = source
        self._sink = sink
        self.field_names = []

        for key, value in kwargs.iteritems():
            setattr(self, key, value)
            self.field_names.append(key)

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
        base_info = "%d--%d" % (self.source, self.sink)
        for field in self.field_names:
            base_info += "\n\t%s: %s" % (field, str(getattr(self, field)))
        return base_info


class LazyEdge(Edge):
    def __init__(self, row, nodes_table=None, auto_update=True):
        self.reserved = {'_row', '_nodes_table', 'auto_update', '_source_node', '_sink_node'}
        self._row = row
        self._nodes_table = nodes_table
        self.auto_update = auto_update
        self._source_node = None
        self._sink_node = None

    def _set_item(self, item, value):
        self._row[item] = value
        if self.auto_update:
            self.update()

    def __getattr__(self, item):
        if item == 'reserved' or item in self.reserved:
            return object.__getattribute__(self, item)
        try:
            return self._row[item]
        except KeyError:
            raise AttributeError("Attribute not supported (%s)" % str(item))

    def __setattr__(self, key, value):
        if key == 'reserved' or key in self.reserved:
            super(LazyEdge, self).__setattr__(key, value)
        else:
            self._row[key] = value
            if self.auto_update:
                self.update()

    def update(self):
        self._row.update()

    @property
    def source(self):
        return self._row['source']

    @property
    def sink(self):
        return self._row['sink']

    @property
    def source_node(self):
        if self._nodes_table is None:
            raise RuntimeError("Must set the _nodes_table attribute before calling this method!")

        if self._source_node is None:
            source_row = self._nodes_table[self.source]
            return LazyNode(source_row)
        return self._source_node

    @property
    def sink_node(self):
        if self._nodes_table is None:
            raise RuntimeError("Must set the _nodes_table attribute before calling this method!")

        if self._sink_node is None:
            sink_row = self._nodes_table[self.sink]
            return LazyNode(sink_row)
        return self._sink_node

    def __repr__(self):
        base_info = "%d--%d" % (self.source, self.sink)
        return base_info


class RegionPairs(Maskable, MetaContainer, RegionsTable):
    """
    Class for working with data associated with pairs of regions.

    Generally, a RegionPairs object has two components:

    - Nodes or regions: (Non-overlapping) genomic regions
      obtained by splitting the genome into distinct pieces.
      See also :class:`~GenomicRegion` and :class:`~RegionsTable`

    - Edges or contacts: Pairs of genomic regions. See also
      :class:`~Edge`
    """

    _classid = 'REGIONPAIRS'

    class EntryDescription(t.IsDescription):
        source = t.Int32Col(pos=0)
        sink = t.Int32Col(pos=1)

    class EdgeIter(object):
        def __init__(self, this, _iter=None):
            self.this = this
            if _iter is None:
                self.iter = iter(this._edges)
            else:
                self.iter = iter(_iter)
            self.row_conversion_args = list()
            self.row_conversion_kwargs = dict()
            self.only_intrachromosomal = False
            self.regions_dict = None

        def __getitem__(self, item):
            res = self.this._edges[item]

            if isinstance(res, np.ndarray):
                edges = []
                for edge in res:
                    edges.append(self.this._row_to_edge(edge, *self.row_conversion_args, **self.row_conversion_kwargs))
                return edges
            else:
                edge = self.this._row_to_edge(res, *self.row_conversion_args, **self.row_conversion_kwargs)
                return edge

        def __iter__(self):
            if self.only_intrachromosomal:
                self.regions_dict = self.this.regions_dict
            return self

        def __call__(self, *args, **kwargs):
            if 'only_intrachromosomal' in kwargs:
                self.only_intrachromosomal = kwargs['only_intrachromosomal']
                del kwargs['only_intrachromosomal']
            self.row_conversion_args = args
            self.row_conversion_kwargs = kwargs
            return iter(self)

        def next(self):
            row = self.iter.next()
            if self.only_intrachromosomal:
                while self.regions_dict[row['source']].chromosome != self.regions_dict[row['sink']].chromosome:
                    row = self.iter.next()
            return self.this._row_to_edge(row, *self.row_conversion_args, **self.row_conversion_kwargs)

        def __len__(self):
            return len(self.this._edges)

    def __init__(self, file_name=None, mode='a', additional_fields=None, tmpdir=None,
                 _table_name_nodes='nodes', _table_name_edges='edges'):

        """
        Initialize a :class:`~RegionPairs` object.

        :param file_name: Path to a save file
        :param mode: File mode to open underlying file
        :param additional_fields: Additional fields (in PyTables notation) associated with
                                  edge data, e.g. {'weight': tables.Float32Col()}
        :param _table_name_nodes: (Internal) name of the HDF5 node for regions
        :param _table_name_edges: (Internal) name of the HDF5 node for edges
        """

        # private variables
        self._max_node_ix = -1

        if file_name is not None:
            file_name = os.path.expanduser(file_name)

        # initialize inherited objects
        RegionsTable.__init__(self, file_name=file_name, _table_name_regions=_table_name_nodes,
                              mode=mode, tmpdir=tmpdir)
        Maskable.__init__(self, self.file)
        MetaContainer.__init__(self, self.file)

        # create edge table
        if _table_name_edges in self.file.root:
            self._edges = self.file.get_node('/', _table_name_edges)
        else:
            basic_fields = self._get_field_dict(additional_fields=additional_fields)

            self._edges = MaskedTable(self.file.root, _table_name_edges, basic_fields)

        # index edge table
        create_col_index(self._edges.cols.source)
        create_col_index(self._edges.cols.sink)

        # update field names
        self._source_field_ix = 0
        self._sink_field_ix = 0
        self.field_names = []
        self._field_dict = self._edges.coldescrs
        for i, name in enumerate(self._edges.colnames):
            if not name.startswith("_"):
                self.field_names.append(name)
            if name == 'source':
                self._source_field_ix = i
            if name == 'sink':
                self._sink_field_ix = i

    def _get_field_dict(self, additional_fields=None):
        basic_fields = RegionMatrixTable.EntryDescription().columns.copy()
        if additional_fields is not None:
            if not isinstance(additional_fields, dict) and issubclass(additional_fields, t.IsDescription):
                # IsDescription subclass case
                additional_fields = additional_fields.columns

            current = len(basic_fields)
            for key, value in sorted(additional_fields.iteritems(), key=lambda x: x[1]._v_pos):
                if key not in basic_fields:
                    if value._v_pos is not None:
                        value._v_pos = current
                        current += 1
                    basic_fields[key] = value
        return basic_fields

    def add_node(self, node, flush=True):
        """
        Add a :class:`~Node` or :class:`~GenomicRegion`.

        :param node: :class:`~Node` or :class:`~GenomicRegion`,
                     see :func:`~RegionsTable.add_region` for details
        :param flush: Write data to file immediately after import.
        """
        return self.add_region(node, flush)

    def add_edge(self, edge, check_nodes_exist=True, flush=True, replace=False, row=None):
        """
        Add an edge to this object.

        :param edge: :class:`~Edge`, dict with at least the
                     attributes source and sink, optionally weight,
                     or a list of length 2 (source, sink) or 3
                     (source, sink, weight).
        :param check_nodes_exist: Make sure that there are nodes
                                  that match source and sink indexes
        :param flush: Write data to file immediately after import
        :param replace: If row is provided, replace values in existing edge with the ones in edge
        :param row: PyTables row object representing an edge. If provided, edge will be used to
                    modify existing row.
        """
        source = None
        sink = None

        # object
        is_object = True
        try:
            source = edge.source
            sink = edge.sink
        except AttributeError:
            is_object = False

        # dictionary
        is_dict = False
        if not is_object:
            is_dict = True
            try:
                source = edge['source']
                sink = edge['sink']
            except TypeError:
                is_dict = False

        # list
        is_list = False
        if not is_object and not is_dict:
            is_list = True
            try:
                source = edge[self._source_field_ix]
                sink = edge[self._sink_field_ix]
            except TypeError:
                is_list = False

        if source is None and sink is None:
            raise ValueError("Edge type not recognised (%s)" % str(type(edge)))

        if check_nodes_exist:
            n_regions = len(self._regions)
            if source >= n_regions or sink >= n_regions:
                raise ValueError("Node index exceeds number of nodes in object")

        if is_object:
            new_edge = self._edge_from_object(edge)
        elif is_dict:
            new_edge = self._edge_from_dict(edge)
        elif is_list:
            new_edge = self._edge_from_list(edge)
        else:
            raise ValueError("Edge type not recognised (%s)" % str(type(edge)))

        self._add_edge(new_edge, row=row, replace=replace)

        if flush:
            self.flush()

    def _add_edge(self, edge, row, replace=False):
        source, sink = edge.source, edge.sink
        if source > sink:
            source, sink = sink, source

        update = True
        if row is None:
            update = False
            row = self._edges.row
        row['source'] = source
        row['sink'] = sink
        for name in self.field_names:
            if not name == 'source' and not name == 'sink':
                try:
                    value = getattr(edge, name)
                    if replace or not update:
                        row[name] = value
                    else:
                        row[name] += value
                except AttributeError:
                    pass
        if update:
            row.update()
        else:
            row.append()

    def _edge_from_object(self, edge):
        return edge

    def _edge_from_dict(self, edge):
        source, sink = edge['source'], edge['sink']

        attributes = dict()
        for name, value in edge.iteritems():
            if not name == 'source' and not name == 'sink':
                 attributes[name] = value

        return Edge(source, sink, **attributes)

    def _edge_from_list(self, edge):
        source, sink = edge[self._source_field_ix], edge[self._sink_field_ix]

        attributes = dict()
        for i, name in enumerate(self.field_names):
            if not name == 'source' and not name == 'sink':
                try:
                    attributes[name] = edge[i]
                except IndexError:
                    break

        return Edge(source, sink, **attributes)

    def add_nodes(self, nodes):
        """
        Bulk-add nodes from a list.

        :param nodes: List (or iterator) of nodes. See
                      :func:`~RegionMatrixTable.add_node`
                      for details.
        """
        self.add_regions(nodes)

    def add_edges(self, edges):
        """
        Bulk-add edges from a list.

        :param edges: List (or iterator) of edges. See
                      :func:`~RegionMatrixTable.add_edge`
                      for details
        """
        for edge in edges:
            self.add_edge(edge, flush=False)
        self.flush(flush_nodes=False)

    def flush(self, flush_nodes=True, flush_edges=True, update_index=True):
        """
        Write data to file and flush buffers.

        :param flush_nodes: Flush nodes tables
        :param flush_edges: Flush edges table
        :param update_index: Update mask indices in edges table
        """
        if flush_nodes:
            self._regions.flush()

        if flush_edges:
            self._edges.flush(update_index=update_index)

    def edge_subset(self, key=slice(0, None, None), lazy=False, auto_update=True,
                    only_intrachromosomal=False):
        """
        Get a subset of edges.

        :param key: Possible key types are:

                    Region types

                    - Node: Only the ix of this node will be used for
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
        :param lazy: Enable lazy loading of edge attributes
        :param auto_update: Automatically update edge attributes on change
        :return: generator (:class:`~Edge`)
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

        row_ranges = list(self._get_node_ix_ranges(nodes_ix_row))
        col_ranges = list(self._get_node_ix_ranges(nodes_ix_col))

        # fill matrix with weights
        for row_range in row_ranges:
            for col_range in col_ranges:
                for edge_row in self._edge_row_range(row_range[0], row_range[1],
                                                     col_range[0], col_range[1],
                                                     only_intrachromosomal=only_intrachromosomal):
                    yield self._row_to_edge(edge_row, lazy=lazy, auto_update=auto_update)

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

    def _edge_row_range(self, source_start, source_end, sink_start, sink_end, only_intrachromosomal=False):
        condition = "(source > %d) & (source < %d) & (sink > %d) & (sink < %d)"
        condition1 = condition % (source_start-1, source_end+1, sink_start-1, sink_end+1)
        condition2 = condition % (sink_start-1, sink_end+1, source_start-1, source_end+1)

        if source_start > sink_start:
            condition1, condition2 = condition2, condition1

        regions_dict = None
        if only_intrachromosomal:
            regions_dict = self.regions_dict

        overlap = range_overlap(source_start, source_end, sink_start, sink_end)

        for edge_row in self._edges.where(condition1):
            if (only_intrachromosomal and
                    regions_dict[edge_row['source']].chromosome != regions_dict[edge_row['sink']].chromosome):
                continue
            yield edge_row

        for edge_row in self._edges.where(condition2):
            if overlap is not None:
                if (overlap[0] <= edge_row['source'] <= overlap[1]) and (overlap[0] <= edge_row['sink'] <= overlap[1]):
                    continue

            if (only_intrachromosomal and
                    regions_dict[edge_row['source']].chromosome != regions_dict[edge_row['sink']].chromosome):
                continue
            yield edge_row

    def _get_node_ix_ranges(self, nodes_ix=None):
        if not isinstance(nodes_ix, list):
            nodes_ix = [nodes_ix]

        # get range generator
        return ranges(nodes_ix)

    def _getitem_nodes(self, key, as_index=False):
        # 'chr1:1234:56789'
        if isinstance(key, str):
            key = GenomicRegion.from_string(key)

        # Node('chr1', 1234, 56789, ix=0)
        if isinstance(key, Node):
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
                region_nodes = [self._row_to_region(row) for row in self._regions.where(condition)]

            return region_nodes

        # 1:453
        if isinstance(key, slice):
            if as_index:
                return [row['ix'] for row in self._regions.iterrows(key.start, key.stop, key.step)]
            else:
                return [self._row_to_region(row) for row in self._regions.iterrows(key.start, key.stop, key.step)]

        # 432
        if isinstance(key, int):
            row = self._regions[key]
            if as_index:
                return row['ix']
            else:
                return self._row_to_region(row)

        # [item1, item2, item3]
        all_nodes_ix = []
        for item in key:
            nodes_ix = self._getitem_nodes(item, as_index=as_index)
            if isinstance(nodes_ix, list):
                all_nodes_ix += nodes_ix
            else:
                all_nodes_ix.append(nodes_ix)
        return all_nodes_ix

    def _row_to_node(self, row, lazy=False):
        if lazy:
            return LazyNode(row)
        return Node(chromosome=row["chromosome"], start=row["start"],
                    end=row["end"], ix=row["ix"])

    def get_node(self, key):
        """
        Get a single node by key.

        :param key: For possible key types see :func:`~RegionMatrixTable.__getitem__`
        :return: A :class:`~Node` matching key
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

        :param key: For possible key types see :func:`~RegionMatrixTable.__getitem__`
        :return: A list of :class:`~Node` objects matching key
        """
        return self._getitem_nodes(key)

    def _row_to_edge(self, row, lazy=False, auto_update=True):
        if not lazy:
            source = row["source"]
            sink = row["sink"]
            d = dict()
            for field in self.field_names:
                if field != 'source' and field != 'sink':
                    d[field] = row[field]

            source_node_row = self._regions[source]
            source_node = self._row_to_node(source_node_row)
            sink_node_row = self._regions[sink]
            sink_node = self._row_to_node(sink_node_row)
            return Edge(source_node, sink_node, **d)
        else:
            return LazyEdge(row, self._regions, auto_update=auto_update)

    def get_edge(self, ix, lazy=False):
        """
        Get an edge from this object's edge list.

        :param ix: integer
        :param lazy: Use lazy loading of object attributes. Do not
                     use lazy objects outside of loop iterations!
        :return:
        """
        return self._row_to_edge(self._edges[ix], lazy=lazy)

    def nodes(self):
        """
        Iterator over this object's nodes/regions.

        See :func:`~RegionsTable.regions` for details.
        :return: Iterator over :class:`~GenomicRegions`
        """
        return self._nodes_iter()

    def _nodes_iter(self):
        return self.regions()

    @property
    def edges(self):
        """
        Iterate over :class:`~Edge` objects.

        :return: Iterator over :class:`~Edge`
        """
        return self._edges_iter()

    def _edges_iter(self):
        return RegionPairs.EdgeIter(self)

    def _is_sorted(self, sortby):
        column = getattr(self._edges.cols, sortby)
        if (column.index is None or
                not column.index.is_csi):
            return False
        return True

    def edges_sorted(self, sortby, reverse=False, *args, **kwargs):
        """
        Iterate over edges sorted by a specific column.

        :param sortby: Name of column to sort over
        :param reverse: Iterate in reverse order
        :return: EdgeIter iterator
        """
        # ensure sorting on qname_ix column
        column = getattr(self._edges.cols, sortby)

        if not self._is_sorted(sortby):
            try:
                logging.info("Sorting %s..." % sortby)
                if not column.is_indexed:
                    column.create_csindex()
                elif not column.index.is_csi:
                    column.reindex()
            except t.exceptions.FileModeError:
                raise RuntimeError("This object is not sorted by requested column! "
                                   "Cannot sort manually, because file is in read-only mode.")
        if reverse:
            step = -1
        else:
            step = None
        edge_iter = RegionMatrixTable.EdgeIter(self, _iter=self._edges.itersorted(sortby, step=step))
        return edge_iter(*args, **kwargs)

    def __iter__(self):
        return self.edges

    def __len__(self):
        return len(self._edges)


class AccessOptimisedRegionPairs(RegionPairs):
    """
    Extends :class:`~RegionPairs` with a backend that partitions edges into chromosomes.

    This partitioning should greatly speed up edge and matrix queries for large Hi-C data sets,
    such as high-resolution (<=10kb) human Hi-C maps.

    Iterating over sorted edges performance is somewhat reduced due to the fact that we have to
    integrate tens to hundreds of tables in the sorting.
    """

    _classid = 'ACCESSOPTIMISEDREGIONPAIRS'

    class EdgeIter(RegionPairs.EdgeIter):
        """
        Class providing iterator functionality to a :class:`~RegionPairs` object.
        """
        def __init__(self, this, _iter=None):
            RegionPairs.EdgeIter.__init__(self, this, _iter=_iter)
            self.iter = _iter
            self.interchromosomal = True
            self.intrachromosomal = True

        def __getitem__(self, item):
            if isinstance(item, int):
                return self.this.get_edge(item, intrachromosomal=self.intrachromosomal,
                                          interchromosomal=self.interchromosomal)
            elif isinstance(item, slice):
                edges = []
                start = 0 if item.start is None else item.start
                stop = len(self.this.edges) if item.stop is None else item.stop
                step = 1 if item.step is None else item.step
                if step != 1:
                    raise ValueError("Step sizes != 1 not currently supported in slices. %s" % str(step))

                l = 0
                for edge_table in self.this._edge_table_iter(intrachromosomal=self.intrachromosomal,
                                                             interchromosomal=self.interchromosomal):
                    # not yet in range
                    if start >= l + len(edge_table):
                        l += len(edge_table)
                        continue
                    # over range - can stop here
                    if stop < l:
                        break

                    # in range, get edges
                    r = (max(0, start - l), min(len(edge_table), stop - l))
                    print(r)
                    res = edge_table[r[0]:r[1]]
                    for edge in res:
                        edges.append(self.this._row_to_edge(edge,
                                                            *self.row_conversion_args,
                                                            **self.row_conversion_kwargs))
                    l += len(edge_table)
                return edges

        def __iter__(self):
            return self

        def __call__(self, *args, **kwargs):
            if 'only_intrachromosomal' in kwargs:
                self.interchromosomal = False
                self.intrachromosomal = True
                del kwargs['only_intrachromosomal']
            if 'intrachromosomal' in kwargs:
                self.intrachromosomal = kwargs['intrachromosomal']
                del kwargs['intrachromosomal']
            if 'interchromosomal' in kwargs:
                self.interchromosomal = kwargs['interchromosomal']
                del kwargs['interchromosomal']
            self.row_conversion_args = args
            self.row_conversion_kwargs = kwargs
            return iter(self)

        def next(self):
            if self.iter is None:
                self.iter = self.this._edge_row_iter(intrachromosomal=self.intrachromosomal,
                                                     interchromosomal=self.interchromosomal)
            row = self.iter.next()
            return self.this._row_to_edge(row, *self.row_conversion_args, **self.row_conversion_kwargs)

        def __len__(self):
            return len(self.this)

    def __init__(self, file_name=None, mode='a', tmpdir=None, additional_fields=None,
                 _table_name_nodes='nodes', _table_name_edges='edges'):
        # private variables
        self._max_node_ix = -1

        if file_name is not None:
            file_name = os.path.expanduser(file_name)

        # initialize inherited objects
        RegionsTable.__init__(self, file_name=file_name, _table_name_regions=_table_name_nodes,
                              mode=mode, tmpdir=tmpdir)
        Maskable.__init__(self, self.file)
        MetaContainer.__init__(self, self.file)

        # create edge table
        self._field_dict = None
        self.field_names = None
        self._edge_table_dict = dict()
        self._source_field_ix = 0
        self._sink_field_ix = 1

        # existing one
        if _table_name_edges in self.file.root:
            self._edges = self.file.get_node('/', _table_name_edges)
            for edge_table in self._edges._f_iter_nodes():
                if self._field_dict is None:
                    self._field_dict = edge_table.coldescrs
                    self._update_field_names(edge_table=edge_table)
                source_partition = edge_table.attrs['source_partition']
                sink_partition = edge_table.attrs['sink_partition']
                self._edge_table_dict[(source_partition, sink_partition)] = edge_table
        else:
            # create edge table definition
            self._field_dict = self._get_field_dict(additional_fields=additional_fields)

            self._edges = self.file.create_group('/', _table_name_edges)
            # will always have 0-0 partition
            edge_table = self._create_edge_table(0, 0)
            self._update_field_names(edge_table=edge_table)

        # update partitions
        self._update_partitions()

    def _update_field_names(self, edge_table=None):
        """
        Set internal object variables related to edge table field names.
        """
        if edge_table is None:
            for et in self._edges._f_iter_nodes():
                edge_table = et
                break

        if edge_table is None:
            return

        # update field names
        self._source_field_ix = 0
        self._sink_field_ix = 0
        self.field_names = []
        for i, name in enumerate(edge_table.colnames):
            if not name.startswith("_"):
                self.field_names.append(name)
            if name == 'source':
                self._source_field_ix = i
            if name == 'sink':
                self._sink_field_ix = i

    def _update_partitions(self):
        """
        Update the list of partition break points (split by chromosome)
        """
        self.partitions = []
        previous_chromosome = None
        for i, region in enumerate(self.regions(lazy=True)):
            if region.chromosome != previous_chromosome and previous_chromosome is not None:
                self.partitions.append(i)
            previous_chromosome = region.chromosome

    def flush(self, flush_nodes=True, flush_edges=True, update_index=True):
        """
        Write data to file and flush buffers.

        :param flush_nodes: Flush nodes tables
        :param flush_edges: Flush edges table
        :param update_index: Update mask indices in edges table
        """
        if flush_nodes:
            self._regions.flush()
            # update partitions
            self._update_partitions()

        if flush_edges:
            for edge_table in self._edges._f_iter_nodes():
                edge_table.flush(update_index=update_index)

    def _get_field_dict(self, additional_fields=None):
        """
        Generate a dictionary of PyTables fields to create edge table.

        Save fields dict in object variable - we need to generate a lot of tables.
        """
        if self._field_dict is not None:
            return self._field_dict
        return RegionPairs._get_field_dict(self, additional_fields=additional_fields)

    def _get_partition_ix(self, region_ix):
        """
        Bisect the partition table to get the partition index for a region index.
        """
        return bisect_right(self.partitions, region_ix)

    def _create_edge_table(self, source_partition, sink_partition):
        """
        Create and register an edge table for a partition combination.
        """
        edge_table = MaskedTable(self._edges,
                                 'chrpair_' + str(source_partition) + '_' + str(sink_partition),
                                 self._field_dict)
        edge_table.attrs['source_partition'] = source_partition
        edge_table.attrs['sink_partition'] = sink_partition

        # index
        create_col_index(edge_table.cols.source)
        create_col_index(edge_table.cols.sink)

        self._edge_table_dict[(source_partition, sink_partition)] = edge_table
        return edge_table

    def _get_edge_table(self, source, sink):
        """
        Return an edge table for this particular region index combination.
        """
        if source > sink:
            source, sink = sink, source

        source_partition = self._get_partition_ix(source)
        sink_partition = self._get_partition_ix(sink)

        if (source_partition, sink_partition) in self._edge_table_dict:
            return self._edge_table_dict[(source_partition, sink_partition)]

        edge_table = self._create_edge_table(source_partition, sink_partition)

        return edge_table

    def _edge_table_iter(self, intrachromosomal=True, interchromosomal=True):
        """
        Iterate over internal edge tables.

        :param intrachromosomal: If true, include intra-chromosomal edge tables
        :param interchromosomal: If true, include inter-chromosomal edge tables
        :return: Edge table iterator
        """
        # intra-chromosomal
        if intrachromosomal:
            for i in xrange(len(self.partitions) + 1):
                if (i, i) in self._edge_table_dict:
                    yield self._edge_table_dict[(i, i)]

        # inter-chromosomal
        if interchromosomal:
            for i in xrange(len(self.partitions) + 1):
                for j in xrange(i + 1, len(self.partitions) + 1):
                    if (i, j) in self._edge_table_dict:
                        yield self._edge_table_dict[(i, j)]

    def _edge_from_list(self, edge):
        source, sink = edge[self._source_field_ix], edge[self._sink_field_ix]
        self._get_edge_table(source, sink)

        return RegionPairs._edge_from_list(self, edge)

    def _add_edge(self, edge, row, replace=False):
        """
        Add an edge to an internal edge table.
        """
        source, sink = edge.source, edge.sink
        if source > sink:
            source, sink = sink, source

        update = True
        if row is None:
            update = False
            table = self._get_edge_table(source, sink)
            row = table.row
        row['source'] = source
        row['sink'] = sink
        for name in self.field_names:
            if not name == 'source' and not name == 'sink':
                try:
                    value = getattr(edge, name)
                    if replace or not update:
                        row[name] = value
                    else:
                        row[name] += value
                except AttributeError:
                    pass
        if update:
            row.update()
        else:
            row.append()

    def get_edge(self, item, intrachromosomal=True, interchromosomal=True,
                 *row_conversion_args, **row_conversion_kwargs):
        """
        Get an edge by index.

        :param intrachromosomal: If true, include intra-chromosomal edge tables in index count
        :param interchromosomal: If true, include inter-chromosomal edge tables in index count
        :param row_conversion_args: Arguments passed to :func:`RegionPairs._row_to_edge`
        :param row_conversion_args: Keyword arguments passed to :func:`RegionPairs._row_to_edge`
        :return: :class:`~Edge`
        """
        l = 0
        for edge_table in self._edge_table_iter(intrachromosomal=intrachromosomal,
                                                interchromosomal=interchromosomal):
            if l <= item < l + len(edge_table):
                res = edge_table[item - l]
                return self._row_to_edge(res, *row_conversion_args, **row_conversion_kwargs)
            l += len(edge_table)
        raise IndexError("index out of range (%d)" % item)

    @property
    def edges(self):
        """
        Iterate over :class:`~Edge` objects.

        :param lazy: Enable lazy loading of edge attributes,
                     only works in the loop iteration this
                     edge is accessed.
        :return: Iterator over :class:`~Edge`
        """
        return self._edges_iter()

    def _edges_iter(self):
        return AccessOptimisedRegionPairs.EdgeIter(self)

    def _edge_row_iter(self, intrachromosomal=True, interchromosomal=True):
        """
        Yield rows in edge tables, ordered by partition.
        """
        for edge_table in self._edge_table_iter(intrachromosomal=intrachromosomal, interchromosomal=interchromosomal):
            for row in edge_table:
                yield row

    def _partition_ix_range(self, start, stop):
        """
        Get a range of partitions with start and stop indices per partition from global start and stop region indices.

        :param start: Region start index
        :param stop: Region stop index
        :return: tuple, where first element is
        """
        start_partition = self._get_partition_ix(start)
        stop_partition = self._get_partition_ix(stop)

        def _is_start_of_partition(start_ix, partition_ix):
            if partition_ix == 0:
                if start_ix == 0:
                    return True
            else:
                if start_ix == self.partitions[partition_ix - 1]:
                    return True
            return False

        def _is_end_of_partition(stop_ix, partition_ix):
            if partition_ix == len(self.partitions):
                if stop_ix == len(self.regions)-1:
                    return True
            else:
                if stop_ix == self.partitions[partition_ix]-1:
                    return True
            return False

        if start_partition == stop_partition:
            complete = _is_start_of_partition(start, start_partition) and _is_end_of_partition(stop, stop_partition)
            return [(start, stop, start_partition, complete)]

        partition_ranges = []
        start_range_complete = _is_start_of_partition(start, start_partition)
        start_range = (start, self.partitions[start_partition] - 1, start_partition, start_range_complete)
        partition_ranges.append(start_range)

        for i in xrange(start_partition + 1, stop_partition):
            partition_ranges.append((self.partitions[i-1], self.partitions[i]-1, i, True))

        stop_range_complete = _is_end_of_partition(stop, stop_partition)
        stop_range = (self.partitions[stop_partition - 1], stop, stop_partition, stop_range_complete)
        partition_ranges.append(stop_range)
        return partition_ranges

    def _edge_row_range(self, source_start, source_end, sink_start, sink_end, only_intrachromosomal=False):
        """
        Iterate over a range of rows in this object's edge tables.

        Rows are selected based on region indices of interacting regions.
        """
        source_partition_ranges = self._partition_ix_range(source_start, source_end)
        sink_partition_ranges = self._partition_ix_range(sink_start, sink_end)

        covered = set()
        for source_partition_range in source_partition_ranges:
            source_start, source_end, source_partition, source_complete = source_partition_range

            for sink_partition_range in sink_partition_ranges:
                sink_start, sink_stop, sink_partition, sink_complete = sink_partition_range

                if only_intrachromosomal and source_partition != sink_partition:
                    continue

                if sink_partition < source_partition:
                    key = (sink_partition, source_partition)
                else:
                    key = (source_partition, sink_partition)

                if key in covered:
                    continue
                else:
                    covered.add(key)

                if key in self._edge_table_dict:
                    table = self._edge_table_dict[key]

                    # entire partition is requested, no need for where query
                    if source_complete and sink_complete:
                        for edge_row in table:
                            yield edge_row
                    else:
                        condition = "(source > %d) & (source < %d) & (sink > %d) & (sink < %d)"
                        condition1 = condition % (source_start - 1, source_end + 1, sink_start - 1, sink_end + 1)
                        condition2 = condition % (sink_start - 1, sink_end + 1, source_start - 1, source_end + 1)

                        if source_start > sink_start:
                            condition1, condition2 = condition2, condition1

                        overlap = range_overlap(source_start, source_end, sink_start, sink_end)

                        for edge_row in table.where(condition1):
                            yield edge_row

                        for edge_row in table.where(condition2):
                            if overlap is not None:
                                if (overlap[0] <= edge_row['source'] <= overlap[1]) and (overlap[0] <= edge_row['sink'] <= overlap[1]):
                                    continue

                            yield edge_row

    def _is_sorted(self, sortby):
        """
        For each edge table, check if it is sorted.
        """
        for edge_table in self._edge_table_iter():
            column = getattr(edge_table.cols, sortby)
            if (column.index is None or
                    not column.index.is_csi):
                return False
        return True

    def _edge_row_iter_sorted(self, sortby, step=None, intrachromosomal=True, interchromosomal=True):
        """
        Yield rows in edge tables, ordered by partition.
        """
        table_iterators = []
        for edge_table in self._edge_table_iter(intrachromosomal=intrachromosomal, interchromosomal=interchromosomal):
            table_iterators.append(iter(edge_table.itersorted(sortby, step=step)))

        rows = []
        for i, table_iterator in enumerate(table_iterators):
            try:
                row = table_iterator.next()
                rows.append(row)
            except StopIteration:
                del table_iterators[i]

        while len(table_iterators) > 0:
            # find current minimum or maximum
            current = None
            current_ix = None
            if step is None or step >= 0:
                for i, row in enumerate(rows):
                    if current is None or row[sortby] < current:
                        current = row[sortby]
                        current_ix = i
            else:
                for i, row in enumerate(rows):
                    if current is None or row[sortby] > current:
                        current = row[sortby]
                        current_ix = i

            yield rows[current_ix]

            try:
                rows[current_ix] = table_iterators[current_ix].next()
            except StopIteration:
                del table_iterators[current_ix]
                del rows[current_ix]

    def edges_sorted(self, sortby, reverse=False, *args, **kwargs):
        """
        Iterate over edges sorted by *sortby*.
        """
        for edge_table in self._edge_table_iter():
            # ensure sorting on sortby column
            column = getattr(edge_table.cols, sortby)

            if not self._is_sorted(sortby):
                try:
                    logging.info("Sorting %s..." % sortby)
                    if not column.is_indexed:
                        column.create_csindex()
                    elif not column.index.is_csi:
                        column.reindex()
                except t.exceptions.FileModeError:
                    raise RuntimeError("This object is not sorted by requested column! "
                                       "Cannot sort manually, because file is in read-only mode.")
        if reverse:
            step = -1
        else:
            step = None
        edge_iter = AccessOptimisedRegionPairs.EdgeIter(self, _iter=self._edge_row_iter_sorted(sortby, step=step))
        return edge_iter(*args, **kwargs)

    def __len__(self):
        l = 0
        for edge_table in self._edge_table_iter():
            l += len(edge_table)
        return l

    def __iter__(self):
        return self.edges


class RegionMatrixTable(RegionPairs):
    """
    Class for working with matrix-based data.

    Generally, a RegionMatrix object has two components:

    - Nodes or regions: (Non-overlapping) genomic regions
      obtained by splitting the genome into distinct pieces.
      See also :class:`~GenomicRegion` and :class:`~RegionsTable`

    - Edges or contacts: Pairs of genomic regions with optionally
      associated weight or contact strength. See also
      :class:`~Edge`

    This is a memory-efficient implementation of a matrix data
    container. Internally, this is achieved by saving entries
    of the matrix in sparse notation, i.e. in a list of
    non-zero contacts.

    Its bracket-notation access behaves like a numpy
    array and handles data retrieval and assignment in matrix-
    fashion, e.g. m[1:3] would return rows 1 and 2 of
    the matrix m (0-based index). However, the bracket
    notation can also handle :class:`~GenomicRegion` descriptor
    strings, i.e. m['chr1','chr5'] will extract the inter-
    chromosomal matrix between chromosomes 1 and 5 only.

    Examples:

    .. code:: python

        m = RegionMatrix(file_name="/path/to/save/file")

        # load genomic regions
        genome = Genome.from_folder("/path/to/fasta/folder")
        regions = genome.get_regions("HindIII")
        m.add_regions(regions)

        # load edges
        edges = []
        edges.append(Edge(source=10, sink=23, weight=3)
        edges.append(Edge(source=8, sink=9, weight=57)
        # ...
        m.add_edges(edges)
    """

    _classid = 'REGIONMATRIXTABLE'

    def __init__(self, file_name=None, mode='a', additional_fields=None, tmpdir=None,
                 default_field=None, _table_name_nodes='nodes', _table_name_edges='edges'):

        """
        Initialize a :class:`~RegionMatrixTable` object.

        :param file_name: Path to a save file
        :param mode: File mode to open underlying file
        :param _table_name_nodes: (Internal) name of the HDF5 node for regions
        :param _table_name_edges: (Internal) name of the HDF5 node for edges
        """

        # private variables
        self.default_field = default_field
        RegionPairs.__init__(self, file_name=file_name, mode=mode, additional_fields=additional_fields, tmpdir=tmpdir,
                             _table_name_nodes=_table_name_nodes, _table_name_edges=_table_name_edges)

        if default_field is None:
            for field_name in self._edges.colnames:
                if not field_name.startswith("_") and field_name != "source" and field_name != "sink":
                    self.default_field = field_name
                    break

    def _flush_edge_buffer(self, e_buffer, replace=False, update_index=True,
                           clean_zero=True, default_column=None):
        if default_column is None:
            default_column = self.default_field
        # update current rows
        for row in self._edges:
            key = (row["source"], row["sink"])

            if key in e_buffer:
                value = e_buffer[key]
                # it is a weight
                try:
                    if replace:
                        row[default_column] = float(value)
                    else:
                        row[default_column] += float(value)
                    row.update()
                except TypeError:
                    self.add_edge(value, check_nodes_exist=False, flush=False, replace=replace, row=row)
                del e_buffer[key]
        self.flush(update_index=False)

        # flush remaining buffer
        for source, sink in e_buffer:
            key = (source, sink)
            value = e_buffer[key]
            try:
                v = float(value)
                if v == 0:
                    continue
                new_edge = self._edge_from_dict({'source': source, 'sink': sink, default_column: v})
                self.add_edge(new_edge, check_nodes_exist=False, flush=False)
            except TypeError:
                self.add_edge(value, check_nodes_exist=False, flush=False)

        self.flush(update_index=True)
        if clean_zero:
            self._remove_zero_edges(update_index=update_index, weight_column=default_column)

    def __getitem__(self, key):
        return self.as_matrix(key)

    def as_matrix(self, key=slice(0, None, None), values_from=None, mask_missing=False, impute_missing=False):
        """
        Get a chunk of the matrix.

        :param key: Possible key types are:

                    Region types

                    - Node: Only the ix of this node will be used for
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
        :param values_from: Determines which column will be used to populate
                            the matrix. Default is 'weight'.
        :param mask_missing: if True, will mask missing/unmappable contacts
        :param impute_missing: if True, will average missing contacts
        :return: :class:`RegionMatrix`
        """
        if values_from is None:
            values_from = self.default_field

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

        row_ranges = list(self._get_node_ix_ranges(nodes_ix_row))
        col_ranges = list(self._get_node_ix_ranges(nodes_ix_col))

        m = self._get_matrix(row_ranges, col_ranges, weight_column=values_from)

        # select the correct output format
        # empty result: matrix
        rm = None
        if m.shape[0] == 0 and m.shape[1] == 0:
            rm = RegionMatrix(m, col_regions=[], row_regions=[])
        # both selectors are lists: matrix
        elif isinstance(nodes_ix_row, list) and isinstance(nodes_ix_col, list):
            rm = RegionMatrix(m, col_regions=nodes_col, row_regions=nodes_row)
        # row selector is list: vector
        elif isinstance(nodes_ix_row, list):
            rm = RegionMatrix(m[:, 0], col_regions=[nodes_col], row_regions=nodes_row)
        # column selector is list: vector
        elif isinstance(nodes_ix_col, list):
            rm = RegionMatrix(m[0, :], col_regions=nodes_col, row_regions=[nodes_row])

        if rm is not None:
            if mask_missing or impute_missing:
                mappable = self.mappable()
                mask = np.zeros(m.shape, dtype=bool)
                current_row = 0
                for row_range in row_ranges:
                    for i in xrange(row_range[0], row_range[1]+1):
                        if not mappable[i]:
                            mask[current_row] = True
                        current_row += 1

                current_col = 0
                for col_range in col_ranges:
                    for i in xrange(col_range[0], col_range[1]+1):
                        if not mappable[i]:
                            mask[:, current_col] = True
                        current_col += 1
                masked_rm = np.ma.MaskedArray(rm, mask=mask)

                if impute_missing:
                    return self._impute_missing_contacts(masked_rm)
                return masked_rm
            else:
                return rm

        # both must be indexes
        return m[0, 0]

    def _impute_missing_contacts(self, hic_matrix=None, _expected_contacts=None):
        """
        Impute missing contacts in a Hi-C matrix.

        :param hic_matrix: a :class:`~HicMatrix` object
        :param _expected_contacts: An ExpectedContacts object for this matrix
        :return: the input matrix with imputed values
        """
        if not hasattr(hic_matrix, "mask"):
            raise ValueError("hic_matrix must be a numpy masked array!")

        # here to avoid circular dependency
        from kaic.architecture.hic_architecture import ExpectedContacts

        if _expected_contacts is not None:
            ex = _expected_contacts
            close_ex = False
        else:
            ex = ExpectedContacts(self, smooth=True)
            close_ex = True

        intra_expected = ex.intra_expected()
        inter_expected = ex.inter_expected()

        for i in xrange(hic_matrix.shape[0]):
            row_region = hic_matrix.row_regions[i]
            for j in xrange(hic_matrix.shape[1]):
                col_region = hic_matrix.col_regions[j]

                if hic_matrix.mask[i, j]:
                    if row_region.chromosome == col_region.chromosome:
                        d = abs(row_region.ix-col_region.ix)
                        hic_matrix[i, j] = intra_expected[d]
                    else:
                        hic_matrix[i, j] = inter_expected

        if close_ex:
            ex.close()

        return hic_matrix

    def _get_matrix(self, row_ranges, col_ranges, weight_column=None):
        if weight_column is None:
            weight_column = self.default_field

        n_rows = 0
        for row_range in row_ranges:
            n_rows += row_range[1]-row_range[0]+1

        n_cols = 0
        for col_range in col_ranges:
            n_cols += col_range[1]-col_range[0]+1

        # create empty matrix
        m = np.zeros((n_rows, n_cols))

        # fill matrix with weights
        row_offset = 0
        for row_range in row_ranges:
            n_rows_sub = row_range[1] - row_range[0] + 1
            col_offset = 0
            for col_range in col_ranges:
                n_cols_sub = col_range[1] - col_range[0] + 1

                for edge_row in self._edge_row_range(row_range[0], row_range[1], col_range[0], col_range[1]):
                    source = edge_row['source']
                    sink = edge_row['sink']
                    weight = edge_row[weight_column]

                    ir = source - row_range[0] + row_offset
                    jr = sink - col_range[0] + col_offset
                    if 0 <= ir < m.shape[0] and 0 <= jr < m.shape[1]:
                        m[ir, jr] = weight

                    ir = sink - row_range[0] + row_offset
                    jr = source - col_range[0] + col_offset
                    if 0 <= ir < m.shape[0] and 0 <= jr < m.shape[1]:
                        m[ir, jr] = weight

                col_offset += n_cols_sub
            row_offset += n_rows_sub

        return m

    def as_data_frame(self, key, weight_column=None):
        """
        Get a pandas data frame by key.

        For key types see :func:`~RegionMatrixTable.__getitem__`.

        :param key: For key types see :func:`~RegionMatrixTable.__getitem__`.
        :param weight_column: Determines which column populates the DF
        :return: Pandas data frame, row and column labels are
                 corresponding node start positions
        """
        if weight_column is None:
            weight_column = self.default_field

        nodes_ix_row, nodes_ix_col = self._get_nodes_from_key(key, as_index=True)
        nodes_row, nodes_col = self._get_nodes_from_key(key, as_index=False)

        row_ranges = list(self._get_node_ix_ranges(nodes_ix_row))
        col_ranges = list(self._get_node_ix_ranges(nodes_ix_col))

        m = self._get_matrix(row_ranges, col_ranges, weight_column=weight_column)
        labels_row = []
        for node in nodes_row:
            labels_row.append(node.start)
        labels_col = []
        for node in nodes_col:
            labels_col.append(node.start)
        df = p.DataFrame(m, index=labels_row, columns=labels_col)

        return df

    def __setitem__(self, key, item):
        self.set_matrix(key, item, clean_zero=True)

    def set_matrix(self, key, item, clean_zero=False, values_to=None):
        """
        Set a chunk of the matrix.

        :param key: Possible key types are:

                    Region types

                    - Node: Only the ix of this node will be used for
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
                    e.g.: 'chr1:0-1000, chr4:2300-3000' will set the entries
                    of the relevant regions between chromosomes 1 and 4.
        :param item: matrix to replace existing values
        :param clean_zero: Remove edges where 'values_to' colum is zero
        :param values_to: Determines which column is replaced by the provided
                          item. Default: weight
        """
        if values_to is None:
            values_to = self.default_field

        nodes_ix_row, nodes_ix_col = self._get_nodes_from_key(key, as_index=True)
        self._set_matrix(item, nodes_ix_row, nodes_ix_col, clean_zero=clean_zero, weight_column=values_to)

    def _set_matrix(self, item, nodes_ix_row=None, nodes_ix_col=None,
                    clean_zero=False, weight_column=None):
        if weight_column is None:
            weight_column = self.default_field

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

        self._flush_edge_buffer(replacement_edges, replace=True, clean_zero=clean_zero, default_column=weight_column)

    def _update_edge_weight(self, source, sink, weight, add=False, flush=True, weight_column=None):
        if weight_column is None:
            weight_column = self.default_field

        if source > sink:
            source, sink = sink, source

        value_set = False
        for row in self._edges.where("(source == %d) & (sink == %d)" % (source, sink)):
            original = 0
            if add:
                original = row[weight_column]
            row[weight_column] = weight + original
            row.update()
            value_set = True
            if flush:
                self.flush()
        if not value_set:
            self.add_edge(Edge(source=source, sink=sink, weight=weight), flush=flush)

    def _remove_zero_edges(self, flush=True, update_index=True, weight_column=None):
        if weight_column is None:
            weight_column = self.default_field

        zero_edge_ix = []
        ix = 0
        for row in self._edges.iterrows():
            if row[weight_column] == 0:
                zero_edge_ix.append(ix)
            ix += 1

        for ix in reversed(zero_edge_ix):
            self._edges.remove_row(ix)

        if flush:
            self.flush(update_index=update_index)

    def marginals(self, weight_column=None):
        """
        Get the marginals vector of this Hic matrix.
        """
        if weight_column is None:
            weight_column = self.default_field

        # prepare marginals dict
        marginals = np.zeros(len(self.regions()), float)

        for edge in self.edges(lazy=True):
            marginals[edge.source] += getattr(edge, weight_column)
            if edge.source != edge.sink:
                marginals[edge.sink] += getattr(edge, weight_column)

        return marginals

    def mappable(self):
        """
        Get the mappability vector of this matrix.
        """
        # prepare marginals dict
        mappable = np.zeros(len(self.regions()), dtype=bool)

        for edge in self.edges(lazy=True):
            mappable[edge.source] = True
            if edge.source != edge.sink:
                mappable[edge.sink] = True

        return mappable

    def scaling_factor(self, matrix, weight_column=None):
        """
        Compute the scaling factor to another matrix.

        Calculates the ratio between the number of contacts in
        this Hic object to the number of contacts in another
        Hic object.

        :param matrix: A :class:`~Hic` object
        :param weight_column: Name of the column to calculate the scaling factor on
        :return: float
        """
        if weight_column is None:
            weight_column = self.default_field

        logging.info("Calculating scaling factor...")
        m1_sum = 0.0
        for edge in self.edges(lazy=True):
            m1_sum += getattr(edge, weight_column)
        m2_sum = 0.0
        for edge in matrix.edges(lazy=True):
            m2_sum += getattr(edge, weight_column)
        scaling_factor = m1_sum / m2_sum
        logging.debug("Scaling factor: %f" % scaling_factor)
        return scaling_factor

    def get_combined_matrix(self, matrix, key=None, scaling_factor=None, weight_column=None):
        """
        Return a :class:`~HicMatrix` where values above the diagonal
        are from this object and values below the diagonal are from
        another :class:`~Hic` object.

        "Above the diagonal" refers to the diagonal of the complete
        Hic object, not the diagonal of the returned matrix.

        :param matrix: Another :class:`~Hic` object
        :param key: A matrix selector. Use tuple to selct row and
                    columns, also see __getitem__
        :param scaling_factor: Factor to scale the hic values. If None,
                               will be computed using
                               :func:`~Hic.scaling_factor`.
        :param weight_column: Name of the column used to combine matrices into one
        :return: :class:`~HicMatrix`
        """
        if key is None:
            key = slice(0, None, None)

        if scaling_factor is None:
            scaling_factor = self.scaling_factor(matrix, weight_column=weight_column)

        m_top = self.as_matrix(key=key, values_from=weight_column)

        # find diagonal
        row_region = m_top.row_regions[0]
        matching_index = None
        for i, col_region in enumerate(m_top.col_regions):
            if col_region == row_region:
                matching_index = i

        if matching_index is None:
            col_region = m_top.col_regions[0]
            for i, row_region in enumerate(m_top.row_regions):
                if col_region == row_region:
                    matching_index = -1 * i

        if matching_index is None:
            return m_top

        # replace diagonal
        m_bottom = matrix.as_matrix(key=key, values_from=weight_column) * scaling_factor
        top_indices = np.triu_indices(m_top.shape[0], matching_index, m_top.shape[1])
        m_bottom[top_indices] = m_top[top_indices]
        return m_bottom


class AccessOptimisedRegionMatrixTable(RegionMatrixTable, AccessOptimisedRegionPairs):
    """
    Class with faster access to matrix data, based on :class:`~AccessOptimisedRegionPairs`.
    """

    _classid = 'ACCESSOPTIMISEDREGIONMATRIXTABLE'

    def __init__(self, file_name=None, mode='a', tmpdir=None, additional_fields=None,
                 default_field=None, _table_name_nodes='nodes', _table_name_edges='edges'):
        """
        Initialize a :class:`~AccessOptimisedRegionMatrixTable` object.

        :param file_name: Path to a save file
        :param mode: File mode to open underlying file
        :param _table_name_nodes: (Internal) name of the HDF5 node for regions
        :param _table_name_edges: (Internal) name of the HDF5 node for edges
        """

        # private variables
        self.default_field = default_field
        AccessOptimisedRegionPairs.__init__(self, file_name=file_name, mode=mode, additional_fields=additional_fields,
                                            tmpdir=tmpdir, _table_name_nodes=_table_name_nodes,
                                            _table_name_edges=_table_name_edges)

        if default_field is None:
            self.default_field = self.field_names[2]

    def _flush_edge_buffer(self, e_buffer, replace=False, update_index=True,
                           clean_zero=True, default_column=None):
        if default_column is None:
            default_column = self.default_field

        # re-arrange edge buffer
        partition_e_buffer = defaultdict(dict)
        for key, edge in e_buffer.iteritems():
            source_partition = self._get_partition_ix(key[0])
            sink_partition = self._get_partition_ix(key[1])
            partition_e_buffer[(source_partition, sink_partition)][key] = e_buffer[key]

        # update current rows
        for partition_key, e_buffer in partition_e_buffer.iteritems():
            if partition_key in self._edge_table_dict:
                edge_table = self._edge_table_dict[partition_key]

                for row in edge_table:
                    key = (row["source"], row["sink"])

                    if key in e_buffer:
                        value = e_buffer[key]
                        # it is a weight
                        try:
                            if replace:
                                row[default_column] = float(value)
                            else:
                                row[default_column] += float(value)
                            row.update()
                        except TypeError:
                            self.add_edge(value, check_nodes_exist=False, flush=False, replace=replace, row=row)
                        del e_buffer[key]
                self.flush(update_index=False)

            # flush remaining buffer
            for source, sink in e_buffer.iterkeys():
                key = (source, sink)
                value = e_buffer[key]
                try:
                    v = float(value)
                    if v == 0:
                        continue
                    new_edge = self._edge_from_dict({'source': source, 'sink': sink, default_column: v})
                    self.add_edge(new_edge, check_nodes_exist=False, flush=False)
                except TypeError:
                    self.add_edge(value, check_nodes_exist=False, flush=False)

        self.flush(update_index=True)
        if clean_zero:
            self._remove_zero_edges(update_index=update_index, weight_column=default_column)

    def _remove_zero_edges(self, flush=True, update_index=True, weight_column=None):
        if weight_column is None:
            weight_column = self.default_field

        for edge_table in self._edge_table_iter():
            zero_edge_ix = []
            ix = 0
            for row in edge_table.iterrows():
                if row[weight_column] == 0:
                    zero_edge_ix.append(ix)
                ix += 1

            for ix in reversed(zero_edge_ix):
                edge_table.remove_row(ix)

        if flush:
            self.flush(update_index=update_index)


class Hic(RegionMatrixTable):
    """
    Class for working with Hi-C data.

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

    _classid = 'HIC'

    class HicRegionAnnotationDescription(t.IsDescription):
        bias = t.Float32Col(pos=0, dflt=1)

    def __init__(self, data=None, file_name=None, mode='a', tmpdir=None,
                 _table_name_nodes='nodes', _table_name_edges='edges',
                 _table_name_node_annotations='node_annot'):

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

        RegionMatrixTable.__init__(self, additional_fields={'weight': t.Float64Col(pos=0)},
                                   file_name=file_name, mode=mode, tmpdir=tmpdir,
                                   default_field='weight',
                                   _table_name_nodes=_table_name_nodes,
                                   _table_name_edges=_table_name_edges)

        if _table_name_node_annotations in self.file.root:
            self._node_annotations = self.file.get_node('/', _table_name_node_annotations)
        elif mode not in ('r', 'r+'):
            self._node_annotations = t.Table(self.file.root, _table_name_node_annotations,
                                             Hic.HicRegionAnnotationDescription)
            self._node_annotations.flush()
        else:
            # compatibility with existing objects
            self._node_annotations = None

        # add data
        self._add_data(data)

    def _add_data(self, data):
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

    def load_read_fragment_pairs(self, pairs, excluded_filters=(), _max_buffer_size=5000000):
        """
        Load data from :class:`~kaic.construct.seq.FragmentMappedReadPairs`.

        This method automatically sums up reads mapping to the same
        fragment pairs and creates exactly one edge per fragment pair.

        :param pairs: A :class:`~kaic.construct.seq.FragmentMappedReadPairs`
                      object.
        :param excluded_filters: Filters to ignore when loading the data
        :param _max_buffer_size: Number of edges kept in buffer before
                                 writing to Table.
        """
        # add regions
        if len(self._regions) != 0:
            raise RuntimeError("When importing from read pairs you MUST start from an empty data set!")

        self.add_regions(pairs.regions())

        l = len(pairs)

        with RareUpdateProgressBar(max_value=l) as pb:
            edge_buffer = {}
            for i, pair in enumerate(pairs.pairs(lazy=True, excluded_filters=excluded_filters)):
                source = pair.left.fragment.ix
                sink = pair.right.fragment.ix
                if source > sink:
                    source, sink = sink, source
                key = (source, sink)
                if key not in edge_buffer:
                    edge_buffer[key] = 0
                edge_buffer[key] += 1

                if len(edge_buffer) > _max_buffer_size:
                    logging.info("Flushing buffer")
                    self._flush_edge_buffer(edge_buffer, replace=False, update_index=False)
                    edge_buffer = {}
                try:
                    pb.update(i)
                except ValueError:
                    pass
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
            with RareUpdateProgressBar(max_value=len(hic.regions)) as pb:
                for i, region in enumerate(hic.regions()):
                    self.add_region(region, flush=False)
                    pb.update(i)
            self.flush()

            with RareUpdateProgressBar(max_value=len(hic.edges)) as pb:
                for i, edge in enumerate(hic.edges()):
                    self.add_edge(edge, check_nodes_exist=False, flush=False)
                    pb.update(i)
            self.flush()
            self.bias_vector(hic.bias_vector())
        # if already have nodes in this HiC object...
        else:
            logging.info("Binning Hi-C contacts")
            # create region "overlap map"
            overlap_map = _get_overlap_map(hic.regions(), self.regions())

            edge_buffer = defaultdict(int)
            with RareUpdateProgressBar(max_value=len(hic.edges)) as pb:
                for i, old_edge in enumerate(hic.edges()):
                    old_source = old_edge.source
                    old_sink = old_edge.sink
                    old_weight = old_edge.weight
                    new_edges = _edges_by_overlap_method([old_source, old_sink, old_weight], overlap_map)

                    for new_edge in new_edges:
                        key_pair = (new_edge[0], new_edge[1])
                        edge_buffer[key_pair] += new_edge[2]

                    if len(edge_buffer) > _edge_buffer_size:
                        self._flush_edge_buffer(edge_buffer, replace=False, update_index=False)
                        edge_buffer = defaultdict(int)
                    pb.update(i)
            self._flush_edge_buffer(edge_buffer)

    def copy(self, file_name, tmpdir=None):
        return Hic(data=self, file_name=file_name, tmpdir=tmpdir, mode='w')

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
            chromosome_list.append(Chromosome(name=chromosome, length=self.chromosome_lens[chromosome]))

        genome = Genome(chromosomes=chromosome_list)
        hic = self.__class__(file_name=file_name, mode='w')
        regions = genome.get_regions(bin_size)
        hic.add_regions(regions)
        regions.close()
        hic.load_from_hic(self)

        return hic

    @classmethod
    def from_hic(cls, hics, file_name=None, tmpdir=None, only_intrachromosomal=False):
        if isinstance(hics, Hic):
            hics = [hics]

        logging.info("Checking if regions are identical")
        identical = True
        for i in xrange(1, len(hics)):
            if len(hics[i].regions) != len(hics[0].regions):
                identical = False
                break

            for self_region, hic_region in zip(hics[0].regions, hics[i].regions):
                if self_region.chromosome != hic_region.chromosome:
                    identical = False
                    break
                if self_region.start != hic_region.start:
                    identical = False
                    break
                if self_region.end != hic_region.end:
                    identical = False
                    break

        if not identical:
            raise ValueError("Regions must be identical in both Hic objects to merge!")

        merged_hic = cls(file_name=file_name, tmpdir=tmpdir, mode='w')
        for region in hics[0].regions:
            merged_hic.add_region(region, flush=False)
        merged_hic.flush()

        chromosomes = hics[0].chromosomes()
        for i in xrange(len(chromosomes)):
            r2 = xrange(i, i + 1) if only_intrachromosomal else xrange(i, len(chromosomes))
            for j in r2:
                logging.info("Chromosomes: {}-{}".format(chromosomes[i], chromosomes[j]))
                edges = dict()
                for hic in hics:
                    for edge in hic.edge_subset(key=(chromosomes[i], chromosomes[j])):
                        key = (edge.source, edge.sink)
                        if key not in edges:
                            edges[key] = edge
                        else:
                            edges[key].weight += edge.weight

                for edge in edges.itervalues():
                    merged_hic.add_edge(edge, flush=False)
        merged_hic.flush()

        return merged_hic


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
        chrms = {hl.genome.chrmStartsBinCont[i]: hl.genome.chrmLabels[i] for i in xrange(0, len(hl.genome.chrmLabels))}
        chromosome = ''
        for i in xrange(0, len(hl.genome.posBinCont)):
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

            for i in xrange(0, data.shape[0]):
                iNode = i+chr1StartBin
                start = i
                if chr1 != chr2:
                    start = 0
                for j in xrange(start, data.shape[1]):
                    jNode = j+chr2StartBin

                    if data[i, j] != 0:
                        hic.add_edge([iNode, jNode, data[i, j]], flush=False)
        hic.flush(flush_nodes=False)

        return hic

    def _merge(self, hic, _edge_buffer_size=5000000):
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

            l = len(hic.regions)

            with RareUpdateProgressBar(max_value=l) as pb:
                for i, region in enumerate(hic.regions):
                    ix = self._get_region_ix(region)
                    if ix is None:
                        ix = self.add_region([region.chromosome, region.start, region.end], flush=False)
                    ix_conversion[region.ix] = ix
                    pb.update(i)
                self._regions.flush()

        # merge edges
        self.log_info("Merging contacts...")
        edge_buffer = {}
        l = len(hic.edges)
        with RareUpdateProgressBar(max_value=l) as pb:
            for i, merge_edge in enumerate(hic.edges):
                merge_source = ix_conversion[merge_edge.source]
                merge_sink = ix_conversion[merge_edge.sink]
                merge_weight = merge_edge.weight

                if merge_source > merge_sink:
                    merge_source, merge_sink = merge_sink, merge_source

                edge_buffer[(merge_source, merge_sink)] = merge_weight

                pb.update(i)

                if len(edge_buffer) > _edge_buffer_size:
                    logging.info("Flushing buffer...")
                    self._flush_edge_buffer(edge_buffer, replace=False, update_index=False)
                    edge_buffer = {}
        logging.info("Final flush...")
        self._flush_edge_buffer(edge_buffer, replace=False, update_index=False)

    def merge(self, hic_or_hics, _edge_buffer_size=5000000):
        """
        Merge this object with other :class:`~Hic` objects.

        First merges genomic regions, then merges edges.
        It is strongly advised that the genomic regions in
        both objects are the same, although this method will attempt to
        "translate" regions from one object to the other if
        this is not the case.

        :param hic_or_hics: :class:`~Hic` object or a list
                            of :class:`~Hic` objects to be
                            merged into this one
        """
        import traceback
        if isinstance(hic_or_hics, Hic):
            hic = hic_or_hics
            try:
                self._merge(hic, _edge_buffer_size=_edge_buffer_size)
            except Exception as e:
                hic.__exit__(e, e.message, traceback.format_exc())
            else:
                hic.__exit__(None, None, None)
        else:
            try:
                for hic in hic_or_hics:
                    logging.info("Merging {}".format(hic.file_name))
                    try:
                        self._merge(hic, _edge_buffer_size=_edge_buffer_size)
                    except Exception as e:
                        hic.__exit__(e, e.message, traceback.format_exc())
                    else:
                        hic.__exit__(None, None, None)
            except TypeError:
                logging.info('{} is not a Hic object or an iterable'.format(hic_or_hics))

        logging.info("Removing zero edges")
        self._remove_zero_edges(update_index=True)

    def flush(self, flush_nodes=True, flush_edges=True, update_index=True):
        """
        Write data to file and flush buffers.

        :param flush_nodes: Flush nodes tables
        :param flush_edges: Flush edges table
        :param update_index: Update mask indices in edges table
        """
        RegionMatrixTable.flush(self, flush_nodes=flush_nodes, flush_edges=flush_edges, update_index=update_index)
        self._node_annotations.flush()

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

    def filter_low_coverage_regions(self, rel_cutoff=None, cutoff=None, queue=False):
        """
        Convenience function that applies a :class:`~LowCoverageFilter`.

        The cutoff can be provided in two ways:
        1. As an absolute threshold. Regions with contact count below this
        absolute threshold are filtered
        2. As a fraction relative to the median contact count of all regions.

        If both is supplied, whichever threshold is lower will be selected.

        If no parameter is supplied, rel_cutoff will be chosen as 0.1.

        :param rel_cutoff: A cutoff as a fraction (0-1) of the median contact count of all
                           regions.
        :param cutoff: A cutoff in absolute contact counts (can be float) below
                       which regions are considered "low coverage"
        :param queue: If True, filter will be queued and can be executed
                      along with other queued filters using
                      run_queued_filters
        """
        if cutoff is None and rel_cutoff is None:
            rel_cutoff = 0.1

        mask = self.add_mask_description('low_coverage',
            'Mask low coverage regions in the Hic matrix '
            '(absolute cutoff {:.4}, relative cutoff {:.1%}'.format(float(cutoff) if cutoff else 0., float(rel_cutoff) if rel_cutoff else 0.))

        low_coverage_filter = LowCoverageFilter(self, rel_cutoff=rel_cutoff, cutoff=cutoff, mask=mask)
        self.filter(low_coverage_filter, queue)
    
    def bias_vector(self, vector=None):
        """
        Get or set the bias vector of this Hic matrix.

        :param vector: Numpy float vector. If provided, sets the
                       the bias vector of the object.
        """

        if vector is not None:
            if len(vector) != len(self.regions()):
                raise ValueError("Bias vector must be the same length as number of regions (%d)" % len(self.regions()))

            # overwrite biases
            if len(self._node_annotations) == len(vector):
                for i, row in enumerate(self._node_annotations):
                    row['bias'] = vector[i]
                    row.update()
            # create new biases
            else:
                row = self._node_annotations.row
                for value in vector:
                    row['bias'] = value
                    row.append()
            self._node_annotations.flush()
            return vector

        vector = np.ones(len(self.regions()))
        if len(self._node_annotations) > 0:
            for i, row in enumerate(self._node_annotations):
                vector[i] = row['bias']

        return vector

    def mappable_regions(self):
        marginals = self.marginals()
        mappable = defaultdict(int)
        for i, region in enumerate(self.regions()):
            if marginals[i] > 0:
                mappable[region.chromosome] += 1
        return mappable

    def possible_contacts(self, _mappable=None):
        if _mappable is None:
            _mappable = self.mappable_regions()

        # calculate possible combinations
        intra_possible = 0
        inter_possible = 0
        chromosomes = _mappable.keys()
        for i in xrange(len(chromosomes)):
            chromosome1 = chromosomes[i]
            n1 = _mappable[chromosome1]
            intra_possible += n1**2/2 + n1/2
            for j in xrange(i+1, len(chromosomes)):
                chromosome2 = chromosomes[j]
                n2 = _mappable[chromosome2]
                inter_possible += n1*n2

        return intra_possible, inter_possible

    @property
    def architecture(self):
        import kaic.architecture.hic_architecture as ha
        return ha.HicArchitecture(self)


class AccessOptimisedHic(Hic, AccessOptimisedRegionMatrixTable):
    """
    Class with faster access to matrix data, based on :class:`~AccessOptimisedRegionPairs`.
    """

    _classid = 'ACCESSOPTIMISEDHIC'

    def __init__(self, data=None, file_name=None, mode='a', tmpdir=None,
                 _table_name_nodes='nodes', _table_name_edges='edges',
                 _table_name_node_annotations='node_annot'):
        """
        Initialize a :class:`~AccessOptimisedHic` object.

        :param data: Can be the path to an XML file denoting a Hic object,
                     another Hic object, a :class:`~FragmentMappedReadPairs`
                     object, or a path to a save file. In the latter case,
                     this parameter may replace file_name, but only if
                     file_name is None.
        :param file_name: Path to a save file
        :param _table_name_nodes: (Internal) name of the HDF5 node for regions
        :param _table_name_edges: (Internal) name of the HDF5 node for edges
        """

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

        AccessOptimisedRegionMatrixTable.__init__(self, additional_fields={'weight': t.Float64Col(pos=0)},
                                                  file_name=file_name, mode=mode, tmpdir=tmpdir,
                                                  default_field='weight',
                                                  _table_name_nodes=_table_name_nodes,
                                                  _table_name_edges=_table_name_edges)

        if _table_name_node_annotations in self.file.root:
            self._node_annotations = self.file.get_node('/', _table_name_node_annotations)
        elif mode not in ('r', 'r+'):
            self._node_annotations = t.Table(self.file.root, _table_name_node_annotations,
                                             Hic.HicRegionAnnotationDescription)
            self._node_annotations.flush()
        else:
            # compatibility with existing objects
            self._node_annotations = None

        # add data
        self._add_data(data)

    def flush(self, flush_nodes=True, flush_edges=True, update_index=True):
        """
        Write data to file and flush buffers.

        :param flush_nodes: Flush nodes tables
        :param flush_edges: Flush edges table
        :param update_index: Update mask indices in edges table
        """
        AccessOptimisedRegionMatrixTable.flush(self, flush_nodes=flush_nodes,
                                               flush_edges=flush_edges, update_index=update_index)
        self._node_annotations.flush()

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
            for edge_table in self._edge_table_iter():
                edge_table.filter(edge_filter, _logging=log_progress)
        else:
            for edge_table in self._edge_table_iter():
                edge_table.queue_filter(edge_filter)

    def run_queued_filters(self, log_progress=False):
        """
        Run queued filters.

        :param log_progress: If true, process iterating through all edges
                             will be continuously reported.
        """
        for edge_table in self._edge_table_iter():
            edge_table.run_queued_filters(_logging=log_progress)


def load_hic(file_name, mode='r', tmpdir=None, _edge_table_name='edges'):
    f = t.open_file(file_name, mode='r')
    n = f.get_node('/' + _edge_table_name)
    if isinstance(n, MaskedTable):
        hic_class = Hic
    elif isinstance(n, t.group.Group):
        hic_class = AccessOptimisedHic
    else:
        raise ValueError("%s is not a valid Hi-C object file" % file_name)

    f.close()
    return hic_class(file_name=file_name, mode=mode, tmpdir=tmpdir)


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
    chosen at 10% of the mean contact count of all regions.
    """
    def __init__(self, hic_object, cutoff=None, rel_cutoff=None, mask=None):
        """
        Initialize filter with these settings.

        The cutoff can be provided in two ways:
        1. As an absolute threshold. Regions with contact count below this
        absolute threshold are filtered
        2. As a fraction relative to the median contact count of all regions.

        If both is supplied, whichever threshold is lower will be selected.

        :param hic_object: The :class:`~Hic` object that this
                           filter will be called on. Needed for
                           contact count calculation.
        :param rel_cutoff: A cutoff as a fraction (0-1) of the median contact count of all
                           regions.
        :param cutoff: A cutoff in absolute contact counts (can be float) below
                       which regions are considered "low coverage"
        :param mask: Optional Mask object describing the mask
                     that is applied to filtered edges.
        """
        HicEdgeFilter.__init__(self, mask=mask)

        self._marginals = hic_object.marginals()
        if cutoff is None and rel_cutoff is None:
            raise ValueError("Either rel_cutoff or cutoff must be given")
        cutoff = min(cutoff if cutoff else float("inf"),
                     self.calculate_cutoffs(rel_cutoff)[0] if rel_cutoff else float("inf"))
        logging.info("Final absolute cutoff threshold is {:.4}".format(float(cutoff)))

        self._regions_to_mask = set()
        for i, contacts in enumerate(self._marginals):
            if contacts < cutoff:
                self._regions_to_mask.add(i)
        logging.info("Selected a total of {} ({:.1%}) regions to be masked".format(
            len(self._regions_to_mask), len(self._regions_to_mask)/len(hic_object.regions)))

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


class RegionMatrix(np.ndarray):
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
            try:
                col_key = slice(0, len(self.col_regions), 1)
            except TypeError:
                col_key = None
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

        try:
            col_regions = self.col_regions[col_key]
        except TypeError:
            col_regions = None

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
            return slice(start, stop+1, 1)
        return key

    def __reduce__(self):
        # Get the parent's __reduce__ tuple
        pickled_state = super(RegionMatrix, self).__reduce__()
        # Create our own tuple to pass to __setstate__
        new_state = pickled_state[2] + (pickle.dumps(self.row_regions), pickle.dumps(self.col_regions))
        # Return a tuple that replaces the parent's __setstate__ tuple with our own
        return pickled_state[0], pickled_state[1], new_state

    def __setstate__(self, state):
        self.row_regions = pickle.loads(state[-2])
        self.col_regions = pickle.loads(state[-1])
        # Call the parent's __setstate__ with the other tuple elements.
        super(RegionMatrix, self).__setstate__(state[0:-2])


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
                return Node(ix=ix, chromosome=chromosome, start=start, end=end)
            
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
                return Edge(source=source, sink=sink, weight=weight)
            
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
