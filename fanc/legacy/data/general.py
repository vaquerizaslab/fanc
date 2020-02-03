"""
Provide basic convenience data types.

Most classes are based on pytables and hdf5 dictionaries,
allowing on-disk storage and therefore processing of large
files. Other features include indexing and querying.
"""

from __future__ import division
from fanc.config import config
import tables as t
from tables.nodes import filenode
from fanc.tools.files import create_or_open_pytables_file, tmp_file_name
from fanc.tools.general import RareUpdateProgressBar, create_col_index
import os
from tables.exceptions import NoSuchNodeError
from abc import ABCMeta, abstractmethod
import shutil
import warnings
from collections import defaultdict
from .registry import class_id_dict, class_name_dict
import tempfile
from future.utils import with_metaclass, string_types
from builtins import object
import logging
logger = logging.getLogger(__name__)

_filter = t.Filters(complib="blosc", complevel=2, shuffle=True)


class MetaFileBased(type):
    """
    Metaclass that exists to register classes by name in the module.
    Taken pretty much verbatim from the same system in PyTables.
    """
    def __init__(cls, name, bases, dict_):
        super(MetaFileBased, cls).__init__(name, bases, dict_)

        # Always register into class name dictionary.
        class_name_dict[cls.__name__] = cls

        # Register into class identifier dictionary only if the class
        # has an identifier and it is different from its parents'.
        cid = getattr(cls, '_classid', None)
        if cid is not None:
            for base in bases:
                pcid = getattr(base, '_classid', None)
                if pcid == cid:
                    break
            else:
                class_id_dict[cid] = cls


class FileBased(with_metaclass(MetaFileBased, object)):
    _classid = 'FILEBASED'

    def __init__(self, file_name=None, mode='a', tmpdir=None,
                 _meta_group='meta_information'):

        if hasattr(self, 'tmp_file_name'):
            self.tmp_file_name = getattr(self, 'tmp_file_name')
        else:
            self.tmp_file_name = None

        # open file or keep in memory
        if hasattr(self, 'file'):
            if not isinstance(self.file, t.file.File):
                raise ValueError("'file' attribute already exists, but is no pytables File")
        else:
            self.file = None
            self.file_name = file_name

            if tmpdir is None or (isinstance(tmpdir, bool) and not tmpdir):
                self.tmp_file_name = None
                self._init_file(file_name, mode)
            else:
                logger.info("Working in temporary directory...")
                if isinstance(tmpdir, bool):
                    tmpdir = tempfile.gettempdir()
                else:
                    tmpdir = os.path.expanduser(tmpdir)
                self.tmp_file_name = tmp_file_name(tmpdir, prefix='tmp_fanc', extension='h5')
                logger.info("Temporary output file: {}".format(self.tmp_file_name))
                if mode in ('r+', 'r'):
                    shutil.copyfile(file_name, self.tmp_file_name)
                elif mode == 'a' and file_name is not None and os.path.isfile(file_name):
                    shutil.copyfile(file_name, self.tmp_file_name)
                self._init_file(self.tmp_file_name, mode)

        if not hasattr(self, 'meta'):
            self._meta_group_name = _meta_group
            self.meta = None
            self._init_meta()
            self._update_classid()

    def _init_meta(self):
        class MetaAccess(object):
            def __init__(self, meta_attributes=None):
                self._meta_attributes = meta_attributes

            def __getitem__(self, item):
                if self._meta_attributes is None:
                    raise KeyError("No such key: {}", item)
                return self._meta_attributes[item]

            def __getattr__(self, item):
                if item == '_meta_attributes':
                    return object.__getattribute__(self, item)
                if self._meta_attributes is None:
                    raise AttributeError("No such attribute: {}", item)
                return getattr(self._meta_attributes, item)

            def __setattr__(self, key, value):
                if key == '_meta_attributes':
                    object.__setattr__(self, key, value)
                    return

                if self._meta_attributes is None:
                    raise t.FileModeError("File not writable, attribute cannot be set!")
                self._meta_attributes[key] = value

            def __setitem__(self, key, value):
                if self._meta_attributes is None:
                    raise t.FileModeError("File not writable, item cannot be set!")
                setattr(self._meta_attributes, key, value)

            def __contains__(self, item):
                if self._meta_attributes is None:
                    return False
                try:
                    _ = self._meta_attributes[item]
                    return True
                except KeyError:
                    return False

            def keys(self):
                return vars(self._meta_attributes).keys()

            def values(self):
                return vars(self._meta_attributes).values()

        # existing?
        try:
            meta_group = self.file.get_node('/' + self._meta_group_name)
        except t.NoSuchNodeError:
            try:
                meta_group = self.file.create_group('/', self._meta_group_name)
            except t.FileModeError:
                logger.debug("File not open for writing, not creating meta group.")
                self.meta = MetaAccess()
                return

        try:
            meta_node = self.file.get_node(meta_group, 'meta_node')
        except t.NoSuchNodeError:
            try:
                meta_node = filenode.new_node(self.file, where=meta_group, name='meta_node')
            except t.FileModeError:
                logger.debug("File not open for writing, not creating meta node.")
                self.meta = MetaAccess()
                return

        self.meta = MetaAccess(meta_node.attrs)

    def _update_classid(self, force=False):
        if '_classid' not in self.meta or force:
            try:
                self.meta._classid = self._classid
            except (AttributeError, t.FileModeError):
                pass

    def close(self, copy_tmp=True, remove_tmp=True):
        if not self.file.isopen:
            warnings.warn("File {} is already closed!".format(self.file.filename))
            return
        file_mode = self.file.mode
        self.file.close()
        if self.tmp_file_name is not None:
            if copy_tmp and self.file_name is not None and file_mode not in ('r', 'r+'):
                logger.info("Moving temporary output file to destination {}".format(self.file_name))
                shutil.copyfile(self.tmp_file_name, self.file_name)
            if remove_tmp:
                os.remove(self.tmp_file_name)

    def _init_file(self, file_name, mode):
        if file_name is None:
            self.file = create_or_open_pytables_file()
        elif type(file_name) == str:
            self.file = create_or_open_pytables_file(file_name, mode=mode)
        elif isinstance(file_name, t.file.File):
            self.file = file_name
        else:
            raise TypeError("file_name is not a recognisable type")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        report_exception = exc_type is None
        self.close(copy_tmp=report_exception)
        return report_exception


class FileGroup(FileBased):

    _classid = 'FILEGROUP'

    def __init__(self, group, file_name=None, mode='a', tmpdir=None):
        FileBased.__init__(self, file_name=file_name, mode=mode, tmpdir=tmpdir)

        try:
            group_node = self.file.get_node("/" + group)
            if not isinstance(group_node, t.group.Group):
                raise TypeError("%s is not a group, but %s" % (group, str(type(group_node))))
            self._group = group_node
        except t.NoSuchNodeError:
            self._group = self.file.create_group('/', group)


class TableObject(object):
    def __getitem__(self, key):
        try:
            return self.__getattribute__(key)
        except AttributeError:
            raise IndexError("No item " + str(key) + " in object")


class Mask(object):
    """
    Class providing Mask details.
    """
    
    def __init__(self, name=None, description='', ix=None):
        if ix is not None:
            self.ix = ix
        else:
            self.ix = 0
        self.name = name.decode() if isinstance(name, bytes) else name
        self.description = description.decode() if isinstance(description, bytes) else description

    def __repr__(self):
        return "%d. %s: %s" % (self.ix, self.name, self.description)


class Maskable(FileBased):
    """
    Class that adds masking functionality to tables.
    
    It can be desirable to hide a subset of data in 
    a table rather than deleting it outright. This is 
    especially useful in pytables tables, as they are
    not designed for frequent deletions and may 
    suffer from performance problems as a consequence.
    The process of hiding data is referred to as 
    'masking'.
    
    This class creates a pytables Table that contains
    meta-information about masks, i.e. mask IDs, mask
    names, and mask descriptions. Since it is backed
    by pytables, Mask descriptions will be saved to 
    the corresponding hdf5 file automatically. This 
    class additionally provides tools to create and 
    retrieve Mask descriptions.
    
    Note that this class only provides meta-
    information for Masks. Masking can technically 
    performed without Maskable, but then the reasons
    for masking will not be saved.
    
    The actual masking is performed using a MaskFilter
    on a MaskedTable object using the 
    MaskedTable.filter(MaskFilter) function. 
    """

    _classid = 'MASKABLE'
    
    class MaskDescription(t.IsDescription):
        ix = t.Int16Col(pos=0)
        name = t.StringCol(50, pos=1)
        description = t.StringCol(255, pos=2)

    def __init__(self, data=None, file_name=None, table_name="mask", mode='a', tmpdir=None):
        """
        Enable recording of masking in pytables-backed object.
        
        This constructor is built with considerable flexibility
        in the ways it can be called. The 'data' parameter
        serves here as a multi-purpose container that simplifies
        calling the constructor with just one argument.
        
        Args:
            data (...):         multi-purpose container, simplifies one-argument
                                constructor calls. Can stand in for a pytables 
                                Table or File, or a path.
                (None):         An existing mask description table is expected in 
                                self._mask. If none is found, the constructor will 
                                look for a pytables File in self.file and attempt 
                                to create the mask description table in there.
                (str):          An exisiting pytable File at this location will be
                                opened in append mode or one will be created.
                (tables.file.File):
                                A mask description table will be created in this
                                pytables File.
                (tables.table.Table):
                                This pytables Table will be used as mask 
                                description table (if it has the necessary fields
                                ix, name, description).
            table_name (str):   name of the mask table that is created in
                                pytables file, does not usually need to be 
                                modified
        """
        # parse potential unnamed argument
        if data is not None:
            # data is file name
            if type(data) is str or isinstance(data, t.file.File):
                if file_name is None:
                    file_name = data
                    data = None
            elif type(data) == t.table.Table:
                self._set_mask_table(data)
        else:
            if hasattr(self, '_mask'):
                self._set_mask_table(self._mask)

        FileBased.__init__(self, file_name, tmpdir=tmpdir, mode=mode)
                
        if (not hasattr(self, '_mask') or self._mask is None) and self.file is not None:
            try:
                self._mask = self.file.get_node('/' + table_name)
            except NoSuchNodeError:
                self._mask = self.file.create_table("/", table_name, Maskable.MaskDescription)
                row = self._mask.row
                row['ix'] = 0
                row['name'] = 'default'
                row['description'] = 'Default mask'
                row.append()
                self._mask.flush()
        
    def _set_mask_table(self, table):
        if type(table) == t.table.Table:
            if ('ix' not in table.colnames or
                    'name' not in table.colnames or
                    'description' not in table.colnames):
                raise ValueError("Object already has a mask table, \
                                  but it does not have all the necessary \
                                  columns (ix, name, description)")
            self._mask = table
        else:
            raise ValueError("Table is not a MaskDescription table")

    def add_mask_description(self, name, description):
        """
        Add a mask description to the _mask table and return its ID.
        

        Args:
            name (str):        name of the mask
            description (str): description of the mask
            
        Returns:
            int:    id of the mask
        """
        ix = len(self._mask)
        row = self._mask.row
        row['ix'] = ix
        row['name'] = name
        row['description'] = description
        row.append()
        self._mask.flush()

        return Mask(ix=ix, name=name, description=description)

    def masks(self):
        this = self

        class MaskIter(object):
            def __init__(self):
                self.iter = iter(this._mask)

            def __iter__(self):
                return self

            def __next__(self):
                row = next(self.iter)
                return Maskable._row_to_mask(row)

        return MaskIter()

    @staticmethod
    def _row_to_mask(row):
        return Mask(name=row['name'], ix=row['ix'], description=row['description'])

    def get_binary_mask_from_masks(self, masks):
        o = []
        for m in masks:
            if isinstance(m, string_types):
                o.append(2**self.get_mask(m).ix)
            elif isinstance(m, MaskFilter):
                o.append(2**m.mask_ix)
            elif isinstance(m, Mask):
                o.append(2**m.ix)
            elif isinstance(m, int):
                o.append(2**m)
            else:
                raise ValueError('Can only get binary mask from mask names, indexes and MaskFilter instances')
        return sum(o)

    def get_mask(self, key):
        """
        Search _mask table for key and return Mask.
        
        Args:
            key (str): search by mask name
            key (int): search by mask ID
            
        Returns:
            Mask
        """
        
        if isinstance(key, int):
            for row in self._mask.where("ix == %d" % key):
                return Maskable._row_to_mask(row)
        else:
            for row in self._mask.where("name == b'%s'" % str(key)):
                return Maskable._row_to_mask(row)
        return KeyError("Unrecognised key type")
    
    def get_masks(self, ix):
        """
        Extract mask IDs encoded in parameter and return masks.
        
        IDs are powers of 2, so a single int field in the table can hold
        multiple masks by simply adding up the IDs. Similar principle to
        UNIX chmod (although that uses base 8)
        
        Args:
            ix (int): integer that is the sum of powers of 2. Note that this value
                      is not necessarily itself a power of 2.
            
        Returns:
            list (Mask): list of Masks extracted from ix
        """
        
        key_copy = ix

        last_ix = int(self._mask[-1][0])
        
        masks = []
        while last_ix > 0:
            if key_copy - 2**last_ix >= 0:
                mask = self.get_mask(last_ix)
                if mask is not None:
                    masks.append(mask)
                key_copy -= 2**last_ix
            last_ix -= 1
        if key_copy > 0:
            masks.append(self.get_mask(0))
        
        return list(reversed(masks))

    def _mask_statistics_table(self, table, include_unmasked=True):
        masks = defaultdict(int)
        if include_unmasked:
            masks['unmasked'] = 0

        if 'mask_stats' in table.attrs:
            stats = table.attrs['mask_stats']
        else:
            stats = table.mask_stats()

        for mask_bit, count in stats.items():
            row_masks = self.get_masks(mask_bit)

            found_masks = False
            for mask in row_masks:
                found_masks = True
                masks[mask.name] += count

            if not found_masks and include_unmasked:
                masks['unmasked'] += count
        return masks

    def _mask_statistics_group(self, group, include_unmasked=True):
        masks = defaultdict(int)

        for table in group:
            for mask, count in self._mask_statistics_table(table, include_unmasked=include_unmasked).items():
                masks[mask] += count

        return masks

    def mask_statistics(self, table, include_unmasked=True):
        if isinstance(table, t.Group):
            return self._mask_statistics_group(table, include_unmasked=include_unmasked)
        elif isinstance(table, MaskedTable):
            return self._mask_statistics_table(table, include_unmasked=include_unmasked)
        else:
            raise ValueError("First arg must be PyTable Group or MaskedTable!")


class MaskedTableView(object):
    def __init__(self, masked_table, it=None, excluded_masks=0):
        self.masked_table = masked_table
        self.iter = it if it else self.masked_table._iter_visible_and_masked()
        self.excluded_mask_ix = excluded_masks

    def __iter__(self):
        return self

    def __next__(self):
        row = next(self.iter)
        # bit-shift magic! Go @alexis!
        # a is a subset of b if and only if a | b == b.
        # If this condition is satisfied for each byte, return TRUE. Otherwise return FALSE
        while row[self.masked_table._mask_field] | self.excluded_mask_ix != self.excluded_mask_ix:
            row = next(self.iter)
        return row


class MaskedTable(t.Table):
    """
    Wrapper that adds masking functionality to a pytables table. 
    
            
    MaskedTable is simply a wrapper around a
    pytables Table with certain columns. It
    adds masking functionality to table rows.
        
    It can be desirable to hide a subset of data in 
    a table rather than deleting it outright. This is 
    especially useful in pytables tables, as they are
    not designed for frequent deletions and may 
    suffer from performance problems as a consequence.
    
    This class extends the functionality of classes
    that work on the basis of a pytables table. It
    allows the masking of rows and provides table iterators
    and len() functions that ignore masked rows, and 
    therefore give the illusion of deleted rows.
    Masked rows, however, can be inspected separately
    at any point.
    
    The minimum requirement for a table to be masked
    is the inclusion of a 'mask' column. For full
    iterator and len() functionality, the table also
    needs to include a 'mask_ix' column. Checks for
    this are performed when specifying the table to
    be masked.
    """
    
    _c_classid = 'MASKEDTABLE'
    
    def __init__(self, parentnode, name, description=None,
                 title="", filters=None, expectedrows=None,
                 chunkshape=None, byteorder=None, _log=False,
                 mask_field='_mask', mask_index_field='_mask_ix',
                 ignore_reserved_fields=False):
        """
        Pytables Table extension to provide masking functionality.
        """
        
        # set instance variables
        self._queued_filters = []
        self._mask_field = mask_field
        self._mask_index_field = mask_index_field

        if description is not None:
            # try converting description to dict
            try:
                description = description.columns
            except AttributeError:
                pass
            
            # fill in fields required for masking
            if isinstance(description, dict):
                masked_description = description.copy()
            else:
                raise ValueError("Unrecognised description type (%s)" % str(type(description)))
            
            # check that reserved keys are not used
            if mask_field in masked_description:
                if not ignore_reserved_fields:
                    raise ValueError("{} field is reserved in MaskedTable!".format(mask_field))
            else:
                masked_description[mask_field] = t.Int32Col()

            if mask_index_field in masked_description:
                if not ignore_reserved_fields:
                    raise ValueError("{} field is reserved in MaskedTable!".format(mask_index_field))
            else:
                masked_description[mask_index_field] = t.Int64Col()
        else:
            masked_description = None
                
        t.Table.__init__(self, parentnode, name,
                         description=masked_description, title=title,
                         filters=filters,
                         expectedrows=expectedrows,
                         _log=_log,
                         chunkshape=chunkshape,
                         byteorder=byteorder)
        
        self.enable_mask_index()

    def disable_mask_index(self):
        mask_ix_col = getattr(self.cols, self._mask_index_field)
        mask_ix_col.remove_index()

    def enable_mask_index(self):
        mask_ix_col = getattr(self.cols, self._mask_index_field)
        create_col_index(mask_ix_col)

    def flush(self, update_index=False, log_progress=True):
        """
        Flush buffered rows.
        
        Also updates the mask index, if requested.
        """
        self._flush(update_index, log_progress=log_progress)
    
    def _flush(self, update_index=False, log_progress=True):
        # commit any previous changes
        t.Table.flush(self)

        if update_index:
            self._update_ix(log_progress=log_progress)
            # force flush of index if
            # autoindex is disabled
            if not self.autoindex:
                self.flush_rows_to_index()
            # commit index changes
            t.Table.flush(self)

    def iterrows(self, start=None, stop=None, step=None, excluded_masks=0):
        it = t.Table.iterrows(self, start, stop, step)
        return MaskedTableView(self, it, excluded_masks=excluded_masks)

    def itersorted(self, sortby, checkCSI=False,
                   start=None, stop=None, step=None, excluded_masks=0):
        it = t.Table.itersorted(self, sortby, checkCSI=checkCSI,
                                start=start, stop=stop, step=step)
        return MaskedTableView(self, it, excluded_masks=excluded_masks)

    def _iter_visible_and_masked(self):
        """
        Return an iterator over all rows, including masked ones.
        """
        return t.Table.iterrows(self)

    def masked_rows(self):
        """
        Return an iterator over masked rows.
        """
        this = self
        it = self._iter_visible_and_masked()

        class MaskedRows(MaskedTableView):
            def __init__(self, masked_table, it_):
                super(MaskedRows, self).__init__(masked_table, it_)

            def __next__(self):
                row = next(self.iter)
                while row[this._mask_field] == 0:
                    row = next(self.iter)
                return row

            def __getitem__(self, key):
                if isinstance(key, int):
                    if key >= 0:
                        key = -1*key - 1
                        res = [x.fetch_all_fields() for x in
                               super(MaskedTable, this).where("%s == %d" % (this._mask_index_field, key))]
                        if len(res) == 1:
                            return res[0]
                        if len(res) == 0:
                            raise IndexError("Index %d out of bounds" % key)
                        raise RuntimeError("Duplicate row for key %d" % key)
                    else:
                        l = this._visible_len()
                        return self[l+key]
                else:
                    raise KeyError('Cannot retrieve row with key ' + str(key))

        return MaskedRows(self, it)

    def __getitem__(self, key):
        return self._get_visible_item(key)
    
    def _get_visible_item(self, key):
        if type(key) == slice:
            res = []
            # set sensible defaults
            start = key.start
            if start is None:
                start = 0
            stop = key.stop
            if stop is None:
                stop = len(self)
            step = key.step
            if step is None:
                step = 1
                  
            for i in range(start, stop, step):
                res.append(self._get_visible_item(i))
            return res
        else:
            try:
                # this could be a numpy int
                try:
                    key = key.item()
                except AttributeError:
                    pass
                
                key = int(key)
                
                if key >= 0:
                    res = [x.fetch_all_fields() for x in self.where("%s == %d" % (self._mask_index_field, key))]
                    if len(res) == 1:
                        return res[0]
                    if len(res) == 0:
                        raise IndexError("Index %d out of bounds" % key)
                    raise RuntimeError("Duplicate row for key %d" % key)
                else:
                    l = self._visible_len()
                    if l+key >= 0:
                        return self._get_visible_item(l+key)
                    raise KeyError('Cannot retrieve row with key %s' % str(key))
            except ValueError as e:
                raise KeyError('Cannot retrieve row with key %s (%s)' % (str(key), str(e)))
    
    def _original_getitem(self, key):
        return t.Table.__getitem__(self, key)
    
    def __len__(self):
        """
        Return the 'perceived' length of the masked table.
          
        If the table has masked rows, these will not be counted.
        """
          
        return self._visible_len()
    
    def _visible_len(self):
        if 'masked_length' not in self.attrs or self.attrs['masked_length'] == -1:
            return sum(1 for _ in iter(self.where("%s >= 0" % self._mask_index_field)))
        return int(self.attrs['masked_length'])
    
    def _original_len(self):
        return t.Table.__len__(self)
    
    # new index update method
    def _update_ix(self, log_progress=not config.hide_progressbars):
        """
        Update the row indexes of the Table.
        
        Should be run after masking. Will assign auto-
        incrementing integers (from 0) to the mask index
        field of each row in the table if it is not 
        masked, -1 otherwise.
        """

        if log_progress:
            logger.debug("Updating mask indices")

        stats = defaultdict(int)

        l = self._original_len()
        with RareUpdateProgressBar(max_value=l, silent=not log_progress) as pb:
            ix = 0
            masked_ix = -1
            for i, row in enumerate(self._iter_visible_and_masked()):
                m = row[self._mask_field]
                stats[m] += 1

                if m > 0:
                    row[self._mask_index_field] = masked_ix
                    masked_ix -= 1
                else:
                    row[self._mask_index_field] = ix
                    ix += 1
                row.update()
                pb.update(i)
            self.attrs['masked_length'] = ix
        try:
            self.attrs['mask_stats'] = stats
        except t.FileModeError:
            pass

    def mask_stats(self):
        stats = defaultdict(int)
        for row in self._iter_visible_and_masked():
            stats[row[self._mask_field]] += 1

        try:
            self.attrs['mask_stats'] = stats
        except t.FileModeError:
            pass

        return stats
    
    def _get_masks(self, binary_mask):
        def bits(n):
            while n:
                b = n & (~n+1)
                yield b
                n ^= b
        return list(bits(binary_mask))

    def _row_masks(self, row):
        return self._get_masks(row[self._mask_field])

    def _has_mask(self, row, mask):
        return mask in self._row_masks(row)

    def filter(self, mask_filter, _logging=not config.hide_progressbars):
        """
        Run a MaskFilter on this table.
        
        This functions calls the MaskFilter.valid function on
        every row and masks them if the function returns False.
        After running the filter, the table index is updated
        to match only unmasked rows.

        :param mask_filter: A :class:`~MaskFilter` object
        :param _logging: Print progress to stderr
        """

        total = 0
        ix = 0
        mask_ix = -1

        # statistics
        stats = defaultdict(int)
        with RareUpdateProgressBar(max_value=self._original_len(), silent=not _logging) as pb:
            for i, row in enumerate(self._iter_visible_and_masked()):
                total += 1

                if not mask_filter.valid(row):
                    row[self._mask_field] += 2**mask_filter.mask_ix

                stats[row[self._mask_field]] += 1

                # update index
                if row[self._mask_field] > 0:
                    row[self._mask_index_field] = mask_ix
                    mask_ix -= 1
                else:
                    row[self._mask_index_field] = ix
                    ix += 1
                row.update()
                pb.update(i)
        self.attrs['masked_length'] = ix

        if _logging:
            logger.info("Total: %d. Filtered: %d" % (total, -1*(mask_ix-1)))

        try:
            self.attrs['mask_stats'] = stats
        except t.FileModeError:
            pass

        self.flush(update_index=False)

        return stats

    def queue_filter(self, filter_definition):
        """
        Add a MaskFilter to filter queue.
        
        Queued filters can be run at a later time using
        the run_queued_filters function.
        """
        self._queued_filters.append(filter_definition)
        
    def run_queued_filters(self, _logging=not config.hide_progressbars):
        """
        Run queued MaskFilters.
        
        MaskFilters can be queued using the
        queue_filter function.

        :param _logging: If True, prints log to stderr
        """
        ix = 0
        mask_ix = -1
        total = 0

        stats = defaultdict(int)
        with RareUpdateProgressBar(max_value=self._original_len(), silent=not _logging) as pb:
            for i, row in enumerate(self._iter_visible_and_masked()):
                total += 1
                for f in self._queued_filters:
                    if not f.valid(row):
                        row[self._mask_field] += 2**f.mask_ix

                stats[row[self._mask_field]] += 1

                # update index
                if row[self._mask_field] > 0:
                    row[self._mask_index_field] = mask_ix
                    mask_ix -= 1
                else:
                    row[self._mask_index_field] = ix
                    ix += 1
                row.update()
                pb.update(i)

        self.attrs['masked_length'] = ix

        if _logging:
            logger.info("Total: %d. Filtered: %d" % (total, -1*(mask_ix-1)))

        try:
            self.attrs['mask_stats'] = stats
        except t.FileModeError:
            pass

        self.flush(update_index=False)
        return stats

    def where(self, condition, condvars=None,
              start=None, stop=None, step=None):
        condition = "(" + condition + ") & (%s >= 0)" % self._mask_index_field
        return super(MaskedTable, self).where(condition, condvars=None,
                                              start=None, stop=None, step=None)


class MaskFilter(with_metaclass(ABCMeta, object)):
    """
    Abstract class that defines a filter for MaskedTable.
    """

    def __init__(self, mask=None, mask_ix=0, mask_name='default', mask_description="Default mask."):
        """
        Create a MaskFilter.
        
        Sets values for mask index, name, and description
        to be used by MaskedTable.filter function.
        
        Args:
            mask (Mask, int):    if Mask, the attributes mask_ix,
                                 mask_name, and mask_description
                                 will be copied here.
                                 if int, mask_name and mask_description
                                 will be set to defaults, while
                                 mask_ix will be set to this value.
            mask_ix (int):       Index of the mask to be applied.
                                 Defaults to 0
            mask_name (str):     Name of the mask to be applied.
                                 Defaults to 'default'
            mask_description (str):
                                 Description of the mask to be applied.
                                 Defaults to 'Default mask'
        """
        
        if mask is not None:
            if isinstance(mask, Mask):
                self.mask_ix = mask.ix
                self.mask_name = mask.name
                self.mask_description = mask.description
            else:
                self.mask_ix = mask
                self.mask_name = mask_name
                self.mask_description = mask_description
        else:
            self.mask_ix = mask_ix
            self.mask_name = mask_name
            self.mask_description = mask_description
            
        if not type(self.mask_ix) == int:
            raise ValueError("mask_ix must be an integer!")
    
    @abstractmethod
    def valid(self, row):
        """
        Test if a row is valid according to this filter.
        
        Args:
            row (tables.tableextension.Row):
                A pytables Table Row
        
        Returns:
            bool: True if row is valid, False otherwise
        """
        pass
