import itertools
import random
import collections
import progressbar
import tables as t
import re
import os
import errno
from builtins import object
from datetime import datetime


class Map(dict):
    """
    Example:
    m = Map({'first_name': 'Foo'}, last_name='Bar', age=99, sports=['Baz'])
    """
    def __init__(self, *args, **kwargs):
        super(Map, self).__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.items():
                    self[k] = v

        if kwargs:
            for k, v in kwargs.items():
                self[k] = v

    def __getattr__(self, attr):
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(Map, self).__setitem__(key, value)
        self.__dict__.update({key: value})

    def __delattr__(self, item):
        self.__delitem__(item)

    def __delitem__(self, key):
        super(Map, self).__delitem__(key)
        del self.__dict__[key]


class CachedIterator(object):
    def __init__(self, it, cache_size=1):
        self.it = iter(it)
        self.current = None
        self.cache = collections.deque(maxlen=cache_size)
        self.cache_offset = 0

    def __iter__(self):
        return self

    def __next__(self):
        try:
            if self.cache_offset:
                self.cache_offset -= 1
                if self.cache_offset:
                    return self.cache[-self.cache_offset]
                return self.current
            self.cache.append(self.current)
            self.current = next(self.it)
            return self.current
        except StopIteration:
            return None

    def prev(self):
        try:
            self.cache_offset += 1
            return self.cache[-self.cache_offset]
        except IndexError as e:
            self.cache_offset -= 1
            raise e


def ranges(i):
    for _, b in itertools.groupby(enumerate(i), lambda xy: xy[1] - xy[0]):
        b = list(b)
        yield b[0][1], b[-1][1]


def range_overlap(source_start, source_end, sink_start, sink_end):
    overlap = (max(source_start, sink_start), min(source_end, sink_end))
    if overlap[0] <= overlap[1]:
        return overlap
    return None


def bit_flags_from_int(number, base=2):
    # find largest power of 2
    exponent = -1
    while base**(1+exponent) <= number:
        exponent += 1
    
    # number does not contain any power of 2 (0)
    if exponent == -1:
        return set([])
    
    # find all powers
    bits = []
    while exponent >= 0:
        power = base**exponent
        if number - power >= 0:
            number -= power
            bits.append(exponent)
        exponent -= 1
    
    return set(bits)


def distribute_integer(value, divisor, _shuffle=True):
    a = [int(value/divisor)] * divisor
    remaining = value - sum(a)
    i = 0
    while remaining:
        a[i] += 1
        i += 1
        remaining -= 1
    if _shuffle:
        random.shuffle(a)
    return a


def create_col_index(col):
    """
    Create index for Pytables column in a
    safe manner.
    If index already exists, does nothing.
    If index is corrupt, recreates index.
    """
    try:
        col.create_index()
    except ValueError:
        # Index exists
        pass
    except t.NodeError:
        col.remove_index()
        col.create_index()


def to_slice(l):
    if len(l) == 0:
        return slice(0, 0)
    if len(l) == 1:
        return slice(l[0], l[0]+1)
    if len(l) == 2:
        return slice(l[0], l[1]+1, l[1]-l[0])

    d = l[1] - l[0]
    for i in range(len(l)-1):
        if l[i+1] - l[i] != d:
            raise ValueError("Step size between elements varies, cannot make slice")
    return slice(l[0], l[-1]+1, d)


def mkdir(dir_name):
    dir_name = os.path.expanduser(dir_name)

    try:
        os.makedirs(dir_name)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    if not dir_name.endswith('/'):
        dir_name += '/'

    return dir_name


def which(program):
    """
    Check if executable exists in PATH
    :param program: executable name or path to executable
    :return: full path if found, None if not
    """
    import os

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def human_format(num, precision=0):
    """
    Format a number as a string, suffixing letter for 1000 (K), 100000 (M), ...
    :param num: any number larger than zero
    :param precision: number of positions after decimal point
    :return: string representing the number
    """
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    # add more suffixes if you need them
    return '{:.{prec}f}{}'.format(num, ['', 'K', 'M', 'G', 'T', 'P'][magnitude], prec=precision)


def natural_sort(l):
    def convert(text):
        return int(text) if text.isdigit() else text.lower()
    return sorted(l, key=lambda key: [convert(c) for c in re.split('([0-9]+)', key)])


def natural_cmp(pa, pb):
    i = 0
    j = 0

    try:
        while i < len(pa) and j < len(pb):
            if pa[i].isdigit() and pb[j].isdigit():
                # digit comparison
                while pa[i] == '0':
                    i += 1
                while pb[j] == '0':
                    j += 1
                while pa[i].isdigit() and pb[j].isdigit() and pa[i] == pb[j]:
                    i += 1
                    j += 1
                if pa[i].isdigit() and pb[j].isdigit():
                    k = 0
                    try:
                        while pa[i+k].isdigit() and pb[j+k].isdigit():
                            k += 1
                    except IndexError:
                        if i+k < len(pa):
                            return 1
                        if j+k < len(pb):
                            return -1
                        # both ran over index
                        return int(pa[i]) - int(pb[j])
                    return 1 if pa[i+k].isdigit() else -1 if pb[j+k].isdigit() else int(pa[i]) - int(pb[j])
                elif pa[i].isdigit():
                    return 1
                elif pb[j].isdigit():
                    return -1
                elif i != j:
                    return 1 if i < j else -1
            else:
                # string comparison
                if pa[i] != pb[j]:
                    return -1 if pa[i] < pb[j] else 1
                i += 1
                j += 1
    except IndexError:
        pass
    return -1 if i < len(pa) else 1 if j < len(pb) else 0


class RareUpdateProgressBar(progressbar.ProgressBar):
    def __init__(self, min_value=0, max_value=None, widgets=None,
                 left_justify=True, initial_value=0, poll_interval=None,
                 percent_update_interval=1, silent=False, prefix=None,
                 **kwargs):
        self.manual_poll_interval = True
        if poll_interval is None:
            self.manual_poll_interval = False

        if silent:
            kwargs['fd'] = open(os.devnull, "w")

        if widgets is None:
            widgets = [progressbar.widgets.Percentage(), ' (', progressbar.widgets.SimpleProgress(), ') ',
                       progressbar.widgets.Bar(), ' ', progressbar.widgets.Timer(), ' ',
                       progressbar.widgets.AdaptiveETA()]

        if prefix is not None:
            widgets = [prefix + ' '] + widgets

        progressbar.ProgressBar.__init__(self, min_value=min_value, max_value=max_value, widgets=widgets,
                                         left_justify=left_justify, initial_value=initial_value,
                                         poll_interval=poll_interval, **kwargs)

        if self.max_value is not None:
            self.one_percent = self.max_value*(1.0*percent_update_interval)/100
        else:
            self.one_percent = 1

        self.silent = silent

    def _needs_update(self):
        """
        Returns whether the ProgressBar should redraw the line.
        """
        if self.value > self.next_update or self.end_time:
            self.next_update += self.one_percent
            return True

        elif self.manual_poll_interval:
            delta = datetime.now() - self.last_update_time
            return delta > self.poll_interval
