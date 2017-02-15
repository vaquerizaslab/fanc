import itertools
import random
import collections
import progressbar
import tables as t
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


class RareUpdateProgressBar(progressbar.ProgressBar):
    def __init__(self, min_value=0, max_value=None, widgets=None,
                 left_justify=True, initial_value=0, poll_interval=None,
                 percent_update_interval=1, silent=False,
                 **kwargs):
        self.manual_poll_interval = True
        if poll_interval is None:
            self.manual_poll_interval = False

        if silent:
            kwargs['fd'] = open(os.devnull, "w")

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
