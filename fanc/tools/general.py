import itertools
try:
    from itertools import izip as zip
except ImportError:
    pass
import random
import collections
import progressbar
import tables as t
import re
import os
import errno
from builtins import object
from datetime import datetime
from Bio import Restriction
from Bio.Seq import reverse_complement
from future.utils import string_types
import pysam
import threading
import multiprocessing as mp
import warnings

random.seed(0)

sam_reference_consumers = {0, 2, 3, 7, 8}  # M, D, N, =, X
sam_query_consumers = {0, 1, 4, 7, 8}  # M, I, S, =, X


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


def find_alignment_match_positions(alignment, longest=False):
    reference_pos = alignment.reference_start
    query_pos = 0
    reference_matches = []
    query_matches = []
    longest_match_ix = None
    longest_match_length = 0
    try:
        for operation, n in alignment.cigartuples:
            # longest actual match
            if operation == 0 or operation == 7:
                match_start_ref, match_end_ref = reference_pos, reference_pos + n
                match_start_query, match_end_query = query_pos, query_pos + n
                reference_matches.append((match_start_ref, match_end_ref))
                query_matches.append((match_start_query, match_end_query))
                if longest_match_length < n:
                    longest_match_ix = len(reference_matches) - 1
                    longest_match_length = n

            if operation in sam_reference_consumers:
                reference_pos += n

            if operation in sam_query_consumers:
                query_pos += n

        if longest:
            return reference_matches[longest_match_ix], query_matches[longest_match_ix]
        else:
            return reference_matches, query_matches
    except TypeError:
        if alignment.is_unmapped:
            warnings.warn("Cannot extract match positions for {} because it is unmapped!".format(alignment.qname))
        else:
            warnings.warn("Cannot extract match positions for {} (unknown reason)!".format(alignment.qname))
        return (None, None), (None, None)


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


def mkdir(*dir_name):
    dir_name = os.path.expanduser(os.path.join(*dir_name))

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


def human_format(num, precision=2, lowercase=False):
    """
    Format a number as a string, suffixing letter for 1000 (K), 100000 (M), ...
    :param num: any number larger than zero
    :param precision: number of positions after decimal point
    :param lowercase: return lowercase suffix letters if True
    :return: string representing the number
    """
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0

    num = round(num, precision)

    if num - int(num) == 0:
        num = int(num)
    num_str = '{}{}'.format(num, ['', 'k', 'M', 'G', 'T', 'P'][magnitude])
    return num_str if not lowercase else num_str.lower()


def str_to_int(num_string, decimal_separator='.', thousand_separator=','):
    try:
        num_string = num_string.replace(thousand_separator, '').lower()
    except AttributeError:
        return num_string

    try:
        return int(num_string)
    except ValueError:
        i = 0
        while i < len(num_string) and (num_string[i].isdigit() or num_string[i] == decimal_separator):
            i += 1

        try:
            number = float(num_string[:i])
            suffix = num_string[i:]

            multipliers = {
                'gb': 1000000000,
                'mb': 1000000,
                'kb': 1000,
                'bp': 1,
                'g': 1000000000,
                'm': 1000000,
                'k': 1000,
                'b': 1,
            }

            return int(number * multipliers[suffix])
        except (KeyError, ValueError):
            raise ValueError("Cannot convert '{}' to integer!".format(num_string))


def natural_sort(l):
    def convert(text):
        return int(text) if text.isdigit() else text.lower()
    return sorted(l, key=lambda key: [convert(c) for c in re.split('([0-9]+)', key)])


def natural_cmp(pa, pb):
    i = 0
    j = 0

    if pa[-2] == '/' == pb[-2]:
        pa = pa[:-2]
        pb = pb[:-2]

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
    return 1 if i < len(pa) else -1 if j < len(pb) else 0


def ligation_site_pattern(restriction_enzyme):
    if isinstance(restriction_enzyme, string_types):
        if "^" in restriction_enzyme and "_" in restriction_enzyme:
            cut_pattern = restriction_enzyme
        else:
            restriction_enzyme = getattr(Restriction, restriction_enzyme)
            cut_pattern = restriction_enzyme.elucidate()
    else:
        cut_pattern = restriction_enzyme.elucidate()

    left_side = []
    right_side = []
    for character in cut_pattern:
        if not (len(left_side) > 0 and left_side[-1] == '_'):
            if not character == '^':
                left_side.append(character)

        if character == '^' or len(right_side) > 0:
            if not character == '_':
                right_side.append(character)

    left_side_re = "".join(left_side[:-1])
    right_side_re = "".join(right_side[1:])

    ll = len(left_side_re)
    lr = len(right_side_re)

    forward = left_side_re + right_side_re
    reverse = reverse_complement(forward)

    forward = forward.replace('N', '[ACGT]')
    reverse = reverse.replace('N', '[ACGT]')

    return forward, ll, reverse, lr


def split_at_ligation_junction(sequence, pattern):
    patterns = []
    if isinstance(pattern, string_types):
        patterns.append(ligation_site_pattern(pattern))
    elif isinstance(pattern, tuple) or isinstance(pattern, list):
        if (isinstance(pattern, tuple) and len(pattern) == 4 and
                isinstance(pattern[1], int) and isinstance(pattern[3], int)):
            patterns.append(pattern)
        else:
            for p in pattern:
                if isinstance(p, tuple) and len(p) == 4:
                    patterns.append(p)
                else:
                    patterns.append(ligation_site_pattern(p))

    if isinstance(sequence, bytes):
        sequence = sequence.decode()

    hits = []
    for pattern in patterns:
        forward, lf, reverse, lr = pattern
        for m in re.finditer(forward, sequence, re.IGNORECASE):
            hits.append(m.start() + lf)

        if forward != reverse:
            for m in re.finditer(reverse, sequence, re.IGNORECASE):
                hits.append(m.start() + lr)

    sub_sequences = []
    previous_hit = 0
    for hit in sorted(hits):
        sub_sequences.append(sequence[previous_hit:hit])
        previous_hit = hit
    sub_sequences.append(sequence[previous_hit:])

    return sub_sequences


def add_dict(x, y):
    return {k: x.get(k, 0) + y.get(k, 0) for k in set(x) | set(y)}


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
        self.next_update = min_value

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


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


class WorkerMonitor(object):
    def __init__(self, value=0, manager=None):
        if manager is None:
            manager = mp.Manager()

        self.counter_lock = manager.Lock()
        self.worker_lock = manager.Lock()

        with self.counter_lock:
            self.val = value

        with self.worker_lock:
            self.worker_states = dict()

    def increment(self):
        with self.counter_lock:
            self.val += 1

    def value(self):
        with self.counter_lock:
            return self.val

    def set_worker_busy(self, worker_uuid):
        with self.worker_lock:
            self.worker_states[worker_uuid] = True

    def set_worker_idle(self, worker_uuid):
        with self.worker_lock:
            self.worker_states[worker_uuid] = False

    def get_worker_state(self, worker_uuid):
        with self.worker_lock:
            return self.worker_states[worker_uuid]

    def workers_idle(self):
        with self.worker_lock:
            for busy in self.worker_states.values():
                if busy:
                    return False
            return True


def get_sam_mapper(sam_file):
    try:
        if isinstance(sam_file, pysam.AlignmentFile):
            return sam_file.header['PG'][0]['ID']
        else:
            with pysam.AlignmentFile(sam_file) as sam:
                return sam.header['PG'][0]['ID']
    except Exception:
        return False
