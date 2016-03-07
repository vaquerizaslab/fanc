'''
Created on Aug 28, 2015

@author: kkruse1
'''

import itertools
import random
import collections


class CachedIterator:
    def __init__(self, it, cache_size=1):
        self.it = iter(it)
        self.current = None
        self.cache = collections.deque(maxlen=cache_size)
        self.cache_offset = 0

    def __iter__(self):
        return self

    def next(self):
        try:
            if self.cache_offset:
                self.cache_offset -= 1
                if self.cache_offset:
                    return self.cache[-self.cache_offset]
                return self.current
            self.cache.append(self.current)
            self.current = self.it.next()
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
    for _, b in itertools.groupby(enumerate(i), lambda (x, y): y - x):
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
            number = number - power
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
