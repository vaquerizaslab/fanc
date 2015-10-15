'''
Created on Aug 28, 2015

@author: kkruse1
'''

import itertools
import random

def ranges(i):
    for _, b in itertools.groupby(enumerate(i), lambda (x, y): y - x):
        b = list(b)
        yield b[0][1], b[-1][1]
        

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
    