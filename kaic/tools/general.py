'''
Created on Aug 28, 2015

@author: kkruse1
'''

import itertools

def ranges(i):
    for _, b in itertools.groupby(enumerate(i), lambda (x, y): y - x):
        b = list(b)
        yield b[0][1], b[-1][1]