cdef extern from  "natural_cmp.c":
    int strnum_cmp(const char *_a, const char *_b)

def natural_cmp(a, b):
    assert (isinstance(a, bytes))
    assert (isinstance(b, bytes))
    return strnum_cmp(a, b)
