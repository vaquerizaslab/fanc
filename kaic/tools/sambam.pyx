cdef extern from  "natural_cmp.c":
    int strnum_cmp(const char *_a, const char *_b)

slash_int = b'/'[0]

def natural_cmp(a, b):
    assert (isinstance(a, bytes))
    assert (isinstance(b, bytes))

    if a[-2] == b[-2] == slash_int:
        return strnum_cmp(a[:-2], b[:-2])
    else:
        return strnum_cmp(a, b)
