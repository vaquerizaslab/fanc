cdef extern from  "natural_cmp.c":
    int strnum_cmp(const char *_a, const char *_b)

slash_int = b'/'[0]

def natural_cmp(a, b, ignore_pe_suffix=True):
    assert (isinstance(a, bytes))
    assert (isinstance(b, bytes))

    if ignore_pe_suffix and len(a) > 1 and len(b) > 1 and a[-2] == b[-2] == slash_int:
        return strnum_cmp(a[:-2], b[:-2])
    else:
        return strnum_cmp(a, b)
