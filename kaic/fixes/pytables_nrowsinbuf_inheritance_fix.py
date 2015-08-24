'''
Monkey-patch _calc_nrowsinbuf method to account for
classes inheriting Table

@author: kaukrise
'''

from tables import Leaf
from tables.exceptions import PerformanceWarning
import warnings

def monkeypatch_method(cls):
    def decorator(func):
        setattr(cls, func.__name__, func)
        return func
    return decorator


@monkeypatch_method(Leaf)
def _calc_nrowsinbuf(self):
    """Calculate the number of rows that fits on a PyTables buffer."""

    params = self._v_file.params
    # Compute the nrowsinbuf
    rowsize = self.rowsize
    buffersize = params['IO_BUFFER_SIZE']
    if rowsize != 0:
        nrowsinbuf = buffersize // rowsize
        if (self.__class__.__name__ == "Table" or
            "Table" in [base.__name__ for base in self.__class__.__bases__]):
            # The number of rows in buffer needs to be an exact multiple of
            # chunkshape[0] for queries using indexed columns.
            # Fixes #319 and probably #409 too.
            nrowsinbuf -= nrowsinbuf % self.chunkshape[0]
    else:
        nrowsinbuf = 1

    # tableextension.pyx performs an assertion
    # to make sure nrowsinbuf is greater than or
    # equal to the chunksize.
    # See gh-206 and gh-238
    if (self.chunkshape is not None and (self.__class__.__name__ == "Table" or
        "Table" in [base.__name__ for base in self.__class__.__bases__])):
        if nrowsinbuf < self.chunkshape[0]:
            nrowsinbuf = self.chunkshape[0]

    # Safeguard against row sizes being extremely large
    if nrowsinbuf == 0:
        nrowsinbuf = 1
        # If rowsize is too large, issue a Performance warning
        maxrowsize = params['BUFFER_TIMES'] * buffersize
        if rowsize > maxrowsize:
            warnings.warn("""\
The Leaf ``%s`` is exceeding the maximum recommended rowsize (%d bytes);
be ready to see PyTables asking for *lots* of memory and possibly slow
I/O.  You may want to reduce the rowsize by trimming the value of
dimensions that are orthogonal (and preferably close) to the *main*
dimension of this leave.  Alternatively, in case you have specified a
very small/large chunksize, you may want to increase/decrease it."""
                          % (self._v_pathname, maxrowsize),
                          PerformanceWarning)
    return nrowsinbuf