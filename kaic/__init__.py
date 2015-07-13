"""

KaiC
====

Provides
    1. Classes for working with Hi-C data
    2. Classes for working with tabular data

"""

import logging
logging.basicConfig(level=logging.INFO)

from .data.genomic import HicBasic, HicNode, HicEdge, Genome, Chromosome
from .data.general import Table 

