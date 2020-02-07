import genomic_regions as gr
import numpy as np
import logging

logger = logging.getLogger(__name__)


def load_contacts(contacts_file_name, sep=None, ix_converter=None):
    m = None
    try:  # numpy binary format
        m = np.load(contacts_file_name)
    except (IOError, ValueError):  # not an .npy file

        # check number of fields in file
        with open(contacts_file_name, 'r') as f:
            line = f.readline()
            while line.startswith('#'):
                line = f.readline()
            line = line.rstrip()
            n_fields = len(line.split(sep))

        if n_fields > 3:  # square matrix format
            m = np.loadtxt(contacts_file_name)

        if m is not None:
            for i, row in enumerate(m):
                for j, val in enumerate(row):
                    yield i, j, val
        else:
            with open(contacts_file_name, 'r') as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    line = line.rstrip()
                    fields = line.split(sep)
                    if ix_converter is None:
                        source, sink, weight = int(fields[0]), int(fields[1]), float(fields[2])
                    else:
                        source = ix_converter[fields[0]]
                        sink = ix_converter[fields[1]]
                        weight = float(fields[2])

                    yield source, sink, weight


def load_regions(file_name, sep=None):
    regions = []
    ix_converter = None
    with open(file_name, 'r') as f:
        for i, line in enumerate(f):
            line = line.rstrip()
            fields = line.split(sep)
            if len(fields) > 2:
                chromosome = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                ix = i
                if len(fields) > 3 and fields[3] != '.':  # HicPro
                    if ix_converter is None:
                        ix_converter = dict()
                    if fields[3] in ix_converter:
                        raise ValueError("Name column in region BED must "
                                         "only contain unique values! ({})".format(fields[3]))
                    ix_converter[fields[3]] = ix
                regions.append(gr.GenomicRegion(chromosome=chromosome, start=start, end=end, ix=ix))
    return regions, ix_converter
