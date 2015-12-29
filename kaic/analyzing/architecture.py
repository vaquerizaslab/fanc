from kaic.data.genomic import Hic
import numpy as np
import logging
log = logging.getLogger(__name__)
log.setLevel(10)

def insulation_index(hic, d):
    chr_bins = hic.chromosome_bins
    n = len(hic.regions())
    log.debug("Fetching matrix")
    hic_matrix = hic[:,:]
    ins_matrix = np.empty(n)
    log.debug("Starting processing")
    for i, r in enumerate(hic.regions()):
        if i - chr_bins[r.chromosome][0] < d:
            ins_matrix[i] = np.nan
            continue
        if chr_bins[r.chromosome][1] - i < d:
            ins_matrix[i] = np.nan
            continue
        ins_matrix[i] = np.ma.sum(hic_matrix[i: i + d, i - d:i])
    return ins_matrix

def contact_band(hic, d1, d2, use_oe_ratio=False):
    chr_bins = hic.chromosome_bins
    n = len(hic.regions())
    log.debug("Fetching matrix")
    hic_matrix = hic[:,:]
    band = np.empty(n)
    if use_oe_ratio:
        raise NotImplementedError("oe not implemented yet :(")
    log.debug("Starting processing")
    for i, r in enumerate(hic.regions()):
        if i - chr_bins[r.chromosome][0] < d2:
            band[i] = np.nan
            continue
        if chr_bins[r.chromosome][1] - i < d2:
            band[i] = np.nan
            continue
        band[i] = np.ma.sum(hic_matrix[i - d2:i - d1, i + d1:i + d2])
    return band

def observed_expected_ratio(hic, per_chromosome=True):
    pass
