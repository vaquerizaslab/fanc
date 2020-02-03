def cis_trans_ratio(hic, normalise=False):
    """
    Calculate the cis/trans ratio for a Hic object.

    :param hic: :class:`~fanc.Hic` object
    :param normalise: If True, will normalise ratio to the possible number of cis/trans contacts
                      in this genome. Makes ratio comparable across different genomes
    :return: tuple (ratio, cis, trans, factor)
    """
    cis = 0
    trans = 0
    regions_dict = hic.regions_dict
    for edge in hic.edges(lazy=True, norm=False):
        if regions_dict[edge.source].chromosome == regions_dict[edge.sink].chromosome:
            cis += edge.weight
        else:
            trans += edge.weight
    if not normalise:
        return cis / (cis + trans), cis, trans, 1.0

    intra_total, chromosome_intra_total, inter_total = hic.possible_contacts()

    f = intra_total / inter_total

    return cis / (cis + trans * f), cis, trans, f

