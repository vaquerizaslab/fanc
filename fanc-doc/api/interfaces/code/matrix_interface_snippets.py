import fanc
import numpy as np

# start snippet check
hic = fanc.load("examples/output/hic/binned/fanc_example_500kb.hic")
isinstance(hic, fanc.matrix.RegionMatrixContainer)  # True if interface supported
# end snippet check

# start snippet matrix whole genome
m = hic.matrix()
# end snippet matrix whole genome

# start snippet matrix type
type(m)  # fanc.matrix.RegionMatrix
isinstance(m, np.ndarray)  # True
# end snippet matrix type

# start snippet matrix no mask
m_unmasked = hic.matrix(mask=False)
# end snippet matrix no mask

# start snippet matrix subset
m_chr19 = hic.matrix(('chr19', 'chr19'))
# end snippet matrix subset

# start snippet matrix row col
m_inter1 = hic.matrix(('chr18', 'chr19'))
m_inter1.shape  # (157, 119)
m_inter2 = hic.matrix(('chr19', 'chr18'))
m_inter2.shape  # (119, 157)
# end snippet matrix row col

# start snippet matrix marginals
m_chr19.shape  # (119, 119)
marginals = np.sum(m_chr19, axis=0)
marginals.shape  # (119,)
marginals[:5]
# [1.0000000074007467, 0.9999999585562779,
# 1.0000000102533806, 0.999999987196381, 1.0000000140165086]
# end snippet matrix marginals

# start snippet matrix regions
m_inter1.row_regions
# [chr18:1-500000,
#  chr18:500001-1000000,
#  chr18:1000001-1500000,
#  chr18:1500001-2000000,
# ...
m_inter1.col_regions
# [chr19:1-500000,
#  chr19:500001-1000000,
#  chr19:1000001-1500000,
#  chr19:1500001-2000000,
# ...
# end snippet matrix regions

# start snippet matrix region matrix subset
# subset by index
m_chr19_sub1 = m_chr19[0:3, 0:3]
m_chr19_sub1.row_regions
# [chr19:1-500000, chr19:500001-1000000, chr19:1000001-1500000]
m_chr19_sub1.col_regions
# [chr19:1-500000, chr19:500001-1000000, chr19:1000001-1500000]

# subset by region interval
m_chr19_sub2 = m_chr19['chr19:2mb-3mb', 'chr19:500kb-1mb']
m_chr19_sub2.row_regions
# [chr19:1500001-2000000, chr19:2000001-2500000, chr19:2500001-3000000]
m_chr19_sub2.col_regions
# end snippet matrix region matrix subset
# [chr19:1-500000, chr19:500001-1000000]


# start snippet matrix no norm
m_chr19_uncorr = hic.matrix(('chr19', 'chr19'), norm=False)
# end snippet matrix no norm

# start snippet matrix oe
m_chr19_oe = hic.matrix(('chr19', 'chr19'), oe=True)
# end snippet matrix oe

# start snippet matrix log oe
m_chr19_log_oe = hic.matrix(('chr19', 'chr19'), oe=True, log=True)
# end snippet matrix log oe


