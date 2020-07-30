# start snippet comparisons setup
import fanc
import fanc.plotting as fancplot

import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

hic_esc = fanc.load("output/hic/binned/fanc_example_100kb.hic")
hic_cn = fanc.load("output/hic/binned/fanc_example_100kb.hic")

ins_esc_bed = fanc.load("output/hic/binned/fanc_example_100kb.hic")
ins_cn_bed = fanc.load("output/hic/binned/fanc_example_100kb.hic")


# end snippet comparisons setup
import matplotlib.pyplot as plt

# start snippet comparisons compare
diff_esc_cn = fanc.DifferenceMatrix.from_matrices(hic_esc, hic_cn,
                                                  file_name="architecture/comparisons/esc_vs_cn.diff")
fc_esc_cn = fanc.FoldChangeMatrix.from_matrices(hic_esc, hic_cn, log_cmp=True,
                                                file_name="architecture/comparisons/esc_vs_cn.fc")
# end snippet comparisons compare

# start snippet comparisons plot
p_esc = fancplot.TriangularMatrixPlot(hic_esc, vmin=0, vmax=0.05)
p_cn = fancplot.TriangularMatrixPlot(hic_cn, vmin=0, vmax=0.05)
p_diff = fancplot.TriangularMatrixPlot(diff_esc_cn, vmin=-0.03, vmax=0.03)
p_fc = fancplot.TriangularMatrixPlot(fc_esc_cn, vmin=-1, vmax=1)

gf = fancplot.GenomicFigure([p_esc, p_cn, p_diff, p_fc])
fig, axes = gf.plot("chr1:167.9mb-168.7mb")
# end snippet comparisons plot
fig.savefig("../docsrc/api/analyse/images/comparisons_matrices.png")
plt.close(fig)

# start snippet comparisons custom
from fanc.architecture.comparisons import ComparisonMatrix


class CustomComparisonMatrix(ComparisonMatrix):
    _classid = 'CUSTOMCOMPARISONMATRIX'

    def __init__(self, *args, **kwargs):
        ComparisonMatrix.__init__(self, *args, **kwargs)

    def compare(self, weight1, weight2):
        if weight1 < weight2:
            return -1
        if weight1 > weight2:
            return 1
        return 0
# end snippet comparisons custom


# start snippet scores compare
diff_ins_esc_cn = fanc.DifferenceRegions.from_regions(ins_esc_bed, ins_cn_bed)
fc_ins_esc_cn = fanc.FoldChangeRegions.from_regions(ins_esc_bed, ins_cn_bed)
# end snippet scores compare

# start snippet comparisons plot
p_orig = fancplot.LinePlot([ins_esc_bed, ins_cn_bed], vmin=-1, vmax=1,
                           colors=['darkturquoise', 'orange'],
                           style='mid')
p_diff = fancplot.LinePlot(diff_ins_esc_cn, vmin=-0.7, vmax=0.7,
                           colors=['green'],
                           style='mid')
p_fc = fancplot.LinePlot(fc_ins_esc_cn, vmin=-1, vmax=1,
                         colors=['grey'],
                         style='mid')
gf = fancplot.GenomicFigure([p_orig, p_diff, p_fc])
fig, axes = gf.plot("chr1:167.9mb-168.7mb")
# end snippet comparisons plot
fig.savefig("../docsrc/api/analyse/images/comparisons_regions.png")
plt.close(fig)