# start snippet comparisons setup
import fanc
import fanc.plotting as fancplot

import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

hic_esc = fanc.load("architecture/other-hic/bonev2017_esc_10kb.chr1.hic")
hic_cn = fanc.load("architecture/other-hic/bonev2017_cn_10kb.chr1.hic")

ins_esc = fanc.load("architecture/domains/bonev2017_esc_10kb.chr1.insulation")
ins_cn = fanc.load("architecture/domains/bonev2017_cn_10kb.chr1.insulation")

ins_esc_100kb = ins_esc.score_regions(100000)
ins_cn_100kb = ins_cn.score_regions(100000)
# end snippet comparisons setup
import matplotlib.pyplot as plt

# start snippet comparisons compare
diff_esc_cn = fanc.DifferenceMatrix.from_matrices(hic_esc, hic_cn,
                                                  file_name="architecture/comparisons/esc_vs_cn.diff")
fc_esc_cn = fanc.FoldChangeMatrix.from_matrices(hic_esc, hic_cn, log_cmp=True,
                                                file_name="architecture/comparisons/esc_vs_cn.fc")
# end snippet comparisons compare

# start snippet comparisons plot
p_esc = fancplot.TriangularMatrixPlot(hic_esc, vmin=0, vmax=0.05,
                                      max_dist=400000)
p_cn = fancplot.TriangularMatrixPlot(hic_cn, vmin=0, vmax=0.05,
                                     max_dist=400000)
p_diff = fancplot.TriangularMatrixPlot(diff_esc_cn, vmin=-0.02, vmax=0.02,
                                       colormap='bwr', max_dist=400000)
p_fc = fancplot.TriangularMatrixPlot(fc_esc_cn, vmin=-1.5, vmax=1.5,
                                     colormap='PiYG', max_dist=400000)

gf = fancplot.GenomicFigure([p_esc, p_cn, p_diff, p_fc], ticks_last=True)
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


# start snippet regions compare
diff_ins_esc_cn_100kb = fanc.DifferenceRegions.from_regions(ins_esc_100kb, ins_cn_100kb)
# end snippet regions compare

# start snippet regions plot
p_orig = fancplot.LinePlot([ins_esc_100kb, ins_cn_100kb], ylim=(-1, 1),
                           colors=['darkturquoise', 'orange'],
                           style='mid', fill=False)
p_diff = fancplot.LinePlot(diff_ins_esc_cn_100kb, ylim=(-1., 1.),
                           colors=['aquamarine'],
                           style='mid')
gf = fancplot.GenomicFigure([p_orig, p_diff], ticks_last=True)
fig, axes = gf.plot("chr1:167.9mb-168.7mb")

axes[0].set_ylabel("Insulation\nscore")
axes[1].set_ylabel("Insulation\ndifference")
# end snippet regions plot
fig.savefig("../docsrc/api/analyse/images/comparisons_regions.png")
plt.close(fig)


# start snippet scores compare
diff_ins_esc_cn = fanc.DifferenceScores.from_scores(ins_esc, ins_cn)
# end snippet scores compare

# start snippet scores plot
p_esc = fancplot.GenomicVectorArrayPlot(ins_esc, colormap='RdBu_r',
                                        vmin=-0.5, vmax=0.5,
                                        genomic_format=True, y_scale='log')
p_cn = fancplot.GenomicVectorArrayPlot(ins_cn, colormap='RdBu_r',
                                       vmin=-0.5, vmax=0.5,
                                       genomic_format=True, y_scale='log')
p_diff = fancplot.GenomicVectorArrayPlot(diff_ins_esc_cn, colormap='PiYG',
                                         vmin=-0.5, vmax=0.5,
                                         genomic_format=True, y_scale='log')

gf = fancplot.GenomicFigure([p_esc, p_cn, p_diff], ticks_last=True)
fig, axes = gf.plot("chr1:167.9mb-168.7mb")

axes[0].set_ylabel("Insulation\nscore ESC")
axes[1].set_ylabel("Insulation\nscore CN")
axes[2].set_ylabel("Insulation\ndifference")
# end snippet scores plot
fig.savefig("../docsrc/api/analyse/images/comparisons_scores.png")
plt.close(fig)