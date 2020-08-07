# start snippet pca setup
import fanc
import fanc.plotting as fancplot
from fanc.architecture.comparisons import hic_pca

lowc_hindiii_100k = fanc.load("architecture/other-hic/lowc_hindiii_100k_1mb.hic")
lowc_hindiii_5M = fanc.load("architecture/other-hic/lowc_hindiii_5M_1mb.hic")
lowc_mboi_1M = fanc.load("architecture/other-hic/lowc_mboi_1M_1mb.hic")
lowc_mboi_100k = fanc.load("architecture/other-hic/lowc_mboi_100k_1mb.hic")
lowc_mboi_50k = fanc.load("architecture/other-hic/lowc_mboi_50k_1mb.hic")
# end snippet pca setup
import matplotlib.pyplot as plt


# start snippet pca run
pca_info, pca_result = hic_pca(lowc_hindiii_5M, lowc_hindiii_100k,
                               lowc_mboi_1M, lowc_mboi_50k,
                               lowc_mboi_100k, ignore_zeros=True,
                               scale=False, region='chr19', sample_size=100000,
                               strategy='variance')
# end snippet pca run

# start snippet pca plot
fig, ax = fancplot.pca_plot(pca_result, variance=pca_info,
                            names=["HindIII 5M", "HindIII 100k",
                                   "MboI 1M", "MboI 50k", "MboI 100k"])
# end snippet pca plot
fig.savefig('../docsrc/api/analyse/images/pca_default.png')
plt.close(fig)


# start snippet pca adjust
fig, ax = fancplot.pca_plot(pca_result, variance=pca_info,
                            names=["HindIII 5M", "HindIII 100k",
                                   "MboI 1M", "MboI 50k", "MboI 100k"],
                            colors=["orange", "orange", "cyan", "cyan", "cyan"],
                            markers=["*", "o", "*", "o", "o"])
# end snippet pca adjust
fig.savefig('../docsrc/api/analyse/images/pca_adjust.png')
plt.close(fig)


# start snippet pca ev
fig, ax = fancplot.pca_plot(pca_result, variance=pca_info,
                            eigenvectors=(1, 2),
                            names=["HindIII 5M", "HindIII 100k",
                                   "MboI 1M", "MboI 50k", "MboI 100k"],
                            colors=["orange", "orange", "cyan", "cyan", "cyan"],
                            markers=["*", "o", "*", "o", "o"])
# end snippet pca ev
fig.savefig('../docsrc/api/analyse/images/pca_ev.png')
plt.close(fig)
