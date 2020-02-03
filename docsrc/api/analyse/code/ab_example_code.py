# start snippet ab setup
import fanc
import fanc.plotting as fancplot
import matplotlib.pyplot as plt

hic_1mb = fanc.load("output/hic/binned/fanc_example_1mb.hic")
# end snippet ab setup


# start snippet ab matrix
ab = fanc.ABCompartmentMatrix.from_hic(hic_1mb)
# end snippet ab matrix

# start snippet ab subset
ab_chr18 = ab.matrix(('chr18', 'chr18'))
# end snippet ab subset

# start snippet ab fancplot-correlation
fig, ax = plt.subplots()
mp = fancplot.SquareMatrixPlot(ab, ax=ax,
                           norm='lin', colormap='RdBu_r',
                           vmin=-1, vmax=1,
                           draw_minor_ticks=False)
mp.plot('chr18')
plt.show()
# end snippet ab fancplot-correlation
fig.savefig('../fanc-doc/api/analyse/images/ab_1mb_correlation.png')


# start snippet ab ev
ev = ab.eigenvector()
# end snippet ab ev

# start snippet ab gc-ev
gc_ev = ab.eigenvector(genome='hg19_chr18_19.fa', force=True)
# end snippet ab gc-ev


# start snippet ab plot-ev
fig, ax = plt.subplots()
lp = fancplot.LinePlot(ab)
lp.plot('chr18')
plt.show()
# end snippet ab plot-ev
fig.savefig('../fanc-doc/api/analyse/images/ab_1mb_ev.png')

