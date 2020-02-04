import os

# start snippet fancplot import
import fanc
import fanc.plotting as fancplot
# end snippet fancplot import

# start snippet fancplot matplotlib
import matplotlib.pyplot as plt
# end snippet fancplot matplotlib

#fanc_base = '..'
fanc_base = '/Users/kkruse/dev/fanc'

# start snippet fancplot load hic
hic = fanc.load("output/hic/binned/fanc_example_50kb.hic")
# end snippet fancplot load hic

# start snippet fancplot triangular string
hp = fancplot.TriangularMatrixPlot(hic, vmax=0.05)
hp.plot('chr18:6mb-10mb')
hp.show()
# end snippet fancplot triangular string

hp.save(os.path.join(fanc_base, 'docsrc/api/plot/images/plot_triangular.png'), figsize=(6, 3))


# start snippet fancplot plt axes
fig, ax = plt.subplots()
hp = fancplot.TriangularMatrixPlot(hic, vmax=0.05, ax=ax)
hp.plot('chr18:6mb-10mb')
ax.set_xticks([7500000, 8500000])
ax.set_xticklabels(['customise', 'everything!'])
hp.show()
# end snippet fancplot plt axes

hp.save(os.path.join(fanc_base, 'docsrc/api/plot/images/plot_triangular_plt.png'), figsize=(6, 3))