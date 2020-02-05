import os

# start snippet fancplot import
import fanc
import fanc.plotting as fancplot
# end snippet fancplot import

# start snippet fancplot matplotlib
import matplotlib.pyplot as plt
# end snippet fancplot matplotlib

fanc_base = '..'
#fanc_base = '/Users/kkruse/dev/fanc'

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


# start snippet fancplot square string
hp = fancplot.SquareMatrixPlot(hic)
hp.plot('chr18:6mb-10mb')
hp.show()
# end snippet fancplot square string

hp.save(os.path.join(fanc_base, 'docsrc/api/plot/images/plot_square.png'), figsize=(6, 6))


# start snippet fancplot square vmax
hp = fancplot.SquareMatrixPlot(hic, vmax=0.05)
hp.plot('chr18:6mb-10mb')
hp.show()
# end snippet fancplot square vmax

hp.save(os.path.join(fanc_base, 'docsrc/api/plot/images/plot_square_vmax.png'), figsize=(6, 6))


# start snippet fancplot square cmap
hp = fancplot.SquareMatrixPlot(hic, vmax=0.05, colormap='white_red')
hp.plot('chr18:6mb-10mb')
hp.show()
# end snippet fancplot square cmap

hp.save(os.path.join(fanc_base, 'docsrc/api/plot/images/plot_square_cmap.png'), figsize=(6, 6))


# start snippet fancplot square log
hp = fancplot.SquareMatrixPlot(hic, vmax=0.05, norm='log')
hp.plot('chr18:6mb-10mb')
hp.show()
# end snippet fancplot square log

hp.save(os.path.join(fanc_base, 'docsrc/api/plot/images/plot_square_log.png'), figsize=(6, 6))

# start snippet fancplot square uncorrected
hp = fancplot.SquareMatrixPlot(hic, vmax=30, matrix_norm=False)
hp.plot('chr18:6mb-10mb')
hp.show()
# end snippet fancplot square uncorrected

hp.save(os.path.join(fanc_base, 'docsrc/api/plot/images/plot_square_uncorrected.png'), figsize=(6, 6))


# start snippet fancplot square oe
hp = fancplot.SquareMatrixPlot(hic, vmin=-2, vmax=2, oe=True, log=True,
                               colormap='RdBu_r')
hp.plot('chr18:6mb-10mb')
hp.show()
# end snippet fancplot square oe

hp.save(os.path.join(fanc_base, 'docsrc/api/plot/images/plot_square_oe.png'), figsize=(6, 6))


fig, ax = plt.subplots(figsize=(6, 3))
# start snippet fancplot triangular maxdist
hp = fancplot.TriangularMatrixPlot(hic, vmax=0.05, max_dist='2mb')
hp.plot('chr18:6mb-10mb')
hp.show()
# end snippet fancplot triangular maxdist

hp.save(os.path.join(fanc_base, 'docsrc/api/plot/images/plot_triangular_maxdist.png'), figsize=(6, 3))


fig, ax = plt.subplots(figsize=(10, 10))
# start snippet fancplot mirror
# first create two triangular plots
top_plot = fancplot.TriangularMatrixPlot(hic, vmax=0.05, max_dist='2mb')
bottom_plot = fancplot.TriangularMatrixPlot(hic, vmin=-2, vmax=2, oe=True, log=True, colormap='RdBu_r', max_dist='2mb')
# then merge them
hp = fancplot.MirrorMatrixPlot(top_plot, bottom_plot)
hp.plot('chr18:6mb-10mb')
hp.show()
# end snippet fancplot mirror

hp.save(os.path.join(fanc_base, 'docsrc/api/plot/images/plot_mirror.png'))