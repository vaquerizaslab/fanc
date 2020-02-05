import os
import fanc
import fanc.plotting as fancplot

# start snippet fancplot matplotlib
import matplotlib.pyplot as plt
# end snippet fancplot matplotlib

#fanc_base = '..'
fanc_base = '/Users/kkruse/dev/fanc'

# start snippet fancplot load bed
insulation_scores_1mb = fanc.load("architecture/domains/fanc_example_100kb.insulation_1mb.bed")
insulation_scores_2mb = fanc.load("architecture/domains/fanc_example_100kb.insulation_2mb.bed")
boundaries_1mb = fanc.load("architecture/domains/fanc_example_100kb.insulation_boundaries_1mb.bed")
boundaries_2mb = fanc.load("architecture/domains/fanc_example_100kb.insulation_boundaries_2mb.bed")
# end snippet fancplot load bed

# start snippet fancplot line fill
hp = fancplot.LinePlot(insulation_scores_1mb)
hp.plot('chr18:6mb-10mb')
hp.show()
# end snippet fancplot line fill

hp.save(os.path.join(fanc_base, 'docsrc/api/plot/images/plot_line.png'))

# start snippet fancplot line nofill
hp = fancplot.LinePlot(insulation_scores_1mb, fill=False)
hp.plot('chr18:6mb-10mb')
hp.show()
# end snippet fancplot line nofill

hp.save(os.path.join(fanc_base, 'docsrc/api/plot/images/plot_line_nofill.png'))

# start snippet fancplot line mid
hp = fancplot.LinePlot(insulation_scores_1mb, style='mid')
hp.plot('chr18:6mb-10mb')
hp.show()
# end snippet fancplot line mid

hp.save(os.path.join(fanc_base, 'docsrc/api/plot/images/plot_line_mid.png'))

# start snippet fancplot line col
hp = fancplot.LinePlot(insulation_scores_1mb, style='mid', colors='cyan',
                       plot_kwargs={'alpha': 0.5})
hp.plot('chr18:6mb-10mb')
hp.show()
# end snippet fancplot line col

hp.save(os.path.join(fanc_base, 'docsrc/api/plot/images/plot_line_col.png'))

# start snippet fancplot line ylim
hp = fancplot.LinePlot(insulation_scores_1mb, style='mid', colors='cyan',
                       plot_kwargs={'alpha': 0.5}, ylim=(-1, 1))
hp.plot('chr18:6mb-10mb')
hp.show()
# end snippet fancplot line ylim

hp.save(os.path.join(fanc_base, 'docsrc/api/plot/images/plot_line_ylim.png'))

# start snippet fancplot line multi
hp = fancplot.LinePlot((insulation_scores_1mb, insulation_scores_2mb),
                       style='mid', colors=('cyan', 'magenta'),
                       plot_kwargs={'alpha': 0.5}, ylim=(-1, 1),
                       labels=('1mb', '2mb'))
hp.plot('chr18:6mb-10mb')
hp.show()
# end snippet fancplot line multi

hp.save(os.path.join(fanc_base, 'docsrc/api/plot/images/plot_line_multi.png'))

# start snippet fancplot bar basic
hp = fancplot.BarPlot((insulation_scores_1mb, insulation_scores_2mb),
                      colors=('cyan', 'magenta'),
                      plot_kwargs={'alpha': 0.5}, ylim=(-1, 1),
                      labels=('1mb', '2mb'))
hp.plot('chr18:6mb-10mb')
hp.show()
# end snippet fancplot bar basic

hp.save(os.path.join(fanc_base, 'docsrc/api/plot/images/plot_bar.png'))

# start snippet fancplot bar boundaries
hp = fancplot.BarPlot(boundaries_1mb, colors='cyan')
hp.plot('chr18:20mb-30mb')
hp.show()
# end snippet fancplot bar boundaries

hp.save(os.path.join(fanc_base, 'docsrc/api/plot/images/plot_bar_boundaries.png'))