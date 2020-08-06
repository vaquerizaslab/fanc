# start snippet aggregate setup
import fanc
import fanc.plotting as fancplot

import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

hic_100kb = fanc.load("output/hic/binned/fanc_example_100kb.hic")
hic_rao = fanc.load("architecture/loops/rao2014.chr11_77400000_78600000.hic")

boundaries = fanc.load("architecture/domains/fanc_example_100kb.insulation_boundaries_score_0.7_1mb.bed")
tads = fanc.load("architecture/domains/gm12878_tads.bed")
loops = fanc.load("architecture/loops/rao2014.chr11_77400000_78600000.loops_no_singlets.bedpe")
# end snippet aggregate setup
import matplotlib.pyplot as plt

# start snippet aggregate cooler
hic_100kb = fanc.load("architecture/other-hic/fanc_example.mcool@100kb")
# end snippet aggregate cooler

# start snippet aggregate juicer
hic_100kb = fanc.load("architecture/other-hic/fanc_example.juicer.hic@100kb")
# end snippet aggregate juicer

# start snippet aggregate fixed
fixed_am = fanc.AggregateMatrix.from_center(hic_100kb, boundaries.regions,
                                            window=5000000)
# end snippet aggregate fixed

fig, ax = plt.subplots()
# start snippet aggregate plotfixed
ax = fancplot.aggregate_plot(fixed_am, vmin=-1, vmax=1,
                             labels=['-2.5Mb', 'boundary', '+2.5Mb'])
# end snippet aggregate plotfixed
fig.savefig("../docsrc/api/analyse/images/aggregate_fixed.png")
plt.close(fig)


# start snippet aggregate variable
variable_am = fanc.AggregateMatrix.from_regions(hic_100kb, tads.regions)
# end snippet aggregate variable

fig, ax = plt.subplots()
# start snippet aggregate plotvariable
ax = fancplot.aggregate_plot(variable_am, vmin=-0.5, vmax=0.5,
                             relative_label_locations=[0.33, 0.66])
# end snippet aggregate plotvariable
fig.savefig("../docsrc/api/analyse/images/aggregate_variable.png")
plt.close(fig)


# start snippet aggregate rescale
rescale_am = fanc.AggregateMatrix.from_regions(hic_100kb, tads.regions,
                                               rescale=True, log=False)
# end snippet aggregate rescale

fig, ax = plt.subplots()
# start snippet aggregate plotrescale
ax = fancplot.aggregate_plot(rescale_am, vmax=0.045,
                             relative_label_locations=[0.33, 0.66],
                             colormap='germany')
# end snippet aggregate plotrescale
fig.savefig("../docsrc/api/analyse/images/aggregate_rescale.png")
plt.close(fig)


# start snippet aggregate loops
loops_am = fanc.AggregateMatrix.from_center_pairs(hic_rao, loops.region_pairs())
# end snippet aggregate loops

fig, ax = plt.subplots()
# start snippet aggregate plotloops
ax = fancplot.aggregate_plot(loops_am, vmin=-0.5, vmax=0.5,
                             relative_label_locations=[0.5],
                             labels=['loop anchor'])
# end snippet aggregate plotloops
fig.savefig("../docsrc/api/analyse/images/aggregate_loops.png")
plt.close(fig)

# start snippet aggregate matrix
m = fixed_am.matrix()
type(m)  # numpy.ndarray
# end snippet aggregate matrix

# start snippet aggregate components
for component_matrix in fixed_am.components():
    shape = component_matrix.shape  # 51, 51
    # ...
# end snippet aggregate components

