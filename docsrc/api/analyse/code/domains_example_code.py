# start snippet domains setup
import fanc
import fanc.plotting as fancplot
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

hic_100kb = fanc.load("output/hic/binned/fanc_example_100kb.hic")
# end snippet domains setup
import matplotlib.pyplot as plt

# start snippet alternative cooler
hic_100kb = fanc.load("architecture/other-hic/fanc_example.mcool@100kb")
# end snippet alternative cooler

# start snippet alternative juicer
hic_100kb = fanc.load("architecture/other-hic/fanc_example.juicer.hic@100kb")
# end snippet alternative juicer

# start snippet insulation basic
insulation = fanc.InsulationScores.from_hic(hic_100kb,
                                            [5000000, 1000000, 1500000,
                                             2000000, 2500000, 3000000,
                                             3500000, 4000000],
                                            file_name="architecture/domains/fanc_example_100kb.insulation")
# end snippet insulation basic

# start snippet insulation scores
scores = insulation.scores(1000000)
# end snippet insulation scores

# start snippet insulation regions
for region in insulation.regions('chr18'):
    score = region.insulation_1000000
    # do something with score
    # ...
# end snippet insulation regions

# start snippet insulation score_regions
insulation_1mb = insulation.score_regions(1000000)
for r in insulation_1mb.regions:
    score = r.score
# end snippet insulation score_regions

# start snippet insulation export
insulation.to_bed("architecture/domains/fanc_example_100kb.insulation_1mb.bed", 1000000)
# end snippet insulation export

fig, ax = plt.subplots(figsize=(6, 3))
# start snippet insulation multiplot
p = fancplot.GenomicVectorArrayPlot(insulation, colormap='RdBu_r', vmin=-1, vmax=1,
                                    genomic_format=True)
p.plot('chr18:18mb-28mb')
# end snippet insulation multiplot
fig.savefig("../docsrc/api/analyse/images/domains_multi.png")
plt.close(fig)

fig, ax = plt.subplots()
# start snippet insulation singleplot
p = fancplot.LinePlot(insulation, ylim=(-1, 1), colors=['darkturquoise'], style="mid",
                      attribute="insulation_1000000")
p.plot('chr18:18mb-28mb')
# end snippet insulation singleplot
fig.savefig("../docsrc/api/analyse/images/domains_single.png")
plt.close(fig)


# start snippet boundaries run
boundaries = fanc.Boundaries.from_insulation_score(insulation, window_size=1000000)
# end snippet boundaries run

# start snippet boundaries regions
for boundary in boundaries.regions:
    score = boundary.score
# end snippet boundaries regions

# start snippet boundaries plot
ph = fancplot.TriangularMatrixPlot(hic_100kb, max_dist=5000000, vmin=0, vmax=0.05)
pb = fancplot.BarPlot(boundaries)
f = fancplot.GenomicFigure([ph, pb])
fig, axes = f.plot('chr18:18mb-28mb')
# end snippet boundaries plot
fig.savefig("../docsrc/api/analyse/images/boundaries.png")
plt.close(fig)


# start snippet boundaries fp
min_score = 1.0
robust_boundaries = []
for boundary in boundaries.regions:
    score = boundary.score
    if score >= min_score:
        robust_boundaries.append(boundary)
# end snippet boundaries fp

# start snippet boundaries minscore
robust_boundaries = fanc.Boundaries.from_insulation_score(insulation, window_size=1000000,
                                                          min_score=1.0)
# end snippet boundaries minscore


# start snippet directionality basic
directionality = fanc.DirectionalityIndexes.from_hic(hic_100kb,
                                                     [1000000, 1500000,
                                                      2000000, 2500000, 3000000,
                                                      3500000, 4000000],
                                                     file_name="architecture/domains/fanc_example_100kb.di")
# end snippet directionality basic

fig, ax = plt.subplots(figsize=(6, 3))
# start snippet directionality multiplot
p = fancplot.GenomicVectorArrayPlot(directionality, colormap='RdBu_r', vmin=-0.1, vmax=0.1,
                                    genomic_format=True)
p.plot('chr18:18mb-28mb')
# end snippet directionality multiplot
fig.savefig("../docsrc/api/analyse/images/directionality_multi.png")
plt.close(fig)

fig, ax = plt.subplots()
# start snippet directionality singleplot
p = fancplot.LinePlot(directionality, ylim=(-0.1, 0.1), colors=['darkturquoise'], style="mid",
                      attribute="directionality_4000000")
p.plot('chr18:18mb-28mb')
# end snippet directionality singleplot
fig.savefig("../docsrc/api/analyse/images/directionality_single.png")
plt.close(fig)