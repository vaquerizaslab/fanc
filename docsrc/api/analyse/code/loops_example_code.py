# start snippet loops setup
import fanc
import fanc.peaks
import fanc.plotting as fancplot

import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

#hic_10kb = fanc.load("architecture/other-hic/bonev2017_esc_10kb.chr1.hic")
hic_sample = fanc.load("architecture/loops/rao2014.chr11_77400000_78600000.hic")
# end snippet loops setup
import matplotlib.pyplot as plt


# start snippet loops annotate
loop_caller = fanc.RaoPeakCaller()
loop_info = loop_caller.call_peaks(hic_sample,
                                   file_name="architecture/loops/rao2014.chr11_77400000_78600000.loop_info")
# end snippet loops annotate


# start snippet loops iter
for edge in loop_info.peaks(lazy=True):
    print(edge.oe_d, edge.fdr_d)
# end snippet loops iter


# start snippet loops matrix
# regular Hi-C
p_hic = fancplot.TriangularMatrixPlot(loop_info, vmin=0, vmax=0.05,
                                      max_dist=600000)
# Donut O/E values
p_oe_d = fancplot.TriangularMatrixPlot(loop_info, vmin=0, vmax=4,
                                       weight_field='oe_d', colormap='Reds',
                                       max_dist=600000)
# Donut FDR values
p_fdr_d = fancplot.TriangularMatrixPlot(loop_info, vmin=0, vmax=0.1,
                                        weight_field='fdr_d', colormap='binary_r',
                                        max_dist=600000)

gf = fancplot.GenomicFigure([p_hic, p_oe_d, p_fdr_d], ticks_last=True)
fig, axes = gf.plot('chr11:77.4mb-78.6mb')
# end snippet loops matrix
fig.savefig("../docsrc/api/analyse/images/loops_annotate.png")
plt.close(fig)


# start snippet loops filters
# filter on O/E
enrichment_filter = fanc.peaks.EnrichmentPeakFilter(
        enrichment_ll_cutoff=1.75,
        enrichment_d_cutoff=3,
        enrichment_h_cutoff=2,
        enrichment_v_cutoff=2)
loop_info.filter(enrichment_filter, queue=True)

# filter on FDR
fdr_filter = fanc.peaks.FdrPeakFilter(
    fdr_ll_cutoff=0.05,
    fdr_d_cutoff=0.05,
    fdr_h_cutoff=0.05,
    fdr_v_cutoff=0.05)
loop_info.filter(fdr_filter, queue=True)

# filter on mappability
mappability_filter = fanc.peaks.MappabilityPeakFilter(
    mappability_ll_cutoff=0.7,
    mappability_d_cutoff=0.7,
    mappability_h_cutoff=0.7,
    mappability_v_cutoff=0.7)
loop_info.filter(mappability_filter, queue=True)

# filter on distance from diagonal
distance_filter = fanc.peaks.DistancePeakFilter(5)
loop_info.filter(distance_filter, queue=True)

# filter on minimum number of observed values
observed_filter = fanc.peaks.ObservedPeakFilter(10)
loop_info.filter(observed_filter, queue=True)

loop_info.run_queued_filters()
# end snippet loops filters


# start snippet loops remaining
# regular Hi-C
p_hic = fancplot.TriangularMatrixPlot(loop_info, vmin=0, vmax=0.05,
                                      max_dist=600000)
# Donut O/E values
p_oe_d = fancplot.TriangularMatrixPlot(loop_info, vmin=0, vmax=4,
                                       weight_field='oe_d', colormap='Reds',
                                       max_dist=600000)
# Donut FDR values
p_fdr_d = fancplot.TriangularMatrixPlot(loop_info, vmin=0, vmax=0.1,
                                        weight_field='fdr_d', colormap='binary_r',
                                        max_dist=600000)

gf = fancplot.GenomicFigure([p_hic, p_oe_d, p_fdr_d], ticks_last=True)
fig, axes = gf.plot('chr11:77.4mb-78.6mb')
# end snippet loops remaining
fig.savefig("../docsrc/api/analyse/images/loops_remaining.png")
plt.close(fig)


# start snippet loops merge
merged_peaks = loop_info.merged_peaks()
# end snippet loops merge

# start snippet loops singlets
singlet_filter = fanc.peaks.RaoMergedPeakFilter(cutoff=-1)
merged_peaks.filter(singlet_filter)
# end snippet loops singlets


# start snippet loops export
# regular Hi-C
merged_peaks.to_bedpe("architecture/loops/rao2014.chr11_77400000_78600000.bedpe")
# chr11	77807143	77842857	chr11	77947143	77982857	.	0.0
# end snippet loops export