"""
Provide plotting functions for genomic data types.

Many common data types used in genomics research are supported. Including, but
not limited to, :class:`Hi-C <fanc.plotting.HicPlot>`,
:class:`bed <fanc.plotting.GenomicFeaturePlot>`,
:class:`bigwig <fanc.plotting.BigWigPlot>`
and :class:`gene (GTF) <fanc.plotting.GenePlot>` file visualization.
The basic idea is that figures can be composed of multiple panels which are
arranged vertically and share a common x-axis representing genomic coordinates.

Each panel is created separately and then combined into a single figure.
For example, when analyzing Hi-C data it is often interesting to correlate
features in the Hi-C map with ChIP-seq tracks. In that case, one would first
create a :class:`~fanc.plotting.HicPlot` object, which visualizes
Hi-C data, and then a :class:`~fanc.plotting.BigWigPlot` object, which
can plot bigwig files that are used during ChIP-seq analysis. Finally, the two
objects are used to create a :class:`~fanc.plotting.GenomicFigure`.

Example
-------

FAN-C comes with a few example datasets that can be used to explore the basic
functionalities of the plotting module. The paths in this example are relative
to the top-level fanc directory where the setup.py file is located.

.. note::
    The paths to the example datasets can be accessed easily using the
    example_data dictionary:

    .. code:: python

        import fanc
        print(fanc.example_data)

    Yields a dict with the user-specific path to the datasets:

    ::

        {'chip_bedgraph': '/Users/example/fanc/fanc/test/data/test_plotting/CTCF_ChIP_FE_chr11_77-80Mb_mouse_embryo_fibroblasts.bedgraph.gz',
         'chip_bigwig': '/Users/example/fanc/fanc/test/data/test_plotting/CTCF_ChIP_FE_chr11_77-80Mb_mouse_embryo_fibroblasts.bigwig',
         'chip_peak_bed': '/Users/example/fanc/fanc/test/data/test_plotting/CTCF_ChIP_FE_chr11_77-80Mb_mouse_embryo_fibroblasts.peaks.bed.gz',
         'gene_gtf': '/Users/example/fanc/fanc/test/data/test_plotting/genes_mm10_chr11_77-80Mb.gtf.gz',
         'hic': '/Users/example/fanc/fanc/test/data/test_network/rao2014.chr11_77400000_78600000.hic'}

.. code:: python

    import fanc.plotting as fancplot

    # Create Hic plot
    hplot = fancplot.HicPlot("fanc/test/data/test_network/rao2014.chr11_77400000_78600000.hic")

    # Create plot showing some CTCF ChIP-seq data from ENCODE
    bplot = fancplot.BigWigPlot("fanc/test/data/test_plotting/CTCF_ChIP_FE_chr11_77-80Mb_mouse_embryo_fibroblasts.bigwig",
                             title="CTCF ChIP", ylabel="fold enrichment")

    # Create plot of all genes in the region. The squash option in combination
    # with the group_by option causes all exons of each gene to be merged.
    # This is useful if the number of alternative transcripts is overwhelming
    # the plot.
    gplot = fancplot.GenePlot("fanc/test/data/test_plotting/genes_mm10_chr11_77-80Mb.gtf.gz",
                           group_by="gene_name", squash=True, show_labels=False,
                           title="mm10 genes")

    # The created plots are used to generate a figure by passing them as a list
    # to the GenomicFigure constructor. The order in which they are passed
    # determines the order of panels in the figure.
    gfig = fancplot.GenomicFigure([hplot, bplot, gplot])

    # Plot a specific region of the genome
    fig, axes = gfig.plot("chr11:77400000-78600000")

    # Open plot in an interactive window
    fig.show()

    Example rendering of the above code saved using
    ``fig.savefig("example.png", dpi=100)``.

Editing figure and axes
-----------------------

The :meth:`GenomicFigure.plot() <fanc.plotting.GenomicFigure.plot>`
function returns standard matplotlib Figure and a list of Axes instances that
can be further adjusted using standard matplotlib methods. The matplotlib axes
instance associated with each plot is also accesible from the "ax" property of
each plot.

.. warning:: The Axes instances of the plot should only be edited after
    :meth:`GenomicFigure.plot() <fanc.plotting.GenomicFigure.plot>`
    has been called. Otherwise any changes that were made may be overwritten
    when the plot() method is called.

For example, to add a bit of annotating text at a specific location in the
BigWigPlot, the example above can be edited as follows:

.. code:: python

    fig, axes = gfig.plot("chr11:77400000-78600000")
    bplot.ax.text(77850000, 60, "Interesting peak!")
    fig.show()

The coordinates in the Axes are data coordinates, the x-axis is genomic
coordinates on the current chromosome and the y-axis in this case the
fold-enrichment of the bigwig track.

Basic Plot types and options
----------------------------

.. note::
    An explanation of each plot class and the parameters that it supports can be
    accessed by suffixing a question mark (in Ipython/Jupyter) or calling the
    help() function:

    .. code::

        import fanc.plotting as fancplot
        fancplot.BigWigPlot? # Ipython/Jupyter
        help(fancplot.BigWigPlot) # standard python


A few basic parameters such as a title and the aspect ratio are available for
all plot classes. The aspect ratio parameter is a floating point number between
0 and 1 that determines the height of the plot. A value of 1 results in a square
plot, .5 represents a plot that is half as high as it is wide.

The :class:`~fanc.plotting.GenomicFigure` provides a few convenience
parameters. Setting ``ticks_last=True`` for example removes tick labels from all
panels but the last one which makes the overall plot more compact.

Independent x-axis and inverting x-axis
'''''''''''''''''''''''''''''''''''''''

By default, the x-axis of all plots in a figure are linked, meaning that all
plots display exactly the same region. In some situations in can be helpful
to plot a multiple regions, such as when features in syntenic regions across
multiple species need to be compared. Since syntenic regions can be on the +
strand in some species and on the - strande, sometimes the x-axis also needs to
be inverted to maintain correct orientation.

In this situation the ``independent-x`` option should be set in the
:class:`~fanc.plotting.GenomicFigure`. As a result, the
:meth:`GenomicFigure.plot() <fanc.plotting.GenomicFigure.plot>` method
no longer expects a single region as argument, but a list of regions equal to
the number of plots in the figure. We can modify the example above to illustrate
this point:

.. code:: python

    gfig = fancplot.GenomicFigure([hplot, bplot, gplot], independent_x=True)
    fig, axes = gfig.plot(["chr11:77400000-78600000", "chr11:77500000-78600000",
                           "chr11:77200000-78600000"])
    fig.show()

Synchronize y-axis limits for multiple datasets
'''''''''''''''''''''''''''''''''''''''''''''''

Sometimes it can be useful to synchronize y-axis limits across multiple datasets
in order to compare their signals. Many plot types support plotting multiple
datasets in the same panel, making comparisons easy:

.. code:: python

    bplot = fancplot.BigWigPlot(["dataset1.bigwig", "dataset2.bigwig"])

Alternatively it is possible to synchronize y-axis limits across panels using
:class:`~fanc.plotting.LimitGroup` instances. In this example the y-axis limits
are shared between bplot1 and bplot2 and seperately between bplot3 and bplot4:

.. code:: python

    ygroup = fancplot.LimitGroup()
    ygroup2 = fancplot.LimitGroup(limit=(None, 10))
    bplot1 = fancplot.BigWigPlot("dataset1.bigwig", ylim=ygroup)
    bplot2 = fancplot.BigWigPlot("dataset2.bigwig", ylim=ygroup)
    bplot3 = fancplot.BigWigPlot("dataset3.bigwig", ylim=ygroup2)
    bplot4 = fancplot.BigWigPlot("dataset4.bigwig", ylim=ygroup2)

It is also possible to constrain the y-axis limits. Passing ``limit=(None, 10)``
to the constructor constrains the upper limit to a maximum of 10 while leaving
the lower limit unconstrained.

Programmatic plotting using loops
---------------------------------

fanc plotting is ideally suited for programmatic generation of many plots or
dynamic assembly of multiple datasets in a single figure. In this example three
Hi-C datasets are visualized in a single figure:

.. code:: python

    import fanc.plotting as fancplot
    hic_datasets = ["my_data1.hic", "my_data2.hic", my_data3.hic"]
    hic_plots = [fancplot.HicPlot(h, max_dist=500000) for h in hic_datasets]
    gfig = fancplot.GenomicFigure(hic_plots)
    regions = ["chr11:77400000-78600000", "chr11:1100000-13600000"]
    for r in regions:
        fig, axes = gfig.plot(r)
        fig.savefig("plot_region_{}.png".format(r.replace(":", "_")))

"""

from fanc import config
from fanc.plotting.hic_plotter import HicPlot, HicPlot2D, HicComparisonPlot2D, \
    HicSlicePlot, HicPeakPlot, TriangularMatrixPlot, SquareMatrixPlot, SplitMatrixPlot
from fanc.plotting.plotter import VerticalSplitPlot, GenomicVectorArrayPlot, GenomicFeaturePlot, GenomicRegionsPlot, \
    GenomicFeatureScorePlot, GenomicFigure, BigWigPlot, GenePlot, LinePlot, \
    FeatureLayerPlot, GenomicDataFramePlot, HighlightAnnotation, RegionsValuesPlot, \
    BarPlot, Virtual4CPlot, MirrorMatrixPlot
from fanc.plotting.helpers import append_axes, absolute_wspace_hspace, SymmetricNorm, \
                                  style_ticks_whitegrid, LimitGroup
from fanc.plotting.statistics import *
from fanc.plotting.colormaps import *
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import logging

logger = logging.getLogger(__name__)

if config['pdf_font_as_text']:
    logger.debug("Using text for PDFs instead of paths")
    matplotlib.rcParams['pdf.fonttype'] = 42

sns.set_style("ticks")

plt.register_cmap(name='RdBuWhitespace_r', cmap=fc_cmap)
plt.register_cmap(name='germany', cmap=germany_cmap)
plt.register_cmap(name='germany_r', cmap=germany_r)
plt.register_cmap(name='white_red', cmap=white_red)
plt.register_cmap(name='white_red_r', cmap=white_red_r)
