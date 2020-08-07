.. _fanc-aggregate:


#######################
Hi-C aggregate analysis
#######################

.. note::

    The following examples use the matrix files in FAN-C format. If you want to try the same
    commands using Juicer ``.hic`` files, replace ``output/hic/binned/fanc_example_100kb.hic``
    with ``architecture/other-hic/fanc_example.juicer.hic@100kb``. If you want to work with
    Cooler files in this tutorial, use ``architecture/other-hic/fanc_example.mcool@100kb``.
    The results will be minimally different due to the "zooming" and balancing applied by
    each package.

It can be very informative to view the average Hi-C matrix for a set of regions, rather
than the Hi-C matrix at each individual region. This can help you in identifying common
structural features of these regions or over-/under-representation of contacts in the
vicinity.

Here are examples of TAD and loop aggregate plots from our recent preprint
(`Kruse et al. (2019) <https://www.biorxiv.org/content/10.1101/523712v1>`_):

.. image:: images/aggregate_examples_biorxiv.png

You can easily create your own aggregate plots using ``fanc aggregate``.

.. argparse::
   :module: fanc.commands.fanc_commands
   :func: aggregate_parser
   :prog: fanc aggregate
   :nodescription:
   :nodefault:

You can provide ``fanc aggregate`` with a list of genomic regions in any of the common
region-based formats (BED, GFF, BigWig, ...) or with a list of genomic region pairs in
BEDPE format. For lists of regions, the aggregate matrix will be located at the Hi-C
matrix diagonal. For pairs of regions, matrix subsets can be anywhere in the genome.

************************************
Aggregate over variable size regions
************************************

By default, if you provide ``fanc aggregate`` with a list of regions, it will extract
the square Hi-C sub-matrices along the diagonal for each region and interpolate them
to match the width set by ``--pixels`` (90 by default). It will then calculate the
average value for each pixel, which then form the aggregate matrix.

Let's try this on TADs called using the arrowhead algorithm (`Rao and Huntley et al.,
2014 <http://dx.doi.org/10.1016/j.cell.2014.11.021>`_). ``fanc aggregate`` will ignore
all regions in the file that are not present in the Hi-C matrix. In our example Hic file,
that is everything outside of chromosomes 18 and 19:

.. literalinclude:: code/aggregate_example_code
    :language: bash
    :start-after: start snippet aggregate basic
    :end-before: end snippet aggregate basic

This command only produces an AggregateMatrix file (``fanc_example_100kb.agg``), which
is useful for further usage with FAN-C, but not easily readable. To extract the aggregate
matrix in txt format, simply add ``-m`` and to plot it just use ``-p``:

.. literalinclude:: code/aggregate_example_code
    :language: bash
    :start-after: start snippet aggregate plot
    :end-before: end snippet aggregate plot

The resulting plot looks a bit strange:

.. image:: images/fanc_example_100kb.agg.png


Important note: if your input regions have variable sizes, as assumed in this section,
the resulting aggregate matrix with default settings is highly misleading, since smaller
regions will have larger average signal in the interpolated matrix due to being closer to
the diagonal. You can easily correct for this effect using O/E matrices instead of the
regular Hi-C matrix. Simply set the ``-e`` flag for this. ``-e`` works very well with
log2-transformed data (``-l``). Let's see how this changes your matrix:

.. literalinclude:: code/aggregate_example_code
    :language: bash
    :start-after: start snippet aggregate oe
    :end-before: end snippet aggregate oe

.. image:: images/fanc_example_100kb_oe.agg.png

This still does not look like much of a TAD, but we can add a little more context by
expanding the plotting region relative to the region size using ``-r``:

.. literalinclude:: code/aggregate_example_code
    :language: bash
    :start-after: start snippet aggregate expand
    :end-before: end snippet aggregate expand

.. image:: images/fanc_example_100kb_oe_large.agg.png


That plot depicts a region that is 3x the size of the TAD located in its center and
already looks like we would expect: High signal in the center, especially at the TAD
corner, where the corner loops are typically located.

We can further apply an exponential rescaling (``--rescale``) of the data to make this
look more like a Hi-C matrix:

.. literalinclude:: code/aggregate_example_code
    :language: bash
    :start-after: start snippet aggregate rescale
    :end-before: end snippet aggregate rescale

.. image:: images/fanc_example_100kb_oe_large_res.png

Here, we are not log-transforming the data and we are setting the saturation of the
pixel values at 0.045 using ``--vmax``.


**************
Aggregate TADs
**************

For both the log2(O/E) and rescaled versions of the aggregate matrices, there are
preset flags you can use called ``--tads`` and ``--tads-imakaev``, respectively. The
latter is named after the first author of the publication that first used rescaled
aggregate matrices in this fashion
(`Flyamer, Gassler, and Imakaev et al., 2017 <http://www.nature.com/doifinder/10.1038/nature21711>`_). In
the above example, you can simply run

.. literalinclude:: code/aggregate_example_code
    :language: bash
    :start-after: start snippet aggregate tads
    :end-before: end snippet aggregate tads


*******************
Fixed-width regions
*******************

Sometimes, you may want to use a fixed window surrounding a set of features in the
aggregate analysis, such as TAD boundaries. ``fanc aggregate`` provides the ``-w``
option to plot the aggregate Hi-C matrix in a window os size w around the center
of each region in the list provided.

.. literalinclude:: code/aggregate_example_code
    :language: bash
    :start-after: start snippet aggregate fixed
    :end-before: end snippet aggregate fixed

.. image:: images/fanc_example_100kb_boundaries.agg.png


You can see the relatively faint "average boundary" in the centre of the plot. When using
O/E and log2-transformed matrices, this becomes much more obvious:

.. literalinclude:: code/aggregate_example_code
    :language: bash
    :start-after: start snippet aggregate oefixed
    :end-before: end snippet aggregate oefixed

.. image:: images/fanc_example_100kb_boundaries_oe.agg.png


You can change the viewpoint to other positions within a region, such as the 5' end,
using the ``-v`` option.


****************************************
Loops and other pairwise genomic regions
****************************************

When you have loop calls or other pairwise genomic regions in BEDPE format, you can use
``fanc aggregate`` to make aggregate loop plots. The preset for this is ``--loops``.

.. literalinclude:: code/aggregate_example_code
    :language: bash
    :start-after: start snippet aggregate loops
    :end-before: end snippet aggregate loops

.. image:: images/rao2014.chr11_77400000_78600000.loops_no_singlets.agg.png

Control the size of the plot using the ``--pixels`` argument.