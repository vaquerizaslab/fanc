.. _fanc-loops:


############
Loop calling
############


Loops frequently form between two genomic regions, and are visible in the Hi-C matrix as
patches of increased contact intensity:

.. literalinclude:: code/loops_example_code
    :language: bash
    :start-after: start snippet loops example
    :end-before: end snippet loops example

.. image:: images/rao2014.chr11_77400000_78600000.png

We can use ``fanc loops`` to call loops in Hi-C matrices using the HICCUPS algorithm
(`Rao and Huntley et al., 2014 <http://dx.doi.org/10.1016/j.cell.2014.11.021>`_).
Please refer to the original paper for details on the algorithm, specifically the
different types of local neighborhoods defined to make loop calling robust.

.. argparse::
   :module: fanc.commands.fanc_commands
   :func: loops_parser
   :prog: fanc loops
   :nodescription:
   :nodefault:


**********************************
Annotating pixels for loop calling
**********************************

The first step in HICCUPS consists of annotating each pixel with various measures
related to their loop probability. The most important ones are

- Enrichment over expected values in the local neighborhood
- FDR of the local enrichment
- Mappability of the local neighborhood

.. literalinclude:: code/loops_example_code
    :language: bash
    :start-after: start snippet loops annotate
    :end-before: end snippet loops annotate

When run like this, ``fanc loops`` does not actually call any loops, but merely returns
a matrix object where every pixel is annotated with the above properties. Importantly,
as this is the most computationally expensive step, you are strongly advised to choose
a large number of threads using the ``-t`` option. Even better, if you have access to
a computational cluster running Sun/Oracle Grid Engine, ``fanc loops`` can automatically
submit annotation jobs to the cluster if you set the ``--sge`` flag. The ``-t`` option
then specifies the number of jobs allowed to run in parallel instead of the number of
local threads used for multiprocessing.

By default, ``fanc loops`` assumes a loop size of 25kb. This determines the area around
a pixel that is not included in the local neighborhood calculations. If this is chosen
too small, the neighborhood will lie within the peak region, and enrichments are going
to be lower. If this is chosen too big, the neighborhood will no longer be local. If you
have reason to believe your loops size differes from the default, you can set it explicitly
with ``-p``.

Similarly, the width of the neighborhood is determined as ``p + 3`` by default. If you want
to in- or decrease the neighborhood width, use the ``-w`` parameter. You should know,
however, that this is just a starting value, and the neighbborhood width might be increased
on demand internally.

Finally, you can control the size of the submatrices sent to each thread using the
``--batch-size`` parameter. The default, 200, should suit most purposes, but if your
individual jobs are taking too long, you should reduce this number.

We can now use the output object with annotated pixels for downstream processing.

**************************
Filtering annotated pixels
**************************

We need to apply filters to the annotated pixel object that remove all pixels with
a low probability of being a loop. These filters typically consist of enrichment filters,
FDR filters, and mappability filters. Additionally, there are filters for minimum distance
between regions, and the minimum number of unnormalised valid pairs in a pixel.

You can either set a global enrichment filter that acts on all neighborhoods
using ``-e``, or choose individual thresholds for each local neighborhood with
``--enrichment-donut``, ``--enrichment-vertical``, ``--enrichment-horizontal``,
and ``--enrichment-lower-left``. You usually want to set  at least the
``--enrichment-donut`` cutoff to something like 2.

For FDR values, also called q-values, you can set a global filter using ``-q``. Control
the filtering of individual neighborhoods using ``--fdr-donut``, ``--fdr-vertical``,
``--fdr-horizontal``, and ``--fdr-lower-left``. Typical values for each neighborhood
are around 0.1.

Mappability filters act on pixels where a certain fraction of pixels in their local
neighborhoods is unmappable. To set a global mappability cutoff for all neighborhoods,
use the ``-m`` option. Again, local neighborhood mappability filters can be fine-tuned
using the ``--mappability-donut``, ``--mappability-vertical``, ``--mappability-horizontal``,
and ``--mappability-lower-left`` options.

It is generally a good idea to filter on the minimum distance between regions to
consider forming a loop, as a lot of false positive loops will be close to the diagonal.
You can use the ``-d <b>`` parameter to set a threshold on the minimum distance, where
*b* is expressed in number of bins.

In addition, we highly recommend applying a filter on the minimum number of valid pairs
in a pixel (``-o``), so that false-positive loops due to noise are avoided.

Finally, we have included the filter applied by Rao and Huntley et al. in their original
HICCUPS algorithm as a convenient preset ``--rh-filter``. It only retains peaks that are
at least 2-fold enriched over either the donut or lower-left neighborhood, at least
1.5-fold enriched over the horizontal and vertical neighborhoods, at least 1.75-fold
enriched over both the donut and lower-left neighborhood, and have an FDR <= 0.1 in
every neighborhood.

An example command could look like this:

.. literalinclude:: code/loops_example_code
    :language: bash
    :start-after: start snippet loops filter
    :end-before: end snippet loops filter

This filters the vast majority of pixels in the matrix.


************************************
Merging unfiltered pixels into loops
************************************

Pixels that pass all filtering steps are good candidates for loops. Often, these pixels
appear in clusters, which we merge/join in this step. Pixels that do not form a cluster
are generally false-positives, so we filter them using ``--remove-singlets``.

.. literalinclude:: code/loops_example_code
    :language: bash
    :start-after: start snippet loops merge
    :end-before: end snippet loops merge

******************
Exporting to BEDPE
******************

Finally, we can export all the merged loops to BEDPE using ``-b``:

.. literalinclude:: code/loops_example_code
    :language: bash
    :start-after: start snippet loops export
    :end-before: end snippet loops export

