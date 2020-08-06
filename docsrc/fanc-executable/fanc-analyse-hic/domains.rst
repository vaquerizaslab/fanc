.. _fanc-domains:


####################
Hi-C domain analysis
####################

.. note::

    The following examples use the matrix files in FAN-C format. If you want to try the same
    commands using Juicer ``.hic`` files, replace ``output/hic/binned/fanc_example_100kb.hic``
    with ``architecture/other-hic/fanc_example.juicer.hic@100kb``. You will also need to adjust the
    `-vmax` value in the triangular matrix plot to `50`. If you want to work with
    Cooler files in this tutorial, use ``architecture/other-hic/fanc_example.mcool@100kb``.
    The results will be minimally different due to the "zooming" and balancing applied by
    each package.

Like compartments, topologically associating domains, or TADs, for a fundamental level of
genome organisation.

.. literalinclude:: code/domains_example_code
    :language: bash
    :start-after: start snippet domains basic
    :end-before: end snippet domains basic

.. image:: images/fanc_example_100kb_tads.png


FAN-C provides multiple "scores" that are designed to find the boundaries between domains.


****************
Insulation Score
****************

The insulation score (`Crane et al. 2015 <http://www.nature.com/doifinder/10.1038/nature14450>`_)
adds up contacts in a sliding window align the Hi-C matrix diagonal.

.. image:: images/fanc_example_100kb_tads_insulation_score_example.png

Regions with low score are "insulating", i.e. regions between domains. Regions with high scores
are most likely found inside domains.

Use ``fanc insulation`` to calculate the insulation score from the command line:

.. argparse::
   :module: fanc.commands.fanc_commands
   :func: insulation_parser
   :prog: fanc insulation
   :nodescription:
   :nodefault:

=======
Example
=======

``fanc insulation`` is typically used to calculate insulation scores with multiple window
sizes at the same time, as a single window size might be prone to local matrix differences:

.. literalinclude:: code/domains_example_code
    :language: bash
    :start-after: start snippet domains window
    :end-before: end snippet domains window

Window sizes are chosen using the ``-w`` parameter.

We can easily plot all insulation scores at the same time using ``fancplot``:

.. literalinclude:: code/domains_example_code
    :language: bash
    :start-after: start snippet domains scores
    :end-before: end snippet domains scores

.. image:: images/fanc_example_50kb_tads_insulation.png


==============
Output formats
==============

The default output format has maximum compatibility within FAN-C, but other tools like
genome browsers won't be able to read it. Export of insulation scores to a different
format is simple. Just choose ``-output-format bed`` to export to BED file,
``-output-format gff`` to export to GFF format, or ``-output-format bigwig`` to export
to an indexed BigWig file. The latter is a binary file, and hence not readable as
text format, but is the fastest for accessing subsets of scores.

When using an output format other than the default, the second positional argument
(output) becomes an output file prefix. It is appended by the window size (in
abbreviated form, i.e. 1000000 becomes "1mb") and gets the file ending of the chosen
format. Example:

.. literalinclude:: code/domains_example_code
    :language: bash
    :start-after: start snippet domains bed
    :end-before: end snippet domains bed

This produces the output files:

.. code::

    architecture/domains/fanc_example_100kb.insulation_1.5mb.bed
    architecture/domains/fanc_example_100kb.insulation_1mb.bed
    architecture/domains/fanc_example_100kb.insulation_2.5mb.bed
    architecture/domains/fanc_example_100kb.insulation_2mb.bed
    architecture/domains/fanc_example_100kb.insulation_3.5mb.bed
    architecture/domains/fanc_example_100kb.insulation_3mb.bed
    architecture/domains/fanc_example_100kb.insulation_4mb.bed

Of course, you can also simply convert existing insulation scores to another format
without having to recalculate everything. Simply run:

.. literalinclude:: code/domains_example_code
    :language: bash
    :start-after: start snippet domains simplebed
    :end-before: end snippet domains simplebed

and the insulation scores for all window sizes in the object will be converted to BED
files using the input file name as prefix. If you only want to convert specific window
sizes, use the ``-w`` parameter. To find out which window sizes are available in a
previously calculated scores object, simply run ``fanc insulation`` without any
parameters:

.. literalinclude:: code/domains_example_code
    :language: bash
    :start-after: start snippet domains info
    :end-before: end snippet domains info

This prints:

.. code::

    Window sizes available in object:
    1mb 1.5mb 2mb 2.5mb 3mb 3.5mb 4mb


You can plot scores from one or more window sizes using the ``line`` plot in ``fancplot``:

.. literalinclude:: code/domains_example_code
    :language: bash
    :start-after: start snippet domains line
    :end-before: end snippet domains line


.. image:: images/fanc_example_50kb_tads_insulation_1mb.png


=============
Normalisation
=============

By default, ``fanc insulation`` will normalise the insulation scores to the chromosomal
average and the log-transform them. You can get raw, untransformed scores using ``-N`` and
``--L``, respectively. If you want to normalise the scores, but to a smaller region on the
chromosome (to take into account local variability in insulation), you can choose the
normalisation window size with ``--normalisation-window``. The window is specified in bins.

Normally, ``fanc insulation`` will use the arythmetic mean of the chromosomal scores to
normalise. This has the effect that scores upon log2-transformation are not perfectly
centred around 0. To remedy this, you can use the geometric mean instead, with the ``-g``
option.

When you are working with matrices that are already log2-transformed, you may want to use
the ``-s`` option to normalise the scores by subtracting, instead of dividing the chromosomal
average.

If you have a lot of outliers and sharp score changes, you may use a trimmed mean
to calculate average scores with ``--trim-mean <f>``, which will ignore the top and bottom
fraction *f* of scores for calculating the average.


=====================
Impute missing values
=====================

In the above examples, you will notice the region on the left that is unmappable in the Hi-C
matrix. In the insulation score calculation, if the insulation window is covered by more than
50% of unmappable regions, the score will be NaN. ``fanc insulation`` offers the option to
impute the unmappable values from the expected values of the chromosome using ``--impute``.

.. literalinclude:: code/domains_example_code
    :language: bash
    :start-after: start snippet domains impute
    :end-before: end snippet domains impute


This will result in score without NaN (at least in the center of chromosomes), but can also
be misleading if the region of interest happens to lie in an unmappable region. Therefore use
this capability with caution!

.. literalinclude:: code/domains_example_code
    :language: bash
    :start-after: start snippet domains plotimpute
    :end-before: end snippet domains plotimpute

.. image:: images/fanc_example_50kb_tads_insulation_imputed.png


**************************************
Insulating boundaries (TAD boundaries)
**************************************

Regions in the genome where the insulation score reaches a local minimum represent the region
between two self-interacting domains, or TADs. You can use ``fanc boundaries`` to identify these
regions:

.. argparse::
   :module: fanc.commands.fanc_commands
   :func: boundaries_parser
   :prog: fanc boundaries
   :nodescription:
   :nodefault:

When we run ``fanc boundaries`` on the above example using 1mb and 2mb as the window sizes:

.. literalinclude:: code/domains_example_code
    :language: bash
    :start-after: start snippet domains boundaries
    :end-before: end snippet domains boundaries

We get two output files with all insulation score minima and associated scores (the depth of
the minimum compared to the two neighboring maxima):

.. code::

    fanc_example_100kb.insulation_boundaries_1mb.bed
    fanc_example_100kb.insulation_boundaries_2mb.bed

Let's plot the boundaries from the 1mb scores:

.. literalinclude:: code/domains_example_code
    :language: bash
    :start-after: start snippet domains plotboundaries
    :end-before: end snippet domains plotboundaries

.. image:: images/fanc_example_50kb_tads_insulation_1mb_boundaries.png

As you can see, lower minima get higher scores. By default, ``fanc boundaries`` outputs all
minima, but you may set a threshold using ``--min-score <s>`` to report only boundaries with
scores greater than *s*.

.. literalinclude:: code/domains_example_code
    :language: bash
    :start-after: start snippet domains minscore
    :end-before: end snippet domains minscore

By default, ``fanc boundaries`` will return minima as matrix bins. However, since the boundary
calls rely on a smoothed insulation score track, it can attempt to identify the boundary location
with sub-bin resolution. Use ``-x`` to try this, but be aware that this is not precise.


.. _directionality-index:

********************
Directionality Index
********************

The directionality index (`Dixon et al. 2012 <http://www.nature.com/doifinder/10.1038/nature11082>`_)
measures the bias in contact frequency up- and downstream of an Hi-C region. When inside TADs,
this measure tends towards zero, as interactions in either direction are equally frequent. However,
when approaching a TAD boundary this measure changes drastically, as one direction will remain
inside the TAD, where there is a high contact intensity, whereas the other direction will lie in
a low intensity region outside the TAD.

Use ``fanc directionality`` to calculate the directionality index from the command line:

.. argparse::
   :module: fanc.commands.fanc_commands
   :func: directionality_parser
   :prog: fanc directionality
   :nodescription:
   :nodefault:

=======
Example
=======

``fanc directionality`` is very similar in syntax to ``fanc insulation``.
It is typically used to calculate directionality indexes with multiple window
sizes at the same time, as a single window size might be prone to local matrix differences:

.. literalinclude:: code/domains_example_code
    :language: bash
    :start-after: start snippet domains directionality
    :end-before: end snippet domains directionality

Window sizes are chosen using the ``-w`` parameter.

We can easily plot all directionality indexes at the same time using ``fancplot``:

.. literalinclude:: code/domains_example_code
    :language: bash
    :start-after: start snippet domains plotdirectionality
    :end-before: end snippet domains plotdirectionality


.. image:: images/fanc_example_50kb_tads_directionality.png

To export the directionality index to other genomic formats using ``fanc directionality``
follow the instructions as for ``fanc insulation``.


*********************
A note on TAD calling
*********************

There are a lot of tools available for calling TADs in Hi-C matrices, including one that
we have written called `TADtool <https://github.com/vaquerizaslab/tadtool>`_. However,
and this is a point we are also making with TADtool specifically, TAD calling algorithms
often depend critically on their input parameters, and different TAD callers can lead to
very different results. We are therefore currently not bundling a TAD calling tool with
FAN-C, and refer the user to one of the many available tools for TAD calling that offer
a wide range of features.