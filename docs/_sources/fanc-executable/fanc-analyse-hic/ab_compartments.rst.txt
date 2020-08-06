.. _fanc-ab-compartments:

#######################
AB compartment analysis
#######################

.. note::

    The following examples use the matrix files in FAN-C format. If you want to try the same
    commands using Juicer ``.hic`` files, replace ``output/hic/binned/fanc_example_1mb.hic``
    with ``architecture/other-hic/fanc_example.juicer.hic@1mb``. If you want to work with
    Cooler files in this tutorial, use ``architecture/other-hic/fanc_example.mcool@1mb``.
    The results will be minimally different due to the "zooming" and balancing applied by
    each package.

Regions in a Hi-C matrix can generally be assigned to either the active or the inactive
compartment, also called 'A' and 'B' compartments, respectively.

.. argparse::
   :module: fanc.commands.fanc_commands
   :func: compartments_parser
   :prog: fanc compartments
   :nodescription:
   :nodefault:


******************
Correlation matrix
******************

Compartments are derived from a correlation matrix, in which each entry i, j corresponds
to the Pearson correlation between row i and column j of the (Hi-C) matrix.

The ``fanc compartments`` command can produce a correlation matrix (AB compartment) object
from a FAN-C matrix file. Example:

.. literalinclude:: code/ab_compartments_example_code
    :language: bash
    :start-after: start snippet compartments basic
    :end-before: end snippet compartments basic

By default, compartment matrices are calculated on a per-chromosome basis, since each
chromosome might be normalised differently. To force the AB calculation on the whole genome,
use the ``-w`` option.

We can quickly plot the correlation matrix using the ``fancplot`` command:

.. literalinclude:: code/ab_compartments_example_code
    :language: bash
    :start-after: start snippet compartments plot
    :end-before: end snippet compartments plot

.. image:: images/fanc_example_1mb.ab.png


**************
AB Eigenvector
**************

The eigenvector of the correlation matrix is used to derive compartment type and strength
for each matrix bin. Generally, regions with positive values are assigned the 'A',
regions with negative values the 'B' compartment. To calculate the eigenvector in addition
to the correlation matrix, simply add the ``-v <file_name>`` option to the previous command:

.. literalinclude:: code/ab_compartments_example_code
    :language: bash
    :start-after: start snippet compartments ev
    :end-before: end snippet compartments ev

If ``architecture/compartments/fanc_example_1mb.ab`` already exists, it will not be recalculated
but the matrix is loaded from file. You can use the ``-f`` option to overwrite the existing file
in any case.

Alternatively, if you have already calculated the correlation matrix previously, simply run

.. literalinclude:: code/ab_compartments_example_code
    :language: bash
    :start-after: start snippet compartments pre
    :end-before: end snippet compartments pre

Note that I have not specified an output file for the correlation matrix in this case.

We can plot the eigenvector using ``fancplot``:

.. literalinclude:: code/ab_compartments_example_code
    :language: bash
    :start-after: start snippet compartments evplot
    :end-before: end snippet compartments evplot

.. image:: images/fanc_example_1mb.ab_and_ev.png

``fanc compartments`` outputs the first eigenvector by default. In some cases it might be
useful to choose a different eigenvector (sometimes the first eigenvector identifies
chromosomal arms rather than compartments). To change the eiegnvector use the ``--i`` option,
e.g. ``-i 2`` for the second instead of the first correlation matrix eigenvector.

The sign of the eigenvector does not necessarily correspond well to the A or B compartment.
Often, the eigenvector is "flipped" (inverted signs on its entries). Mathematically, if x
is an eigenvector, so is -x. You can use external information to "orient" the eigenvector,
so that it most likely corresponds to the active and inactive compartments. Specifically,
you can supply a FASTA file with the genomic sequence to ``fanc compartments`` using the
``-g <fasta_file>`` argument. This is generally a good idea: ``fanc compartments`` then
calculates the average GC content of regions with positive and those with negative eigenvector
entries. As GC content has previously been shown to correlate well with compartmentalisation,
the eigenvector is oriented in such a way that negative entries correspond to 'B' (low GC
content) and positive entries to 'A' (high GC content).

.. literalinclude:: code/ab_compartments_example_code
    :language: bash
    :start-after: start snippet compartments gc
    :end-before: end snippet compartments gc


**********
AB domains
**********

Consecutive matrix bins with the same eigenvector sign are considered part of a "domain".
You can use the ``-d <domain_file>`` option to write the AB domains to a BED file:

.. literalinclude:: code/ab_compartments_example_code
    :language: bash
    :start-after: start snippet compartments domain
    :end-before: end snippet compartments domain

Again, it is recommended to supply the genome for a better assignment of 'A' and 'B' to
each domain:

.. literalinclude:: code/ab_compartments_example_code
    :language: bash
    :start-after: start snippet compartments gcdomain
    :end-before: end snippet compartments gcdomain

The output contains the domain type in the "name" field and the average eigenvector entry
values of all bins in the domain in the "score" field.

.. code::

    chr18	1	1000000	B	-0.051104836165905	.
    chr18	1000001	2000000	B	-0.08958033099770546	.
    chr18	2000001	3000000	B	-0.08728895584742229	.
    chr18	3000001	4000000	A	-0.04081033542752266	.
    chr18	4000001	5000000	B	-0.060093003511428836	.
    chr18	5000001	6000000	B	-0.06160238136847814	.
    chr18	6000001	7000000	A	-0.04591525399259159	.
    chr18	7000001	8000000	A	-0.028244723100215197	.
    chr18	8000001	9000000	A	-0.010580475545591779	.
    chr18	9000001	10000000	A	0.005653513595461845	.
    ...


**********************
AB enrichment profiles
**********************

The compartmentalisation of a genome can be visualised in an enrichment profile plot.
This will divide bins into percentiles using the associated eigenvector values, and then
the average observed/expected (O/E) values of contacts in each pair of percentile bins
is calculated and plotted.

Use the ``-e <plot_file>`` option to generate this plot. You can also use the
``-m <matrix_file>`` option to output the enrichment matrix values to file for further
analysis.

.. literalinclude:: code/ab_compartments_example_code
    :language: bash
    :start-after: start snippet compartments enrichment
    :end-before: end snippet compartments enrichment

Note that for the enrichment plot you will also need to provide the Hi-C matrix file for
enrichment calculation. If the AB compartment matrix already exists at the output location,
it will not be overwritten unless you are using the ``-f`` option.

As with the above commands, it is highly recommended to supply a genome file to properly
orient the eigenvector:

.. literalinclude:: code/ab_compartments_example_code
    :language: bash
    :start-after: start snippet compartments gcenrichment
    :end-before: end snippet compartments gcenrichment

This is what the plot looks like for the example:

.. image:: images/fanc_example_1mb.ab_profile.png

You can customise the enrichment analysis using additional parameters. By default, the
percentiles for eigenvector binning are chosen at 20, 40, 60, 80, and 100. To choose a
finer binning, for example, you can use ``-p 10 20 30 40 50 60 70 80 90 100``. The "0"
percentile is always added automatically. The eigenvector values are not necessarily
centred around 0, which can lead to some unexpected enrichment plots. To force the central
enrichment matrix bin to 0, and perform separate percentile calculations for values < 0 and
>= 0, use the ``-s 0`` option. Note, however, that this will lead to differences in the
number of bins plotted on the left and right side of the matrix.

.. image:: images/fanc_example_1mb.ab_profile_sym.png

To format the plot, you can choose a different colormap with ``-c <cmap>``. See
`here <https://matplotlib.org/examples/color/colormaps_reference.html>`_ for colormap
options. You can set the colorbar scale with ``--enrichment-min <n>`` and
``--enrichment-max <n>``, where n is -1 and 1, respectively.