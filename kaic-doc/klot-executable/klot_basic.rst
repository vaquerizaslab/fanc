.. _klot-executable:

===========
Basic usage
===========

``klot`` is modular and multiple plots can be combined into a single figure.
Each plot can be configured individually, with a number of customisation options.

There are two modes: by default, an interactive plot opens in a new window,
which allows synchronised scrolling, zooming and panning for plots that support it.
By specifying an output file, the plot is saved to file directly.

********
Overview
********

.. argparse::
   :module: kaic.commands.kaic_commands
   :func: kaic_parser
   :prog: kaic
   :nodescription:

******************************************
Setting up the figure and plotting regions
******************************************

The main argument of ``klot`` is a region specification. You can list one or more
regions by region selector (of the form <chromosome>:<start>-<end>), or you can
provide the path to a file with region information (BED, GFF, Kai-C
object). You can also mix the two. Regions will be plotted in the order they are
listed.

By default, each provided region is exactly the area that is plotted. It is, however,
also possible to merely center the plot on the provided region, and to set a fixed
plotting window using the ``-w`` option. This is especially useful when providing
regions from another analysis, such as ChIP-seq peaks or insulation boundaries,
and allows you to quickly survey other genomic features in their immediate surrounding.

You can also customise the figure properties by additional arguments, controlling
figure name, proportions, and spacing between panels.

*************
Adding panels
*************

After setting up the figure and plotting region(s), you can start adding plots to the
figure. On the command line, the ``-p`` or ``--plot`` argument initiates a new plotting
section. ``klot`` will try to choose a plot type automatically according to the type of
data you provide, but you can always choose exactly which plot you want by using the
``-t`` option.

A basic call to ``klot``, plotting a 2 Megabase region on chromosome 11 of a Hi-C matrix
in classic 2D view, could look like this:

.. code:: bash

    klot chr11:68000000-70000000 -p -t hic2d /path/to/hic_file

To add a BigWig track of ChIP-seq peaks, simply do

.. code:: bash

    klot chr11:68000000-70000000 -p -t hic2d /path/to/hic_file -p -t bigwig /path/to/bigwig_file

which will add a BigWig plot underneath the Hi-C plot.

To list all possible plot types, type:

.. code:: bash

    klot -p -h

which will print:

.. argparse::
   :module: kaic.commands.klot_commands
   :func: type_parser
   :prog: klot


To get more information on a specific plot, simply type:

.. code:: bash

    klot -p -t <plot_type> -h

For example, ``klot -p -t hic -h`` will print the help text for the Hi-C triangle plot:

.. argparse::
   :module: kaic.commands.klot_commands
   :func: hic_parser
   :prog: klot


