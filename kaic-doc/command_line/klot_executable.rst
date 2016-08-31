====
klot
====

Kai-C provides a dedicated plotting executable (``klot``) for plotting Hi-C matrices and other kinds of genomic data.

``klot`` is modular and multiple plots can be combined into a single figure. Each plot can be configured individually,
with a number of customisation options.

There are two modes: by default, an interactive plot opens in a new window, which allows synchronised scrolling, zooming
and panning for plots that support it. By specifying an output file, the plot is saved to file directly.

Here is ``klot`` s help function:

.. code:: bash

    usage:
                    klot [<klot global parameters>] <region> [<region> ...]
                          --plot <plot data file(s)> [<plot parameters>] [...]

                    Run klot --plot -t <plot type> -h for help on a specific subplot.


    klot plotting tool for kaic

    positional arguments:
      regions               List of region selectors (<chr>:<start>-<end>) or
                            files with region information (BED, GTF, ...).

    optional arguments:
      -h, --help            show this help message and exit
      -o OUTPUT, --output OUTPUT
                            Suppresses interactive plotting window and redirects
                            plot to file. Specify path to file when plotting a
                            single region, and path to a folder for plotting
                            multiple regions.
      -s SCRIPT, --script SCRIPT
                            Use a script file to define plot.
      -p PLOT, --plot PLOT  New plot, type will be chosen automatically by file
                            type, unless '-t' is provided.
      -n NAME, --name NAME  Plot name to be used as prefix when plotting multiple
                            regions. Is ignored for single region and interactive
                            plot.
      --height HEIGHT       Height of the figure in inches.
      --width WIDTH         Width of the figure in inches.
      -w WINDOW_SIZE, --window-size WINDOW_SIZE
                            Plotting region size in base pairs. If provided, the
                            actual size of the given region is ignored and instead
                            a region <chromosome>:<region center - window
                            size/2>-<region cener + window size/2> will be
                            plotted.
      -vs HSPACE, --vertical-space HSPACE
                            Vertical distance between plots in fraction of figure.


Setting up the figure and plotting regions
------------------------------------------

The main argument of ``klot`` is a region specification. You can list one or more regions by region selector (of the
form <chromosome>:<start>-<end>), or you can provide the path to a file with region information (BED, GFF, Kai-C
object). You can also mix the two. Regions will be plotted in the order they are listed.

By default, each provided region is exactly the area that is plotted. It is, however, also possible
to merely center the plot on the provided region, and to set a fixed plotting window using the ``-w`` option.
This is especially useful when providing regions from another analysis, such as ChIP-seq peaks or insulation
boundaries, and allows you to quickly survey other genomic features in their immediate surrounding.

You can also customise the figure properties by additional arguments, controlling figure name, proportions, and spacing
between panels.


Adding panels
-------------

After setting up the figure and plotting region(s), you can start adding plots to the figure. On the command line,
the ``-p`` or ``--plot`` argument initiates a new plotting section. ``klot`` will try to choose a plot type
automatically according to the type of data you provide, but you can always choose exactly which plot you want by
using the ``-t`` option.

A basic call to ``klot``, plotting a 2 Megabase region on chromosome 11 of a Hi-C matrix in classic 2D view,
could look like this:

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

.. code:: bash

    usage: klot [-t TYPE] [data [data ...]]

    klot subplot identifier

    positional arguments:
      data                  Data to be plotted in subplot.

    optional arguments:
      -t TYPE, --type TYPE  Manually specify subplot type. Options:
                            hic          Hi-C plot, cropped triangle style
                            hic2d        Hi-C plot, matrix style
                            hicsplit     Hi-C vs Hi-C plot, split matrix
                            hicvhic      Hi-C vs Hi-C plot, matrices "mirrored"
                            fc           Fold-change plot, cropped triangle style
                            hicvfc       Hi-C vs fold-change plot, matrices "mirrored"
                            array        Array "flame" plot (e.g. insulation index)
                            region       Bar plot with region score (e.g. BED)
                            line         Line plot with values per region
                            bigwig       Plot BigWig files

To get more information on a specific plot, simply type:

.. code:: bash

    klot -p -t <plot_type> -h

For example, ``klot -p -t hic -h`` will print the help text for the Hi-C triangle plot:

.. code:: bash

    usage: klot [-h] [--aspect-ratio ASPECT_RATIO] [--title TITLE]
                [--fix-chromosome] [-vmin VMIN] [-vmax VMAX] [-d MAX_DIST] [-l]
                [-r] [-c COLORMAP] [-C]
                hic

    Hi-C plot.

    positional arguments:
      hic                   Hi-C object.

    optional arguments:
      -h, --help            show this help message and exit
      --aspect-ratio ASPECT_RATIO
                            Aspect ratio of this panel. Default is determined by
                            figure type (usually 1.0).
      --title TITLE         Title of this plot.
      --fix-chromosome      Fix chromosome identifier for this plot (add or remove
                            'chr' as required)
      -vmin VMIN, --minimum-value VMIN
                            Minimum value assigned the first color in the
                            colorbar.
      -vmax VMAX, --maximum-value VMAX
                            Maximum value assigned the last color in the colorbar.
      -d MAX_DIST, --maximum-distance MAX_DIST
                            Maximum distance between two points after which
                            triangle will be truncated.
      -l, --log             Log-transform heatmap values
      -r, --range-slider    Add vmax/vmin slider to plot
      -c COLORMAP, --colormap COLORMAP
                            Matplotlib colormap (default: viridis)
      -C, --no-colorbar     Do not show colorbar in plot


Plot types
----------

TODO: List all plot types with sample images.
