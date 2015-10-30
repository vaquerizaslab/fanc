Getting started
===============

.. contents::
   :depth: 2

kaic offers access to its Hi-C processing pipeline on multiple levels, including
a high-level executable and a low-level Python 2.7 API.

Installation
------------

You can download the kaic source code from GitLab by cloning its repository.

.. code:: bash

  git clone

kaic depends on several python packages which should mostly be installed
automatically during the installation process.

However, one core dependency is not installed automatically:
`PyTables <https://github.com/PyTables/PyTables>`_ is currently on version 3.2.2
and there is a critical bug affecting kaic performance that will be fixed in
version 3.3. To circumvent the bug, it is highly recommended to install the
PyTables development version.

PyTables depends on the `HDF5 <https://www.hdfgroup.org/HDF5/>`_ library. To
install the latest version it is easiest to use Homebrew (on OS X) or Linuxbrew.

.. code:: bash

  brew install hdf5

Or you can install
`from source <https://www.hdfgroup.org/HDF5/release/obtain5.html>`_.

To ensure this version of HDF5 is used during the installation of PyTables, set
the HDF5_DIR environment variable to the path of your HDF5 installation.

Then install PyTables:

.. code:: bash

  pip install git+https://github.com/PyTables/PyTables.git

or if you do not have root access to the Python installation:

.. code:: bash

  pip install --user git+https://github.com/PyTables/PyTables.git

Finally, from the root directory of kaic do

.. code:: bash

  pip install .

or

.. code:: bash

  pip install --user .



Command-line executable
-----------------------

kaic provides a high-level executable that can perform most kaic functions. Here
is its help screen:

.. code:: bash

    usage: kaic <command> [options]

        Commands:

        --- Mapping
        iterative_mapping  Iteratively map a FASTQ file to a Bowtie 2 index

        --- Reads
        load_reads         Load a SAM/BAM file into a Reads object
        filter_reads       Filter a Reads object

        --- Pairs
        reads_to_pairs     Convert a Reads object into a Pairs object
        filter_pairs       Filter a Pairs object

        --- Hic
        pairs_to_hic       Convert a pairs object into a Hic object
        merge_hic          Merge multiple Hic objects
        bin_hic            Bin a Hic object into same-size regions
        correct_hic        Correct a Hic object for biases

        --- Plotting
        plot_ligation_err  Plot the ligation error of a Pairs object
        plot_hic_corr      Plot the correlation of two Hic objects
        plot_hic_matrix    Plot a Hic matrix
        plot_diff          Plot the difference between two Hic matrices
        plot_ratio         Plot the ratio between two Hic matrices

        --- Other
        hiclib_to_kaic     Convert a hiclib object to a Hic object

    Run kaic <command> -h for help on a specific command.


Mapping
~~~~~~~


A kaic pipeline generally expects to start from mapped reads in a SAM/BAM file. For convenience, however, kaic offers a
single mapping command: iterative_mapping.

iterative_mapping
_________________

It currently requires the `hiclib package <https://bitbucket.org/mirnylab/hiclib>`_ to be installed.

.. code:: bash

    usage: kaic iterative_mapping [-h] [-b BOWTIE] [-bo BOWTIE2_OPTIONS]
                              [-m MIN_LENGTH] [-s STEP] [-t THREADS]
                              [-tmp TMP] [-nj] [-nc]
                              input index output

    Iteratively map a FASTQ file to a Bowtie 2 index

    positional arguments:
      input                 Input FASTQ files (comma-separated)
      index                 Bowtie 2 index (include prefix)
      output                Output files (must be same number as input files)

    optional arguments:
      -h, --help            show this help message and exit
      -b BOWTIE, --bowtie BOWTIE
                            Bowtie 2 executable path (will check PATH variable by
                            default)
      -bo BOWTIE2_OPTIONS, --bowtie2-options BOWTIE2_OPTIONS
                            Bowtie 2 command line options, default: --very-
                            sensitive --no-unal
      -m MIN_LENGTH, --min-length MIN_LENGTH
                            Minimum sequence length to attempt the mapping
      -s STEP, --step STEP  Step size to increase mapped sequence length
      -t THREADS, --threads THREADS
                            Number of threads
      -tmp TMP, --temp-dir TMP
                            Temporary directory
      -nj, --no-join        Do not join partial SAM files into a single file
      -nc, --no-clean       Do not delete partial SAM files (*.sam.\d)

iterative_mapping will truncate reads in a FASTQ file to a shorter size (given by MIN_LENGTH) and map them to a Bowtie 2
reference genome. Reads that do not map uniquely to the reference are collected, extended by STEP, and mapped again to
the reference genome. This process is repeated until non-uniquely mapping reads are extended to their full length.
Finally, all uniquely mapping reads are joined in a single SAM file.

Example use:

.. code:: bash

    kaic iterative_mapping /path/to/some.fastq /path/to/bowtie/index/prefix /path/to/output.sam -m 50 -s 5


Reads
~~~~~

Reads objects represent a list of mapped reads. kaic offers functionality to load reads from a SAM/BAM file and to
filter reads based on several mapping criteria.

load_reads
__________

This command loads reads from a SAM file along with all their mapping properties.

.. code:: bash

    usage: kaic load_reads [-h] input output

    Load a SAM/BAM file into a Reads object

    positional arguments:
      input       Input SAM file
      output      Output file

The result is a Reads object, by convention these should have the ``.reads`` extension.

Example use:

.. code:: bash

    kaic load_reads /path/to/some.sam /path/to/output.reads


filter_reads
____________

This command can be used to filter reads in a Reads object that do not pass certain criteria.

.. code:: bash

    usage: kaic filter_reads [-h] [-m] [-u] [-us] [-q QUALITY] [-s STATS]
                         input [output]

    Filter a Reads object

    positional arguments:
      input                 Input Reads file
      output                Output Reads file. If not provided will filter
                            existing file directly.

    optional arguments:
      -h, --help            show this help message and exit
      -m, --mapped          Filter unmapped reads
      -u, --unique          Filter reads that map multiple times (with a lower
                            score)
      -us, --unique-strict  Strictly filter reads that map multiple times (XS tag)
      -q QUALITY, --quality QUALITY
                            Cutoff for the minimum mapping quality of a read
      -s STATS, --stats STATS
                            Path for saving stats pdf

The ``-m`` option filters out all unmapped reads. The ``-u`` option filter reads with duplicate alignments of the same
quality to the reference genome, while ``-us`` filters reads if they have duplicate alignments regardless of quality.
With ``-q QUALITY`` it is possible to filter reads with a mapping quality lower than ``QUALITY``.

By adding the ``-s STATS`` option it is possible to get a PDF overview of the filtering process in a simple bar chart:

.. image:: images/reads.stats.png

Example use:

.. code:: bash

    kaic filter_reads /path/to/original.reads /path/to/filtered.reads -m -us -q 30 -s /path/to/stats.pdf


Pairs
~~~~~

A Pairs object represents pairs of mapped reads that have been assigned to regions in a reference genome. Typically,
regions are restriction fragments, which mark the lowest achievable resolution in a Hi-C experiment.

reads_to_pairs
______________

This command converts two (paired) Reads objects to a Pairs object by first identifying the genomic region each read
falls in, and then saving matching pairs of reads. It requires a reference sequence in FASTA format and the name of the
restriction enzyme used in the experiment.

.. code:: bash

    usage: kaic reads_to_pairs [-h] reads1 reads2 genome restriction_enzyme output

    Convert a Reads object into a Pairs object

    positional arguments:
      reads1              First half of input reads
      reads2              Second half of input reads
      genome              Can be an HDF5 Genome object, a FASTA file, a folder
                          with FASTA files, or a comma-separated list of FASTA
                          files.
      restriction_enzyme  Restriction enzyme used in the experiment, e.g. HindIII
      output              Output file for mapped pairs

The ``genome`` parameter is very flexible in its usage: To ensure that the regions in the final Hic object occur in the
desired order, it is recommended to use a comma-separated string with the paths of FASTA files with each chromosome
reference sequence.

Example:

.. code:: bash

    kaic reads_to_pairs /path/to/first.reads /path/to/second.reads /path/to/chr1.fa,/path/to/chr2.fa HindIII /path/to/output.pairs


filter_pairs
____________

Similar to ``filter_reads``, this command filters pairs of mapped reads in a Pairs object.

.. code:: bash

    usage: kaic filter_pairs [-h] [-i INWARD] [-o OUTWARD] [-r REDIST] [-s STATS]
                         input [output]

    Filter a Pairs object

    positional arguments:
      input                 Input FragmentMappedPairs file
      output                Output FragmentMappedPairs file. If not provided will
                            filter input file in place.

    optional arguments:
      -h, --help            show this help message and exit
      -i INWARD, --inward INWARD
                            Minimum distance for inward-facing read pairs
      -o OUTWARD, --outward OUTWARD
                            Minimum distance for outward-facing read pairs
      -r REDIST, --re-distance REDIST
                            Maximum distance for a read to the nearest restriction
                            site
      -s STATS, --stats STATS
                            Path for saving stats pdf

The ``-i`` option can be used to filter *inward-facing* read pairs, while ``-o`` filter *outward-facing* reads at a
certain distance (see `Jin et al. 2013 <http://www.nature.com/nature/journal/v503/n7475/full/nature12644.html>`_).
``-r`` filters pairs where at least one read maps more than a certain distance to the nearest restriction site.

Example use:

.. code:: bash

    kaic filter_pairs /path/to/original.pairs /path/to/filtered.pairs -i 10000 -o 25000 -s /path/to/stats.pdf


Hic
~~~

The Hic object represents a Hi-C matrix. This includes both variable-region matrices, such as those based on restriction
fragments, and equi-distant regions, such as binned Hi-C matrices. It handles common tasks, such as binning or merging
Hic objects, and can be used to `plot <Plotting>`_ Hi-C data in a variety of ways.


pairs_to_hic
____________

This command converts a Pairs object into a Hic object by summing up pairs with the same fragments and using that as
a weight (or contact count). The regions defined in the pairs object are transferred to the new object without changes,
i.e. the order of regions as defined in the `reads_to_pairs`_ command will be the order of regions along the axes of the
Hi-C matrix.

.. code:: bash

    usage: kaic pairs_to_hic [-h] pairs hic

    Convert a pairs object into a Hic object

    positional arguments:
      pairs       Input FragmentMappedReadPairs file
      hic         Output path for Hic file

Example:

.. code:: bash

    kaic pairs_to_hic /path/to/my.pairs /path/to/new.hic


merge_hic
_________

Merges multiple Hic objects into one. The command will try to merge smartly, i.e. it should even work in cases where the
genomic regions differ between objects (for example when merging a chr1 with a chr2 matrix). In a first step, regions
will be merged and regions that exist in both matrices will be assigned new indices. In the second step, contacts will
be merged.

.. code::bash

    usage: kaic merge_hic [-h] hic [hic ...] output

    Merge multiple Hic objects

    positional arguments:
      hic         Input Hic files
      output      Output binned Hic object

Example:

.. code:: bash

    kaic merge_hic /path/to/old_1.hic /path/to/old2.hic /path/to/old3.hic /path/to/merged.hic


bin_hic
_______

This command bins regions in the genome into same-size chunks. The default strategy to distribute reads in the case of
old regions overlapping two or more regions in the binned Hic object is given by
`Rao et al. (2014) <http://www.cell.com/abstract/S0092-8674%2814%2901497-4>`_. Please note that, due to the nature of
the binning strategy, it is very likely that the last region in the genome is shorter than the requested bin size.

.. code:: bash

    usage: kaic bin_hic [-h] hic output bin_size

    Bin a Hic object into same-size regions

    positional arguments:
      hic         Input Hic file
      output      Output binned Hic object
      bin_size    Bin size in base pairs

Example to bin an existing object at 50kb resolution:

.. code:: bash

    kaic bin_hic /path/to/old.hic /path/to/binned.hic 50000


correct_hic
___________

You can use this command to correct Hic matrices using matrix balancing. By default, it uses the efficient matrix
balancing approach by `Knight and Ruiz (2012) <http://imajna.oxfordjournals.org/content/33/3/1029>`_, but providing the
``-i`` option switches to the iterative ICE method by
`Imakaev et al. (2012) <http://www.nature.com/nmeth/journal/v9/n10/full/nmeth.2148.html?WT.ec_id=NMETH-201210>`_.

.. code:: bash

    usage: kaic correct_hic [-h] [-i] [-c] input [output]

    Correct a Hic object for biases

    positional arguments:
      input             Input Hic file
      output            Output Hic file. If not provided will filter existing file
                        in place.

    optional arguments:
      -h, --help        show this help message and exit
      -i, --ice         Use ICE iterative correction instead of Knight matrix
                        balancing
      -c, --chromosome  Correct intra-chromosomal data individually, ignore inter-
                        chromosomal data

Sometimes it is not wanted to correct the entire matrix in one go, for example due to computer memory constraints or
the quality of inter-chromosomal data. In this case the ``-c`` option will cause the command to correct each
intra-chromosomal sub-matrix individually, leaving the inter-chromosomal data untouched.

Example use:

.. code:: bash

    kaic correct_hic /path/to/uncorrected.hic /path/to/corrected.hic


Plotting
~~~~~~~~

kaic provides a growing list of plotting commands to quickly assess the data at hand.

plot_ligation_err
_________________

Plot the ligation error of mapped read pairs in a Pairs object. For an explanation of the different types of read pairs
see `Jin et al. (2013) <http://www.nature.com/nature/journal/v503/n7475/full/nature12644.html>`_. The point at which
the red and blue curves converge toward the dotted line can be used as a rough guideline for cutoffs in the
`filter_pairs`_ command.

.. code:: bash

    usage: kaic plot_ligation_err [-h] [-p POINTS] input [output]

    Plot the ligation error of a Pairs object

    positional arguments:
      input                 Input FragmentMappedPairs file
      output                Output pdf

    optional arguments:
      -h, --help            show this help message and exit
      -p POINTS, --points POINTS
                            Data points that make up one increment of the x axis.
                            More=smoother=less detail.

``-p POINTS`` can be used to control the smoothing of the curve, but generally the auto-selected value provides a good
balance between smooting and detail.

Example:

.. code:: bash

    kaic plot_ligation_err /path/to/my.pairs /path/to/error.pdf


plot_hic_matrix
_______________

Plot the Hi-C matrix represented by a Hic object. By default, the command tries to pick the