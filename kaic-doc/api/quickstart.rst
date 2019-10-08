.. _quickstart:

=============================
Quickstart with the Kai-C API
=============================

After you have installed Kai-C (see :ref:`kaic_installation`), you can import the module
from a Python console or script:

.. code::

    import kaic

The following tutorials will assume that you have loaded the :code:`kaic` module in this
manner.

****************
Loading datasets
****************

Any analysis typically begins with loading datasets into your workspace. Kai-C tries to make this
as simple as possible with the :func:`~kaic.__init__.load` function.
If you already have processed Hi-C files, either from the :code:`kaic` command line
application (see :ref:`kaic-auto` or :ref:`kaic-modular`), or from a compatible Hi-C application
(:code:`.cool` or :code:`.mcool` from `Cooler <https://github.com/mirnylab/cooler>`_ or :code:`.hic` from
`Juicer <https://github.com/aidenlab/juicer>`_), simply load them into your workspace
using :func:`~kaic.__init__.load` - no need to specify the type of file you are loading:

.. code::

     data = kaic.load("/path/to/file.hic")

When dealing with multi-resolution Hi-C files such as :code:`.mcool` from Cooler or :code:`.hic` from
Juicer, you can load a specific resolution using the :code:`@` notation:

.. code::

     data = kaic.load("/path/to/file.mcool@25000")

:func:`~kaic.__init__.load` is not limited to Hi-C files, but also works on any other file
produced with Kai-C, such as :class:`~kaic.pairs.RegionPairs` files, :class:`~kaic.regions.Genome`,
analysis results like :class:`~kaic.architecture.comparisons.FoldChangeMatrix` and generally most
other Kai-C files.

:func:`~kaic.__init__.load` even works on most of the common file formats for genomic
datasets, such as BED, GFF, BigWig, Tabix, BEDPE and more. Try it out on your dataset of choice and
chances are :func:`~kaic.__init__.load` can handle it. And if it does not, consider raising
an issue on `Github <https://github.com/vaquerizaslab/kaic/issues>`_ to ask for support.

Internally, :func:`~kaic.__init__.load` finds a suitable class for the type of data
in the supplied file, and opens the file using that class. For example, the result of

.. code::

    hic = kaic.Hic("output/hic/binned/kaic_example_1mb.hic")

is equivalent to

.. code::

    hic = kaic.load("output/hic/binned/kaic_example_1mb.hic")

with the big advantage that you don't need to worry about remembering class names or
their location within the Kai-C module hierarchy. In both cases, the type of the
returned object is :code:`kaic.hic.Hic`:

.. code::

    # check the type fo the object
    type(hic)  # kaic.hic.Hic

Here are a few more examples:

.. code::

    cool = kaic.load("test.cool")
    type(cool)  # kaic.compatibility.cooler.CoolerHic

    juicer = kaic.load("test_juicer.hic")
    type(juicer)  # kaic.compatibility.juicer.JuicerHic

    fragments = kaic.load("hg19_chr18_19_re_fragments.bed")
    type(fragments)  # genomic_regions.regions.Bed

    bam = kaic.load("test.bam")
    type(bam)  # pysam.libcalignmentfile.AlignmentFile

    ins = kaic.load("architecture/domains/kaic_example_100kb.insulation")
    type(ins)  # kaic.architecture.domains.InsulationScores

    # and many other data types

The next section will discuss :ref:`common_interfaces` that make working with genomic data
in general and Kai-C objects specifically straightforward and simple.
