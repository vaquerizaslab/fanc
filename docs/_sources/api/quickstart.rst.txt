.. _quickstart:

=============================
Quickstart with the FAN-C API
=============================

After you have installed FAN-C (see :ref:`fanc_installation`), you can import the module
from a Python console or script:

.. code::

    import fanc

The following tutorials will assume that you have loaded the :code:`fanc` module in this
manner.

****************
Loading datasets
****************

Any analysis typically begins with loading datasets into your workspace. FAN-C tries to make this
as simple as possible with the :func:`~fanc.__init__.load` function.
If you already have processed Hi-C files, either from the :code:`fanc` command line
application (see :ref:`fanc-auto` or :ref:`fanc-modular`), or from a compatible Hi-C application
(:code:`.cool` or :code:`.mcool` from `Cooler <https://github.com/mirnylab/cooler>`_ or :code:`.hic` from
`Juicer <https://github.com/aidenlab/juicer>`_), simply load them into your workspace
using :func:`~fanc.__init__.load` - no need to specify the type of file you are loading:

.. code::

     data = fanc.load("/path/to/file.hic")

When dealing with multi-resolution Hi-C files such as :code:`.mcool` from Cooler or :code:`.hic` from
Juicer, you can load a specific resolution using the :code:`@` notation:

.. code::

     data = fanc.load("/path/to/file.mcool@25000")

:func:`~fanc.__init__.load` is not limited to Hi-C files, but also works on any other file
produced with FAN-C, such as :class:`~fanc.pairs.RegionPairs` files, :class:`~fanc.regions.Genome`,
analysis results like :class:`~fanc.architecture.comparisons.FoldChangeMatrix` and generally most
other FAN-C files.

:func:`~fanc.__init__.load` even works on most of the common file formats for genomic
datasets, such as BED, GFF, BigWig, Tabix, BEDPE and more. Try it out on your dataset of choice and
chances are :func:`~fanc.__init__.load` can handle it. And if it does not, consider raising
an issue on `Github <https://github.com/vaquerizaslab/fanc/issues>`_ to ask for support.

Internally, :func:`~fanc.__init__.load` finds a suitable class for the type of data
in the supplied file, and opens the file using that class. For example, the result of

.. code::

    hic = fanc.Hic("output/hic/binned/fanc_example_1mb.hic")

is equivalent to

.. code::

    hic = fanc.load("output/hic/binned/fanc_example_1mb.hic")

with the big advantage that you don't need to worry about remembering class names or
their location within the FAN-C module hierarchy. In both cases, the type of the
returned object is :code:`fanc.hic.Hic`:

.. code::

    # check the type fo the object
    type(hic)  # fanc.hic.Hic

Here are a few more examples:

.. code::

    cool = fanc.load("test.cool")
    type(cool)  # fanc.compatibility.cooler.CoolerHic

    juicer = fanc.load("test_juicer.hic")
    type(juicer)  # fanc.compatibility.juicer.JuicerHic

    fragments = fanc.load("hg19_chr18_19_re_fragments.bed")
    type(fragments)  # genomic_regions.regions.Bed

    bam = fanc.load("test.bam")
    type(bam)  # pysam.libcalignmentfile.AlignmentFile

    ins = fanc.load("architecture/domains/fanc_example_100kb.insulation")
    type(ins)  # fanc.architecture.domains.InsulationScores

    # and many other data types

The next section will discuss :ref:`common_interfaces` that make working with genomic data
in general and FAN-C objects specifically straightforward and simple.

*******
Logging
*******

FAN-C uses the ``logging`` module. In a ``python`` session, use a statement like

.. literalinclude:: generate/code/generate_example_code.py
    :language: python
    :start-after: start snippet logging
    :end-before: end snippet logging

to enable basic console logging. Have a look the
`logging documentation <https://docs.python.org/3/library/logging.html>`_
for more information on log levels and handlers.