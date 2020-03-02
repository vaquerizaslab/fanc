.. _api_map:

===================
Mapping FASTQ files
===================

To map FASTQ files using the API, we will use the :mod:`fanc.map` module.

*******
Mappers
*******

First, we need to decide which mapper to use, and under what circumstances an unaligned
read should be truncated and re-mapped to the genome. The :class:`~fanc.map.SimpleBowtie2Mapper`,
for example, only attempts to align a read once and does no iterative mapping. The
:class:`~fanc.map.Bowtie2Mapper` resubmits unaligned reads and reads with a low score,
There are equivalent mappers for BWA named
:class:`~fanc.map.SimpleBWAMapper` and :class:`~fanc.map.BWAMapper`, but here we
choose the :class:`~fanc.map.Bowtie2Mapper`. It requires only the path of the corresponding
``bowtie2`` index:

.. literalinclude:: code/generate_example_code.py
    :language: python
    :start-after: start snippet mapper
    :end-before: end snippet mapper

The ``threads`` parameter controls how many threads are given to each ``bowtie2-align``
process. By using a `.bam` ending, the output is converted to a BAM file at the end of
mapping automatically.


*****************
Iterative mapping
*****************

Now we can use :func:`fanc.map.iterative_mapping` to start the actual mapping process:

.. literalinclude:: code/generate_example_code.py
    :language: python
    :start-after: start snippet iterative mapping
    :end-before: end snippet iterative mapping

Note that we are calling iterative mapping twice, independently for each FASTQ file, as
appropriate for a Hi-C experiment. ``min_size`` determines the minimum size of a truncated
read after which it will be discarded, while ``step_size`` determines the truncation amount.

.. note::

   When downloading FASTQ files from SRA using SRAtools, e.g. with `fastq-dump`, do not
   use the ``-I / --readids`` option, which appends ``.1`` or ``.2`` to the read name. This
   interferes with the sorting and read pairing step in FAN-C. **Read names of the two mates
   must be identical**.

By providing a ``restriction_enzyme`` name, we enable ligation junction splitting. Each read
will be scanned for a predicted ligation junction of the provided restriction enzyme and if
one is encountered, it will be split at the junction before alignment. This can greatly increase
alignment rates, especially for longer reads.

***************
SAM/BAM sorting
***************

After the mapping is complete, we can sort the files by read name for further processing.

.. literalinclude:: code/generate_example_code.py
    :language: python
    :start-after: start snippet sort sam
    :end-before: end snippet sort sam

The above command replaces the original file with the sorted version. You can use the
``output_file`` parameter to output to a different file, if you prefer to keep the unsorted
version.
