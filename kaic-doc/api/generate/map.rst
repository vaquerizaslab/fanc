===================
Mapping FASTQ files
===================

To map FASTQ files using the API, we will use the :mod:`kaic.map` module.

*******
Mappers
*******

First, we need to decide which mapper to use, and under what circumstances an unaligned
read should be truncated and re-mapped to the genome. The :class:`~kaic.map.SimpleBowtie2Mapper`,
for example, only attempts to align a read once and does no iterative mapping. The
:class:`~kaic.map.Bowtie2Mapper` resubmits unaligned reads and reads with a low score,
There are equivalent mappers for BWA named
:class:`~kaic.map.SimpleBWAMapper` and :class:`~kaic.map.BWAMapper`, but here we
choose the :class:`~kaic.map.Bowtie2Mapper`. It requires only the path of the corresponding
``bowtie2`` index:

.. literalinclude:: code/generate_example_code.py
    :language: python
    :start-after: start snippet mapper
    :end-before: end snippet mapper

The ``threads`` parameter controls how many threads are given to each ``bowtie2-align``
process.

Now we can use :func:`kaic.map.iterative_mapping` to start the actual mapping process:

.. literalinclude:: code/generate_example_code.py
    :language: python
    :start-after: start snippet iterative mapping
    :end-before: end snippet iterative mapping

Note that we are calling iterative mapping twice, independently for each FASTQ file, as
appropriate for a Hi-C experiment. ``min_size`` determines the minimum size of a truncated
read after which it will be discarded, while ``step_size`` determines the truncation amount.

By providing a ``restriction_enzyme`` name, we enable ligation junction splitting. Each read
will be scanned for a predicted ligation junction of the provided restriction enzyme and if
one is encountered, it will be split at the junction before alignment. This can greatly increase
alignment rates, especially for longer reads.

After the mapping is complete, we can sort the files by read name for further processing.

.. literalinclude:: code/generate_example_code.py
    :language: python
    :start-after: start snippet sort sam
    :end-before: end snippet sort sam

The above command replaces the original file with the sorted version. You can use the
``output_file`` parameter to output to a different file