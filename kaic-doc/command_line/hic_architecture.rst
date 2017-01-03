.. _kaic-analysis:

#################
Hi-C architecture
#################

Kai-C provides a lot of analysis methods that have become standard in the field. In this section, we will use the
Hi-C matrix generated from the example data in the ``kaic/test/examples/`` folder on our `GitHub page
<http://www.github.com/vaquerizaslab/kaic>`_. If you want to follow the examples on this page, please first generate
Hi-C matrices using one of these tutorials: :ref:`example-kaic-auto` or :ref:`modular-analysis`.


***************************************
Methods to find TADs and TAD boundaries
***************************************

Topologically associating domains were a central finding of Hi-C experiments in higher eukaryotes. Kai-C contains
two of the most popular methods for identifying TADs: the insulation index and the directionality index.


Insulation index
================

