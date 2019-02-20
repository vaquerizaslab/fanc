.. _kaic-comparisons:

############################
Matrix and score comparisons
############################

Kai-C provides a central utility named ``kaic compare`` to compare Hi-C matrices
and measures derived from them between two conditions, such as different cell types,
treatments, etc.

.. argparse::
   :module: kaic.commands.kaic_commands
   :func: compare_parser
   :prog: kaic compare
   :nodescription:

*******
Example
*******
