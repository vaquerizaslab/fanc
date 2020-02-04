.. _fanc-pca:

############
PCA analysis
############

When working with multiple Hi-C libraries, it is often useful to assess the variability
between replicates and samples from different conditions. This can provide valuable
information about potential experimental biases and whether samples from different
replicates can be safely merged.

.. argparse::
   :module: fanc.commands.fanc_commands
   :func: pca_parser
   :prog: fanc pca
   :nodescription:
   :nodefault:

*******
Example
*******
As an example, we ran a PCA analysis on the 1mb resolution mESC Hi-C matrices from our
`Low-C paper <https://doi.org/10.1038/s41467-018-06961-0>`_ using different restriction
enzymes (MboI and HindIII), as well as different input cell numbers.

.. code::

    fanc pca -n "HindIII 100k" "HindIII 5M" "MboI 100k" "MboI 1M" "MboI 50k" \
             -Z -s 100000 -r chr19 -p architecture/pca/lowc.pca.png \
             architecture/other-hic/lowc_hindiii_100k_1mb.hic \
             architecture/other-hic/lowc_hindiii_5M_1mb.hic \
             architecture/other-hic/lowc_mboi_100k_1mb.hic \
             architecture/other-hic/lowc_mboi_1M_1mb.hic \
             architecture/other-hic/lowc_mboi_50k_1mb.hic \
             architecture/pca/lowc.pca

The result looks like this, where you can clearly see the division between the different
restriction enzymes:

.. image:: images/lowc.pca.png

*******
Filters
*******

By default, PCA is run on the whole genome. In the example above, we have restricted the
analysis to chromosome 19 using the ``-r chr19`` argument. ``-Z`` instructs
``fanc pca`` to use only non-zero matrix entries for the PCA - this can help mitigate
the effect of very weak contacts on the variability.

In the example, we are limiting the number of contacts used for the PCA to 100,000
using ``-s 100000``. The type of contacts in the sample are chose according to the
``--strategy`` option. By default, the contacts with the largest variability across
the genome are chosen, but you can also use ``--strategy fold-change`` to choose
contacts with the largest fold-change or ``--strategy passthrough`` to make no
prior selection of contacts.

If you only want to include contacts up to (or above a) a certain distance, you can
specify that distance using the ``--max-distance`` (or ``min-distance``) option.

Finally, ``fanc pca`` offers two filters designed to remove uninformative contacts
before PCA. The first, ``--expected-filter <f>``, removes all contacts in which all
samples have a signal below *f* x *E(d)*, where E is the expected value function depending
on the distance d. The second, ``--background-filter <f>`` removes contacts in which
all samples have an *f* -fold enrichment over the expected inter-chromosomal contacts.
