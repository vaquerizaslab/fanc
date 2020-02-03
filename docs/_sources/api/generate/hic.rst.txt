.. _api_hic:

##############################################
Generate, bin, filter, and correct Hic objects
##############################################

After we have obtained a filtered :class:`~fanc.pairs.ReadPairs` object, we can easily
convert it into a :class:`~fanc.hic.Hic` matrix.

.. literalinclude:: code/generate_example_code.py
    :language: python
    :start-after: start snippet hic convert
    :end-before: end snippet hic convert

Note that this only uses valid pairs for Hi-C object creation and creates **fragment-level**
Hi-C matrix, using the same fragment definitions as in the original :class:`~fanc.pairs.ReadPairs`
object.

.. note::

    If you have multiple Hi-C objects for the same sample, for example from different replicates or
    sequencing runs, we recommend merging them at this stage:

    .. code:: python

       from fanc.hic import Hic
       merged_hic = Hic.merge([hic_rep1, hic_rep2]])

The fragment-level Hi-C objects can very easily be binned:

.. literalinclude:: code/generate_example_code.py
    :language: python
    :start-after: start snippet hic bin
    :end-before: end snippet hic bin

Just like :class:`~fanc.pairs.ReadPairs`, :class:`~fanc.hic.Hic` can also be filtered. Currently,
the two available filters are :class:`~fanc.hic.LowCoverageFilter` and
:class:`~fanc.edge.DiagonalFilter`:

.. literalinclude:: code/generate_example_code.py
    :language: python
    :start-after: start snippet hic filter
    :end-before: end snippet hic filter

- :class:`~fanc.hic.LowCoverageFilter` removes all Hi-C "edges" (pixels) in regions that have
  low coverage. The cutoff can either be expressed in absolute numbers of pairs (``cutoff``) or
  as a fraction of the median coverage of all regions (``rel_cutoff``). We highly recommend
  filtering for low coverage before matrix balancing
- :class:`~fanc.edge.DiagonalFilter` removes all edges at the diagonal. Some researcher have
  achieved better normalisation results by setting the matrix diagonal to 0 before normalisation.

Finally, we can normalise the matrix using matrix balancing. You have the choice of either
Knight-Ruiz matrix balancing (KR, :class:`~fanc.hic.kr_balancing`) or iterative correction (ICE,
:class:`~fanc.hic.ice_balancing`). KR balancing is typically faster, but also consumes a lot
more memory, especially for high matrix resolutions. ICE is slow, and might not always converge
to a solution, but is much more memory efficient in its implementation.

.. literalinclude:: code/generate_example_code.py
    :language: python
    :start-after: start snippet hic balance
    :end-before: end snippet hic balance

With KR balancing, you can choose to restore the original valid pairs count using
``restore_coverage=True``. If this is set to ``False``, edge weights represent ligation
probabilities rather than counts. If you set ``whole_matrix=True``, the balancing will be done
on the whole genome matrix. The default is to do it on a per-chromosome basis.

The bias vector will be stored in the matrix, and is automatically applied when querying the object.