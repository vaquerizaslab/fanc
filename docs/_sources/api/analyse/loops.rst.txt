.. _api_loops:


############
Loop calling
############

To follow this tutorial, download the FAN-C example data, for example through our
`Keeper library <https://keeper.mpdl.mpg.de/d/147906745b634c779ed3/>`_, and set up
your Python session like this:

.. literalinclude:: code/loops_example_code.py
    :language: python
    :start-after: start snippet loops setup
    :end-before: end snippet loops setup

Also have a look at the command line documentation at :ref:`fanc-loops` for the
command-line approach to loop calling in FAN-C.

Loop calling in FAN-C is a 3-step process:

1. Annotate pixels with local enrichment, probability, and mappability information
2. Filter pixels based on suitable cutoffs to identify pixels that are part of a loop
3. Merge nearby pixels into a single loop call


***************
Annotate pixels
***************

Start by annotating each pixel with loop information using the
:class:`~fanc.peaks.RaoPeakCaller` and :func:`~fanc.peaks.RaoPeakCaller.call_peaks` function:

.. literalinclude:: code/loops_example_code.py
    :language: python
    :start-after: start snippet loops annotate
    :end-before: end snippet loops annotate

This creates a :class:`~fanc.peaks.RaoPeakInfo` object, which contains a bunch of information
on each pixel. Here are the most important ones:

- ``weight`` and ``uncorrected`` are the normalised and raw contact strength, respectively
- ``e_d``, ``e_h``, ``e_v``, and ``e_ll`` are the expected values calculated from the
  donut, horizontal, vertical, and lower left neighborhoods, repectively, of the pixel
- ``oe_d``, ``oe_h``, ``oe_v``, ``oe_ll`` are the observed/expected (O/E) values
  calculated from the donut, horizontal, vertical, and lower left neighborhoods, repectively,
  of the pixel
- ``fdr_d``, ``fdr_h``, ``fdr_v``, ``fdr_ll`` are the false discovery rates (FDRs)
  calculated from the donut, horizontal, vertical, and lower left neighborhoods, repectively,
  of the pixel
- ``mappability_d``, ``mappability_h``, ``mappability_v``, ``mappability_ll`` are the mappabilities
  of the donut, horizontal, vertical, and lower left neighborhoods, repectively,
  of the pixel. This is a value between 0 and 1 stating how many pixels in each neighborhood
  have valid weights

By default, most intra-chromosomal pixels will receive this information, but you have some
influence over the computation with additional parameters. For example, loops tend to occur
away from the diagonal, so you may want to exclude pixels near the diagonal using ``min_locus_dist``,
whihc is expressed in bins. You can also change the values for ``p`` and ``w_init`` (see original
publication for details) - otherwise sensible defaults are chosen. ``min_mappable_fraction``
controls which pixels are excluded if their local neighbourhoods are below a certain mappability
fraction, which is 0.7 by default, i.e. 70% of each neighbourhood must be valid.

Finally, loops calling is very resource intensive and needs to be heavily parallelised. The number
of parallel processes is controlled with ``n_processes``. By default, this runs loop calling
locally. However, we highly recommend setting ``cluster=True`` if you are in a cluster environment
that supports DRMAA (e.g. SGE, Slurm). You may need to set the ``DRMAA_LIBRARY_PATH`` environment
variable to the location of your ``libdrmaa.so`` in your shell in order for this to work.
:func:`~fanc.peaks.RaoPeakCaller.call_peaks` must then be called from a head node in order to be
able to submit jobs to the cluster. Read more about the interface with DRMAA in the
`gridmap package <https://github.com/pygridtools/gridmap>`_.

To access the new pixel information, you can use the :func:`~fanc.peaks.RaoPeakInfo.peaks` iterator:

.. literalinclude:: code/loops_example_code.py
    :language: python
    :start-after: start snippet loops iter
    :end-before: end snippet loops iter

We can also plot it like a regular matrix:

.. literalinclude:: code/loops_example_code.py
    :language: python
    :start-after: start snippet loops matrix
    :end-before: end snippet loops matrix

.. image:: images/loops_annotate.png

Note how we are controlling which attribute is plotted with ``weight_field``.

These kinds of plots can be extremely useful to choose appropriate filtering thresholds for
each of these attributes, which we will see in the next section.

*************
Filter pixels
*************

Based on the annotation, we can try to filter out pixels that are not part of loops. This can
be done with filters, specifically

- :class:`~fanc.peaks.EnrichmentPeakFilter` filters pixels below a minimum O/E in specific
  neighborhoods
- :class:`~fanc.peaks.FdrPeakFilter` filters pixels with an FDR higher than specified in each
  neighborhood
- :class:`~fanc.peaks.MappabilityPeakFilter` filters pixels with a minimum mappable fraction
  below the specified cutoffs
- :class:`~fanc.peaks.ObservedPeakFilter` filters pixels that do not have a minimum number of
  raw reads mapping to them
- :class:`~fanc.peaks.DistancePeakFilter` filters pixels that are closer than a certain number
  of bins from the diagonal

Each of these is instantiated and then added to the object using
:func:`~fanc.peaks.RaoPeakInfo.filter`. When applying multiple filters it is recommended to
first queue them using ``queue=True`` and to run all of them together with
:func:`~fanc.peaks.RaoPeakInfo.run_queued_filters`:

.. literalinclude:: code/loops_example_code.py
    :language: python
    :start-after: start snippet loops filter
    :end-before: end snippet loops filter

.. image:: images/loops_remaining.png

We could also tweak the above cutoffs to further remove noisy pixels. The next step is to merge
pixels belonging to the same loop into loop calls.

************
Merge pixels
************

To merge the remaining pixels into loops use

.. literalinclude:: code/loops_example_code.py
    :language: python
    :start-after: start snippet loops merge
    :end-before: end snippet loops merge

``merged_peaks`` now contains putative loop calls. However, there will still be a number of
false-positive loops in this object, which typically consist of a single enriched pixel, likely
due to noise. Real loops generally consist of multiple pixels. We can remove the "singlets"
using :class:`~fanc.peaks.RaoMergedPeakFilter`:

.. literalinclude:: code/loops_example_code.py
    :language: python
    :start-after: start snippet loops singlet
    :end-before: end snippet loops singlet

A cutoff here of ``-1`` remove all singlets. In Rao, Huntley et al. (2014) they only remove
singlets if their qvalue sum is larger than ``0.02``, which you can use instead of ``-1``
and particularly strong singlets will still be kept in the data.

For our example dataset this leaves only a single loop, which we can export with
:func:`~fanc.peaks.PeakInfo.to_bedpe`:

.. literalinclude:: code/loops_example_code.py
    :language: python
    :start-after: start snippet loops export
    :end-before: end snippet loops export
