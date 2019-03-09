.. _plots:

.. currentmodule:: multi_locus_analysis

Plotting utilities
==================

Helper functions used to generate the examples in the :ref:`Analytical Theory
<theory>` tutorial are contained in the
:mod:`multi_locus_analysis.plotting` module. Please see the :ref:`API docs
<api>` for details.

Historically, this package's author has been able to make most plots related
to multi locus analysis with a couple of short lines of :mod:`matplotlib` code.

Histograms in Python are not the most elegant, but
:mod:`multi_locus_analysis.finite_window` has several methods for
representing histograms fairly as smooth (or step) functions that can simply be
plotted.

Documentation on how to get around several common issues encountered when plotting (in particular, using
colors, colorbars, normalization and facets, can be found in the `docs <bruno_util.rtfd.io>`_ for
the external module :mod:`bruno_util.plotting`.

All generally-applicable plotting routines developed within this codebase
have been moved to :mod:`bruno_util`, along with their documentation.

Displacement Distributions
--------------------------

The displacement distribution, :math:`f_{V_i^\delta}(x)` describes how far
(:math:`x`) a particle has moved after a given time :math:`\delta`.

This is a distribution for every value of :math:`\delta`. Plotting one curve for
every value of :math:`\delta` is unwieldy:

.. note::

    There will be a plot here, but not made yet.

.. .. image::

..     figures/all_disps_hist.png

We provide :func:`multi_locus_analysis.plotting.make_all_disps_hist`
for quickly plotting subsets of the histograms, comparing to analytical
theory, etc. See :ref:`analytical-disps-hist` for mathematical details.

.. note::

    There will be a plot here, but not made yet.

.. .. image::

..     figures/all_disps_hist_subsample.png

Velocity Correlations
^^^^^^^^^^^^^^^^^^^^^

Similarly to the displacements distribution, the velocity correlation
function :math:`<V_{ij}^{\delta}(t)\cdot{}V_{ij}^{\delta}(t)>` is a function for
each value of :math:`\delta`, and so is cumbersome to plot.

We provide :func:`multi_locus_analysis.plotting.cvv_plot_sized` to aid in
subsampling, weighting the plotted curves visually based on their statistical
certainty, and other useful tricks, like comparing to the analytical theory of
polymer diffusion as described in our :ref:`analytical theory tutorial
<analytical-velocity-correlation>`.
