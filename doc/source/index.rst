Welcome to multi_locus_analysis's documentation!
================================================

MLA is a collection of utilities for analyzing multi-locus particle-tracking
data.

.. For a full tutorial on how to reproduce the figures from (Newman, Beltran,
   and Calhoon et. at.,  *in preparation*), see :ref:`Yeast Meiosis <burgess>`.

If you're not sure what you're looking for, try perusing :ref:`the plotting
guide <plots>` to get an idea of what this module can do!

To compute MSDs, velocity correlations, and other time-dependent moments of the
trajectory, see :ref:`Trajectory Statistics <stats>`. For
computing waiting times—given a trajectory that can take multiple different
states (such as "paired" and "unpaired")—see :ref:`Waiting Times
<finite_window>`. For examples of common plots (especially
for viewing the output of :mod:`multi_locus_analysis.stats`) and their
interpretations, see :ref:`Plotting <plots>`. For examples of comparing
the analytical theory of (polymer) diffusions to the data, see :ref:`Analytical
Theory <theory>`.

.. toctree::
    :maxdepth: 1

    Trajectory Statistics <stats>
    Waiting Times <finite_window>
    Plotting <plots>
    Analytical Theory <theory>
    API reference <api>

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
