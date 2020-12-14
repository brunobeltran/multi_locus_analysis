Welcome to multi_locus_analysis's documentation!
================================================

MLA is a collection of utilities for analyzing multi-locus particle-tracking
data.

Please see one of our various tutorials to get started analyzing diffusing
particle trajectories, first passage times, MSDs, correlations, etc. and
comparing these to established analytical diffusion theory.

.. TODO: For a full tutorial on how to reproduce the figures from (Newman, Beltran,
   and Calhoon et. at.,  *in preparation*), see :ref:`Yeast Meiosis <burgess>`.

.. TODO: For a full tutorial on how to reproduce the figures from (Beltran,
    MacPherson, Spakowitz, *in preparation*), see :ref:`the Waiting Times
    tutorial <finite_window>

To compute MSDs, velocity correlations, and other time-dependent moments of the
trajectory, see :ref:`Trajectory Statistics <stats>`. For
computing waiting times—given a trajectory that can take multiple different
states (such as "paired" and "unpaired")—see :ref:`Waiting Times
<finite_window>`. For examples of comparing
the analytical theory of (polymer) diffusions to the data, see :ref:`Analytical
Theory <theory>`. Custom plotting code used in the :ref:`theory <theory>`
section is documented in our :ref:`plotting guide <plots>`, but for actual
example plots, see the other tutorials.

For a full walkthrough of how this code was used to analyze tens of thousands
of homologous locus trajectories across dozens of experimental conditions, see
the our :ref:`Analyzing homologous loci example <homolog_example>`.

.. toctree::
    :maxdepth: 1

    Trajectory Statistics <stats>
    Waiting Times <finite-window>
    Analytical Theory <theory>
    Plotting <plots>
    API reference <api>
    Homologs Example <homolog_example>

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
