.. _homolog_example:

.. currentmodule:: multi_locus_analysis.examples.burgess

Analyzing homologous loci trajectories
======================================

In collaboration with the `Burgess lab
<https://smburgess.faculty.ucdavis.edu/>`_ at UC Davis, we have compiled an
example data set that involves measurements of pairs of homologous loci
diffusing throughout the nucleii of *S. cerevisiae* cells at different stages of
meiosis I.

In the following sections, we lay out the process of

1. extracting the nuclear radius, diffusivity, and other parameters from this
   kind of data.
2. analyzing the time it takes the loci to colocalize, and how long they stay
   colocalized.
3. comparing these results to analytical theory from the `wlcsim
   <wlcsim.rtfd.io>`_ module.

In order for our plots to have the same styling as in our paper, we first set
our matplotlib style:

.. plot::
    :nofigs:
    :context: close-figs

    >>> from multi_locus_analysis.examples.burgess.styles import *
    >>> use_pnas_style()


Description of the data
-----------------------

Our study used yeast strains containing chromosomes carrying FROS tags,
comprised of chromosomally-integrated *tet* operator arrays of 112 repeats bound
by fluorescent TetR-GFP protein.  Operators were inserted at
either the *URA3* locus---which is on the short arm of chr.~V near the
centromere, or the *LYS2* locus---which is in the center of the long arm
of chr. |nbsp| II.

On the population level, these loci are known to start off largely colocalized at
the end of G0, then separate at the start of meiosis, before becoming
colocalized again during prophase |nbsp| I:

.. figure:: _static/homologs/Fig1.svg
    :alt: Figure showing (a) prophase stages, (b) colocalization progression, and (c) tagged locus genomic locations

    A schematic of the relative timing of the chromosome events of meiosis in
    SK1 strains of budding yeast,
    :cite:`padmore1991,weiner1994,cha2000,tesse2003,brar2009,borner2004,peoples2002`.
    (a) Chromosomes in pre-meiotic cells arrested in G0 are in the Rabl
    configuration with centromeres tethered to the nuclear periphery. (b) Early
    to mid prophase is marked by DSB formation and the initiation of synapsis.
    (c) Late prophase is marked by the end-to-end alignment of homologs by the
    synaptonemal complex. (d) Fraction cells over time that demonstrate
    colocalization of the *URA3* locus and completion of meiosis I (MI). The
    x-axis measures the time :math:`T_i` (:math:`i` hours) after induction of
    sporulation that the cells in question were prepared for imaging.
    Pre-meiotic colocalization is lost during DNA replication and is restored
    during meiotic prophase, culminating in the full-length alignment of
    homologs joined by the synaptonemal complex (SC). Soon afterwards, cells
    begin to complete meiosis I (MI). (e) The relative positions along the
    chromosome of our tagged loci are shown. These loci were chosen to probe the
    dependence of colocalization on centromere proximity.

We wish to uncover what the forces are pulling these loci together. In order to
do so, we must first establish a baseline model for what the diffusion of these
loci would look like in the absense of any force.

Confined Rouse polymer with linkages
------------------------------------

The Rouse model is the most basic possible model for a diffusing polymer, and
applies even to semiflexible polymers such as DNA at long enough length and
time scales (significantly longer than the persistence length).

TODO: describe the model more in detail, especially linkages, and relevant
parameters.

Running the Simulation
^^^^^^^^^^^^^^^^^^^^^^

The code used to run our simulations (including the final parameters chosen) can
be reproduced by running:

.. code:: bash

    >>> python -m multi_locus_analysis.examples.burgess.simulation

In short, this calls `~.simulation.run_homolog_param_scan`, which parallelizes
the running of `~.simulation.run_interior_sim` for various values of the mean
linkage density :math:`mu`. These calls are wrapped in the
`~.simulation.save_interior_sim` function to ensure that each individual
simulation is saved to its own sub-directory (making this parallelization
thread-safe, and even multi-node safe). The output of this command will be a
directory structure that looks like:

.. code:: bash

    >>> tree homolog-sim | head
    homolog-sim
    ├── homolog-sim.tower13.1
    │   ├── all_beads.csv
    │   └── params.pkl
    ├── homolog-sim.tower13.10
    │   ├── all_beads.csv
    │   └── params.pkl
    .
    .
    .

where each folder is named ``"{base_name}/{base_name}.{hostname}.{i}"``, and the
code guarantees a unique folder name ``{base_name}.{hostname}.{i}`` per
simulation. The ``params.pkl`` file holds the parameters that were passed to
`wlcsim.bd.homolog.rouse`. Each ``all_beads.csv`` file contains a Pandas
dataframe, a representative example might look like:

.. code:: python

    >>> all_beads.head()
        t  bead         X1   ...         Z2  is_loop  is_tether    FP
    0  0.0     0   0.000000  ...  23.677611        0          0  0.04
    1  0.0     1  -0.357303  ...  27.119501        0          0  0.04
    2  0.0     2 -14.284399  ...  12.135141        0          0  0.04
    3  0.0     3 -29.879801  ...  16.638270        0          0  0.04
    4  0.0     4 -71.418055  ... -61.567177        0          0  0.04

where ``'t'`` is the simulation time, ``'bead'`` is the bead index, ``'X{i}'``
is the :math:`x`-coordinate of polymer :math:`i` (since we're simulating
a homologous pair, :math:`i\in{1, 2}`), ``'is_loop'`` is true if that bead is a
linkage, ``'is_tether'`` is true if that bead is tethered to the nuclear
envelope, and ``'FP'`` is the fraction of beads that are linkages on average, in
this simulation run.

Parameterization of Rouse model
-------------------------------

Importing Raw Data
^^^^^^^^^^^^^^^^^^

First, we import the raw data from the examples module. For a complete
description of each column, see the documentation in
`.multi_locus_analysis.examples.burgess`.

.. plot::
    :nofigs:
    :context:

    >>> from multi_locus_analysis.examples import burgess
    >>> burgess.df[['X', 'Y', 'Z']].head()

    locus genotype exp.rep meiosis cell frame spot
    HET5  WT       2       t0      1    1     1     2.13328  3.19992  7.75
                                        2     1     2.53327  1.99995  7.75
                                        3     1     2.66660  2.39994  7.50
                                        4     1     2.79993  2.53327  7.50
                                        5     1     2.93326  2.39994  7.50

Justification for the Kuhn length
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO: copy in justifications for both 15nm and 50nm.

Determining nuclear radius
^^^^^^^^^^^^^^^^^^^^^^^^^^

The budding yeast nucleus is typically described as being a sphere with radius
of around :math:`1\mu{}m`, :cite:`phillips2012`. However, the nucleur size can
change drastically between strains :cite:`berger2008`, and between growth
conditions :cite:`jorgensen2007`. In addition, the nucleus grows two-fold, along with
the cell, between the G1 and S phases :cite:`jorgensen2007`.

In addition, the nucleus is almost never percisely spherical, however, recent
measurements have shown that the nucleus is typically approximately spherical
throughout the cell cycle (with a mild elongation at the start of mitosis) once
spherical abberations from the typical imaging setup are corrected for
:cite:`wang2016`. To our knowledge, no reliable measurements of this same type
exist for cells entering meiosis, so we will instead use the convex hull of the
volume explored by our tagged loci to estimate the nuclear radius (more
precisely, to set a lower bound on this radius).


.. plot::
    :context: close-figs

    >>> from multi_locus_analysis.stats import convex_hull
    >>> fig, ax = plt.subplots(constrained_layout=True, figsize=(col_width, golden_ratio*col_width))
    >>> chull_volume = burgess.df \
    >>>     .groupby(burgess.cell_cols) \
    >>>     .apply(convex_hull, xcol='X', ycol='Y', zcol='Z', volume=True)
    >>> volume_to_r = lambda V: np.power(3/4/np.pi*V, 1/3)
    >>> chull_volume.loc['HET5', 'WT', :, 't5'] \
    >>>             .apply(volume_to_r) \
    >>>             .hist(ax=ax)
    >>> ax.set_xlabel(r'Nuclear Radius ($\mu{}m$)')
    >>> ax.set_ylabel('Count')


Phenomenological MSCD Correction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Unlike in our analytical theory, the

Determining diffusivity
^^^^^^^^^^^^^^^^^^^^^^^

TODO: fit initial power-law slope of WT URA t3 to analytical theory to determine
D.

Example "cells"
---------------

Given all of these parameters, we can now use our model to visualize how
heterogeneity between individual cells affects the single-cell MSCD curves, and
compare that to what we see in the experiments. Consider the following example
cells:

.. plot::
    :context: close-figs

    >>> from multi_locus_analysis.examples.burgess import plotting as mplt
    >>> cells = [homolog.generate_poisson_homologs(4, burgess.chrv_size_bp)
    >>>          for i in range(5)]
    >>> mplt.draw_cells(cells)


.. plot::
    :context: close-figs

    >>> t = np.arange(30, 1501, 30)
    >>> plateaus = [homolog.mscd_plateau(links) for links in cells]
    >>> i = np.argsort(plateaus)
    >>> for i, cell in enumerate(cells):
    >>>     pass
    #TODO: change everything (draw_cells) to um


.. include:: determining-diffusivity.rst

Experimental waiting times
--------------------------

.. include:: experimental-waiting-times.rst


.. bibliography:: homologs.bib
