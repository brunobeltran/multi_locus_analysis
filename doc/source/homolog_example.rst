.. _homolog_example:

.. currentmodule:: multi_locus_analysis.examples.burgess

Analysing homologous loci trajectories
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
3. comparing these results to analytical theory from the :ref:`wlcsim` module.

Description of the data
-----------------------

Our study used yeast strains containing chromosomes carrying FROS tags,
comprised of chromosomally-integrated *tet* operator arrays of 112 repeats bound
by fluorescent TetR-GFP protein.  Operators were inserted at
either the *URA3* locus---which is on the short arm of chr.~V near the
centromere, or the *LYS2* locus---which is in the center of the long arm
of chr. |nbsp| II.

On the population leve, these loci are known to start off largely colocalized at
the end of G0, then separate at the start of meiosis, before becoming
colocalized again during prophase |nbsp| I:

.. raw:: html

    <object data="/homologs/Fig1.svg" type="image/svg+xml"></object>

We wish to uncover what the forces are pulling these loci together. In order to
do so, we must first establish a baseline model for what the diffusion of these
loci would look like in the absense of any force.


Parameterization of Rouse model
-------------------------------

The Rouse model is the most basic possible model for a diffusing polymer, and
applies even to semiflexible polymers such as DNA at long enough length and
time scales (significantly longer than the persistence length).

Justification for the Kuhn length
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Determining diffusivity
^^^^^^^^^^^^^^^^^^^^^^^

.. include:: determining-diffusivity.rst


Experimental waiting times
--------------------------

.. include:: experimental-waiting-times.rst




