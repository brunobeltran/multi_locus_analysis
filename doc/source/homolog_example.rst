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

.. figure:: _static/homologs/Fig1.svg
    :alt: Figure showing (a) prophase stages, (b) colocalization progression, and (c) tagged locus genomic locations

    A schematic of the relative timing of the chromosome events of meiosis in SK1
    strains of budding yeast.


.. {padmore1991,weiner1994,
.. cha2000,tesse2003,brar2009,borner2004,peoples2002}. (a) Chromosomes in
.. pre-meiotic cells arrested in G0 are in the Rabl configuration with centromeres
.. tethered to the nuclear periphery. (b) Early to mid prophase is marked by DSB
.. formation and the initiation of synapsis. (c) Late prophase is marked by the
.. end-to-end alignment of homologs by the synaptonemal complex. (d) Fraction cells
.. over time that demonstrate colocalization of the \textit{URA3} locus and
.. completion of meiosis I (MI). The x-axis measures the time $T_i$ ($i$ hours)
.. after induction of sporulation that the cells in question were prepared for
.. imaging. Pre-meiotic colocalization is lost during DNA replication and is
.. restored during meiotic prophase, culminating in the full-length alignment of
.. homologs joined by the synaptonemal complex (SC). Soon afterwards, cells begin
.. to complete meiosis I (MI). (e) The relative positions along the chromosome of
.. our tagged loci are shown. These loci were chosen to probe the dependence of
.. colocalization on centromere proximity.


.. .. raw:: html

..     <object data="_static/homologs/Fig1.svg" type="image/svg+xml"></object>

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

TODO: copy in justifications for both 15nm and 50nm.

Determining diffusivity
^^^^^^^^^^^^^^^^^^^^^^^

TODO: fit initial power-law slope of WT URA t3 to analytical theory to determine
D.

.. include:: determining-diffusivity.rst

Determining nuclear radius
^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO: copy in from Jupyter notebook here.

Experimental waiting times
--------------------------

.. include:: experimental-waiting-times.rst
