.. _theory:

.. currentmodule:: multi_locus_analysis

Analytical theory
=================

Here, we show some best practices for comparing MSDs, displacement
distributions, and velocity correlations to the theory of diffusing molecules.

The literature contains many good reviews on these topics, so we cover
only some often overlooked pitfalls.

First, we discuss the :ref:`effects of confinement <analytical-msds>` on
MSDs, especially on extracting diffusivity and the subdiffusivity
coefficient (i.e. Hurst
index).

Then, we discuss the effects of heterogeneity between cells on the :ref:`gaussianity
of the displacement distribution <analytical-disps-hist>` and what that looks like for real data.

Finally, we discuss how to use the velocity cross correlation of two
particles to :ref:`extract information <analytical-velocity-correlation>` about whether those two particles are connected
in a drift-independent manner.

.. _analytical-msds:

MSDs
^^^^

The crowded cellular environment is often subdiffusive, and across bacteria, in
some yeast and some higher mammalian cells, it has been shown that this
diffusion is well described by a fractional brownian motion (for example, see
`(Weber et. al., PRL, 2010)
<https://web.stanford.edu/~ajspakow/publications/WST-ecoli.pdf>`_).

A particle undergoing fractional brownian motion has a probability distribution
given by the fractional diffusion equation

.. math::

    \partial_t^\alpha p(x, t) = D*\partial_x^2 p(x,t)

where we use the usual Caputo definition of the fractional derivative
:math:`\partial_t^\alpha`. This corresponds to a power-law memory kernel,
:math:`K(t-t_0) = \frac{(2-\alpha)(1-\alpha)}{|t - t_0|^\alpha}`, where

.. math::

    \xi \int_0^t K(t - t_0) \frac{d X(t_0)}{dt_0} dt_0 = \ldots{}

replaces the constant viscosity of the regular Langevin equation

.. math::

    \xi \frac{d X(t)}{dt} = \ldots{}

Intuitively, this means that there is a viscoelastic force that tries to return
the particle to it's original location (by opposing the historical velocity of
the particle). The strength of this force falls off as :math:`t^{-\alpha}`.

.. note::

    The continuous time random walk (CTRW) can be modeled---more specifically,
    the random waiting times of the CTRW model can be included, regardless of
    whether a power-law viscoelastic memory kernel is assumed---by replacing the
    2nd spatial derivative in our diffusion equation with a fraction derivative.

    In what follows, we consider the diffusion to be well described by a simple
    fractional brownian motion, with no random waiting times between collisions.

The viscoelastic memory kernel is parameterized by the single parameter
:math:`\alpha`, which describes how much the system wishes to return
(elastically) to a previous configuration. More precisely, the MSD of a particle
diffusing according to the equations above can be written as

.. math::

    \frac{3k_BT}{\xi}\frac{\sin(\alpha\pi)}{\pi(1-\alpha/2)(1-\alpha)\alpha} t^\alpha

In other words, in log-log space, the MSD will be a line with slope
:math:`\alpha`.

If we instead consider the diffusion of a polymer instead of a point, we get a
new equation of motion.

TODO: add details on polymer MSD here.

This means that for a point on the polymer, the observed MSD in log-log space
will be a line with slope :math:`\alpha/2` for short time scales. This persists
until the terminal relaxation time :math:`t_R = [N^2b^2\xi/(k_BT)]^{1/\alpha}`
of the polymer, which is the time scale at which the polymer as a whole can be
said to be diffusing (in other words, the time scale at which the diffusion of
the center of mass of the polymer begins to dominate over the internal
fluctuations of the polymer). At times scales longer than the terminal
relaxation time, the MSD will be a line with slope :math:`\alpha`, and the
polymer as a whole will diffuse as if experiencing an effective viscosity of
:math:`N\xi`.

TODO: add plot of rouse_mid_msd showing two regimes of MSD.

So in experimental MSDs, where the chromosomal loci are confined to a given
spatial region (whether that be a chromosome territory or the nucleus as a
whole), we expect to see the MSD curve cut off at a different location depending
on how the confinement time scale compares to the terminal relaxation time of
the chromosome (or chromosomes, if they are linked together).

.. _analytical-disps-hist:

Displacement Distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO: document gaussianity vs laplace distributions of dispalcements here.

.. _analytical-velocity-correlation:

Velocity Cross-Correlation
^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to measure whether (and how) stress is communicated between the loci
being measured, we can look at their velocity cross correlation functions.

As shown in `(Lampo & Spakowitz,
Biophysical Journal, 2016)
<https://www.sciencedirect.com/science/article/pii/S0006349515047049>`_, the
velocity cross-correlation of two loci on a Rouse polymer can be computed
exactly.




