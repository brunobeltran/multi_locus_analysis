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

.. math::

    C_{vv}^{(\delta)}(t, n_1, n_2) = \left\langle \vec{V}^{(\delta)}(n_1, t)
    \cdot \vec{V}^{(\delta)}(n_2, 0) \right\rangle

after a cosine transform in :math:`n`, and expanding in the fraction diffusion
eigenbasis :math:`E_{\alpha,1}` (the Mittag-Leffler functions), we get

.. math::

    C_{vv}^{(\delta)}(t, n_1, n_2) = \frac{3 k_B T}{N \xi \Gamma(3 - \alpha)
    \Gamma(1 + \alpha)} \frac{2 C_v^{(\delta)}(t)}{C_v^{(\delta)}(0)}

.. math::

    + \sum_{p=1}^\infty \frac{3k_BT}{k_p}\phi_p(n_1)\phi_p(n_2) \left(
    -E_{\alpha,1}( -\frac{k_p|t+\delta|^\alpha}{N\xi\Gamma(3 - \alpha)} )
    \right.

.. math::

    \left. + 2E_{\alpha,1}( -\frac{k_p|t|^\alpha}{N\xi\Gamma(3 - \alpha)} )
    -E_{\alpha,1}( -\frac{k_p|t-\delta|^\alpha}{N\xi\Gamma(3 - \alpha)} )
    \right)

The first term captures the motion of the center of mass of the polymer, while
the terms in the sum capture the motion due to different "modes" of the polymer.

The dependence on the location along the polymer is captured by these Rouse
modes, as seen below

.. plot::
    :context: close-figs

    >>> from multi_locus_analysis import analytical
    >>> ps = np.arange(1000)
    >>> modes = analytical.rouse_mode(ps, 0.5, 1)
    >>> plt.scatter(ps, modes)
    >>> modes = analytical.rouse_mode(ps, 0.51, 1)
    >>> plt.scatter(ps, modes)
    >>> modes = analytical.rouse_mode(ps, 0.91, 1)
    >>> plt.scatter(ps, modes)
    >>> plt.xlabel('p')
    >>> plt.ylabel('$\phi_p(n/N)$')
    >>> plt.legend(['n = 0.5', 'n = 0.51', 'n = 0.91'])

The coefficients :math:`k_p = \frac{3\pi^2k_BT}{Nb^2} p^2` then determine how
strongly to weight each mode as a function of :math:`t` and :math:`\delta`.

Because :math:`E_{\alpha,1}(x)\to 0` as :math:`x\to-\infty`, and since
:math:`k_p \propto p^2`, the terms of the summation will be increasingly small.

Their absolute values are also incredibly small, which is strange, considering
that they must eventually influence the numerically large factor
:math:`2AC^{(\delta)}_v(t)/C^{(\delta)}_v(0)`.

We can instead express the velocity correlation in terms of Meijer G-functions

.. math::

    G^{m,n}_{p,q} \left( \left. \begin{matrix} a_1, \dots, a_n ; a_{n+1} \dots a_p \\ b_1, \dots, b_m ; b_{m+1} \dots b_q \end{matrix}\; \right| \; z ; r \right) = \frac{1}{2 \pi i} \int_L \frac{\prod_{j=1}^m \Gamma(b_j+s) \prod_{j=1}^n\Gamma(1-a_j-s)} {\prod_{j=n+1}^{p}\Gamma(a_j+s) \prod_{j=m+1}^q \Gamma(1-b_j-s)} z^{-s/r} ds

in short, we have (for the case :math:`\alpha = 1`)

.. math::
    C_{vv}^{(\delta)}(t, \Delta n) = \frac{3k_BT}{\delta^2\sqrt{\xi k}}

.. math::
    \times \left\{ |t - \delta|^{1/2} G_{1,2}^{2,0}\left[\left. \frac{|\Delta{}n|^2 \xi}{4k} |t - \delta|^{-1} \right|^{3/2}_{0,1/2}\right] \right.

.. math::
    + |t + \delta|^{1/2} G_{1,2}^{2,0}\left[ \frac{|\Delta{}n|^2 \xi}{4k} |t + \delta|^{-1} |^{3/2}_{0,1/2}\right]

.. math::
    \left. -2|t|^{1/2} G_{1,2}^{2,0}\left[ \frac{|\Delta{}n|^2 \xi}{4k} |t|^{-1} |^{3/2}_{0,1/2}\right] \right\}
