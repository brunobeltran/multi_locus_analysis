.. _theory:

.. currentmodule:: multi_locus_analysis

Analytical theory
=================

Here, we show some best practices for comparing MSDs, displacement
distributions, and velocity correlations to the theory of diffusing molecules.

The literature contains many good reviews on these topics, so we cover
only some often overlooked pitfalls.

First, we review the most-used :ref:`model <model-overview>` of viscoelastic
diffusion (for both free molecules and large polymers).

The first common pitfall we address is the :ref:`effect of confinement
<analytical-msds>` on MSDs, especially on extracting diffusivity and the
subdiffusivity coefficient (i.e. Hurst index).

Then, we discuss the effects of heterogeneity between measurements on the
:ref:`gaussianity of the displacement distribution <analytical-disps-hist>` and
what that looks like for real data.

Finally, we discuss how to use the velocity cross correlation of two particles
to :ref:`extract information <analytical-velocity-correlation>` about whether
those two particles are connected (such as by being on the same polymer) using a
drift-independent approach.

.. _model-overview:

Modeling viscoelastic diffusion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The crowded cellular environment is often subdiffusive, and across bacteria, in
some yeast and some higher mammalian cells, it has been shown that this
diffusion is well described by a fractional Brownian motion (for example, see
`(Weber et. al., PRL, 2010)
<https://web.stanford.edu/~ajspakow/publications/WST-ecoli.pdf>`_).

A particle undergoing fractional Brownian motion can be modeled by a fractional
Langevin equation (ignoring inertial terms)

.. math::

    \partial_t^\alpha X(t) = F_B(t)

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

The right-hand side of this equation is called the "Brownian force" and is a
gaussian process with mean zero and covariance determined by the
fluctuation-dissipation theorem

.. math::

    \langle F_B(t) F_B(t') \rangle = \xi k_B T K(t - t') I.

The viscoelastic memory kernel is parameterized by the single parameter
:math:`\alpha`, which describes how much the system wishes to return
(elastically) to a previous configuration.

The simplest way to incorporate the effects of being on a polymer is to use the
Rouse model for polymer dynamics. This model is a formal way to describe how a
so-called `"Gaussian chain" <https://en.wikipedia.org/wiki/Ideal_chain>`_ (i.e.
the limit as the number of beads goes to infinity of a model where the polymer
is composed of beads :math:`X_i(t)` connected by springs) behaves when
undergoing thermal fluctuations.

.. math::

    \partial_t^\alpha X(n,t) = \frac{3k_B T}{b^2} \partial_n^2 X(n, t) + F_B(t)

where :math:`b` is called the "Kuhn length" of the Gaussian chain, and describes
(for a fixed time) how changing the length of the polymer changes the volume
occupied by the polymer.

While the Gaussian chain (and by proxy, the Rouse model) is often unrealistic
for isolated, short polymer chains, it can be shown mathematically (for example,
see the section on Gaussian chains in Doi and Edwards) that every polymer has a
Kuhn length, and behaves identically to a Gaussian chain on long enough length
scales. Therefore, for the case of two loci on chromatin, separated by more than
a couple of kilobases, we argue that the Rouse model is more than sufficient.

.. note::

    Notice that this equation is very reminiscent of the
    Fokker-Planck-Kolmogorov equation for continuous time random walk
    (CTRW).

    A the probability distribution for the position of a free particle
    undergoing a CTRW where the step sizes are Gaussian and the waiting times
    between steps are distributed like :math:`\hat{\phi}(s) = 1/(1 +
    1/\hat{\eta}(s))` (where the hat denotes a Laplace transform), can be
    modeled via

    .. math::

        \partial_t p(x, t) = \partial_t \int_0^t \eta(t - t') \partial_x^2 p(x,
        t') dt'

    In the case where :math:`\eta` is a power law (like :math:`t^{1-\alpha}`),
    this reduces to our Langevin equation above, with the Brownian force
    removed.

    These equations are fundamentally different however, because in our
    model the fluctuation dissipation theorem connects the Brownian force
    :math:`F_B` to the viscosity. This means the velocity correlation for a fBm
    particle is given by :math:`\langle V(t) V(0) \rangle \propto -t^\alpha`,
    whereas since the CTRW takes independent steps by construction, and so has
    :math:`\langle V(t) V(0) \rangle = 0` almost everywhere.

    In what follows, we consider the diffusion to be well described by a simple
    fractional Brownian motion, with no random waiting times between collisions.
    The combination of these two models ((f)Bm and a CTRW) is discussed to some
    degree in the Supplemental Materials of (Beltran and Kannan, PRL 2019), and
    represents well a random walk with heavy-tailed "defects".


.. _analytical-msds:

MSDs
^^^^

The MSD of a particle diffusing according to the equations above can be written
as

.. math::

    \frac{3k_BT}{\xi}\frac{\sin(\alpha\pi)}{\pi(1-\alpha/2)(1-\alpha)\alpha} t^\alpha

In other words, in log-log space, the MSD will be a line with slope
:math:`\alpha`.

.. note::

    While this is true even for arbitrarily short times in our
    fractional diffusion model (due to its fractal nature), for real particles
    diffusing in actual viscoelastic materials (like beads in a polymer gel), as
    we probe shorter and shorter times, there should be a time scale on which
    the elastic forces of the material are not felt. At these size scales in
    real data we will often see an MSD that looks like a regular diffusion (a
    line with slope :math:`\alpha = 1` in log-log space). At even shorter
    times (usually too short to measure), the particle will eventually exhibit
    ballistic behavior, as it travels processively between collisions between
    the molecules that surround it.

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

.. plot::
    :context: close-figs

    >>> from multi_locus_analysis import analytical
    >>> alpha = 0.8; N = 100; b = 1
    >>> t = np.logspace(-6, 6, 100)
    >>> msd = analytical.rouse_mid_msd(t, alpha, b, N, num_modes=1000)
    >>> plt.loglog(t[t>1], 3/N*np.sin(alpha*np.pi)/(
    >>>     np.pi*(1-alpha/2)*(1-alpha)*alpha)*np.power(t[t>1], alpha),
    >>>     '-.', label=r'$\propto t^\alpha$')
    >>> plt.loglog(t, np.power(t, alpha/2), '-.', label=r'$t^{\alpha/2}$')
    >>> plt.loglog(t, msd, 'k', label='MSD(t)')
    >>> plt.xlabel('time (AU)')
    >>> plt.ylabel(r'MSD ($\mu{}m^2$)')
    >>> plt.legend()

Note that the drop-off on the left-hand side of the MSD plot is simply due to
the finite number of Rouse modes included when numerically evaluating the MSD.
The Rouse polymer is a "fractal" model, so the proper solution to the equations
of motion presented above has an MSD whose short-time behavior scales as
:math:`t^{\alpha/2}` indefinitely. The log of the number of modes included
(optional arg `num_modes` in the call) corresponds exactly to the number of
orders of magnitude that will have the correct scaling in the calculated MSD.

However, this drop-off is also instructive, because this fractal behavior is exactly
the most unrealistic part of the Rouse model. In real systems, at short enough
time scales, a monomer on the polymer will also transition back to the
background slope of :math:`\propto\alpha`. This means that if you measure an actual
locus on a polymer, say a particular locus of a small plasmid in a mammalian
cell (where :math:`\alpha\approx 1`). Then we will see a curve that looks just
like the black curve above. The long time crossover (from :math:`t^{\alpha/2}`
to :math:`t^\alpha`) will measure the terminal relaxation time of the plasmid
(the time scale after which diffusion of a locus on the plasmid is dominated by
diffusion of the plasmid as a whole). The short crossover time (from :math:`t^\alpha` to :math:`t^{\alpha/2}`) will correspond to the time scale at which thermal forces from the surrounding medium begin to be dominated by the elastic restoring force from the polymer chain.


.. note::

    The fractional Brownian motion itself is also a fractal model, and as such
    has the similar shortcomings to the Rouse model. Namely, at long enough time
    scales, if confinement isn't a factor, it is likely that the MSD will
    transition eventually to a slope of :math:`\alpha=1` for any realistic
    system. This may also be true for short enough times. In a sense, what one
    can say is that the fBm Rouse model is really telling us is that for the
    time scales described above, whatever the medium value of :math:`\alpha` is,
    being on a polymer simply halves it. Thus it is always important to compare
    polymer MSDs to the MSD of a free particle in the same medium, since many
    complex mediums will not have a single characteristic value of
    :math:`\alpha`, but instead have an MSD whose slope depends on the time or
    length scale (e.g. due to defects or traps of a particular size scale).
    If comparison to a free particle is not possible experimentally, there are
    also other ways to estimate the terminal relaxation time of the polymer if
    this is not possible (see velocity cross correlation discussion below).

So in experimental MSDs, where the chromosomal loci are confined to a given
spatial region (whether that be a chromosome territory or the nucleus as a
whole), we expect to see the MSD curve cut off at a different location depending
on how the confinement time scale compares to the terminal relaxation time of
the chromosome (or chromosomes, if they are linked together).

TODO: make plot showing three options (plateau before polymer regime, before
terminal relaxation time, and after terminal relaxation time).

As can be seen from the plot above, even for this simplest possible model of
Rouse polymer motion, simply fitting MSD curves with straight lines in log-log
space can lead to extremely misleading results, due to the non-linear nature of
the MSD curve itself. Therefore, it is important to establish the confinement
diameter, and terminal relaxation time of the polymer at a minimum, before
trying to fit a raw MSD to extract, for example, a value for :math:`\alpha`.

.. _analytical-disps-hist:

Displacement Distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO: document gaussianity vs laplace distributions of dispalcements here.

The `paper of Lampo et al
<http://stanford.edu/~ajspakow/publications/LSBWS-celldiff.pdf>`_ is extremely
valuable for understanding how heterogeneity in diffusivities between different
cells can lead to a population that appears to have a non-guassian displacements
distribution even though each underlying cell may itself have gaussian
displacements.

A more general overview of the combined effects of intra-cellular and
inter-cellular heterogeneity can be found in the `later paper by Stylianidou and
Lampo <https://journals.aps.org/pre/abstract/10.1103/PhysRevE.97.062410>_`.

This section will eventually summarize how these issues appear in practice,
where the distribution of diffusivities is often not as clean as in the data
that Lampo et al present.

.. _analytical-velocity-correlation:

Velocity Cross-Correlation
^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to measure whether (and how) stress is communicated between the loci
being measured, we can look at their velocity cross correlation functions.

The velocity autocorrelation of a free particle undergoing fractional Brownian
motion will serve as a useful comparison throughout what follows in order to
understand limiting behavior of the more complicated, polymer case. Such a
particle has a velocity autocorrelation function given by

.. math::

    \frac{C_v^{(\delta)}(t)}{C_v^{(\delta)}(0)} = \frac{ | t - \delta|^\alpha -
    2|t|^\alpha + | t + \delta|^\alpha}{2\delta^\alpha}

where we should notice that :math:`C_v^{(\delta)}(0)` is just the MSD as a
function of :math:`\delta`. We plot this function below for a few values of
:math:`\alpha`.

.. plot::
    :context: close-figs

    >>> td = np.linspace(0, 4, 401)
    >>> betas = np.concatenate([np.linspace(0, 1, 6), np.array([1.5, 2])])
    >>> betas[0] += 0.0001
    >>> cmap = plt.get_cmap('viridis')
    >>> cnorm = mpl.colors.Normalize(vmin=0, vmax=1.5)
    >>> for beta in betas:
    >>>     label = f'$\\beta = {beta:0.1}$' if beta < 1 else f'$\\beta = {beta}$'
    >>>     plt.plot(td, analytical.frac_discrete_cv_normalized(td, 1, beta),
    >>>              label=label, c=cmap(cnorm(beta)))
    >>> plt.xlabel(r'$t/\delta$')
    >>> plt.ylabel(r'$C_v^{(\delta)}(t)$')
    >>> plt.legend(loc='upper right')

Notice that as :math:`\alpha` goes to 1, the velocity correlation goes to that
of a regular Brownian motion, :math:`\max\{1 - t/\delta, 0\}`. As :math:`\alpha`
goes to 2, the velocity correlation becomes identically 1 (as expected for
ballistic motion). As :math:`\alpha` goes to 0, the velocity correlation is
identically zero, with unit spikes at zero and one in the positive and negative
directions, respecitvely (as expected for totally arrested motion).

As shown in `(Lampo & Spakowitz, Biophysical Journal, 2016)
<https://www.sciencedirect.com/science/article/pii/S0006349515047049>`_, the
velocity cross-correlation of two loci on a Rouse polymer can be computed
exactly.

.. math::

    C_{vv}^{(\delta)}(t, n_1, n_2) = \left\langle \vec{V}^{(\delta)}(n_1, t)
    \cdot \vec{V}^{(\delta)}(n_2, 0) \right\rangle

after a cosine transform in :math:`n`, and expanding in the eigenbasis of the fractional derivative, :math:`E_{\alpha,1}` (the Mittag-Leffler functions), we get

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

.. The dependence on the location along the polymer is captured by these Rouse
.. modes, as seen below

.. .. plot::
..     :context: close-figs

..     >>> from multi_locus_analysis import analytical
..     >>> ps = np.arange(1000)
..     >>> modes = analytical.rouse_mode(ps, 0.5, 1)
..     >>> plt.scatter(ps, modes)
..     >>> modes = analytical.rouse_mode(ps, 0.51, 1)
..     >>> plt.scatter(ps, modes)
..     >>> modes = analytical.rouse_mode(ps, 0.91, 1)
..     >>> plt.scatter(ps, modes)
..     >>> plt.xlabel('p')
..     >>> plt.ylabel('$\phi_p(n/N)$')
..     >>> plt.legend(['n = 0.5', 'n = 0.51', 'n = 0.91'])

The coefficients :math:`k_p = \frac{3\pi^2k_BT}{Nb^2} p^2` then determine how
strongly to weight each mode as a function of :math:`t` and :math:`\delta`.

Because :math:`E_{\alpha,1}(x)\to 0` as :math:`x\to-\infty`, and since
:math:`k_p \propto p^2`, the summation converges fairly quickly. However, using
the numerical values for :math:`\phi_p(n_i)` tends to cause numerical
instabilities. Instead, we consider several limiting cases where properties of
the cosine function allow us to simplify the expression.

Firstly, the polymer locus's autocorrelation function, i.e. :math:`n_1 = n_2`
is just

.. math::

    TODO

The unnormalized version of this equation can be
found in Eq. 33 of `(Weber & Spakowitz
Physical Review E, 2010)
<https://journals.aps.org/pre/pdf/10.1103/PhysRevE.82.011913>`_.

Infinite Polymer Limit
----------------------

In the simple, but typical case where :math:`\alpha = 1`, the infinite polymer
case becomes significantly simpler than the general case of a finite polymer.

We can express the velocity correlation of two loci on an infinite polymer in a
regular viscous medium in terms of the Meijer G-functions, which are defined as
the Mellin-Barnes integrals

.. math::

    G^{m,n}_{p,q} \left( \left. \begin{matrix} a_1, \dots, a_n ; a_{n+1} \dots a_p \\ b_1, \dots, b_m ; b_{m+1} \dots b_q \end{matrix}\; \right| \; z ; r \right) = \frac{1}{2 \pi i} \int_L \frac{\prod_{j=1}^m \Gamma(b_j+s) \prod_{j=1}^n\Gamma(1-a_j-s)} {\prod_{j=n+1}^{p}\Gamma(a_j+s) \prod_{j=m+1}^q \Gamma(1-b_j-s)} z^{-s/r} ds

For two loci separated by a distance :math:`\Delta n`, the correlation function
is simply given by

.. math::
    C_{vv}^{(\delta)}(t, \Delta n) = \frac{3k_BT}{\delta^2\sqrt{\xi k}}

.. math::
    \times \left\{ | t - \delta|^{1/2} G_{1,2}^{2,0}\left[\left. \frac{| \Delta{}n|^2 \xi}{4k} | t - \delta|^{-1} \right|^{3/2}_{0,1/2}\right] \right.

.. math::
    + | t + \delta|^{1/2} G_{1,2}^{2,0}\left[ \frac{| \Delta{}n|^2 \xi}{4k} | t + \delta|^{-1} | ^{3/2}_{0,1/2}\right]

.. math::
    \left. -2|t|^{1/2} G_{1,2}^{2,0}\left[ \frac{| \Delta{}n|^2 \xi}{4k} | t|^{-1} | ^{3/2}_{0,1/2}\right] \right\}.


This function can be evaluated exactly using any numerical library that has
routines for evaluating Meijer G-functions. In particular, in Python, this
functionality is implemented in the well-tested package `mpmath`.

Velocity Cross-Correlation "Trick"
----------------------------------

It is simple enough to notice that


