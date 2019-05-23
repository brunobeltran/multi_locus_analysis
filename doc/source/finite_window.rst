.. _finite_window:

.. currentmodule:: multi_locus_analysis.finite_window

Finite window correction
========================

Measurements of "waiting times" or "survival times" taken within a finite
time interval require special statistical treatment to account for the bias
towards short measurements introduced by the measurement itself. For
mathematical details, see the manuscript in preparation `(Beltran et. al.,
*in preparation*)
<https://www.overleaf.com/project/5c6c76d88d00f95294959995>`_.

Otherwise, the process for these corrections in explained in the
module's doctring

:mod:`multi_locus_analysis.finite_window`
-----------------------------------------
    .. automodule:: multi_locus_analysis.finite_window
        :noindex:

Tutorial
--------

Generating AB Trajectories
^^^^^^^^^^^^^^^^^^^^^^^^^^

We use the following two functions to generate example data:

.. autofunction:: multi_locus_analysis.finite_window.ab_window_fast
    :noindex:

.. autofunction:: multi_locus_analysis.finite_window.ab_window
    :noindex:


In short, the former is useful if you can compute the means of the waiting time
distributions, and the latter is more general (but *much* slower).


.. plot::
    :nofigs:
    :context: close-figs

    >>> import scipy.stats
    >>> from multi_locus_analysis import finite_window as fw
    >>> trajs = fw.ab_window_fast([scipy.stats.expon(scale=1).rvs,
    >>>     scipy.stats.expon(scale=3).rvs], [1, 3], 5, 1000)
    >>> trajs.head()
    replicate  state  start_time  end_time  window_start  window_end
    0          0      0   -0.008954  0.641606             0           5
    1          0      1    0.641606  2.593285             0           5
    2          0      0    2.593285  4.291466             0           5
    3          0      1    4.291466  9.386009             0           5
    4          1      0   -1.117854  0.455587             0           5

from each trajectory, we can extract the wait times:

.. plot::
    :nofigs:
    :context:

    >>> waits = trajs.groupby('replicate').apply(fw.traj_to_waits)
    >>> waits.head()
    state_changes_to_wait_timecate  state  start_time  end_time  window_start  window_end  wait_time  window_size       wait_type
    replicate rank_order
    0         0                   0      0    0.000000  0.641606             0           5   0.641606            5   left exterior
            1                   0      1    0.641606  2.593285             0           5   1.951679            5        interior
            2                   0      0    2.593285  4.291466             0           5   1.698181            5        interior
            3                   0      1    4.291466  5.000000             0           5   0.708534            5  right exterior
    1         4                   1      0    0.000000  0.455587             0           5   0.455587            5   left exterior

Rescaling Exact Waiting Times
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

    There will eventually be documentation here about how to use ecdf_window to
    exactly compute the histogram of waiting times from `waits`. For now,
    however, we focus on the more difficult case (so that it can be reviewed at
    lab meeting tomorrow) of movies.

As you can see, even if we only take "interior" times, the empirical CDF of the
data does not match the distribution that we put in (`expon(scale=1)` in this
case, from above).

.. plot::
    :context:

    >>> plt.figure(figsize=[4, 3])
    >>> interior_0 = waits.loc[
    >>>     (waits['wait_type'] == 'interior') & (waits['state'] == 0),
    >>>     'wait_time'
    >>> ].values
    >>> x, cdf = fw.ecdf(interior_0)
    >>> plt.plot(x, cdf, label='Empirical CDF, state 1')
    >>> plt.plot(x, scipy.stats.expon().cdf(x), 'k-.', label='Actual CDF')
    >>> plt.xlabel('t')
    >>> plt.ylabel(r'$P(\mathrm{wait} <= t)$')

We can correct this discrepancy by simply using
:func:`multi_locus_analysis.finite_window.ecdf_windowed`.

.. plot::
    :context: close-figs

    >>> plt.figure(figsize=[4,3])
    >>> x, cdf = fw.ecdf_windowed(interior_0, 5)
    >>> plt.plot(x, cdf*scipy.stats.expon().cdf(5), label='Empirical CDF, state 1')
    >>> plt.plot(x, scipy.stats.expon().cdf(x), 'k-.', label='Actual CDF')
    >>> plt.xlabel('t')
    >>> plt.ylabel(r'$P(\mathrm{wait} <= t)$')


A Caveat for Discrete Movies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

While we would ideally have an exact measurement of when transitions between
`A` and `B` states happen, it is more often the case that we have a "movie" of
sorts: where we measure the state of the system at a fixed set of times.

This only provides us with upper and lower bounds for the actual waiting time
that we're trying to observe. For example, consider the trajectory depicted
below.

::

        A A A A A B B B B B B B A A A A B
        | - - - | - - - | - - - | - - - |  ...
        |                               |
    t = 0s      1s      2s      3s      4s

This trajectory, when measured at the discrete times shown, would look like

>>> pd.Series({0: 'A', 1: 'A', 2: 'B', 3: 'A', 4: 'B'}).head()
    0    A
    1    A
    2    B
    3    A
    4    B
    dtype: object

Naively, if you only had this movie in front of you with no knowledge of the
actual underlying state change times, it might seem to suggest that there was an
exterior-censored "A" of length 2, one each interior censored times of length 1,
and one exterior-censored "B" time of length 1. However, by looking at the true
trajctory above, we see that the first "A" wait was much shorter than 2s, and
the first "B" wait was much longer than 1s, whereas the last "A" wait just
happened to match up with our prediction of 1s.

Because our normalization factor depends non-linearly on the observed waiting
time, one might guess that simply using the "naive" times might cause bias. We
will show that this is the case by generating some artificial movies ourselves.

Generating Discrete Trajectories (Movies)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:mod:`multi_locus_analysis.finite_window` includes a convenience method for
generating "movies" from the output of the `AB_window*` functions.

.. autofunction:: state_changes_to_trajectory
    :noindex:

This function has an alias for convenience
(:func:`multi_locus_analysis.finite_window.traj_to_movie`).

>>> movies = trajs.groupby('replicate').apply(traj_to_movie,
times=np.linspace(0, 1, 10))
>>> movies.head()
t          0.0  0.5  1.0  1.5  2.0  2.5  3.0  3.5  4.0  4.5  5.0
replicate
0            0    0    1    1    1    1    0    0    0    1    1
1            0    1    1    1    1    1    1    1    1    1    1
2            1    1    1    1    0    1    1    1    0    1    1
3            0    0    1    1    1    1    1    1    1    1    0
4            0    0    0    1    1    1    1    1    1    1    0

We can get a long-form DataFrame by simply doing

>>> movies.T.unstack()
>>> movies.name = 'state' # name resulting column from unstack
>>> movies.head()
replicate  t
0          0.0    0
           0.5    0
           1.0    1
           1.5    1
           2.0    1
Name: state, dtype: int64

As is clear from the following plot, the data being effectively
discretized creates a bias in the tail of the distribution, even when the
times are corrected with our method.

.. TODO : finish

