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


>>> from multi_locus_analysis import finite_window as fw
>>> trajs = fw.ab_window_fast([scipy.stats.expon(3).rvs, scipy.stats.expon(1).rvs],
>>>                           [3, 1], 1, num_replicates=1000)
>>> trajs.head()
   replicate  state  start_time  end_time  window_start  window_end
0          0      0   -2.078675  1.556861             0           1
1          1      0   -2.082373  0.968224             0           1
2          1      1    0.968224  2.603759             0           1
3          2      0   -2.079407  0.980203             0           1
4          2      1    0.980203  2.030800             0           1

from each trajectory, we can extract the wait times:

>>> waits = trajs.groupby('replicate').apply(fw.state_changes_to_wait_times)

Rescaling Exact Waiting Times
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

    There will eventually be documentation here about how to use ecdf_window to
    exactly compute the histogram of waiting times from `waits`. For now,
    however, we focus on the more difficult case (so that it can be reviewed at
    lab meeting tomorrow) of movies.

A Caveat for Discrete Movies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

While we would ideally have an exact measurement of when transitions between
`A` and `B` states happen, it is more often the case that we have a "movie" of
sorts: where we measure the state of the system at a fixed set of times.

This only provides us with upper and lower bounds for the actual waiting time
that we're trying to observe. For example

::

        A A A A A B B B B B B B A A A A B
        | - - - | - - - | - - - | - - - |  ...
        |                               |
    t = 0s      1s      2s      3s      4s

This trajectory, when measured, would look like

>>> pd.Series({0: 'A', 1: 'A', 2: 'B', 3: 'A', 4: 'B'}).head()
    0    A
    1    A
    2    B
    3    A
    4    B
    dtype: object

Naively, this would seem to suggest that there was an exterior-censored "A" of
length 2, one each interior censored times of length 1, and one
exterior-censored "B" time of length 1. However, by looking at the true
trajctory above, we see that the first "A" wait was much shorter than 2s, the
first "B" wait was much longer than 1s, and the last "A" wait just happened to
match up with our prediction of 1s.

Because our normalization factor depends non-linearly on the observed waiting
time, one might guess that simply using the "naive" times might cause bias. We
will show that this is the case by generating some artificial movies ourselves.

Generating Discrete Trajectories (Movies)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:mod:`multi_locus_analysis.finite_window` includes a convenience method for
generating "movies" from the output of the `AB_window*` functions.

.. autofunction:: state_changes_to_trajectory
    :noindex:

This function has an alias for convenience (`fw.traj_to_movie`).

>>> movie = trajs.groupby('replicate').apply(traj_to_movie,
times=np.linspace(0, 1, 10))
>>> movie.head()


