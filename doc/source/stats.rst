.. _stats:

.. currentmodule:: multi_locus_analysis

Trajectory Statistics
=====================

All utilities for calculating moments of a trajectory (i.e. statistics of
particle motion), are contained within the :mod:`multi_locus_analysis.stats`
module.

This statement is made more precise in that module's docstring:

:mod:`multi_locus_analysis.stats`
---------------------------------
    .. automodule:: multi_locus_analysis.stats
        :noindex:

Tutorial
--------

.. tip::

    Like the rest of :mod:`multi_locus_analysis`, these functions are optimized
    for use with long-form :mod:`pandas` DataFrames. The following tutorial
    assumes that your data set is small enough for approximately :math:`n^3`
    data points to fit in memory. In order to take up only linear space, simply
    compose the desired functions instead of applying them one-by-one (which
    keeps around unnecessary intermediates).

Example data
^^^^^^^^^^^^

You can load in the data we will use in this tutorial as follows:

.. plot::
    :nofigs:
    :context: close-figs

    >>> import multi_locus_analysis as mla
    >>> from multi_locus_analysis.examples import burgess
    >>> df = mla.examples.burgess.df
    >>> df.head()
                                                      X    Y    Z foci    t
    locus genotype exp.rep meiosis cell frame spot
    HET5  WT       2       t0      1    1     1     1.6  2.4  3.1  unp    0
                                        2     1     1.9  1.5  3.1  unp   30
                                        3     1     2.0  1.8  3.0  unp   60
                                        4     1     2.1  1.9  3.0  unp   90
                                        5     1     2.2  1.8  3.0  unp  120

As is typical in this kind of data, there are some number of columns holding
the actual measurements, and some number which are simply record-keeping
columns. In the Burgess data a uniquely-identifying-subset of the
record-keepign columns have been pre-emptively made into an index.

We will be using named collections of columns throughout this tutorial to group
together parts of the DataFrame as needed (for example, to get each
continuous trajectory we use :code:`traj_cols`). Please see
:mod:`multi_locus_analysis.examples.burgess` for detailed exaplanations of the
information in each column. For now, we need only point out that each unique
combination of ``['locus', 'genotype', 'exp.rep', 'meiosis', 'cell']`` defines
a movie, which contains two spots diffusing relative to each other: these are
homologous loci at different times in meiosis. To save some typing, we
extract the groupings we will use later:

.. plot::
    :nofigs:
    :context:

    >>> cell_cols = burgess.cell_cols
    >>> traj_cols = burgess.traj_cols
    >>> frame_cols = burgess.frame_cols

Different tasks will be easier or harder with the index set to identify a
single spot at a single time, or a single time point with different columns
for each spot. We can use the :func:`multi_locus_analysis.pivot_loci` helper to
move from one representation to the other

    >>> mla.pivot_loci(df, pivot_cols=['X', 'Y', 'Z']).head()
                                              foci    t   X1   Y1   Z1   X2   Y2   Z2
    locus genotype exp.rep meiosis cell frame
    HET5  WT       2       t0      1    1      unp    0  1.6  2.4  3.1  2.1  1.5  3.1
                                        2      unp   30  1.9  1.5  3.1  1.9  2.5  3.1
                                        3      unp   60  2.0  1.8  3.0  1.5  2.5  3.4
                                        4      unp   90  2.1  1.9  3.0  1.4  2.2  3.4
                                        5      unp  120  2.2  1.8  3.0  1.5  2.4  3.4

Velocities to MS(C)Ds
^^^^^^^^^^^^^^^^^^^^^

With almost any set of trajectories, we first extract all the velocities that we
are interested in. If we have drift-free, individual trajectories, this can be
done simply as

    >>> # subsampling makes computation easy for demo purposes
    >>> # group by remaining columns that uniquely id a trajectory
    >>> vels = df.loc['URA3', 'WT', 9].groupby(traj_cols[3:]).apply(
    >>>     mla.stats.pos_to_all_vel, xcol='X', ycol='Y', zcol='Z', framecol='t')
    >>> vels.head()
                                 tf   vx   vy   vz
    meiosis cell spot ti delta
    t0      2    1    0  0        0  0.0  0.0  0.0
                         30      30  0.0 -0.1  0.3
                         60      60  0.0  0.0  0.3
                         90      90  0.1 -0.1  0.2
                         120    120  0.2 -0.2  0.3

If there are sufficient trajectories at each time point, then sometimes it is
possible to first subtract out the drift first

>>> raise NotImplementedError('Pandas shape error in following code')
>>> for i in ['X', 'Y', 'Z']:
>>>     df['dc'+i] = df[i] - df.groupby(frame_cols)[i].mean()
>>> # groupby remaining cols that uniquely id a trajectory
>>> dc_vels = df.loc['HET5', 'WT', 2].groupby(traj_cols[3:]).apply(
>>>     mla.stats.pos_to_all_vel, xcol='dcX', ycol='dcY', zcol='dcZ',
>>>     framecol='t')
>>> dc_vels.head()

However, in our example data, we have two trajectories, so it makes
more sense to simply work directly with the distance between them as our
drift-free measurement.

.. plot::
    :nofigs:
    :context:

    >>> df = mla.pivot_loci(df, pivot_cols=['X', 'Y', 'Z'])
    >>> for i in ['X', 'Y', 'Z']:
    >>>     df[i+'12'] = df[i+'2'] - df[i+'1']
    >>> df.head()

We can then get all the :math:`V^\delta_{12}(t_k)` by simply using these new
:math:`X_{12}(t_k)` columns:

>>> d_vels = df.groupby(cell_cols).apply(mla.stats.pos_to_all_vel, xcol='dX', ycol='dY', zcol='dZ', framecol='t')

By default, :func:`pos_to_all_vel` computes all velocities. If all
non-overlapping velocities are required instead, such as when computing standard
error in a time-averaged MS(c)D measurement, pass the
:code:`force_independence` kwarg.

.. plot::
    :nofigs:
    :context:

    >>> # first subsample so that this runs quickly on a single machine
    >>> df = df.loc['URA3', 'WT', 9].copy() # one experiment
    >>> remaining_index = cell_cols[3:]
    >>> d_vels = df.groupby(remaining_index).apply(mla.stats.pos_to_all_vel, \
    >>>     xcol='X12', ycol='Y12', zcol='Z12', framecol='t', \
    >>>     force_independence=True)
    >>> d_vels['v'] = np.sqrt(d_vels['vx']**2 + d_vels['vy']**2 + d_vels['vz']**2)
    >>> mscds_per_stage = d_vels.groupby(['meiosis', 'delta'])['v'].agg(['mean', 'std', 'count'])

The MSD of the distance between two homolog pairs is often called the MSCD in
the meiosis literature (standing for Mean Squared Change in Displacement).
These MSCDs can be simply plotted via

.. plot::
    :context:

    >>> plt.figure(figsize=[4, 3])
    >>> for label, data in mscds_per_stage.groupby('meiosis'):
    >>>     data = data.reset_index()
    >>>     data = data[data['delta'] > 0] # for log-log scale prettiness
    >>>     # SEM bars
    >>>     plt.errorbar(data['delta']/60, data['mean'], \
    >>>                  data['std']/np.sqrt(data['count']), \
    >>>                  label=label)
    >>> plt.yscale('log')
    >>> plt.xscale('log')
    >>> plt.ylim([0.01, 1])
    >>> plt.xlabel(r'time lag, $\delta$ (min)')
    >>> plt.ylabel(r'MSCD, $<V_{12}^\delta(t)>$ ($\mu{}m^2/s$)')
    >>> plt.legend()


Velocity correlations
^^^^^^^^^^^^^^^^^^^^^

In order to calculate the velocity autocorrelation, we include the
function :func:`all_vel_to_corr`, which can be used to get numbers that can be groupby'd and averaged to get the velocity autocorrelation.

However, doing this naively on any reasonably-sized dataset will result in way
too large of a DataFrame to fit in RAM. This is because this function is simply
calculating :math:`V^\delta(t + t_k) \cdot V^\delta(t_k)` for every single
:math:`\delta` and :math:`t_k`. While this works well for small datasets (like
single movies), even clever use of Pandas will often make applying this
function to a large dataset spicy. Therefore, we provide convenience functions
:func:`vels_to_cvvs_by_hand` and :func:`vvc_stats_by_hand` to first write out
all these products and then take the appropriate averages (respectively).

If we wanted to use them on the above velocities, we would simply run:

.. warning::

    The following creates a temporary file of approximately 2.1GB.

.. code-block:: python

    >>> mla.stats.vels_to_cvvs_by_hand(d_vels, ['meiosis'], 'vvs.csv',
    >>>                                max_t_over_delta=4)
    >>> cvv_stats = mla.stats.vvc_stats_by_hand('vvs.csv', ['meiosis'])
    >>> cvv_stats = mla.stats.cvv_by_hand_make_usable(cvv_stats, ['meiosis'])
    >>> cvv_stats.set_index(['meiosis', 'delta', 't'], inplace=True)
                             sqs      sum     ...      cvv_normed  ste_normed
    meiosis delta t                           ...
    t0      30    0    8787.4250  1427.42     ...        1.000000    0.065668
                  30   6635.5702  -570.96     ...       -0.421770    0.060173
                  60   6503.1084   -88.80     ...       -0.067135    0.060968
                  90   6000.1380   -76.52     ...       -0.058923    0.059647
                  120  5884.4852   191.24     ...        0.150432    0.060341

We can then plot the velocity correlation curves (for example at meiotic
stage ``t0``) by simply running

    >>> from multi_locus_analysis import plotting
    >>> plotting.cvv_plot_sized(cvv.loc['t0'].reset_index(), data_deltas=30*np.arange(1, 5))

.. note::

    There should be a plot here but it has not been inserted yet pending come
    organizational decisions.

where we ignore larger values of ``delta`` because the amount of data available
for averaging decreases to a single point as ``delta`` approaches the length of
the observation.

Notice how the correlation decays as $\delta$ increases. For details on why
this means that the correlation is actually increasing as delta increases,
please see the paper in preparation (Newman, Beltran, and Calhoon et al).

Notes
^^^^^

For a listing of all functions exposed
by the :mod:`stats` module, see the :ref:`API reference <api>` section.

