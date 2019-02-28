.. _stats:

.. currentmodule:: multi_locus_analysis

Trajectory Statistics
=====================

All utilities for calculating moments of a trajectory (i.e. statistics of
particle motion), are contained within the :mod:`multi_locus_analysis.stats`
module.

.. sidebar:: :mod:`multi_locus_analysis.stats` docstring

    .. automodule:: multi_locus_analysis.stats
        :noindex:

Tutorial
--------

.. tip::

    Like the rest of :mod:`multi_locus_analysis`, these functions are optimized
    for use with long-form :mod:`pandas` DataFrames. The following tutorial
    assumes that your data set is small enough for approximately :math:`n^4`
    data points to fit in memory. In order to take up only linear space, simply
    compose the desired functions instead of applying them one-by-one (which
    keeps around unnecessary intermediates).

Example data
^^^^^^^^^^^^

You can load in the data we will use in this tutorial as follows:

    >>> import multi_locus_analysis as mla
    >>> from mla.examples import burgess
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

    >>> frame_cols = mla.examples.burgess.frame_cols
    >>> traj_cols = mla.examples.burgess.traj_cols

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
^^^^^^^^^^^^^^^^^^

With almost any set of trajectories, we first extract all the velocities that we
are interested in. If we have drift-free, individual trajectories, this can be
done simply as

>>> vels = df.groupby(traj_cols).apply(mla.stats.pos_to_all_vel, xcol='X', ycol='Y', zcol='Z', frame_col='t')

If there are sufficient trajectories at each time point, then sometimes it is
possible to first subtract out the drift first

>>> df.set_index(frame_cols, inplace=True)
>>> for i in ['X', 'Y', 'Z']:
>>>     df['dc'+i] = df[i] - df.groupby(frame_cols)[i].mean()
>>> dc_vels = df.groupby(traj_cols).apply(mla.stats.pos_to_all_vel, xcol='dcX', ycol='dcY', zcol='dcZ', frame_col='frame id')

However, in our example data, we have two trajectories, so it makes
more sense to simply work directly with the distance between them as our
drift-free measurement.

>>> df = mla.pivot_loci(df, pivot_cols=['X', 'Y', 'Z'])
>>> for i in ['X', 'Y', 'Z']:
>>>     df['d'+i] = df[i+'2'] - df[i+'1']
>>> d_vels = df.groupby(cell_cols).apply(mla.stats.pos_to_all_vel, xcol='dX', ycol='dY', zcol='dZ', framecol='t')

By default, :func:`pos_to_all_vel` computes all velocities. If all
non-overlapping velocities are required instead, such as when computing standard
error in a time-averaged MS(c)D measurement, pass the
:code:`force_independence` kwarg.

>>> d_vels = df.groupby(cell_cols).apply(mla.stats.pos_to_all_vel, \
        xcol='dX', ycol='dY', zcol='dZ', framecol='t', \
        force_independence=True)
>>> mscds_per_stage = d_vels.groupby(['locus', 'genotype', 'meiosis',
'delta']).agg(['mean', 'std', 'count'])

So named because the MSD of the distance between two homolog pairs is often
called the MSCD in the meiosis literature (for Mean Squared Change in
Displacement). These can be simply plotted via

>>> for label, data in mscds_per_stage(['locus', 'genotype', 'meiosis',
'delta']):
        data = data.reset_index()
        plt.errorbar(data['delta'], data['mean'], data['std']/np.sqrt(data['count']))


Velocity correlations
^^^^^^^^^^^^^^^^^^^^^

TODO: finish documentation


Notes
^^^^^

For a listing of all functions exposed
by the :mod:`stats` module, see the :ref:`API reference <api>` section.

