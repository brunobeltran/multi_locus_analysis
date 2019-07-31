"""Show various versions of the MS(C)D plots for the Burgess data, do some
simple analysis to extract effective diffusion coefficients, confinement
radii, etc.

The differences between MSD and MSCD, as well as those caused by splitting
the trajectories as opposed to treating bound and unbound as different
states of the same trajectory are explored.


There are several reasonable definitions for what a "trajectory" is in this
data.

    1) spot==1 positions and spot==2 positions, including the time when they
    spend being the "same spot", while they are not broken by internal NaN.
    2) spot==1 positions for a single cell, while they are not broken by an
    internal NaN. spot==2 positions only invalid when spots are paired, but
    assuming it makes sense to calculate velocities when at an intermediate
    time the spots where paired.
    3) spot==1 positions for a single cell, while they are not broken by an
    internal NaN. spot==2 positions only valid in between times they are
    paired, assuming that it makes no sense to calculate velocities between
    two unpaired times when there is an intermediate "pair" time.
    4) spot==1 positions and spot==2 positions, only when they are not bound
    together. separately, spot==1&2 positions.

(1) is easily extracted from the default state of burgess.df, e.g. in
burgess.make_all_intermediates, where the MSCD is saved in msds_movies_file.

(2) requires breaking up spot==2 trajectories by "wait_id", defined in add_wait_id

(3) requires simply setting spot==2 to NaN when paired, not breaking up by
wait_id

(4) requires simply taking "pair" and "unp" times separately into the
calculations. (after breaking up by wait_id).

These different interpretations lead to several different MS(C)D's::

    a) (from (1)) MSCDs after groupby(cell_cols)
    b) (from (2)) MSCDs after groupby(cell_cols), dropping data from 'pair' times.
    c) (from (3)) (MSCDs after groupby(cell_cols, + ['wait_id']). this ignores
    'pair' times implicitly.
    d) (from (1)) MSDs each spot after groupby(cell_cols)
    e) (from (2)) MSDs of spot1 after groupby(cell_cols) and of spot2 ignoring
    "pair" times
    f) (from (3)) MSDs of spot1 after groupby(cell_cols) and of spot2 after
    groupby(cell_cols + ['wait_id'])
    g) (from (4)) MSDs of spot1 and spot2 after groupby(cell_cols+['wait_id']),
    ignoring spot2 whenever foci=='pair'. separately: the MSDs of the spot
    during "paired" times.

In our estimation, which of these makes sense depends on how we interpret what
is happening when the loci are paired. If there is a non-markovianity induced
by pairing, such as if the loci are fixed together by some external force at
that time, then only (c) and (g) likely make any sense.

on the other hand, if the diffusion of each locus does not change when they are
paired (such as would be the case if they were not binding, and only coming
into proximity), then (a/b) and (d/e) make sense, depending on whether our
resolution is high enough that we can safely assign a distance "0" to the loci
for the times at which we see them "paired" (otherwise, we will be skewing our
MSDs to the small side).
"""

from . import *

import scipy.special

from pathlib import Path

def msd(df, mscd=True, include_z=False, traj_group=cell_cols,
        groups=condition_cols, vel_file=None):
    """Catch-all function to compute the various versions of the MSD curves
    that you might want to investigate. In particular, see Notes below for how
    to compute all (a)-(g) versions of the MS(C)Ds described in the module
    docs.

    Parameters
    ----------
    df : pd.DataFrame
        the burgess data
    mscd : bool
        whether or not to calculate MSCDs instead of MSDs
    include_z : bool
        whether or not to use the "z" dimension in the data, which is a little
        weird since the nuclei are not spherical
    traj_group : List[str]
        columns by which to group the data in order to extract one trajectory
        per group. this allows calculation of different types of MSDs, as
        described below.
    groups : List[str]
        columns by which to perform the final groupby to get one MSD curve per
        group. this simply allows you to calculate e.g. one MSD per cell, or
        one MSD per experiment, etc. etc.

    Returns
    -------
    msds : pd.DataFrame
        dataframe with 'delta' and :param:`groups` as the index and ['mean',
        'std', 'count'] as the columns, corresponding to the time-and-ensemble
        averaged quantities for the MSD as a function of 'delta' for each group
        in :param:`groups`.
    """
    if traj_group is None:
        traj_group = cell_cols if mscd else cell_cols+['spot']
    x = 'X'; y = 'Y'; z = 'Z' if include_z else None
    if mscd:
        x = 'd'+x; y = 'd'+y;
        if include_z:
            z = 'd'+z
    all_vel = df \
            .groupby(traj_group) \
            .apply(pos_to_all_vel, xcol=x, ycol=y, zcol=z, framecol='t')
    absv2 = np.power(all_vel['vx'], 2)
    absv2 += np.power(all_vel['vy'], 2)
    if include_z:
        absv2 += np.power(all_vel['vz'], 2)
    all_vel['abs(v)'] = np.sqrt(absv2)
    if vel_file:
        all_vel.to_csv(vel_file)
    if 'delta' not in groups:
        groups = groups + ['delta']
    msds = all_vel.groupby(groups)['abs(v)'].agg(['mean', 'std', 'count'])
    msds['ste'] = msds['std']/np.sqrt(msds['count']-1)
    msds['ste_norm'] = msds['ste']*np.sqrt(2/(msds['count']-1)) \
            *scipy.special.gamma(msds['count']/2) \
            /scipy.special.gamma((msds['count']-1)/2)
    return msds

msd_args = {
    'dvel': {'df': df_flat, 'mscd': True},
    'dvel_unp': {'df': df_flat[df_flat['foci'] == 'unp'], 'mscd': True},
    'dvel_by_wait': {'df': df_flat[df_flat['foci'] == 'unp'],
                     'mscd': True,
                     'traj_group': cell_cols + ['wait_id']},
    'vel_double_counted': {'df': df, 'mscd': False},
    # 'vel_no_double': {'df': df},
    # 'vel_spot2_by_wait': {'df': df},
    'vel_unp_by_wait': {'df': df[df['foci'] == 'unp'],
                        'traj_group': cell_cols + ['spot', 'wait_id'],
                        'mscd': False},
    'vel_pair_by_wait': {'df': df[df['foci'] == 'pair'],
                         'traj_group': cell_cols + ['spot', 'wait_id'],
                         'mscd': False},
}
def precompute_msds(prefix=burgess_dir, force_redo=False):
    """Precomputes a bunch of different "MS(C)Ds"

    The various types of MSDs described in the module documentation can be
    computed as follows (with include_z and groups as you see fit)::

        (a) msd(df_flat, mscd=True)
        (b) msd(df_flat[df_flat['foci'] == 'unp'], mscd=True)
        (c) msd(df_flat[df_flat['foci'] == 'unp'], mscd=True, traj_group=cell_cols+['wait_id'])
        (d) msd(df)
        (e) msd(df[(df['foci'] == 'unp') | (df['spot'] == 1)])
        (f) df.loc[df['spot'] == 1, 'wait_id'] = 0
            msd(df[(df['foci'] == 'unp') | (df['spot'] == 1)],
            traj_group=cell_cols+['spot'+'wait_id'])
        (g) msd(df[df['foci'] == 'unp'], traj_group=cell_cols+['spot'+'wait_id'])
            msd(df[(df['foci'] == 'pair') & (df['spot'] == 1)], traj_group=cell_cols+['wait_id'])

    We refer to the velocities calculated from each of the above subsets of
    the Burgess data by easy-to-remember names::

        (a) dvel
        (b) dvel_unp
        (c) dvel_by_wait
        (d) vel_double_counted
        (e) vel_no_double: NOT CURRENTLY COMPUTED
        (f) vel_spot2_by_wait: NOT CURRENTLY COMPUTED
        (g) vel_unp_by_wait
            vel_pair_by_wait

    """
    for name, kwargs in msd_args.items():
        kwargs['vel_file'] = prefix / Path(name + '.csv')
        msd_file = prefix / Path('msds_' + name + '.csv')
        if not Path(msd_file).exists() or force_redo:
            msds = msd(**kwargs)
            msds.to_csv(msd_file)

