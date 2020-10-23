"""
Wrapper functions to (pre)compute MS(C)Ds, compute diffusivity, nuc_radius.

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

(2) requires breaking up spot==2 trajectories by "wait_id", defined in
``add_wait_id``

(3) requires simply setting spot==2 to NaN when paired, not breaking up by
wait_id

(4) requires simply taking "pair" and "unp" times separately into the
calculations. (after breaking up by wait_id).

These different interpretations lead to several different MS(C)D's::

    a) (from (1)) MSCDs after groupby(cell_cols)
    b) (from (2)) MSCDs after groupby(cell_cols), but ignoring any frames
    during which the loci are 'pair'ed.
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

Instead of showing them all in the final paper, we have chosen to simply use
(b) and (d), and to make sure to process our simulation output in the same way,
so that we can phenomenologically correct for the bias that this causes.
"""
from . import cell_cols, condition_cols, burgess_dir
import multi_locus_analysis.examples.burgess as burgess  # NOQA
from multi_locus_analysis.stats import pos_to_all_vel

import scipy.special
import numpy as np

from pathlib import Path


def msd(df, mscd=True, include_z=True, traj_group=cell_cols,
        groups=condition_cols, vel_file=None, dims=None, **kwargs):
    """
    Catch-all function to compute the various versions of the MSD.

    Parameters
    ----------
    df : pd.DataFrame
        The burgess data.
    mscd : bool
        Whether or not to calculate MSCDs instead of MSDs (effective whether or
        not to prefix intermediate and output dimension column names with 'd').
    include_z : bool
        Whether or not to use the "z" dimension in the data, which is a little
        weird since the nuclei are not spherical
    traj_group : List[str]
        Columns by which to group the data in order to extract one trajectory
        per group. this allows calculation of different types of MSDs, as
        described below.
    groups : List[str]
        Columns by which to perform the final groupby to get one MSD curve per
        group. this simply allows you to calculate e.g. one MSD per cell, or
        one MSD per experiment, etc. etc.
    vel_file : Optional[str]
        A file name for the "vel" file if it is to be saved. This holds the
        intermediate output of ``pos_to_all_vel`` before it is averaged to
        extract the MS(C)Ds.
    dims : List[Optional[str]]
        Manually specify the names of the columns to be used for each dimension
        of position. A None specifies to ignore that column.
    **kwargs : Dict[str, object]
        forwarded to mla.pos_to_all_vel

    Returns
    -------
    msds : pd.DataFrame
        dataframe with 'delta' and *groups* as the index and ['mean',
        'std', 'count'] as the columns, corresponding to the time-and-ensemble
        averaged quantities for the MSD as a function of 'delta' for each group
        in *groups*.

    """
    if traj_group is None:
        traj_group = cell_cols if mscd else cell_cols+['spot']
    if dims is None:
        dims = ['X', 'Y', 'Z'] if include_z else ['X', 'Y', None]
    if mscd:
        dims = ['d' + x if x is not None else None for x in dims]
    all_vel = df \
        .groupby(traj_group) \
        .apply(pos_to_all_vel, xcol=dims[0], ycol=dims[1], zcol=dims[2],
               framecol='t', **kwargs)
    # output of pos_to_all_vel always has these cols by default
    vdims = ['vx', 'vy', 'vz']
    vdims = [v for v in vdims if v in all_vel.columns]
    # accumulate vx^2 + vy^2 + vz^2 as directed by dims
    absv2 = np.zeros_like(all_vel[vdims[0]])  # tf col always exists
    for vdim in vdims:
        absv2 += np.power(all_vel[vdim], 2)
    all_vel['abs(v)'] = np.sqrt(absv2)
    if vel_file:
        all_vel.to_csv(vel_file)
    if 'delta' not in groups:
        groups = groups + ['delta']
    msds = all_vel.groupby(groups)['abs(v)'].agg(['mean', 'std', 'count'])
    msds['ste'] = msds['std']/np.sqrt(msds['count']-1)
    msds['ste_norm'] = msds['ste']*np.sqrt(2/(msds['count']-1)) \
        * scipy.special.gamma(msds['count']/2) \
        / scipy.special.gamma((msds['count']-1)/2)
    return msds


preset_msd_args_ = """{
            'dvel': {'df': burgess.df_flat, 'mscd': True},
            'dvel_unp': {
                'df': burgess.df_flat[burgess.df_flat['foci'] == 'unp'],
                'mscd': True
            },
            'dvel_unp_by_wait': {
                'df': burgess.df_flat[burgess.df_flat['foci'] == 'unp'],
                'mscd': True,
                'traj_group': cell_cols + ['wait_id']
            },
            'dvel_unp_by_wait_na': {
                'df': burgess.df_flat[burgess.df_flat['foci'] == 'unp'],
                'mscd': True, 'traj_group': cell_cols + ['wait_id', 'na_id']
            },
            'vel': {'df': burgess.df.xs(1, level='spot'), 'mscd': False},
            'vel_double_counted': {'df': burgess.df, 'mscd': False},
            'vel_spot2': {'df': burgess.df.xs(2, level='spot'), 'mscd': False},
            'vel_unp_by_wait': {'df': burgess.df[burgess.df['foci'] == 'unp'],
                                'traj_group': cell_cols + ['spot', 'wait_id'],
                                'mscd': False},
            'vel_pair_by_wait': {'df': burgess.df[burgess.df['foci'] == 'pair'],
                                'traj_group': cell_cols + ['spot', 'wait_id'],
                                'mscd': False},
        }
"""
preset_msd_args = eval(preset_msd_args_)


def precompute_msds(prefix=burgess_dir, force_redo=False, **kwargs):
    """
    Precomputes a bunch of different "MS(C)Ds".

    Parameters
    ----------
    prefix : path_like
        Folder name in with to save output (default is burgess_dir)
    force_redo : bool
        Whether or not to redo the calculation if a file with the requested
        name already exists.
    **kwargs
        Forwarded to msds.msd

    Notes
    -----
    The default is to compute msds using the following 'name': argument pairs:

    .. code::

    """
    for name, msd_args in preset_msd_args.items():
        msd_args['vel_file'] = prefix / Path(name + '.csv')
        msd_file = prefix / Path('msds_' + name + '.csv')
        if not Path(msd_file).exists() or force_redo:
            msds = msd(**msd_args, **kwargs)
            msds.to_csv(msd_file)


precompute_msds.__doc__ += '    ' + preset_msd_args_
