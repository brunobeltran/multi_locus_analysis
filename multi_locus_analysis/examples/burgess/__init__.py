r"""``burgess`` Data Set

Various movies (many cells per movie) of yeast cells undergoing meiosis. In
each cell, two loci are tagged (in the same color). Various mutants and stages
of meiosis were imaged.

Data interface
^^^^^^^^^^^^^^

``df``
    The data.

``cell_cols``
    The columsn to groupby to get each unique "cell" (i.e. each pair of
    trajectories :math:`X_1(t_k)` and :math:`X_2(t_k)`.

``traj_cols``
    The columns to groupby to get each trajectory (one particle at a time).

``frame_cols``
    The columns to groupby to get each frame taken (including both particles).

``spot_cols``
    The columns to groupby to get localization (one spot at one time).

Data columns
^^^^^^^^^^^^

``locus``
    a designator of which locus was tagged. ``HET5`` corresponds to a
    heterozygous cross of the ``URA3`` and ``LYS2`` tags.

``genotype``
    ``WT`` for wildtype or ``SP`` for :math:`\Delta`\ *spo11*.

``exp.rep``
    an unique integer for each experimental replicate (only unique if
    all other ``movie_cols`` are specified.

``meiosis``
    the stage of progression through meiosis. movies were taken by spotting
    cells onto a slide every thirty minutes. the times are labelled ``t#``,
    where the number nominally corresponds to the number of hours since the
    cells were transferred to sporulation media, but don't take it very
    seriously.

``cell``
    unique identifier for the different cells in a given movie.

``frame``
    frame counter for each movie

``t``
    number of seconds since beginning of movie. since only 1/30s frame
    rates were used, this is just 30 times the ``frame`` column.

``X``
    x-coordinate of a loci

``Y``
    y-coordinate of a loci

``Z``
    z-coordinate of a loci
"""
# required for importing data
from ...dataframes import pivot_loci
# for processing data
from ...stats import pos_to_all_vel, vels_to_cvvs_by_hand, vvc_stats_by_hand, cvv_by_hand_make_usable
from ...finite_window import discrete_trajectory_to_wait_times
from ...fitting import get_best_fit_fixed_beta

import pandas as pd
import numpy as np

from pathlib import Path

burgess_dir = Path(__file__).resolve().parent

movie_cols = ['locus', 'genotype', 'exp.rep', 'meiosis']
cell_cols = movie_cols + ['cell']
frame_cols = cell_cols + ['frame']
traj_cols = cell_cols + ['spot']
spot_cols = cell_cols + ['frame', 'spot']


df = pd.read_csv(burgess_dir / Path('xyz_conf_okaycells9exp.csv'))

def add_foci(df):
    foci1 = (np.isfinite(df.X1) & np.isfinite(df.Y1) & np.isfinite(df.Z1))
    foci2 = (np.isfinite(df.X2) & np.isfinite(df.Y2) & np.isfinite(df.Z2))
    notfoci2 = ~((np.isfinite(df.X2) | np.isfinite(df.Y2) | np.isfinite(df.Z2)))
    paired = foci1 & notfoci2
    unpaired = foci1 & foci2
    foci_col = df.observation.copy()
    foci_col[paired] = 'pair'
    foci_col[unpaired] = 'unp'
    foci_col[~(paired | unpaired)] = np.nan
    df['foci'] = foci_col
    return df

def replace_na(df):
    # apparently this doesn't work
    # df.loc[np.isnan(df['X2']), ['X2', 'Y2', 'Z2']]
    # so instead
    for i in ['X', 'Y', 'Z']:
        df.loc[np.isnan(df[i+'2']), i+'2'] = df.loc[np.isnan(df[i+'2']), i+'1']
    df.dropna(inplace=True)
    return df

# munge the raw data provided by Trent from the Burgess lab into the format our
# code expects
# df = df[df['observation'] == 'Okay'] # already done by trent for this file
df = add_foci(df)
del df['observation']
del df['desk']
cols = list(df.columns)
cols[5] = 'frame'
cols[6] = 't'
df.columns = cols
df = replace_na(df)
df.set_index(frame_cols, inplace=True)
df = pivot_loci(df, pivot_cols=['X', 'Y', 'Z'])


def make_all_intermediates(prefix=burgess_dir, force_redo=False):
    prefix = Path(prefix)

    df2 = pivot_loci(df, pivot_cols=['X', 'Y', 'Z'])
    for X in ['X', 'Y', 'Z']:
        df2['d'+X] = df2[X+'2'] - df2[X+'1']

    all_vel3_file = prefix / Path('all_vel3.csv')
    if all_vel3_file.exists() and not force_redo:
        all_vel3 = pd.read_csv(all_vel3_file)
    else:
        # happens instantaneously
        all_vel3 = df2 \
                .groupby(cell_cols) \
                .apply(pos_to_all_vel, xcol='dX', ycol='dY', zcol='dZ', framecol='frame')
        all_vel3['abs3(v)'] = np.sqrt(np.power(all_vel3['vx'], 2)
                                    + np.power(all_vel3['vy'], 2)
                                    + np.power(all_vel3['vz'], 2))
        all_vel3['abs2(v)'] = np.sqrt(np.power(all_vel3['vx'], 2)
                                    + np.power(all_vel3['vy'], 2))
        all_vel3.to_csv(all_vel3_file)
    all_vvc3_file = prefix / Path('all_vvc3.csv')
    # this file is only needed to generate teh cvv_stats, and shoul dnot be
    # loaded in, as it takes up huge amounts of space
    if not all_vvc3_file.exists() or force_redo:
        vels_to_cvvs_by_hand(all_vel3, cell_cols, all_vvc3_file,
                                 dzcol=None, max_t_over_delta=4)
    cvv_stats_file = prefix / Path('cvv_stats.csv')
    if cvv_stats_file.exists() and not force_redo:
        cvv_stats = pd.read_csv(cvv_stats_file)
    else:
        cvv_stats = vvc_stats_by_hand(all_vvc3_file, movie_cols)
        cvv_stats = cvv_by_hand_make_usable(cvv_stats, movie_cols)
        cvv_stats.to_csv(cvv_stats_file)
    cvv_fits_file = prefix / Path('cvv_fits.csv')
    if cvv_fits_file.exists() and not force_redo:
        cvv_fits = pd.read_csv(cvv_fits_file)
    else:
        # cvv_stats['t'] *= 30
        # cvv_stats['delta'] *= 30
        cvv_fits = cvv_stats.groupby(movie_cols).apply(get_best_fit_fixed_beta, bounds=([0.1, 1], [2, 173]))
        cvv_fits = cvv_fits.apply(pd.Series)
        cvv_fits.to_csv(cvv_fits_file)
    waitdf_file = prefix / Path('waitdf.csv')
    if waitdf_file.exists() and not force_redo:
        waitdf = pd.read_csv(waitdf_file)
    else:
        waitdf = df2.groupby(cell_cols).apply(discrete_trajectory_to_wait_times, t_col='t', state_col='foci')
        waitdf.to_csv(waitdf_file)
    msds_movies_file = prefix / Path('msds_movies.csv')
    if msds_movies_file.exists() and not force_redo:
        msds_movies = pd.read_csv(msds_movies_file)
    else:
        msds_movies = all_vel3.groupby(movie_cols + ['delta'])['abs3(v)'].agg(['mean', 'std', 'count'])
        msds_movies.to_csv(msds_movies_file)
    return all_vel3, cvv_fits, waitdf, msds_movies
