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

condition_cols = ['locus', 'genotype', 'meiosis']
movie_cols = ['locus', 'genotype', 'exp.rep', 'meiosis']
cell_cols = movie_cols + ['cell']

frame_cols = cell_cols + ['frame']
traj_cols = cell_cols + ['spot']
spot_cols = cell_cols + ['frame', 'spot']

df_xyz = pd.read_csv(burgess_dir / Path('xyz_conf_okaycells9exp.csv'))

def add_foci(df):
    """Extract a column labeling whether the loci are paired at each frame.

    Data comes in as raw trajectories. "Paired" trajectories (where the
    second loci cannot be measured because it is coincident with the first) are
    labeled only by the fact that the second loci has NaN values.

    NOTE: when both X1 and X2 are NaN, this is simply a bad frame, these should
    remain NaN.
    """
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
    """Assuming we have add_foci'd, we don't *need* to artificially set the
    second trajectory's values to NaN, so undo that here."""
    # apparently this doesn't work
    # df.loc[np.isnan(df['X2']), ['X2', 'Y2', 'Z2']]
    # so instead
    for i in ['X', 'Y', 'Z']:
        df.loc[np.isnan(df[i+'2']), i+'2'] = df.loc[np.isnan(df[i+'2']), i+'1']
    return df

def breakup_by_na(traj):
    """Take a Burgess trajectory, and create a new column "na_id" that
    uniquely tracks continuous chunks of trajectory where there are no
    NAN's.

    This should be applied to the "flat" dataframe, without breaking up the
    trajectories by which spot it is. This is because if one spot was
    incorrectly measured, we have effectively zero confidence that the other
    one was.
    """
    # first make sure that the following code makes sense. we use the
    # assumption that all frames are included in the data, with frames that
    # failed to yield a good measurement being labeled simply by a row of all
    # NaN in the data. So we first check that no "frames" are "left out"
    assert(np.all(np.diff(traj.reset_index()['frame'].values) == 1))
    # now find the NaN rows, if they exist. the following array will be
    # true if we should break to a "new" na_id at that row (some redundancy,
    # since strings of consecutive NAN's are "broken" at each row)
    break_on = np.any(traj.isna(), axis=1)
    # unique ID
    traj['na_id'] = np.cumsum(break_on)
    return traj

def add_wait_id(traj):
    """Take a single Burguess "cell", and add a column that uniquely tracks the
    individual stretches of time over which that cell has "pair" or "unp" loci.

    Assumes that wait_id is not changed by internal NaN's. If you wish to break
    over both wait_id and NaN, simply use both this function and breakup_by_na,
    and groupby both _id's simultaneously.
    """
    fnum = (traj['foci'] == 'pair').astype(int)
    # 0 or 1 depending on whether we should start a new wait_id at that row
    break_on = np.insert(np.abs(np.diff(fnum.values)), 0, 0)
    traj['wait_id'] = np.cumsum(break_on)
    return traj


def munge_data(df):
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
    del cols
    df = replace_na(df)
    df.set_index(frame_cols, inplace=True)
    df = df.groupby(cell_cols).apply(breakup_by_na)
    df = df.groupby(cell_cols).apply(add_wait_id)
    df_flat = df
    df = pivot_loci(df, pivot_cols=['X', 'Y', 'Z'])
    for X in ['X', 'Y', 'Z']:
        df_flat['d'+X] = df_flat[X+'2'] - df_flat[X+'1']
    return df, df_flat

df_file = burgess_dir / Path('df.csv')
df_flat_file = burgess_dir / Path('df_flat.csv')
if not (df_file.exists() and df_flat_file.exists()):
    df, df_flat = munge_data(df_xyz)
    df.to_csv(df_file)
    df_flat.to_csv(df_flat_file)
else:
    df = pd.read_csv(df_file)
    df_flat = pd.read_csv(df_flat_file)

