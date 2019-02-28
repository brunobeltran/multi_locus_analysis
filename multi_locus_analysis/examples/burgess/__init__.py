r"""``burgess`` Data Set

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
from ... import pivot_loci

import pandas as pd
import numpy as np

from pathlib import Path

burgess_dir = Path(__file__).resolve().parent

cell_cols = ['locus', 'genotype', 'exp.rep', 'meiosis', 'cell']
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


