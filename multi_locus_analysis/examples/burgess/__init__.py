r"""``burgess`` Data Set

Various movies (many cells per movie) of yeast cells undergoing meiosis. In
each cell, two loci are tagged (in the same color). Various mutants and stages
of meiosis were imaged.

Data interface
^^^^^^^^^^^^^^

``df``
    The data, one row per (locus, frame).

``df_flat``
    The data, one row per frame, both loci.

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

import pandas as pd
import numpy as np

from pathlib import Path

# known biological constants, in bp
bp_per_nuc = 147  # 146 or 147, depending on the crystal
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2493396/
cerevisiae_nrl = 165  # nucleosome repeat length in bp
cerevisiae_linker = cerevisiae_nrl - bp_per_nuc
chrii_size_bp = 813184
chrv_size_bp = 576874
um_per_bp = 0.00034

# measurement parameters
t_data = np.arange(0, 1501, 30).astype(float)

location_lys_bp = np.mean([469748, 473926])
location_cen2_bp = np.mean([238207, 238323])
lys_locus_frac = location_lys_bp / chrii_size_bp
chrii_centromere_frac = location_cen2_bp / chrii_size_bp

location_ura_bp = np.mean([116167, 116970])
location_cen5_bp = np.mean([151987, 152104])
ura_locus_frac = location_ura_bp/chrv_size_bp
chrv_centromere_frac = location_cen5_bp/chrv_size_bp

# bare DNA wlc parameters (found to not work well)
kuhn_length_wlc = 0.05  # 50nm
chrv_size_wlc_um = 0.34 * chrv_size_bp / 1000  # um
chrv_size_wlc_kuhn = chrv_size_wlc_um / kuhn_length_wlc
location_ura_wlc_um = location_ura_bp*(
    chrv_size_wlc_um / chrv_size_bp
)
location_cen5_wlc_um = location_cen5_bp*(
    chrv_size_wlc_um / chrv_size_bp
)

# effective nucleosome chain parameters (found to work much better)
if cerevisiae_linker == 18:
    kuhn_length_nuc_chain = 0.01776  # um, from Beltran et al, PRL 2019
chrv_num_nuc = chrv_size_bp / cerevisiae_nrl
chrv_total_linker = chrv_num_nuc * cerevisiae_linker
chrv_linker_um = chrv_total_linker * um_per_bp
chrv_size_nuc_chain_kuhn = chrv_linker_um / kuhn_length_nuc_chain
location_ura_nuc_chain_um = location_ura_bp*(
    chrv_linker_um / chrv_size_bp
)
location_cen5_nuc_chain_um = location_cen5_bp*(
    chrv_linker_um / chrv_size_bp
)

chrii_num_nuc = chrii_size_bp / cerevisiae_nrl
chrii_total_linker = chrii_num_nuc * cerevisiae_linker
chrii_linker_um = chrii_total_linker * um_per_bp
chrii_size_nuc_chain_kuhn = chrii_linker_um / kuhn_length_nuc_chain
location_lys_nuc_chain_um = location_lys_bp*(
    chrii_linker_um / chrii_size_bp
)
location_cen2_nuc_chain_um = location_cen2_bp*(
    chrii_linker_um / chrii_size_bp
)

# backcompat
chrii_size_nuc_chain_um = chrii_linker_um
chrv_size_nuc_chain_um = chrv_linker_um

# derived parameters
# Average of het5 msd convex hull distribution. Also matches well the
# higher of the two wt, ura3, t3 "modes" of the mscd plateau distribution
nuc_radius_um = 1.3
sim_nuc_radius_um = 1
old_sim_D = 20  # um^2/s, old value used for existing sims, 10/2020
sim_D = 3.4872926794253383  # see discussion in "determining-diffusivity" docs

# finally, the actual data
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
    notfoci2 = ~(np.isfinite(df.X2) | np.isfinite(df.Y2) | np.isfinite(df.Z2))
    paired = foci1 & notfoci2
    unpaired = foci1 & foci2
    foci_col = df.observation.copy()
    foci_col[paired] = 'pair'
    foci_col[unpaired] = 'unp'
    foci_col[~(paired | unpaired)] = np.nan
    df['foci'] = foci_col
    return df


def pixels_to_units(df):
    """
    Keep track of necessary conversions into real units.

    For now all data was collected with Z-stacks where pixels in the X,Y plane
    have a "real" width of 0.13333 um/pixels, whereal the z-stacks are spaced
    at 0.25um intervals. The df_xyz file has units of "pixels/10"

    In the future, Trent will want to add code here to make sure that other
    movies with different pixel sizes can be compared directly to old
    experiments.
    """
    x_pix_to_um = 0.13333
    y_pix_to_um = 0.13333
    z_pix_to_um = 0.25
    df['X1'] *= 10*x_pix_to_um
    df['Y1'] *= 10*y_pix_to_um
    df['Z1'] *= 10*z_pix_to_um
    df['X2'] *= 10*x_pix_to_um
    df['Y2'] *= 10*y_pix_to_um
    df['Z2'] *= 10*z_pix_to_um


def replace_na(df):
    """
    Undo Trent's decision to make pair times have one locus as X=NaN.

    Assuming we have add_foci'd, we don't *need* to artificially set the
    second trajectory's values to NaN, so undo that here.
    """
    # apparently this doesn't work
    # df.loc[np.isnan(df['X2']), ['X2', 'Y2', 'Z2']]
    # so instead
    for i in ['X', 'Y', 'Z']:
        df.loc[np.isnan(df[i+'2']), i+'2'] = df.loc[np.isnan(df[i+'2']), i+'1']
    return df


def breakup_by_na(traj):
    """
    Split up a trajectory each time it contains a NaN.

    Takes a Burgess trajectory, and create a new column "na_id" that
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
    """
    Split a trajectory each time a colocalization event occurs.

    Takes a single Burguess "cell", and add a column that uniquely tracks the
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
    """
    Munge raw data provided by Trent from the Burgess lab into desired format.

    Trent already does the equivalent of df = df[df['observation'] == 'Okay'],
    by manually classifying trajectories, largely...
    """
    df = add_foci(df)
    del df['observation']
    del df['desk']
    cols = list(df.columns)
    cols[5] = 'frame'
    cols[6] = 't'
    df.columns = cols
    del cols
    pixels_to_units(df)
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
    df.sort_index(inplace=True)
    df.to_csv(df_file)
    df_flat.sort_index(inplace=True)
    df_flat.to_csv(df_flat_file)
else:
    df = pd.read_csv(df_file)
    df.set_index(frame_cols + ['spot'], inplace=True)
    df_flat = pd.read_csv(df_flat_file)
    df_flat.set_index(frame_cols, inplace=True)


# API
from . import analysis  # NOQA
from . import msds  # NOQA
from . import plotting  # NOQA
from . import simulation  # NOQA
from . import workflow  # NOQA
