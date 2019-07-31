"""Functions that outline the actual workflow of the paper."""

from . import *
from .msds import *

def make_all_intermediates(prefix=burgess_dir, force_redo=False):
    prefix = Path(prefix)

    # first write out all useful interpretations of the "velocity" of our
    # particle(s)
    precompute_msds(prefix, force_redo)

    # now compute the velocity correlation of the one that corresponds to the
    # most interpretable result (2D, vel of distance between particles,
    # starting a new "trajectory" every time the particles touch to eliminate
    # bias towards dX being zero)
    all_vel = pd.read_csv(burgess_dir / Path('dvel_by_wait.csv'))
    # this file is only needed to generate the cvv_stats, and should not be
    # loaded in, as it takes up huge amounts of space
    all_vvc_file = prefix / Path('all_vvc.csv')
    if not all_vvc_file.exists() or force_redo:
        vels_to_cvvs_by_hand(all_vel, cell_cols, all_vvc_file,
                             dzcol=None, max_t_over_delta=4)
    # file name to save actual velocity correlations
    cvv_stats_file = prefix / Path('cvv_stats.csv')
    if cvv_stats_file.exists() and not force_redo:
        cvv_stats = pd.read_csv(cvv_stats_file)
    else:
        cvv_stats = vvc_stats_by_hand(all_vvc_file, movie_cols)
        cvv_stats = cvv_by_hand_make_usable(cvv_stats, movie_cols)
        cvv_stats.to_csv(cvv_stats_file)
    cvv_fits_file = prefix / Path('cvv_fits.csv')
    if cvv_fits_file.exists() and not force_redo:
        cvv_fits = pd.read_csv(cvv_fits_file)
    else:
        # cvv_stats['t'] *= 30
        # cvv_stats['delta'] *= 30
        cvv_fits = cvv_stats.groupby(movie_cols).apply(get_best_fit_fixed_beta, p0=(1, 250), bounds=([0.1, 30], [2, 1500]))
        cvv_fits = cvv_fits.apply(pd.Series)
        cvv_fits.to_csv(cvv_fits_file)
    waitdf_file = prefix / Path('waitdf.csv')
    if waitdf_file.exists() and not force_redo:
        waitdf = pd.read_csv(waitdf_file)
    else:
        waitdf = df_flat.groupby(cell_cols).apply(discrete_trajectory_to_wait_times, t_col='t', state_col='foci')
        waitdf.to_csv(waitdf_file)

    return all_vel, cvv_fits, waitdf, msds_movies

