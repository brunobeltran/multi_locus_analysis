from multi_locus_analysis import finite_window as fw
import pandas as pd
import numpy as np


def is_equal_or_both_nan(x, y):
    return np.all((x == y) | (np.isnan(x) & np.isnan(y)))

def test_traj_to_movie():
    trajs = np.array([
        [0, 1, -0.100,  0.100, 0, 1], # easy case
        [0, 0,  0.100,  1.100, 0, 1],
        [1, 0, -0.100,  1.000, 0, 1], # switch exactly at window_end
        [2, 0,  0.000,  0.101, 0, 1],
        [2, 1,  0.101,  0.109, 0, 1], # should be a "missed" time
        [2, 0,  0.109,  0.200, 0, 1],
        [2, 1,  0.200,  1.100, 0, 1],
        [3, 1,  0.000,  0.199, 0, 1],
        [3, 0,  0.199,  0.200, 0, 1], # several "aligned" in a row
        [3, 1,  0.200,  0.300, 0, 1],
        [3, 0,  0.300,  1.000, 0, 1],
        [4, 1,  0.000,  0.099, 0, 1], # starts exactly at window_start
        [4, 0,  0.099,  1.200, 0, 1],
        [5, 1, -0.100,  1.100, 0, 1], # no transitions, both cases
        [6, 0, -0.100,  1.100, 0, 1]])
    trajs = pd.DataFrame(trajs)
    trajs.columns = ['replicate', 'state', 'start_time', 'end_time',
                     'window_start', 'window_end']

    # first ignore end_time
    movies = trajs.groupby('replicate').apply(fw.traj_to_movie, times=np.linspace(0, 1, 11))
    expected_movies = np.array([
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])

    assert(np.all(expected_movies == movies.values))

    # now including end_time
    movies = trajs.groupby('replicate').apply(fw.traj_to_movie,
        times=np.linspace(0, 1, 11), end_times_col='end_time')
    expected_movies = np.array([
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, np.nan],
        [0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, np.nan],
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])

    assert(is_equal_or_both_nan(expected_movies, movies.values))

