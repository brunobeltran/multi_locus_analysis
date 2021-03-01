import numpy as np
import pandas as pd
import pytest
from scipy.stats import beta

from multi_locus_analysis import finite_window as fw


# for Beta distributions, making the window size a little less than 1
# guarantees there's not too many, and also not too few waits per
# trajectory, so that the effect we want to show is visible
example_window = 0.8
wait_vars = [beta(5, 2), beta(2, 2)]


def is_close_or_both_nan(x, y):
    return np.isclose(x, y) | (np.isnan(x) & np.isnan(y))


def is_equal_or_both_nan(x, y):
    return np.all((x == y) | (np.isnan(x) & np.isnan(y)))


def test_old_v_new_sim_to_frame_obs(N_traj=1000):
    # the old method doesn't work in the multi-window case. We've tested the
    # new method by eye, but no automated tests yet
    sim = fw.ab_window(
        [var.rvs for var in wait_vars],
        window_size=example_window,
        offset=-100*np.sum([var.mean() for var in wait_vars]),
        num_replicates=N_traj,
        states=['A', 'B'],
    )
    traj_cols = ['replicate']

    # first use old method for getting discrete frames, very involved
    T = sim.window_end.max()
    movie_frame_t = np.linspace(0, T, 21)
    movies = sim.groupby(traj_cols).apply(
        fw.traj_to_movie,
        times=movie_frame_t
    )
    movies = movies.T.unstack()
    movies.name = 'state'
    old_obs = pd.DataFrame(movies) \
        .reset_index() \
        .groupby(traj_cols) \
        .apply(fw.movie_to_waits)

    # now new method, much happier
    frames = fw.munging.simulation_to_frame_times(sim, movie_frame_t, traj_cols=traj_cols)
    new_obs = fw.sim_to_obs(frames.reset_index(), traj_cols=traj_cols)

    # now make the index match for easy comparison
    new_obs = new_obs.reset_index()
    new_obs['rank_order'] -= 1
    new_obs = new_obs.set_index(traj_cols + ['rank_order'])

    assert np.all(np.isclose(
        0, np.abs(new_obs['wait_time'] - old_obs['wait_time']) > 0
    ))


def test_old_v_new_obs_to_sim(N_traj=1000):
    # test the full multi-window case
    het_trajs = [
        fw.ab_window(
            [var.rvs for var in wait_vars],
            window_size=window,
            offset=-100*np.sum([var.mean() for var in wait_vars]),
            num_replicates=N_traj,
            states=['A', 'B'],
        )
        for window in np.array([1/2, 1, 2])*example_window
    ]
    sim = pd.concat(het_trajs, ignore_index=True)
    traj_cols = ['window_end', 'replicate']
    multi_T_waits = fw.sim_to_obs(sim, traj_cols=traj_cols)
    with pytest.deprecated_call():
        multi_T_waits_old = sim \
            .groupby(traj_cols) \
            .apply(fw.traj_to_waits)
    for col in traj_cols:
        del multi_T_waits_old[col]

    assert np.all(
        multi_T_waits.reset_index()['rank_order']
        == multi_T_waits_old.reset_index()['rank_order']
    )
    assert np.all(np.isclose(
        multi_T_waits['wait_time'], multi_T_waits_old['wait_time']
    ))
    assert np.all(multi_T_waits['wait_type'] == multi_T_waits_old['wait_type'])
    assert np.all(multi_T_waits['state'] == multi_T_waits_old['state'])


def test_traj_to_movie():
    trajs = np.array([
        [0, 1, -0.100,  0.100, 0, 1],  # easy case
        [0, 0,  0.100,  1.100, 0, 1],
        [1, 0, -0.100,  1.000, 0, 1],  # switch exactly at window_end
        [2, 0,  0.000,  0.101, 0, 1],
        [2, 1,  0.101,  0.109, 0, 1],  # should be a "missed" time
        [2, 0,  0.109,  0.200, 0, 1],
        [2, 1,  0.200,  1.100, 0, 1],
        [3, 1,  0.000,  0.199, 0, 1],
        [3, 0,  0.199,  0.200, 0, 1],  # several "aligned" in a row
        [3, 1,  0.200,  0.300, 0, 1],
        [3, 0,  0.300,  1.000, 0, 1],
        [4, 1,  0.000,  0.099, 0, 1],  # starts exactly at window_start
        [4, 0,  0.099,  1.200, 0, 1],
        [5, 1, -0.100,  1.100, 0, 1],  # no transitions, both cases
        [6, 0, -0.100,  1.100, 0, 1]])
    trajs = pd.DataFrame(trajs)
    trajs.columns = ['replicate', 'state', 'start_time', 'end_time',
                     'window_start', 'window_end']

    # first ignore end_time
    movies = trajs.groupby('replicate').apply(
        fw.traj_to_movie, times=np.linspace(0, 1, 11))
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
    movies = trajs.groupby('replicate').apply(
        fw.traj_to_movie,
        times=np.linspace(0, 1, 11), end_times_col='end_time'
    )
    expected_movies = np.array([
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])

    assert(is_equal_or_both_nan(expected_movies, movies.values))


def test_movie_to_waits():
    movies = pd.DataFrame([], columns=['replicate', 't', 'state'])
    movies['replicate'] = np.concatenate([
        np.ones(10), 2*np.ones(10), 3*np.ones(10)
    ]).astype(int)
    movies['t'] = np.concatenate(3*[np.arange(10)]).astype(float)
    movies['state'] = 0
    movies.loc[0, 'state'] = 1
    movies.loc[9, 'state'] = 1
    movies.loc[20:30, 'state'] = 1

    waits = movies.groupby('replicate').apply(fw.movie_to_waits, t_col='t',
                                              state_col='state')
    expected_waits = np.array([
        [1., 0., 0., 1., 1., 1., 0., 1., 'left exterior', 9.],
        [1., 1., 1., 9., 8., 0., 7., 9., 'interior', 9.],
        [1., 2., 9., 9., 0., 1., 0., 1., 'right exterior', 9.],
        [2., 0., 0., 9., 9., 0., 9., 9., 'full exterior', 9.],
        [3., 0., 0., 9., 9., 1., 9., 9., 'full exterior', 9.]
    ], dtype=object)

    assert(np.all(waits.reset_index().values == expected_waits))
