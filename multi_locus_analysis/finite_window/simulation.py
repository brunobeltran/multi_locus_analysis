from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd

import bruno_util
import bruno_util.random
from .. import stats as _mla_stats


@bruno_util.random.strong_default_seed
def ab_window_fast(rands, means, window_size, num_replicates=1, states=[0, 1],
                   seed=None, random_state=None):
    """WARNING: BUGGY! Needs to use t*f(t) for first time. Doesn't.

    Simulate a two-state system switching between states A and B.

    In addition to functions that can generate random waiting times for each
    state, this "fast" version of the code requires the
    average waiting times are are means[0], means[1], respectively.

    .. warning::

        apparently, :func:`bruno_util.random.strong_default_seed` is broken (or
        this function is) because passing a seed does not make the output
        reproducible.

    Parameters
    ----------
    rands : (2,) List[scipy.stats.rv_continuous]
        One of the random variables defined in :mod:`scipy.stats`.
        Alternatively, any callable that takes `random_state` and `size`
        kwargs. `random_state` should accept a :class:`np.random.RandomState`
        seed. `size` will be a tuple specifying output shape of random number
        array requested.
    means : (2,) array_like
        average waiting times for each of the states
    window_size : float
        the width of the window over which the observation takes place
    num_replicates : int
        number of times to run the simulation, default to 1
    states : (2,) array_like
        the "names" of each state, default to [0,1]
    seed : np.random.RandomState
        state to start the simulation with

    Returns
    -------
    df : pd.DataFrame
        The start/end times of each waiting time simulated. This data frame has
        `columns=['replicate', 'state', 'start_time', 'end_time',
        'window_start', 'window_end']`.

    Notes
    -----
    Consider the waiting time intersecting the left boundary of the
    observation window. The left boundary will be a uniform fraction of
    the way through this wait time. This can easily be seen in the case of
    finite-variance wait times using CLT and starting the switching process
    arbitrarily far left of the window of observation, or in general be imposed
    by requiring time-homogeneity of the experiment.

    We use this fact here to speed up correct simulation of time-homogenous
    windows by directly simulating only the waiting times within the windows
    instead of also simulating a long run of "pre-equilibrating" waiting times
    some offset before the window, as in :func:`ab_window`.
    """
    raise NotImplementedError("Currently buggy. Need to fix start time "
                              "generation, see docstring.")
    # np.concatenate can't handle concatenating nothing into nothing, so...
    if num_replicates <= 0:
        return pd.DataFrame(columns=['replicate', 'state', 'start_time',
                                     'end_time', 'window_start', 'window_end'])
    pool_size_guess = int(num_replicates*1.5*window_size/(means[0] + means[1]))
    pool_size_guess = max(pool_size_guess, 16)
    rands = [
        bruno_util.random.make_pool(rand, pool_size_guess,
                                    random_state=random_state)
        for rand in rands
    ]

    state_names = np.array(states)
    start_times = []
    end_times = []
    states = []
    for i in range(num_replicates):
        start_times.append([])
        end_times.append([])
        states.append([])

        prob_start_0 = means[0]/(means[0] + means[1])
        if np.random.random_sample() < prob_start_0:
            sim_state = 0
        else:
            sim_state = 1
        # fraction of way through left exterior-censored waiting time the left
        # window boundary lands
        u = np.random.random_sample()
        # full wait time, left exterior-censored
        r = rands[sim_state]()
        # true start
        t0 = -u*r
        # remainder of wait time not spend left of window
        t = (1 - u)*r

        start_times[-1].append(t0)
        states[-1].append(sim_state)
        sim_state = 1 - sim_state  # switch between 0 and 1
        while t <= window_size:
            end_times[-1].append(t)
            start_times[-1].append(t)
            states[-1].append(sim_state)
            t += rands[sim_state]()
            sim_state = 1 - sim_state  # switch between 0 and 1
        end_times[-1].append(t)
    start_times = np.concatenate([np.array(ts) for ts in start_times])
    end_times = np.concatenate([np.array(ts) for ts in end_times])
    replicate_ids = np.concatenate([i*np.ones_like(ts) for i, ts in enumerate(states)])
    states = np.concatenate([np.array(ss) for ss in states])
    states = state_names[states]
    df = pd.DataFrame.from_dict({'replicate': replicate_ids, 'state': states,
                                 'start_time': start_times, 'end_time': end_times})
    df['window_start'] = 0
    df['window_end'] = window_size
    return df


@bruno_util.random.strong_default_seed
def ab_window(rands, window_size, offset, num_replicates=1, states=[0, 1],
              seed=None, random_state=None):
    r"""
    Simulate an asynchronous two-state system from time 0 to `window_size`.

    Similar to :func:`multi_locus_analysis.finite_window.ab_window_fast`, but
    designed to work when the means of the distributions being used are hard to
    calculate.

    Simulate asynchronicity by starting the simulation in a uniformly random
    state at a time :math:`-t_\infty` (a large negative number).

    .. note::

        This number must be specified in the `offset` parameter and if it is
        not much larger than the means of the waiting times being used, the
        asynchronicity approximation will be very poor.

    The simulation only records between times 0 and window_size.

    Parameters
    ----------
    rands : (2,) List[scipy.stats.rv_continuous]
        Callable that takes "random_state" and "size" kwargs that accept
        a np.random.RandomState seed and a tuple specifying array sizes, resp.
    window_size : float
        The width of the window over which the observation takes place
    offset : float
        The (negative) time at which to start (in state 0) in order to
        equilibrate the simulation state by time t=0.
    states : (2,) array_like
        the "names" of each state, default to [0,1]
    num_replicates : int
        Number of times to run the simulation, default 1.
    seed : Optional[int]
        Random seed to start the simulation with
    random_state : np.random.RandomState
        Random state to start the simulation with. Preempts the seed argument.

    Returns
    -------
    df : pd.DataFrame
        The start/end times of each waiting time simulated. This data frame has
        columns=['replicate', 'state', 'start_time', 'end_time',
        'window_start', 'window_end'].
    """
    # np.concatenate can't handle concatenating nothing into nothing, so we
    if num_replicates <= 0:
        return pd.DataFrame(columns=['replicate', 'state', 'start_time',
                                     'end_time', 'window_start', 'window_end'])
    # accept relative or absolute offset
    if offset > 0:
        offset = -offset

    pool_size_guess = int(np.abs(offset)*num_replicates)
    rands = [bruno_util.random.make_pool(
        rand, pool_size_guess, random_state=random_state
    ) for rand in rands]
    # set aside state names to convert from 0,1 all at once later
    state_names = np.array(states)

    start_times = []
    end_times = []
    states = []
    for i in range(num_replicates):
        start_times.append([])
        end_times.append([])
        states.append([])
        t = offset
        # holds most recently used sim_state
        sim_state = 0
        r = rands[sim_state]()
        while t + r < 0:
            t += r
            sim_state = 1 - sim_state
            r = rands[sim_state]()
        states[i].append(sim_state)
        start_times[i].append(t)
        while t + r <= window_size:
            t += r
            end_times[i].append(t)
            sim_state = 1 - sim_state
            r = rands[sim_state]()
            states[i].append(sim_state)
            start_times[i].append(t)
        end_times[i].append(t + r)
    start_times = np.concatenate([np.array(ts) for ts in start_times])
    end_times = np.concatenate([np.array(ts) for ts in end_times])
    replicate_ids = np.concatenate([
        i*np.ones_like(ts) for i, ts in enumerate(states)
    ])
    states = np.concatenate([np.array(ss) for ss in states])
    states = state_names[states]
    df = pd.DataFrame.from_dict({
        'replicate': replicate_ids, 'state': states, 'start_time': start_times,
        'end_time': end_times
    })
    df['window_start'] = 0
    df['window_end'] = window_size
    return df

def _boot_final_est(N_var_T):
    import multi_locus_analysis.finite_window as fw
    import multi_locus_analysis.plotting.finite_window as fplt
    N_traj_per_boot, var_pair, err_t = N_var_T
    T = np.max(err_t)
    var_pair = {name: var for name, var in var_pair}
    sim = ab_window(
        [var.rvs for _, var in var_pair.items()],
        offset=-100*np.sum([var.mean() for _, var in var_pair.items()]),
        window_size=T,
        num_replicates=N_traj_per_boot,
        states=[name for name in var_pair]
    )
    obs = fw.sim_to_obs(sim)
    res = {}
    for name, var in var_pair.items():
        exterior = fplt._ext_from_obs(obs, name)
        interior, _ = fplt._int_win_from_obs(obs, name)
        bin_centers, final_est = fw.ecdf_combined(exterior, interior, T)
        res[name] = var.cdf(err_t) - np.interp(err_t, bin_centers, final_est)
    return res


def _boot_int_cdf(N_var_T):
    import multi_locus_analysis.finite_window as fw
    import multi_locus_analysis.plotting.finite_window as fplt
    N_traj_per_boot, var_pair, err_t = N_var_T
    T = np.max(err_t)
    var_pair = {name: var for name, var in var_pair}
    sim = ab_window(
        [var.rvs for _, var in var_pair.items()],
        offset=-100*np.sum([var.mean() for _, var in var_pair.items()]),
        window_size=T,
        num_replicates=N_traj_per_boot,
        states=[name for name in var_pair]
    )
    obs = fw.sim_to_obs(sim)
    res = {}
    for name, var in var_pair.items():
        exterior = fplt._ext_from_obs(obs, name)
        interior, _ = fplt._int_win_from_obs(obs, name)
        x, cdf = fw.ecdf_windowed(interior, T)
        res[name] = var.cdf(err_t)/var.cdf(T) - np.interp(err_t, x, cdf)
    return res


def bootstrap_int_error(var_pair, T, N_boot=1000, N_traj_per_boot=1000,
                        N_err_t=101):

    err_t = np.linspace(0, T, 101)
    derr_t = np.diff(err_t)
    err_ave = {var.name: np.zeros_like(err_t) for var in var_pair}
    err_std = {var.name: np.zeros_like(err_t) for var in var_pair}
    err_l2 = {var.name: 0 for var in var_pair}

    def accumulate_ave_std(res):
        accumulate_ave_std.N += 1
        for name in res:
            delta = res[name] - err_ave[name]
            err_ave[name] += delta/accumulate_ave_std.N
            err_std[name] += delta*(res[name] - err_ave[name])
            # separately, keep running average of l2 error
            res_mid = (res[name][1:]**2 + res[name][:-1]**2) / 2
            err_l2[name] += np.sum(np.sqrt(res_mid * derr_t))
    accumulate_ave_std.N = 0

    pickleable_var_pair = tuple((var.name, var.rv) for var in var_pair)

    with Pool(processes=cpu_count()) as pool:
        lazy_res = [
            pool.apply_async(
                _boot_int_cdf,
                ((N_traj_per_boot, pickleable_var_pair, err_t), )
            )
            for i in range(N_boot)
        ]
        for res in lazy_res:
            accumulate_ave_std(res.get())
    # for i in range(N_boot):
    #     accumulate_ave_std(_boot_int_i((N_traj_per_boot, pickleable_var_pair,
    #                                     err_t)))
    for name in err_std:
        err_std[name] /= accumulate_ave_std.N - 1

    return err_t, err_ave, err_std, err_l2
