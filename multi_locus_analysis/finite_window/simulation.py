from multiprocessing import Pool, cpu_count

import lifelines
import numpy as np
import pandas as pd
import scipy.stats
from scipy.stats import expon

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
    l2_ave = {var.name: 0 for var in var_pair}
    l2_std = {var.name: 0 for var in var_pair}
    linf_ave = {var.name: 0 for var in var_pair}
    linf_std = {var.name: 0 for var in var_pair}

    def accumulate_ave_std(res):
        accumulate_ave_std.N += 1
        for name in res:
            delta = res[name] - err_ave[name]
            err_ave[name] += delta/accumulate_ave_std.N
            err_std[name] += delta*(res[name] - err_ave[name])
            # separately, keep running average of l2 error
            res_mid = (res[name][1:]**2 + res[name][:-1]**2) / 2
            l2 = np.sum(np.sqrt(res_mid * derr_t))
            delta = l2 - l2_ave[name]
            l2_ave[name] += delta/accumulate_ave_std.N
            l2_std[name] += delta*(l2 - l2_ave[name])
            # separately, keep running average of l^\inf error
            linf = np.max(np.abs(res[name]))
            delta = linf - linf_ave[name]
            linf_ave[name] += delta/accumulate_ave_std.N
            linf_std[name] += delta*(linf - linf_ave[name])

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

    return err_t, err_ave, err_std, l2_ave, l2_std, linf_ave, linf_std


def _mean_from_exp_cdf(x, cdf):
    y = np.log(1 - cdf)
    i = np.isfinite(x) & np.isfinite(y)
    return -1/scipy.stats.linregress(x[i], y[i]).slope


def _example_lambda_fit(V_T_N):
    import multi_locus_analysis.finite_window as fw
    import multi_locus_analysis.plotting.finite_window as fplt
    lambdas, T, N_traj = V_T_N
    var_pair = [
        fplt.Variable(expon(scale=lam), name=f"Exp({lam})")
        for lam in lambdas
    ]
    sim = fw.ab_window(
        [var.rvs for var in var_pair],
        offset=-100*np.sum([var.mean() for var in var_pair]),
        window_size=T,
        num_replicates=N_traj,
        states=[var.name for var in var_pair]
    )
    obs = fw.sim_to_obs(sim)

    mean_est = fw.average_lifetime(obs)
    true_mean = {var.name: var.mean() for var in var_pair}
    naive_slope_est = {}
    correct_slope_est = {}
    kaplan_slope_est = {}
    uncensored_baseline = {}
    for var in var_pair:
        # naive
        interior, windows = fplt._int_win_from_obs(obs, var.name)
        try:
            x_int, cdf_int = fw.ecdf_windowed(interior, windows)
            naive_slope_est[var.name] = _mean_from_exp_cdf(x_int, cdf_int)
        except:
            naive_slope_est[var.name] = np.nan
        # corrected
        exterior = fplt._ext_from_obs(obs, var.name)
        try:
            bin_centers, final_cdf = fw.ecdf_combined(exterior, interior, T)
            correct_slope_est[var.name] = _mean_from_exp_cdf(bin_centers, final_cdf)
        except:
            correct_slope_est[var.name] = np.nan
        # kaplan
        times = np.concatenate([interior, exterior])
        is_interior = np.concatenate(
            [np.ones_like(interior), np.zeros_like(exterior)]
        ).astype(bool)
        try:
            kmf = lifelines.KaplanMeierFitter() \
                    .fit(times, event_observed=is_interior)
            x_kap = kmf.cumulative_density_.index.values
            cdf_kap = kmf.cumulative_density_.values.flatten()
            kaplan_slope_est[var.name] = _mean_from_exp_cdf(x_kap, cdf_kap)
        except:
            kaplan_slope_est[var.name] = np.nan
        # uncensored baseline
        num_obs = len(interior)
        try:
            x_unc, cdf_unc = _mla_stats.ecdf(var.rvs(size=(num_obs,)),
                                            pad_left_at_x=0)
            uncensored_baseline[var.name] = _mean_from_exp_cdf(x_unc, cdf_unc)
        except:
            uncensored_baseline[var.name] = np.nan
    df = pd.concat(map(pd.Series,
        [true_mean, correct_slope_est, naive_slope_est,
        mean_est, kaplan_slope_est, uncensored_baseline]
    ), axis=1)
    df.columns = ['true', 'corrected', 'naive', 'count-based', 'kaplan',
                  'uncensored']
    return df


# def _choose_N_traj(means, min_obs=100, T=1):
#     # if mean is *too* large, then we have to beware of seeing enough total
#     # exterior times
#     exterior_per_traj = 1/np.sum(means) * T
#     required_for_ext = min_obs / exterior_per_traj
#     interior_per_traj =

def bootstrap_exp_fit_error(N_boot=1000, N_traj=1000, N_lam=41, T=1):
    lambdas_uniq = np.logspace(-2, 2, N_lam)
    lambdas = np.random.choice(lambdas_uniq, size=(N_boot,2))
    Ts = T*np.ones((N_boot,))
    N_trajs = (N_traj*np.ones((N_boot,))).astype(int)
    # # another approach might be to just sample less the ridiculous ones
    # N_boot = np.linspace(min_boot, max_boot, N_lam).astype(int)
    # N_tot = np.sum(N_boot)
    # lambdas = np.zeros((N_tot,))
    # j = 0
    # for i, N in enumerate(N_boot):
    #     lambdas[j:j+N] = lambdas_uniq[i]
    #     j += N
    V_T_Ns = list(zip(lambdas, Ts, N_trajs))
    with Pool(processes=cpu_count()) as pool:
        df = pd.concat(pool.map(_example_lambda_fit, V_T_Ns))
    df['N_traj'] = N_traj
    return df


def _alpha_from_cdf(x, cdf, xmin):
    xlog = np.log10(x)
    ylog = np.log10(1 - cdf)
    i = (x > xmin) & np.isfinite(xlog) & np.isfinite(ylog)
    # -slope, and +1 because cdf not pdf
    return 1 - scipy.stats.linregress(xlog[i], ylog[i]).slope


def _example_pareto_alpha(V_T_N):
    import multi_locus_analysis.finite_window as fw
    import multi_locus_analysis.plotting.finite_window as fplt

    # unpack parameters first
    (betas, xmin), T, N_traj = V_T_N
    var_pair = [
        fplt.Variable(scipy.stats.pareto(beta, scale=xmin),
                      name=f'Pareto({beta:0.3g})')
        for beta in betas
    ]
    # run one simulation
    sim = fw.ab_window(
        [var.rvs for var in var_pair],
        offset=-100*np.sum([var.mean() for var in var_pair]),
        window_size=T,
        num_replicates=N_traj,
        states=[var.name for var in var_pair]
    )
    obs = fw.sim_to_obs(sim)

    # now extract alpha several different ways
    true_alpha = {var.name: var.args[0] + 1 for var in var_pair}
    mle_interior_est = {}
    mle_uncensored_baseline = {}
    fit_interior = {}
    fit_corrected = {}
    fit_kaplan = {}
    fit_uncensored_baseline = {}
    for var in var_pair:
        # mle, interior
        try:
            interior, windows = fplt._int_win_from_obs(obs, var.name)
            num_obs = len(interior)
            mle_interior_est[var.name] = _mla_stats.power_law_slope_mle(
                interior, xmin, num_obs)
        except:
            mle_interior_est[var.name] = np.nan
        # fit, interior
        try:
            x_int, cdf_int = fw.ecdf_windowed(interior, windows)
            fit_interior[var.name] = _alpha_from_cdf(x_int, cdf_int, xmin)
        except:
            fit_interior[var.name] = np.nan
        # fit, corrected
        try:
            exterior = fplt._ext_from_obs(obs, var.name)
            bin_centers, final_cdf = fw.ecdf_combined(exterior, interior, T)
            fit_corrected[var.name] = _alpha_from_cdf(bin_centers, final_cdf, xmin)
        except:
            fit_corrected[var.name] = np.nan
        # fit, kaplan
        try:
            times = np.concatenate([interior, exterior])
            is_interior = np.concatenate(
                [np.ones_like(interior), np.zeros_like(exterior)]
            ).astype(bool)
            kmf = lifelines.KaplanMeierFitter() \
                .fit(times, event_observed=is_interior)
            x_kap = kmf.cumulative_density_.index.values
            cdf_kap = kmf.cumulative_density_.values.flatten()
            fit_kaplan[var.name] = _alpha_from_cdf(x_kap, cdf_kap, xmin)
        except:
            fit_kaplan[var.name] = np.nan
        # mle, uncensored baseline
        try:
            uncensored_obs = var.rvs(size=(num_obs,))
            mle_uncensored_baseline[var.name] = _mla_stats.power_law_slope_mle(
                uncensored_obs, xmin, num_obs)
        except:
            mle_uncensored_baseline[var.name] = np.nan
        # fit, uncensored baseline
        try:
            x_unc, cdf_unc = _mla_stats.ecdf(uncensored_obs, pad_left_at_x=0)
            fit_uncensored_baseline[var.name] = \
                _alpha_from_cdf(x_unc, cdf_unc, xmin)
        except:
            fit_uncensored_baseline[var.name] = np.nan
    df = pd.concat(
        map(
            pd.Series,
            [true_alpha, mle_interior_est, mle_uncensored_baseline,
            fit_interior, fit_corrected, fit_kaplan,
            fit_uncensored_baseline]
        ), axis=1
    )
    df.columns = ['true', 'mle-interior', 'mle-uncensored', 'fit-interior',
                 'fit-corrected', 'fit-kaplan', 'fit-uncensored']
    return df


def bootstrap_alpha_fit_error(N_boot=1000, N_traj=1000, xmin=0.1, T=1):
    betas_uniq = np.linspace(1, 2, 21)[1:]
    # # incantation to repeat each betas_uniq N_boot times
    # betas = np.tile(betas_uniq, (N_boot, 1)).T.flatten()
    betas = np.random.choice(betas_uniq, size=(N_boot,2))
    Ts = T*np.ones((N_boot,))
    N_trajs = (N_traj*np.ones_like(Ts)).astype(int)
    params = zip(betas, xmin*np.ones_like(Ts))
    V_T_Ns = list(zip(params, Ts, N_trajs))
    with Pool(processes=cpu_count()) as pool:
        df = pd.concat(pool.map(_example_pareto_alpha, V_T_Ns))
    df['N_traj'] = N_traj
    return df
