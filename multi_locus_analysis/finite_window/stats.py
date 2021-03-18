import warnings

import numpy as np
import pandas as pd
import scipy.optimize

from ..stats import ecdf


def average_lifetime(obs, traj_cols=['replicate']):
    """
    Estimate the true means of each state of the process.

    Doesn't require removing any censoring, because the process is stationary.
    """
    true_ends = obs[~np.isclose(obs['end_time'], obs['window_size'])]
    # the number of start times of each state are highly
    # correlated estimators of the total number of "start times"
    # of the full process that cycles. we combine them all to
    # average out edge effects, but don't gain much from it
    ave_num_ends = true_ends.groupby('state')['end_time'] \
        .count() \
        .mean()
    total_window = obs \
        .groupby(traj_cols)['window_size'].first() \
        .sum()
    total_mean_est = total_window / ave_num_ends
    # # this also works, but is obviously less robust
    # start_state = obs.groupby(traj_cols)['state'].first()
    # start_counts = start_state.value_counts()
    # mean_est = {
    #     name: total_mean_est*(count/np.sum(start_counts))
    #     for name, count in start_counts.items()
    # }
    total_state_times = obs.groupby('state')['wait_time'].sum()
    mean_est = {
        name: total_mean_est*(time/np.sum(total_state_times))
        for name, time in total_state_times.items()
    }
    return pd.Series(mean_est)


def window_sf(obs, traj_cols):
    window_sizes = obs.groupby(traj_cols)['window_size'].first().values
    # now sorted
    window_sizes, window_cdf = ecdf(window_sizes, pad_left_at_x=0)
    return 1 - window_cdf


def ecdf_combined(exterior, interior, window_sizes, ext_bins='auto', **kwargs):

    all_times, cdf_int, cdf_ext, Z_X, F_T = ecdf_ext_int(
        exterior, interior, window_sizes, **kwargs
    )

    y, t_bins = np.histogram(exterior, bins=ext_bins)
    bin_centers = (t_bins[:-1] + t_bins[1:]) / 2
    Z_hist = np.sum(y*np.diff(t_bins))
    # exterior estimator
    ext_est = (1 - y/Z_hist/Z_X)
    # standard error of mean, Normal approximation.
    ext_var = y / (Z_hist*Z_X)**2
    ext_weight = 1/ext_var
    # interior estimator
    int_est = np.interp(bin_centers, all_times, cdf_int*F_T)
    # Dvoretzky-Kiefer-Wolfowitz lower bound on the error at
    # confidence level alpha
    alpha = 0.32  # to match 1 sigma (std dev)
    N_int = len(interior)
    int_var = 1/(2*N_int) * np.log(2/alpha)
    int_weight = 1/int_var

    # final estimate is just the weighted average
    final_est = (
        int_est*int_weight + ext_est*ext_weight
    ) / (
        int_weight + ext_weight
    )
    return bin_centers, final_est


def _ccdf_int(x_int, cdf_int):
    if not np.isclose(cdf_int[0], 0):
        raise ValueError("Need to integrate interior CDF, but values start >0")
    ccdf_int = np.zeros_like(cdf_int)
    ccdf_int[1:] = np.cumsum(cdf_int[1:] * np.diff(x_int))
    return ccdf_int


def ecdf_ext_int(exterior, interior, window_sizes, window_sf=None,
                 times_allowed=None, pad_left_at_x=0):
    if len(exterior) < 5 or len(interior) < 5:
        raise ValueError("Need enough exterior *and* interior times to "
                         "perform correction.")

    x_ext, cdf_ext = ecdf(exterior, times_allowed, pad_left_at_x=pad_left_at_x)
    x_int, cdf_int = ecdf_windowed(
        interior, window_sizes, window_sf, times_allowed,
        pad_left_at_x=pad_left_at_x,
    )
    # integral of CDF w.r.t. t
    ccdf_int = _ccdf_int(x_int, cdf_int)

    if times_allowed is None:
        times_allowed = np.sort(np.concatenate((x_int, x_ext)))
        cdf_ext = np.interp(times_allowed, x_ext, cdf_ext)
        ccdf_int = np.interp(times_allowed, x_int, ccdf_int)
        cdf_int = np.interp(times_allowed, x_int, cdf_int)

    def err_f(ab, t, integrated_cdf, exterior_cdf):
        return np.linalg.norm(
            ab[0] * t + ab[1] * integrated_cdf - exterior_cdf
        )

    opt_out = scipy.optimize.minimize(
        err_f,
        x0=[1, -1],
        args=(times_allowed, ccdf_int, cdf_ext),
        bounds=((0, np.inf), (-np.inf, 0)),
    )
    if not opt_out.success:
        raise ValueError("Unable to compute F_X(T)!")
    Z_X, b = opt_out.x
    F_T = -b/Z_X
    return times_allowed, cdf_int, cdf_ext, Z_X, F_T


def ecdf_windowed(
        times, window_sizes, window_sf=None, times_allowed=None,
        auto_pad_left=None, pad_left_at_x=0, normalize=True,
        skip_times_allowed_check=False):
    """
    Empirical cumulative distribution for windowed observations.

    Parameters
    ----------
    times : (N,) array_like
        "Interior" waiting times.
    window_sizes : float or (N,) array_like
        The window size used. If a single value is passed, the window size is
        assumed to be constant.
    times_allowed : (M,) array_like
        Unique values that the data can take. Mostly useful for adding
        eCDF values at locations where data could or should have been observed
        but none was recorded (e.g. if a movie was taken with a given framerate
        but not all possible window lengths were observed).
    auto_pad_left : bool
        Deprecated. It makes more sense to default to left padding at zero for
        a renewal process. If left False, the data will not have a data value
        at the point where the eCDF equals zero. Use mean inter-data spacing to
        automatically generate an aesthetically reasonable such point. You
        *must* pass ``pad_left_at_x=False`` manually for this to work as
        expected.
    pad_left_at_x : float, default: 0
        Same as ``auto_pad_left``, but specify the point at which to add
        the leftmost point.
    window_sf : (M,) array_like of float
        For each unique window size in *window_sizes*, the number of
        trajectories with *at least* that window size. If not specified, it is
        assumed that each unique value of window size correponds to a unique
        trajectory. For the case of constant window size, this option is
        ignored.

    Returns
    -------
    x : (M,) array_like
        The values at which the eCDF was computed. By default
        :code:`np.sort(np.unique(y))`.
    cdf : (M,) array_like
        Values of the eCDF at each x.

    Notes
    -----
    If using ``times_allowed``, the *pad_left* parameters are redundant.
    """
    times = np.array(times)
    window_sizes = np.array(window_sizes)
    if times_allowed is not None and not skip_times_allowed_check:
        # again, not many resources in numpy for dealing with sorted arrays, so
        # we can't just take advantage of that and do it later...
        close_to_allowed = np.any(
            np.isclose(times[None, :], times_allowed[:, None]),
            axis=0  # one for each in times
        )
        if not np.all(close_to_allowed):
            # you'd have to modify the core algorithm to allow things to not
            # compare is_close
            raise ValueError("times_allowed passed does not contain all "
                             "elements of times passed (not np.is_close).")
    # unique doesn't take advantage of sorting anyway, so just do this here
    times_allowed = np.unique(times) \
        if times_allowed is None else np.array(times_allowed)
    times_allowed.sort()
    # flags only needed for multiple window sizes
    ignore_window_sf = False
    # allow providing single window size
    if window_sizes.size == 1:
        # and skip the extra pass where we calculate window_sf in that case
        ignore_window_sf = True
        window_sizes = window_sizes*np.ones_like(times)
    # algorithm is basically to sort then count, as in mla.stats.ecdf, but with
    # weights
    i = np.argsort(times)
    times = times[i]
    window_sizes = window_sizes[i]
    if auto_pad_left:
        dx = np.mean(np.diff(times_allowed))
        times_allowed = np.insert(times_allowed, 0, times_allowed[0] - dx)
    elif pad_left_at_x is not None:
        # remove this test for now, as weibull/pareto variables can generate
        # VERY small waits
        if not times_allowed[0] == pad_left_at_x:
            if times_allowed[0] < pad_left_at_x:
                warnings.warn('pad_left_at_x not left of x in ecdf_windowed! '
                              'Ignoring...')
            else:
                times_allowed = np.insert(times_allowed, 0, pad_left_at_x)

    # base weights due to interior censoring
    weights = (window_sizes - times)
    if not ignore_window_sf:
        if window_sf is None:
            # get fraction of windows that are at *least* of each width
            traj_window, window_ecdf = ecdf(window_sizes, pad_left_at_x=0)
            window_sf = 1 - window_ecdf
        else:
            # make sure to pad_left_at_x=0
            traj_window = np.insert(np.unique(window_sizes), 0, 0)
            if not np.isclose(window_sf[0], 1.0):
                raise ValueError("window_sf must be pad_left_at_x=0.")
        # for each observed time, we can get number of windows in which it can
        # have been observed
        # minus 1 because of how searchsorted returns indices
        window_i = np.searchsorted(traj_window, times) - 1
        # scale by fraction of trajectories where each time was observable
        weights = weights*window_sf[window_i]
    # TODO probably numba compile just this for loop
    i = 0
    num_obs = len(times)
    # weight contributed by each time. need to simply combine weights due to
    # "isclose" times to get final cdf output.
    full_cdf = np.cumsum(1/weights)
    # we loop over the allowed times instead of full_cdf itself so that we can
    # catch "empty" buckets, introduced either by times_allowed or pad_left*
    cdf = np.zeros(times_allowed.shape, dtype=np.dtype('float'))
    for cdf_i, t in enumerate(times_allowed):
        # if the times_allowed extend beyond any observations, just pad out the
        # result in a vectorized way
        if i >= num_obs - 1:
            cdf[cdf_i:] = full_cdf[-1]
            break
        # this triggers both if t is smaller than any observation (i == 0), or
        # if this particular cdf_i should have the same value as the previous
        # one (because there are no observations in this bucket)
        if times[i] > t and not np.isclose(times[i], t):
            if i == 0:
                # we're already zero-initialized
                continue
            # otherwise, just use the previous value, since there's guaranteed
            # to be at least one if we've already advanced i
            cdf[cdf_i] = cdf[cdf_i - 1]
        # otherwise, there is some amount of weight in this bucket. advance i
        # until we've counted all the weight
        while i + 1 < num_obs and np.isclose(times[i+1], t):
            i += 1
        cdf[cdf_i] = full_cdf[i]
    if normalize:
        cdf = cdf/cdf[-1]
    # some heavy-tailed dists can generate VERY small times, (isclose to 0)
    if pad_left_at_x == 0:
        cdf[0] = 0
    return times_allowed, cdf


def ecdf_simple(waits, T, pad_left_at_x=0):
    """cdf of interior times (ts > 0) observed in window of size T"""
    ts, counts = np.unique(waits, return_counts=True)
    i = np.argsort(ts)
    ts = ts[i]
    counts = counts[i]
    weights = T/(T - ts)
    return np.insert(ts, 0, pad_left_at_x), \
        np.insert(np.cumsum(counts*weights), 0, 0)/len(weights)
