import warnings

import numpy as np
import scipy.optimize

from ..stats import ecdf


def ecdf_ext_int(exterior, interior, window_sizes, window_sf=None,
                 times_allowed=None, pad_left_at_x=0):
    x_ext, cdf_ext = ecdf(exterior, times_allowed, pad_left_at_x=pad_left_at_x)
    x_int, cdf_int = ecdf_windowed(
        interior, window_sizes, times_allowed,
        pad_left_at_x=pad_left_at_x,
        window_sf=window_sf
    )
    # now compute integral of CDF w.r.t. t
    ccdf_int = np.zeros_like(cdf_int)
    ccdf_int[1:] = np.cumsum(cdf_int[1:] * np.diff(x_int))

    if times_allowed is None:
        times_allowed = np.sort(np.concatenate((x_int, x_ext)))
        cdf_ext = np.interp(times_allowed, x_ext, cdf_ext)
        ccdf_int = np.interp(times_allowed, x_int, ccdf_int)

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
    Z_X, F_T = opt_out.x

    return times_allowed, cdf_int*(-F_T/Z_X)


def ecdf_windowed(
        times, window_sizes, window_sf=None, times_allowed=None,
        auto_pad_left=None, pad_left_at_x=0, normalize=True):
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
    # set up y (values to hist) ymax (windows) and x (bins)
    # well, the eCDF equivalents, not actually a hist...
    y = times
    ymax = window_sizes
    y = np.array(y)
    ymax = np.array(ymax)
    # variables only needed for multiple window sizes
    ignore_window_sf = False
    uniq_ymax = None
    # allow providing single window size
    if ymax.size == 1:
        ignore_window_sf = True
        ymax = ymax*np.ones_like(y)
    i = np.argsort(y)
    y = y[i]
    ymax = ymax[i]
    if times_allowed is not None:
        x = np.unique(times_allowed)
    else:
        x = np.unique(y)
    x.sort()
    if auto_pad_left:
        dx = np.mean(np.diff(x))
        x = np.insert(x, 0, x[0] - dx)
    elif pad_left_at_x is not None:
        if x[0] <= pad_left_at_x:
            warnings.warn('pad_left_at_x not left of x in ecdf_windowed! '
                          'Ignoring...')
        else:
            x = np.insert(x, 0, pad_left_at_x)

    num_obs = len(y)
    cdf = np.zeros(x.shape, dtype=np.dtype('float'))
    if not ignore_window_sf and window_sf is None:
        # get fraction of windows that are at *least* of each width
        uniq_ymax, window_sf = ecdf(ymax, pad_left_at_x=0)
    weights = (ymax - y)
    if not ignore_window_sf:
        # for each observed time, we can get number of windows in which it can
        # have been observed
        if uniq_ymax is None:
            # don't forget to pad_left_at_x=0
            uniq_ymax = np.insert(np.unique(ymax), 0, 0)
        # minus 1 because of how searchsorted returns indices
        window_i = np.searchsorted(uniq_ymax, y) - 1
        frac_trajs_observable = window_sf[window_i]
        weights = weights*frac_trajs_observable
    full_cdf = np.cumsum(1/weights)  # before repeats removed
    i = 0
    for xi, xx in enumerate(x):
        while i + 1 < num_obs and np.isclose(y[i+1], xx):
            i += 1
        cdf[xi] = full_cdf[i]
    if normalize:
        cdf = cdf/full_cdf[-1]
    return x, cdf


def ecdf_simple(waits, T, pad_left_at_x=0):
    """cdf of interior times (ts > 0) observed in window of size T"""
    ts, counts = np.unique(waits, return_counts=True)
    i = np.argsort(ts)
    ts = ts[i]
    counts = counts[i]
    weights = T/(T - ts)
    return np.insert(ts, 0, pad_left_at_x), \
        np.insert(np.cumsum(counts*weights), 0, 0)/len(weights)
