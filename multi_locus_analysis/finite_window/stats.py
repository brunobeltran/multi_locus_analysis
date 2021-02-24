import numpy as np

from ..stats import ecdf


def ecdf_windowed(
        times, window_sizes, times_allowed=None, auto_pad_left=None,
        pad_left_at_x=None, window_cumulant=None, normalize=True):
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
        If left False, the data will not have a data value at the point where
        the eCDF equals zero. Use mean inter-data spacing to automatically
        generate an aesthetically reasonable such point.
    pad_left_at_x : bool
        Same as ``auto_pad_left``, but specify the point at which to add
        the leftmost point.
    window_cumulant : (M,) array_like of float
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
    y = times
    ymax = window_sizes
    y = np.array(y)
    ymax = np.array(ymax)
    # variables only needed for multiple window sizes
    ignore_window_cumulant = False
    uniq_ymax = None
    # allow providing single window size
    if ymax.size == 1:
        ignore_window_cumulant = True
        ymax = ymax*np.ones_like(y)
    i = np.argsort(y)
    y = y[i]
    ymax = ymax[i]
    if times_allowed is not None:
        x = np.unique(times_allowed)
    else:
        x = np.unique(y)
    x.sort()
    num_obs = len(y)
    cdf = np.zeros(x.shape, dtype=np.dtype('float'))
    if not ignore_window_cumulant and window_cumulant is None:
        # get fraction of windows that are at *least* of each width
        uniq_ymax, window_cumulant = ecdf(ymax, pad_left_at_x=0)
    weights = (ymax - y)
    if not ignore_window_cumulant:
        window_frac_at_least = 1 - window_cumulant
        # for each observed time, we can get number of windows in which it can
        # have been observed
        if uniq_ymax is None:
            # don't forget to pad_left_at_x=0
            uniq_ymax = np.insert(np.unique(ymax), 0, 0)
        # minus 1 because of how searchsorted returns indices
        window_i = np.searchsorted(uniq_ymax, y) - 1
        frac_trajs_observable = window_frac_at_least[window_i]
        weights = weights*frac_trajs_observable
    full_cdf = np.cumsum(1/weights)  # before repeats removed
    i = 0
    for xi, xx in enumerate(x):
        while i + 1 < num_obs and np.isclose(y[i+1], xx):
            i += 1
        cdf[xi] = full_cdf[i]
    if normalize:
        cdf = cdf/full_cdf[-1]
    if auto_pad_left:
        dx = np.mean(np.diff(x))
        x = np.insert(x, 0, x[0] - dx)
        cdf = np.insert(cdf, 0, 0)
    elif pad_left_at_x is not None:
        if x[0] == pad_left_at_x:
            raise ValueError('pad_left_at_x value already exists in x!')
        x = np.insert(x, 0, pad_left_at_x)
        cdf = np.insert(cdf, 0, 0)
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
