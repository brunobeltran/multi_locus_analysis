"""
Tools for analysis of distributional statistics of trajectories.

Anything that's not a moment of the trajectory that you want to analyze, you
typically have to look at distributionally (e.g. displacement distributions,
waiting time distributions, etc.). This module contains code that
facilitates working with probability and cumulative distribution functions,
bootstrapping, etc., as well as various miscellaneous distributional tools like
normality tests.
"""
from functools import partial
import warnings

import numpy as np
import pandas as pd
import scipy
import scipy.stats
from scipy.signal import savgol_filter
import statsmodels.stats.proportion as binom

from .moments import pos_to_all_vel
from ..util import mark_code_only


def power_law_slope_mle(x, xmin, N=None):
    if N is None:
        N = len(x)
    return 1 + N / np.sum(np.log(x/xmin))


def ecdf(y, y_allowed=None, auto_pad_left=False, pad_left_at_x=None):
    """Compute empirical cumulative distribution function (eCDF) from data.

    Parameters
    ----------
    y : (N,) array_like
        Values of the data.
    y_allowed : (M,) array_like
        Unique values that the data can take. Mostly useful for adding
        eCDF values at locations where data could or should have been observed
        but none was recorded.
    auto_pad_left : bool
        If left False, the data will not have a data value at the point where
        the eCDF equals zero. Use mean inter-data spacing to automatically
        generate an aesthetically reasonable such point.
    pad_left_at_x : bool
        If ``auto_pad_left`` is False, you may explicitly specify the value at
        which to add the leftmost extra point.

    Returns
    -------
    x : (M,) array_like
        The values at which the eCDF was computed. By default
        :code:`np.sort(np.unique(y))`.
    cdf : (M,) array_like
        Values of the eCDF at each x.

    Notes
    -----
    If using ``y_allowed``, the *pad_left* parameters are redundant, and should
    typically be left False/None.
    """
    y = np.array(y)
    y.sort()
    if y_allowed is not None:
        x = np.unique(y_allowed)
    else:
        x = np.unique(y)
    x.sort()
    if auto_pad_left:
        dx = np.mean(np.diff(x))
        x = np.insert(x, 0, x[0] - dx)
    elif pad_left_at_x is not None:
        if not np.isclose(pad_left_at_x, x[0]):
            if x[0] < pad_left_at_x:
                warnings.warn('pad_left_at_x not left of x in ecdf_windowed! '
                              'Ignoring...')
            else:
                x = np.insert(x, 0, pad_left_at_x)
    num_obs = len(y)
    cdf = np.zeros(x.shape, dtype=np.dtype('float'))
    i = 0
    for xi, xx in enumerate(x):
        while i < num_obs and y[i] <= xx:
            i += 1
        cdf[xi] = i/float(num_obs)
    return x, cdf


def _double_up(x):
    """[1,2,3] to [1,1,2,2,3,3]"""
    return np.tile(x, (2, 1)).T.flatten()


def bars_given_hist(y, bins):
    return _double_up(bins)[1:-1], _double_up(y)


def bars_given_discrete_cdf(x, cdf):
    """like bars_given_cdf for when you've used ecdf's times_allowed arg."""
    x_mid = (x[1:] + x[:-1]) / 2
    real_cdf = cdf[:-1]
    X, Y = bars_given_cdf(x_mid, real_cdf)
    return np.insert(X, 0, [0, 0]), np.insert(Y, 0, [0, 0])


def bars_given_cdf(x, cdf):
    """takes x, cdf from cdf_exact* functions and makes a plottable histogram
    by tracing out the PDF. Works well for CDFs that come from observations on
    a fixed grid, and not well for continuous observations. (i.e.
    discrete_trajectory_to_wait_times output will work well, but not
    state_changes_to_wait_times)."""
    if np.any(np.diff(x) == 0):
        raise ValueError('Values in x repeated!')
    pmf = np.diff(cdf)/np.diff(x)
    # tested by inspection
    return _double_up(x)[1:-1], _double_up(pmf)


def smooth_pdf(x, cdf, bw_method=None):
    """Takes x, cdf from cdf_exact* and returns a kernel density estimator that
    can be evaluated at any X to get an estimate of pdf(X).

    use bw_method to specify the way scipy.stats.gaussian_kde should determine
    the bandwidth of the gaussian. """
    Dcdf = np.diff(cdf)
    Dx = x[1:]
    return scipy.stats.gaussian_kde(Dx, weights=Dcdf,
                                    bw_method=bw_method)


def simultaneous_confint_from_cdf(alpha, n_samples, x, cdf):
    return binom.multinomial_proportions_confint(
        n_samples*np.diff(cdf), alpha=alpha, method='sison-glaz'
    )


def pointwise_confint_from_cdf(alpha, n_samples, x, cdf, bonferroni=True):
    n_bins = len(x) - 1
    if bonferroni:
        alpha = alpha/n_bins
    # pmf = np.diff(cdf)*n_samples
    # WARNING: not clear that this code ever worked, what is it supposed to be
    # doing?
    return binom.proportion_confint(cdf, n_samples, alpha)


def bars_given_confint(x, confint):
    """takes x, confint from cdf_exact*, binom.multinomial_proportions_confint
    (respectively), and rescales confint to correctly fit around the output of
    bars_given_cdf in an aesthetic way.
    """
    # binom speaks in counts, but bars_given_cdf speaks in density, so we need
    # to translate to densities by scaling by dx
    if np.any(np.diff(x) == 0):
        raise ValueError('Values in x repeated!')
    confint = confint/np.tile(np.diff(x), (2, 1)).T
    conf_bar_y = np.stack((
        _double_up(confint[:, 0]),
        _double_up(confint[:, 1])
    ))
    return _double_up(x)[1:-1], conf_bar_y


def sample_from_cdf(n, x, cdf):
    """Takes a sample count, cdf in the form x,cdf, like from output of
    cdf_exact* functions [i.e. pairs of (x, P(X<=x))]. Samples from the
    empirical distribution function at the maximum "x" resolution allowed by x.

    Returns fraction of the resampled data that fell into each bin. In other
    words, it returns pmf, as if one had done:

    >>> pmf, x = np.histogram(samples, bins=x)
    """
    r = np.random.rand(n)
    i = np.searchsorted(cdf, r)
    pmf = np.bincount(i)[1:]
    pmf = pmf/np.sum(pmf)
    extra_zeros = len(x) - len(pmf) - 1
    return np.concatenate([pmf, np.zeros((extra_zeros,))])


def bootstrapped_pmf_confint(n_samples, alpha, x, cdf, num_bootstraps=1000,
                             bonferroni=True):
    r"""Given an empirical cdf (x, cdf), this function generates bootstrapped
    error bars that represent, pointwise, the area that a second observation of
    n_samples (need not equal the number of samples used to generate (x, cdf))
    would lie between with probability 1-alpha if it had a true CDF given by
    the (continuous, linear interpolation of the) empircal CDF.

    Parameters
    ----------
    n_samples : int
        How many samples the secondary measurement has. This is the number of
        data points drawn in each bootstrap iteration.
    alpha : float \in [0,1]
        1 - confidence level desired
    x : (N,) array_like
        Values at which the empirical CDF was measured
    cdf : (N,) array_like
        Values of the empirical CDF
    num_bootstraps : (optional) int
        Number of bootstrapping interations to perform. WARNING: scales the
        memory required for now.
    bonferroni : (optional) bool
        Whether to scale alpha based on the number of bins so that the plot can
        be used to visually assert pointwise statistical significance at the
        requested alpha.

    Returns
    -------
    confint : (2,N-1) array_like
        Upper and lower bounds of the confidence interval calculated.
    """
    n_bins = len(x) - 1
    if bonferroni:
        alpha = alpha/n_bins
    bootstrapped_pmfs = np.zeros((n_bins, num_bootstraps))
    for i in range(num_bootstraps):
        bootstrapped_pmfs[:, i] = sample_from_cdf(n_samples, x, cdf)
    # axis arg is index of array's shape that will be "deleted"
    # we pass two alphas, this new dimension will be prepended
    return np.quantile(bootstrapped_pmfs, [alpha, 1-alpha],
                       overwrite_input=True, axis=1)


def bootstrapped_pmf_confint_bars(n_samples, x, cdf, num_bootstraps=1000):
    """same as non-bars version, but returns the actual samples as pmf bars,
    ready to plot."""
    n_bins = len(x) - 1
    bootstrapped_pmfs = np.zeros((n_bins, num_bootstraps))
    for i in range(num_bootstraps):
        bootstrapped_pmfs[:, i] = sample_from_cdf(n_samples, x, cdf)
    bootstrapped_cdfs = np.cumsum(bootstrapped_pmfs, axis=0)
    bootstrapped_cdfs = np.concatenate(
        [np.zeros((1, num_bootstraps)), bootstrapped_cdfs]
    )
    pmf_bars = np.zeros((2*n_bins, num_bootstraps))
    for i in range(num_bootstraps):
        _, pmf_bars[:, i] = bars_given_cdf(x, bootstrapped_cdfs[:, i])
    x_bars, _ = bars_given_cdf(x, bootstrapped_cdfs[:, 0])
    return x_bars, pmf_bars


def bootstrapped_pmf_from_waits_(
        n_samples, num_bootstraps, times, window_sizes, times_allowed,
        progress_bar=False, **kwargs):
    """Does shared work of creating non-bar pmfs, used by other
    bootstrapped_pmf_from_waits_* functions."""
    from ..finite_window import ecdf_windowed
    # the most common bootstrap thing is to draw samples your own size
    if n_samples is None:
        n_samples = len(times)
    n_bins = len(times_allowed) - 1
    num_waits = len(times)
    bootstrapped_pmfs = np.zeros((n_bins, num_bootstraps))
    for i in range(num_bootstraps):
        if progress_bar:
            print('\r{x:02d}%\r'.format(x=int(float(i)/num_bootstraps*100)),
                  end='')
        samples = np.floor(np.random.rand(n_samples)*num_waits).astype(int)
        t = times[samples]
        w = window_sizes[samples]
        _, cdf = ecdf_windowed(t, w, times_allowed, **kwargs)
        bootstrapped_pmfs[:, i] = np.diff(cdf)
    return bootstrapped_pmfs


def bootstrapped_pmf_from_waits(
        times, window_sizes, times_allowed, n_samples=None, alpha=0.05,
        num_bootstraps=1000, bonferroni=True, progress_bar=False, **kwargs):
    """Takes n_samples, num_bootstraps (# iterations), and calculates the pmf
    of the data num_bootsraps times using n_samples-sized samples drawn with
    replacement from the wait_times/windows that were passed.

    The internal call to cdf_exact_given_windows_quinn needs you to use the
    times_allowed argument, but all kwargs are forwarded to that function just
    in case."""
    n_bins = len(times_allowed) - 1
    if bonferroni:
        alpha = alpha/n_bins
    bootstrapped_pmfs = bootstrapped_pmf_from_waits_(
        n_samples, num_bootstraps, times, window_sizes, times_allowed,
        progress_bar=progress_bar, **kwargs)
    return np.quantile(bootstrapped_pmfs, [alpha, 1-alpha],
                       overwrite_input=True, axis=1)


def bootstrapped_pmf_from_waits_bars(
        times, window_sizes, times_allowed, n_samples=None,
        num_bootstraps=1000, **kwargs):
    """Same as bootstrapped_pmf_from_waits, but returns bars_given_cdf, ready to
    plot results."""
    n_bins = len(times_allowed) - 1
    bootstrapped_pmfs = bootstrapped_pmf_from_waits_(
        n_samples, num_bootstraps, times, window_sizes, times_allowed, **kwargs
    )
    bootstrapped_cdfs = np.cumsum(bootstrapped_pmfs, axis=0)
    bootstrapped_cdfs = np.concatenate([
        np.zeros((1, num_bootstraps)), bootstrapped_cdfs
    ])
    pmf_bars = np.zeros((2*n_bins, num_bootstraps))
    for i in range(num_bootstraps):
        _, pmf_bars[:, i] = bars_given_cdf(
            times_allowed, bootstrapped_cdfs[:, i]
        )
    x_bars, _ = bars_given_cdf(times_allowed, bootstrapped_cdfs[:, 0])
    return x_bars, pmf_bars


@mark_code_only
def scale_and_test_normality(vx):
    """Should be applied to a vector of velocities, vx"""
    vx = vx[np.isfinite(vx)]
    vx -= np.mean(vx)
    std = np.std(vx)
    if std > 0:
        vx /= np.std(vx)
    num_disps = vx.size
    ksstat, ks_pval = scipy.stats.kstest(vx, 'norm', mode='asymp')
    ngp = np.mean(np.power(vx, 4)) \
        / (3*np.power(np.mean(np.power(vx, 2)), 2)) - 1
    try:
        shapiro_stat, shapiro_pval = scipy.stats.shapiro(vx)
    except Exception:  # shut up linter. don't remember what this raises
        shapiro_stat = shapiro_pval = np.nan
    is_shapiro_pval_accurate = num_disps < 5000  # from scipy docs, v 19.1
    anderson_stat, crit_vals, crit_levels = scipy.stats.anderson(vx)
    if np.isfinite(anderson_stat):
        crit_vals = np.concatenate(([0], crit_vals, [np.inf]))
        crit_alphas = np.concatenate(([100], crit_levels, [0]))/100
        ihigh = np.where(crit_vals < anderson_stat)[0][0]
        ilow = ihigh + 1
    else:
        ihigh = ilow = 0
        crit_vals = [np.nan]
        crit_alphas = [np.nan]
    return {'NGP': ngp, 'KS Statistic': ksstat, 'KS p-value': ks_pval,
            'Shapiro-Wilk Statistic': shapiro_stat,
            'Shapiro-Wilk p-value': shapiro_pval,
            'Shapiro-Wilk p-value is accurate': is_shapiro_pval_accurate,
            'Anderson-Darling Statistic': anderson_stat,
            'Anderson-Darling alpha lower bound': crit_alphas[ilow],
            'Anderson-Darling alpha upper bound': crit_alphas[ihigh],
            'Anderson-Darling upper crit value': crit_vals[ihigh],
            'Anderson-Darling lower crit value': crit_vals[ilow]
            }


@mark_code_only
def add_savgol(traj, window_size, order=3):
    s = 'savgol' + str(window_size) + '_'
    if len(traj['x']) <= window_size:
        traj[s + 'x'] = np.nan
        traj[s + 'y'] = np.nan
        traj[s + 'x_fluct'] = np.nan
        traj[s + 'y_fluct'] = np.nan
    else:
        traj[s + 'x'] = savgol_filter(traj['x'], window_size, order)
        traj[s + 'y'] = savgol_filter(traj['y'], window_size, order)
        traj[s + 'x_fluct'] = traj['x'] - traj[s + 'x']
        traj[s + 'y_fluct'] = traj['y'] - traj[s + 'y']
    return traj


@mark_code_only
def add_savgol_window(data, window_sizes):
    for window_size in window_sizes:
        data = data.groupby(
            ['experiment', 'movie name', 'molecule id']
        ).apply(partial(add_savgol, window_size=window_size))
    return data


@mark_code_only
def get_savgol_vels(data, window_sizes):
    # we're gonna do the same thing twice, basically get vels for each window
    # size
    savgol_vels = []
    for window_size in window_sizes:
        fluct_vels = data.groupby(
            ['experiment', 'movie name', 'molecule id']
        ).apply(partial(
            pos_to_all_vel, delta=1,
            xcol='savgol'+str(window_size)+'_x',
            ycol='savgol'+str(window_size)+'_y')
        )
        fluct_vels['window_size'] = window_size
        savgol_vels.append(fluct_vels)
    savgol_vels = [vel.reset_index().set_index(
        ['experiment', 'movie name', 'molecule id', 'window_size', 'ti',
         'delta']
    ) for vel in savgol_vels]
    savgol_vels = pd.concat(savgol_vels)
# now for _fluct columns
    corrected_vels = []
    for window_size in window_sizes:
        fluct_vels = data.groupby(
            ['experiment', 'movie name', 'molecule id']
        ).apply(partial(
            pos_to_all_vel, delta=1,
            xcol='savgol'+str(window_size)+'_x_fluct',
            ycol='savgol'+str(window_size)+'_y_fluct'
        ))
        fluct_vels['window_size'] = window_size
        corrected_vels.append(fluct_vels)
    corrected_vels = [vel.reset_index().set_index(
        ['experiment', 'movie name', 'molecule id', 'window_size', 'ti',
         'delta']
    ) for vel in corrected_vels]
    return savgol_vels, pd.concat(corrected_vels)


@mark_code_only
def get_savgol_msds(savgol_data, window_sizes):
    msds = []
    for window_size in window_sizes:
        savgol_disps = savgol_data.groupby(
            ['experiment', 'movie name', 'molecule id']
        ).apply(partial(
            pos_to_all_vel, absolute_time=True,
            xcol='savgol'+str(window_size)+'_x_fluct',
            ycol='savgol'+str(window_size)+'_y_fluct'
        ))
        sq_disps = pd.DataFrame(index=savgol_disps.index)
        sq_disps['sq_disp'] = np.power(savgol_disps['vx'], 2) \
            + np.power(savgol_disps['vy'], 2)
        sq_disps['delta_abs'] = savgol_disps['delta_abs'].round(3)
        savgol_msds = sq_disps \
            .groupby(['delta_abs'])['sq_disp'] \
            .agg(['mean', 'std', 'count'])
        savgol_msds['window_size'] = window_size
        savgol_msds = savgol_msds.reset_index().set_index(['window_size',
                                                           'delta_abs'])
        msds.append(savgol_msds)
    return pd.concat(msds)
