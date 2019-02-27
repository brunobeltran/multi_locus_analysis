""""Interior Windowing" Statistical Corrections

This module is designed to facilitate working with Markov processes which have
been observed over finite windows of time.

Suppose you wish to observe a process that switches between states A and B, and
observe the process over a window of time [0, T].

There will be two types of "observed" times, which we we call "interior times"
(those when the entire time spent in the state is observed) and "exterior
times" (the lengths spent in the state that we observed at the start/end point
of our window).

For example, suppose you are observing two fluorescent loci that are either in
contact (A) or not (B). If you measure for 5s, and the loci begin in contact,
come apart at time t=2s, and rejoin at time t=4s,
    A A A A A A B B B B B B A A A A
    | - - | - - | - - | - - | - - | -- | -- | ...
    |                             |
t = 0s    1s    2s    3s    4s    5s

we would say that T = 5s, and you measured one "interior" time of length 2s,
and two exterior times, one of length 2s and one of length 1s (at the start and
end of the window, respectively).

When histogramming these times, we must apply a statistical correction to the
final histogram weights due to the effects of the finite observation window.

The distribution of the exterior times $f_E(t)$ is exactly equal to the survival
function of the actual distribution $S(t) = 1 - \int_0^t f(s) ds$ normalized to
equal one over the observation interval.
No functions are currently included
that leverage this, since extracting information from the survival function is
likely only worth it if a large fraction of your observations are exterior
times.

On the other hand, the interior times distribution is given by
$f_I(t) = f(t)(T - t)$ for $t\in[0,T]$. In order to plot the actual shape of
$f(t)$ with this biasing removed, we provide the cdf_exact_given_windows
function below.

A typical workflow is, given an array of interior times, :code:`t`, and an array
of the window sizes each time was observed within, :code:`w`, is to extract the
CDF exactly, then optionally convert that to a PDF to display a regular
histogram::

    >>> x, cdf = cdf_exact_given_windows(t, w, pad_left_at_x=0)
    >>> xp, pdf = bars_given_cdf(x, cdf)
    >>> confint = simultaneous_confint_from_cdf(0.05, len(t), x, cdf)
    >>> xc, confs = bars_given_confint(x, confint)
    >>> plt.plot(xp, pdf, 'k', xc, confs, 'r-.')

The remainder of the functions herein are related to extracting confidence
intervals for these distributional estimates.

For details of the derivation, see the Biophysical Journal paper in preparation
by Beltran and Spakowitz.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import scipy
from scipy import stats
from scipy import interpolate
import statsmodels
import statsmodels.stats.proportion as binom
from lifelines import KaplanMeierFitter

import bruno_util
import bruno_util.random
from bruno_util import pandas_util


###############{{{
# simulation code

@bruno_util.random.strong_default_seed
def shortcut_finite_window_sim_2(rands, means, window_size, num_replicates=1,
                                 seed=None, random_state=None):
    """Simulate a two-state system switching between states 0 and 1 with
    waiting times that can be sampled by rands[0] and rands[1], and whose
    average waiting times are are means[0], means[1], respectively.

    Parameters
    ----------
    rands : (2,) List[scipy.stats.rv_continuous]
        Callable that takes "random_state" and "size" kwargs that accept
        a np.random.RandomState seed and a tuple specifying array sizes, resp.
    means : (2,) array_like
        average waiting times for each of the states
    window_size : float
        the width of the window over which the observation takes place
    num_replicates : int
        number of times to run the simulation
    seed : np.random.RandomState
        state to start the simulation with

    Returns
    -------
    df : pd.DataFrame
        The start/end times of each waiting time simulated. This data frame has
        columns=['replicate', 'state', 'start_time', 'end_time',
        'window_start', 'window_end'].

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
    some offset before the window, as in explicit_finite_window_sim_2.
    """
    # np.concatenate can't handle concatenating nothing into nothing, so we
    if num_replicates <= 0:
        return pd.DataFrame(columns=['replicate', 'state', 'start_time',
                                     'end_time', 'window_start', 'window_end'])
    pool_size_guess = int(num_replicates*1.5*window_size/(means[0] + means[1]))
    pool_size_guess = max(pool_size_guess, 16)
    rands = [bruno_util.random.make_pool(rand, pool_size_guess, random_state=random_state)
             for rand in rands]

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
        sim_state = 1 - sim_state # switch between 0 and 1
        while t < window_size:
            end_times[-1].append(t)
            start_times[-1].append(t)
            states[-1].append(sim_state)
            t += rands[sim_state]()
            sim_state = 1 - sim_state # switch between 0 and 1
        end_times[-1].append(t)
    start_times = np.concatenate([np.array(ts) for ts in start_times])
    end_times = np.concatenate([np.array(ts) for ts in end_times])
    replicate_ids = np.concatenate([i*np.ones_like(ts) for i, ts in enumerate(states)])
    states = np.concatenate([np.array(ss) for ss in states])
    df = pd.DataFrame.from_dict({'replicate': replicate_ids, 'state': states,
                                 'start_time': start_times, 'end_time': end_times})
    df['window_start'] = 0
    df['window_end'] = window_size
    return df

@bruno_util.random.strong_default_seed
def explicit_finite_window_sim_2(rands, window_size, offset, num_replicates=1,
                        seed=None, random_state=None):
    """Simulate a two-state system switching between states 0 and 1 with
    waiting times that can be sampled by rands[0] and rands[1], and whose
    average waiting times are are means[0], means[1], respectively. Emulate
    time-homogeneity of the finite observation window by starting the
    simulation at some specificed offset before time zero, then only recording
    between times 0 and window_size.

    Parameters
    ----------
    rands : (2,) List[scipy.stats.rv_continuous]
        Callable that takes "random_state" and "size" kwargs that accept
        a np.random.RandomState seed and a tuple specifying array sizes, resp.
    window_size : float
        The width of the window over which the observation takes place
    offset : float
        The (negative) time at which to start (in state 0) in order to equilibrate teh
        simulation state by time t=0.
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
    rands = [bruno_util.random.make_pool(rand, pool_size_guess, random_state=state)
             for rand in rands]

    start_times = []
    end_times = []
    states = []
    wait_type = []
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
            sim_state += 1
            r = rands[sim_state]()
        states[i].append(sim_state)
        start_times[i].append(t)
        while t + r < window_size:
            t += r
            end_times[i].append(t)
            sim_state += 1
            r = rands[sim_state]()
            states[i].append(sim_state)
            start_times[i].append(t)
        end_times[i].append(t + r)
    start_times = np.concatenate([np.array(ts) for ts in start_times])
    end_times = np.concatenate([np.array(ts) for ts in end_times])
    replicate_ids = np.concatenate([i*np.ones_like(ts) for i, ts in enumerate(states)])
    states = np.concatenate([np.array(ss) for ss in states])
    df = pd.DataFrame.from_dict({'replicate': replicate_ids, 'state': states,
                                 'start_time': start_times, 'end_time': end_times})
    df['window_start'] = 0
    df['window_end'] = window_size
    return df

# end simulation code
###############}}}

###############{{{
# keep in dataframe functions

def state_changes_to_wait_times(data, start_times_col, state_col):
    """Converts a Series of state change times into a dataframe containing each
    wait time, with its start, end, rank order, and the state it's leaving.

    This function deals with "continuous" wait times, as in not measured at
    discrete time points (on a grid), so the wait times it returns are
    exact."""
    pass

def discrete_trajectory_to_wait_times(data, time_column, states_column):
    """Converts a discrete trajectory to a dataframe containing each wait time,
    with its start, end, rank order, and the state it's leaving.

    Discrete here means that the state of the system was observed at finite
    time points (on a lattice in time), as opposed to a system where the exact
    times of transitions between states are known.

    Because a discrete trajectory only bounds the wait times, and does not
    determine their exact lengths (as a continuous trajectory might),
    additional columns are included that explictly bound the wait times, in
    addition to returning the "natural" estimate.

    Parameters
    ----------
    data : pd.DataFrame
        should have at least states_column and time_column columns, and already
        be groupby'd so that there's only one "trajectory" within the
        DataFrame. One row should correspond to an observation at a particular
        time point.
    time_column : string
        the name of the column containing the time of each time point
    states_column : string
        the name of the column containing the state for each time point

    Returns
    -------
    wait_df : pd.DataFrame
        columns are ['wait_time', 'start_time', 'end_time', 'wait_state',
        'wait_type', 'wait_bound_low', 'wait_bound_high'], where
        [wait,end,start]_time columns are self explanatory, wait_state is the
        value of the states_column during that waiting time, and wait_type is
        one of 'interior', 'left exterior', 'right exterior', 'full exterior',
        depending on what kind of waiting time was observed. See
        discrete_trajectory_to_wait_times_ for detailed explanation of these
        categories. The 'wait_bound_low/high' columns contain the
        minimum/maximum possible value of the wait time (resp.), given the
        observations.

        The default index is named "rank_order", since it tracks the order,
        zero-indexed in whihc the wait times occured.

    .. _discrete_trajectory_to_wait_times:

    Notes
    -----
    the following types of wait times are of interest to us

    1) *interior* censored times: whenever you are observing a switching
    process for a finite amount of time, any waiting time you observe the
    entirety of is called "interior" censored

    2) *left exterior* censored times: whenever the waiting time you observe
    started before you began observation (it overlaps the "left" side of your
    interval of observation)

    3) *right exterior* censored times: same as above, but overlappign the
    "right" side of your interval of observation.

    4) *full exterior* censored times: whenever you observe the existence of a
    single, particular state throughout your entire window of observation.
    """

    wait_columns = ['wait_time', 'start_time', 'end_time', 'wait_state', 'wait_type']
    states = data[states_column].values
    times = data[time_column].values
    num_measurements = len(data)

    # now iterate through valid part of trajectory to establish wait times
    start_times = []
    end_times = []
    wait_state = []
    wait_type = []
    k0 = 0 # index at which current state began
    state = states[i]
    state_has_changed = False
    for k in range(num_measurements):
        # if no state change, continue
        if states[k] == state:
            continue
        # otherwise, store change
        start_times.append(times[k0])
        end_times.append(times[k])
        wait_state.append(state)
        # if this is the first state change, we store it separately
        if not state_has_changed:
            wait_type.append('left exterior')
            state_has_changed = True
        # otherwise, a normal state change
        else:
            wait_type.append('interior')
        # either way, state has changed
        state = states[k]
        k0 = k
    # also store the time spent in final state
    start_times.append(times[k0])
    end_times.append(times[k])
    # state type also stored specially
    wait_state.append(state)
    if not state_has_changed:
        wait_type.append('full exterior')
    else:
        wait_type.append('right exterior')
    start_times = np.array(start_times)
    end_times = np.array(end_times)
    wait_times = end_times - start_times
    df = pd.DataFrame({'start_time': start_times, 'end_time': end_times,
                       'wait_time': wait_times, 'wait_state': wait_state})
    df.index.name = 'rank_order'
    df['window_size'] = valid_time_window
    return df

# end keep in dataframe functions
###############}}}

###############{{{
# useful aliases for random variables

def abs_norm_rvs_factory(*args, **kwargs):
    n = scipy.stats.norm(*args, **kwargs)
    def rvs(*args, **kwargs):
        return np.abs(n.rvs(*args, **kwargs))
    return rvs

# end useful aliases for random variables
###############}}}


###############{{{
# statistical window corrections
def cdf_exact(y, y_allowed=None, auto_pad_left=False, pad_left_at_x=None):
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
        Same as :param:`auto_pad_left`, but specify the point at which to add
        the leftmost point.

    Returns
    -------
    x : (M,) array_like
        The values at which the eCDF was computed. By default
        :code:`np.sort(np.unique(y))`.
    cdf : (M,) array_like
        Values of the eCDF at each x.

    Notes
    -----
    If using :param:`times_allowed`, the *pad_left* parameters are redundant.
    """
    y = np.array(y)
    y.sort()
    if y_allowed is not None:
        x = np.unique(y_allowed)
    else:
        x = np.unique(y)
    x.sort()
    num_obs = len(y)
    cdf = np.zeros(x.shape, dtype=np.dtype('float'))
    i = 0
    for xi, xx in enumerate(x):
        while i < num_obs and y[i] <= xx:
            i += 1
        cdf[xi] = i/float(num_obs)
    if auto_pad_left == True:
        dx = np.mean(np.diff(x))
        x = np.insert(x, 0, x[0] - dx)
        # x = np.append(x, x[-1] + dx)
    elif pad_left_at_x is not None:
        x = np.insert(x, 0, pad_left_at_x)
        cdf = np.insert(cdf, 0, 0)
    return x, cdf

exact_cdf = cdf_exact

def cdf_exact_given_windows(times, window_sizes, times_allowed=None,
        auto_pad_left=None, pad_left_at_x=None, window_func=None):
    """Compute empirical cumulative distribution function (eCDF) from data
    taken within a finite observation interval.

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
        but none was recorded (i.e. if a movie was taken with a given framerate
        but not all possible window lengths were observed.
    auto_pad_left : bool
        If left False, the data will not have a data value at the point where
        the eCDF equals zero. Use mean inter-data spacing to automatically
        generate an aesthetically reasonable such point.
    pad_left_at_x : bool
        Same as :param:`auto_pad_left`, but specify the point at which to add
        the leftmost point.
    window_func : function
        The function used to reweight the observations. The default is
        equivalent to
        :code:`lambda yi: np.sum(np.heaviside(ymax - yi, 0)*(ymax - yi)))`.


    Returns
    -------
    x : (M,) array_like
        The values at which the eCDF was computed. By default
        :code:`np.sort(np.unique(y))`.
    cdf : (M,) array_like
        Values of the eCDF at each x.

    Notes
    -----
    If using :param:`times_allowed`, the *pad_left* parameters are redundant.
    """
    y = times
    ymax = window_sizes
    y = np.array(y)
    ymax = np.array(ymax)
    # allow providing single window size
    if ymax.size == 1:
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
    if window_func is None:
        window_func = np.vectorize(lambda yi: np.sum(np.heaviside(ymax - yi, 0)*(ymax - yi)))
    weights = window_func(y)
    full_cdf = np.cumsum(1/weights) # before repeats removed
    i = 0
    for xi, xx in enumerate(x):
        while i + 1 < num_obs and y[i+1] == xx:
            i += 1
        cdf[xi] = full_cdf[i]/full_cdf[-1]
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

def double_up(x):
    """[1,2,3] to [1,1,2,2,3,3]"""
    return np.tile(x, (2, 1)).T.flatten()

def bars_given_cdf(x, cdf):
    """takes x, cdf from cdf_exact* functions and makes a plottable histogram
    by tracing out the PDF. """
    if np.any(np.diff(x) == 0):
        raise ValueError('Values in x repeated!')
    pmf = np.diff(cdf)/np.diff(x)
    # tested by inspection
    return double_up(x)[1:-1], double_up(pmf)

def simultaneous_confint_from_cdf(alpha, n_samples, x, cdf):
    return binom.multinomial_proportions_confint(n_samples*np.diff(cdf),
            alpha=alpha, method='sison-glaz')

def pointwise_confint_from_cdf(alpha, n_samples, x, cdf, bonferroni=True):
    n_bins = len(x) - 1
    if bonferroni:
        alpha = alpha/n_bins
    pmf = np.diff(cdf)*n_samples
    return binom.proportion_confint(countsdf['pair'], n_samples, alpha)

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
    conf_bar_y = np.stack((double_up(confint[:,0]), double_up(confint[:,1])))
    return double_up(x)[1:-1], conf_bar_y

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
    """Given an empirical cdf (x, cdf), this function generates bootstrapped
    error bars that represent, pointwise, the area that a second observation of
    n_samples (need not equal the number of samples used to generate (x, cdf))
    would lie between with probability 1-alpha if it had a true CDF given by the
    (continuous, linear interpolation of the) empircal CDF.

    Inputs
    ------
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

    Output
    ------
    confint : (2,N-1) array_like
        Upper and lower bounds of the confidence interval calculated.
    """
    n_bins = len(x) - 1
    if bonferroni:
        alpha = alpha/n_bins
    bootstrapped_pmfs = np.zeros((n_bins, num_bootstraps))
    for i in range(num_bootstraps):
        bootstrapped_pmfs[:,i] = sample_from_cdf(n_samples, x, cdf)
    # axis arg is index of array's shape that will be "deleted"
    # we pass two alphas, this new dimension will be prepended
    return np.quantile(bootstrapped_pmfs, [alpha, 1-alpha], overwrite_input=True, axis=1)

def bootstrapped_pmf_confint_bars(n_samples, x, cdf, num_bootstraps=1000):
    """same as non-bars version, but returns the actual samples as pmf bars,
    ready to plot."""
    n_bins = len(x) - 1
    bootstrapped_pmfs = np.zeros((n_bins, num_bootstraps))
    for i in range(num_bootstraps):
        bootstrapped_pmfs[:,i] = sample_from_cdf(n_samples, x, cdf)
    bootstrapped_cdfs = np.cumsum(bootstrapped_pmfs, axis=0)
    bootstrapped_cdfs = np.concatenate([np.zeros((1, num_bootstraps)), bootstrapped_cdfs])
    pmf_bars = np.zeros((2*n_bins, num_bootstraps))
    for i in range(num_bootstraps):
        _, pmf_bars[:,i] = bars_given_cdf(x, bootstrapped_cdfs[:,i])
    x_bars, _ = bars_given_cdf(x, bootstrapped_cdfs[:,0])
    return x_bars, pmf_bars

def bootstrapped_pmf_from_waits_(n_samples, num_bootstraps, times, window_sizes,
                              times_allowed, progress_bar=False, **kwargs):
    """Does shared work of creating non-bar pmfs, used by other
    bootstrapped_pmf_from_waits_* functions."""
    # the most common bootstrap thing is to draw samples your own size
    if n_samples is None:
        n_samples = len(times)
    n_bins = len(times_allowed) - 1
    num_waits = len(times)
    bootstrapped_pmfs = np.zeros((n_bins, num_bootstraps))
    for i in range(num_bootstraps):
        if progress_bar:
            print('\r{x:02d}%\r'.format(x=int(float(i)/num_bootstraps*100)), end='')
        samples = np.floor(np.random.rand(n_samples)*num_waits).astype(int)
        t = times[samples]
        w = window_sizes[samples]
        _, cdf = cdf_exact_given_windows_quinn(t, w, times_allowed, **kwargs)
        bootstrapped_pmfs[:,i] = np.diff(cdf)
    return bootstrapped_pmfs

def bootstrapped_pmf_from_waits(times, window_sizes,
        times_allowed, n_samples=None, alpha=0.05, num_bootstraps=1000,
        bonferroni=True, progress_bar=False, **kwargs):
    """Takes n_samples, num_bootstraps (# iterations), and calculates the pmf
    of the data num_bootsraps times using n_samples-sized samples drawn with
    replacement from the wait_times/windows that were passed.

    The internal call to cdf_exact_given_windows_quinn needs you to use the
    times_allowed argument, but all kwargs are forwarded to that function just
    in case."""
    n_bins = len(times_allowed) - 1
    if bonferroni:
        alpha = alpha/n_bins
    bootstrapped_pmfs = bootstrapped_pmf_from_waits_(n_samples, num_bootstraps,
            times, window_sizes, times_allowed, progress_bar=progress_bar, **kwargs)
    return np.quantile(bootstrapped_pmfs, [alpha, 1-alpha], overwrite_input=True, axis=1)

def bootstrapped_pmf_from_waits_bars(times, window_sizes, times_allowed,
        n_samples=None, num_bootstraps=1000, **kwargs):
    """Same as bootstrapped_pmf_from_waits, but returns bars_given_cdf, ready to
    plot results."""
    n_bins = len(times_allowed) - 1
    bootstrapped_pmfs = bootstrapped_pmf_from_waits_(n_samples, num_bootstraps,
            times, window_sizes, times_allowed, **kwargs)
    bootstrapped_cdfs = np.cumsum(bootstrapped_pmfs, axis=0)
    bootstrapped_cdfs = np.concatenate([np.zeros((1, num_bootstraps)), bootstrapped_cdfs])
    pmf_bars = np.zeros((2*n_bins, num_bootstraps))
    for i in range(num_bootstraps):
        _, pmf_bars[:,i] = bars_given_cdf(x, bootstrapped_cdfs[:,i])
    x_bars, _ = bars_given_cdf(x, bootstrapped_cdfs[:,0])
    return x_bars, pmf_bars

# end statistical window corrections
###############}}}