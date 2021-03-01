import warnings

import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype

censor_t = CategoricalDtype(
    categories=['interior', 'left exterior', 'right exterior',
                'full exterior'],
    ordered=False
)


def state_changes_to_wait_times(traj):
    """
    DEPRECATED: in favor of :func:`simulations_to_observations`.

    Converts the output of :func:`ab_window_fast` into a
    :class:`pd.DataFrame` containing each wait time, with its start, end, rank
    order, and the state it's leaving.

    This function deals with "continuous" wait times, as in not measured at
    discrete time points (on a grid), so the wait times it returns are
    exact."""
    warnings.warn('Use sim_to_obs instead!', DeprecationWarning)
    waits = traj.copy()
    num_waits = len(waits)
    waits['rank_order'] = np.arange(num_waits) + 1
    waits.set_index('rank_order', inplace=True)
    if num_waits == 0:
        return waits
    waits['wait_time'] = waits['end_time'] - waits['start_time']
    waits.loc[1, 'wait_time'] \
        = waits.loc[1, 'end_time'] - waits.loc[1, 'window_start']
    waits.loc[num_waits, 'wait_time'] = \
        waits.loc[num_waits, 'window_end'] \
        - waits.loc[num_waits, 'start_time']
    waits['window_size'] = waits['window_end'] - waits['window_start']
    waits['wait_type'] = 'interior'
    if num_waits == 1:
        waits['wait_type'] = 'full exterior'
        waits['wait_time'] = waits['window_size'].copy()
    else:
        waits.loc[1, 'wait_type'] = 'left exterior'
        waits.loc[num_waits, 'wait_type'] = 'right exterior'
    return waits


def traj_to_waits(*args, **kwargs):
    """Alias of :func:`state_changes_to_wait_times`"""
    return state_changes_to_wait_times(*args, **kwargs)


def simulations_to_observations(simulations, traj_cols=['replicate']):
    """
    Vectorized extracting of "observed" wait times.
    """
    sim = simulations
    obs = sim[traj_cols + ['state']].copy()
    # observed times are only within the interval
    obs['start_time'] = np.max(sim[['start_time', 'window_start']].to_numpy(),
                               axis=1)
    obs['end_time'] = np.min(sim[['end_time', 'window_end']].to_numpy(),
                             axis=1)
    # allows vectorized computation of all wait times
    obs['wait_time'] = obs['end_time'] - obs['start_time']
    obs['window_size'] = sim['window_end'] - sim['window_start']

    # vectorized-ly extract rank order/num_waits to classify each wait time
    obs['rank_order'] = obs.groupby(traj_cols)['start_time'].rank().astype(int)
    obs.set_index(traj_cols, inplace=True)
    obs['num_waits'] = obs.groupby(traj_cols)['rank_order'].max()

    # default to assuming its an interior time
    obs['wait_type'] = 'interior'
    # don't think there's a way to make this categorical at "assign" time
    obs['wait_type'] = obs['wait_type'].astype(censor_t)
    # the first and last ones are left/right exterior
    obs.loc[obs['rank_order'] == 1, 'wait_type'] = 'left exterior'
    obs.loc[
        obs['rank_order'] == obs['num_waits'],
        'wait_type'
    ] = 'right exterior'
    # and isolated wait times are doubly exterior
    obs.loc[obs['num_waits'] == 1, 'wait_type'] = 'full exterior'

    return obs.set_index('rank_order', append=True)


sim_to_obs = simulations_to_observations


def state_changes_to_movie_frames(
        traj, times, state_col='state', start_times_col='start_time',
        end_times_col='end_time', wait_type_col='wait_type'):
    """
    Convert state changes into discrete-time observations of state.

    Takes a Series of state change times into a Series containing
    observations at the times requested. The times become the index.

    Parameters
    ----------
    times : (N,) array_like
        times at which to "measure" what state we're in to make the new
        trajectories.
    traj : pd.DataFrame
        should have *state_col* and *start_times_col* columns. the values of
        *state_col* will be copied over verbatim.
    state_col : str, default: 'state'
        name of column containing the state being transitioned out of for each
        measurement in *traj*.
    start_times_col : str, default: 'start_times'
        name of column containing times at which *traj* changed state
    end_times_col : (optional) str
        by default, the function assumes that times after the last provided
        state transition time are in the same state. if passed, this column is
        used to determine at what time the last state "finished". Times after
        this will be labeled as NaN. Analagously to how start times are
        treated, if the end time *exactly* matches the "transition" time,
        assume this is an "exterior" measurement.

    Returns
    -------
    movie : pd.Series
        Series defining the "movie" with frames taken at `times` that
        simply measures what state `traj` is in at each frame. index is
        `times`, `state_col` is used to name the Series.

    Notes
    -----
    A start time means that if we observe at that time, the state transition
    will have already happened (right-continuity), except in the case when the
    transition happens at *exactly* the window end time. This is confusing in
    words, but simple to see in an example (see the example below).

    Examples
    --------
    For the DataFrame

        >>> df = pd.DataFrame([['A',  -1, 0.1], ['B', 0.1, 1.0]],
        >>>     columns=['state', 'start_time', 'end_time'])

    the discretization into tenths of seconds would give

        >>> state_changes_to_movie_frames(df, times=np.linspace(0, 1, 11),
        >>>     end_times_col='end_time')
        t
        0.0      A
        0.1      B
        0.2      B
        0.3      B
        0.4      B
        0.5      B
        0.6      B
        0.7      B
        0.8      B
        0.9      B
        1.0      B
        Name: state, dtype: object

    Notice in particular how at 0.1, the state is already 'B'. However at time
    1.0 the state is not already "unknown". This is what is meant by the Notes
    section above.

    If the `end_times_col` argument is omitted, then the last observed state is
    assumed to continue for all `times` requested from then on:

        >>> state_changes_to_movie_frames(df, times=np.linspace(0, 2, 11))
        t
        0.0    A
        0.2    B
        0.4    B
        0.6    B
        0.8    B
        1.0    B
        1.2    B
        1.4    B
        1.6    B
        1.8    B
        2.0    B
        Name: state, dtype: object


    """
    if len(traj) <= 0:
        raise ValueError('Need a non-empty trajectory, otherwise how will I '
                         'know what the possible states are?')
    times = np.sort(np.array(times))
    traj.sort_values(start_times_col)
    if times[0] < traj[start_times_col].iloc[0]:
        raise ValueError('Requested a time before any measurements were made!')
    # initialize in a random state to get dtype correct
    movie = pd.Series(traj[state_col].iloc[0], index=times)
    i = 0
    i_max = len(traj)
    for time in times:
        while i < i_max and traj[start_times_col].iloc[i] <= time:
            i += 1
        movie.loc[time] = traj[state_col].iloc[i-1]
    if end_times_col is not None:
        last_observed_time = traj[end_times_col].iloc[-1]
        if last_observed_time < times[-1]:
            first_bad_time = np.argmax(last_observed_time < times)
            movie.loc[times[first_bad_time]:] = np.nan
    movie.name = state_col
    movie.index.name = 't'
    return movie


def simulation_to_frame_times(simulation, t, traj_cols=['replicate']):
    """
    Vectorized conversion into frames, using np.searchsorted.

    Getting what "frames" each transition happened on is pretty easy. Consider
    the example trajectory with ``window_start=1.2``, ``window_end=2.9``, and
    transition times:

        >>> transitions = [0.9, 1.5, 1.7, 2.5, 3.1]

    with possible frames:

        >>> frames = [0, 1, 2, 3, 4, 5, 6]

    It may seem that by definition, we can only observe the transitions between
    *window_start* and *window_end*, however, this is not actually true here,
    because the stationary, "continuous" Markov renewal process that we've
    generated using the `fw.simulation` functions is valid from the "true"
    start of the left exterior time to the "true" end of the right exterior
    time, even if we normally discard these times for the purposes of
    simulating censored statistics. So the discrete movie's window size will
    just be the number of frames that are in the interval defined by
    ``[min(transitions), max(transitions)]``.

    Similarly, if instead ``window_end=3.05``, then nothing actually changes,
    because either way in a discrete movie we directly observe the state at
    ``t=3``.

    We've got nice strats for extracting the first and last possible frame in a
    vectorized way. So the algorithm is to first do that. Then, we discard
    transition times that are beyond our observation frames in either
    direction, and searchsort the remaining ones into the frames. In the
    above, that's:

        >>> t_i = [1, 2, 2, 3, 4]

    or equivalently:

        >>> start_t_i = [1, 2, 2, 3]

    and

        >>> end_t_i = [2, 2, 3, 4]

    Now we simply notice that the *last* start time that searchsorts to being
    before a given frame is going to dictate the state that we observe at that
    frame. So if the states after each ``start_t_i`` were:

        >>> states = ['A', 'B', 'A', 'B']

    then we will observe states ['A', 'A', 'B'] at frames [1, 2, 3], with no
    observation at frames 0 or 4, 5 and 6.

    Finally, sometimes there will be multiple state switches between two
    frames, and yet by chance we end up in the same state as when we started.
    We must combine these.

    To get the correct output for *every* frame, we would simply have to notice
    that sometimes a waiting time observation will cover many frames, and add
    one last step where we simply fill frames not specifically marked by a
    start time with whatever the state was at the previous frame. We don't do
    this right now because my thesis is due tomorrow.
    """
    sim = simulation

    # First, get the discrete window

    left_lext = sim.groupby(traj_cols)['start_time'].min()
    start_i = np.searchsorted(t, left_lext)
    # start times that searchsort past the end means that the entire trajectory
    # occured *after* all frames were over, so it's unobservable
    t_rnan = np.append(t, np.nan)
    leftmost_frame = pd.Series(t_rnan[start_i], index=left_lext.index)
    right_rext = sim.groupby(traj_cols)['end_time'].max()
    # normally would be " - 1", but we insert nan at start of t instead
    # to automagically catch the case that the rightmost wait end searchsorts
    # left of *all* frames, in which case the trajectory is unobservable
    end_i = np.searchsorted(t, right_rext)
    t_lnan = np.insert(t, 0, np.nan)
    rightmost_frame = pd.Series(t_lnan[end_i], index=right_rext.index)

    # Second, build the frames df

    # exclude useless waits
    sim = sim[sim['start_time'] < t[-1]]
    sim = sim[sim['end_time'] > t[0]]
    states = sim[traj_cols + ['state']].copy()
    start_times = sim['start_time'].values
    states['start_i'] = np.searchsorted(t, start_times)
    frames = states.groupby(traj_cols + ['start_i'])['state'].last()
    frames = frames.reset_index().set_index(traj_cols)
    frames['start_time'] = t[frames['start_i'].values]
    frames['window_start'] = leftmost_frame
    frames['window_end'] = rightmost_frame
    # Before adding the end times we need to remove neighboring wait
    # times that end up corresponding to the same state.
    # HACK: create a special index that tracks which unique trajectory each
    # wait time belongs to, then blindly diff state col, and use the new
    # traj_id col to tell when a lack of diff is due to a redundant wait time
    # (same state as prev one) and when it's just due to a new trajectory.
    traj_id = frames.groupby(traj_cols)['start_time'].first() * 0
    traj_id.iloc[:] = np.arange(len(traj_id))
    frames['traj_id'] = traj_id
    # Need to add a single element start to get size match. First wait time by
    # definition can't be redundant, so add traj_id which will never matches so
    # it can't get marked as _is_redundant.
    frames['traj_id2'] = np.insert(frames['traj_id'].values[:-1], 0, -1)
    states = frames['state'].unique()
    state_id = np.arange(len(states))
    state_map = {s: state_id[i] for i, s in enumerate(states)}
    frames['state_id'] = frames['state'].replace(state_map)
    # insert, so the *second* of the redundant guys is deleted
    state_diff = np.insert(np.diff(frames['state_id']), 0, 1)
    redundant = (state_diff == 0) & (frames['traj_id'] == frames['traj_id2'])
    frames = frames.loc[~redundant]
    # Use the same hack as for the redundant guys, but we need the nan on the
    # other side this time, to match what we'll do to the end time col
    traj_id = frames.groupby(traj_cols)['start_time'].first() * 0
    traj_id.iloc[:] = np.arange(len(traj_id))
    frames['traj_id'] = traj_id
    frames['traj_id2'] = np.append(frames['traj_id'].values[1:], np.nan)
    frames['end_time'] = np.append(frames['start_time'].values[1:], np.nan)
    is_last = frames['traj_id'] != frames['traj_id2']
    frames.loc[is_last, 'end_time'] = frames.loc[is_last, 'window_end']
    # get rid of hacky temporary columns
    del frames['traj_id']
    del frames['traj_id2']

    return frames

def traj_to_movie(*args, **kwargs):
    "Alias of :func:`.state_changes_to_movie_frames`."
    return state_changes_to_movie_frames(*args, **kwargs)


def discrete_trajectory_to_wait_times(data, t_col='t', state_col='state'):
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
        columns are ['wait_time', 'start_time', 'end_time', 'state',
        'wait_type', 'min_waits', 'max_waits'], where [wait,end,start]_time
        columns are self explanatory, state is the value of the states_column
        during that waiting time, and wait_type is one of 'interior', 'left
        exterior', 'right exterior', 'full exterior', depending on what kind of
        waiting time was observed. See the `Notes` section below for detailed
        explanation of these categories. The 'min/max_waits' columns contain
        the minimum/maximum possible value of the wait time (resp.), given the
        observations.

        The default index is named "rank_order", since it tracks the order
        (zero-indexed) in which the wait times occured.

    Notes
    -----

    the following types of wait times are of interest to us

    1) *interior* censored times: whenever you are observing a switching
    process for a finite amount of time, any waiting time you observe the
    entirety of is called "interior" censored

    2) *left exterior* censored times: whenever the waiting time you observe
    started before you began observation (it overlaps the "left" side of your
    interval of observation)

    3) *right exterior* censored times: same as above, but overlapping the
    "right" side of your interval of observation.

    4) *full exterior* censored times: whenever you observe the existence of a
    single, particular state throughout your entire window of observation.
    """

    states = data[state_col].values
    times = data[t_col].values
    num_measurements = len(data)

    # now iterate through valid part of trajectory to establish wait times
    start_times = []
    end_times = []
    earliest_st = []  # bounds on start time
    latest_st = []
    earliest_et = []  # bounds on end time
    latest_et = []
    wait_state = []
    wait_type = []
    k0 = 0  # index at which current state began
    state = states[k0]
    state_has_changed = False
    for k in range(num_measurements):
        # if no state change, continue
        if states[k] == state:
            continue
        # otherwise, store change
        start_times.append(times[k0])
        end_times.append(times[k])
        wait_state.append(state)
        # bounds on true wait time value
        if k0 == 0:  # left exterior times have exactly determined "start"
            earliest_st.append(times[k0])
        else:
            earliest_st.append(times[k0-1])
        latest_st.append(times[k0])
        earliest_et.append(times[k-1])
        latest_et.append(times[k])
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
    if k0 == 0:  # full exterior times also have exactly determined "start"
        earliest_st.append(times[k0])
    else:
        earliest_st.append(times[k0-1])
    latest_st.append(times[k0])
    # right/full exterior times have exactly determined "end"
    earliest_et.append(times[k])
    latest_et.append(times[k])
    # state type stored specially for final state
    wait_state.append(state)
    if not state_has_changed:
        wait_type.append('full exterior')
    else:
        wait_type.append('right exterior')
    start_times = np.array(start_times)
    end_times = np.array(end_times)
    wait_times = end_times - start_times
    min_waits = np.array(earliest_et) - np.array(latest_st)
    max_waits = np.array(latest_et) - np.array(earliest_st)
    df = pd.DataFrame({'start_time': start_times, 'end_time': end_times,
                       'wait_time': wait_times, 'state': wait_state,
                       'min_waits': min_waits, 'max_waits': max_waits,
                       'wait_type': wait_type})
    df.index.name = 'rank_order'
    df['window_size'] = times[-1] - times[0]
    return df


def movie_to_waits(*args, **kwargs):
    """Alias of :func:`.discrete_trajectory_to_wait_times`"""
    return discrete_trajectory_to_wait_times(*args, **kwargs)
