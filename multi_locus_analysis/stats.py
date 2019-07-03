r"""For computing statistics of particle motions

When dealing with trajectory data, all the common readouts can be thought of as
moments of the trajectory with different time lags or discrete derivatives
applied.

More precisely, suppose you have trajectory data :math:`X_i(t_k)`, then
we define the discrete velocity (the first discrete derivative of the particle
as :math:`V_i^\delta(t_k) = X_i(t_k + \delta) - X_i(t_k)`
and the discrete relative velocities
:math:`V_{ij}^\delta(t_k) = X_{ij}(t_k + \delta) - X_{ij}(t_k)`
where :math:`X_{ij}` is the distance between the points :math:`i` and :math:`j`.
This module simply helps quickly calculate moments of the :math:`V^\delta(t)`.

For example, the ensemble MSD is simply :math:`E_i[(V_i^\delta(0))^2]` as a
function of :math:`\delta`, where :math:`E_\chi[\cdot]` indicates the average is
taken over the variable :math:`\chi`. The time-averaged MSD for a single
particle on the other hand is given by :math:`E_t[(V_i^\delta(t))^2]`.

This module provides a simple set of functions for easily extracting the moment
of any product of any derivatives with any combination of time lags of a set of
trajectories, :math:`X_i(t_k)`.
"""
import numpy as np
import pandas as pd
import scipy
from scipy.io import loadmat

import gc
from pathlib import Path

###############{{{
# Call on DataFrame directly

def re_agg_mean(df, group_cols, mean_col='mean', count_col='count'):
    def ave_wrap(df):
        return np.average(df[mean_col], weights=df[count_col])
    return df.groupby(group_cols).apply(ave_wrap)

def autocorr_1(df, xcol, tcol=None, **kwargs):
    XX = corr_prod(df, xcol=xcol, tcol=tcol, **kwargs)
    return XX.groupby(['t'])['XX'].agg(['mean', 'std', 'count'])

def corr_prod(df, xcol, tcol=None, max_dt=None, max_t_over_delta=4,
              allowed_dts=None, t_step=None):
    r"""Calculates the time averaged autocorraltion correlation over a column,
    with an optional time column.

    The time average can be written
    .. math::

        C_X(t) = E_\tau[X_{\tau + t}X_\tau]
    """
    x = df[xcol]
    if tcol is not None:
        t = df[tcol]
    else:
        t = np.arange(len(x))
    # C_x(i-j) = dx[i]*dx[j] so all n^2 products
    XX = x[None,:]*x[:,None]
    dts = t[None,:] - t[:,None]
    tau = t[None,:] - np.zeros_like(t[:,None])
    good_ix = dts >= 0
    if max_dt:
        good_ix &= dts < max_dt
    if max_t_over_delta and 'delta' in df:
        good_ix &= dts/df['delta'].iloc[0] <= max_t_over_delta
    if allowed_dts:
        good_ix &= np.isin(dts, allowed_dts)
    arrs = [pd.Series(arr[good_ix]) for arr in [dts, XX, tau]]
    XX = pd.concat(arrs, axis=1)
    XX.columns = ['t', 'XX', 'tau']
    XX.set_index(['t', 'tau'], inplace=True)
    return XX

# end Call on DataFrame directly
###############}}}


###############{{{
# Apply to groupby'd DataFrame

def traj_to_msds(traj, xcol='x', ycol='y', framecol='frame id'):
    """Data should be groupby'd trajectories.
    Takes x,y coordinate columns and the "time" column, expects integer index
    for time column."""
    x = traj[xcol]
    y = traj[ycol]
    ft = traj[framecol]
    dfts = ft[None,:] - ft[:,None]
    dxs = x[None,:] - x[:,None]
    dys = y[None,:] - y[:,None]
    dX = np.power(dxs, 2) + np.power(dys, 2)
    old_index = np.array(traj.index)
    old_index_start = old_index[:,None] - np.zeros((len(traj.index),))[None,:].astype(int)
    old_index_end = old_index[None,:] - np.zeros((len(traj.index),))[:,None].astype(int)
    pos_ix = dfts > 0
    displacements = pd.DataFrame({'displacements': dX[pos_ix],
                               'old_index_end': old_index_end[pos_ix],
                               'old_index_start': old_index_start[pos_ix],
                               'delta': dfts[pos_ix]})
    return displacements

def pos_to_all_vel(trace, xcol='x', ycol='y', zcol=None, framecol='i',
                   delta=None, delta_max=None, deltas=None,
                   absolute_time=None, force_independence=False):
    """For getting displacement distributions for all possible deltas
    simultaneously.

    >>> all_vel = df.groupby(bp.track_columns).apply(lambda df:
                mla.pos_to_all_vel(df, xcol='dX', ycol='dY', framecol='tp.n'
        ))

    """
    # framecol could feasibly be an index column of the dataframe, so we should
    # reset index just in case
    trace = trace.reset_index()
    t = trace[framecol]
    x = trace[xcol]
    y = trace[ycol]
    if zcol is not None:
        z = trace[zcol]
    ti = t[:,None] - np.zeros_like(t)[None,:]
    tf = t[None,:] - np.zeros_like(t)[:,None]
    dt = tf - ti
    vx = x[None,:] - x[:,None]
    vy = y[None,:] - y[:,None]
    if zcol is not None:
        vz = z[None,:] - z[:,None]
    good_dt = dt >= 0
    if force_independence:
        i = np.arange(len(ti))
        overlap_id = i[None,:] % (i[None,:] - i[:,None])
        good_dt &= overlap_id == 0
    if delta:
        good_dt &= dt == delta
    if delta_max:
        good_dt &= dt < delta_max
    if deltas is not None:
        good_dt &= np.isin(dt, deltas)
    all_vel = {'ti': ti[good_dt], 'delta': dt[good_dt], 'tf': tf[good_dt],
               'vx': vx[good_dt], 'vy': vy[good_dt]}
    if zcol is not None:
        all_vel['vz'] = vz[good_dt]
    all_vel = pd.DataFrame(all_vel)
    if absolute_time:
        frames_to_sec = trace['frames_to_sec'].unique()
        if len(frames_to_sec) != 1:
            raise ValueError("You asked me to convert to absolute time, but there's not one unique frame_to_sec value")
        all_vel['delta_abs'] = all_vel['delta']*frames_to_sec[0]
    all_vel.set_index(['ti', 'delta'], inplace=True)
    return all_vel

def vel_to_pos(vel, dxcol='vx', dycol='vy', framecol='tf'):
    t = vel[framecol].as_matrix()
    vx = vel[dxcol].as_matrix()
    vy = vel[dycol].as_matrix()
    sorted_ix = np.argsort(t)
    df = pd.DataFrame({'x': np.cumsum(vx[sorted_ix]),
                       'y': np.cumsum(vy[sorted_ix])},
                      index=t[sorted_ix])
    df.index.name = 't'
    return df

def all_vel_to_corr(vel, dxcol='vx', dycol='vy', framecol='tf',
        max_dt=None, max_t_over_delta=4, allowed_dts=None):
    """
    >>> all_vel.reset_index(level=pandas_util.multiindex_col_ix(all_vel, 'ti'), inplace=True)
    """
    t = vel[framecol]
    dx = vel[dxcol]
    dy = vel[dycol]
    # cvv[i,j] = dx[i]*dx[j] + dy[i]*dy[j], so all n^2 dot products
    cvv = dx[None,:]*dx[:,None] + dy[None,:]*dy[:,None]
    tau = t[None,:] - np.zeros_like(t[:,None])
    dts = t[None,:] - t[:,None]
    good_ix = dts >= 0
    if max_dt:
        good_ix &= dts < max_dt
    if max_t_over_delta and 'delta' in vel:
        good_ix &= dts/vel['delta'].iloc[0] <= max_t_over_delta
    if allowed_dts:
        good_ix &= np.isin(dts, allowed_dts)
    arrs = [pd.Series(arr[good_ix]) for arr in [dts, cvv, tau]]
    cvvs = pd.concat(arrs, axis=1)
    cvvs.columns = ['t', 'cvv', 'tau']
    cvvs.set_index(['t', 'tau'], inplace=True)
    return cvvs

def moments(df, cols=None, ns=[0,1,2]):
    """Take a DataFrame and optionally a list of columns of interest and for
    each n in ns, calculate the nth moment of each column, where the 0th moment
    in defined for convenience to be the count.

    Returns a series with multi-index with two levels. Level 0 contains the
    requested columns and level 1 contains the requested moment numbers (ns)."""
    if cols:
        df = df[cols]
    def moment(n):
        if n == 0:
            return lambda x: np.nansum(np.isfinite(x))
        elif n > 0:
            return lambda x: np.nanmean(np.power(x, n))
    # rows have cols, columns have ns
    df = pd.DataFrame({n: df.apply(moment(n), axis=0) for n in ns})
    df.rename_axis('Moment', axis='columns', inplace=True)
    return df.T.unstack()

# end Apply to groupby'd DataFrame
###############}}}

###############{{{
# "combining" functions, post groupby/apply

def combine_moments(data):
    """Takes a column from output of df.groupby([...]).apply(moments) and
    combines the values to get overall moments for full dataset.

    Example:

    .. code-block:: python

        grouped_stats = df.groupby(['experiment']).apply(lp.moments)
        overall_stats = grouped_stats.groupby(level=0, axis=1).apply(lp.combine_moments).unstack()

    Example:

    .. code-block:: python

        def grouper(g):
            return g.groupby(level=0, axis=1).apply(lp.combine_moments).unstack()
        grouped_stats[new_col] = make_interesting_value(grouped_stats)
        overall_stats = df.groupby(new_col).apply(grouper)
    """
    # get top level column multi-index name if it exists
    if data.columns.names[0] != 'Moment':
        name = data.columns.levels[0][data.columns.labels[0][0]]
        data = data[name]
    total = np.nansum(data[0])
    weights = data[0]/total
    combined_vals = data.apply(lambda x: np.nansum(weights*x))
    combined_vals[0] = total
    return combined_vals

def cvv_by_hand_make_usable(cvv_stats, group_cols):
    """add columns to vvc_stats_by_hand output that we typically want. (cvv,
    std, ste, cvv_normed, ste_normed). requires the group-by columns used to do
    the normalization via msds using groupby"""
    cvv_stats['cvv'] = cvv_stats['sum']/cvv_stats['cnt']
    cvv_stats['std'] = np.sqrt(cvv_stats['sqs']/cvv_stats['cnt']
                                - np.power(cvv_stats['cvv'], 2))
    cvv_stats['ste'] = cvv_stats['std']/np.sqrt(cvv_stats['cnt'])
    def subtract_t0(cvvs):
        cvv0 = cvvs[cvvs['t'] == 0]['cvv']
        if len(cvv0) > 1:
            raise ValueError('Ambiguous t0, did you groupby the wrong thing?')
        cvvs['cvv_normed'] =  cvvs['cvv']/cvv0.iloc[0]
        cvvs['ste_normed'] = cvvs['ste']/cvv0.iloc[0]
        return cvvs
    cvv_stats = cvv_stats.groupby(group_cols + ['delta']).apply(subtract_t0)
    return cvv_stats

# "combining" functions, post groupby/apply
###############}}}


###############{{{
# Doing Pandas stuff manually cause memory leaks

def groupby_apply_efficient(df, group_cols, apply_fun, n_between_gc=50,
                            print_progress=False, apply_to_cols=None):
    """An attempt to make groupby not leak memory. Kinda worked at one point
    (late 2016)."""
    output_df = []
    for i, (group_vals, group) in enumerate(df.groupby(group_cols)):
        if apply_to_cols:
            group = group[apply_to_cols]
        if i % n_between_gc == 0:
            gc.collect()
            if print_progress:
                print('.')
        group_df = apply_fun(group).reset_index()
        for j, col in enumerate(group_cols):
            group_df[col] = group_vals[j]
        output_df.append(group_df)
    return pd.concat(output_df, ignore_index=True)

def vels_to_cvvs_by_hand(vels, groups, file_name, framecol='ti', dxcol='vx',
                            dycol='vy', dzcol='vz', max_t_over_delta=None, max_dt=None,
                            allowed_dts=None, deltas_of_interest=None,
                            include_group=True):
    """should be passed a velocities dataframe (as from pos_to_all_vel.

    groups vels by groups, optionally ignores group labels to make a flat
    file with all vvc information with repeats of t,delta pairs as
    appropriate."""
    # plausibly, vels has some of groups in its index, so we will want to make
    # them columns, espcially delta and framecol
    vels = vels.reset_index()
    if deltas_of_interest is None:
        deltas_of_interest = vels.delta.unique()
    # delta == 0, then all velocities equal zero, wasting computation time
    deltas_of_interest = deltas_of_interest[deltas_of_interest > 0]
    if 'delta' not in groups:
        groups.append('delta')
    with open(file_name, 'w') as f:
        f.write('tau,t,cvv')
        for group in groups:
            f.write(',' + str(group))
        f.write('\n')
        for label, vel in vels[
            np.isin(vels.delta, deltas_of_interest)
        ].groupby(groups):
            t = vel[framecol]
            dx = vel[dxcol]
            dy = vel[dycol]
            if dzcol is not None:
                dz = vel[dzcol]
                cvv = dx[None,:]*dx[:,None] + dy[None,:]*dy[:,None] + dz[None,:]*dz[:,None]
            # cvv[i,j] = dx[i]*dx[j] + dy[i]*dy[j], so all n^2 dot products
            else:
                cvv = dx[None,:]*dx[:,None] + dy[None,:]*dy[:,None]
            tau = t[None,:] - np.zeros_like(t[:,None])
            dts = t[None,:] - t[:,None]
            good_ix = dts >= 0
            if max_dt:
                good_ix &= dts < max_dt
            if max_t_over_delta and 'delta' in vel:
                good_ix &= dts/vel['delta'].iloc[0] <= max_t_over_delta
            if allowed_dts:
                good_ix &= np.isin(dts, allowed_dts)
            # if tmax_or_dtmax:
            #     good_ix &= ((dts < tmax_or_deltamax[0])
            #             | (dts/vel['delta'].iloc[0] <= tmax_or_deltamax[1]))
            arrs = [pd.Series(arr[good_ix]) for arr in [tau, dts, cvv]]
            cvvs = pd.concat(arrs, axis=1)
            cvvs.columns = ['tau', 't', 'cvv']
            cvvs.dropna(inplace=True)
            if include_group:
                for i, group in enumerate(groups):
                    cvvs[group] = label[i]
            cvvs.to_csv(f, index=False, header=False)

def vvc_stats_by_hand(file_name, groups, print_interval=None, skip_lines=None):
    """Calculates moments using the output of vels_to_cvvs_by_hand."""
    grouped_sqs = {}
    grouped_sum = {}
    grouped_cnt = {}
    group_ixs = []
    group_names = []
    with open(file_name) as f:
        for i,line in enumerate(f):
            # get group names on first go round
            if i == 0:
                columns = line.rstrip().split(',')
                for j, column in enumerate(columns):
                    if column == 'cvv':
                        cvv_ix = j
                    if column in groups or column == 'delta' or column == 't':
                        group_ixs.append(j)
                        group_names.append(column)
                group_ixs = np.array(group_ixs)
                group_names = group_names
            # <= accounts for first line being header
            if skip_lines and i <= skip_lines:
                continue
            if print_interval and i % print_interval == 0:
                print(i)
            vals = np.array(line.rstrip().split(','))
            try:
                cvv = float(vals[cvv_ix])
            except ValueError:
#                 cvv = np.nan
                continue
            group = tuple(vals[group_ixs])
            if group not in grouped_sqs:
                grouped_sqs[group] = 0
                grouped_sum[group] = 0
                grouped_cnt[group] = 0
            grouped_sqs[group] += cvv*cvv
            grouped_sum[group] += cvv
            grouped_cnt[group] += 1
    cvvsqs = pd.DataFrame.from_dict(grouped_sqs, orient='index')
    cvv_stats = cvvsqs.reset_index()['index'].apply(pd.Series)
    cvv_stats.columns = group_names
    cvv_stats['sqs'] = cvvsqs[0].values
    cvvsum = pd.DataFrame.from_dict(grouped_sum, orient='index')
    cvv_stats['sum'] = cvvsum[0].values
    cvvcnt = pd.DataFrame.from_dict(grouped_cnt, orient='index')
    cvv_stats['cnt'] = cvvcnt[0].values
    cvv_stats['t'] = pd.to_numeric(cvv_stats['t'])
    cvv_stats['delta'] = pd.to_numeric(cvv_stats['delta'])
    return cvv_stats

# end Doing Pandas stuff manually cause memory leaks
###############}}}


###############{{{
# TODO: not yet made general, but useful statistics
def scale_and_test_normality(vx):
    """Should be applied to a vector of velocities, vx"""
    vx = vx[np.isfinite(vx)]
    vx -= np.mean(vx)
    std = np.std(vx)
    if std > 0:
        vx /= np.std(vx)
    num_disps = vx.size
    ksstat, ks_pval = scipy.stats.kstest(vx, 'norm', mode='asymp')
    ngp = np.mean(np.power(vx, 4))/(3*np.power(np.mean(np.power(vx, 2)), 2)) - 1
    try:
        shapiro_stat, shapiro_pval = scipy.stats.shapiro(vx)
    except:
        shapiro_stat = shapiro_pval = np.nan
    is_shapiro_pval_accurate = num_disps < 5000 # from scipy docs, v 19.1
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
            'Shapiro-Wilk Statistic': shapiro_stat, 'Shapiro-Wilk p-value':
            shapiro_pval, 'Shapiro-Wilk p-value is accurate':
            is_shapiro_pval_accurate, 'Anderson-Darling Statistic':
            anderson_stat, 'Anderson-Darling alpha lower bound': crit_alphas[ilow],
            'Anderson-Darling alpha upper bound': crit_alphas[ihigh],
            'Anderson-Darling upper crit value': crit_vals[ihigh],
            'Anderson-Darling lower crit value': crit_vals[ilow]
            }

def add_savgol(traj, window_size, order=3):
    if len(traj['x']) <= window_size:
        traj['savgol'+str(window_size)+'_x'] = np.nan
        traj['savgol'+str(window_size)+'_y'] = np.nan
        traj['savgol'+str(window_size)+'_x_fluct'] = np.nan
        traj['savgol'+str(window_size)+'_y_fluct'] = np.nan
    else:
        traj['savgol'+str(window_size)+'_x'] = savgol_filter(traj['x'], window_size, order)
        traj['savgol'+str(window_size)+'_y'] = savgol_filter(traj['y'], window_size, order)
        traj['savgol'+str(window_size)+'_x_fluct'] = traj['x'] - traj['savgol'+str(window_size)+'_x']
        traj['savgol'+str(window_size)+'_y_fluct'] = traj['y'] - traj['savgol'+str(window_size)+'_y']
    return traj

def add_savgol_window(data, window_sizes):
    for window_size in window_sizes:
        data = data.groupby(['experiment', 'movie name', 'molecule id']
                            ).apply(partial(add_savgol, window_size=window_size))
    return data

def get_savgol_vels(data, window_sizes):
# we're gonna do the same thing twice, basically get vels for each window size
    savgol_vels = []
    for window_size in window_sizes:
        fluct_vels = data.groupby(['experiment', 'movie name', 'molecule id']).apply(partial(pos_to_all_vel, delta=1, xcol='savgol'+str(window_size)+'_x', ycol='savgol'+str(window_size)+'_y'))
        fluct_vels['window_size'] = window_size
        savgol_vels.append(fluct_vels)
    savgol_vels = [vel.reset_index().set_index(['experiment', 'movie name', 'molecule id', 'window_size', 'ti', 'delta']) for vel in savgol_vels]
    savgol_vels = pd.concat(savgol_vels)
# now for _fluct columns
    corrected_vels = []
    for window_size in window_sizes:
        fluct_vels = data.groupby(['experiment', 'movie name', 'molecule id']).apply(partial(pos_to_all_vel, delta=1, xcol='savgol'+str(window_size)+'_x_fluct', ycol='savgol'+str(window_size)+'_y_fluct'))
        fluct_vels['window_size'] = window_size
        corrected_vels.append(fluct_vels)
    corrected_vels = [vel.reset_index().set_index(['experiment', 'movie name', 'molecule id', 'window_size', 'ti', 'delta']) for vel in corrected_vels]
    return savgol_vels, pd.concat(corrected_vels)

def get_savgol_msds(savgol_data, window_sizes):
    msds = []
    for window_size in window_sizes:
        savgol_disps = savgol_data.groupby(['experiment', 'movie name', 'molecule id']).apply(partial(pos_to_all_vel, absolute_time=True, xcol='savgol'+str(window_size)+'_x_fluct', ycol='savgol'+str(window_size)+'_y_fluct'))
        sq_disps = pd.DataFrame(index=savgol_disps.index)
        sq_disps['sq_disp'] = np.power(savgol_disps['vx'], 2) +np.power(savgol_disps['vy'], 2)
        sq_disps['delta_abs'] = savgol_disps['delta_abs'].round(3)
        savgol_msds = sq_disps.groupby(['delta_abs'])['sq_disp'].agg(['mean', 'std', 'count'])
        savgol_msds['window_size'] = window_size
        savgol_msds = savgol_msds.reset_index().set_index(['window_size', 'delta_abs'])
        msds.append(savgol_msds)
    return pd.concat(msds)

# end TODO: not yet made general, but useful statistics
###############}}}
