"""For plotting multi_locus_analysis results"""
from .analytical import vvc_normalized_theory
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def cvv_plot_sized(cvvs, analytical_deltas=[], delta_col='delta', t_col='t',
                   cvv_col='cvv_normed', max_t_over_delta=4, data_deltas=None,
                   A=1, beta=0.5, tDeltaN=None, cmap_name='viridis', fig=None,
                   alpha_map=None, size_map=None, theory_linewidth=2,
                   include_errorbar=True, include_lines=False,
                   include_points=False, data_line_alpha=1, data_linewidth=1,
                   BASE_DOT_SIZE=10000, **kwargs):
    """One pretty version of the Velocity-Velocity correlation plots for the
    experimental data, with some theory overlaid."""
    # set up data/plot
    if fig is None:
        fig = plt.figure(figsize=(10,10))
    if len(analytical_deltas) > 0 and tDeltaN is None:
        raise ValueError("Specify tDeltaN (and A/beta) if you want to plot analytical curves.")
    if data_deltas is None:
        data_deltas = np.sort(cvvs[delta_col].unique())
    # prune data if requested
    good_ix = np.isin(cvvs[delta_col], data_deltas)
    t_over_delta = cvvs[t_col]/cvvs[delta_col]
    good_ix = good_ix & (t_over_delta <= max_t_over_delta)
    cvvs = cvvs[good_ix]
    # "aesthetic" marker size and alpha values to make sparser stuff
    # (i.e. small delta) equally visible
    max_delta = np.max(cvvs[delta_col])
    cmap = plt.get_cmap(cmap_name)
    cnorm = mpl.colors.Normalize(vmin=0, vmax=max_delta)
    if alpha_map is None:
        alpha_map = lambda x: 0.4 + 0.3*(1-x/max_delta)**2
    if size_map is None:
        MAGIC_NUM = BASE_DOT_SIZE # chosen to make aesthetic dot sizes
        size_map = lambda x: MAGIC_NUM/x
    colors = cmap(cvvs[delta_col]/np.max(cvvs[delta_col]))
    colors[:,3] = alpha_map(cvvs[delta_col])
    if include_points:
        plt.scatter(cvvs[t_col]/cvvs[delta_col],
                    cvvs[cvv_col]/A,
                    color=colors,
                    s=size_map(cvvs[delta_col]), **kwargs)
    if include_errorbar:
        for delta,data in cvvs.groupby(delta_col):
            td = data['t']/data['delta']
            i = np.argsort(td)
            plt.errorbar(td.iloc[i], data['cvv_normed'].iloc[i],
                    data['ste_normed'].iloc[i], c=cmap(cnorm(delta)),
                         linestyle='',
                         alpha=data_line_alpha, **kwargs)
    if include_lines:
        for delta,data in cvvs.groupby(delta_col):
            color = cmap(delta/max_delta)[0:3] + (data_line_alpha, )
            x = data['t']/delta
            y = np.power(delta, 2 - beta)*data[cvv_col]/A
            i = np.argsort(x)
            plt.plot(x.iloc[i], y.iloc[i], c=color,
                     linewidth=data_linewidth, **kwargs)
    plt.xlim([0, max_t_over_delta])
    sm = mpl.cm.ScalarMappable(cnorm, cmap)
    sm.set_array([])
    cbar = plt.colorbar(sm)
    cbar.set_label(r'$\delta$')
    plt.xlabel(r'$t/\delta$')
    plt.ylabel(r'Normalized Autocorrelation\n$<V_{ij}^\delta(t) \cdot{} V_{ij}^\delta(0)>/<V_{ij}^\delta(0)>$')

    # plot theoretical fit curves
    t = np.arange(0.0, 4, 0.001)
    # t = time lag of the correlation
    # delta = step size used to calculate the velocities
    # beta = alpha/2
    # A = "diffusivity" i.e. MSD(delta) = A*delta^alpha
    # tDeltaN = predicted stress correlation time
    for delta in analytical_deltas:
        plt.plot(t, vvc_normalized_theory(t, delta=delta, beta=beta, A=A, tDeltaN=tDeltaN),
                 color=cmap(cnorm(delta)), linewidth=theory_linewidth)
    #     x = np.unique(cvvs[cvvs[:,0] == delta, 1])/delta
    #     x = x[x <= 4]
    #     y = testfunc2(x, delta=delta, beta=beta, A=A, tDeltaN=tdeltaN)
    #     plt.scatter(x, y, color=color, marker='x')


def make_all_disps_hist(displacements, centering="mean,std",
                        cmap=None, reverse=False, alpha=1,
                        cmap_log=False, factor_by=None,
                        xlim=None, ylim=None,
                        include_theory_data=False, include_theory_exact=True,
                        vxcol='vx', vycol='vy', max_plots=np.inf,
                        laplace=True, normal=True, axs=None,
                        no_tick_labels=False, cbar=True,
                        xbins=None, frames_to_sec=None, yscale='log',
                        omit_title=True):
    """Make a manual factor plot of displacement histograms.

    Parameters
    ----------
    factor_by : List<str>
        list of level numbers into the displacements heirarchical index telling
        which set of levels to factor on.
    centering : str
        controls how to center the data so that it fits on a single plot. one
        of "none", "mean", "mean,std" for nothing, mean subtraction, and mean
        subtraction+dividing by std deviation (respectively).
    reverse : bool
        True plots smallest deltas on top, False the opposite
    """
    if xlim is None and centering != None:
        xlim = [-7.5, 7.5] if yscale == 'log' else [-4, 4]
    if ylim is None and yscale == 'log':
        ylim=[0.00005, 1]
    if ylim is None and yscale == 'linear':
        ylim=[0, 0.5]
    frames_multiplier = 1 if frames_to_sec is None else frames_to_sec
    def scale_only(disps):
        X = disps[vxcol]
        Xmean = np.nanmean(X)
        X -= Xmean
        X /= np.nanstd(X)
        X += Xmean
        return X
    def center_and_scale(disps):
        X = disps[vxcol]
        X -= np.nanmean(X)
        X /= np.nanstd(X)
        return X
    def no_centering(disps):
        return disps[vxcol]
    make_disp_plottable = {"none": no_centering, "mean,std": center_and_scale,
                           "std": scale_only}[centering]
    if cmap is None:
        cmap = plt.get_cmap('viridis')
    if factor_by is None:
        # factor by experiment name by default
        factor_by = displacements.index.names.index('experiment')
    if xbins is None:
        xbins = np.linspace(-6, 6, 100)
    if axs is None:
        axs = []
    ys = []
    deltas = []
    # for experiment in displacements.index.levels[displacements.index.names.index('experiment')]:
    #     to_plot = displacements.loc[experiment].reset_index()
    for num_plots, (group_name, to_plot) in enumerate(displacements.groupby(level=factor_by)):
        if num_plots >= max_plots:
            break
        # to_plot = to_plot.reset_index()
        max_dt = to_plot['delta_abs'].max()
        min_dt = to_plot['delta_abs'].min()
        if cmap_log is True:
            cnorm = mpl.colors.LogNorm(vmin=min_dt, vmax=max_dt)
        else:
            cnorm = mpl.colors.Normalize(vmin=min_dt, vmax=max_dt)
        if len(axs) <= num_plots:
            fig, ax = plt.subplots()
            axs.append(ax)
        else:
            ax = axs[num_plots]
            plt.sca(ax)
        if reverse:
            hard_reverse = [(dt, disps) for dt, disps in to_plot.groupby('delta_abs')]
            maybe_reverse = reversed(hard_reverse)
        else:
            maybe_reverse = to_plot.groupby('delta_abs')
        for dt, disps in maybe_reverse:
            X = make_disp_plottable(disps)
            X = X.loc[np.isfinite(X)]
            if not np.any(np.isfinite(X)) or len(X) < 3:
                continue
            delta = dt
            y, _, _ = ax.hist(X, bins=xbins, **{"histtype": "step",
                    "color": cmap(cnorm(dt)), "normed": True,
                    "linewidth": 1, "alpha": alpha})
            deltas.append(delta)
            ys.append(y)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=cnorm)
        # fake up the array of the scalar mappable. Urgh...
        sm._A = []
        if cbar:
            cbar = plt.colorbar(sm)
            if frames_to_sec is None:
                cbar.ax.set_ylabel(r'$\delta$ (frames)')
            else:
                cbar.ax.set_ylabel(r'$\delta$ (s)')
        if include_theory_data:
            if normal:
                ax.hist(np.random.standard_normal((to_plot.groupby('delta_abs').count().max()[vxcol],1)),
                        bins=xbins, **{"histtype": "step", "linewidth": 2, "color": "k",
                        "alpha": 1, "normed": True})
            if laplace:
                ax.hist(np.random.laplace(size=(to_plot.groupby('delta_abs').count().max()[vxcol],1)),
                        bins=xbins, **{"histtype": "step", "linewidth": 2, "color": "k",
                        "alpha": 1, "normed": True})
        if include_theory_exact:
            xbins[0] *= 10
            xbins[-1]*= 10
            if normal:
                ax.plot(xbins, scipy.stats.norm.pdf(xbins), 'k-.')
            if laplace:
                ax.plot(xbins, scipy.stats.laplace.pdf(xbins), 'k-.')
        plt.xlabel(r'$v_x$')
        if xlim is not None:
            plt.xlim(xlim)
        if ylim is not None:
            plt.ylim(ylim)
        plt.yscale(yscale)
        title_str = str(group_name).replace(r'_', r'\_') # escape LaTeX as necessary
        title_str += ' Displacement Distribution'
        if include_theory_data or include_theory_exact:
            title_str += '\n'
            if laplace:
                title_str += 'Laplace '
                if normal:
                    title_str += '+ '
                else:
                    title_str += 'Distribution '
            if normal:
                title_str += 'Standard Normal '
            title_str += 'for scale'
        if not omit_title:
            plt.title(title_str)
        if no_tick_labels:
            # labels = [item.get_text() for item in ax.get_xticklabels()]
            # empty_string_labels = ['']*len(labels)
            # ax.set_xticklabels(empty_string_labels)
            labels = [item.get_text() for item in ax.get_yticklabels()]
            empty_string_labels = ['']*len(labels)
            ax.set_yticklabels(empty_string_labels)
    return axs, deltas, ys
