import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def plot_colored_line(t, x, y, cmap='viridis', linewidth=3, ax=None,
                      colorbar=True):
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    lc = mpl.collections.LineCollection(segments, cmap=plt.get_cmap(cmap))
    lc.set_array(t)
    lc.set_linewidth(linewidth)

    if ax is None:
        fig, ax = plt.subplots()
    ax.add_collection(lc)
    ax.set_xlim(get_lim(x))
    ax.set_ylim(get_lim(y))

    if colorbar:
        cnorm = mpl.colors.Normalize(vmin=np.min(t), vmax=np.max(t))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=cnorm)
        sm._A = t
        cax = plt.gca()
        plt.sca(ax)
        cbar = plt.colorbar(sm)
        plt.sca(cax)
    return ax

def get_lim(x, margin=0.1):
    min = np.min(x)
    max = np.max(x)
    dx = max - min
    return [min - margin*dx, max + margin*dx]

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
    Inputs:
        factor_by - list of level numbers into the displacements heirarchical index
                    telling which set of levels to factor on.
        centering - controls howto center the data so that it fits on a signle
                    plot
        reverse - True plots smallest deltas on top, False the opposite
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
                cbar.ax.set_ylabel('$\delta$ (frames)')
            else:
                cbar.ax.set_ylabel('$\delta$ (s)')
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
        plt.xlabel('$v_x$')
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
