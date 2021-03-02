"""
Code specialized for plotting test runs of the finite_window code.
"""
from .. import finite_window as fw
from .. import stats

import lifelines
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, PathPatch, Rectangle
from matplotlib.legend_handler import HandlerBase
from matplotlib.path import Path
from matplotlib.transforms import Affine2D
import numpy as np
import pandas as pd
import scipy
import seaborn as sns

import bruno_util.plotting as bplt
figure_size = bplt.use_cell_style(mpl.rcParams)

km_color = sns.color_palette('colorblind')[0]
interior_color = sns.color_palette('colorblind')[1]
interior_linestyle = '-.'


# make cool legends where stair lines and straight lines are treated
# differently
stair_path = Path([
    [-1, -1], [-1/3, -1], [-1/3, 1], [1/3, 1], [1/3, -1], [1, -1], [1, 1],
    [5/3, 1]
])
bbox = stair_path.get_extents()
make_unit = Affine2D() \
        .translate(-bbox.x0, -bbox.y0) \
        .scale(1/bbox.width, 1/bbox.height)
stair_path = make_unit.transform_path(stair_path)


def _int_win_from_obs(obs, state):
    interior = obs.loc[
        (obs['state'] == state) & (obs['wait_type'] == 'interior'),
        ['wait_time', 'window_size']
    ]
    return interior.wait_time.values, interior.window_size.values


def _ext_from_obs(obs, state):
    exterior = obs.loc[
            (obs['state'] == state)
            & (obs['wait_type'] != 'interior')
            & (obs['wait_type'] != 'full exterior'),
            ['wait_time']
    ]
    return exterior.wait_time.values


def _get_axes(var_pair, name=None, size=None):
    if size is None:
        if name is None:
            raise ValueError("Must specify figure size or figure_size[] key!")
        size = figure_size[name]
    fig = plt.figure(figsize=size, constrained_layout=True)
    axs = fig.subplot_mosaic([[var.name for var in var_pair]])
    return fig, axs


def _init_entries(var_pair):
    """
    Make base for smart legend handler accumulator.

    We make our column-per-variable legends by using a fake patch to store the
    label that we want for the top of the column, then manually shifting the
    label left to center it over the column after all the other handles have
    been drawn.

    Add handles to each of the ``legend_entries[var.name] : List[Patch]``.
    """
    return {
        var.name: [FakePatch(alpha=0, label=var.pretty_name)]
        for var in var_pair
    }


def _legend_from_entries(ax, legend_entries, label_hspace=-4.8, **kwargs):
    handles = [h for _, patches in legend_entries.items() for h in patches]
    legend = ax.legend(handles=handles, ncol=2, **kwargs, handler_map={
        FakePatch: FakeHandlerPatch(), StairLine: StairLineHandler()
    })

    # hack to left-align my "fake" legend column "titles"
    # first remove unpredictable "packing" whitespace
    for vpack in legend._legend_handle_box.get_children():
        for hpack in vpack.get_children()[:1]:
            hpack.get_children()[0].set_width(0)
    legend_texts = legend.get_texts()
    # then scoot horizontally by label_hspace points
    for i, handle in enumerate(legend.legendHandles):
        if isinstance(handle, FakeRectangle):
            legend_texts[i].set_position(
                (label_hspace/72*mpl.rcParams['figure.dpi'], 0)
            )


class FakePatch(Patch):
    """Identify legend handles for special treatment."""
    pass


class FakeRectangle(Rectangle):
    """Identify post-handler legendHandles for special treatment."""
    pass


class FakeHandlerPatch(HandlerBase):
    """
    Handler `.FakePatch` instances. Ensure `.FakeRectangle` used.
    """
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        p = FakeRectangle(xy=(-xdescent, -ydescent),
                          width=width, height=height)
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)
        return [p]


class StairLine(object):
    def __init__(self, line):
        self.line = line

    def get_label(self):
        return self.line._label

    def get_color(self):
        return self.line.get_color()

    def get_linewidth(self):
        return self.line.get_linewidth()

    def get_alpha(self):
        return self.line.get_alpha()


class StairLineHandler(object):

    def __init__(self, **kwargs):
        self.path_patch_kw = kwargs

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        put_in_leg = Affine2D().translate(x0, y0).scale(width, height)
        leg_patch = put_in_leg.transform_path(stair_path)
        patch = PathPatch(leg_patch, color=orig_handle.get_color(),
                          alpha=orig_handle.get_alpha(),
                          transform=handlebox.get_transform(),
                          lw=orig_handle.get_linewidth(),
                          fill=False,
                          **self.path_patch_kw)
        handlebox.add_artist(patch)
        return [patch]


class Variable:
    """
    Wrap scipy rv's so they're easy to plot.
    """

    # which kwargs to init have defaults that change every time you instantiate
    _cyclers = {
        'linestyle': ['-', '--', '-.', ':'],
        'color': sns.color_palette('colorblind')
    }
    _cyclers_i = {cycler: 0 for cycler in _cyclers}

    def scaled_ylim(self, T):
        if self._scaled_ylim is not None:
            return self._scaled_ylim
        t = np.linspace(0, T, 101)
        return 1.1/self.rv.cdf(T) * np.max(self.rv.pdf(t))

    @staticmethod
    def _get_rv_name(rv):
        name = rv.dist.name
        if len(rv.args) > 0:
            name += '(' + str(rv.args[0])
            for arg in rv.args[1:]:
                name += ', ' + arg
            name += ')'
        return name[0].upper() + name[1:]

    def __init__(self, rv, **kwargs):
        """
        Parameters
        ----------
        rv : scipy.stats.rv_frozen
            The random variable to be sampled.
        name : Optional[str]
            The name to be used as the dataframe index for this variable (can
            be autogenerated from *rv*).
        pretty_name : Optional[str]
            The name to be used when creating plot labels. Same as *name* if
            not specified.
        linestyle : Optional[matplotlib.lines.LineStyle-like]
            The linestyle to be used to identify this variable (if one is
            used). Cycles through ``['-', '--', '-.', ':']``.
        color : Optional[matplotlib.color.Color-like]
            The color to be used to identify this variable (if one is used).
            Cycles through ``sns.color_palette('colorblind')``.
        """
        self.rv = rv

        if "name" not in kwargs:
            kwargs['name'] = Variable._get_rv_name(rv)
        if "pretty_name" not in kwargs:
            kwargs["pretty_name"] = kwargs['name']
        self._scaled_ylim = kwargs.pop("scaled_ylim", None)

        for cycler, options in Variable._cyclers.items():
            if cycler not in kwargs:
                i = Variable._cyclers_i[cycler]
                kwargs[cycler] = options[i]
                Variable._cyclers_i[cycler] += 1
                Variable._cyclers_i[cycler] %= len(options)

        self.__dict__.update({
            key: getattr(rv, key) for key in dir(rv)
        })
        self.__dict__.update(kwargs)


# so the user doesn't have to actually know how to construct these to use the
# code, we make reasonable defaults. marking name as "None" will be hard-coded
# to rely instead on the dataframe's contents below
_default_vars = [
    Variable(rv=None, name=None, linestyle=':',
             color=sns.color_palette('colorblind')[2]),
    Variable(rv=None, name=None, linestyle='--',
             color=sns.color_palette('colorblind')[3]),
]


def compare_interior_kaplan(obs, var_pair, rescale_kaplan=False,
                            rescale_interior=False):
    """
    Interior vs kaplan est for `multi_locus_analysis.finite_window.ab_window`.

    Compare the Kaplan-Meier estimator to the empirical distribution function
    (eCDF) of interior times of data generated using the
    `multi_locus_analysis.finite_window.ab_window` or
    `multi_locus_analysis.finite_window.ab_window_fast` functions.
    """
    kmfs = {}
    for name, state in obs.groupby('state'):
        times = state['wait_time'].values
        not_censored = (state['wait_type'] == 'interior').values
        kmfs[name] = lifelines.KaplanMeierFitter().fit(
            times, event_observed=not_censored,
            label=r'Meier-Kaplan Estimator, $\pm$95% conf int'
        )

    fig, axs = _get_axes(
        var_pair,
        name='two-by-half column, four legend entries above'
    )
    T = obs.window_size.max()
    for var in var_pair:
        ax = axs[var.name]

        # extract KM CDF fit
        tk = kmfs[var.name].cumulative_density_.index.values
        kmf = kmfs[var.name].cumulative_density_.values
        # and confidence intervals
        low, high = kmfs[var.name] \
            .confidence_interval_cumulative_density_.values.T
        Z = kmf[-1] / var.cdf(T) if rescale_kaplan else 1
        km_l = ax.plot(tk, kmf/Z, color=km_color, label='Kaplan-Meier')[0]
        ax.fill_between(tk, low/Z, high/Z, color=km_color, alpha=0.4)

        # plot actual distribution
        t = np.linspace(0, T, 101)
        analytical_l,  = ax.plot(
            t, var.cdf(t), color='k', label='Actual CDF'
        )

        # now compute the empirical distribution of the "interior" times
        interior, _ = _int_win_from_obs(obs, var.name)
        x, cdf = fw.ecdf(interior, pad_left_at_x=0)

        Z = 1/var.cdf(x[-1]) if rescale_interior else 1
        interior_l, = ax.plot(
            x, cdf / Z, c=var.color, ls=interior_linestyle,
            label='"Interior" eCDF'
        )

        # prettify the plot
        ax.set_xlim([0, T])
        ax.set_ylim([0, 1])
        ax.set_xlabel('time')
        ax.set_ylabel(r'Cumulative probability')

        ax.legend(
            title=var.pretty_name,
            handles=[interior_l, km_l, analytical_l],
            # align bottom of legend 2% ax height above axis, filling full axis
            # width
            bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
            ncol=1, mode="expand", borderaxespad=0.
        )
    return fig


def exterior_hist(obs, var_pair, ext_bins=51):
    fig, ax = plt.subplots(
        figsize=figure_size['full column'],
        constrained_layout=True
    )
    legend_entries = _init_entries(var_pair)
    T = obs.window_size.max()
    t = np.linspace(0, T, 101)
    for var in var_pair:
        exterior = _ext_from_obs(obs, var.name)

        y, t_bins = np.histogram(exterior, bins=ext_bins, density=1)
        X, Y = stats.bars_given_hist(y, t_bins)
        line, = ax.plot(X, Y, c=var.color, label='Exterior Histogram')
        legend_entries[var.name].append(line)

        scale = T - scipy.integrate.quad(var.cdf, 0, T)[0]
        line, = ax.plot(t, var.sf(t)/scale, c='k', ls=var.linestyle,
                        label='Rescaled\nsurvival function')
        legend_entries[var.name].append(line)
    ax.set_xlabel('time')
    ax.set_ylabel(r'$P(X_\mathrm{exterior} = t)$')
    ax.set_xlim([0, T])
    ax.set_ylim([0, ax.get_ylim()[1]])
    _legend_from_entries(ax, legend_entries)


def corrected_interior_pdf(obs, var_pair, traj_cols=['replicate'],
                           int_bins=51):
    fig, ax = plt.subplots(
        figsize=figure_size['full column'],
        constrained_layout=True
    )
    legend_entries = _init_entries(var_pair)

    traj_windows = obs.groupby(traj_cols)['window_size'].first().values
    # now sorted
    traj_windows, window_cdf = stats.ecdf(traj_windows, pad_left_at_x=0)
    window_sf = 1 - window_cdf

    T = obs.window_size.max()
    t = np.linspace(0, T, 101)
    for var in var_pair:
        interior, window_sizes = _int_win_from_obs(obs, var.name)
        window_i = np.searchsorted(traj_windows, interior) - 1
        correction = 1/window_sf[window_i]/(window_sizes - interior)
        y, t_bins = np.histogram(interior, weights=correction, bins=int_bins,
                                 density=1)
        X, Y = stats.bars_given_hist(y, t_bins)

        line, = ax.plot(t, var.pdf(t), ls=var.linestyle,
                        c='0.5', label=r'$f_X(t)$')
        legend_entries[var.name].append(line)
        line, = ax.plot(X, Y, c=var.color, label='Corrected\nInterior PDF')
        legend_entries[var.name].append(line)
        line, = ax.plot(t, var.pdf(t)/var.cdf(T),
                        c='k', ls=var.linestyle, label=r'$f_X(t)/F_X(T)$')
        legend_entries[var.name].append(line)
    ax.set_xlabel('time')
    ax.set_ylabel(r'$P(X_\mathrm{interior} = t)$')
    ax.set_xlim([0, T])
    ylim = np.max([var.scaled_ylim(T) for var in var_pair])
    ax.set_ylim([0, ylim])
    _legend_from_entries(ax, legend_entries)
    return fig


def corrected_interior_ecdf(obs, var_pair):
    fig, ax = plt.subplots(
        figsize=figure_size['full column'],
        constrained_layout=True
    )
    legend_entries = _init_entries(var_pair)
    T = obs.window_size.max()
    t = np.linspace(0, T, 100)
    for var in var_pair:
        interior, window_sizes = _int_win_from_obs(obs, var.name)
        x, cdf = fw.ecdf_windowed(interior, window_sizes, pad_left_at_x=0)

        line, = ax.plot(t, var.cdf(t), ls=var.linestyle,
                        c='0.5', label=r'$F_X(t)$')
        legend_entries[var.name].append(line)
        line, = ax.plot(x, cdf, c=var.color, label='Corrected\nInterior CDF')
        legend_entries[var.name].append(line)
        line, = ax.plot(t, var.cdf(t)/var.cdf(T), c='k', ls=var.linestyle,
                        label=r'$F_X(t)/F_X(T)$')
        legend_entries[var.name].append(line)
    ax.set_xlabel('time')
    ax.set_ylabel(r'$P(X_\mathrm{interior} \leq t)$')
    ax.set_xlim([0, T])
    ax.set_ylim([0, 1])
    _legend_from_entries(ax, legend_entries)


def discrete_exterior_est(obs, var_pair, times_allowed):
    fig, ax = plt.subplots(
        figsize=figure_size['full column'],
        constrained_layout=True
    )
    legend_entries = _init_entries(var_pair)
    T = obs.window_size.max()
    t = np.linspace(0, T, 101)
    for var in var_pair:
        exterior = _ext_from_obs(obs, var.name)
        ext_bins = times_allowed[1:] - np.diff(times_allowed)/2
        y, t_bins = np.histogram(exterior, bins=ext_bins, density=1)
        X, Y = stats.bars_given_hist(y, t_bins)
        X -= np.diff(times_allowed)/2
        line, = ax.plot(X, Y, c=var.color, label='Exterior Histogram')
        legend_entries[var.name].append(line)

        scale = T - scipy.integrate.quad(var.cdf, 0, T)[0]
        line, = ax.plot(t, var.sf(t)/scale, c='k', ls=var.linestyle,
                        label='Rescaled\nsurvival function')
        legend_entries[var.name].append(line)
    ax.set_xlabel('time')
    ax.set_ylabel(r'$P(X_\mathrm{exterior} = t)$')
    ax.set_xlim([0, T])
    ax.set_ylim([0, ax.get_ylim()[1]])
    _legend_from_entries(ax, legend_entries)


def discrete_bars_from_cdf(obs, var_pair, times_allowed):
    fig, ax = plt.subplots(
        figsize=figure_size['full column'],
        constrained_layout=True
    )
    legend_entries = _init_entries(var_pair)
    T = obs.window_size.max()
    t = np.linspace(0, T, 100)
    for var in var_pair:
        interior, window_sizes = _int_win_from_obs(obs, var.name)
        x, cdf = fw.ecdf_windowed(interior, window_sizes,
                                  times_allowed=times_allowed)

        X, Y = stats.bars_given_discrete_cdf(x, cdf)
        line, = ax.plot(t, var.pdf(t), ls=var.linestyle,
                        c='0.5', label=r'$F_X(t)$')
        legend_entries[var.name].append(line)
        line, = ax.plot(X, Y, c=var.color, label='Corrected\nInterior PDF')
        legend_entries[var.name].append(line)
        line, = ax.plot(t, var.pdf(t)/var.cdf(T), c='k', ls=var.linestyle,
                        label=r'$F_X(t)/F_X(T)$')
        legend_entries[var.name].append(line)
    ax.set_xlabel('time')
    ax.set_ylabel(r'$P(X_\mathrm{interior} = t)$')
    ax.set_xlim([0, T])
    ylim = np.max([var.scaled_ylim(T) for var in var_pair])
    ax.set_ylim([0, ylim])
    _legend_from_entries(ax, legend_entries)


def scaling_normalizers(obs, var_pair):
    fig, ax = plt.subplots(
        figsize=figure_size['full column'],
        constrained_layout=True
    )

    legend_entries = _init_entries(var_pair)
    T = obs.window_size.max()
    cdf_int_to_ext_cdf = {}
    for var in var_pair:
        exterior = _ext_from_obs(obs, var.name)
        x_ext, cdf_ext = fw.ecdf(exterior.wait_time.values, pad_left_at_x=0)
        interior, window_sizes = _int_win_from_obs(obs, var.name)
        x_int, cdf_int = fw.ecdf_windowed(interior, window_sizes,
                                          pad_left_at_x=0)
        # now compute integral of CDF w.r.t.
        ccdf_int = fw.stats._ccdf_int(x_int, cdf_int)

        t_all = np.sort(np.concatenate((x_int, x_ext)))
        resampled_int_cdf = np.interp(t_all, x_int, ccdf_int)
        resampled_ext_cdf = np.interp(t_all, x_ext, cdf_ext)

        def err_f(ab, t, integrated_cdf, exterior_cdf):
            return np.linalg.norm(
                ab[0] * t + ab[1] * integrated_cdf - exterior_cdf
            )

        opt_out = scipy.optimize.minimize(
            err_f,
            x0=[1, -1],
            args=(t_all, resampled_int_cdf, resampled_ext_cdf),
            bounds=((0, np.inf), (-np.inf, 0)),
        )
        if not opt_out.success:
            raise ValueError("Unable to compute F_X(T)!")
        a, b = opt_out.x
        cdf_int_to_ext_cdf[var.name] = {'a': a, 'b': b}

        line, = ax.plot(x_ext, cdf_ext, c=0.7*np.array(var.color), ls='--',
                        lw=1.5,
                        label=r'$\hat{F}_\mathrm{ext}(t)$')
        legend_entries[var.name].append(line)

        line, = ax.plot(x_int, ccdf_int, c=var.color, ls=':',
                        label=r"$\mathcal{I}[\hat{F}_X](t)$")
        legend_entries[var.name].append(line)

        line, = ax.plot(x_int, a*x_int + b*ccdf_int, c=var.color, ls='-', lw=1,
                        label=f"${a:02.1f} - {-b:02.1f}"
                              r"\mathcal{I}[\hat{F}_X](t)$")
        legend_entries[var.name].append(line)
    ax.set_xlabel('time')
    ax.set_ylabel(r'$P(X_\mathrm{exterior} \leq t)$')
    ax.set_xlim([0, T])
    ax.set_ylim([0, 1])
    _legend_from_entries(ax, legend_entries, loc='upper left')

    return fig


def fully_corrected_interior_pdf(obs, var_pair, traj_cols=['replicate'],
                                 int_bins=51):
    fig, ax = plt.subplots(
        figsize=figure_size['full column'],
        constrained_layout=True
    )
    legend_entries = _init_entries(var_pair)

    traj_windows = obs.groupby(traj_cols)['window_size'].first().values
    # now sorted
    traj_windows, window_cdf = stats.ecdf(traj_windows, pad_left_at_x=0)
    window_sf = 1 - window_cdf

    T = obs.window_size.max()
    t = np.linspace(0, T, 101)
    for var in var_pair:
        interior, window_sizes = _int_win_from_obs(obs, var.name)
        exterior = _ext_from_obs(obs, var.name)
        _, _, _, _, F_T = fw.ecdf_ext_int(
            exterior, interior, window_sizes, window_sf=window_sf
        )
        window_i = np.searchsorted(traj_windows, interior) - 1
        correction = 1/window_sf[window_i]/(window_sizes - interior)
        y, t_bins = np.histogram(interior, weights=correction, bins=int_bins,
                                 density=1)
        X, Y = stats.bars_given_hist(y, t_bins)

        line, = ax.plot(t, var.pdf(t), ls=var.linestyle,
                        c='k', label='True $f_X(t)$')
        legend_entries[var.name].append(line)
        line, = ax.plot(X, Y*F_T, c=var.color,
                        label='Fully Corrected\nInterior PDF')
        legend_entries[var.name].append(line)
    ax.set_xlabel('time')
    ax.set_ylabel(r'$P(X_\mathrm{interior} = t)$')
    ax.set_xlim([0, T])
    ylim = np.max([var.scaled_ylim(T) for var in var_pair])
    ax.set_ylim([0, ylim])
    _legend_from_entries(ax, legend_entries)

    return fig


def int_ext_cdf_comparison(obs, var_pair, traj_cols=['replicate'],
                           ext_bins='auto', **kwargs):
    fig, ax = plt.subplots(
        figsize=figure_size['full column'],
        constrained_layout=True
    )
    legend_entries = _init_entries(var_pair)
    window_sf = fw.window_sf(obs, traj_cols)
    T = obs.window_size.max()
    t = np.linspace(0, T, 101)
    for var in var_pair:
        exterior = _ext_from_obs(obs, var.name)
        interior, window_sizes = _int_win_from_obs(obs, var.name)
        all_times, cdf_int, cdf_ext, Z_X, F_T = fw.ecdf_ext_int(
            exterior, interior, window_sizes, window_sf=window_sf,
            **kwargs
        )
        y, t_bins = np.histogram(exterior, bins=ext_bins, density=1)
        X, Y = stats.bars_given_hist(y, t_bins)

        line, = ax.plot(t, var.cdf(t), ls=var.linestyle,
                        c='k', label='True $F_X(t)$')
        legend_entries[var.name].append(line)
        line, = ax.plot(all_times, cdf_int*F_T, c=var.color, alpha=0.8,
                        label='Interior CDF estimate')
        legend_entries[var.name].append(line)
        line, = ax.plot(X, 1 - Y/Z_X, c=var.color, alpha=0.8,
                        label='Exterior CDF estimate')
        legend_entries[var.name].append(StairLine(line))
    ax.set_xlabel('time')
    ax.set_ylabel(r'$P(X_\mathrm{interior} \leq t)$')
    ax.set_xlim([0, T])
    ax.set_ylim([0, 1])
    _legend_from_entries(ax, legend_entries)

    return fig


def ext_est_std(obs, var_pair, traj_cols=['replicate'], ext_bins='auto'):
    """Plot error bars on the final exterior estimator."""
    fig, ax = plt.subplots(
            figsize=figure_size['full column'],
            constrained_layout=True
        )
    legend_entries = _init_entries(var_pair)
    window_sf = fw.window_sf(obs, traj_cols)
    T = obs.window_size.max()
    t = np.linspace(0, T, 101)
    for var in var_pair:
        exterior = _ext_from_obs(obs, var.name)
        interior, window_sizes = _int_win_from_obs(obs, var.name)
        all_times, cdf_int, cdf_ext, Z_X, F_T = fw.ecdf_ext_int(
            exterior, interior, window_sizes, window_sf=window_sf
        )
        y, t_bins = np.histogram(exterior, bins=ext_bins)
        bin_centers = (t_bins[:-1] + t_bins[1:]) / 2
        Z_hist = np.sum(y*np.diff(t_bins))
        hist_vals = (1 - y/Z_hist/Z_X)
        # standard error of mean, Normal approximation.
        hist_sigma = np.sqrt(y)/Z_hist/Z_X

        line, = ax.plot(t, var.cdf(t), ls=var.linestyle,
                        c='k', label='True $F_X(t)$')
        legend_entries[var.name].append(line)
        line, = ax.plot(bin_centers, hist_vals, c=var.color,
                        label='Exterior CDF estimate')
        legend_entries[var.name].append(line)
        fill = ax.fill_between(
            bin_centers,
            hist_vals + 2*hist_sigma,
            hist_vals - 2*hist_sigma,
            facecolor=var.color,
            edgecolor=None,
            alpha=0.5,
            label=r'$\pm2\sigma$'
        )
        legend_entries[var.name].append(fill)
    ax.set_xlabel('time')
    ax.set_ylabel(r'$P(X_\mathrm{interior} \leq t)$')
    ax.set_xlim([0, T])
    ax.set_ylim([-0.1, 1.1])
    _legend_from_entries(ax, legend_entries)

    return fig


def final_cdf_comparison(obs, var_pair, traj_cols=['replicate'],
                         ext_bins='auto'):
    fig, ax = plt.subplots(
            figsize=figure_size['full column'],
            constrained_layout=True
        )

    legend_entries = _init_entries(var_pair)
    window_sf = fw.window_sf(obs, traj_cols)
    T = obs.window_size.max()
    t = np.linspace(0, T, 101)
    for var in var_pair:
        exterior = _ext_from_obs(obs, var.name)
        interior, window_sizes = _int_win_from_obs(obs, var.name)
        bin_centers, final_est = fw.ecdf_combined(
            exterior, interior, window_sizes, window_sf=window_sf
        )

        line, = ax.plot(t, var.cdf(t), ls=var.linestyle,
                        c='k', label='True $F_X(t)$')
        legend_entries[var.name].append(line)
        line, = ax.plot(bin_centers, final_est, c=var.color,
                        label='Combined estimator', alpha=0.8)
        legend_entries[var.name].append(line)
    ax.set_xlabel('time')
    ax.set_ylabel(r'$P(X_\mathrm{interior} \leq t)$')
    ax.set_xlim([0, T])
    ax.set_ylim([-0.1, 1.1])
    _legend_from_entries(ax, legend_entries)


def _multi_t_waits(var_pair, window_size, N_traj):
    # generate same distribution, but with three different window sizes
    het_trajs = [
        fw.ab_window(
            [var.rvs for var in var_pair],
            window_size=window,
            offset=-100*np.sum([var.mean() for var in var_pair]),
            num_replicates=N_traj,
            states=[var.name for var in var_pair]
        )
        for window in np.array([1/2, 1, 2])*window_size
    ]
    het_trajs = pd.concat(het_trajs, ignore_index=True)

    multi_T_waits = fw.sim_to_obs(
        het_trajs, traj_cols=['window_end', 'replicate']
    )
    return multi_T_waits


def _incorrect_multi_window_pdf(obs, var_pair, int_bins=51):
    fig, ax = plt.subplots(
        figsize=figure_size['full column'],
        constrained_layout=True
    )

    legend_entries = _init_entries(var_pair)
    T = obs.window_size.max()
    t = np.linspace(0, T, 101)
    for var in var_pair:
        interior, window_sizes = _int_win_from_obs(obs, var.name)
        correction = 1/(window_sizes - interior)
        y, t_bins = np.histogram(interior, weights=correction, bins=int_bins,
                                 density=1)
        X, Y = stats.bars_given_hist(y, t_bins)

        line, = ax.plot(t, var.pdf(t)/var.cdf(T), ls=var.linestyle,
                        c='k', label='$f_X(t)/F_X(T)$')
        legend_entries[var.name].append(line)
        line, = ax.plot(X, Y, c=var.color, label='Corrected\nInterior PDF')
        legend_entries[var.name].append(line)

    ax.set_xlabel('time')
    ax.set_ylabel(r'$P(X_\mathrm{interior} = t)$')
    ax.set_xlim([0, T])
    ylim = np.max([var.scaled_ylim(T) for var in var_pair])
    ax.set_ylim([0, ylim])
    _legend_from_entries(ax, legend_entries)

    return fig


def multi_window_demo(var_pair, window_size, N_traj=30_000):
    multi_T_waits = _multi_t_waits(var_pair, window_size, N_traj)
    fig1 = _incorrect_multi_window_pdf(multi_T_waits, var_pair)
    fig2 = corrected_interior_pdf(multi_T_waits, var_pair,
                                  traj_cols=['window_end', 'replicate'])
    return fig1, fig2
