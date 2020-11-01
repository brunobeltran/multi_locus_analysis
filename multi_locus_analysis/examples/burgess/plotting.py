"""Some example plots with the existing code."""
from . import location_ura_bp, location_cen5_bp, chrv_size_bp
from .styles import col_width, golden_ratio

import bruno_util.plotting as bplt
from bruno_util.plotting import cmap_from_list
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import numpy as np
import pandas as pd

from pathlib import Path


def mscds_by_genotype():
    mscd_file = Path('./msds_dvel_by_wait.csv')
    if not mscd_file.exists():
        from .msds import precompute_msds
        precompute_msds()

    mscd_ubound_only = pd.read_csv(mscd_file)

    cmap = cmap_from_list(mscd_ubound_only.reset_index()['genotype'].unique())
    for label, d in mscd_ubound_only.groupby(['locus', 'genotype', 'meiosis']):
        d = d.reset_index()
        d = d[d['delta'] > 0]
        plt.errorbar(d['delta'], d['mean'], d['ste'], c=cmap(label[1]),
                     label=str(label[1]))
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()


def mscds_by_meiotic_stage():
    mscd_file = Path('./msds_dvel_by_wait.csv')
    if not mscd_file.exists():
        from .msds import precompute_msds
        precompute_msds()

    mscd_ubound_only = pd.read_csv(mscd_file)

    def meiosis_to_time(m):
        return -1 if m == 'ta' else int(m[1:])

    cmap = cmap_from_list(mscd_ubound_only
                          .reset_index()['meiosis']
                          .apply(meiosis_to_time)
                          .unique())
    for label, d in mscd_ubound_only.groupby(['locus', 'genotype', 'meiosis']):
        d = d.reset_index()
        d = d[d['delta'] > 0]
        plt.errorbar(d['delta'], d['mean'], d['ste'],
                     c=cmap(meiosis_to_time(label[2])),
                     label=str(label[2]))
    plt.xscale('log')
    plt.yscale('log')
    sm = mpl.cm.ScalarMappable(cmap='viridis',
                               norm=mpl.colors.Normalize(vmin=0, vmax=5))
    sm.set_array([])
    plt.colorbar(sm)
    # plt.legend()


def per_cell_msd(msd, cond=None, skip=None, curve_count=10, tri_x0=None,
                 tri_xhalf=None, **kwargs):
    plt.figure()
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('t (s)')
    plt.ylabel(r'MSD ($\mu{}m^2$)')
    curve_cols = ['exp.rep', 'cell']
    if skip is None:
        total_curves = len(msd.groupby(curve_cols).first())
        skip = int(total_curves/curve_count)
    if cond is not None:
        plt.title('Per-Cell MSDs for Condition: ' + str(cond))
    for i, (cell_id, data) in enumerate(msd.groupby(curve_cols)):
        if not i % skip == 0:
            continue
        data = data.reset_index()
        data = data[(data['delta'] > 0) & (data['mean'] > 0)]
        data = data.sort_values('delta')
        plt.errorbar(data['delta'], data['mean'], data['ste'], **kwargs)
    #     plt.plot(data['delta'], data['mean'])
    if tri_x0 is None:
        tri_x0 = [30, 1.05]
    if tri_xhalf is None:
        tri_xhalf = [200, 0.67]
    bplt.draw_power_law_triangle(alpha=0, x0=tri_x0, width=0.7,
                                 orientation='down', x0_logscale=False,
                                 label=r'$\alpha=0$')
    bplt.draw_power_law_triangle(alpha=0.5, x0=tri_xhalf, width=0.7,
                                 orientation='down', x0_logscale=False,
                                 label=r'$\alpha=0.5$')


def draw_cells(cells, min_y=0.05, max_y=0.95, label_loc=location_ura_bp,
               cen_loc=location_cen5_bp, chr_size=chrv_size_bp,
               label_colors=None, ax=None):
    r"""
    Render the model of homologous chromosomes with linkages.

    Parameters
    ----------
    cells : array of float array
        Each element should be a list of the link positions between the
        homologous chromosomes for a given cell. The array corresponding to
        each cell *MUST* be sorted.
    min_y : float
        Bottom of chromosome line in figure coordinates.
    max_y : float
        Top of chromosome line in figure coordinates.
    label_loc : float
        Location of label.
    cen_loc : float
        Centromere location.
    chr_size : float
        Total length of chromosome.
    label_colors : Sequence[Color-like], optional
        If provided, the cells will have colored labels, Cell 1, Cell 2, ...
    """
    def chr_coords(s):
        """Map from [0, 1] to locaion on plot."""
        return max_y - (max_y - min_y)*s
    # rescale linkages to [0, 1]
    cells = [np.array(links) / chr_size for links in cells]
    n_cells = len(cells)
    # and all relevant locations
    locus_frac = label_loc / chr_size
    centromere_frac = cen_loc / chr_size
    if ax is None:
        # fill entire figure with invisible axes to draw in
        fig = plt.figure(figsize=(col_width, col_width/golden_ratio))
        ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    # center each of N "cells" directly between N+1 fenceposts spanning [0, 1]
    n_fences = n_cells + 1
    fence_posts = np.linspace(0, 1, n_fences)
    width_per_cell = np.diff(fence_posts)[0]
    cell_centers = (fence_posts[1:] + fence_posts[:-1]) / 2
    # (1/2) times the spacing between centers of two chromosomes in each "cell"
    width_to_chr_center = width_per_cell / 5
    chr_width = 15
    # only works with mixed backends, where 72"PX"/in is always true, otherwise
    # you need to do something like:
    # transAxes.inverted().transform(dpi_scale_trans.transform([1/72, 1/72])
    pt_to_ax = ax.transAxes.inverted().transform(
        ax.get_figure().dpi_scale_trans.transform([1/72, 1/72])
    )
    for i, x in enumerate(cell_centers):
        for dx in [width_to_chr_center, -width_to_chr_center]:
            cap_radius_ax = chr_width/2 * pt_to_ax[1]
            # draw the chromosomes
            ax.plot(
                [[x + dx, x + dx], [x + dx, x + dx]],
                [[chr_coords(0), chr_coords(centromere_frac) - cap_radius_ax],
                [chr_coords(centromere_frac) + cap_radius_ax, chr_coords(1)]],
                transform=ax.transAxes, linewidth=chr_width,
                solid_capstyle='round', color=[50/255, 50/255, 50/255]
            )
            ax.plot(
                [[x + dx, x + dx], [x + dx, x + dx]],
                [[chr_coords(0), chr_coords(centromere_frac) - cap_radius_ax],
                [chr_coords(centromere_frac) + cap_radius_ax, chr_coords(1)]],
                transform=ax.transAxes, linewidth=chr_width-2,
                solid_capstyle='round', color=[197/255, 151/255, 143/255]
            )
            # draw the centromere black dot
            ax.scatter([x + dx], [chr_coords(centromere_frac)],
                    zorder=10, transform=ax.transAxes, s=200, color='k')
            # draw the label, green star
            ax.scatter([x + dx], [chr_coords(locus_frac)],
                    zorder=15, transform=ax.transAxes, s=500, color='g',
                    marker='*', edgecolors='k')
        for linkage in cells[i]:
            ax.plot([x - width_to_chr_center, x + width_to_chr_center],
                    2*[chr_coords(linkage)],
                    color=(0, 0, 1), transform=ax.transAxes,
                    linewidth=5, solid_capstyle='round')
        num_linkages = len(cells[i])
        j = np.searchsorted(cells[i], locus_frac)
        closest_links = []
        if j != 0:
            closest_links.append(cells[i][j - 1])
        if j != num_linkages:
            closest_links.append(cells[i][j])
        closest_links = np.array(closest_links)
        if len(closest_links) > 0:
            linewidths = 1.2*np.ones_like(closest_links)
            closestest_link = np.argmin(np.abs(closest_links - locus_frac))
            linewidths[closestest_link] = 3.5
        for k, linkage in enumerate(closest_links):
            ax.plot([x - width_to_chr_center, x - width_to_chr_center,
                    x + width_to_chr_center, x + width_to_chr_center],
                    [chr_coords(locus_frac), chr_coords(linkage),
                    chr_coords(linkage), chr_coords(locus_frac)],
                    color=(1, 1, 1), transform=ax.transAxes,
                    linewidth=linewidths[k], linestyle='--',
                    dash_capstyle='butt', zorder=100)
        if label_colors:
            # add extra height above chromosomes to account for rounded end
            # caps and then some extra
            ax.transAxes
            ax.text(x, max_y, f'Cell {i}\n', ha='center',
                    va='bottom', color=label_colors[i],
                    transform=ax.transAxes,
                    fontsize=mpl.rcParams['axes.titlesize'])
    return ax

