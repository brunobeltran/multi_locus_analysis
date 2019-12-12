"""Some example plots with the existing code."""
from . import *
import bruno_util.plotting as bplt
from bruno_util.plotting import cmap_from_list

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

from pathlib import Path

def mscds_by_genotype():
    mscd_file = Path('./msds_dvel_by_wait.csv')
    if not mscd_file.exists():
        from .msds import precompute_msds
        precompute_msds()

    mscd_unbound_only = pd.read_csv(mscd_file)

    cmap = cmap_from_list(mscd_unbound_only.reset_index()['genotype'].unique())
    for label, d in mscd_unbound_only.groupby(['locus', 'genotype', 'meiosis']):
        d = d.reset_index()
        d = d[d['delta'] > 0]
        plt.errorbar(d['delta'], d['mean'], d['ste'], c=cmap(label[1]), label=str(label[1]))
    plt.yscale('log'); plt.xscale('log')
    plt.legend()

def mscds_by_meiotic_stage():
    mscd_file = Path('./msds_dvel_by_wait.csv')
    if not mscd_file.exists():
        from .msds import precompute_msds
        precompute_msds()

    mscd_unbound_only = pd.read_csv(mscd_file)


    meiosis_to_time = lambda m: -1 if m == 'ta' else int(m[1:])
    cmap = cmap_from_list(mscd_unbound_only.reset_index()['meiosis'].apply(meiosis_to_time).unique())
    for label, d in mscd_unbound_only.groupby(['locus', 'genotype', 'meiosis']):
        d = d.reset_index()
        d = d[d['delta'] > 0]
        plt.errorbar(d['delta'], d['mean'], d['ste'], c=cmap(meiosis_to_time(label[2])), label=str(label[2]))
    plt.yscale('log'); plt.xscale('log')
    sm = mpl.cm.ScalarMappable(cmap='viridis', norm=mpl.colors.Normalize(vmin=0, vmax=5))
    sm.set_array([])
    plt.colorbar(sm)
    # plt.legend()

def per_cell_msd(msd, cond=None, skip=None, curve_count=10, tri_x0=None,
                 tri_xhalf=None, **kwargs):
    plt.figure()
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('t (s)')
    plt.ylabel('MSD ($\mu{}m^2$)')
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
    bplt.draw_power_law_triangle(alpha=0, x0=tri_x0, width=0.7, orientation='down', x0_logscale=False, label=r'$\alpha=0$')
    bplt.draw_power_law_triangle(alpha=0.5, x0=tri_xhalf, width=0.7, orientation='down', x0_logscale=False, label=r'$\alpha=0.5$')
