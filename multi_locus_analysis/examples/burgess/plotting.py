"""Some example plots with the existing code."""
from . import *
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


