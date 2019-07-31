"""Miscellaneous analysis of the Burgess data not worth saving anywhere
else."""

from . import *

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import scipy

def alpha_fit():
    """We notice from bplot.mscds_by_meiotic_stage that some of the t0 times
    actually spend a reasonable time not confined. We attempt to use these to
    measure :math:`alpha`.

    First we extract the relevant mscd's (we're looking for the two "slowest"
    trajetories::

        >>> m = pd.read_csv('./msds_dvel_by_wait.csv')
        >>> m0 = m[(m['delta'] == 210) & (m['meiosis'] == 't0')]
        >>> m0
                locus genotype meiosis  delta      mean       std  count       ste
            7     HET5       WT      t0    210  0.806006  0.553395   2119  0.012022
            357   LYS2       SP      t0    210  0.452109  0.508022   1723  0.012239
            707   LYS2       WT      t0    210  0.510926  0.469129   1251  0.013264
            1029  URA3       SP      t0    210  0.599170  0.427937   1430  0.011316
            1395  URA3       WT      t0    210  0.639609  0.472704   1633  0.011698

    So we conclude that it's the LYS2, t0 trajectories that are so "slow". Fit
    their MSCD curves in the linear regime (0, 800)-ish by eye.

    CONCLUSION: :math:`alpha` is something like 0.2.
    """
    mscd_unbound_only = pd.read_csv(burgess_dir / Path('msds_dvel_by_wait.csv'))
    mscd_unbound_only.set_index(condition_cols, inplace=True)
    mslow = mscd_unbound_only.loc['LYS2', :, 't0'].copy()
    mslow = mslow[mslow['delta'] > 0]
    mslow = mslow[mslow['delta'] < 800]
    sp = mslow.loc['LYS2', 'SP', 't0']
    wt = mslow.loc['LYS2', 'WT', 't0']
    return [scipy.stats.linregress(np.log10(d['delta']), np.log10(d['mean']))
            for d in [sp, wt]]

def get_confinement_values_by_hand():
    """MSDs in this dataset are too incomplete/noisy to be very confident of
    any automatic method to call confinement levels. this function allows you
    to specify them by hand."""
    confinements = {}
    cmap = cmap_from_list(mscd_unbound_only.reset_index()['meiosis'].unique())
    for label, d in mscd_unbound_only.groupby(['locus', 'genotype', 'meiosis']):
        plt.cla()
        plt.errorbar(d.reset_index()['delta'], d['mean'], d['ste'], c=cmap(label[2]), label=str(label))
        plt.yscale('log')
        plt.xscale('log')
        plt.legend()
        ax = plt.gca()
        ax.yaxis.set_major_locator(MultipleLocator(0.1))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%0.01f'))
        ax.yaxis.set_minor_locator(MultipleLocator(0.01))
        ax.yaxis.set_minor_formatter(FormatStrFormatter(''))
        plt.grid(True)
        plt.pause(0.01)
        ctype = input('confinement type: ')
        cval = float(input('approximage conf level'))
        confinements[label] = {'type': str(ctype), 'R': cval}

