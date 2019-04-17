"""For computing analytical results relevant to diffusing loci"""
import numpy as np
import scipy
from scipy.signal import savgol_filter, savgol_coeffs

from functools import lru_cache
from pathlib import Path

def vc(t, delta, beta):
    """velocity correlation of locus on rouse polymer. beta = alpha/2."""
    return ( np.power(np.abs(t - delta), beta)
           + np.power(np.abs(t + delta), beta)
           - 2*np.power(np.abs(t), beta)
           )/( 2*np.power(delta, beta) )

def vcp(t, delta, beta):
    return ( np.power(np.abs(t + delta), beta)
           - np.power(np.abs(t), beta) )*np.power(delta, beta)

# precomupted velocity cross-correlation for rouse polymer
# from Lampo et al, was pre-computed on a grid...
deltas_ = np.linspace(-3, 3, 25) # placeholders for corresponding values in logspace
alphas_ = np.linspace(0.25, 1, 31) # steps of 0.025
tOverDeltas_ = np.linspace(0, 5, 501) # steps of 0.01
vvcf_table_ = np.reshape(np.loadtxt(Path(__file__).parent / Path('vvcf_table.csv'),
        delimiter=','), (31, 25, 501))
def calc_vel_corr_fixed_(tOverDelta, deltaOverTDeltaN, alpha):
    # this performs interpolation in logspace for "delta"/deltaOverTDeltaN
    deltaOverTDeltaN = np.log10(deltaOverTDeltaN)
    return scipy.interpolate.interpn((alphas_, deltas_, tOverDeltas_),
                                     vvcf_table_,
                                     (alpha, deltaOverTDeltaN, tOverDelta))
calc_vel_corr_fixed_.vvcf = None
calc_vel_corr_fixed = np.vectorize(calc_vel_corr_fixed_)

def vvc_rescaled_theory(t, delta, beta, A, tDeltaN):
    """velocity cross correlation of two points on rouse polymer."""
    return 2*A*np.power(delta, beta)*(vc(t*delta, delta, beta) - calc_vel_corr_fixed(t, delta/tDeltaN, 2*beta))
