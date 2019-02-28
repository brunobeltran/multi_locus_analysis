"""For computing analytical results relevant to diffusing loci"""
import numpy as np
from scipy.signal import savgol_filter, savgol_coeffs

from functools import lru_cache
from pathlib import Path

def frac_msd_confined_fit(t, D, Cr, t0, deltaT=1, alpha=0.5):
    """puts parameters in more convenient order for fixing known ones with
    kwargs. deltaT := t1 - t0 you probably want to specify before hand
    instead of letting it scan both t0 and t1."""
    return frac_msd_confined(t, alpha, D, Cr, t0, t0+deltaT)

def frac_msd_confined(t, alpha, D, Cr, t0, t1):
    r"""Approximation of MSD of a fractional brownian motion in confinement.

    for short times, MSD in log-log space will be like alpha*t + D

    for long times, MSD will be constant (at confinement radius)

    here, we arbitrarily specify crossover times t0 and t1, in between these
    the MSD is assumed to be a cubic spline that interpolates between the short
    and long time behavior.

    The cubic function s(t) satisfying

    .. math::

        s(t_0) = \alpha*t + D

        s(t_1) = C_r

        s'(t_0) = \alpha

        s'(t_1) = 0

    with the continuity condition :math:`s''(t_1) = 0` is given by the equation

    .. math::
        s(t) = a (t - k - t_1)(t - t_1)^2 + C_r

    where

    .. math::
        k = \frac{\alpha(t_0 - t_1)(2t_0 + t_1) + 3(D - C_r)(t_0 - t_1)}
                 {\alpha (t_0 + t_1) + 2(D - C_r)}

        a = -\frac{\alpha(t_0 - t_1) + 2(D - C_r)}
                  {(t_0 - t_1)^3}

    agh the formula for a above is wrong, for now just use

    .. math::
        a = -\frac{\alpha}{(t_0 - t_1)(3t_0 - 3t_1 - 2k)}


    warning
    -------
    this function is pretty suboptimal, since only enforcing
    :math:`C^{(2)}` at the right endpoitn of the spline means that the curve
    will go slightly above the max Cr before coming back down.

    in principle, we should recalculate these coefficients to make a better
    fit, but I have lab meeting right now.
    """
    t = np.array(t)
    result = np.nan*np.zeros_like(t)
    result[t < t0] = alpha*t[t<t0] + D
    result[t > t1] = Cr
    if np.any(np.isnan(result)):
        knum = alpha*(t0 - t1)*(2*t0 + t1) + 3*(D - Cr)*(t0 - t1)
        kdenom = alpha*(t0 + t1) + 2*(D - Cr)
        k = knum/kdenom
        a = alpha/((t0 - t1)*(3*t0 - 3*t1 - 2*k))
        t = t[np.isnan(result)]
        result[np.isnan(result)] = a*(t - k - t1)*(t - t1)**2 + Cr
    return result

def vc(t, delta, beta):
    """velocity correlation of locus on rouse polymer. beta = alpha/2."""
    return ( np.power(np.abs(t - delta), beta)
           + np.power(np.abs(t + delta), beta)
           - 2*np.power(np.abs(t), beta)
           )/( 2*np.power(delta, beta) )

def vcp(t, delta, beta):
    return ( np.power(np.abs(t + delta), beta)
           - np.power(np.abs(t), beta) )*np.power(delta, beta)

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

def haar_coeffs_w(n):
    return np.sign(np.arange(-n, n+1))/n/(n+1)

@lru_cache(maxsize=128, typed=False)
def savgol_coeffs_w(n, order=3):
    return savgol_coeffs(2*n+1, order, deriv=1, use='dot')

@lru_cache(maxsize=128, typed=False)
def savgol_coeffs_what(n, order=3):
    return -np.cumsum(savgol_coeffs_w(n, order))

def wavelet_c_from_w_hat(what):
    window_size = len(what)
    n = (window_size-1)/2

    raise NotImplementedError('TODO')

def savgol_A(k, n):
# equation 10, lena's BPJ 2014
    A = 0
    for i in range(-n, k+n-2 + 1):
        A += np.power(savgol_chat(i,k,n), 2)
    return A

def savgol_B(k, n):
# equation 10, lena's BPJ 2014
    B = 0
    for i in range(-n, k+n-1 + 1):
        B += np.power(savgol_c(i,k,n), 2)
    return B/2

def savgol_c(i, k, n):
    if i == -n:
        return -savgol_chat(i, k, n)
    elif i <= k + n - 2:
        return savgol_chat(i-1,k,n) - savgol_chat(i,k,n)
    elif i == k + n -1:
        return savgol_chat(i-1,k,n)
    else:
        raise ValueError('-n <= i <= k+n-1 required')

def savgol_chat(i,k,n):
    hi = 1 if (0 <= i < k) else 0
    chat = hi
    for j in range(max(-n, i-k+1), min(n-1,i) + 1):
        chat -= savgol_what(j,n)
    return chat

@lru_cache(maxsize=128, typed=False)
def savgol_what(j, n, order=3):
    return savgol_coeffs_what(n, order)[n + j]
