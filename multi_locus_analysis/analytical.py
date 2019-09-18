"""For computing analytical results relevant to diffusing loci"""
from bruno_util.mittag_leffler import ml as mittag_leffler

import numpy as np
import scipy
from scipy.special import gamma
from scipy.signal import savgol_filter, savgol_coeffs
from numba import jit
import mpmath

from functools import lru_cache
from pathlib import Path

def frac_msd(t, alpha, kbT=1, xi=1):
    """MSD of fractionally diffusing free particle.

    Weber, Phys Rev E, 2010 (Eq 10)"""
    return 3*kbT/xi*np.sin(alpha*np.pi)*t**alpha \
        /(np.pi*(1-alpha/2)*(1-alpha)*alpha)

def frac_cv(t, alpha, kbT=1, xi=1):
    """Velocity autocorrelation of a fractionally-diffusing particle.
    Weber, Phys Rev E, 2010 (Eq 32)"""
    return -(3*kbT/xi)*np.sin(alpha*np.pi)/(np.pi*(2-alpha))*np.abs(np.power(t, alpha-2))

def frac_discrete_cv(t, delta, alpha, kbT=1, xi=1):
    """Discrete velocity autocorrelation of a fractionally-diffusing particle.
    Weber, Phys Rev E, 2010 (Eq 33)"""
    t = np.atleast_1d(t)
    delta = np.atleast_1d(delta)
    # too many divide by t's to rewrite to avoid these warnings in less than
    # 5min, not worth it.
    with np.errstate(divide='ignore', invalid='ignore'):
        eta = delta/t
        t = t + np.zeros_like(eta) # fix to full size if not already
        delta = delta + np.zeros_like(eta) # fix to full size if not already
        cv_delta_t = frac_cv(t, alpha, kbT, xi)/(eta*eta*alpha*(1 - alpha))
        cv_delta_t[eta<=1] = cv_delta_t[eta<=1] \
            *(2 - np.power(1 - eta[eta<=1], alpha) - np.power(1 + eta[eta<=1], alpha))
        cv_delta_t[eta>1] = cv_delta_t[eta>1] \
            *(2 + np.power(eta[eta>1] - 1, alpha) - np.power(1 + eta[eta>1], alpha)) \
            + frac_msd(delta[eta>1] - t[eta>1], alpha, kbT, xi)/delta[eta>1]/delta[eta>1]
        cv_delta_t[t == 0] = frac_msd(delta[t == 0], alpha, kbT, xi)/np.power(delta[t == 0], 2)
    return cv_delta_t

def frac_discrete_cv_normalized(t, delta, alpha):
    """Normalized discrete velocity autocorrelation of a fractionally-diffusing
    particle. Should be equivalent to

        frac_discrete_cv(t, delta, 1, 1)/frac_discrete_cv(0, delta, 1, 1)

    Lampo, BPJ, 2016 (Eq 5)"""
    return (np.power(np.abs(t - delta), alpha)
        - 2*np.power(np.abs(t), alpha)
        + np.power(np.abs(t + delta), alpha)
        )/(2*np.power(delta, alpha))

@jit
def rouse_mode(p, n, N=1):
    """Eigenbasis for Rouse model.

    Indexed by p, depends only on position n/N along the polymer of length N.
    N=1 by default.

    Weber, Phys Rev E, 2010 (Eq 14)"""
    p = np.atleast_1d(p)
    phi = np.sqrt(2)*np.cos(p*np.pi*n/N)
    phi[p == 0] = 1
    return phi

def rouse_mode_coef(p, b, N, kbT=1):
    """k_p: Weber Phys Rev E 2010, after Eq. 18."""
    # alternate: k*pi**2/N * p**2, i.e. k = 3kbT/b**2
    return 3*np.pi**2*kbT/(N*b**2)*p**2

def rouse_mode_corr(p, t, alpha, b, N, kbT=1, xi=1):
    """Weber Phys Rev E 2010, Eq. 21."""
    kp = rouse_mode_coef(p, b, N, kbT)
    return (3*kbT/kp)*mittag_leffler((-kp*t**alpha)/(N*xi*gamma(3-alpha)), alpha, 1)

def simple_rouse_mid_msd(t, b, N, kbT=1, xi=1, num_modes=1000):
    """
    modified from Weber Phys Rev E 2010, Eq. 24.
    """
    rouse_corr = 0
    for p in range(1, num_modes+1):
        k2p = rouse_mode_coef(2*p, b, N, kbT)
        rouse_corr += 12*kbT/k2p*(1 - np.exp(-k2p*t/(N*xi)))
    return rouse_corr + 6*kbT/xi/N*t

def rouse_mid_msd(t, alpha, b, N, kbT=1, xi=1, num_modes=1000):
    """
    Weber Phys Rev E 2010, Eq. 24.
    """
    if alpha == 1:
        return simple_rouse_mid_msd(t, b, N, kbT, xi, num_modes)
    rouse_corr = 0
    for p in range(1, num_modes+1):
        k2p = rouse_mode_coef(2*p, b, N, kbT)
        rouse_corr += 12*kbT/k2p*(1 - mittag_leffler(-k2p*t**alpha/(N*xi*gamma(3-alpha)),
                alpha, 1))
    return rouse_corr + frac_msd(t, alpha, kbT, xi)/N

def tR(alpha, b, N, kbT=1, xi=1):
    """Lampo et al, BPJ, 2016, eq 8"""
    return np.power(N*N*b*b*xi/kbT, 1/alpha)

def tDeltaN(n1, n2, alpha, b, kbT, xi):
    """Lampo et al, BPJ, 2016, eq 11"""
    delN = np.abs(n2 - n1)
    return np.power(delN*delN*b*b*xi/kbT, 1/alpha)

def rouse_cv_mid(t, alpha, b, N, kbT=1, xi=1, min_modes=1000):
    """Velocity autocorrelation of midpoint of a rouse polymer.

    Weber Phys Rev E 2010, Eq. 33."""
    t = np.atleast_1d(t)
    psum = np.zeros_like(t)
    for p in range(1, min_modes+1):
        k2p = rouse_mode_coef(2*p, b, N, kbT)
        coef = -k2p/(N*xi*gamma(3-alpha))
        psum += mittag_leffler(coef*np.power(t, alpha), alpha, alpha-1)
    gam = 1
    return frac_cv(t, alpha, kbT, xi)*(1 + 2*gam*(alpha - 1)*psum)

def rouse_cvv_ep(t, delta, p, alpha, b, N, kbT=1, xi=1):
    """Term in parenthesis in Lampo, BPJ, 2016 Eq. 10."""
    t = np.atleast_1d(t)
    delta = np.atleast_1d(delta)
    p = np.atleast_1d(p)
    tpdelta = np.power(np.abs(t + delta), alpha)
    tmdelta = np.power(np.abs(t - delta), alpha)
    ta = np.power(np.abs(t), alpha)
    kp = rouse_mode_coef(p, b, N, kbT)
    z = -kp/(N*xi*gamma(3 - alpha))
    return -mittag_leffler(z*tpdelta, alpha, 1) \
        + 2*mittag_leffler(z*ta, alpha, 1) \
        - mittag_leffler(z*tmdelta, alpha, 1)

def rouse_large_cvv_g(t, delta, deltaN, b, kbT=1, xi=1):
    """Cvv^delta(t) for infinite polymer.

    Lampo, BPJ, 2016 Eq. 16."""
    k = 3*kbT/b**2
    ndmap = lambda G, arr: np.array(list(map(G, arr)))
    G = lambda x: float(mpmath.meijerg([[],[3/2]], [[0,1/2],[]], x))
    gtmd = ndmap(G, np.power(deltaN, 2)*xi/(4*k*np.abs(t-delta)))
    gtpd = ndmap(G, np.power(deltaN, 2)*xi/(4*k*np.abs(t+delta)))
    gt = ndmap(G, np.power(deltaN, 2)*xi/(4*k*np.abs(t)))
    return 3*kbT/(np.power(delta, 2)*np.sqrt(xi*k))* (
        np.power(np.abs(t - delta), 1/2)*gtmd
      + np.power(np.abs(t + delta), 1/2)*gtpd
      - 2*np.power(np.abs(t), 1/2)*gt
    )


