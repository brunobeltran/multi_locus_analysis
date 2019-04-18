"""For fitting analytical theory to multi_locus_analysis results"""
from .analytical import vvc_normalized_theory

import numpy as np
import scipy
import scipy.optimize

def fit_full_fix_beta(tdelta, A, tDeltaN, beta=0.5):
    """Fit velocity cross-correlation to the case of two loci on a polymer,
    freely diffusing, but attached to each other (so beta = alpha/2 = 1/2)."""
    t = tdelta[0]
    delta = tdelta[1]
    return vvc_normalized_theory(t/delta, delta, beta, A, tDeltaN)

def get_best_fit_fixed_beta(df, t_col='t', delta_col='delta',
                            cvv_col='cvv_normed',
                            ste_col='ste', p0=[1, 10], bounds=([0, 1.5], [10, 1000]),
                            counts_col=None, hack=False):
    """beta = alpha/2"""
    x = np.stack((df[t_col].values, df[delta_col].values))
    y = df[cvv_col].values
    sigma = df[ste_col].values
    good_ix = np.isfinite(x[0]) & np.isfinite(x[1]) & np.isfinite(y) \
            & (x[1] > 0) & np.isfinite(sigma)
    # in two batches to prevent "divide by zero" warnings
    x = x[:,good_ix]
    y = y[good_ix]
    good_ix = x[0]/x[1] <= 4 # within bounds we can use
    if counts_col:
        good_ix &= df[counts_col].values > 3
    x = x[:,good_ix]
    y = y[good_ix]
    # sigma = sigma[good_ix]
    # precomputed values for VCC function are
    # log10(delta/tDeltaN) \in np.linspace(-3, 3, 25)
    # alpha \in np.linspace(0.25, 1, 31) # steps of 0.025
    # t/delta \in np.linspace(0, 5, 501) # steps of 0.01
    # since the smallest delta we have is "30s" (1 frame), the largest tDeltaN
    # allowed is 30_000.
    # the largest delta in this exp is 1500s (50 frames), but I doubt we have enough
    # accuracy to go below a single frame (i.e. tDeltaN = 30)
    # A bounds empirically determined
    (A, tDeltaN), pcov = scipy.optimize.curve_fit(fit_full_fix_beta, xdata=x,
            ydata=y, p0=p0, bounds=bounds)
    # (A, tDeltaN), pcov = scipy.optimize.curve_fit(fit_full_fix_beta, xdata=x,
    #         ydata=y, p0=p0, sigma=sigma, bounds=([0, 1.5], [10, 1000]))
    if hack:
        print(df['rspot_id'].iloc[0])
        print(df['anum_id'].iloc[0])
        print(df['dec_id'].iloc[0])
        print({'A': A, 'tDeltaN': tDeltaN, 'pcov': pcov})
    return {'A': A, 'tDeltaN': tDeltaN, 'pcov': pcov}

