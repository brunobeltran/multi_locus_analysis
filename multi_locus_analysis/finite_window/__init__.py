r"""For correcting interior- and exterior-censored data

This module is designed to facilitate working with Markov processes which have
been observed over finite windows of time.

Suppose you wish to observe a process that switches between states A and B, and
observe the process over a window of time :math:`[0, T]`.

There will be two types of "observed" times, which we we call "interior times"
(those when the entire time spent in the state is observed) and "exterior
times" (the lengths spent in the state that we observed at the start/end point
of our window).

For example, suppose you are observing two fluorescent loci that are either in
contact (A) or not (B). If you measure for 5s, and the loci begin in contact,
come apart at time t=2s, and rejoin at time t=4s,

::

        A A A A A A B B B B B B A A A A
        | - - | - - | - - | - - | - - | -- | -- | ...
        |                             |
    t = 0s    1s    2s    3s    4s    5s


we would say that T = 5s, and you measured one "interior" time of length 2s,
and two exterior times, one of length 2s and one of length 1s (at the start and
end of the window, respectively).

When histogramming these times, we must apply a statistical correction to the
final histogram weights due to the effects of the finite observation window.

The distribution of the exterior times :math:`f_E(t)` is exactly equal to the
survival function of the actual distribution :math:`S(t) = 1 - \int_0^t f(s)
ds` normalized to equal one over the observation interval.
No functions are currently included that leverage this, since extracting
information from the survival function is likely only worth it if a large
fraction of your observations are exterior times.

On the other hand, the interior times distribution is given by :math:`f_I(t) =
f(t)(T - t)` for :math:`t\in[0,T]`. In order to plot the actual shape of
:math:`f(t)` with this biasing removed, we provide the cdf_exact_given_windows
function below.

A typical workflow is, given an array of interior times, :code:`t`, and an
array of the window sizes each time was observed within, :code:`w`, is to
extract the CDF exactly, then optionally convert that to a PDF to display a
regular histogram

.. code-block:: python

    >>> x, cdf = cdf_exact_given_windows(t, w, pad_left_at_x=0)
    >>> xp, pdf = bars_given_cdf(x, cdf)
    >>> confint = simultaneous_confint_from_cdf(0.05, len(t), x, cdf)
    >>> xc, confs = bars_given_confint(x, confint)
    >>> plt.plot(xp, pdf, 'k', xc, confs, 'r-.')

The remainder of the functions herein are related to extracting confidence
intervals for these distributional estimates.

For details of the derivation, see the Biophysical Journal paper in preparation
by Beltran and Spakowitz.
"""
from .simulation import *
from .munging import *
from .stats import *
