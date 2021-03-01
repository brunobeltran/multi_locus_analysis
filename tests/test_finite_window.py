import numpy as np

import multi_locus_analysis.finite_window as fw
import multi_locus_analysis as mla


def test_ecdf_windowed_vs_ecdf():
    times = [1, 2, 1, 1, 2, 2, 3]
    window_sizes = 4e10
    # in the large window size limit, these should match
    tw, cdf_w = fw.ecdf_windowed(times, window_sizes)
    t, cdf = mla.stats.ecdf(times, pad_left_at_x=0)
    assert np.all(tw == t)
    assert np.all(np.isclose(cdf_w, cdf))
    # with specific times_allowed
    times_allowed = np.arange(5)
    tw, cdf_w = fw.ecdf_windowed(times, window_sizes,
                                 times_allowed=times_allowed)
    assert np.all(tw == times_allowed)
    t, cdf = mla.stats.ecdf(times, times_allowed, pad_left_at_x=0)
    assert np.all(np.isclose(cdf_w, cdf))
