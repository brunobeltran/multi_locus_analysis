import numpy as np
import scipy.spatial

from .moments import *
from .distributions import *


def convex_hull(df, xcol='x', ycol='y', zcol=None, tcol='t', allowed_t=None,
                max_t=None, volume=False, area=False):
    """
    Compute the convex hull of a trajectory.

    For a DataFrame containing a trajectory with (X,[Y,[Z]]) values in `xcol`,
    `ycol`, and `zcol` (respectively), use scipy.spatial.ConvexHull (uses
    "QHull" under the hood) to calculate the convex hull (including
    area/volume), optionally only looking at certain times along the
    trajectory.

    2D by default (xcol='x', ycol='y').
    """
    cols = [x for x in [xcol, ycol, zcol] if x is not None]
    good_ix = ~np.isnan(df[xcol])
    if allowed_t is not None:
        good_ix = good_ix & np.isin(df[tcol], allowed_t)
    if max_t is not None:
        good_ix = good_ix & (df[tcol] <= max_t)
    points = df[cols].values[good_ix, :]
    # ensure we have enough points to even construct a single simplex given the
    # dimension that we're interested in, otherwise return 0 (the correct
    # "measure" of the reduced dimensional simplex may not be zero, but this is
    # ignored, so that all returned "n-volumes" are measured for the same n)
    min_points = 4 if zcol else 3
    if len(df[cols].drop_duplicates()) < min_points:
        return np.nan
    try:
        qhull = scipy.spatial.ConvexHull(points)
    except scipy.spatial.qhull.QhullError:
        # likely co-planarity/co-linearity
        return np.nan
    if volume:
        return qhull.volume
    elif area:
        return qhull.area
    else:
        return qhull
