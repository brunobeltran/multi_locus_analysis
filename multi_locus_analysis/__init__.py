"""Module for fitting multi-particle stuff."""
from .stats import *
from .plotting import *

import re

__version__ = "0.0.1"


def pivot_loci(df, pivot_cols=['x', 'y', 'z'], spot_col='spot'):
    """Move between "long" and "short" forms for the spot id column.

    Simply put, we want to be able to transform between the following two
    dataframes:

    ```
                                                   X1   Y1   Z1   X2   Y2   Z2     t foci
        locus genotype exp.rep meiosis cell frame
        HET5  WT       2       t0      1    1      1.6  2.4  3.1  2.1  1.5  3.1    0  unp
                                            2      1.9  1.5  3.1  1.9  2.5  3.1   30  unp
                                            3      2.0  1.8  3.0  1.5  2.5  3.4   60  unp
                                            4      2.1  1.9  3.0  1.4  2.2  3.4   90  unp
                                            5      2.2  1.8  3.0  1.5  2.4  3.4  120  unp
    ```

    and

    ```
                                                        X    Y    Z      t foci
        locus genotype exp.rep meiosis cell frame spot
        HET5  WT       2       t0      1    1     1     1.6  2.4  3.1    0  unp
                                            2     1     1.9  1.5  3.1   30  unp
                                            3     1     2.0  1.8  3.0   60  unp
                                            ...
                                            1     2     2.1  1.5  3.1    0  unp
                                            2     2     1.9  2.5  3.1   30  unp
                                            3     2     1.5  2.5  3.4   60  unp
    ```

    This function can infer which direction to pivot. Because of this, I have
    found using this function much more convenient (and a smaller cognitive
    load) than using a multiindex for the column names and using e.g. pd.unstack
    and friends.

    Parameters
    ----------
    pivot_cols : List<str>
        The names of the columns over which to pivot (without their numerical
        suffixes, these will be inferred).
    spot_col : str
        The name of the column that holds (or will hold) the spot id.

    Returns
    -------
    df : pd.DataFrame
        The pivot-ed DataFrame.
    """
    cols = list(df.columns)
    rs = [re.compile(col+'([0-9]+)') for col in pivot_cols]
    pivot_cols = [col for col in cols if any([r.match(col) for r in rs])]
    # if we are creating the numbered columns
    if spot_col in cols and len(pivot_cols) == 0:
        return df
    # if we are creating a spot column from numbered columns
    elif spot_col not in cols and len(pivot_cols) > 0:
        extra_cols = list(set(cols) - set(pivot_cols))
        # loci id => [existing column names]
        spot_cols = {}
        for col in pivot_cols:
            for r in rs:
                if r.match(col):
                    spot_id = int(r.match(col).groups()[0])
                    if spot_id in spot_cols:
                        spot_cols[spot_id].append(col)
                    else:
                        spot_cols[spot_id] = [col]
        spot_dfs = {}
        for spot_id in spot_cols:
            # trasnform loci id => [existing columns names]
            # to loci id => df with only columns from that spot id
            spot_dfs[spot_id] = df[spot_cols[spot_id]].copy()
            spot_dfs[spot_id].columns = [col[:-len(str(spot_id))]
                    for col in spot_dfs[spot_id].columns]
            # copy over non-index, non-pivot columns as-is
            for col in extra_cols:
                spot_dfs[spot_id][col] = df[col]
            # now add the correct spot value to a new column for each new
            # dataframe we've created
            spot_dfs[spot_id][spot_col] = spot_id
        df = pd.concat(list(spot_dfs.values()))
        df = df.set_index(spot_col, append=True)
    else:
        raise ValueError('''Could not determine which way to pivot.
Either your pivot_cols must exist as numbered columns or your spot_col column
should exist, but not neither or both.''')
    return df

