"""Utilities for massaging DataFrames"""

import re
import json
import pandas as pd
import numpy as np

def array_from_numpy_string(s):
    """For unserializing np.array's after DataFrame.to_csv

    Notes
    -----

    As of 2019-03-13, looks like instead of commas separating the numpy array
    elements, *two* spaces are printed. We just replace these spaces with
    commas after removing all other spaces and load the string as json."""

    s = re.sub(r' *\[ +', '[', s)
    s = re.sub(' ([-0-9 ])', '\\1,', s)
    s = re.sub('\n', ',', s)
    return np.array(json.loads(s))

def pivot_loci(df, pivot_cols=['x', 'y', 'z'], spot_col='spot'):
    """Move between "long" and "short" forms for the spot id column.

    Simply put, we want to be able to transform between the following two
    dataframes::

        Condensed form
                                                   X1   Y1   Z1   X2   Y2   Z2     t foci
        locus genotype exp.rep meiosis cell frame
        HET5  WT       2       t0      1    1      1.6  2.4  3.1  2.1  1.5  3.1    0  unp
                                            2      1.9  1.5  3.1  1.9  2.5  3.1   30  unp
                                            3      2.0  1.8  3.0  1.5  2.5  3.4   60  unp
                                            4      2.1  1.9  3.0  1.4  2.2  3.4   90  unp
                                            5      2.2  1.8  3.0  1.5  2.4  3.4  120  unp

    and
    ::

        "Long" form
                                                        X    Y    Z      t foci
        locus genotype exp.rep meiosis cell frame spot
        HET5  WT       2       t0      1    1     1     1.6  2.4  3.1    0  unp
                                            2     1     1.9  1.5  3.1   30  unp
                                            3     1     2.0  1.8  3.0   60  unp
                                            ...
                                            1     2     2.1  1.5  3.1    0  unp
                                            2     2     1.9  2.5  3.1   30  unp
                                            3     2     1.5  2.5  3.4   60  unp

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
    cols_to_pivot = [col for col in cols if any([r.match(col) for r in rs])]
    # if we are creating the numbered columns
    if spot_col in df.index.names and len(cols_to_pivot) == 0 \
    and all(col in df.columns for col in pivot_cols):
        extra_cols = list(set(cols) - set(pivot_cols))
        def rename_cols(data):
            data = data.copy()
            spot_id = str(data.index.get_level_values(spot_col)[0])
            for col in df.columns:
                if col in pivot_cols:
                    data[col+spot_id] = data[col].copy()
                del data[col]
            data.index = data.index.droplevel(spot_col)
            return data
        # make one dataframe per spot_id, with just that data in it
        dfs = [rename_cols(data) for _, data in df.groupby(spot_col)]
        # now copy out the non-traj data (that wasn't split apart in prev line)
        df = df[extra_cols].copy()
        # flatten across the "spot_id" dimension by removing that index and
        # then deleting redundant index values (this assumes that extra_cols
        # are always constant as you change spot_id), if they are not you
        # should be including them in the "traj" cols. Otherwise, only the
        # "spot==1" value is kept here
        df.index = df.index.droplevel(spot_col)
        df = df[~df.index.duplicated(keep='first')]
        # now add back in the columns we manually split apart
        for data in dfs:
            df[data.columns] = data
    # if we are creating a spot column from numbered columns
    elif spot_col not in df.index.names and len(cols_to_pivot) > 0:
        extra_cols = list(set(cols) - set(cols_to_pivot))
        # loci id => [existing column names]
        spot_cols = {}
        for col in cols_to_pivot:
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


