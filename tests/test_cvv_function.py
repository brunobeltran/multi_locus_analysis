from multi_locus_analysis import analytical

import pytest
import itertools
import pandas as pd
import numpy as np

import multiprocessing



@pytest.mark.skip(reason='Not implemented')
def test_cvv_produces_msd():
    pass

def f(p):
    return analytical.rouse_cvv(*p)

@pytest.mark.skip(reason='Test should only run by hand...lengthy.')
def test_cvv_versus_tom():
    # deltaOverTDeltaNs = np.logspace(-3, 3, 25)
    # deltaOverTDeltaNs = np.power(10, analytical.deltas_)
    # alphas = np.linspace(0.25, 1, 31) # steps of 0.025
    # alphas = analytical.alphas_
    # tOverDeltas_ = np.linspace(0, 5, 501) # steps of 0.01
    # tOverDeltas = analytical.tOverDeltas_
    alphas = [0.25, 0.5, 0.75, 1]
    deltaOverTDeltaNs = np.logspace(-3, 3, 5)
    tOverDeltas = np.linspace(0, 5, 51)
    params = np.stack(list(itertools.product(deltaOverTDeltaNs, alphas)))
    full_params = list(map(lambda p: analytical.un_rouse_nondim(tOverDeltas, *p), params))
    with multiprocessing.Pool(31) as p:
        cvvs = p.map(f, full_params)
    df = pd.DataFrame(np.stack(cvvs))
    df['deltaOverTDeltaN'] = params[:,0]
    df['alpha'] = params[:,1]
    df = df.set_index(['alpha','deltaOverTDeltaN'])
    df.columns = tOverDeltas
    return df
