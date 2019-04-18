from multi_locus_analysis import analytical

import pytest
import itertools



@pytest.mark.skip(reason='Not implemented')
def test_cvv_produces_msd():
    pass

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
    full_params = map(lambda p: analytical.un_rouse_nondim(tOverDeltas, *p), params)
    with multiprocessing.Pool(31) as p:
        f = lambda p: analytical.rouse_cvv(*p)
        cvvs = p.map(f, full_params)
    return cvvs
