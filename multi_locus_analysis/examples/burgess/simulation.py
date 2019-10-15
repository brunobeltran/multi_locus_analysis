"""run some BD simulations using our lab's `wlcsim` module to compare to the
data."""
import os
import multiprocessing
import socket # for getting hostname
from pathlib import Path
import datetime

import numpy as np
import pandas as pd

from wlcsim.bd import rouse
from pscan import Scan


def run_interior_sim(FP):
    """for various different "connectivities" (fraction of beads "homolog
    paired") run 100 or so BD simulations each to get statistics on looping
    probabilities that are comparable to the experiment."""

    # FP = 0.1 # fraction loci homolog paired
    N = int(1e2+1); L = 17475; R = 1000; b = 15; D = 2e7 # ChrV
    dt = rouse.recommended_dt(N, L, b, D)
    Aex = 100; # strength of force confining beads within nucleus
    tether_list = np.array([]).astype(int) # no tethered loci

    # now define the time grid on which to run the BD, and specify which of those times to save

    Nt = 1e8; # total number of equi-space time steps, about 2500s
    t = np.arange(0, Nt*dt, dt) # take Nt time steps
    i30 = np.argmin(np.abs(t - 30)) # index at "30s"
    t_save = t[0::i30] # save every 30s

    # now run the simulation
    tether_list, loop_list, X = rouse.rouse_homologs(N, FP, L, b, D, Aex, R, R, R, t, t_save, tether_list)

    X1, X2 = rouse.split_homologs_X(X, N, loop_list)

    # make bead #, t, and binary "looped"/"tethered" arrays to same size as X1/X2
    bead_id = np.arange(N, dtype=int)
    bead_id = bead_id[None,:]
    bead_id = np.tile(bead_id, (len(t_save), 1))

    t = t_save[:,None]
    t = np.tile(t, (1, N))

    is_looped = np.zeros((N,)).astype(int)
    if len(loop_list) > 1:
        is_looped[loop_list[1:,0]] = 1
    is_looped = is_looped[None,:]
    is_looped = np.tile(is_looped, (len(t_save), 1))

    is_tethered = np.zeros((N,)).astype(int)
    if len(tether_list) > 0:
        is_tethered[thether_list] = 1
    is_tethered = is_tethered[None,:]
    is_tethered = np.tile(is_tethered, (len(t_save), 1))

    # finally, put all these together into a DataFrame
    df = np.array([arr.flatten() for arr in [t, bead_id, X1[:,:,0], X1[:,:,1],
            X1[:,:,2], X2[:,:,0], X2[:,:,1], X2[:,:,2], is_looped, is_tethered]])
    df = pd.DataFrame(df.T, columns=['t', 'bead', 'X1', 'Y1', 'Z1', 'X2', 'Y2', 'Z2',
            'is_loop', 'is_tether'])
    df['FP'] = FP
    for col in ['bead', 'is_loop', 'is_tether']:
        df[col] = df[col].astype(int)
    return df

def save_interior_sim(p):
    FP = p['FP']
    tether_list = p['tether_list']
    base_dir = p['output_dir']

    hostname = socket.gethostname()

    # first call dibs on a unique output directory
    sim_num = 0
    while True:
        # as long as the OS is sane, only one thread can successfully mkdir
        try:
            # periods aren't allowed in hostnames, so use as delimiter
            sim_dir = Path(base_dir) / Path('homolog-sim.' + str(hostname) + '.' + str(sim_num))
            os.makedirs(sim_dir, mode=0o755)
            break
        except OSError:
            sim_num = sim_num + 1

    df = run_interior_sim(FP)
    df.to_csv(sim_dir / Path('all_beads.csv'))

def get_bead_df(dir):
    pass

if __name__ == '__main__':
    # save each run of this script in a unique directory
    base_name = './homolog-sim/no-tether'
    run_num = 0
    while True:
        # as long as the OS is sane, only one thread can successfully mkdir
        try:
            # periods aren't allowed in hostnames, so use as delimiter
            base_dir = base_name + '.' + str(run_num)
            os.makedirs(base_dir, mode=0o755)
            break
        except OSError:
            run_num = run_num + 1

    # set up parameters (FP, tether lists, output directories)
    params = {'FP': [0.2],
              'tether_list': [np.array([]).astype(int)],
              'output_dir': [base_dir]}
    scan = Scan(params)
    scan.add_count(lambda p: 100)

    # set up multiprocessing
    num_cores = multiprocessing.cpu_count() - 1
    p = multiprocessing.Pool(num_cores)

    # now run simulations, one per core until complete
    script_name = os.path.basename(__file__)
    print(script_name + ': Running scan!')
    for params in p.imap_unordered(save_interior_sim, scan.params(), chunksize=1):
        print(script_name + ": " + datetime.datetime.now().isoformat()
              + ": completed run with params: " + str(params))
    print(script_name + ': Completed scan!')

    # for params in scan.params():
    #     save_interior_sim(params)
    #     print(script_name + ": " + datetime.datetime.now().isoformat()
    #           + ": completed run with params: " + str(params))
