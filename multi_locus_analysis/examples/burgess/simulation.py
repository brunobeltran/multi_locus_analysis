"""run some BD simulations using our lab's `wlcsim` module to compare to the
data."""
import datetime
import multiprocessing
import os
from pathlib import Path
import pickle
import socket # for getting hostname
import time
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from bruno_util import numpy as bnp
from pscan import Scan
from wlcsim.bd import rouse
from wlcsim.bd import homolog

from . import sim_D, chrv_size_nuc_chain_um, kuhn_length_nuc_chain, nuc_radius_um
from ... import finite_window as fw


def get_sim_t(Nt, Nlin, overlap, N, L, b, D, preequilibrate=True):
    dt = rouse.recommended_dt(N, L, b, D)
    # terminal relaxation time, i.e. how long to diffuse before equilibrated
    # the factor of two compared to a linear chain accounts for the two
    # polymers being linked
    if preequilibrate:
        t_R = rouse.terminal_relaxation(N, 2*L, b, D)
    else:
        t_R = 0
    t = np.arange(0, t_R + Nt*dt, dt)
    min_save_i = np.where(t >= t_R)[0][0]
    # for saved values of time_params, this contains ~30,60,90,etc
    t_i, i_i = bnp.loglinsample(Nt, Nlin, overlap)
    t_i += min_save_i
    t_save = t[t_i]
    return t, t_save, t_i, i_i


def run_interior_sim(FP, time_params, rouse_params):
    """for various different "connectivities" (fraction of beads "homolog
    paired") run 100 or so BD simulations each to get statistics on looping
    probabilities that are comparable to the experiment."""
    rouse_params = rouse_params.copy()
    rouse_params['FP'] = FP

    tether_list = np.array([]).astype(int)  # no tethered loci for now

    # now define the time grid on which to run the BD, and specify which of those times to save
    dt_params = {k: rouse_params[k] for k in ['N', 'L', 'b', 'D']}
    t, t_save, _, _ = get_sim_t(**time_params, **dt_params)

    # now run the simulation
    tether_list, loop_list, X = homolog.rouse(
        t=t, t_save=t_save, tether_list=tether_list, **rouse_params
    )

    N = rouse_params.pop('N')
    X1, X2 = homolog.split_homologs_X(X, N, loop_list)

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
    base_dir = p['output_dir']
    time_params = p['time_params']
    rouse_params = p['rouse_params']

    # read in existing parameters from the current base_dir if they exist, in
    # order to check if the parameters that we're using match the directory
    # that we're in...to prevent a bunch of simulations with hella different
    # parameters.
    params_dir = Path(base_dir) / Path('params')
    param_file = params_dir / Path('shared.csv')
    # FP is saved in the output dataframe of run_interior_sim, not here
    all_params = pd.Series(dict(time_params, **rouse_params))
    try:
        os.makedirs(params_dir, mode=0o755)
        # if the above doesn't raise, nobody's established the parameters yet
        all_params.to_csv(param_file, header=False)
        shared_params = all_params.copy()
    except OSError:
        if not param_file.exists():
            # somebody might have started at the same time as us, given them a
            # second to write out the params file
            time.sleep(1)
            if not param_file.exists():
                raise OSError(f"Param file not found: {param_file}.")
        shared_params = pd.read_csv(param_file, index_col=0, header=None,
                                    squeeze=True)

    hostname = socket.gethostname()
    # now call dibs on a unique output directory
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

    # make sure to warn the user about which specific subdirectory will have
    # weird parameters if we are in fact doing something weird like that
    if not np.all(shared_params == all_params):
        warnings.warn(f'Params of simulation in new directory ({sim_dir}) do '
                      'not match previous simulations in base directory '
                      f'({param_file}), be careful!')
    # run and save simulation
    print(f"Starting run: {sim_dir} at time: "
          + datetime.datetime.now().isoformat())
    start = time.process_time()
    df = run_interior_sim(FP, time_params, rouse_params)
    elapsed_time = time.process_time() - start
    print(f"Completed run: {sim_dir} at time: "
          + datetime.datetime.now().isoformat()
          + f"\nTime elapsed: {elapsed_time}")

    df.to_csv(sim_dir / Path('all_beads.csv'), index=False)
    all_params.to_csv(sim_dir / Path('params.csv'), header=False)


_ura3_bead_out_of_101 = 20
def get_bead_df(base_dir, bead_id=_ura3_bead_out_of_101, force_redo=False,
                index_col=None):
    """Extracts all ura3 data from simulations in a set of directories."""
    base_dir = Path(base_dir)
    dfs = []
    for sim in base_dir.glob('homolog-sim.*'):
        sim_file = sim / Path('all_beads.csv')
        bead_file = sim / Path(f'bead_{bead_id}.csv')
        if bead_file.exists() and not force_redo:
            bead = pd.read_csv(bead_file)
            dfs.append(bead)
            continue
        # skip failed/in progress simulations
        if not sim_file.exists():
            continue
        df = pd.read_csv(sim_file, index_col=index_col)
        df['sim_name'] = sim
        bead = df[df['bead'] == bead_id].copy()
        t0 = df[df['t'] == df['t'].iloc[0]]
        loops = t0.loc[t0.is_loop == 1, 'bead'].values
        left_n = np.where(loops <= bead_id)[0]
        bead['left_neighbor'] = loops[left_n[-1]] if len(left_n) > 0 else np.nan
        right_n = np.where(loops >= bead_id)[0]
        bead['right_neighbor'] = loops[right_n[0]] if len(right_n) > 0 else np.nan
        bead['min_bead'] = t0['bead'].min()
        bead['max_bead'] = t0['bead'].max()
        dfs.append(bead)
        dfs[-1].to_csv(bead_file, index=False)

    df = pd.concat(dfs, ignore_index=True)
    df = df.set_index(['FP', 'sim_name', 'bead', 't'])
    df = df.sort_index()
    df.to_csv(base_dir / Path(f'all_bead_{bead_id}.csv'))
    return df

def select_exp_times(df, desired_t=np.arange(0, 1501, 30), force_redo=False):
    t_uniq = df['t'].unique()
    min_t = np.min(t_uniq)
    # get indices into t_uniq for time closest to each desired time
    ind = []
    for t in desired_t:
        sim_time = t_uniq - min_t
        ind.append(np.argmin(np.abs(sim_time - t)))
    # the t's to keep
    t = t_uniq[ind]
    if np.any(np.abs(1 - t / np.round(t)) > 0.01):
        raise ValueError("Cannot find simulation times sufficiently similar "
                         "to the experiments.")
    df_exp = df[np.isin(df['t'], t)].copy()
    # if possible, in int index column for time is nice
    if len(np.unique(np.diff(desired_t))) == 1:
        dt = desired_t[1] - desired_t[0]
        df_exp['ti'] = np.round((df_exp['t'] - min_t)/ dt).astype(int)
        # and if we're close enough, just round the time itself as well
        if np.abs(1 - np.round(dt) / dt) < 0.001:
            df_exp['t'] = np.round(dt*df_exp['ti']).astype(int)
    return df_exp

def add_paired_cols(df, paired_distances=None):
    if paired_distances is None:
        # paired_distances = [10, 50, 100, 250, 500, 750, 1000]
        paired_distances = [.250]
    for x in ['X', 'Y', 'Z']:
        df['d' + x] = df[x + '2'] - df[x + '1']
    df['||dX||'] = np.sqrt(
            np.power(df['dX'], 2)
            + np.power(df['dY'], 2)
            + np.power(df['dZ'], 2)
    )
    for dist in paired_distances:
        df['pair' + str(dist)] = df['||dX||'] < dist
    return df

def get_interior_times(df, state_col='pair.25'):  #, TOL=None):
    waitdf = df.groupby(['FP', 'sim_name']).apply(
            fw.discrete_trajectory_to_wait_times, t_col='t', state_col=state_col)
    waitdf = waitdf[waitdf['wait_type'] == 'interior'].copy()
    # # also because we're using floating times, we need to de-duplicate
    # # choose TOL  to be two more decimal points than dt, or so
    # if TOL is None:
    #     TOL = dt/1e2
    # wtimes = np.sort(waitdf['wait_time'].unique().copy())
    # diff = np.append(True, np.diff(wtimes))
    # wtimes = wtimes[diff > TOL]
    # for uniq_time in wtimes:
    #     waitdf.loc[np.isclose(waitdf['wait_time'], uniq_time), 'wait_time'] = uniq_time
    return waitdf

def plot_interior_times(waitdf):
    fig_unpair = plt.figure()
    for FP, data in waitdf.groupby('FP'):
        paired = data[~data['wait_state']]
        try:
            x, cdf = fw.ecdf_windowed(paired['wait_time'].values, paired['window_size'].values, pad_left_at_x=0)
        except:
            continue
        xp, pdf = fw.bars_given_cdf(x, cdf)
        plt.plot(xp, pdf, label='FP = ' + str(FP))
    plt.yscale('log'); plt.xscale('log')
    plt.xlabel('time (s)')
    plt.ylabel('Probaility')
    plt.legend()
    plt.title('Unpaired PDFs')

    fig_pair = plt.figure()
    for FP, data in waitdf.groupby('FP'):
        paired = data[data['wait_state']]
        try:
            x, cdf = fw.ecdf_windowed(paired['wait_time'].values, paired['window_size'].values, pad_left_at_x=0)
        except:
            continue
        xp, pdf = fw.bars_given_cdf(x, cdf)
        plt.plot(xp, pdf, label='FP = ' + str(FP))
    plt.yscale('log'); plt.xscale('log')
    plt.xlabel('time (s)')
    plt.ylabel('Probaility')
    plt.legend()
    plt.title('Paired PDFs')

    return fig_pair, fig_unpair


def run_homolog_param_scan(fp_list=np.linspace(0, 0.1, 11), replicates=25,
                           base_dir=None, material_params=None,
                           time_params=None):
    if base_dir is None:
        # let's us re-run script and just get more replicates
        base_dir = './homolog-sim/'
    base_dir = Path(base_dir)

    if material_params is None:
        material_params = {
             # ChrV
            'N': int(1e2+1), 'L': chrv_size_nuc_chain_um, 'R': nuc_radius_um,
            'b': kuhn_length_nuc_chain, 'D': sim_D,
            'Aex': 100  # strength of force confining beads within nucleus
        }
    if time_params is None:
        # 10 repeats each of FP=[0, 0.02, 0.05, 0.1] on one 32 core node takes about
        # ~1min to run for Nt=1e5,
        # ~14min for Nt=1e6
        # ~3.2hrs for Nt=1e7
        # ~45hrs for Nt=1e8
        # the ChrV polymer has t_R/recommended_dt ~ 1e8, so ~48hrs total
        time_params = {
            'Nt': 1e7,  # total number of equi-space time steps, about 2500s
            'Nlin': 1e3,  # how many time steps per log-spaced group
            'overlap': 0.43,  # chosen to get times similar to experiment
        }
    rouse_params = material_params.copy()
    rouse_params['rx'] = rouse_params['ry'] = rouse_params['rz'] = rouse_params['R']
    del rouse_params['R']

    # set up scan parameters (really just FP, the others are just "constants",
    # but we have to pass them all in a single parameter to get around a bug in
    # multiprocessing (2018-10-01)
    params = {'FP': fp_list,
              'tether_list': [np.array([]).astype(int)],
              'output_dir': [base_dir],
              'time_params': [time_params],
              'rouse_params': [rouse_params]}
    scan = Scan(params)
    scan.add_count(lambda p: replicates)

    # set up multiprocessing
    num_cores = multiprocessing.cpu_count() - 1
    p = multiprocessing.Pool(num_cores)

    # now run simulations, one per core until complete
    script_name = os.path.basename(__file__)
    print(script_name + ': Running scan!')
    for _ in p.imap_unordered(save_interior_sim, scan.params(), chunksize=1):
        continue
    print(script_name + ': Completed scan!')

    # # the single-threaded version of the above code
    # for params in scan.params():
    #     save_interior_sim(params)
    #     print(script_name + ": " + datetime.datetime.now().isoformat()
    #           + ": completed run with params: " + str(params))

if __name__ == '__main__':
    import sys
    base_dir = sys.argv[1] if len(sys.argv) > 1 else None
    run_homolog_param_scan(base_dir=base_dir)
