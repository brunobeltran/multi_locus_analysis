"""run some BD simulations using our lab's `wlcsim` module to compare to the
data."""
import os
import multiprocessing
import socket # for getting hostname
from pathlib import Path
import datetime
import pickle

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from bruno_util import numpy as bnp
from wlcsim.bd import rouse
from pscan import Scan

from ... import finite_window as fw

N = int(1e2+1); L = 17.475; R = 1.3; b = 0.015; D = 2.8 # ChrV
Aex = 100; # strength of force confining beads within nucleus
dt = rouse.recommended_dt(N, L, b, D)
ura3_bead = 20;

# 10 repeats each of FP=[0, 0.02, 0.05, 0.1] on one 32 core node takes about
# ~1min to run for Nt=1e5,
# ~14min for Nt=1e6
# ~3.2hrs for Nt=1e7
# TBD for Nt=1e8
Nt = 1e8  # total number of equi-space time steps, about 2500s
Nlin = 1e3  # how many time steps per log-spaced group
overlap = 0.43  # chosen to get times similar to experiment

def run_interior_sim(FP):
    """for various different "connectivities" (fraction of beads "homolog
    paired") run 100 or so BD simulations each to get statistics on looping
    probabilities that are comparable to the experiment."""

    # FP = 0.1 # fraction loci homolog paired
    tether_list = np.array([]).astype(int) # no tethered loci

    # now define the time grid on which to run the BD, and specify which of those times to save
    t = np.arange(0, Nt*dt, dt)  # take Nt time steps
    t_i, i_i = bnp.loglinsample(Nt, 1e3, 0.43)  # contains ~30,60,90,etc
    # i30 = np.argmin(np.abs(t - 30)) # index at "30s"
    # t_save = t[0::i30] # save every 30s
    t_save = t[t_i]

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

    # run and save simulation
    df = run_interior_sim(FP)
    df.to_csv(sim_dir / Path('all_beads.csv'), index=False)

    # make sure to save parameters that were used
    pd.Series({
        'N': N, 'L': L, 'R': R, 'b': b, 'D': D, 'Aex': Aex, 'dt': dt,
        'Nt': Nt, 'Nlin': Nlin, 'linlog_overlap': overlap
    }).to_csv('params.csv')

def get_bead_df(base_dir, bead_id=ura3_bead, force_redo=False):
    """Extracts all ura3 data from simulations in a set of directories."""
    dfs = []
    for sim in base_dir.glob('homolog-sim.*'):
        sim_file = sim / Path('all_beads.csv')
        if sim_file.exists() and not force_redo:
            bead = pd.read_csv(sim_file)
            dfs.append(bead)
            continue
        df['sim_name'] = sim
        bead = df[df['bead'] == ura3_bead].copy()
        t0 = df[df['t'] == df['t'].iloc[0]]
        loops = t0.loc[t0.is_loop == 1, 'bead'].values
        left_n = np.where(loops <= ura3_bead)[0]
        bead['left_neighbor'] = loops[left_n[-1]] if len(left_n) > 0 else np.nan
        right_n = np.where(loops >= ura3_bead)[0]
        bead['right_neighbor'] = loops[right_n[0]] if len(right_n) > 0 else np.nan
        bead['min_bead'] = t0['bead'].min()
        bead['max_bead'] = t0['bead'].max()
        dfs.append(bead)
        dfs[-1].to_csv(sim_file)

    df = pd.concat(dfs, ignore_index=True)
    df = df.set_index(['FP', 'sim_name', 'bead', 't'])
    df = df.sort_index()
    df = add_paired_cols(df)
    df.to_csv(base_dir / Path(f'all_bead_{bead_id}.csv'))
    return df

def select_exp_times(base_dir, Nt=1e8, bead_id=ura3_bead, t=None,
                     force_redo=False):
    base_dir = Path(base_dir)
    bead_file = base_dir / Path(f'all_bead_{bead_id}.csv')
    if not bead_file.exists() or force_redo:
        df = get_bead_df(base_dir)
    else:
        df = pd.read_csv(bead_file)
    # # find the ones that look like experiment times
    # desired_t = np.arange(0, 1501, 30)
    # sim_t = df['t'].unique()
    # err = np.array([np.min(ti - sim_t) for ti in desired_t])
    # for now instead just use the biggest available step size
    with open(base_dir / 'i_i.pkl', 'rb') as f:
        i_i = pickle.load(f)
    # use the spacing that's biggest while still having enough samples
    i_i_i = np.where([len(i) > 50 for i in i_i])[0][-1]
    if t is None:
        # we've used same dt for all simulations so far
        t = np.arange(0, Nt*dt, dt) # so use global value by default
    df_exp = df[np.isin(df['t'], t[t_i[i_i[i_i_i]]])]
    return df_exp

def add_paired_cols(df, paired_distances=None):
    if paired_distances is None:
        # paired_distances = [10, 50, 100, 250, 500, 750, 1000]
        paired_distances = [250]
    df['dX'] = np.sqrt(
            np.power(df['X2'] - df['X1'], 2)
            + np.power(df['Y2'] - df['Y1'], 2)
            + np.power(df['Z2'] - df['Z1'], 2)
    )
    for dist in paired_distances:
        df['pair' + str(dist)] = df['dX'] < dist
    return df

def get_interior_times(df, state_col='pair250', TOL=None):
    waitdf = df.groupby(['FP', 'sim_name']).apply(
            fw.discrete_trajectory_to_wait_times, t_col='t', state_col=state_col)
    def interior(df):
        return df.reset_index().iloc[1:-1]
    waitdf = waitdf.groupby(['FP', 'sim_name']).apply(interior)
    # for some reason the index gets duplicated (reset_index above makes easier
    # to delete duplicated columns, since del is easier than droplevel)
    # waitdf.index = waitdf.index.droplevel(0)
    del waitdf['FP']
    del waitdf['sim_name']
    # also because we're using floating times, we need to de-duplicate
    # choose TOL  to be two more decimal points than dt, about
    if TOL is None:
        TOL = dt/1e2
    wtimes = np.sort(waitdf['wait_time'].unique().copy())
    diff = np.append(True, np.diff(wtimes))
    wtimes = wtimes[diff > TOL]
    for uniq_time in wtimes:
        waitdf.loc[np.isclose(waitdf['wait_time'], uniq_time), 'wait_time'] = uniq_time
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
                           base_name=None):
    if base_name is None:
        # save each run of this script in a unique directory
        base_name = './homolog-sim/no-tether-more-saves'
    # run_num = 0
    # while True:
    #     # as long as the OS is sane, only one thread can successfully mkdir
    #     try:
    #         # periods aren't allowed in hostnames, so use as delimiter
    #         base_dir = base_name + '.' + str(run_num)
    #         os.makedirs(base_dir, mode=0o755)
    #         break
    #     except OSError:
    #         run_num = run_num + 1

    # jk no, don't do that. much more useful to be able to re-run script and
    # just get more replicates
    base_dir = base_name

    # set up parameters (FP, tether lists, output directories)
    params = {'FP': fp_list,
              'tether_list': [np.array([]).astype(int)],
              'output_dir': [base_dir]}
    scan = Scan(params)
    scan.add_count(lambda p: replicates)

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

if __name__ == '__main__':
    run_homolog_param_scan()
