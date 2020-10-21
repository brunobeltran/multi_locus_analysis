Now that we have estimates for the Kuhn length of the polymer, its total length,
and the average radius of the nucleus, all of the relevant length scales of the
problem are fixed.

It would be difficult to extract any relevant time scales from the single-cell
trajectories, since they are largely plateaus even at the shortest time scales
of our experiment (:math:`t = 30s`). Therefore, in order to determine the time
scale, we can look to the ensemble averaged MSCD curves.

Since all of the length scales of the problem are now fixed, given a model for
linkage formation along the chromosome, we can extract the average linkage
spacing by simply looking at the plateau levels of ensemble averaged MSCD
curves.

While we are only measuring the loci that are more than 250nm apart from each
other in our experiments, this same procedure can be repeated for comparable
simulations to extract how large of a bias it should cause, then invert it.

With the simulation output from running `burgess.simulations.run_interior_sim`
(with post-processing via `burgess.simulations.get_bead_df` and
`burgess.simulations.select_exp_times`) saved in a file (``df_exp.csv`` below),
we can compute the simulation MSCD curves as follows:

.. code:: python

    import pandas as pd
    import multi_locus_analysis as mla
    import numpy as np
    import matplotlib.pyplot as plt

    df = pd.read_csv('doc/source/df_exp.csv')
    df['t'] = (np.round(df['t'] / 30)*30).astype(int)
    df['dX'] = df['X2'] - df['X1']
    df['dY'] = df['Y2'] - df['Y1']
    df['dZ'] = df['Z2'] - df['Z1']

    all_dvel = df.groupby(['FP', 'sim_name']) \
                 .apply(mla.stats.pos_to_all_vel,
                        xcol='dX', ycol='dY', zcol='dZ', framecol='t')

    mscd_fp = dV.groupby(['FP', 'delta']).agg(['mean', 'std', 'count'])

    def plot_mscd(df):
        df = df.reset_index()
        fp = df['FP'].iloc[0]*100
        df = df.sort_values('delta')
        df = df[df['delta'] > 0]
        plt.errorbar(df['delta'], df['mean'], df['std']/np.sqrt(df['count']), c=sim_cmap(sim_cnorm_continuous(fp)))

    mscd_fp.groupby(['FP']).apply(plot_mscd)
    plt.yscale('log')
    plt.xscale('log')
    plt.colorbar(sim_sm_continuous)
    plt.xlabel('time (s)')
    plt.ylabel('Ensemble MSCD ($\mu{}m^2$)')

which will produce a plot like the following:

(TODO, copy in)

.. .. figure::

but to see the effects of extra , we need to redo...
