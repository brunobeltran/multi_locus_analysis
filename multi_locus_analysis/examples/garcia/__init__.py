"""A project led by Myron Child to explore homologous loci around
gastrulation."""

chrom_name = '2L'
locus_name = '38F1'

# very approximate, in base pairs
loc = 20_820_000
chrom_len = 23_520_000
# bp to nanometers
bp_to_nm = 0.34

# guess from nucleosome paper (Beltran & Kannan, in review)
kuhn_length = 40 # (nm)

def plot_baseline(alpha):
    td = np.linspace(0, 3, 91)
    plt.figure()
    plt.plot(td, analytical.vc(td, np.ones_like(td), alpha), 'k-.',
            label=f'$\\beta = {alpha}$')
    alpha = alpha/2
    plt.plot(td, analytical.vc(td, np.ones_like(td), alpha), 'k--',
            label=f'$\\beta = {alpha}$')
    plt.plot([0,3], [0,0], 'k', label=f'uncorrelated')
    plt.ylim([-0.5, 1])
    plt.xlim([0, 3])
    plt.legend()
