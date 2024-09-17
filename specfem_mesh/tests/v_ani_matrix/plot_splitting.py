import numpy as np
import matplotlib.pyplot as plt
from wetools.plotting import setup_we_mpl

setup_we_mpl()
from wetools.plotting import save_figs_to_single_pdf

plot_SEM = True

mdir = './Tromp_1993_model'
ddir = "/matrices/"
t = 'S'
ells = [2]
ns = [6]

nmodes = len(ells)
assert (len(ns) == len(ells))

figures = []

for imode in range(nmodes):

    fig, ax = plt.subplots()
    fig.set_tight_layout(True)

    l = ells[imode]
    n = ns[imode]

    rad = np.loadtxt(f'{mdir}/{ddir}/radial_{n}{t}{l}.txt')[:2 * l + 1, :]
    rad = np.diag(rad)

    m = np.arange(-l, l + 1)
    ax.plot(m, rad)

    if plot_SEM:
        sem = np.loadtxt(f'{mdir}/{ddir}/sem_{n}{t}{l}.txt')[:2 * l + 1, :]
        sem = np.diag(sem)
        ax.plot(m, sem)
        ratio = sem / rad
        print('ratio', ratio)

    ax.legend(['Modes', 'SEM'])

    ax.set_xlabel('m')
    ax.set_ylabel(r'$\delta \omega$')
    ax.set_xlim([-(l + 0.5), l + 0.5])
    ax.set_title(f'{n}{t}{l}', weight='bold')

    figures.append(fig)

#save_figs_to_single_pdf(figures, f'modes_{n}{t}{l}.pdf')
save_figs_to_single_pdf(figures, f'{mdir}/{n}{t}{l}_Tromp_model.pdf')
