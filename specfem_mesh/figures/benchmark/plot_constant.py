import numpy as np
import matplotlib.pyplot as plt
from wetools.plotting import setup_we_mpl, align_y_labels
from broc import broc_map
import matplotlib.gridspec as gridspec

setup_we_mpl()

datadir = './data/constant_ACLNF/'
mode = '6S4'
fig = plt.figure(2, figsize=(4.05,7))
fig.set_tight_layout(True)

spec2 = gridspec.GridSpec(ncols=1, nrows=16, figure=fig)
ax_cbar  = fig.add_subplot(spec2[:2])
ax_mat   = fig.add_subplot(spec2[2:9])
ax_comp  = fig.add_subplot(spec2[9:])


m = np.arange(-4, 5)


# Load the SEM matrix:
sem = np.loadtxt(f"{datadir}/sem_{mode}.txt")[:9]

map = ax_mat.imshow(sem, cmap=broc_map, vmin=-0.001, vmax=0.001)

semdiag = np.diag(sem)

fig.colorbar(cax=ax_cbar, mappable=map, orientation='horizontal')

# Load the radial data
rad = np.diag(np.loadtxt(f"{datadir}/radial_{mode}.txt")[:9])

ax_comp.plot(m, rad, label='Radial')
ax_comp.plot(m, semdiag, ':x', label='GLL Quadrature')


ax_mat.set_xlim([-0.5, 8.5])
ax_comp.set_xlim([-4.3, 4.3])

ax_comp.set_xlabel(r"$m$", weight='bold', fontsize=14)

ax_mat.set_ylabel(r"$m$", weight='bold', fontsize=14)
ax_comp.set_ylabel(r"$\delta\omega$", weight='bold', fontsize=14)

for axi in [ax_mat, ax_comp]:
    axi.yaxis.set_label_coords(-0.3, 0.5)


ax_comp.legend()


ax_cbar.set_xticks([-0.001, 0, 0.001])
ax_cbar.set_title('$\delta \omega$')


ax_mat.spines['top'].set_visible(True)
ax_mat.spines['right'].set_visible(True)
ax_mat.set_xticks([])
ax_mat.set_yticks([])





plt.savefig('constant_ACLNF.pdf', format='pdf')

plt.show()