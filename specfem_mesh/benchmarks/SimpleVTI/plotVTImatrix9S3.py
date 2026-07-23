import numpy as np 
import matplotlib.pyplot as plt 
from broc import broc_map
import matplotlib.gridspec as gridspec

# Plot Figure 9S3 as a matrix and in graph format 
NEX = 176 

fsticks = 12
fslabels= 14

fig= plt.figure(figsize=(12,6))
spec = gridspec.GridSpec(ncols=2, nrows=10, figure=fig)

axmat  = fig.add_subplot(spec[:8, 0])
axcart = fig.add_subplot(spec[:, 1])
axcbar = fig.add_subplot(spec[-1, 0])


N = 11
L = 4
tl1 = 2*L + 1 
name = f"{N}S{L}"

# in mHz
fmz   = np.loadtxt(f"PREMfreqs/{name}.txt", skiprows=1) * 1000
omega = fmz * 2 * np.pi/1000



# Lets convert them to their real amplitudes: by division by freq
sem = np.loadtxt(f"SEM_NEX{NEX}/sem_fast_{name}.txt")[:tl1, :tl1]/(2*omega) 
t95 = np.loadtxt(f"Tromp95/radial_{name}.txt")[:tl1, :tl1]/(2*omega)

m  = np.arange(tl1) - L
mlabs  = m[::2]


# Plot diagonals 
semdiag = sem.diagonal()*1e5
t95diag = t95.diagonal()*1e5

axcart.plot(m, t95diag, 'k', zorder=0, linewidth=2)
axcart.scatter(m, semdiag, marker='x', color='#BB5566', s=60, linewidth=2)

axcart.set_xticks(mlabs)
axcart.set_xticklabels(mlabs, fontsize=fsticks)
axcart.set_yticks([-12, -8 , 4, 0, 4])
axcart.set_yticklabels([-12, -8 , 4, 0, 4], fontsize=fsticks)

axcart.set_xlabel(r"$m$", fontsize=fslabels)
axcart.set_ylabel(r"$\delta \omega $ [x 10$^{\text{-5}}$ rad s${}^{\text{-1}}$]", fontsize=fslabels)

# Plot the actual matrix: 
# We already say x 1e-5 in the cbar label so remove from the cbar by 
sem *= 1e5
vv = np.amax(sem)

im = axmat.imshow(sem, cmap=broc_map, vmin=-vv, vmax=vv)


ticks = np.arange(tl1)[::2]
axmat.set_xticks(ticks)
axmat.set_xticklabels(mlabs, fontsize=fsticks)
axmat.set_yticks(ticks)
axmat.set_yticklabels(mlabs, fontsize=fsticks)


leg = axcart.legend(["Radial eigenfunctions", "SEM"], fontsize=12)
frame = leg.get_frame()
frame.set_edgecolor('k')

axcart.set_ylim([-12, 5])


axcart.spines['bottom'].set_bounds(-4, 4)
axcart.spines[['right', 'top']].set_visible(False)

fig.colorbar(im, cax=axcbar, orientation='horizontal')


axcbar.set_xticks(np.array([-4, -2, 0, 2, 4]) )
axcbar.set_xticklabels(np.array([-4, -2, 0, 2, 4]), fontsize=fsticks)


axcbar.set_xlabel(r"$\delta \omega $ [x 10$^{\text{-5}}$ rad s${}^{\text{-1}}$]", fontsize=fslabels)

fig.savefig(f"{N}S{L}VTI.pdf", format='pdf')
