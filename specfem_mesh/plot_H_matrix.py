import numpy as np 
import matplotlib.pyplot as plt 
from broc import broc_map


n = '6'
t = 'S'
l = '10'

tl1 = 2*int(l) + 1

mats = [f"VTI_sem_fast_{n}{t}{l}.txt", f"sem_fast_{n}{t}{l}.txt"]


fig, ax = plt.subplots(2,2)

minv = 100000
maxv = -100000

for i in range(len(mats)): 
    mat = np.loadtxt(f"output/{mats[i]}")

    print(minv)
    print(maxv)

    minval = np.min(mat)
    maxval = np.max(mat)
    if minval < minv: 
        minv = minval 
    if maxval > maxv:
        maxv = maxval


    print(minv)
    print(maxv)
    print()


maxval = np.max([np.abs(minv), np.abs(maxv)])


for i in range(len(mats)): 
    mat = np.loadtxt(f"output/{mats[i]}")

    real = mat[:tl1, :]
    imag = mat[tl1:, :]


    im = ax[i,0].imshow(real, vmin=-maxval, vmax=maxval, cmap=broc_map)
    im = ax[i,1].imshow(imag, vmin=-maxval, vmax=maxval, cmap=broc_map)

    fs = 7
    ax[i,0].text(x=0.55*tl1, y=1, s="Min:   " + '{:0.1e}'.format(np.min(real)), fontsize=fs)
    ax[i,0].text(x=0.55*tl1, y=2, s="Max:   " + '{:0.1e}'.format(np.max(real)), fontsize=fs)


    ax[i,1].text(x=0.55*tl1, y=1, s="Min:   " + '{:0.1e}'.format(np.min(imag)), fontsize=fs)
    ax[i,1].text(x=0.55*tl1, y=2, s="Max:   " + '{:0.1e}'.format(np.max(imag)), fontsize=fs)


ax[0,0].set_title('Real', weight='bold')
ax[0,1].set_title('Imaginary', weight='bold')


fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)


ax[0,0].set_ylabel(f"VTI model", weight='bold')
ax[1,0].set_ylabel(f"TTI model", weight='bold')

fig.suptitle(f"Mode {n}{t}{l}", weight='bold')

plt.savefig(f'figures/Hmat_{n}{t}{l}.pdf')