import numpy as np
import matplotlib.pyplot as plt
from wetools.plotting import setup_we_mpl; setup_we_mpl()
import matplotlib.gridspec as gridspec
from colour_schemes import hex

HC = hex.Hexes().hex['MedContrast_PT'][::-1]

nexs = ['64', '96', '112', '128', '176']

n = 6
t = "S"
l = 10

tl1 = 2*l + 1
fig = plt.figure(constrained_layout=True)
spec2 = gridspec.GridSpec(ncols=2, nrows=2, figure=fig)
ax_sem   = fig.add_subplot(spec2[0, 0])
ax_fast  = fig.add_subplot(spec2[1, 0])
ax_split = fig.add_subplot(spec2[:, 1])

ax = [ax_sem, ax_fast]


mode = f"{n}{t}{l}"

radial = np.diag(np.loadtxt(f"./radial_{mode}.txt")[:tl1, :])

m = np.arange(-10, 11)



ax_split.plot(m, radial, color='k')

inex = 0
for nex in nexs:
    times = np.loadtxt(f"{mode}_{nex}.txt", skiprows=1)

    sem      =  np.diag(np.loadtxt(f"{nex}/sem_{mode}.txt")[:tl1,:])
    fast     =  np.diag(np.loadtxt(f"{nex}/sem_fast_{mode}.txt")[:tl1,:])

    ax_sem.scatter(np.zeros(tl1)+times[1], 100*abs((sem  - radial)/radial), marker='o', color=HC[inex+1], label=nex)
    ax_fast.scatter(np.zeros(tl1)+times[2], 100*abs((fast - radial)/radial), marker='o', color=HC[inex+1], label=nex)


    ax_split.plot(m, sem, 'x', markersize='8' , color=HC[inex+1], label=nex)

    inex += 1






ax_fast.set_xlabel('Time [s]')
ax_sem.set_xlim([0, 200])
ax_fast.set_xlim([0, 800])
ax_split.set_xlim([-11,11])

for axi in [ax_sem, ax_fast]:
    axi.set_ylabel('% difference from Radial')
    axi.set_ylim([0,6])

for axi in [ax_sem, ax_fast, ax_split]:
    h, l = axi.get_legend_handles_labels()
    ph = [axi.plot([],marker="", ls="")[0]]*1
    handles = ph[:1] + h[::1] + ph[1:] + h[1::1]
    labels = ["  NEX"] + l
    leg = axi.legend(handles, labels, ncol=1)
    for vpack in leg._legend_handle_box.get_children():
        for hpack in vpack.get_children()[:1]:
            hpack.get_children()[0].set_width(0)

    ax_sem.set_title('Diagonal matrix elements', weight='bold', fontsize=11)
    ax_fast.set_title('All matrix elements', weight='bold', fontsize=11)

fig.suptitle(r"Mode ${}_{6}S_{10}$", weight='bold')
plt.show()






plt.show()