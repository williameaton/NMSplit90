import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.patches as patches



modeNs = [2, 5, 6, 7, 8, 21, 7, 9, 3, 9, 9, 11, 11, 13, 13, 13, 13, 15, 15, 18, 18, 20, 21,  27, 25, 21, 16]
modeLs = [3, 3, 3, 4, 5,  7, 5, 2, 2, 3, 4,  4,  5,  1,  2,  3,  6,  3,  4,  3,  4,  1,  6,   2,  2,  8,  7]

nmodes = len(modeLs)

# Load the data from SEM and Tromp 95 method 


fig, ax = plt.subplots(sharex=True)

NEX = 176

SEM = []
T95 = []
PER = []
FmHz = []
FmHzFull = []

names = [] 

clrs = []

for imode in range(nmodes): 

    N = modeNs[imode]
    L = modeLs[imode]
    tl1 = 2*L + 1 
    name = f"{N}S{L}"
    fancyname = fr"${{}}_{{{N}}} {{S}}_{{{L}}}$"

    names += [fancyname]

    # in mHz
    fmz = np.loadtxt(f"PREMfreqs/{name}.txt", skiprows=1) * 1000
    FmHz.append(fmz)

    for i in range(L+1): 
        FmHzFull += [fmz]

    omega = fmz * 2 * np.pi/1000


    # Lets convert them to their real amplitudes: by division by freq
    sem = np.loadtxt(f"SEM_NEX{NEX}/sem_fast_{name}.txt")[:tl1, :tl1].diagonal()/(2*omega) 
    t95 = np.loadtxt(f"Tromp95/radial_{name}.txt")[:tl1, :tl1].diagonal()/(2*omega)

    difference = sem - t95

    perc = 100*np.abs(difference/t95)

    maxperc = np.max(perc)

    for ii in range(L+1): 
        if perc[ii] == maxperc:
            clrs.append('#BB5566')
        else:
            clrs.append('k')


    if np.max(perc) > 2: 
        print(name)

    SEM += list(sem[:L+1])
    T95 += list(t95[:L+1])
    PER += list(perc[:L+1])


fig, ax = plt.subplots(2, sharex=True)

ax[0].axhline(1, color='grey', linestyle='--', alpha=0.5)

xmodes = range(nmodes)
ax[0].scatter(FmHzFull, PER, marker='o', s=4, c=clrs)
ax[0].set_yscale('log')


 
# plot in mHz
ax[1].scatter(FmHzFull, np.array(T95)*1000, marker='o', color='k', s=15)
ax[1].scatter(FmHzFull, np.array(SEM)*1000, marker='x', color='#BB5566', s=15,  linewidth=1)
ax[1].set_yscale('log')



# Sort the frequencies into increasing values: 
ifreq = np.argsort(FmHz)


modenameticks = np.array(names)[ifreq]


ax[1].set_xlabel("Frequency [mHz]")
ax[0].set_ylabel("% SEM misfit\nfrom radial kernels")

ax[1].set_ylabel(r"$\delta\omega$ [mHz]")


# For top ticks of mode names 
ax2 = ax[0].twiny()
ax2.set_xticks(np.array(FmHz)[ifreq])

print("Ticks: ", np.array(FmHz)[ifreq])
print("Ticks: ", modenameticks)

ax2.set_xticklabels(modenameticks, rotation=90, fontsize=9)


ticks = ax2.xaxis.get_major_ticks()

for i, tick in enumerate(ticks):
    if i % 2 == 0:
        tick.tick1line.set_markersize(2)  # long tick
        tick.tick2line.set_markersize(2)
    else:
        tick.tick1line.set_markersize(16)   # short tick
        tick.tick2line.set_markersize(16)

labels = ax2.get_xticklabels()
for i, label in enumerate(labels):
    if i % 2 == 1:
        x, y = label.get_position()
        label.set_position((x, y + 0.13))  # increase for more separation
    else:
        x, y = label.get_position()
        label.set_position((x, y - 0.03))  # increase for more separation

# Ensure same x lims 



xp = 8.93
ax[0].add_patch(patches.Rectangle((xp, 1e-8), 10-xp, 40, linewidth=0, edgecolor='None', facecolor='grey', alpha=0.1))
ax[1].add_patch(patches.Rectangle((xp, 1e-8), 10-xp, 1e1, linewidth=0, edgecolor='None', facecolor='grey', alpha=0.1))


ax[0].set_xlim([1,10])
ax2.set_xlim(ax[0].get_xlim())


leg = ax[1].legend(["Radial kernels", "SEM"], loc='lower right')

frame = leg.get_frame()
frame.set_edgecolor('k')
frame.set_facecolor('white')

fig.savefig("VTIconstant.pdf")

