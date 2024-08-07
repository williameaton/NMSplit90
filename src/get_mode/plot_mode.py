import numpy as np 
import matplotlib.pyplot as plt 
from wetools.plotting import setup_we_mpl

hdrs = {2: ['W', "W'"],
        6: ['U', "U'", 'V', "V'", 'P', "P'"]}

setup_we_mpl()
for eig_file in ['4S2.txt']:

    # Read header and load data
    f = open(eig_file, 'r')
    [n, typ, l] = f.readline().split()
    d = np.loadtxt(eig_file, skiprows=1)

    shp = np.shape(d)
    cols = shp[1]-1
    # Plot
    fig, ax = plt.subplots(1,cols)
    fig.set_tight_layout(True)


    for i in range(cols):
        ax[i].plot(d[:,i+1], d[:,0])
        ax[i].set_title(hdrs[cols][i])

    fig.suptitle(f"{n}{typ}{l}")


# Load from mode.data:
d = np.loadtxt('mode.data')


for i in range(4):
    ax[i].plot(d[:,i+1], d[:,0])

plt.show()