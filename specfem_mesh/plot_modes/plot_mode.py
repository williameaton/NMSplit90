import matplotlib.pyplot as plt
import numpy as np
from wetools.plotting import setup_we_mpl

# Load mode:
ddir = './modes/'

n = 0
mode_type = 'C'
l = 2

d = np.loadtxt(f"{ddir}/{n}{mode_type}{l}.txt", skiprows=1)
rad = d[:,0]

naxes = np.shape(d)[1]-1

setup_we_mpl()
fig, ax = plt.subplots(1,naxes, sharey=True)

for i in range(naxes):
    ax[i].plot(d[:,i+1], rad, 'k')

fig.suptitle(f"{n}{mode_type}{l}")

plt.show()