import numpy as np
import matplotlib.pyplot as plt
from wetools.plotting import setup_we_mpl, subplots_hide_xaxes
import scipy.interpolate as ci


hdrs = {2: ['W', "W'"],
        6: ['U', "U'", 'V', "V'", 'P', "P'"]}

setup_we_mpl()

# Reproduce DT 98 figures:

fig, ax = plt.subplots(1,2)
fig.set_tight_layout(True)

icid = 33

# Read header and load data
fname = f"./6S2.txt"
f = open(fname, 'r')
[n, typ, l] = f.readline().split()
d = np.loadtxt(fname, skiprows=1)

# Plot v
#ax.plot(d[:icid, 3], d[:icid, 0], 'k--')

# Plot u
ax[0].plot(d[:icid,1], d[:icid,0], 'k-x')

ax[0].set_title('6S2')


spl = np.loadtxt('./uspline.txt')
splx = spl[:,0]

# Scipy interpolation
cs = ci.CubicSpline(x=d[:icid,0], y=d[:icid,1], bc_type='not-a-knot')
ueven = cs(splx)
ax[0].plot(ueven,splx, 'xr' )

# Load the uspline from f90


ax[0].plot(spl[:,1], splx, 'xb')

ax[1].scatter(splx, ueven-spl[:,1])


ax[0].set_ylabel('Normalised radius')
ax[0].set_xlabel('U eigenfunction')


ax[1].set_xlabel('Normalised radius')
ax[1].set_ylabel('Difference Scipy spline - my spline')

plt.show()