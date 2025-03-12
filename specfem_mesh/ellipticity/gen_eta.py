# Load in density profile:
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

import scipy
mineos_model = np.loadtxt("../../databases/prem_ani_att_database/model")

r   = mineos_model[:,1]/6371000
rho = mineos_model[:,2]

# Lower sides
disc = np.array([0, 33, 66, 71, 126, 130, 132, 139, 148, 154, 157, 168, 181, 185])


# integrate piecewise:
nsecs = len(disc)-1

intnum = np.zeros(nsecs)
intden = np.zeros(nsecs)

fig, ax = plt.subplots(3)


num = rho*(r**4)
den = rho*(r**2)


# Integrate the densitt multiplied by r to the 4:
cumnum = scipy.integrate.cumtrapz(y=num, x=r, initial=0)
cumden = scipy.integrate.cumtrapz(y=den, x=r, initial=0)


ax[0].plot(r, cumnum)
ax[0].plot(r, r*r*cumden)

eta = 6.25*(1 - cumnum/(r*r*cumden))**2 - 1


# Remove the nan from r = 0
eta[0] = 0

#ax[1].plot(r, eta)
# For the smallest 500 km let the eta be linearly tapered to 0:
v = 25
eta[:v] *= ((eta[v-1]-eta[0])/(r[v-1]- r[0]) )
v = 10
eta[:v] *= ((eta[v-1]-eta[0])/(r[v-1]- r[0]) )
ax[1].axvline(r[v-1])

ax[1].plot(r, eta)
ax[1].set_ylabel('Eta')

eta_r = - eta/r



nradial = len(r)
epsilon = np.zeros(nradial)

for ir in range(nradial):

    # Select r range and integrand range
    integrand = eta_r[ir:]
    rrange    = r[ir:]
    integral = scipy.integrate.trapz(y=integrand, x=rrange)
    epsilon[ir] = np.exp(integral)

# Scale by Ea (e.g. hydrostatic is 1/299.8)
Ea = 1/299.8

epsilon *= Ea

epsilon[0] = epsilon[1]


ax[2].plot(r, epsilon)
ax[2].set_ylabel('Epsilon')


#ax[1].set_ylim([0, 0.6])
ax[2].set_ylim([0.0022, 0.0034])




# We may now interpolate these values: these values
# Load the radial values desired: 
rf90 = np.loadtxt("radialpoints.txt")


etainterp = np.zeros(len(rf90))
epsinterp = np.zeros(len(rf90))

for isec in range(nsecs): 
    rlower = r[disc[isec]]
    rupper = r[disc[isec+1]-1]

    r90mask = np.logical_and(rf90 >=rlower, rf90<=rupper)

    rslice   =  r[disc[isec] : disc[isec+1] ]
    etaslice = eta[disc[isec] : disc[isec+1]]
    epsslice = epsilon[disc[isec] : disc[isec+1]]

    etainterp[r90mask]= CubicSpline(rslice, etaslice)(rf90[r90mask])
    epsinterp[r90mask]= CubicSpline(rslice, epsslice)(rf90[r90mask])





ax[1].plot(rf90, etainterp, 'x')
ax[2].plot(rf90, epsinterp, 'x')

fig.savefig('eta.pdf', format='pdf')

np.savetxt(X=etainterp, fname='eta_interpolated.txt')
np.savetxt(X=epsinterp, fname='epsilon_interpolated.txt')

plt.show()
