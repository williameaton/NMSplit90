import matplotlib.pyplot as plt
from numpy import genfromtxt
import numpy as np
import scipy.optimize as so


from determine_B_spline_coeffs import create_splines
RMAX = 0.19172814314864228

#r = np.linspace(0, 1, npts)
r = np.loadtxt("ACLNF/radial/r_for_radmethod")
npts = len(r)


sprnd = create_splines(r/np.max(r))

# load coefficients:
alpha_c = np.loadtxt(f"spline_coeffs/alpha")
beta_c  = np.loadtxt(f"spline_coeffs/beta")
gamma_c = np.loadtxt(f"spline_coeffs/gamma")

alpha = np.zeros(npts)
beta = np.zeros(npts)
gamma = np.zeros(npts)
for i in range(5):
    alpha  += alpha_c[i] * sprnd[i,:]
    beta   +=  beta_c[i] * sprnd[i,:]
    gamma  += gamma_c[i] * sprnd[i,:]

fig, ax = plt.subplots(1,2, figsize=(10,5), sharey=True)

red    = "#BB5566"
blue   = "#004488"
yellow = "#DDAA33"

ax[0].plot(alpha, r, color=yellow)
ax[0].plot(beta,  r, color=red )
ax[0].plot(gamma, r, color=blue )




# Create ACLNF maintaining isotropic values: 

def prem_vals_at_r(radius):
    # Compute PREM values at this radius:
    # Convert from g/cm^3  --> kg/m3
    #              km/s    --> m/s
    x2  = radius ** 2
    rho = (13.088500000 - 8.838100000 * x2) * 1000
    vp  = (11.262200000 - 6.364000000 * x2) * 1000
    vs  = ( 3.667800000 - 4.447500000 * x2) * 1000
    return rho, vp, vs

def ACLNF(rho, vp, vs):
    # Computes ACLNF from rho, vp, vs assuming eta = 1
    A = rho * vp * vp
    C = rho * vp * vp
    L = rho * vs * vs
    N = rho * vs * vs
    F = A - 2*L
    return [A,C,L,N,F]


# Get the PREM values 
rho, vp, vs = prem_vals_at_r(r)
[Aprem, Cprem, Lprem, Nprem, Fprem] = ACLNF(rho, vp, vs)


Kprem = Aprem - (4/3)*Nprem
Mprem = Nprem


A0 = (20*Mprem[0] + 15*Kprem[0]) / (15 + 3*alpha[0] - 4*gamma[0] + 8*beta[0] )

lam = A0*(alpha + 2*gamma + 6*beta)
sig = A0*(alpha - 4*gamma )

N = Mprem - lam/15
A = (9*Kprem + 12*N - sig)/9
C = Aprem + (A0/15) * (12*alpha - 8*beta + 4*gamma)

L = A0*beta + N 

F = A - 2*N - gamma*A0 

scale = 1e-9
ax[1].plot((A - Aprem)*scale, r, color=red)
ax[1].plot((C - Cprem)*scale, r, color='k')
ax[1].plot((L - Lprem)*scale, r, color="#6699CC")
ax[1].plot((N - Nprem)*scale, r, color=blue)
ax[1].plot((F - Fprem)*scale, r, color=yellow)


# Save the perturbations: 
odir = "ACLNF/radial/"
np.savetxt(odir+"A_0", A - Aprem)
np.savetxt(odir+"C_0", C - Cprem)
np.savetxt(odir+"L_0", L - Lprem)
np.savetxt(odir+"N_0", N - Nprem)
np.savetxt(odir+"F_0", F - Fprem)




loc = 'lower right'

fsleg = 11
leg0 = ax[0].legend([r"$\alpha$", r"$\beta$", r"$\gamma$"], loc=loc, fontsize=fsleg)
leg1 = ax[1].legend([r"$\delta A$", 
              r"$\delta C$", 
              r"$\delta L$",
              r"$\delta N$",
              r"$\delta F$"], loc=loc, fontsize=fsleg)


for legend in [leg0, leg1]:
    frame = legend.get_frame()
    frame.set_edgecolor('k')



ax[0].set_ylim([0,RMAX])

fs = 12
ax[0].set_ylabel("Fractional radius of inner core", fontsize=fs)

ax[0].set_xlabel("Anisotropy values", fontsize=fs)
ax[1].set_xlabel("Love parameter perturbation\nfrom PREM [GPa]", fontsize=fs)


for a in ax: 
    a.spines[['right', 'top']].set_visible(False)


ax[0].set_xlim([-0.03, 0.08])
ax[0].spines['bottom'].set_bounds(-0.02, 0.08)


ax[1].set_xlim([-50, 110])
ax[1].spines['bottom'].set_bounds(-40, 100)

fig.set_tight_layout(True)





fig.savefig("Tromp93_model.pdf")

plt.show()