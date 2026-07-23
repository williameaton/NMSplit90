import matplotlib.pyplot as plt
from numpy import genfromtxt
import numpy as np
import scipy.optimize as so

RMAX = 0.19172814314864228

def create_splines(x, N=4):

    xmax = 1
    s    = len(x)
    nknots = N + 1
    knots = np.linspace(0, xmax, nknots)
    h = knots[1] - knots[0]
    splines = np.zeros((nknots, s))
    xminknots = np.zeros((nknots, s))

    for i in range(nknots):
        xminknots[i,:] = x - knots[i]

    v0 = [[0.25, 0, -6/4, 6/4], [-0.25, 0.75, -0.75, 0.25]]
    v1 = [[-0.5, 0, 3/2, 0], [0.75, -6/4, 0, 1], [-0.25, 0.75, -0.75, 0.25]]
    vj = [[0.25, 0,0,0],[-0.75, 0.75, 0.75, 0.25],[0.75, -6/4, 0, 1],[-0.25, 0.75, -0.75, 0.25]]
    vN_1 = [[0.25, 0, 0, 0],[-0.75, 0.75, 0.75, 0.25],[0.5, -1.5, 0, 1]]
    vN = [[0.25,0,0,0],[-0.25, 0.75, 0.75, 0.25]]

    xmk_h = xminknots/h


    for iknot in range(nknots):
        if iknot == 0:
            masks = []
            masks.append(np.logical_and(x >= knots[0], x <= knots[1]))
            masks.append(np.logical_and(x >= knots[1], x <= knots[2]))
            masks.append(np.logical_and(x >= knots[2], x <= knots[3]))
            for mm in range(2):
                for p in range(4):
                        splines[iknot, masks[mm]] += v0[mm][p]* xmk_h[mm, masks[mm]]**(3-p)

        elif iknot == 1:
            masks = []
            masks.append(np.logical_and(x >= knots[0], x <= knots[1]))
            masks.append(np.logical_and(x >= knots[1], x <= knots[2]))
            masks.append(np.logical_and(x >= knots[2], x <= knots[3]))

            for mm in range(3):
                for p in range(4):
                        splines[iknot, masks[mm]] += v1[mm][p]* xmk_h[mm, masks[mm]]**(3-p)


        elif iknot == N-1:
            masks = []
            masks.append(np.logical_and(x >= knots[N-3], x <= knots[N-2]))
            masks.append(np.logical_and(x >= knots[N-2], x <= knots[N-1]))
            masks.append(np.logical_and(x >= knots[N-1], x <= knots[N]))
            for mm in range(3):
                for p in range(4):
                        splines[iknot, masks[mm]] += vN_1[mm][p]* xmk_h[N+mm-3, masks[mm]]**(3-p)


        elif iknot == N:
            masks = []
            masks.append(np.logical_and(x >= knots[N - 2], x <= knots[N-1 ]))
            masks.append(np.logical_and(x >= knots[N-1],     x <= knots[N]))
            for mm in range(2):
                for p in range(4):
                    splines[iknot, masks[mm]] += vN[mm][p] * xmk_h[N + mm - 2, masks[mm]] ** (3 - p)

        else:

            masks = []
            for h in range(4):
                masks.append(np.logical_and(x >= knots[iknot+h-2], x <= knots[iknot+h-1]))

            for mm in range(4):
                for p in range(4):
                        splines[iknot, masks[mm]] += vj[mm][p] * xmk_h[iknot-2 + mm, masks[mm]]**(3-p)

    return splines


def functional(x, n, data, spl):

    pred = 0
    for i in range(n):
        pred += spl[i,:]*x[i]
    misfit = np.sum((data[:,0] - pred)**2)
    return misfit

def fit_splines_to_plot(csv, ax, N=4):

    gamma = genfromtxt(csv, delimiter=',')

    r = gamma[:, 1]
    sp = create_splines(r, N=N)
    # Get the residual
    res = so.minimize(functional, args=(N + 1, gamma, sp), x0=[0, 0, 0, 0, 0]).x

    #ax[0].plot(gamma[:, 0], gamma[:, 1], 'x')

    print('res is ', res)

    return res




fig, ax = plt.subplots(2, figsize=(4.5,8))

N = 4
iproc = 0

coeff = []
gi = 0

greeks = []

for graph in ['alpha', 'beta', 'gamma']:
    coeff = fit_splines_to_plot(f'abg/{graph}.csv', ax=ax, N=N)

    np.savetxt(fname=f'spline_coeffs/{graph}', X=coeff)

    # Load the unique radii:
    #out_r = np.loadtxt(f'./unq_r/unq_r{iproc}')
    out_r = np.loadtxt(f'./ACLNF/radial/r_for_radmethod')
    n_out_r = len(out_r)
    try:
        assert (np.max(out_r) - RMAX) <= 1e-7
    except: 
        raise ValueError(np.max(out_r), RMAX)

    out_r = out_r/RMAX

    sp_out = create_splines(out_r, N=N)


    opt = 0
    for i in range(N+1):
            opt += sp_out[i,:]*coeff[i]
            #ax[1].plot(sp_out[i,:], out_r)

    ax[0].plot(opt, out_r, '-', label=rf'$\{graph}$')

    greeks.append(opt)

ax[0].legend()



# We now have the three greek symbols (alpha,beta,gamma)

# The assumption here will be that A0 is 1
# N = A = 0
# in this case the C = alpha, L = beta, F = -gamma

A = np.zeros(n_out_r)
C = greeks[0]
L = greeks[1]
N = np.zeros(n_out_r)
F = -greeks[2]

ctr = 0
labs = 'ACLNF'
for y in [A, C, L, N, F]:
    if ctr == 3:
        ls = '--'
    else:
        ls = '-'
    ax[1].plot(y, out_r, ls, label=labs[ctr])
    ctr +=1

ax[1].legend()

moddir = 'ACLNF/radial'
np.savetxt(f"{moddir}/A_{iproc}", X=A)
np.savetxt(f"{moddir}/C_{iproc}", X=C)
np.savetxt(f"{moddir}/L_{iproc}", X=L)
np.savetxt(f"{moddir}/N_{iproc}", X=N)
np.savetxt(f"{moddir}/F_{iproc}", X=F)

for i in range(2):
    ax[i].set_xlabel('% perturbation')
    ax[i].set_ylim([0,1])

ax[0].set_ylabel('Normalised radius')

fig.set_tight_layout(True)


fig.savefig("TROMP_MODEL.pdf", format='pdf')

