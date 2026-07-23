import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.backends.backend_pdf as backend

# Table values: 
# Note that 0S8 is wrong in DT98 and correct in DS79
bDS79 = np.array([14.905, 4.621, 1.834, 0.841, 0.407, 0.181, 0.054, 4.173, 2.633, 1.948, 1.437, 0.873, 0.564, 0.427, 0.349,
         0.668, 0.281, 0.159, 0.340, 1.657, 1.485, -0.050, 0.081, 0.019, 0.055, 0.013, 0.005, 0.102, 0.031, 0.020, 0.017])

# List of modes: 
modeNs = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 6, 8, 8, 9, 11, 11, 13, 13, 18, 18]
modeLs = [2, 3, 4, 5, 6, 7, 8, 2, 3, 4, 5, 6, 7, 8, 9, 3, 4, 5, 6, 1, 2, 3, 1, 5, 3,  4,  5,  1,  2,  3,  4]
Types  = "S"

semianalyB = []
semianalyB1066a = []
SEM_bs = []

nmodes = len(modeNs)

xtickslabels = []

for i in range(1): 

    name = f"{modeNs[i]}S{modeLs[i]}"

    # Load semi analytical bs 
    semianalyB.append(np.loadtxt(f"./abc/{name}_{name}.txt")[1] * 1000.0)
    semianalyB1066a.append(np.loadtxt(f"./1066_abc/{name}_{name}.txt")[1] * 1000.0)

    # Load the Wmatrices from SEM and compute B: 
    # Note that what we are seeing is the total difference of the poynomial of order m
    # i.e. \delta w = w_0 (a + bm + cm^2) so we are getting w_0 b - need to divide by
    # the angular frequency (which needs to be nondimensionalised, not dimensionalised)
    wcom = 1.9438865 * 1e-3

    Wmat = np.diag(np.loadtxt(f"./Wmats/Wmat_{name}_{name}.txt")[:2*modeLs[i] + 1, :])
    SEM_bs.append(1e3*(Wmat[1] - Wmat[0])/wcom)



    xtickslabels.append(fr"${{}}_{{{modeNs[i]}}}S_{{{modeLs[i]}}}$")



semianalyB      = np.array(semianalyB)
semianalyB1066a = np.array(semianalyB1066a)

x = np.arange(1)

bDS79 = bDS79[0]

print(SEM_bs)

fig, ax = plt.subplots(2, figsize=(15, 5))
#ax.scatter(x, semianalyB)
#ax[0].scatter(x, semianalyB1066a)
ax[0].scatter(x, bDS79,  marker='o', fc='None', ec='k')
ax[0].scatter(x, SEM_bs, marker='x', color='r')


ax[1].scatter(x, np.abs((semianalyB1066a - bDS79)), marker='o', fc='k')
ax[1].scatter(x, np.abs((SEM_bs - bDS79)), marker='o', fc='k')



axperc = ax[1].twinx()
axperc.scatter(x, np.abs((semianalyB1066a - bDS79)/bDS79)*100, marker='o', fc='grey')
axperc.scatter(x, np.abs((SEM_bs - bDS79)/bDS79)*100, marker='o', fc='grey')



ax[0].set_xticks(x)
ax[0].set_xticklabels([])
ax[1].set_xticks(x)
ax[1].set_xticklabels(xtickslabels, rotation=45)



ax[1].set_ylabel("Absolute Error")
ax[0].set_ylabel(r"$b$ values")
axperc.set_ylabel("% error from DS79")

# figs = []

# for i in range(nmodes):
#     fig, ax = plt.subplots()
#     x = np.arange(3)
#     y = np.array([ semianalyB[i], semianalyB1066a[i], bDS79[i]])
#     ax.scatter(x,  y, marker='x')

#     #ax.scatter(x, semianalyB, marker='x')

#     maxa = np.max(y) * 1.01
#     ax.set_ylim([0, maxa])

#     figs.append(fig)
#     #fig.savefig(f"{i}.pdf", format='pdf')

# pdf = backend.PdfPages("DahlenSailor_B.pdf")
# for ff in figs: ## will open an empty extra figure :(
#     pdf.savefig( ff )
# pdf.close()
fig.set_tight_layout(True)
fig.savefig("DahlenSailor_B.pdf", format='pdf')