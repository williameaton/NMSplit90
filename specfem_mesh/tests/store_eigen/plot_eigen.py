# Plots a list of eigenfunctions: 
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.backends.backend_pdf


modeNs = [2, 5, 6, 7, 8, 21, 7, 9, 3, 9, 9, 11, 11, 13, 13, 13, 13, 15, 15, 18, 18, 20, 21, 25, 27, 21, 16]
modeLs = [3, 3, 3, 4, 5,  7, 5, 2, 2, 3, 4,  4,  5,  1,  2,  3,  6,  3,  4,  3,  4,  1,  6,  2,  2,  8,  7]

nmodes = len(modeLs)

ddir = './'
figlist = []
for i in range(nmodes):
    fig, ax  = plt.subplots(1,2, sharey=True)

    # Load eigenfuncs: 
    name = f"{modeNs[i]}S{modeLs[i]}"
    print(name)
    e = np.loadtxt(f"{ddir}/{name}_eigens.txt")

    rad  = e[:,0]
    u    = e[:,1]
    v    = e[:,2]
    du   = e[:,3]
    dv   = e[:,4]

    ax[0].plot(u, rad, 'k')
    ax[1].plot(v, rad, 'k')

    ax[1].set_title(name)
    figlist.append(fig)


pdf = matplotlib.backends.backend_pdf.PdfPages("eigens.pdf")
for f in figlist:
    pdf.savefig( f )
pdf.close()