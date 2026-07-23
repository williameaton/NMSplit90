import numpy as np 
import matplotlib.pyplot as plt 
from numpy import genfromtxt
import numpy.polynomial.polynomial as poly 

fig, ax = plt.subplots()


n = 9
l = 3

tl1 = 2*l + 1 

outdir = "/scratch/gpfs/TROMP/we3822/NMSplit90/specfem_mesh/output/"

# Load the SEM matrix: 
name = f"{n}S{l}"
Vani = np.loadtxt(f"{outdir}/sem_fast_9S3.txt")[:tl1, :tl1].diagonal()

# This is the dimensionalised vani but for its contribution to H it is 
# 1/(2omega_0) 
# Load the freq, which for some reason is in Hz! 
fHz = np.loadtxt(f"{outdir}/freqs/{name}.txt", skiprows=1)


fmHz = fHz * 1e3


omega = 2*np.pi*fHz


# Perturbation, in angular freq
H = Vani/(2*omega)


# perturbation in mHz
H_mHz = 1000 * H/(2*np.pi)



m = np.arange(tl1) - l

ax.plot(m, fmHz + H_mHz, 'xk')




# Now get the Tromp data:
# With rotation + ani 
all93 = genfromtxt('9S3_points/poly_all.csv', delimiter=',')
rot93 = genfromtxt('9S3_points/poly_rot.csv', delimiter=',')

# Convert to polynomial: 

pall = poly.Polynomial(np.polyfit(all93[:,0], all93[:,1], deg=4)[::-1] )
prot = poly.Polynomial(np.polyfit(rot93[:,0], rot93[:,1], deg=4)[::-1] )


x = np.linspace(-l, l, 100)



#ax.plot(x, pall(x))
#ax.plot(x, prot(x))
#ax.plot(all93[:,0], all93[:,1], 'x')
#ax.plot(rot93[:,0], rot93[:,1], 'x')


# Vani contribution from Tromp paper
vani_93 = pall(x) - prot(x) + fmHz

coeff = 0.0067
ax.plot(x, vani_93+coeff)


fig.savefig("9S3.pdf", format='pdf')


