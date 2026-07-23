# Looking at debugging the a and c values. 
# A should be easiest
import numpy as np 
import matplotlib.pyplot as plt 

#DT98 14.55 givs it as 
#a = 1/3(...) + 1/(2ω^2) ( v - tau*ω*ω)
# where 

# Load r, rho etc
rad  = np.loadtxt("radialpoints.txt")
rho  = np.loadtxt("rhopoints.txt")/5514.3
epsi = np.loadtxt("epsilon_interpolated.txt")
eta  = np.loadtxt("eta_interpolated.txt")


l = 2
k   = np.sqrt(l*(l+1))
ksq  = k**2 

u = np.loadtxt("eigen_u.txt")
v = np.loadtxt("eigen_v.txt")/k

Tbar   = - 6 * u * v 
Tcarot = u*u + (ksq - 3)*v*v


print(u)
print(v)

exit()


prefix = ksq/(( 2*l + 3)* (2*l -1)) 

integrand = (2/3) * epsi * rho * (Tbar - Tcarot*(eta+3)) * rad*rad 

print(len(integrand), integrand[:10])

tau = prefix * np.trapz(y=integrand, x=rad)


print(-tau/2 * 1000)


fig, ax = plt.subplots(1,5)

ax[0].plot(rho, rad)
ax[1].plot(u, rad)
ax[2].plot(v, rad)

fig.savefig("debug.pdf", format="pdf")