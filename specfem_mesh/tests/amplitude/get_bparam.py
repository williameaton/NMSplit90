import numpy as np 
import matplotlib.pyplot as plt 
import scipy

angfreqs = {'0S3': 2.9435544E-03, 
            '0T3': 3.6799361E-03}
SCALE_T = 930.12663175873513

OMEGA = 2*np.pi/86400 


# Load the Eigenfunctions and rho: 
ddir = 'eigens/'

rho = np.loadtxt(f"{ddir}/rho")
rad = rho[:,0]
rho = rho[:,1] * 5514.3 # original dimension

n = 4
t = 'S'
l = 3

k = (l*(l+1))**0.5

# Let us now work only in DIMENSIONALISED UNITS: 
w0 = angfreqs[f"{n}{t}{l}"]
w0_scaled = w0 * SCALE_T

if t == 'S': 
    U   = np.loadtxt(f"{ddir}/{n}{t}{l}_U")[:,1]
    V   = np.loadtxt(f"{ddir}/{n}{t}{l}_V")[:,1]

    fig, ax = plt.subplots()

    ax.plot(U, rad)
    ax.plot(V, rad)


    # Check the normalisation of the original eigenfunctions: 
    # Commented material shows that the modes are normalised according to mineos 
    """rho = rho/5515
    integrand = rho * ((U*U) + (V*V)) * rad * rad

    integral = scipy.integrate.trapz(y=integrand, x=rad)
    # Finally need to multiply by omega^2 but the omega needs to be scaled by time 
    om_scaled =  angfreqs[f"{n}{t}{l}"] * SCALE_T
    print('Test normalisation: ', om_scaled*om_scaled *integral)"""

    # Need to determine chi: 
    integrand = (rho/5515) *(rad*rad) * (V*V + 2*k*U*V)/ (k**2)
    chi = scipy.integrate.trapz(y=integrand, x=rad) 
    chi *= (w0_scaled*w0_scaled)

elif t=='T': 

    W  = np.loadtxt(f"{ddir}/{n}{t}{l}_W")[:,1]

    # Check if normalised as per minios to 1 
    #integral = scipy.integrate.trapz(y= rad*rad*W*W*rho/5515, x=rad)*(angfreqs[f"{n}{t}{l}"] * SCALE_T)**2


    chi_DT = scipy.integrate.trapz(y= rad*rad*W*W*rho/5515, x=rad)/(k**2) * (w0_scaled*w0_scaled)
    chi    = 1/(k**2)

    print('chi_DT: ', chi_DT)
    print('chi: ',     chi)


b = 1000 * chi * OMEGA/w0  
print('b is then ', b )

