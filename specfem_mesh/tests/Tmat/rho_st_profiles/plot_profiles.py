import numpy as np 
import matplotlib.pyplot as plt 



fig, ax = plt.subplots(1,3, sharey=True, sharex=True)

for i in range(3): 

    s = i+1

    for t in range(-s, s+1): 
        d = np.loadtxt(f'rhost_{s}_{t}')

        coeff = d[0,:]
        rad   = d[1:,0]
        rhost = d[1:,1]

        ax[i].plot(rhost, rad, label=f't = {t}')
    
    ax[i].set_title(f"s = {s}")
    ax[i].legend()

plt.savefig('profiles.pdf', format='pdf')