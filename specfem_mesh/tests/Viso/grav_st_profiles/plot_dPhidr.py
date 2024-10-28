# Check d Phi_st d r against FD: 
import numpy as np 
import matplotlib.pyplot as plt 

fig, ax = plt.subplots(3,3)


for s in range(1,4): 
    for t in range(-s, s+1):

        d = np.loadtxt(f"{s}_{t}.txt")



        ax[s-1, 0].plot(d[:,0], d[:,2], '-')

        # Grad phi from the data
        ax[s-1, 2].plot(d[:,0], d[:,3], '-')

        # FD 
        fd = (d[2:,2] - d[:-2,2])/ (2*(d[1,0]-d[0,0]))
        ax[s-1, 1].plot(d[1:-1,0], fd)


s = 2
t = -1
m = np.loadtxt(f'mesh_{s}_{t}.txt')

v = 10 
ax[s-1, 0].plot(m[::v,0], m[::v,2], 'x', markersize=2)
ax[s-1, 2].plot(m[::v,0], m[::v,3], 'x', markersize=2)


plt.savefig('profiles.pdf', format='pdf')