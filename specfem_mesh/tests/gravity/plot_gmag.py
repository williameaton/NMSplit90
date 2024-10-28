import numpy as np 
import matplotlib.pyplot as plt 


g = np.loadtxt('gmag.txt')

fig, ax = plt.subplots(2)

depth = 6371 - g[:,0]/1000
ax[0].plot( depth , g[:,1], 'k' )

dr = depth[0] - depth[1]


ax[0].set_xlim([0, 6371])
ax[0].set_ylabel('gravity, g [m/s^2]')
ax[0].set_ylim([0, 15])
ax[0].set_xlabel('depth [km] ')


#c = np.loadtxt('check.txt')
#ax[0].scatter( 6371 - c[:,0]/1000 , c[:,1], s=3  )


# Plot grav as fd 
gradg = (g[1:,1] - g[:-1,1])/dr

ax[1].plot(depth[1:] , gradg, 'k' )



plt.savefig('gmag.pdf', format='pdf')
