import numpy as np 
import matplotlib.pyplot as plt 


g = np.loadtxt('gmag.txt')

fig, ax = plt.subplots()

ax.plot( 6371 - g[:,0]/1000 , g[:,1], 'k' )


ax.set_xlim([5200, 6371])
ax.set_ylabel('gravity, g [m/s^2]')
ax.set_ylim([0, 15])
ax.set_xlabel('depth [km] ')


c = np.loadtxt('check.txt')

ax.scatter( 6371 - c[:,0]/1000 , c[:,1], s=3  )


plt.savefig('gmag.pdf', format='pdf')
