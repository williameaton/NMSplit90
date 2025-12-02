import numpy as np 
import matplotlib.pyplot as plt 
import scipy.special as ss 



fig, ax = plt.subplots()

for r in np.array([1221.5, 2891, 6371])/6371:

    # Draw circle: 
    circle_x = []
    circle_y = []
    ellipse_x = []
    ellipse_y = []
    colat = []

    eps = 0.1

    # phi is colatitude
    for phi in np.linspace(0, 2*np.pi, 1000):
        
        circle_x.append(r*np.cos(phi))
        circle_y.append(r*np.sin(phi))

        # Compute colatitude: 
        if phi < np.pi/2: 
            colatitude = np.pi/2 - phi 
        elif np.logical_and(phi >= np.pi/2, phi < np.pi):
            colatitude = phi - np.pi/2
        elif np.logical_and(phi >= np.pi, phi < np.pi*1.5):
            colatitude = phi - np.pi/2
        else: 
            colatitude = ((2*np.pi) - phi) + np.pi/2

        deltad = -(2/3)*r*eps*ss.eval_legendre(2, np.cos(colatitude))

        ellipse_x.append((r+deltad)*np.cos(phi))
        ellipse_y.append((r+deltad)*np.sin(phi))


    ax.plot(circle_x, circle_y, 'k')
    ax.plot(ellipse_x, ellipse_y, 'r')
fig.savefig('circle.pdf', format='pdf')
