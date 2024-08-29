import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as ci
#import pytest

def test_spline():

        hdrs = {2: ['W', "W'"],
                6: ['U', "U'", 'V', "V'", 'P', "P'"]}
        mode_type = 'S'
        N    = 23
        L    = 12
        icid = 33

        # Reproduce DT 98 figures:
        fig, ax = plt.subplots(2,4, sharey=True)

        # Read header and load data for the original eigenfunction
        fname = f"./spline/{N}{mode_type}{L}.txt"
        f = open(fname, 'r')
        [n, typ, l] = f.readline().split()
        # Data is in order [radius, U, U', V, V', P, P']
        d = np.loadtxt(fname, skiprows=1)[:icid, :]

        eig_r = d[:,0]
        assert N == int(n)
        assert L == int(l)
        assert typ == mode_type


        for eig_id in range(1,5):
                eigen = d[:,eig_id]

                # Plot original eigs
                ax[0, eig_id - 1].plot(d[:,eig_id], eig_r, 'k-x', alpha=0.2)

                # Load the spline-interpolated eigenfunctions
                spl = np.loadtxt(f'./spline/spline_{n}{mode_type}{l}.txt')
                spl_radii = spl[:,0]
                spl_f90   = spl[:,eig_id]

                # Scipy interpolation of the eigenfuncs
                cs = ci.CubicSpline(x=eig_r, y=eigen, bc_type='not-a-knot')
                eig_interp_scipy = cs(spl_radii)

                # Compute error between the interpolations:
                err = spl_f90 - eig_interp_scipy

                # Check small errors between Scipy and My spline
                assert (abs(err)  < 5e-6).all()

                # Plot the Scipy interpolated points
                ax[0, eig_id - 1].plot(eig_interp_scipy, spl_radii, 'xr', markersize=1 )

                # Plot the spline points from fortran
                ax[0, eig_id - 1].plot(spl_f90, spl_radii, 'xb', markersize=1)

                ax[0, eig_id - 1].axhline(spl_radii[np.where(spl_f90==np.max(spl_f90))[0][0]], alpha=0.2, color='k')
                ax[1, eig_id - 1].scatter(err, spl_radii)
                ax[0, eig_id - 1].set_title(hdrs[6][eig_id-1])
                ax[1, eig_id - 1].set_xlabel(' % Difference scipy - my spline')

        plt.savefig('test_spline.pdf', format='pdf')

test_spline()