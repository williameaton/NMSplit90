import scipy.special as ss
import numpy as np


def test_ylm():
    # IMPORTANT:
    # scipy uses formula with e^{i m theta} and cos(phi)
    # where as we use e^{i m phi} and cos(theta) hence the ordering
    # of phi, theta is 'backwards' in the function call relative to
    # the scipy documentation
    # Scipy also does not include the (-1)**m component of B.58 in
    # DT98 which we use in our code so add in here

    # l, m, ylm
    d = np.loadtxt("testylm.txt")

    theta = 0.5
    phi   = 0.124

    nnums = d.shape[0]
    for i in range(nnums):
        ylm = (-1)**d[i,1] * ss.sph_harm(d[i,1], d[i,0], phi, theta).real

        err = np.abs((ylm-d[i,-1])/ylm)
        assert err < 2e-6


test_ylm()

