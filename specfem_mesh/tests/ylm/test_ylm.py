import scipy.special as ss
import numpy as np
import os


def test_ylm():
    # IMPORTANT:
    # scipy uses formula with e^{i m theta} and cos(phi)
    # where as we use e^{i m phi} and cos(theta) hence the ordering
    # of phi, theta is 'backwards' in the function call relative to
    # the scipy documentation
    # Scipy also does not include the (-1)**m component of B.58 (DT98 definition of Xlm)
    # But it does inlcude it in the definition of Plm, which DT98 does not
    # Hence, the Ylms should be the same.

    # l, m, ylm
    d = np.loadtxt("./ylm/testylm.txt")

    theta = 0.5
    phi   = 0.124

    nnums = d.shape[0]
    for i in range(nnums):
        ylm = ss.sph_harm(d[i,1], d[i,0], phi, theta).real

        err = np.abs((ylm-d[i,-1])/ylm)
        assert err < 2e-6



def test_values(): 
    # l, m, ylm
    d = np.loadtxt("./ylm/ylm_value_error.txt")

    nnums = d.shape[0]
    for i in range(nnums):

        err_r = d[i,0]
        err_i = d[i,1]

        # Real part
        assert abs(err_r) < 1e-12
        assert abs(err_i) < 1e-12