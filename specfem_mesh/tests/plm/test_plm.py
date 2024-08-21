import numpy as np
import scipy.special as ss
import os

def test_plm():

    d = np.loadtxt("./plm/testplm.txt")

    nnums = d.shape[0]
    for i in range(nnums):
        plm = ss.lpmv(d[i,1], d[i,0], np.cos(0.5))
        err = np.abs((plm-d[i,-1])/plm)

        print(err)
        assert err < 1e-5
