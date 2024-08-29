import numpy as np
import scipy.special as ss
import os

def test_plm():

    d = np.loadtxt("./plm/testplm.txt")

    nnums = d.shape[0]
    for i in range(nnums):
        # NOTE THAT OUR PLMS do not include the -1**m that the scipy does
        plm = ss.lpmv(d[i,1], d[i,0], np.cos(0.5)) / ((-1)**d[i,1])
        err = np.abs((plm-d[i,-1])/plm)

        print(err)
        assert err < 1e-5



def test_values(): 

    d = np.loadtxt("./plm/plm_value_error.txt")

    nnums = len(d)
    for i in range(nnums):
        assert abs(d[i]) < 1e-12