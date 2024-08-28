import numpy as np 

def test_xlm():
    ddir = './xlm/'
    d = np.loadtxt(f"{ddir}/xlm_integral.txt")

    for i in range(d.shape[0]):
        if int(d[i,0]) == int(d[i,1]):
            assert abs(d[i,2] - 1/(2*np.pi) )  < 1e-12
        else:
            assert abs(d[i,2])  < 1e-12


