import numpy as np 


ddir = "./integration/"

def test_ylm_integration(): 

    # When m = 2 should be non-zero result for real part
    with open(f"{ddir}/test_ylm_int_2.txt") as file:
        integral = file.readline().split()
        real = float(integral[0])
        imag = float(integral[1])
        analytical = float(file.readline().split()[0])

    assert abs(real - analytical) < 2e-8
    assert abs(imag) < 1e-16


    # When m = 3 should all be zero
    with open(f"{ddir}/test_ylm_int_3.txt") as file:
        integral = file.readline().split()
        real = float(integral[0])
        imag = float(integral[1])
        analytical = float(file.readline().split()[0])

    assert abs(real) < 1e-16
    assert abs(imag) < 1e-16
    assert abs(analytical) < 1e-16



def test_1_integration(): 
    # Integral of unity over the inner core
    d = np.loadtxt(f"{ddir}/test_1_int.txt")
    assert abs(d[0] - d[1]) < 5e-9



def test_r_integration(): 
    # Integral of unity over the inner core
    d = np.loadtxt(f"{ddir}/test_r_int.txt")
    assert abs(d[0] - d[1]) < 1e-9


