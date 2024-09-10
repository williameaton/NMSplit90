#import pytest
import numpy as np


def test_rotation_matrix():
    ddir = './rotation_matrix/'

    n1   = 0
    t1   = 'S'
    l1   = 2

    n2   = 0
    t2   = 'S'
    l2   = 2

    # Load the matrix for the mode solution:
    Wmode = np.loadtxt(f"{ddir}/semi_analytical_{n1}{t1}{l1}_{n2}{t2}{l2}.txt")
    # Load the matrix for the SEM solution:
    Wsem = np.loadtxt(f"{ddir}/Wmat_{n1}{t1}{l1}_{n2}{t2}{l2}.txt")

    n_mode = np.shape(Wmode)
    n_sem  = np.shape(Wsem)

    # Check lengths are the same - comparing same l value
    assert n_mode == n_sem

    # First n_sem[0]/2 rows are real, then next n_sem[0]/2 are imag
    nfreqs = int(n_sem[0]/2)

    # Check length is odd
    assert (nfreqs % 2) != 0

    half_id = int((nfreqs - 1)/2)

    for row in range(int(n_sem[0])):
        for col in range(n_sem[1]):
            if row != col  or row >= nfreqs or (row == half_id and col == half_id) :
                # This is when m = 0 so should be zero for both
                # Note we cant do division for this one
                assert Wsem[row, col]  == 0
                assert Wmode[row, col] == 0
            else:
                # % Error less than 0.1 % seems ok?
                diff_perc = abs(100 * (Wsem[row, col] - Wmode[row, col])/Wmode[row, col])
                assert diff_perc < 0.1

