#import pytest
import numpy as np


def test_rotation_matrix():
    ddir = './rotation_matrix/'
    l    = 2

    # Load the matrix for the mode solution:
    Wmode = np.loadtxt(f"{ddir}/semi_analytical{l}.txt")
    # Load the matrix for the SEM solution:
    Wsem = np.loadtxt(f"{ddir}/Wmat_{l}.txt")

    n_mode = len(Wmode)
    n_sem  = len(Wmode)

    # Check lengths are the same - comparing same l value
    assert n_mode == n_sem

    # Check length is odd
    assert (n_mode % 2) != 0

    for i in range(n_sem):

        if i == (n_sem-1)/2:
            # This is when m = 0 so should be zero for both
            # Note we cant do division for this one
            assert Wsem[i]  == 0
            assert Wmode[i] == 0
        else:
            # Error less than 0.1 % seems ok?
            assert abs(100 * (Wsem[i] - Wmode[i])/Wmode[i]) < 0.1
