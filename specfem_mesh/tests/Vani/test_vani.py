import pytest
import numpy as np

def test_vani():
    ddir = "./Vani/"
    t = 'S'
    ells = [2]
    ns = [6]
    nmodes = len(ells)

    assert (len(ns) == len(ells))

    for imode in range(nmodes):

        l = ells[imode]
        n = ns[imode]

        rad = np.loadtxt(f'{ddir}/radial_{n}{t}{l}.txt')[:2 * l + 1, :]
        rad = np.diag(rad)

        sem = np.loadtxt(f'{ddir}/sem_{n}{t}{l}.txt')[:2 * l + 1, :]
        sem = np.diag(sem)

        # Compute ratio of solutions
        ratio = sem / rad
        err = abs(ratio - 1)

        assert (err < 0.001).all()


