import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad


def compute_gravitational_acceleration(radius_profile, density_profile, G=6.67430e-11):
    """
    Computes the gravitational acceleration as a function of radius for a given piecewise density profile.

    Parameters:
    radius_profile (array): Radii at which density values are given (must be sorted in ascending order).
    density_profile (array): Corresponding density values at each radius.
    G (float): Gravitational constant (default: 6.67430e-11 m^3/kg/s^2).

    Returns:
    g_r (array): Gravitational acceleration at each radius.
    """

    def mass_enclosed(r, radius_profile, density_profile):
        """Computes mass enclosed within radius r using piecewise integration."""
        mass = 0.0
        for i in range(len(radius_profile) - 1):
            r1, r2 = radius_profile[i], radius_profile[i + 1]
            rho = density_profile[i]
            if r <= r1:
                break
            elif r < r2:
                mass += quad(lambda r_p: 4 * np.pi * r_p ** 2 * rho, r1, r)[0]
                break
            else:
                mass += quad(lambda r_p: 4 * np.pi * r_p ** 2 * rho, r1, r2)[0]
        return mass

    g_r = []
    for r in radius_profile:
        if r == 0:
            g_r.append(0)  # No gravity at the exact center
        else:
            M_r = mass_enclosed(r, radius_profile, density_profile)
            g_r.append(G * M_r / r ** 2)

    return np.array(g_r)




radius_profile  = np.loadtxt("radialpoints.txt")*6371000
density_profile = np.loadtxt("rhopoints.txt")


gravity_profile = compute_gravitational_acceleration(radius_profile, density_profile)
print(gravity_profile)

np.savetxt(X=gravity_profile, fname='gravity.txt')


