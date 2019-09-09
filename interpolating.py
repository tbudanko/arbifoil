import numpy as np
from scipy.integrate import trapz
from matplotlib import pyplot as plt
from math import *


def interpolate(f, t, a):
    """
    Interpolating function
    f = f(t) in t = a
    """
    # Normalization of aument
    if a >= 2*np.pi:
        a -= 2*np.pi
    elif a < 0:
        a += 2*np.pi

    # Index finding
    if a > t[-1] or a < t[0]:
        frac = (a - t[-1])/(t[0] + 2*np.pi - t[-1])
        return (1 - frac)*f[-1] + frac*f[0]
    else:
        i = 0
        while np.round(a, decimals = 6) > np.round(t[i], decimals = 6):
            i += 1
            if i == len(t):
                break

        frac = (a - t[i-1])/(t[i] - t[i-1])
        return (1 - frac)*f[i-1] + frac*f[i]

x = np.linspace(0, 2*np.pi, 100, endpoint = False)
y = np.sin(x)

print(interpolate(y,x,0-0.00001))
