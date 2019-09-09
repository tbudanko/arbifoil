"""Recursive derivative function test script"""

import numpy as np
from scipy.integrate import trapz
from matplotlib import pyplot as plt
from math import *

x = np.linspace(0, 2*np.pi, 1000, endpoint = False)
y = np.sin(x)

def derivative(f, t, a, n):
    """
    n-th order derivative of f with respect to t at a
    """
    # Concatenating to allow periodicity

    # Normalization of argument
    if a >= 2*np.pi:
        a -= 2*np.pi
    elif a < 0:
        a += 2*np.pi

    # Index finding
    i = 0
    while np.round(a, decimals = 6) > np.round(t[i], decimals = 6):
        i += 1

    if abs(t[i]-a) > abs(t[i-1]-a):
        i -= 1

    if n == 0:
        return f[i]

    else:
        if i == len(t)-1:
            derivative2 = derivative(f, t, t[0], n-1)
            derivative1 = derivative(f, t, t[i-1], n-1)
            return (derivative2 - derivative1)/(2*np.pi+t[0]-t[i-1])
        elif i == 0:
            derivative2 = derivative(f, t, t[i+1], n-1)
            derivative1 = derivative(f, t, t[i-1], n-1)
            return (derivative2 - derivative1)/(t[i+1]-t[i-1] + 2*np.pi)
        else:
            derivative2 = derivative(f, t, t[i+1], n-1)
            derivative1 = derivative(f, t, t[i-1], n-1)
            return (derivative2 - derivative1)/(t[i+1]-t[i-1])

print(derivative(y,x,np.pi/2,1))
