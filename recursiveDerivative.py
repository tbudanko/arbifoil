"""Recursive derivative function test script"""

import numpy as np
from scipy.integrate import trapz
from matplotlib import pyplot as plt
from math import *

x = np.linspace(0, 2, 100)
y = np.sin(x)

def der(xd, order):
    i = 0
    while xd > x[i]:
        i += 1
    if abs(x[i]-xd) > abs(x[i-1]-xd):
        i -= 1

    if order == 0:
        return y[i]

    else:
        der2 = der(x[i+1], order-1)
        der1 = der(x[i-1], order-1)
        return (der2 - der1)/(x[i+1]-x[i-1])

print(der(1,6))
