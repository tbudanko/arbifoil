"""
ARBIFOIL
"""

import numpy as np
from scipy.integrate import trapz
from matplotlib import pyplot as plt
from math import *

class foil():

    def __init__(self, datFileName):
        """
        Foil initialization
        - .dat file read into z1 complex array
        - inverse Joukowsky z1 -> z2
        - theodorsen mapping
        """
        # Read .dat file
        datFile = open(datFileName, 'r')
        datPoints = datFile.read()
        datFile.close()

        datPoints = datPoints.split('\n')
        datPoints.pop(0) # remove header line
        datPoints.pop(-1) # remove whitespace
        datPoints.pop(-1) # remove repeated trailing edge

        """ z1 - original airfoil complex plane """
        self.z1 =  np.empty(len(datPoints), dtype = complex)
        for i in range(len(self.z1)):
            coordinates = datPoints[i].split()
            self.z1[i] = float(coordinates[0]) + float(coordinates[1])*1j

        LEindex = 1
        while self.z1[LEindex].imag > 0:
            LEindex += 1

        LEradius = getRadius(self.z1[LEindex-1].real, self.z1[LEindex-1].imag, \
                                  self.z1[LEindex].real, self.z1[LEindex].imag, \
                                  self.z1[LEindex+1].real, self.z1[LEindex+1].imag)

        self.z1 = 4/(1-0.5*LEradius)*self.z1
        self.z1 = self.z1 - self.z1[0].real + 2


        """ z2 - pseudocircle complex plane """
        self.z2, self.psi, self.theta =  self.inverseJoukowsky(self.z1)

        """Theodorsen mapping"""
        self.phi, self.eta = self.theodorsen(0.001)

    def inverseJoukowsky(self, z1):
        """
        Inverse Joukowsky mapping
        Two inverses exist, one is chosen arbitrarily and adhered to
        by smoothly completing the pseudocircle.
        """

        z2 = np.empty(len(z1), dtype = complex)
        z2[0] = 1 + 0j # Trailing edge

        z2_1 = z1[1]/2 + ((z1[1]/2)**2 - 1)**0.5
        z2_2 = z1[1]/2 - ((z1[1]/2)**2 - 1)**0.5

        if z2_1.imag > 0:
            z2[1] = z2_1
        else:
            z2[1] = z2_2

        for i in range(2, len(z2)):
            z2_1 = z1[i]/2 + ((z1[i]/2)**2 - 1)**0.5
            z2_2 = z1[i]/2 - ((z1[i]/2)**2 - 1)**0.5

            v3 = z2[i-1] - z2[i-2]
            v1 = z2_1 - z2[i-1]
            v2 = z2_2 - z2[i-1]

            cos1 = (v1.real*v3.real + v1.imag*v3.imag)/np.abs(v1)/np.abs(v3)
            cos2 = (v2.real*v3.real + v2.imag*v3.imag)/np.abs(v2)/np.abs(v3)

            if cos1 > cos2:
                z2[i] = z2_1
            else:
                z2[i] = z2_2

        # Psi and theta arrays are constructed
        # so that the value of psi can be interpolated.
        psi = np.empty(len(z2))
        theta = np.empty(len(z2))
        for i in range(len(z2)):
            psi[i] = np.log(np.abs(z2[i]))
            if np.angle(z2[i]) < 0:
                theta[i] = np.angle(z2[i]) + 2*np.pi
            else:
                theta[i] = np.angle(z2[i])

        return z2, psi, theta

    def theodorsen(self, resTol):
        """
        Theodorsen mapping
        Takes pseudocircle and returns:
        - eta = eta(phi)
        - psi = psi(phi)
        within specified numerical precision (residual tolerance).
        """
        # Circle discretization
        nPoints = 100
        phi = np.linspace(self.theta[0], self.theta[-1], nPoints)
        delPhi = phi[1] - phi[0] # Step size

        k = 0
        res = 1

        etaOld = np.zeros(nPoints)
        etaNew = np.zeros(nPoints)

        while res > resTol:
            # Evaluation of integral (VIII) from reference 1.
            # using the method presented in the appendix.
            for i in range(nPoints):
                for j in range(nPoints):
                    if i==j:
                        etaNew[i] += 2 * delPhi * derivative(self.psi, self.theta, phi[i]-etaOld[i], 1) * (1 - derivative(etaOld, phi, phi[i], 1))
                        #etaNew[i] += 2 * delPhi * self.dPsi(phi[i]-etaOld[i], 1)
                    else:
                        etaNew[i] += 2 * self.fPsi(phi[j]-etaOld[j]) * \
                                np.log(np.sin(0.5*(phi[j] + 0.5*delPhi - phi[i]))/np.sin(0.5*(phi[j] - 0.5*delPhi - phi[i])))
                etaNew[i] *= -1/(2*np.pi)

            res = np.amax(etaNew - etaOld)
            etaOld = etaNew
        return phi, etaNew

    def fPsi(self, arg):
        """
        Interpolating function
        psi = psi(argument)
        """
        # Expand arrays to full circle
        theta = np.append(self.theta, self.theta[0] + 2*np.pi)
        psi = np.append(self.psi, self.psi[0])

        # Normalization of argument
        if arg >= 2*np.pi:
            arg -= 2*np.pi
        elif arg < 0:
            arg += 2*np.pi

        # Index finding
        i = 0
        while np.round(arg, decimals = 6) > np.round(theta[i], decimals = 6):
            i += 1

        f = (arg - theta[i-1])/(theta[i] - theta[i-1])

        return (1 - f)*psi[i-1] + f*psi[i]



    def dPsi(self, arg, order):
        """
        n-th order derivative of psi with respect to argument
        dpsi^(order) / darg
        """

        # Concatenating to allow periodicity
        theta = np.concatenate((self.theta, self.theta + 2*np.pi))
        psi = np.concatenate((self.psi, self.psi))

        # Normalization of argument
        if arg >= 2*np.pi:
            arg -= 2*np.pi
        elif arg < 0:
            arg += 2*np.pi

        # Index finding
        i = 0
        while np.round(arg, decimals = 6) > np.round(theta[i], decimals = 6):
            i += 1

        if abs(theta[i]-arg) > abs(theta[i-1]-arg):
            i -= 1

        if order == 0:
            return psi[i]

        else:
            dPsi2 = self.dPsi(theta[i+1], order-1)
            dPsi1 = self.dPsi(theta[i-1], order-1)
            return (dPsi2 - dPsi1)/(theta[i+1]-theta[i-1])

    def plot(self):
        """
        Plotting method.
        """

        plt.figure()
        plt.plot(self.z1.real, self.z1.imag, 'b')
        plt.plot(self.z2.real, self.z2.imag, 'r')
        plt.gca().set_aspect('equal')
        plt.grid('True')
        plt.figure()
        plt.subplot(211)
        plt.plot(self.theta, self.psi, 'k')
        plt.grid('True')
        plt.subplot(212)
        plt.plot(self.phi-self.eta, self.eta, 'r')
        plt.grid('True')
        plt.show()

def getRadius(x1, y1, x2, y2, x3, y3):
    """
    Auxiliary function :
    Returns local radius of curvature defined by 3 points .
    """

    A = x1*(y2-y3) - y1*(x2-x3) + x2*y3 - x3*y2
    B = (x1**2+y1**2)*(y3-y2) + (x2**2+y2**2)*(y1-y3) + (x3**2+y3**2)*(y2-y1)
    C = (x1**2+y1**2)*(x2-x3) + (x2**2+y2**2)*(x3-x1) + (x3**2+y3**2)*(x1-x2)
    D = (x1**2+y1**2)*(x3*y2-x2*y3) + (x2**2+y2**2)*(x1*y3-x3*y1) + (x3**2+y3**2)*(x2*y1-x1*y2)

    return sqrt((B**2+C**2-4*A*D)/(4*A**2))

def derivative(f, t, a, n):
    """
    n-th order derivative of f with respect to t at a
    """

    # Concatenating to allow periodicity
    t = np.concatenate((t, t + 2*np.pi))
    f = np.concatenate((f, f))

    # Normalization of aument
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
        derivative2 = derivative(f, t, t[i+1], n-1)
        derivative1 = derivative(f, t, t[i-1], n-1)
        return (derivative2 - derivative1)/(t[i+1]-t[i-1])

test = foil('naca6409.txt')
test.plot()
