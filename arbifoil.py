"""
ARBIFOIL
"""

import numpy as np
from scipy.integrate import trapz
from matplotlib import pyplot as plt
from math import *

class foil():
    """

    """
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
            self.z1[i] = 4.02*float(coordinates[0]) - 2.02 + 4.01*float(coordinates[1])*1j


        """ z2 - pseudocircle complex plane """
        self.z2, self.psi, self.theta =  self.inverseJoukowsky(self.z1)

        """Theodorsen mapping"""
        self.phi, self.eta = self.theodorsen(1)



    def inverseJoukowsky(self, z1):
        """
        Inverse Joukowsky mapping
        Two inverses exist, one is chosen arbitrarily and adhered to by smoothly completing the pseudocircle.
        """

        z2 = np.empty(len(z1), dtype = complex)

        z2_1 = z1[0]/2 + ((z1[0]/2)**2 - 1)**0.5
        z2_2 = z1[0]/2 - ((z1[0]/2)**2 - 1)**0.5

        if abs(z2_1) >= 1:
            z2[0] = z2_1
        else:
            z2[0] = z2_2

        for i in range(1, len(z2)):
            z2_1 = z1[i]/2 + ((z1[i]/2)**2 - 1)**0.5
            z2_2 = z1[i]/2 - ((z1[i]/2)**2 - 1)**0.5

            del1 = abs(z2_1 - z2[i-1])
            del2 = abs(z2_2 - z2[i-1])

            if del1 < del2:
                z2[i] = z2_1
            else:
                z2[i] = z2_2

        # Psi and theta arrays are constructed ...
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


    def theodorsen(self, rtol):
        """
        Theodorsen mapping
        Takes pseudocircle and returns:
        - eta = eta(phi)
        - psi = psi(phi)
        within specified numerical precision (residual tolerance).
        """
        nPoints = 100
        phi = np.linspace(self.theta[0], self.theta[-1], nPoints)

        eta = np.empty(nPoints)
        for i in range(nPoints):

            # Subintegral function
            integral = np.empty(nPoints)
            for j in range(nPoints):
                integral[j] = self.fPsi(phi[j]) / tan(phi[j] - phi[i])

            eta[i] = -1/(2*np.pi) * trapz(integral, phi)

        return phi, eta


    def fPsi(self, arg):
        """
        Interpolating function
        psi = psi(complex argument)
        """

        i = 0
        while arg > self.theta[i]:
            i += 1

        f = (arg - self.theta[i-1])/(self.theta[i] - self.theta[i-1])

        return (1 - f)*self.psi[i-1] + f*self.psi[i]

    def dPsi(self, arg, order):
        """
        Derivative of psi with respect to argument
        dpsi^(order) / darg
        """
        i = 0
        while arg > self.theta[i]:
            i += 1
        if abs(self.theta[i]-arg) > abs(self.theta[i-1]-arg):
            i -= 1

        if order == 0:
            return self.psi[i]

        else:
            dPsi2 = self.dPsi(self.theta[i+1], order-1)
            dPsi1 = self.dPsi(self.theta[i-1], order-1)
            return (dPsi2 - dPsi1)/(self.theta[i+1]-self.theta[i-1])



#           for i in range(nPoints):
#            integralFunc = np.empty(nPoints)
#            for j in range(nPoints):
#                integralFunc[i] = self.psi[i] * cot(0.5*())
#
#            self.eta[i] = -1/(2*pi) * trapz()



    def plot(self):
        """
        Plotting method.
        """

        plt.figure()
        plt.plot(self.z1.real, self.z1.imag, 'b')
        plt.plot(self.z2.real, self.z2.imag, 'r')
        plt.gca().set_aspect('equal')
        plt.grid('True')
        plt.show()

test = foil('naca0012.txt')
test.plot()
