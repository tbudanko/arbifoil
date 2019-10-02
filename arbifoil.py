"""
      ARBIFOIL
Arbitrary airfoil analysis using
Theodorsen's conformal mapping method

Author: Toma Budanko
"""

import numpy as np
from scipy.integrate import trapz
from matplotlib import pyplot as plt
from math import *

class foil():
    """
    Foil computational object.

    References:
    [1] Theodorsen, T.: Theory of Wing Sections of Arbitrary Shape, NACA, 1931.
    [2] Theodorsen, T., Garrick I.E.: General Potential Theory of Arbitrary
        Wing Sections, NACA TR-452, 1934.
    [3] Karamcheti, K.: Principles of ideal-fluid aerodynamics,
        R. E. Krieger Publishing Company, 1980.
    """
    def __init__(self, datFileName):
        """
        Foil initialization function:
        - .dat file read into z1 complex array
        - inverse Joukowsky z1 -> z2
        - Theodorsen mapping iteratively determined
        """
        # Read .dat file
        datFile = open(datFileName, 'r')
        datPoints = datFile.read()
        datFile.close()

        datPoints = datPoints.split('\n')
        datPoints.pop(0) # Remove header
        datPoints = list(filter(None, datPoints))
        datPoints.pop(-1) # Remove repeated trailing edge


        """ z1 - original airfoil complex plane """
        self.z1 =  np.empty(len(datPoints), dtype = complex)
        for i in range(len(self.z1)):
            coordinates = datPoints[i].split()
            self.z1[i] = float(coordinates[0]) + float(coordinates[1])*1j


        LEindex = int(np.where(self.z1.real == np.amin(self.z1.real))[0])

        LEradius = getRadius(self.z1[LEindex-1].real, self.z1[LEindex-1].imag, \
                                  self.z1[LEindex].real, self.z1[LEindex].imag, \
                                  self.z1[LEindex+1].real, self.z1[LEindex+1].imag)


        self.LEindex = LEindex

        self.z1 = 4/(1-0.5*LEradius)*self.z1
        self.z1 = self.z1 - self.z1[0].real + 2
        #self.z1[0] = 2 + 0*1j


        # Chord length in the z1 plane
        self.chord = 2 - self.z1[self.LEindex].real

        """ Inverse Joukowsky mapping """
        z2 = np.empty(len(self.z1), dtype = complex)
        #z2[0] = 1 + 0*1j # Trailing edge
        z2[0] = self.z1[0]/2 + ((self.z1[0]/2)**2 - 1)**0.5

        z2_1 = self.z1[1]/2 + ((self.z1[1]/2)**2 - 1)**0.5
        z2_2 = self.z1[1]/2 - ((self.z1[1]/2)**2 - 1)**0.5

        if z2_1.imag > 0:
            z2[1] = z2_1
        else:
            z2[1] = z2_2

        for i in range(2, len(z2)):
            z2_1 = self.z1[i]/2 + ((self.z1[i]/2)**2 - 1)**0.5
            z2_2 = self.z1[i]/2 - ((self.z1[i]/2)**2 - 1)**0.5

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

        self.z2 = z2
        self.psi = psi
        self.theta = theta

        """ Theodorsen mapping """

        # Discretization
        self.nPoints = 50
        phi = np.linspace(0, 2*np.pi, self.nPoints, endpoint = False)
        delPhi = phi[1] - phi[0] # Step size

        # Residual tolerance
        resTol = 0.001
        res = 1 # Initial residual

        print('Theodorsen mapping\nAirfoil: {}\nDiscretization: {} points\nResidual tolerance: {} rad\n'.format(datFileName, self.nPoints, resTol))

        # Evaluation of integral (VIII) from reference [1].
        # using the method presented in the appendix.
        etaNew = np.zeros(self.nPoints)
        etaOld = np.zeros(self.nPoints)

        print('Iteration\tResidual')

        iter=1
        while res > resTol:
            for i in range(self.nPoints):
                etaCur = 0
                for j in range(self.nPoints):
                    if i==j:
                        etaCur += 2 * delPhi * derivative(self.psi, self.theta, phi[i]-etaOld[i], 1) *\
                                                            (1 - derivative(etaOld, phi, phi[i], 1))
                    else:
                        etaCur += 2 * interpolate(self.psi, self.theta, phi[j]-etaOld[j]) * \
                                np.log(np.sin(0.5*(phi[j] + 0.5*delPhi - phi[i]))/np.sin(0.5*(phi[j] - 0.5*delPhi - phi[i])))

                etaCur *= -1/(2*np.pi)
                etaNew[i] = etaCur

            res = np.amax(np.abs(etaNew - etaOld))

            for i in range(self.nPoints):
                etaOld[i] = etaNew[i]

            print('{}\t{}'.format(iter, res))
            iter += 1


        self.phi = phi
        self.eta = etaNew
        print('Calculation done.')

        # Psi, theta, z2 and z1 arrays are recomputed to correspond
        # element-wise to the phi and eta arrays.
        self.psi = np.array([interpolate(self.psi, self.theta, self.phi[i] - self.eta[i]) for i in range(self.nPoints)])
        self.theta = self.phi - self.eta
        self.z2 = np.exp(self.psi)*(np.cos(self.theta)+1j*np.sin(self.theta))
        self.z1 = self.z2 + 1/self.z2

        # Phi and Psi arrays with repeated first values to close the unit circle
        # for integration using trapezoidal rule.
        phi_temp = np.append(self.phi, self.phi[0]+2*np.pi)
        psi_temp = np.append(self.psi, self.psi[0])

        # Average exponential scaling factor
        self.psi0 = 1/2/np.pi*trapz(psi_temp, phi_temp)

        # z3 array (circle plane)
        self.z3 = exp(self.psi0)*(np.cos(self.phi) + 1j*np.sin(self.phi))

        # Mapping coefficients of the combined Theodorsen-Joukowsky map
        self.nCoeffs = 100 # Number of coefficients a
        self.c = np.empty(self.nCoeffs + 1, dtype = complex)
        self.a = np.empty(self.nCoeffs, dtype = complex)
        k      = np.zeros(self.nCoeffs + 2, dtype = complex)
        k[0]   = 1 # k_0 = 1
        h      = np.zeros(self.nCoeffs + 2, dtype = complex)
        h[0]   = 1 # k_0 = 1

        for n in range(self.nCoeffs+1):
            A = np.exp(self.psi0)**(n+1)/np.pi*trapz(np.array([psi_temp[i] * np.cos((n+1)*phi_temp[i]) for i in range(len(phi_temp))]), phi_temp)
            B = np.exp(self.psi0)**(n+1)/np.pi*trapz(np.array([psi_temp[i] * np.sin((n+1)*phi_temp[i]) for i in range(len(phi_temp))]), phi_temp)
            self.c[n] = A + B*1j

            k[n+1] = sum([k[n-o]*self.c[o]*(o+1)/(n+1) for o in range(n+1)])
            h[n+1] = sum([-h[n-o]*self.c[o]*(o+1)/(n+1) for o in range(n+1)])

        for n in range(self.nCoeffs):
            self.a[n] = k[n+2] + h[n]

        # eta_t - angle of attack at zero lift -> phi(theta=0)
        self.eta_t = interpolate(self.eta, self.theta, 0) # + theta_t = 0

        # Surface derivative of z1 wrt. z3
        # Used for surface velocity and pressure coefficient calculation.
        self.dz1dz3_1 = np.empty(self.nPoints)
        for i in range(self.nPoints):
            self.dz1dz3_1[i] = abs((self.z2[i]**2-1)/(self.z2[i]*self.z3[i])*\
                                (1-1j*derivative(self.psi, self.theta, self.theta[i], 1))/\
                                (1+derivative(self.eta, self.theta, self.theta[i], 1)))

        #self.dz1dz3_2 = np.empty(self.nPoints)
        #for i in range(self.nPoints):
        #    self.dz1dz3_2[i] = sqrt((self.z1[i].imag**2*(1/tan(self.theta[i]))**2+\
        #                        self.z1[i].real**2*tan(self.theta[i])**2)*(1+derivative(self.phi,self.theta,self.theta[i],1)**2))/\
        #                        exp(self.psi0)/(1+derivative(self.eta,self.theta,self.theta[i],1))

        #self.dz1dz3_3 = np.ones(self.nPoints)
        #for n in range(self.nCoeffs):
        #    self.dz1dz3_3 += np.abs(-(n+1)*self.a[n]/self.z3**(n+2))

    def testMap(self):
        """
        Visually test the combined mapping.
        """
        z3_test = np.exp(self.psi0)*(np.cos(self.phi)+1j*np.sin(self.phi))
        z1_test = self.c[0] + z3_test
        for n in range(self.nCoeffs):
            z1_test +=  self.a[n]/(z3_test)**(n+1)

        plt.figure()
        plt.plot(np.append(z1_test.real,z1_test[0].real), np.append(z1_test.imag,z1_test[0].imag), 'b')
        plt.plot(np.append(self.z1.real,self.z1[0].real), np.append(self.z1.imag,self.z1[0].imag), 'r')
        plt.plot(self.z3.real, self.z3.imag)
        plt.gca().set_aspect('equal')
        plt.grid('True')
        plt.show()

        plt.figure()
        plt.plot(self.a)
        plt.show()

    def C_L(self, aoa):
        """
        Lift coefficient at specified AoA
        Reference: [2], p.194, expression (45)
        """
        aoa = aoa/180*np.pi

        Gamma = 4*np.pi*sin(aoa-self.eta_t)

        return 2*Gamma*np.exp(self.psi0)/self.chord

    def CoP(self, aoa):
        """
        Nondimensional location of center of pressure at specified AoA
        Reference: [2], p.194, expression (48)
        """
        aoa = aoa/180*np.pi

        m = np.abs(self.c[0]) # Modulus of c1
        delta = np.angle(self.c[0]) # Argument of c1

        b_squared = np.abs(self.a[0]) # Modulus of a1
        gamma_a1 = np.angle(self.a[0]) # Half argument of a1

        hm = -b_squared*sin(2*(-aoa+gamma_a1))/(2*exp(self.psi0)*sin(aoa-self.eta_t))

        return 0.5-(m*cos(delta-aoa) + hm)/cos(aoa)/self.chord

    def C_p(self, aoa):
        """
        Pressure coefficient along top and bottom airfoil surfaces.
        Reference: [2], p.192, expression(40)
        """
        aoa = aoa/180*np.pi

        cp = np.empty(self.nPoints)
        v = 0
        for i in range(self.nPoints):
            v = abs(cos(aoa)-1j*sin(aoa)-cos(aoa-2*self.phi[i])+1j*sin(aoa-2*self.phi[i])+2*1j*sin(aoa-self.eta_t)*(cos(self.phi[i])-1j*sin(self.phi[i])))/self.dz1dz3_2[i]
            cp[i] = 1 - v**2

        # x coordinate from LE towards TE along Chord
        # x(LE) = 0; x(TE) = 1
        x = np.real(self.z1[:])
        xmin = np.amin(x)
        LEindex = int(np.where(x == xmin)[0])
        #x += abs(xmin)
        #x /= self.chord

        # Split into upper and lower airfoil surfaces
        xUpper = np.flip(x[0:LEindex])
        xLower = x[LEindex:-1]

        cpUpper = np.flip(cp[0:LEindex])
        cpLower = cp[LEindex:-1]

        # Plot
        plt.figure()
        #plt.plot(xUpper, cpUpper, 'b')
        #plt.plot(xLower, cpLower, 'r')
        plt.plot(np.append(x,x[0]), np.append(cp,cp[0]), 'b')
        plt.gca().invert_yaxis()
        plt.grid('True')
        plt.show()

        return xUpper, cpUpper, xLower, cpLower

    def C_M_LE(self, aoa):
        """
        Moment coefficient about leading edge at specified AoA.
        Convention is respected with respect to sign, where a positive Moment
        induces a pitch up.

        Reference: [2], p.194, expression (46)
        """
        aoa = aoa/180*np.pi

        m = np.abs(self.c[0]) # Modulus of c1
        delta = np.angle(self.c[0]) # Argument of c1

        b_squared = np.abs(self.a[0]) # Modulus of a1
        gamma_a1 = np.angle(self.a[0]) # Half argument of a1

        return 4*np.pi*b_squared*sin(2*(aoa + gamma_a1))

    def AC(self):
        """
        Coordinates of the airfoil aerodynamic center.

        Reference: [2], p.195
        """
        m = np.abs(self.c[0]) # Modulus of c1
        delta = np.angle(self.c[0]) # Argument of c1

        b_squared = np.abs(self.a[0]) # Modulus of a1
        gamma_a1 = np.angle(self.a[0]) # Half argument of a1

        x = b_squared/exp(self.psi0)*cos(2*gamma_a1 - self.eta_t)
        y = b_squared/exp(self.psi0)*sin(2*gamma_a1 - self.eta_t)

        return x/self.chord, y/self.chord

    def C_M_AC(self):
        """
        Moment coefficient about aerodynamic center.
        Pitch up moment is positive.

        Reference: [2], p.195, expression (51)
        """
        b_squared = np.abs(self.a[0]) # Modulus of a1
        gamma_a1 = np.angle(self.a[0]) # Half argument of a1

        return -4*np.pi*b_squared/(self.chord**2)*sin(2*(gamma_a1-self.eta_t))

    def plot(self):
        """
        Plotting method.
        """

        plt.figure()
        plt.plot(np.append(self.z1.real,self.z1[0].real), np.append(self.z1.imag,self.z1[0].imag), 'b')
        plt.gca().set_aspect('equal')
        plt.grid('True')
        plt.figure()
        plt.subplot(211)
        plt.plot(self.theta, self.psi, 'k')
        plt.grid('True')
        plt.subplot(212)
        plt.plot(self.theta, self.eta, 'r')
        plt.grid('True')
        plt.show()

def getRadius(x1, y1, x2, y2, x3, y3):
    """
    Auxiliary function :
    Returns local radius of curvature defined by 3 points.
    Used for calculating leading edge radius.
    """

    A = x1*(y2-y3) - y1*(x2-x3) + x2*y3 - x3*y2
    B = (x1**2+y1**2)*(y3-y2) + (x2**2+y2**2)*(y1-y3) + (x3**2+y3**2)*(y2-y1)
    C = (x1**2+y1**2)*(x2-x3) + (x2**2+y2**2)*(x3-x1) + (x3**2+y3**2)*(x1-x2)
    D = (x1**2+y1**2)*(x3*y2-x2*y3) + (x2**2+y2**2)*(x1*y3-x3*y1) + (x3**2+y3**2)*(x2*y1-x1*y2)

    if A != 0:
        return sqrt((B**2+C**2-4*A*D)/(4*A**2))
    else:
        return 0

def interpolate(f, t, a):
    """
    Interpolating function
    f = f(t) in t = a where:
        t - angle array covering the unit circle without repeated first value
        f - array of function values corresponding to t
    """
    # Normalization of argument a
    a = redArg(a)

    if a < t[0]:
        a += 2*np.pi

    if a > t[-1]:
        frac = (a - t[-1])/(t[0] + 2*np.pi - t[-1])
        return (1 - frac)*f[-1] + frac*(f[0])

    else:
        i = 0
        while np.round(a, decimals = 6) > np.round(t[i], decimals = 6):
            i += 1
            if i == len(t):
                break
        frac = (a - t[i-1])/(t[i] - t[i-1])
        return (1 - frac)*f[i-1] + frac*f[i]


def derivative(f, t, a, n):
    """
    n-th order derivative of f with respect to t at a
    """
    # Normalization of argument a
    a = redArg(a)

    if a < t[0]:
        a += 2*np.pi

    if a > t[-1]:
        i = 0
        if abs(t[0]+2*np.pi-a) > abs(t[-1]-a):
            i = -1

    else:
        i = 0
        while np.round(a, decimals = 6) > np.round(t[i], decimals = 6):
            i += 1
            if i == len(t):
                break

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

def redArg(a):
    """
    Reduces argument to [0,2pi] interval.
    """
    if a >= 2*np.pi:
        a -= 2*np.pi
    elif a < 0:
        a += 2*np.pi
    return a


test = foil("clarky.txt")
test.AC()
