import numpy
from primordial.equations import Equations as _Equations

class Equations(_Equations):
    def __init__(self, K, V):
        self.K = K
        self.V = V

    def __call__(self, t, y):
        N, phi, dphi = y
        H2 = self.H2(y)
        H = numpy.sqrt(H2)
        ddphi = -3*H*dphi - self.V.d(phi)
        return [H, dphi, ddphi]

    def H2(self, y):
        N, phi, dphi = y
        H2 = (dphi**2/2 + self.V(phi))/3 - self.K*numpy.exp(-2*N)
        return H2

    def inflating(self, y):
        N, phi, dphi = y
        return self.V(phi) - dphi**2 

    def sol(self, sol):
        sol.N, sol.phi, sol.dphi = sol.y
        sol.H = numpy.sqrt(self.H2(sol.y))
        return sol


def KD_initial_conditions(t0, K, N_p, phi_p): 
    b = K * numpy.exp(-2*N_p)
    return [
            N_p + numpy.log(t0)/3 - 9/14 * b * t0**(4./3),
            phi_p - numpy.sqrt(2./3)*numpy.log(t0) - 27*numpy.sqrt(6)/56*b*t0**(4./3),
            -numpy.sqrt(2./3)/t0 - 9*numpy.sqrt(6)/14 * b * t0**(1./3),
            ]
