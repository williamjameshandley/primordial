import numpy
from primordial.equations.equations import Equations as _Equations

class Equations(_Equations):
    """ Background equations in time 

    Solves bacgkround variables in cosmic time for curved and flat universes
    using the Klein-Gordon and Friedmann equations.

    Independent variable:
        N: e-folds (log a)
    
    Variables:
        phi: inflaton field
        dphi: d/dN (phi) 
    
    """
    def __init__(self, K, potential):
        self.potential = potential
        self.K = K
        self.add_variables(['phi', 'dphi'])

    def __call__(self, N, y):
        """ The derivative function for underlying variables,
            computed using the Klein-Gordon equation """
        H2 = self.H2(N, y)
        dphi = self.dphi(N, y)
        dVdphi = self.dVdphi(N, y)
        dlogH = self.dlogH(N, y)
        ddphi = -(dlogH + 3)* dphi - dVdphi/H2

        dy = numpy.zeros_like(y)
        dy[self['phi']] = dphi
        dy[self['dphi']] = ddphi

        return dy

    def H(self, N, y):
        """ Hubble parameter"""
        return numpy.sqrt(self.H2(N, y))

    def H2(self, N, y):
        """ The square of the Hubble parameter,
            computed using the Friedmann equation """
        V = self.V(N, y)
        dphi = self.dphi(N, y) 

        return (2*V - 6*self.K*numpy.exp(-2*N)) / (6 - dphi**2)

    def dlogH(self, N, y):
        """ d/dN log H """
        dphi = self.dphi(N, y) 
        H2 = self.H2(N, y) 
        return - dphi**2/2 + self.K * numpy.exp(-2*N)/H2

    def N(self, N, y):
        return N

    def V(self, N, y):
        """ Potential """
        return self.potential(self.phi(N, y))

    def dVdphi(self, N, y):
        """ Potential derivative """
        return self.potential.d(self.phi(N, y))

    def inflating(self, N, y):
        """ Inflation diagnostic """
        return self.V(N, y) / (self.V(N, y)/2 - self.K*numpy.exp(-2*N)) - self.dphi(N, y)**2 

    def sol(self, sol):
        """ Post-process solution of solve_ivp """
        sol = super(Equations, self).sol(sol)
        sol.N = sol.t
        sol.H = self.H(sol.N, sol.y)
        return sol

class Inflation_start_initial_conditions(object):
    def __init__(self, N_e, phi_e):
        self.t0 = N_e
        self.phi_e = phi_e

    def __call__(self, equations, y0):
        V = equations.potential(self.phi_e)
        N = self.t0
        y0[equations['phi']] = self.phi_e
        y0[equations['dphi']] = -numpy.sqrt(V / (V/2 - equations.K*numpy.exp(-2*N)))
