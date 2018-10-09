import numpy
from primordial.equations.equations import Equations as _Equations
from scipy.interpolate import interp1d

class Equations(_Equations):
    """ Cosmology equations

    Solves background variables in cosmic time for curved and flat universes
    using the Friedmann equation.

    Independent variable:
        N: efolds
    
    Variables:
        t: cosmic time
    
    """
    def __init__(self, H0, Omega_r, Omega_m, Omega_k, Omega_l):
        self.H0 = H0
        Omega = Omega_r + Omega_m + Omega_k + Omega_l
        self.Omega_m = Omega_m/Omega
        self.Omega_l = Omega_l/Omega
        self.Omega_k = Omega_k/Omega
        self.Omega_r = Omega_r/Omega
        if self.Omega_k==0:
            self.N0 = 0
        else:
            self.N0 = numpy.log(1./H0/numpy.sqrt(numpy.abs(self.Omega_k)))

    def H(self, t, y):
        """ Hubble parameter"""
        return numpy.sqrt(self.H2(t, y))

    def H2(self, t, y):
        """ The square of the Hubble parameter,
            computed using the Friedmann equation """
        N = self.N(t, y)
        return self.H0**2 * (
                  self.Omega_r * numpy.exp(4*(self.N0-N)) 
                + self.Omega_m * numpy.exp(3*(self.N0-N))
                + self.Omega_k * numpy.exp(2*(self.N0-N))
                + self.Omega_l
                )

    def sol(self, sol, **kwargs):
        """ Post-process solution of solve_ivp """
        sol.H = self.interp1d(sol.t, self.H(sol.t, sol.y), **kwargs)
        sol = super(Equations, self).sol(sol, **kwargs)
        return sol
