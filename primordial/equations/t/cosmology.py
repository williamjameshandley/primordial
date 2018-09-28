import numpy
from primordial.equations.equations import Equations as _Equations

class Equations(_Equations):
    """ Cosmology equations in time 

    Solves bacgkround variables in cosmic time for curved and flat universes
    using the Friedmann equation.

    Independent variable:
        t: cosmic time
    
    Variables:
        N: efolds
        phi: inflaton field
        dphi: d (phi) / dt
    
    """
    def __init__(self, H0, Omega_r, Omega_m, Omega_l):
        self.H0 = H0
        self.Omega_r = Omega_r
        self.Omega_m = Omega_m
        self.Omega_l = Omega_l
        self.Omega_k = 1 - (Omega_r + Omega_m + Omega_l)
        self.N0 = numpy.log(1/H0/numpy.sqrt(numpy.abs(self.Omega_k)))
        print(self.N0)

        self.add_variables('N')

    def __call__(self, t, y):
        """ The derivative function for underlying variables,
            computed using the Klein-Gordon equation """
        dy = numpy.zeros_like(y)
        dy[self['N']] = self.H(t, y)
        return dy

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

    def sol(self, sol):
        """ Post-process solution of solve_ivp """
        sol = super(Equations, self).sol(sol)
        sol.H = self.H(sol.t, sol.y)
        N0 = self.N0
        sol.Omega_r = self.Omega_r * numpy.exp(4*(self.N0-sol.N)) / (sol.H/self.H0)**2
        sol.Omega_m = self.Omega_m * numpy.exp(3*(self.N0-sol.N)) / (sol.H/self.H0)**2
        sol.Omega_k = self.Omega_k * numpy.exp(2*(self.N0-sol.N)) / (sol.H/self.H0)**2
        sol.Omega_l = self.Omega_l                                / (sol.H/self.H0)**2
        return sol


class initial_conditions(object):
    def __init__(self, Ni):
        self.t0 = 0
        self.Ni = Ni

    def __call__(self, equations, y0):
        t0 = self.t0
        y0[equations['N']] = self.Ni
