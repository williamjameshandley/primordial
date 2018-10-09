import numpy
from primordial.equations.cosmology import Equations as _Equations

class Equations(_Equations):
    """ Cosmology equations in time 

    Solves background variables in cosmic time for curved and flat universes
    using the Friedmann equation.

    Independent variable:
        N: efolds
    
    Variables:
        t: cosmic time
    
    """
    def __init__(self, H0, Omega_r, Omega_m, Omega_l):
        super(Equations, self).__init__(H0, Omega_r, Omega_m, Omega_l)

        self.add_independent_variable('N')
        self.add_variable('t')

    def __call__(self, N, y):
        """ The derivative function for underlying variables,
            computed using the Klein-Gordon equation """
        dy = numpy.zeros_like(y)
        dy[self['t']] = 1/self.H(N, y)
        return dy



    def sol(self, sol):
        """ Post-process solution of solve_ivp """
        N = sol.t
        sol = super(Equations, self).sol(sol)
        sol.N = N
        return sol


class initial_conditions(object):
    def __init__(self, Ni):
        self.t0 = Ni

    def __call__(self, equations, y0):
        y0[equations['t']] = 0