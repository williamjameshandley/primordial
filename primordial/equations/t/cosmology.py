import numpy
from primordial.equations.cosmology import Equations as _Equations

class Equations(_Equations):
    """ Cosmology equations in time 

    Solves background variables in cosmic time for curved and flat universes
    using the Friedmann equation.

    Independent variable:
        t: cosmic time
    
    Variables:
        N: efolds
    
    """
    def __init__(self, H0, Omega_r, Omega_m, Omega_k, Omega_l):
        super(Equations, self).__init__(H0, Omega_r, Omega_m, Omega_k, Omega_l)

        self.add_independent_variable('t')
        self.add_variable('N')

    def __call__(self, t, y):
        """ The derivative function for underlying variables,
            computed using the Klein-Gordon equation """
        dy = numpy.zeros_like(y)
        dy[self['N']] = self.H(t, y)
        return dy


class initial_conditions(object):
    def __init__(self, Ni):
        self.t0 = 0
        self.Ni = Ni

    def __call__(self, equations, y0):
        t0 = self.t0
        y0[equations['N']] = self.Ni
