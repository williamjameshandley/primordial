""" Base classes for inflationary solvers"""
import numpy
from primordial.equations.equations import Equations as _Equations

class Equations(_Equations):

    def __init__(self, K, potential):
        super(Equations, self).__init__()
        self.potential = potential
        self.K = K

    def H(self, t, y):
        """ Hubble parameter"""
        return numpy.sqrt(self.H2(t, y))

    def V(self, t, y):
        """ Potential """
        return self.potential(self.phi(t, y))

    def dVdphi(self, t, y):
        """ Potential derivative """
        return self.potential.d(self.phi(t, y))

    def sol(self, sol, **kwargs):
        """ Post-process solution of solve_ivp """
        sol.H = self._interp1d(sol.t, self.H(sol.t, sol.y), **kwargs)
        sol = super(Equations, self).sol(sol, **kwargs)
        return sol
