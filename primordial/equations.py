from types import MethodType

class Equations(dict):
    """ Base class for equations.

    Allows one to compute derivatives and derived variables.
    Most of the other classes take 'equations' as an object. 

    Derived classes must define:

    __call__(self, t, y)
        The derivative function for underlying variables

    """
    def __init__(self, potential):
        self.potential = potential

    def sol(self, sol):
        return sol

    def V(self, t, y):
        return self.potential(self.phi(t, y))

    def dVdphi(self, t, y):
        return self.potential.d(self.phi(t, y))

    def variable(self, name):
        def method(self, t, y):
            return y[self[name]]
        self[name] = len(self)
        return MethodType(method, self)
