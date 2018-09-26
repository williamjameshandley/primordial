class Equations(object):
    """ Base class for equations.

    Allows one to compute derivatives and derived variables.
    Most of the other classes take 'equations' as an object. 

    Derived classes must define:

    __call__(self, t, y)
        The derivative function for underlying variables

    """
    _i = {}
    def __init__(self, potential):
        self.potential = potential

    def sol(self, sol):
        return sol

    def V(self, t, y):
        return self.potential(self.phi(t, y))

    def dVdphi(self, t, y):
        return self.potential.d(self.phi(t, y))

    def __getitem__(self, key):
        return self._i[key]
