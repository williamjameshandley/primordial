class Equations(object):
    """ Base class for equations.

    Allows one to compute derivatives and derived variables.
    Most of the other classes take 'equations' as an object. 

    Derived classes must define:

    __call__(self, t, y)
        The derivative function for underlying variables

    """
    def sol(self, sol):
        return sol
