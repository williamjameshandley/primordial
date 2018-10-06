import numpy
from types import MethodType
from scipy.interpolate import interp1d

class Equations(dict):
    """ Base class for equations.

    Allows one to compute derivatives and derived variables.
    Most of the other classes take 'equations' as an object. 

    Derived classes must define:

    __call__(self, t, y)
        The derivative function for underlying variables

    and add variable names using ``add_variable`` or ``add_variables``

    """

    def sol(self, sol, **kwargs):
        kind = kwargs.pop('kind', 'cubic')
        bounds_error = kwargs.pop('bounds_error', False)
        t = sol.t[:]
        for name, i in self.items():
            setattr(sol, name, interp1d(t, sol.y[i], kind=kind, bounds_error=bounds_error))
        return sol

    def add_variable(self, name):
        """ Add variable to equations.
        
        * creates an index for the location of variable in y
        * creates a class method of the same name with signature
          name(self, t, y) that should be used to extract the variable value in
          an index-independent manner.
        """
        self[name] = len(self)
        def method(self, t, y):
            return numpy.array(y)[self[name],...]
        setattr(self, name, MethodType(method, self))

    def add_variables(self, names):
        """ Add multiple variables to the equations """
        for name in names:
            self.add_variable(name)
