import numpy
import scipy.interpolate
from types import MethodType

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
        t, j = numpy.unique(sol.t,return_index=True)
        del sol.t
        for name, i in self.items():
            setattr(sol, name, self.interp1d(t, sol.y[i, j], **kwargs))
        tt = self.independent_variable
        setattr(sol, tt + '_events', sol.pop('t_events'))
        setattr(sol, tt, t)
        return sol

    def interp1d(self, x, y, **kwargs):
        kind = kwargs.pop('kind', 'cubic')
        bounds_error = kwargs.pop('bounds_error', False)
        return scipy.interpolate.interp1d(x, y, kind=kind, bounds_error=bounds_error, **kwargs)

    def add_independent_variable(self, name):
        """ Define the name of the independent variable """
        def method(self, t, y):
            return t
        setattr(self, name, MethodType(method, self))
        self.independent_variable = name

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

    @property
    def variables(self):
        return list(val for val in self.values())
