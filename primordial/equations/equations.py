import numpy
import scipy.interpolate
from types import MethodType


class Equations(object):
    """ Base class for equations.

    Allows one to compute derivatives and derived variables.
    Most of the other classes take 'equations' as an object.

    Attributes
    ----------
    i : dict
        dictionary mapping variable names to indices in the solution vector

    independent_variable : string
        name of independent variable
    """
    def __init__(self):
        self.i = {}

    def __call__(self, t, y):
        """ Vector of derivatives

        Parameters
        ----------
        t : float
            Time coordinate

        y : numpy.array
            Variable values

        Returns
        -------
        dy : numpy.array
            Vector of derivatives
        """
        raise NotImplementedError("Equations class must define __call__")

    def sol(self, sol, **kwargs):
        """ Amend solution from from solve_ivp """
        t, j = numpy.unique(sol.t, return_index=True)
        del sol.t
        for name, i in self.i.items():
            setattr(sol, name, self._interp1d(t, sol.y[i, j], **kwargs))
        tt = self.independent_variable
        setattr(sol, tt + '_events', sol.pop('t_events'))
        setattr(sol, tt, t)
        return sol

    def set_independent_variable(self, name):
        """ Set name of the independent variable

        Parameters
        ----------
        name : str
            Name of the independent variable
        """
        def method(self, t, y):
            return t

        method.__doc__ = """ Hi there """

        setattr(self, name, MethodType(method, self))
        self.independent_variable = name

    def add_variable(self, *args):
        """ Add dependent variables to the equations

        * creates an index for the location of variable in y
        * creates a class method of the same name with signature
          name(self, t, y) that should be used to extract the variable value in
          an index-independent manner.

        Parameters
        ----------
        *args : str
            Name of the dependent variables
        """
        for name in args:
            self._add_variable(name)

    def _add_variable(self, name):
        self.i[name] = len(self.i)

        def method(self, t, y):
            return numpy.array(y)[self.i[name], ...]

        method.__doc__ = """ Retrieve %s from the solution vector

        Arguments
        ---------
        t : float
            Time coordinate

        y : numpy.array
            Variable values

        Returns
        -------
        %s : float
            value of  %s
        """ % (name, name, name)

        setattr(self, name, MethodType(method, self))

    def _interp1d(self, x, y, **kwargs):
        kind = kwargs.pop('kind', 'cubic')
        bounds_error = kwargs.pop('bounds_error', False)
        return scipy.interpolate.interp1d(x, y, kind=kind,
                                          bounds_error=bounds_error, **kwargs)
