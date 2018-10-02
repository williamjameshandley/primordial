import numpy

class Event(object):
    """ Base class for events.

    Gives a more usable wrapper to 

    Derived classes must define:

    __call__(self, t, y)
        Scalar root function for determining event

    Attributes:
    ----------

    equations: Equations
        The equations for computing derived variables.

    direction: [-1, 0, +1]
        The direction of the root finding (if any)

    terminal: bool
        Whether to stop at this root
    """
    def __init__(self, equations, direction=0, terminal=False, value=0):
        self.equations = equations
        self.direction = direction
        self.terminal = terminal
        self.value = value


class Inflation(Event):
    """ Inflation entry/exit """
    def __call__(self, t, y):
        return self.equations.inflating(t, y) - self.value


class Stationary(Event):
    """ Tests if a is positive """
    def __call__(self, t, y):
        return self.equations.H2(t, y) - self.value

class ModeExit(Event):
    """ When mode exits the horizon """
    def __call__(self, t, y):
        return numpy.log(numpy.abs(self.equations.H2(t, y)))/2+self.equations.N(t,y) - numpy.log(self.value)

class UntilN(Event):
    """ Stop at N """
    def __call__(self, t, y):
        return self.equations.N(t, y) - self.value
