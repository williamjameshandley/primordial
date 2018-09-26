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
    def __init__(self, equations, direction=0, terminal=False):
        self.equations = equations
        self.direction = direction
        self.terminal = terminal


class Inflation(Event):
    """ Inflation entry/exit """
    def __call__(self, t, y):
        return self.equations.inflating(t, y)


class Stationary(Event):
    """ Tests if a is positive """
    def __call__(self, t, y):
        return self.equations.H2(t, y)

