import numpy


class Event(object):
    """ Base class for events.

    Gives a more usable wrapper to callable event to be passed to
    `scipy.integrate.solve_ivp`

    Parameters
    ----------
    equations: Equations
        The equations for computing derived variables.

    direction: int [-1, 0, +1], optional, default 0
        The direction of the root finding (if any)

    terminal: bool, optional, default False
        Whether to stop at this root

    value: float, optional, default 0
        Offset to root
    """
    def __init__(self, equations, direction=0, terminal=False, value=0):
        self.equations = equations
        self.direction = direction
        self.terminal = terminal
        self.value = value

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
        root : float
            event occurs when this is zero from a given direction
        """
        raise NotImplementedError("Event class must define __call__")


class Inflation(Event):
    """ Inflation entry/exit """
    def __call__(self, t, y):
        return self.equations.inflating(t, y) - self.value


class Collapse(Event):
    """ Tests if H^2 is positive """
    def __call__(self, t, y):
        return self.equations.H2(t, y) - self.value


class ModeExit(Event):
    """ When mode exits the horizon aH """
    def __call__(self, t, y):
        logH = numpy.log(numpy.abs(self.equations.H2(t, y)))/2.
        N = self.equations.N(t, y)
        return logH + N - numpy.log(self.value)


class UntilN(Event):
    """ Stop at N """
    def __call__(self, t, y):
        return self.equations.N(t, y) - self.value
