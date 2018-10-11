import numpy
from scipy import integrate


def solve(equations, ic, interp1d_kwargs={}, *args, **kwargs):
    """ Solve differential equations

    This is a wrapper around ``scipy.integrate.solve_ivp``, with easier
    reusable objects for the equations and initial conditions.

    Parameters
    ----------
    equations: primordial.equations.Equations
        callable to compute the equations

    ic: initial conditions
        callable to set the initial conditions

    interp1d_kwargs: dict
        kwargs to pass to the interpolation functions

    All other arguments are identical to ``scipy.integrate.solve_ivp``

    Returns
    -------
    solution
        Monkey-patched version of the Bunch type usually returned by solve_ivp
    """
    y0 = numpy.zeros(len(equations.i))
    ic(equations, y0)
    sol = integrate.solve_ivp(equations, (ic.t0, 1e300), y0, *args, **kwargs)
    return equations.sol(sol, **interp1d_kwargs)
