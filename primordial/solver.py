import numpy
from scipy import integrate

def solve(equations, ic, interp1d_kwargs={}, *args, **kwargs):
    y0 = numpy.zeros(len(equations.i))
    ic(equations, y0)
    sol = integrate.solve_ivp(equations, (ic.t0, 1e300), y0, *args, **kwargs)
    return equations.sol(sol, **interp1d_kwargs)

