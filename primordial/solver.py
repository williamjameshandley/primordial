import numpy
from scipy import integrate

def solve(equations, ic, *args, **kwargs):
    y0 = numpy.zeros(len(equations))
    ic(equations, y0)
    sol = integrate.solve_ivp(equations, (ic.t0, 1e300), y0, *args, **kwargs)
    return equations.sol(sol)

