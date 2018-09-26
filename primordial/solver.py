import numpy
from scipy import integrate

class Solver(object):
    method = 'RK45'

    def __init__(self, equations, events=None):
        self.equations = equations
        self.events = events

    def solve(self, ic, t=None):
        y0 = numpy.zeros(len(self.equations))
        ic(self.equations, y0)
        sol = integrate.solve_ivp(self.equations, (ic.t0, 1e300), y0, 
                                  method=self.method, t_eval=t, events=self.events)
        return self.equations.sol(sol)

def solve(equations, ic, t=None, events=None, method='RK45'):
    solver = Solver(equations, events)
    return solver.solve(ic, t)

