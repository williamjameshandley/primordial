from scipy import integrate

class Solver(object):
    method = 'RK45'

    def __init__(self, equations, events=None):
        self.equations = equations
        self.events = events

    def solve(self, t, y0):
        try:
            t0 = t[0]
            t1 = t[-1]
            t_eval = t
        except TypeError:
            t0 = t
            t1 = 1e100
            t_eval = None
        sol = integrate.solve_ivp(self.equations, (t0, t1), y0, 
                                  method=self.method, t_eval=t_eval, events=self.events)
        return self.equations.sol(sol)

    def events(self):
        return None
