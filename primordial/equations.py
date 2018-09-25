import numpy
from scipy import integrate

class Solver(object):
    method = 'RK45'
    def __init__(self, K, m=1):
        self.K = K
        self.m = m

    def equations(self, t, y):
        N, phi, dphi = y
        H2 = self.H2(y)
        if H2 < 0:
            print("negative H2")
        H = numpy.sqrt(H2)
        ddphi = -3*H*dphi - self.dVdphi(phi)
        return [H, dphi, ddphi]

    def events(self):
        return [self.inflation_exit]

    def inflation_exit(self, t, y):
        N, phi, dphi = y
        return self.V(phi) - dphi**2 
    inflation_exit.direction = -1
    inflation_exit.terminal = True

    def stationary(self, t, y):
        N, phi, dphi = y
        return self.H2(y)
    inflation_exit.terminal = True

    def H2(self, y):
        N, phi, dphi = y
        H2 = (dphi**2/2 + self.V(phi))/3 - self.K*numpy.exp(-2*N)
        return H2

    def solve(self, t, *args):
        try:
            t0 = t[0]
            t1 = t[-1]
            t_eval = t
        except TypeError:
            t0 = t
            t1 = 1e100
            t_eval = None
        y0 = self.y0(t0, *args)
        sol = integrate.solve_ivp(self.equations, (t0, t1), y0, 
                                  method=self.method, t_eval=t_eval, events=self.events())
        sol.N, sol.phi, sol.dphi = sol.y
        sol.H = numpy.sqrt(self.H2(sol.y))
        return sol

    def V(self, phi):
        return self.m**2*phi**2/2

    def dVdphi(self, phi):
        return self.m**2*phi

    def y0(self, t0, N_p, phi_p): 
        b = self.K * numpy.exp(-2*N_p)
        return [
                N_p + numpy.log(t0)/3 - 9/14 * b * t0**(4./3),
                phi_p - numpy.sqrt(2./3)*numpy.log(t0) - 27*numpy.sqrt(6)/56*b*t0**(4./3),
                -numpy.sqrt(2./3)/t0 - 9*numpy.sqrt(6)/14 * b * t0**(1./3),
                ]