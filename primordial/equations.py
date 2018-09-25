import numpy
from scipy import integrate

class Solver(object):
    method = 'RK45'
    def __init__(self, K, m=1):
        self.K = K
        self.m = m

    def equations(self, t, y):
        N, phi, dphi = y
        H = self.H(y)
        ddphi = -3*H*dphi - self.dVdphi(phi)
        return [H, dphi, ddphi]

    def events(self):
        return [self.inflation_exit]

    def inflation_exit(self, t, y):
        N, phi, dphi = y
        return self.V(phi) - dphi**2 
    inflation_exit.direction = -1
    inflation_exit.terminal = True

    def H(self, y):
        N, phi, dphi = y
        return numpy.sqrt((dphi**2/2 + self.V(phi))/3 - self.K*numpy.exp(-2*N))

    def solve(self, t0, t, *args):
        y0 = self.y0(t0, *args)
        sol = integrate.solve_ivp(self.equations, (t0,t), y0, 
                                  method=self.method, events=self.events())
        sol.N, sol.phi, sol.dphi = sol.y
        sol.H = self.H(sol.y)
        return sol

    def V(self, phi):
        return self.m**2*phi**2/2

    def dVdphi(self, phi):
        return self.m**2*phi

    def y0(self, t0, N_p, phi_p): 
        b = self.K * numpy.exp(-2*N_p)
        return [
                N_p + numpy.log(t0)/3 - 9/14 * b * t0**(4/3),
                phi_p - numpy.sqrt(2/3)*numpy.log(t0) - 27*numpy.sqrt(6)/56*b*t0**(4/3),
                -numpy.sqrt(2/3)/t0 - 9*numpy.sqrt(6)/14 * b * t0**(1/3),
                ]


solver = Solver(-1)

N_p = 0
phi_p = 23
sol = solver.solve(1e-5, 1e4, 0, 23)
print(sol.H)

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(sol.t,sol.phi)
ax.set_xscale('log')
ax.set_yscale('log')
