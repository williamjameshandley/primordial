import numpy
from primordial.equations.N.inflation import Equations as BackgroundEquations, Inflation_start_initial_conditions as BackgroundInflation_start_initial_conditions

class Equations(BackgroundEquations):
    def __init__(self, K, potential, k):
        super(Equations, self).__init__(K, potential)
        self.k = k
        self.add_variable('R1', 'dR1', 'R2', 'dR2')


    def __call__(self, N, y):
        """ The derivative function for underlying variables,
            computed using the Mukhanov-Sasaki equation equation """

        # Compute background variables
        dy = super(Equations, self).__call__(N, y)

        # Get useful variables
        H2 = self.H2(N, y)
        dphi = self.dphi(N, y)
        ddphi = dy[self.i['dphi']]
        R1 = self.R1(N, y)
        R2 = self.R2(N, y)
        dR1 = self.dR1(N, y)
        dR2 = self.dR2(N, y)

        k = self.k
        K = self.K
        if self.K > 0:
            k2 = k*(k+2)-3*K
        else:
            k2 = k**2-3*K

        E = dphi**2 / 2.
        alpha = E*K + k2
        gamma = (K*(k2-K)*dphi**2/2. - 2*K*k2*ddphi/dphi - K*k2 + k2**2)*numpy.exp(-2*N)/H2
        beta = ((dphi**2*K/2. + k2)*K/H2*numpy.exp(-2*N) + K * (3/2.-dphi**2/4.)*dphi**2  + (2*ddphi/dphi - dphi**2/2. + 3)*k2)

        ddR1 = -(beta*dR1 + gamma*R1)/alpha
        ddR2 = -(beta*dR2 + gamma*R2)/alpha

        dy[self.i['R1']] = dR1
        dy[self.i['dR1']] = ddR1
        dy[self.i['R2']] = dR2
        dy[self.i['dR2']] = ddR2

        return dy


class Inflation_start_initial_conditions(BackgroundInflation_start_initial_conditions):
    def __call__(self, equations, y0):
        super(Inflation_start_initial_conditions, self).__call__(equations, y0)
        y0[equations.i['R1']] = 0
        y0[equations.i['dR1']] = equations.k
        y0[equations.i['R2']] = 1
        y0[equations.i['dR2']] = 0
