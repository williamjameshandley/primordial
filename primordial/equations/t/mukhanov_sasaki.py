import numpy
from primordial.equations.t.inflation import Equations as BackgroundEquations, KD_initial_conditions as BackgroundKD_initial_conditions

class Equations(BackgroundEquations):
    def __init__(self, K, potential, k):
        super(Equations, self).__init__(K, potential)
        self.k = k
        self.add_variables(['R1', 'dR1', 'R2', 'dR2'])


    def __call__(self, t, y):
        """ The derivative function for underlying variables,
            computed using the Mukhanov-Sasaki equation equation """

        # Compute background variables
        dy = super(Equations, self).__call__(t, y)

        # Get useful variables
        N = self.N(t, y)
        H = self.H(t, y)
        dphi = self.dphi(t, y)
        ddphi = dy[self['dphi']]
        R1 = self.R1(t, y)
        R2 = self.R2(t, y)
        dR1 = self.dR1(t, y)
        dR2 = self.dR2(t, y)

        k = self.k
        K = self.K
        if self.K > 0:
            k2 = k*(k+2)-3*K
        else:
            k2 = k**2-3*K

        E = dphi**2 / H**2 / 2
        alpha = E*K + k2
        gamma = (k2**2 - k2*K - E * k2*K - 2*k2*ddphi*K/dphi/H - E*k2 + 2*k2*K**2/H**2*numpy.exp(-2*N))*numpy.exp(-2*N)
        beta = (3*H*k2 + 2*E*H*k2 + 2*k2*ddphi/dphi + 3*E*H*K - 2*k2*K/H*numpy.exp(-2*N))

        ddR1 = -(beta*dR1 + gamma*R1)/alpha
        ddR2 = -(beta*dR2 + gamma*R2)/alpha

        dy[self['R1']] = dR1
        dy[self['dR1']] = ddR1
        dy[self['R2']] = dR2
        dy[self['dR2']] = ddR2

        return dy

class KD_initial_conditions(BackgroundKD_initial_conditions):
    def __call__(self, equations, y0):
        super(KD_initial_conditions, self).__call__(equations, y0)
        y0[equations['R1']] = 0
        y0[equations['dR1']] = equations.k
        y0[equations['R2']] = 1
        y0[equations['dR2']] = 0
