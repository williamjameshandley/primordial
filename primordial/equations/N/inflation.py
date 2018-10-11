import numpy
from primordial.equations.inflation import Equations as _Equations

class Equations(_Equations):
    """ Background equations in time 

    Solves bacgkround variables in cosmic time for curved and flat universes
    using the Klein-Gordon and Friedmann equations.

    Independent variable:
        N: e-folds (log a)
    
    Variables:
        phi: inflaton field
        dphi: d/dN (phi) 
        t: cosmic time
    
    """
    def __init__(self, K, potential):
        super(Equations, self).__init__(K, potential)
        self.set_independent_variable('N')
        self.add_variable('phi', 'dphi', 't')

    def __call__(self, N, y):
        """ The derivative function for underlying variables,
            computed using the Klein-Gordon equation """
        H2 = self.H2(N, y)
        dphi = self.dphi(N, y)
        dVdphi = self.dVdphi(N, y)
        dlogH = self.dlogH(N, y)
        ddphi = -(dlogH + 3)* dphi - dVdphi/H2

        dy = numpy.zeros_like(y)
        dy[self.i['phi']] = dphi
        dy[self.i['dphi']] = ddphi
        dy[self.i['t']] = 1./self.H(N, y)

        return dy

    def H2(self, N, y):
        """ The square of the Hubble parameter,
            computed using the Friedmann equation """
        V = self.V(N, y)
        dphi = self.dphi(N, y) 

        return (2*V - 6*self.K*numpy.exp(-2*N)) / (6. - dphi**2)

    def dlogH(self, N, y):
        """ d/dN log H """
        dphi = self.dphi(N, y) 
        H2 = self.H2(N, y) 
        return - dphi**2/2. + self.K * numpy.exp(-2*N)/H2

    def inflating(self, N, y):
        """ Inflation diagnostic """
        return self.V(N, y) / (self.V(N, y)/2. - self.K*numpy.exp(-2*N)) - self.dphi(N, y)**2 


class Inflation_start_initial_conditions(object):
    def __init__(self, N_e, phi_e):
        self.t0 = N_e
        self.phi_e = phi_e

    def __call__(self, equations, y0):
        V = equations.potential(self.phi_e)
        N = self.t0
        y0[equations.i['phi']] = self.phi_e
        y0[equations.i['dphi']] = -numpy.sqrt(V / (V/2. - equations.K*numpy.exp(-2*N)))
        y0[equations.i['t']] = 0.
