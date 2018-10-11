import numpy
from primordial.equations.inflation import Equations as _Equations

class Equations(_Equations):
    """ Background equations in time 

    Solves bacgkround variables in cosmic time for curved and flat universes
    using the Klein-Gordon and Friedmann equations.

    Independent variable:
        t: cosmic time
    
    Variables:
        N: efolds
        phi: inflaton field
        dphi: d (phi) / dt
    
    """
    def __init__(self, K, potential):
        super(Equations, self).__init__(K, potential)
        self.set_independent_variable('t')
        self.add_variable('N', 'phi', 'dphi')

    def __call__(self, t, y):
        """ The derivative function for underlying variables,
            computed using the Klein-Gordon equation """
        H = self.H(t, y)
        dphi = self.dphi(t, y)
        dVdphi = self.dVdphi(t, y)
        ddphi = -3*H*dphi - dVdphi

        dy = numpy.zeros_like(y)
        dy[self.i['N']] = H
        dy[self.i['phi']] = dphi
        dy[self.i['dphi']] = ddphi

        return dy

    def H2(self, t, y):
        """ The square of the Hubble parameter,
            computed using the Friedmann equation """
        N = self.N(t, y)
        V = self.V(t, y)
        dphi = self.dphi(t, y) 
        N = self.N(t, y) 
        return (dphi**2/2. + V)/3. - self.K*numpy.exp(-2*N)

    def inflating(self, t, y):
        """ Inflation diagnostic """
        return self.V(t, y) - self.dphi(t, y)**2 


class KD_initial_conditions(object):
    def __init__(self, t0, N_p, phi_p):
        self.t0 = t0
        self.N_p = N_p
        self.phi_p = phi_p

    def __call__(self, equations, y0):
        t0 = self.t0
        b = equations.K * numpy.exp(-2*self.N_p)
        y0[equations.i['N']] = self.N_p + numpy.log(t0)/3 - 9./14 * b * t0**(4./3)
        y0[equations.i['phi']] = self.phi_p - numpy.sqrt(2./3)*numpy.log(t0) - 27*numpy.sqrt(6)/56*b*t0**(4./3)
        y0[equations.i['dphi']] = -numpy.sqrt(2./3)/t0 - 9*numpy.sqrt(6)/14 * b * t0**(1./3)


class Inflation_start_initial_conditions(object):
    def __init__(self, N_e, phi_e):
        self.t0 = 0.
        self.N_e = N_e
        self.phi_e = phi_e

    def __call__(self, equations, y0):
        y0[equations.i['N']] = self.N_e
        y0[equations.i['phi']] = self.phi_e
        y0[equations.i['dphi']] = -numpy.sqrt(equations.potential(self.phi_e))
