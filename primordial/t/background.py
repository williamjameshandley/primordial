import numpy
from primordial.equations import Equations as _Equations

def make_variable(name, i):
    def method(self, t, y):
        return y[self[name]]
    i[name] = len(i)
    return method


class Equations(_Equations):
    """ Background equations in time 
    
    Variables:
        N: efolds
        phi: inflaton field
        dphi: d (phi) / dt
    
    """
    N = make_variable('N', _Equations._i)
    phi = make_variable('phi', _Equations._i)
    dphi = make_variable('dphi', _Equations._i)

    def __init__(self, K, potential):
        super(Equations, self).__init__(potential)
        self.K = K

    def __call__(self, t, y):
        """ The derivative function for underlying variables """
        H = self.H(t, y)
        dphi = self.dphi(t, y)
        dVdphi = self.dVdphi(t, y)
        ddphi = -3*H*dphi - dVdphi

        dy = numpy.zeros_like(y)
        dy[self['N']] = H
        dy[self['phi']] = dphi
        dy[self['dphi']] = ddphi

        return dy

    def H(self, t, y):
        return numpy.sqrt(self.H2(t, y))

    def H2(self, t, y):
        """ The square of the Hubble constant """
        N = self.N(t, y)
        V = self.V(t, y)
        dphi = self.dphi(t, y) 
        N = self.N(t, y) 
        return (dphi**2/2 + V)/3 - self.K*numpy.exp(-2*N)

    def inflating(self, t, y):
        """ Inflation diagnostic """
        return self.V(t, y) - self.dphi(t, y)**2 

    def sol(self, sol):
        """ Post-process solution of solve_ivp """
        sol.N, sol.phi, sol.dphi = sol.y
        sol.H = numpy.sqrt(self.H2(sol.t, sol.y))
        return sol

    @property
    def n(self):
        return len(self._i)


def KD_initial_conditions(t0, N_p, phi_p, equations): 
    """ Kinetic dominance initial conditions """
    b = equations.K * numpy.exp(-2*N_p)
    y0 = [numpy.nan]* equations.n
    y0[equations['N']] = N_p + numpy.log(t0)/3 - 9/14 * b * t0**(4./3)
    y0[equations['phi']] = phi_p - numpy.sqrt(2./3)*numpy.log(t0) - 27*numpy.sqrt(6)/56*b*t0**(4./3)
    y0[equations['dphi']] = -numpy.sqrt(2./3)/t0 - 9*numpy.sqrt(6)/14 * b * t0**(1./3)
    return y0
