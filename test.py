import numpy
import matplotlib.pyplot as plt
from primordial.solver import solve
from primordial.equations.N.inflation import Equations, Inflation_start_initial_conditions
from primordial.equations.events import Inflation, Stationary, ModeExit
from primordial.equations.inflation_potentials import ChaoticPotential
from primordial.units import H0_h, kstar

h = 0.6732
Omega_b = 0.022383     / h**2
Omega_m = 0.14314      / h**2
Omega_r = 4.18343e-5   / h**2 # Assumes T = 2.7255, and includes neutrinos
Omega_k_ = -0.1
Omega_l = 1-Omega_m-Omega_r-Omega_k_
H0 = H0_h * h
logA = 3.0448
As_ = 1e-10*numpy.exp(logA)
logAs_ = numpy.log(As_)
ns_ = 0.96605
c = Omega_k_**2/As_



def trans(x, y, K, V):
    #phi_e = phi_min + numpy.log(1+x+numpy.sqrt(1+x**2))
    phi_min = 7
    phi_e = phi_min + x+numpy.sqrt(1+x**2)
    Nmax = numpy.log(2*(kstar**2 + K)/V(phi_e))/2
    Nmin = -numpy.log(V(phi_e)/2)/2
    N_e = Nmin + (Nmax-Nmin)*numpy.exp(y)/(1+numpy.exp(y))
    print(N_e, phi_e)
    return N_e, phi_e


def integrate_from_inflation_start(N_e, phi_e, K, V):
    equations = Equations(K, V)
    events= [ModeExit(equations, direction=+1, value=kstar), Inflation(equations, -1, terminal=True), Stationary(equations, terminal=True)]
    ic = Inflation_start_initial_conditions(N_e, phi_e)
    sol = solve(equations, ic, events=events, atol=1e-8, rtol=1e-8)
    return sol


from scipy.interpolate import interp1d
def compute_pps(sol, K, V):
    if len(sol.t_events[0])==0:
        return numpy.nan, numpy.nan, numpy.nan
        print("Too small")
        phi_star = sol.phi[-1]
    else:
        phi = interp1d(sol.t, sol.phi)
        phi_star = phi(sol.t_events[0])[0]
    eps = (V.d(phi_star)/V(phi_star))**2/2
    eta = V.dd(phi_star)/V(phi_star)
    ns = 1 + 2*eta - 6*eps
    logAs = numpy.log(V(phi_star) / eps / 24 / numpy.pi**2)
    logOmega_k = numpy.log(Omega_r)/2 - numpy.log(H0 * sol.H[-1]) - 2*sol.N[-1] 
    return logAs, ns, logOmega_k


def f(x):
    K = 1
    m = 1
    V = ChaoticPotential(m)
    N_e, phi_e = trans(x[0],x[1], K, V)
    sol = integrate_from_inflation_start(N_e, phi_e, K, V)
    logAs, ns, logOmega_k = compute_pps(sol, K, V)

    return [ns - ns_, 2*logOmega_k - logAs - numpy.log(c)]



from scipy.optimize import root
sol = root(f, [0, 0])

K = 1
m = 1
V = ChaoticPotential(m)
N_e = 2.8091921072150887
phi_e = 15.978200860209245
sol = integrate_from_inflation_start(N_e, phi_e, K, V)
logAs, ns, logOmega_k = compute_pps(sol, K, V)
m = numpy.exp(-(logAs - logAs_)/2)

V = ChaoticPotential(m)
N_e = 2.8091921072150887 - numpy.log(m)
phi_e = 15.978200860209245
sol = integrate_from_inflation_start(N_e, phi_e, K, V)
logAs, ns, logOmega_k = compute_pps(sol, K, V)
print("-----------------")
print(logAs - logAs_)
print(ns - ns_)
print(logOmega_k - numpy.log(-Omega_k_))
