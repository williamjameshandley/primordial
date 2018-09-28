import numpy
import matplotlib.pyplot as plt
from primordial.solver import solve
from primordial.equations.t.inflation import Equations, Inflation_start_initial_conditions
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
ns_ = 0.96605
c = Omega_k_**2*As_

from scipy.interpolate import interp1d
def f(x):
    N_e, phi_e = x
    phi_e = numpy.abs(phi_e)
    print(x)

    V = ChaoticPotential()
    K = 1
    equations = Equations(K, V)
    events= [ModeExit(equations, direction=+1, value=kstar), Inflation(equations, -1, terminal=True)]
    ic = Inflation_start_initial_conditions(N_e, phi_e)

    sol = solve(equations, ic, events=events, atol=1e-8, rtol=1e-8)
    if len(sol.t_events[0])==0:
        print("Too small")
        phi_star = phi_e
    else:
        phi = interp1d(sol.t, sol.phi)
        phi_star = phi(sol.t_events[0])[0]

    eps = (V.d(phi_star)/V(phi_star))**2/2
    eta = V.dd(phi_star)/V(phi_star)
    ns = 1 + 2*eta - 6*eps
    As = V(phi_star) / eps / 24 / numpy.pi**2
    Omega_k = -K / ( H0 * sol.H[-1] * numpy.exp(2*sol.N[-1]) ) * numpy.sqrt(Omega_r)

    return [numpy.log(ns/ns_), numpy.log(Omega_k**2*As) - numpy.log(c)]

from scipy.optimize import root
f([0,6])

sol = root(f, [0, 20])
print(sol)

#m = numpy.exp(sol.x[0])
#N_p = sol.x[1] - 2/3 * numpy.log(m) 
#phi_p = sol.x[2] - numpy.sqrt(2/3) * numpy.log(m) 
#t_p = 1e-5/m
#print("solution", t_p, m,  N_p, phi_p)
#
#m, N_p, phi_p = 7.730640614653029e-06, 7.691660989308948, 23.771111141140796
#t_p = 1
#
#V = ChaoticPotential(m)
#K = 1
#equations = Equations(K, V)
#
#events= [Inflation(equations),                   # Record inflation entry and exit
#        Inflation(equations, -1, terminal=True), # Stop on inflation exit
#        Stationary(equations, terminal=True)]    # Stop if universe stops expanding
#
#ic = KD_initial_conditions(t_p, N_p, phi_p)
#sol = solve(equations, ic, events=[ModeExit(equations, direction=+1, value=kstar), Inflation(equations, -1, terminal=True)], atol=1e-8, rtol=1e-8)
#sol
#fig, ax = plt.subplots()
#ax.plot(sol.t,sol.N)
#ax.set_xscale('log')
