import matplotlib.pyplot as plt
from primordial.solver import solve
from primordial.events import Inflation, Stationary
from primordial.t.background import Equations, KD_initial_conditions
from primordial.potentials import ChaoticPotential
import numpy

fig, ax = plt.subplots(3,sharex=True)
for K in [-1, 0, +1]:
    m = 1
    V = ChaoticPotential(m)
    equations = Equations(K, V)

    events= [Inflation(equations),                   # Record inflation entry and exit
            Inflation(equations, -1, terminal=True), # Stop on inflation exit
            Stationary(equations, terminal=True)]    # Stop if universe stops expanding


    N_p = -1.5
    phi_p = 23
    t_p = 1e-5
    ic = KD_initial_conditions(t_p, N_p, phi_p)
    t = numpy.logspace(-5,10,1e6)

    sol = solve(equations, ic, t, events)

    ax[0].plot(sol.N,sol.phi)
    ax[0].set_ylabel(r'$\phi$')

    ax[1].plot(sol.N,sol.H)
    ax[1].set_yscale('log')
    ax[1].set_ylabel(r'$H$')

    ax[2].plot(sol.N,1/(sol.H*numpy.exp(sol.N)))
    ax[2].set_yscale('log')
    ax[2].set_ylabel(r'$1/aH$')

plt.show()
