import matplotlib.pyplot as plt
from primordial.solver import Solver
from primordial.events import Inflation, Stationary
from primordial.t.background import Equations, KD_initial_conditions
from primordial.potentials import ChaoticPotential
import numpy

fig, ax = plt.subplots(3,sharex=True)
for K in [-1, 0, +1]:
    m = 1
    V = ChaoticPotential(m)
    equations = Equations(K, V)

    events= [Inflation(equations, +1),               # Record inflation entry
            Inflation(equations, -1, terminal=True), # Stop on inflation exit
            Stationary(equations, terminal=True)]    # Stop if universe stops expanding


    solver = Solver(equations, events)

    N_p = -1.5
    phi_p = 23
    t = numpy.logspace(-5,10,1e6)

    y0 = KD_initial_conditions(t[0], N_p, phi_p, equations)

    sol = solver.solve(t, y0)

    ax[0].plot(sol.N,sol.phi)
    ax[0].set_ylabel(r'$\phi$')

    ax[1].plot(sol.N,sol.H)
    ax[1].set_yscale('log')
    ax[1].set_ylabel(r'$H$')

    ax[2].plot(sol.N,1/(sol.H*numpy.exp(sol.N)))
    ax[2].set_yscale('log')
    ax[2].set_ylabel(r'$1/aH$')

plt.show()
