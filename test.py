import matplotlib.pyplot as plt
from primordial.solver import SolverKD_t as Solver
from primordial.potentials import ChaoticPotential
import numpy

fig, ax = plt.subplots(3,sharex=True)
for K in [-1, 0, +1]:
    m = 1
    V = ChaoticPotential(m)
    solver = Solver(K, V)

    N_p = -1.5
    phi_p = 23
    t = numpy.logspace(-5,10,1e6)
    sol = solver.solve(t, N_p, phi_p)
    print(sol.t_inflation_entry, sol.t_inflation_exit)

    ax[0].plot(sol.N,sol.phi)
    ax[0].set_ylabel(r'$\phi$')

    ax[1].plot(sol.N,sol.H)
    ax[1].set_yscale('log')
    ax[1].set_ylabel(r'$H$')

    ax[2].plot(sol.N,1/(sol.H*numpy.exp(sol.N)))
    ax[2].set_yscale('log')
    ax[2].set_ylabel(r'$1/aH$')

plt.show()
