import matplotlib.pyplot as plt
from primordial.equations import Solver
import numpy

fig, ax = plt.subplots(3,sharex=True)
for K in [-1, 0, +1]:
    solver = Solver(K)

    N_p = -1.27
    phi_p = 16
    t = numpy.logspace(-5,10,1e6)
    sol = solver.solve(t, N_p, phi_p)
    print(len(sol.t))

    ax[0].plot(sol.N,sol.phi)
    ax[0].set_ylabel(r'$\phi$')

    ax[1].plot(sol.N,sol.H)
    ax[1].set_yscale('log')
    ax[1].set_ylabel(r'$H$')

    ax[2].plot(sol.N,1/(sol.H*numpy.exp(sol.N)))
    ax[2].set_yscale('log')
    ax[2].set_ylabel(r'$1/aH$')

plt.show()
