import numpy
import matplotlib.pyplot as plt
from primordial.solver import solve
from primordial.equations.inflation_potentials import ChaoticPotential
from primordial.equations.t.mukhanov_sasaki import Equations, KD_initial_conditions
from primordial.equations.events import Inflation, Stationary, ModeExit

fig, axes = plt.subplots(3,sharex=True)
for ax, K in zip(axes, [-1, 0, +1]):
    ax2 = ax.twinx()
    m = 1
    V = ChaoticPotential(m)
    k = 100
    equations = Equations(K, V, k)

    events= [
            Inflation(equations),                    # Record inflation entry and exit
            Stationary(equations, terminal=True),    # Stop if universe stops expanding
            ModeExit(equations, -1, terminal=True, value=1e-1)   # Stop on mode exit
            ]


    N_p = -1.5
    phi_p = 23
    t_p = 1e-5
    ic = KD_initial_conditions(t_p, N_p, phi_p)
    t = numpy.logspace(-5,10,1e6)

    sol = solve(equations, ic, t, events)

    ax.plot(sol.N,sol.R1, 'k-')
    ax2.plot(sol.N,-numpy.log(sol.H*numpy.exp(sol.N)), 'b-')

    ax.set_ylabel('$\mathcal{R}$')
    ax2.set_ylabel('$-\log aH$')

    ax.text(0.9, 0.9, r'$K=%i$' % K, transform=ax.transAxes)

axes[-1].set_xlabel('$N$')
