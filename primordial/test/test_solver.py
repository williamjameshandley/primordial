import numpy
from numpy.testing import assert_allclose, assert_almost_equal
from primordial.solver import solve
from primordial.equations.inflation_potentials import ChaoticPotential
from primordial.equations.t.inflation import Equations, KD_initial_conditions
from primordial.equations.events import Inflation, Stationary

def test_Solver():
    for K in [-1, 0, +1]:
        m = 1
        V = ChaoticPotential(m)
        equations = Equations(K, V)

        events= [Inflation(equations, +1), 
                 Inflation(equations, -1, True), 
                 Stationary(equations, terminal=True)]


        t_p = 1e-5
        N_p = -1.5
        phi_p = 23
        ic = KD_initial_conditions(t_p, N_p, phi_p)

        t = numpy.logspace(-5,10,10**6)
        sol = solve(equations, ic, t_eval=t, events=events)
        sol = solve(equations, ic, events=events)
