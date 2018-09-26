import numpy
from numpy.testing import assert_allclose, assert_almost_equal
from primordial.solver import Solver
from primordial.events import Inflation, Stationary
from primordial.t.background import Equations, KD_initial_conditions
from primordial.potentials import ChaoticPotential

def test_Solver():
    for K in [-1, 0, +1]:
        m = 1
        V = ChaoticPotential(m)
        equations = Equations(K, V)

        events= [Inflation(equations, +1), 
                 Inflation(equations, -1, True), 
                 Stationary(equations, terminal=True)]

        solver = Solver(equations, events)

        N_p = -1.5
        phi_p = 23
        t = numpy.logspace(-5,10,10**6)
        y0 = KD_initial_conditions(t[0], N_p, phi_p, equations)
        sol = solver.solve(t, y0)

        t = 1e-5
        sol = solver.solve(t, y0)
