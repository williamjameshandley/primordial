import numpy
from numpy.testing import assert_allclose, assert_almost_equal
from primordial.solver import SolverKD_t as Solver
from primordial.potentials import ChaoticPotential

def test_Solver():
    for K in [-1, 0, +1]:
        m = 1
        V = ChaoticPotential(m)
        solver = Solver(K, V)
        N_p = -1.5
        phi_p = 23
        t = numpy.logspace(-5,10,10**6)
        sol = solver.solve(t, N_p, phi_p)

        t = 1e-5
        sol = solver.solve(t, N_p, phi_p)
