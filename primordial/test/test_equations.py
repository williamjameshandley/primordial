import numpy
from numpy.testing import assert_allclose, assert_almost_equal
from primordial.equations import Solver

def test_Solver():
    for K in [-1, 0, +1]:
        solver = Solver(K)
        N_p = -1.5
        phi_p = 23
        t = numpy.logspace(-5,10,10**6)
        sol = solver.solve(t, N_p, phi_p)

        t = 1e-5
        sol = solver.solve(t, N_p, phi_p)
