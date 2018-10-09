from numpy.testing import assert_allclose, assert_almost_equal
from primordial.equations.t.inflation import Equations as Eqs_t, Inflation_start_initial_conditions as IC_t
from primordial.equations.N.inflation import Equations as Eqs_N, Inflation_start_initial_conditions as IC_N

from primordial.equations.events import Inflation
from primordial.equations.inflation_potentials import ChaoticPotential
from primordial.solver import solve

def test_inflation():
    atol, rtol = 1e-8, 1e-8

    V = ChaoticPotential(m=1)
    N_0, phi_0 = 0, 17

    for K in [-1, 0, 1]:
        # Solve in N
        equations = Eqs_N(K, V)
        ic = IC_N(N_0, phi_0)
        events = Inflation(equations, direction=-1, terminal=True)
        sol_N = solve(equations, ic, events=events, atol=atol*1e-1, rtol=rtol*1e-1)

        N = sol_N.N[1:-1]
        t = sol_N.t(N)

        # Solve in t
        equations = Eqs_t(K, V)
        ic = IC_t(N_0, phi_0)
        events = Inflation(equations, direction=-1, terminal=True)
        sol_t = solve(equations, ic, events=events, t_eval=t, atol=atol*1e-1, rtol=rtol*1e-1)

        # Check t-N consistency
        assert_allclose(sol_t.N(t), N, atol, rtol)
        assert_allclose(t, sol_N.t(N), atol, rtol)

        # Check phi
        assert_allclose(sol_t.phi(t), sol_N.phi(N), atol, rtol)

        # Check H
        assert_allclose(sol_N.H(N), sol_t.H(t), atol, rtol)
