from numpy.testing import assert_allclose, assert_almost_equal
from primordial.equations.t.inflation import Equations as Eqs_t, Inflation_start_initial_conditions as IC_t
from primordial.equations.N.inflation import Equations as Eqs_N, Inflation_start_initial_conditions as IC_N

from primordial.equations.events import Inflation
from primordial.equations.inflation_potentials import ChaoticPotential
from primordial.solver import solve

def test_t_n():
    m = 1
    K = 0
    V = ChaoticPotential()
    atol = 1e-8
    rtol = 1e-8

    N_0 = 0
    phi_0 = 17
    equations = Eqs_t(K, V)
    events = Inflation(equations, direction=-1, terminal=True)
    ic = IC_t(N_0, phi_0)
    sol_t = solve(equations, ic, events=events, atol=atol, rtol=rtol)

    equations = Eqs_N(K, V)
    events = Inflation(equations, direction=-1, terminal=True)
    ic = IC_N(N_0, phi_0)
    sol_N = solve(equations, ic, events=events, atol=atol, rtol=rtol)

    N = sol_N.N[1:-1]
    N_t = sol_t.N(sol_N.t(N))
    assert_allclose(N, N_t, atol*1e2, rtol*1e2)

    t = sol_t.t[1:-1]
    t_N = sol_N.t(sol_t.N(t))
    assert_allclose(t, t_N, atol*1e2, rtol*1e2)
