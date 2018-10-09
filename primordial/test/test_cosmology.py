from numpy.testing import assert_allclose, assert_almost_equal
from primordial.equations.t.cosmology import Equations as Eqs_t, initial_conditions as IC_t
from primordial.equations.N.cosmology import Equations as Eqs_N, initial_conditions as IC_N

from primordial.equations.events import UntilN
from primordial.solver import solve
from primordial.units import H0_h

def test_cosmology():
    atol, rtol = 1e-8, 1e-8

    h = 0.6732
    Omega_b = 0.022383     / h**2
    Omega_m = 0.14314      / h**2
    Omega_r = 4.18343e-5   / h**2 # Assumes T = 2.7255, and includes neutrinos
    Omega_k = -0.01
    Omega_l = 1-Omega_m-Omega_r-Omega_k
    H0 = H0_h * h
    Ni = -20

    # Solve in N
    equations = Eqs_N(H0, Omega_r, Omega_m, Omega_l)
    ic = IC_N(Ni)
    events = UntilN(equations, terminal=True)
    sol_N = solve(equations, ic, events=events, atol=atol, rtol=rtol)

    N = sol_N.N[1:-1]
    t = sol_N.t(N)
    print(t)

    # Solve in t
    equations = Eqs_t(H0, Omega_r, Omega_m, Omega_l)
    ic = IC_t(Ni)
    events = UntilN(equations, terminal=True)
    sol_t = solve(equations, ic, events=events, atol=atol, rtol=rtol)

    ## Check t-N consistency
    #assert_allclose(sol_t.N(t), N, atol, rtol)
    #assert_allclose(t, sol_N.t(N), atol, rtol)
