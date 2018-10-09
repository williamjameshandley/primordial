from numpy.testing import assert_allclose, assert_almost_equal
from primordial.equations.t.cosmology import Equations as Eqs_t, initial_conditions as IC_t
from primordial.equations.N.cosmology import Equations as Eqs_N, initial_conditions as IC_N

from primordial.equations.events import UntilN
from primordial.solver import solve
from primordial.units import H0_h, tp

def test_cosmology():
    atol, rtol = 1e-8, 1e-8

    h = 0.6732
    Omega_b = 0.022383     / h**2
    Omega_m = 0.14314      / h**2
    Omega_r = 4.18343e-5   / h**2 # Assumes T = 2.7255, and includes neutrinos
    H0 = H0_h * h

    for Omega_k in [-0.1, 0, 0.1]:
        Omega_l = 1-Omega_m-Omega_r-Omega_k

        # Solve in N
        equations = Eqs_N(H0, Omega_r, Omega_m, Omega_k, Omega_l)

        Ni = 70 if Omega_k else -70
        ic = IC_N(Ni)

        Ne = equations.N0 if Omega_k else 0
        events = UntilN(equations, terminal=True, value=Ne)
        sol_N = solve(equations, ic, events=events, atol=atol*1e-1, rtol=rtol*1e-1)

        N = sol_N.N[1:-1]
        t = sol_N.t(N)

        # Solve in t
        equations = Eqs_t(H0, Omega_r, Omega_m, Omega_k, Omega_l)
        ic = IC_t(Ni)
        events = UntilN(equations, terminal=True, value=Ne)
        sol_t = solve(equations, ic, events=events, t_eval=t, atol=atol*1e-1, rtol=rtol*1e-1)

        ## Check t-N consistency
        assert_allclose(sol_t.N(t), N, atol, rtol)
        assert_allclose(t, sol_N.t(N), atol, rtol)

        # Check H
        assert_allclose(sol_N.H(N), sol_t.H(t), atol*10, rtol*10)
