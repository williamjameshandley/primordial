========================================
primordial: inflationary equation solver
========================================
:primordial: inflationary equation solver
:Author: Will Handley
:Version: 0.0.7
:Homepage: https://github.com/williamjameshandley/primordial
:Documentation: http://primordial.readthedocs.io/

.. image:: https://travis-ci.org/williamjameshandley/primordial.svg?branch=master
   :target: https://travis-ci.org/williamjameshandley/primordial
   :alt: Build Status
.. image:: https://codecov.io/gh/williamjameshandley/primordial/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/williamjameshandley/primordial
   :alt: Test Coverage Status
.. image:: https://badge.fury.io/py/primordial.svg
   :target: https://badge.fury.io/py/primordial
   :alt: PyPi location
.. image:: https://readthedocs.org/projects/primordial/badge/?version=latest
   :target: https://primordial.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status


Description
===========

``primordial`` is a python package for solving cosmological inflationary equations.

It is very much in beta stage, and currently being built for research purposes.


Example Usage
=============

Plot Background evolution
-------------------------
.. code:: python

    import numpy
    import matplotlib.pyplot as plt
    from primordial.solver import Solver
    from primordial.equations.inflation_potentials import ChaoticPotential
    from primordial.equations.t.inflation import Equations, KD_initial_conditions
    from primordial.equations.events import Inflation, Stationary
    
    fig, ax = plt.subplots(3,sharex=True)
    for K in [-1, 0, +1]:
         m = 1
         V = ChaoticPotential(m)
         equations = Equations(K, V)

         events= [Inflation(equations),                    # Record inflation entry and exit 
                  Inflation(equations, -1, terminal=True), # Stop on inflation exit
                  Stationary(equations, terminal=True)]    # Stop if universe stops expanding


         solver = Solver(equations, events)

         N_p = -1.5
         phi_p = 23
         t_p = 1e-5
         ic = KD_initial_conditions(t_p, N_p, phi_p)
         t = numpy.logspace(-5,10,1e6)

         sol = solver.solve(ic, t)

         ax[0].plot(sol.N,sol.phi)
         ax[0].set_ylabel(r'$\phi$')

         ax[1].plot(sol.N,sol.H)
         ax[1].set_yscale('log')
         ax[1].set_ylabel(r'$H$')

         ax[2].plot(sol.N,1/(sol.H*numpy.exp(sol.N)))
         ax[2].set_yscale('log')
         ax[2].set_ylabel(r'$1/aH$')
         
    ax[-1].set_xlabel('$N$')

|image0|

Plot mode function evolution
----------------------------
.. code:: python

    import numpy
    import matplotlib.pyplot as plt
    from primordial.solver import solve
    from primordial.equations.inflation_potentials import ChaoticPotential
    from primordial.equations.t.mukhanov_sasaki import Equations, KD_initial_conditions
    from primordial.equations.events import Inflation, Stationary, ModeExit

    fig, axes = plt.subplots(3,sharex=True)
    for ax, K in zip(axes, [-1, 0, +1]):
    K=0
    ax = axes[1]
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

|image1|


.. |image0| image:: https://raw.githubusercontent.com/williamjameshandley/primordial/master/figures/background.png
.. |image1| image:: https://raw.githubusercontent.com/williamjameshandley/primordial/master/figures/ms.png 
