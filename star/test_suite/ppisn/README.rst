.. _ppisn:

*****
ppisn
*****

.. tags:: star, very-massive-star, helium-star, pair-instability, pulsational-pair-instability

This test evolves a :math:`72.5\,M_\odot` helium star at :math:`Z = 0.00142`
through its first pulsational pair-instability pulse and verifies the post-pulse relaxation
(see |Marchant2019|).

.. |Marchant2019| replace:: `Marchant et al. 2019 <https://ui.adsabs.harvard.edu/abs/2019ApJ...882...36M/abstract>`__

Run sequence
============

Running ``./rn`` performs two parts:

1. ``inlist_to_he_dep_header`` creates the helium-star model and evolves it
   until the central helium mass fraction falls below :math:`10^{-3}`. It
   writes ``he_dep.mod`` and refreshes ``standard_he_dep.mod``.
2. ``inlist_pulses_header`` loads ``he_dep.mod`` with velocity enabled,
   follows the first pulse, removes unbound surface ejecta directly, and
   relaxes the bound remnant.

Setting ``MESA_SKIP_OPTIONAL`` skips the first part and copies
``standard_he_dep.mod`` to ``he_dep.mod``. The second part stops 100 days
after the first post-pulse relaxation.

Configuration
=============

``inlist_ppisn`` contains the shared pulse physics. ``inlist_hydro_on`` and
``inlist_hydro_off`` are loaded by ``run_star_extras`` when the model enters
and leaves Riemann hydrodynamics.

The initial composition has ``initial_Y = 0.99858d0`` and no hydrogen.
The pulse calculation uses TDC with ``TDC_include_eturb_in_energy_equation = .false.``,
``approx21_cr60_plus_co56.net``, the deBoer
:math:`{}^{12}\mathrm{C}(\alpha,\gamma){}^{16}\mathrm{O}` rate, and
``mass_fraction_limit_for_Skye = 1d-8``. Level 3 gold tolerances relax to
``1d3`` above the temperature set by ``x_ctrl(23) = 5d9`` during evolution,
not during model relaxation.

The inlist profile interval is 200 models. ``run_star_extras`` retains its
100-model interval and 10-model interval near breakout. The test retains
the 6,000-model safety limit and the required 100-day termination check.

Last-Updated: 2026-09 by EbF
