.. _ppisn:

*****
ppisn
*****

.. tags:: star, very-massive-star, helium-star, pair-instability, pulsational-pair-instability

This test evolves a :math:`72.5\,M_\odot` helium star through its first
pulsational pair-instability pulse and verifies the post-pulse relaxation
(see |Marchant2019|).

.. |Marchant2019| replace:: `Marchant et al. 2019 <https://ui.adsabs.harvard.edu/abs/2019ApJ...882...36M/abstract>`__

Run sequence
============

Running ``./rn`` performs two parts through the standard test-suite
``do_one`` helper:

1. ``inlist_to_he_dep_header`` creates the helium-star model and evolves it
   until the central helium mass fraction falls below :math:`10^{-3}`. It
   writes ``he_dep.mod`` and refreshes ``standard_he_dep.mod``.
2. ``inlist_pulses_header`` loads ``he_dep.mod`` with velocity enabled,
   follows the first pulse, removes unbound surface ejecta directly, and
   relaxes the bound remnant.

Setting ``MESA_SKIP_OPTIONAL`` skips the first part and copies
``standard_he_dep.mod`` to ``he_dep.mod``. The test succeeds 100 days after
the first post-pulse relaxation, requires the corresponding termination
code, and limits the calculation to 15000 models.

Configuration
=============

``inlist_ppisn`` contains the shared pulse physics. ``inlist_hydro_on`` and
``inlist_hydro_off`` are loaded by ``run_star_extras`` when the model enters
and leaves Riemann hydrodynamics. Profiles are written every 200 models.

The pulse calculation uses TDC with turbulent energy in the energy equation,
``approx21_cr60_plus_co56.net``, the deBoer
:math:`{}^{12}\mathrm{C}(\alpha,\gamma){}^{16}\mathrm{O}` rate, and
``mass_fraction_limit_for_Skye = 1d-5``.

Last-Updated: 2026-08 by EbF
