.. _dev_TDC_single_star_ppisn:

**************************
dev_TDC_single_star_ppisn
**************************

.. tags:: star, very-massive-star, helium-star, pair-instability, pulsational-pair-instability

This development case evolves an isolated :math:`72.5\,M_\odot` helium star
from the pre-main sequence through helium depletion and then follows
successive pulsational pair-instability events (see |Marchant2019|).

.. |Marchant2019| replace:: `Marchant et al. 2019 <https://ui.adsabs.harvard.edu/abs/2019ApJ...882...36M/abstract>`__

Run sequence
============

Running ``./rn`` selects each part explicitly with ``MESA_INLIST``:

1. ``inlist_to_he_dep_header`` creates the helium star using the mass and
   metallicity in ``inlist_extra``. It stops at a central helium mass fraction
   of :math:`10^{-3}` and writes ``he_dep.mod``.
2. ``inlist_pulses_header`` loads ``he_dep.mod`` with velocity enabled and
   writes the final state to ``final.mod`` when the development run
   terminates.

Setting ``MESA_SKIP_OPTIONAL`` uses ``standard_he_dep.mod`` instead of
rebuilding the helium-depleted model. Running ``./re`` restarts the newest
photo with the pulse inlist.

Configuration
=============

``inlist_ppisn`` contains the persistent pulse physics. ``run_star_extras``
loads ``inlist_hydro_on`` and ``inlist_hydro_off`` as the model enters and
leaves Riemann hydrodynamics, removes unbound surface ejecta directly, and
continues through successive pulses.

Last-Updated: 2026-08 by EbF
