.. _dev_TDC_ppisn_from_binary:

**************************
dev_TDC_ppisn_from_binary
**************************

.. tags:: star, very-massive-star, helium-star, pair-instability, pulsational-pair-instability

This development case follows pulsational pair instability in the accretor
produced by the sibling ``binary_ppisn_progenitor`` calculation (see
|Marchant2019|). It is a single-star pulse calculation; it does not evolve the
binary or construct the progenitor in this directory.

.. |Marchant2019| replace:: `Marchant et al. 2019 <https://ui.adsabs.harvard.edu/abs/2019ApJ...882...36M/abstract>`__

Model handoff
=============

``inlist_pulses`` loads ``MASS2accretor_final.mod`` from the binary
progenitor case. The included ``accretor_postHe.mod`` is an intermediate
binary-evolution snapshot and is not the default starting model selected by
``./rn``.

Run sequence
============

Running ``./rn`` selects ``inlist_pulses_header`` with ``MESA_INLIST``, loads
the binary product with velocity enabled, and writes ``final.mod`` when the
development run terminates. Running ``./re`` restarts the newest photo with
the same pulse inlist.

``inlist_ppisn`` contains the persistent pulse physics.
``run_star_extras`` switches between ``inlist_hydro_on`` and
``inlist_hydro_off`` and removes unbound surface ejecta directly. The model
continues through successive pulses without the test suite's 100-day success
termination or model-number limit. Profiles are written every 200 models.

Last-Updated: 2026-08 by EbF
