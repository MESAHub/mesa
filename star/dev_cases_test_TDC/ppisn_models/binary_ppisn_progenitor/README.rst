.. _binary_ppisn_progenitor:

************************
binary_ppisn_progenitor
************************

.. tags:: star, binary, massive-star, mass-transfer, pair-instability

This development case constructs the binary-evolution progenitor used by
``dev_TDC_ppisn_from_binary``. It evolves an initially
:math:`100\,M_\odot + 50\,M_\odot` binary with an orbital period of 10 days
and metallicity :math:`Z = 0.00146`.

Evolution and model handoff
===========================

``inlist_binary`` selects ``inlist1`` for the initial donor and ``inlist2``
for the initial accretor. Both stars share the physics in ``inlist_both`` and
``inlist_extra``. At donor helium depletion, ``run_binary_extras`` writes
``donor_postHe.mod`` and ``accretor_postHe.mod``, replaces star 1 with a
:math:`2.6\,M_\odot` point mass, and continues evolving star 2.

The configured final model names are ``MASS1donor_final.mod`` and
``MASS2accretor_final.mod``. The latter is the default input model for
``dev_TDC_ppisn_from_binary``.

Run sequence
============

Running ``./rn`` selects ``inlist_binary`` explicitly with ``MESA_INLIST``.
That file also loads ``inlist_pgbinary``. Running ``./re`` restarts the most
recent complete set of binary and stellar photos with the same inlist.

The model uses TDC with turbulent energy in the energy equation,
``approx21_cr60_plus_co56.net``, the deBoer
:math:`{}^{12}\mathrm{C}(\alpha,\gamma){}^{16}\mathrm{O}` rate, and
``mass_fraction_limit_for_Skye = 1d-5``.

Last-Updated: 2026-08 by EbF
