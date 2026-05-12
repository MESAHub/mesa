.. _20M_pre_ms_to_cc:

*****************
20M_pre_ms_to_cc
*****************

This dev case is copied from ``star/test_suite/20M_pre_ms_to_core_collapse``.
It keeps a main ``inlist`` for the final evolution to core collapse,
starting from ``standard_lgTmax.mod``, with pgstar settings in
``inlist_pgstar``.

The source test suite contains seven parts from pre-main sequence to core
collapse.  This dev case flattens the shared settings, mass/Z settings,
and final core-collapse settings into the main ``inlist``. The pgstar
configuration is kept separately in ``inlist_pgstar``.

This case works with ``split_burn``, not without.

For production science we recommend adopting tighter mesh and timestep controls, such as those suggested in the comments of ``inlist``.

Physical checks
===============

None

Inlists
=======

This dev case has two inlist files.

* ``inlist`` loads ``standard_lgTmax.mod`` and evolves until ``fe_core_infall_limit``.
* ``inlist_pgstar`` contains the pgstar plot settings loaded by ``inlist``.

Last-Updated: 12May2026
