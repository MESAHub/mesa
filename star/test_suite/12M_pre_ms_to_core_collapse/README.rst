.. _12M_pre_ms_to_core_collapse:

***************************
12M_pre_ms_to_core_collapse
***************************

This test suite evolves a solar metalicity 12 |MSun| model from the pre-ms to core collapse.

Physical checks
===============

None

Inlists
=======

This test case has six parts.

* Part 1 (``inlist_make_late_pre_zams``) creates a 12 |Msun|, solar metallicity, pre-main sequence model and evolves it for 100 years.

* Part 2 (``inlist_to_zams``) evolves the model to the zero age main sequence.

* Part 3 (``inlist_to_end_core_he_burn``) takes the model to core helium depletion.

* Part 4 (``inlist_to_end_core_c_burn``) takes the model to core carbon depletion.

* Part 5 (``inlist_to_lgTmax``) evolves the model until the core temperature reaches log T =9.55 (approximately silicon burning)

* Part 6 (``inlist_to_cc``) evolves until core collapse.


Last-Updated: 15Sept2021 by rf

