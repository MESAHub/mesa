.. _20M_pre_ms_to_core_collapse:

***************************
20M_pre_ms_to_core_collapse
***************************

This test suite evolves a solar metalicity 20 |MSun| model from the pre-ms to core collapse.
For bit for bit convergence, we recommendedd to run by using the ./run_all script instead of restarting from models,
see https://github.com/MESAHub/mesa/issues/610.

This test_suite has been tested up to 80 solar masses, up to solar metallicity, with mass loss, and produces reasonable HR-tracks.
Note that for higher masses at solar metallicity, some combination of Pextra_factor, mass-loss, and/or superadiabatic convection reduction (e.g. mlt++)
might be necessary to stabilize the surface and avoid numerical issues. See the 80Msun_zams_to_cc test_suite as an example.

For production science we recommend adopting tighter mesh and timestep controls, such as those suggested in the comments of inlist_common.

Physical checks
===============

None

Inlists
=======

This test case has seven parts.

* Part 1 (``inlist_make_late_pre_zams``) creates a 20 |Msun|, Z=1.42*10^-2 metallicity, pre-main sequence model and evolves it for 100 years.

* Part 2 (``inlist_to_zams``) evolves the model to the zero age main sequence.

* Part 3 (``inlist_to_end_core_he_burn``) takes the model to core helium depletion.

* Part 4 (``inlist_remove_envelope``) removes the remaining hydrogen envelope. (optional)

* Part 5 (``inlist_to_end_core_c_burn``) takes the model to core carbon depletion.

* Part 6 (``inlist_to_lgTmax``) evolves the model until the core temperature reaches log T =9.60 (approximately silicon-shell burning)

* Part 7 (``inlist_to_cc``) evolves until core collapse.


Last-Updated: 18Dec2023 by EbF
