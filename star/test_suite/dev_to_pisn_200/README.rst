.. _dev_to_pisn_200:

***************
dev_to_pisn_200
***************

This test case evolves an initialy 200 |Msun| star from ZAMS untill it undergoes a pair instability supernovae (PISN)


Physical checks
===============

This tracks the central value of (Gamma1 - 4/3) and the point when the pressure weighted integral of Gamma1-4/3 (integral_gamma1) first drops below 0.
This value should be 0 or negative otherwise we do not have a PISN but instead PPISN (`Renzo et al (2020) <https://ui.adsabs.harvard.edu/abs/2020A%26A...640A..56R/abstract>`__)



Inlists
=======


This test case has seven parts.

* Part 1 (``inlist_make_late_pre_zams``) This creates a 200 |Msun| Z=1.6*10:sup:-3 model

* Part 2 (``inlist_to_zams_header``) This evolves the model to just before the MS starts

* Part 3 (``inlist_to_end_core_he_burn``) This takes the model up to the end of core helium burning

* Part 4 (``inlist_remove_envelope_header``) This removes whats left of the hydrogren envelope

* Part 5 (``inlist_to_end_core_c_burn``) This takes the model up to the end of carbon burning

* Part 6 (``inlist_convert``) This prepars the model for a PISN, by swithcing from cell-faced to cell-centered hydrodynamics, turning off rotation, and turning on the AMR mesh controls.

* Part 7 (``inlist_finish``) This evoles the model through the PISN by taking it through carbon burning, explosive oxygen ignition, and stops once the star reaches a positive total energy.



Last-Updated: 23Jun2021 (MESA 21a860) by rjfarmer

