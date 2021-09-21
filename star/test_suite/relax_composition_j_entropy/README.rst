.. _relax_composition_j_entropy:


***************************
relax_composition_j_entropy
***************************

This test calls the routines that relax the composition, angular
momentum and energy of a model to given target values.

Part 1 (``inlist_create_input``) evolves a 10 |Msun| model with a 300
km/s surface rotation rate from ZAMS until the effective temperature
drops below 8000 K and saves the composition, angular momentum and
entropy profiles to files.

Part 2 (``inlist_start_relax_composition_j_entropy``) uses MESA's
``relax_initial_...`` controls to relax the same ZAMS model as part 1
directly to the target profiles, also saved in part 1, and
then checks that the luminosity and effective temperature are within
0.5 per cent of the values in the final model of part 1.

Part 3 (``inlist_relax_composition_j_entropy``) evolves the relaxed
model from part 2 until it finishes core helium burning.

Last-Updated: 2021-09-16 (commit 77c4f41e) by Warrick Ball
