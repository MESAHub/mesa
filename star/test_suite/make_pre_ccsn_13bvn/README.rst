.. _make_pre_ccsn_13bvn:

*******************
make_pre_ccsn_13bvn
*******************

This test suite evolves a solar metalicity 12 |MSun| model from the pre-ms to helium depletion, strips the hydrogen envelope, and then evolves to core collapse.
This test suite is made to be similar (almost identical) to the 13BVN model from MESA IV. (This is a sensitive test_suite)

Physical checks
===============

Tests the functionality of velocity drag for v_flag for damping large surfaces velocities from He shell (N-alpha) flashes during advanced burning. 
Tests energy conservation for regions where v_drag is applied.
Tests the use of OP_split_burn as split burn <4d9 leads to large temperature swings and single zone burning in the core, and split burn >4.5d9 needs smaller timesteps. 
This case illustrates strong coupling between burning and hydro: OP_split_burn lower and higher thresholds can only be exceeded with very tight timesteps.

Inlists
=======

This test case has five parts.

* Part 1 (``inlist_to_zams``) creates a 12 |Msun|, Z=2*10^-2 metallicity, pre-main sequence model and evolves it the zero age main sequence.

* Part 2 (``inlist_before_remove``) takes the model to core helium depletion.

* Part 3 (``inlist_remove``) removes the hydrogen envelope, leaving a stripped, evolved helium core.

* Part 4 (``inlist_after_remove``) replaces the remaining hydrogen in the envelope with helium.

* Part 5 (``inlist_to_pre_si_burn``) evolves the model until log_center_T = 9.45d0, just before Si-core burning.

* Part 6 (``inlist_to_cc``) evolves the model until core collapse (300 km/s infall).


Last-Updated: 24May2023 by EbF

