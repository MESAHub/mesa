.. highlight:: console

**********
Test suite
**********

MESA includes a comprehensive test suite.

Building upon test suite cases
------------------------------

Your first stop when setting up a new problem with MESA should be the
MESA test suite. You will find a wide range of sample cases there.
Looking at the test_suite inlists is a quick way to familiarize yourself
with the set of options relevant to your problem. You may want to copy
an inlist from the test suite to one of your working directories to use
as a starting point for a project of your own.

Each test suite problem lives in a subdirectory of

::

   $MESA_DIR/star/test_suite

and you can find (slightly out-of-date, but still useful) descriptions
of some of the test problems in the ``docs/`` sub-directory of each
test_suite case.

For example, take a look at the "high mass" test case. It starts by
creating a pre-main-sequence model of 100 |Msun| with Z=0.02, and then it
"relaxes" Z down to 1e-5 and the mass up to 110 |Msun| before starting the
evolution. It will take under 200 steps (and a few minutes) to reach a
central X of 0.5. To try it yourself,

::

   cd star/test_suite/high_mass
   ./mk
   ./rn

You can do the same with any of the test_suite cases.

If you want to base your work off of a test_suite case, you should make
a copy the directory and then edit this copy.

::

   cp -r $MESA_DIR/star/test_suite/high_mass my_high_mass

The test_suite examples require a few tweaks in order to be used
"outside" of the test_suite directory. First, you need to edit
make/makefile and delete the line

.. note::
   This is no longer needed and is left only as a reference for previous versions of MESA.
   All instances of MESA_DIR have been removed from all test cases.
   You may still need to adjust some inlist paths if they specify relative paths instead of
   absolute paths.

::

   MESA_DIR = ../../../..

Then, edit the inlist files and delete the line

::

   mesa_dir = '../../..'

Then edit the ``rn`` script and delete the line

::

   MESA_DIR=../../..

These changes ensure that you are using the copy of MESA specified by
the ``$MESA_DIR`` environment variable.

You might also need to adjust filenames of any initial models or other
inlists, if they are specified by a relative path.  (You can simply
change these to be an absolute path.)

The test_suite inlists specify rather strict limits on the number of
steps (by setting ``max_model_number``) and/or retries (by setting
``max_number_retries``). You likely want to delete these limits.

The MESA test_suite problems also have non-standard run_star_extras,
including routines that check the runtime of the example. If these annoy
you, they can be pruned by hand.

Tools such as Bill Wolf's
`mesa-cli <http://wmwolf.github.io/mesa_cli/>`__ can automate some of
these steps.


Star tests
----------

:ref:`1.3M_ms_high_Z`
^^^^^^^^^^^^^^^^^^^^^

The test checks the evolution of metal-rich low-mass stars by evolving 
a 1.3 |Msun|, metal-rich Z=0.04 model from the pre-main sequence to core hydrogen depletion.

:ref:`1.4M_ms_op_mono`
^^^^^^^^^^^^^^^^^^^^^^

The test checks the functionality of OP mono opacities. 
The test vehicle is a 1.4 |Msun| solar metallicity model.

:ref:`1.5M_with_diffusion`
^^^^^^^^^^^^^^^^^^^^^^^^^^

The test checks the functionality of element diffusion.
The test vehicle is a 1.5 |Msun| solar metallicity model.

:ref:`15M_dynamo`
^^^^^^^^^^^^^^^^^

The test checks the functionality of element rotation in a 15 |Msun| solar metallicity model.

:ref:`16M_conv_premix`
^^^^^^^^^^^^^^^^^^^^^^

This test suite example re-creates the 16-solar mass main-sequence
evolution with the inclusion of convective premixing (using the Ledoux
criterion), as detailed in Section 5.3 of the MESA V instrument paper
(Paxton et al 2019).

:ref:`16M_predictive_mix`
^^^^^^^^^^^^^^^^^^^^^^^^^

This test suite example re-creates the 16-solar mass main-sequence
evolution with the inclusion of predictive mixing (using the Ledoux
criterion), as detailed in Section 2 of the MESA IV instrument paper
(Paxton et al 2018).

:ref:`1M_pre_ms_to_wd`
^^^^^^^^^^^^^^^^^^^^^^

This test case checks the evolution of a 1 |Msun|, Z=0.02 metallicity from the pre-main sequence to a white dwarf.

:ref:`1M_thermohaline`
^^^^^^^^^^^^^^^^^^^^^^

The test checks thermohaline mixing in a rotating, 1 |Msun|, Z=0.02 metallicity model.

:ref:`12M_pre_ms_to_core_collapse`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This test suite evolves a 12 |MSun| model from the pre-ms to core collapse.

:ref:`20M_pre_ms_to_core_collapse`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This test suite evolves a low metalicity 20 |MSun| model from the pre-ms to core collapse.


:ref:`20M_z2m2_high_rotation`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This test case checks the evolution of a strongly rotating,
Omega/Omega_crit = 0.75, 20 |Msun|, Z=0.02 metallicity model from the
pre-main sequence to the end of core helium burning.


:ref:`5M_cepheid_blue_loop`
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This test case checks that the evolution of a 5 |Msun|, metal-poor Z = 0.008, helium-enriched Y=0.256 model
executes a blue-loop in the HR diagram and crosses the classical Cepheid instability strip boundaries three times.

:ref:`7M_prems_to_AGB`
^^^^^^^^^^^^^^^^^^^^^^

This test case checks that the evolution of a 7 |Msun|, metal-poor Z = 0.001, model reaches the AGB.


:ref:`accreted_material_j`
^^^^^^^^^^^^^^^^^^^^^^^^^^

This test suite example checks the accretion of material and angular momentum onto a 20 |Msun| model.

:ref:`adjust_net`
^^^^^^^^^^^^^^^^^

This test suite example checks the functionality of the adaptive nuclear reaction network.

:ref:`c13_pocket`
^^^^^^^^^^^^^^^^^

This test evolves a 2.0 |Msun| star through one thermal pulse on the
asymptotic giant branch (AGB) and illustrates third dredge up and the
formation of a :math:`^{13}{\rm C}` pocket.


:ref:`carbon_kh`
^^^^^^^^^^^^^^^^

This test suite case evolves a stellar model with a pure carbon
composition as it Kelvin-Helmholtz contracts.  It provides a
convergence example for the different forms of the energy equation.


:ref:`cburn_inward`
^^^^^^^^^^^^^^^^^^^

This test suite example checks the inward propagation of a carbon burning front in a 7.5 |Msun| model.


:ref:`ccsn_IIp`
^^^^^^^^^^^^^^^

This test suite example builds a Type IIp supernova model, including Rayleigh-Taylor Instability mixing, for subsquent use in STELLA.


:ref:`check_pulse_atm`
^^^^^^^^^^^^^^^^^^^^^^

This test checks that the atmosphere structure written to the
pulsation output closely matches what is expected for the
:math:`T(\tau)` relation specified by ``atm_T_tau_relation``.

   
:ref:`conductive_flame`
^^^^^^^^^^^^^^^^^^^^^^^

This test case models a conductively-propagated deflagration wave
("flame") in a high-density, degenerate carbon-oxygen mixture.  It
also provides an example for use of the ``other_build_initial_model``
and ``other_surface_PT`` hooks.


:ref:`conserve_angular_momentum`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This test suite example checks angular momentum conservation from the zero age main-sequence to the formation of a helium core in 1.0 |Msun|, Z=0.02 metallicity, model.


:ref:`conv_core_cpm`
^^^^^^^^^^^^^^^^^^^^

This test case evolves a 1.5 |Msun| star part of the way through
the main sequence with CPM enabled and checks that its convective
core has grown to an appropriate mass coordinate.


:ref:`custom_colors`
^^^^^^^^^^^^^^^^^^^^

This test suite example shows how to use user-defined color filter and extinction files.

:ref:`custom_rates`
^^^^^^^^^^^^^^^^^^^

This test suite case checks the use of custom nuclear reaction rates in an accreting 0.3 |Msun| helium white dwarf model.


:ref:`diffusion_smoothness`
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This test suite case checks that element diffusion produces a sufficiently smooth Brunt profile.

:ref:`extended_convective_penetration`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This test case checks the implementation of the extended convective penetration prescription for core boundary mixing.


:ref:`gyre_in_mesa_bcep`
^^^^^^^^^^^^^^^^^^^^^^^^

This test case checks the implementation of GYRE in MESA for a 12 |Msun|, Z=0.02 metallicity, model evolving from the zero-age main sequence to core hydrogen depletion;
a beta Cephei stellar model.


:ref:`gyre_in_mesa_envelope`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This test case checks the implementation of GYRE in MESA for the envelope of a 12 |Msun|, Z=0.02 metallicity, model.

:ref:`gyre_in_mesa_ms`
^^^^^^^^^^^^^^^^^^^^^^

This test case checks the implementation of GYRE in MESA for a 1 |Msun|, Z=0.02 metallicity, model evolving from the zero-age main sequence to core hydrogen depletion.

:ref:`gyre_in_mesa_rsg`
^^^^^^^^^^^^^^^^^^^^^^^

This test case checks the implementation of GYRE in MESA for a 21 |Msun|, Z=0.02 metallicity, model in the red supergiant regime.

:ref:`gyre_in_mesa_spb`
^^^^^^^^^^^^^^^^^^^^^^^

This test case checks the implementation of GYRE in MESA for a 5 |Msun|, Z=0.02 metallicity, model evolving from the zero-age main sequence to core hydrogen depletion; 
a slowly pulsating B-type star (SPB) stellar model.

:ref:`gyre_in_mesa_wd`
^^^^^^^^^^^^^^^^^^^^^^

This test case checks the implementation of GYRE in MESA for a cooling 0.85 |Msun| white dwarf model.

:ref:`hb_2M`
^^^^^^^^^^^^

This test case shows a 2 |Msun| stellar model evolving
on the horizontal branch (HB) through core helium burning.

:ref:`high_mass`
^^^^^^^^^^^^^^^^

This test case checks the evolution of a 300 |Msun|, Z = 1e-4 metallicity, model through core hydrogen depletion.

:ref:`high_z`
^^^^^^^^^^^^^

This test case checks the capability of evolving high metallicity models through core helium depletion with a 7 |Msun|, Z=0.07 metallicity model.

:ref:`hot_cool_wind`
^^^^^^^^^^^^^^^^^^^^

This test case checks the cool wind, hot wind capability by evolving a 7 |Msun|, Z=0.02 metallicity model from the zero-age main sequence to core helium depletion.

:ref:`hse_riemann`
^^^^^^^^^^^^^^^^^^

This test case checks Riemann HLLC solver can hold an envelope model in hydrostatic equilibrium.

:ref:`irradiated_planet`
^^^^^^^^^^^^^^^^^^^^^^^^

This test case checks the evolution of an ~1 Mjup model after the surface has been irradiated.

:ref:`low_z`
^^^^^^^^^^^^

This test case checks the evolutions of a 0.8 |Msun|, Z=1e-4 metallicity model from the pre-main sequence to core hydrogen depletion.


:ref:`magnetic_braking`
^^^^^^^^^^^^^^^^^^^^^^^

This test case involves the calculation of the spin down caused by a
large-scale magnetic field in a massive star model.


:ref:`make_brown_dwarf`
^^^^^^^^^^^^^^^^^^^^^^^

This test case checks the creation of a 1.05 Mjup, Z=1e-4 metallicity model and its subsequent evolution for 20 billion years.


:ref:`make_co_wd`
^^^^^^^^^^^^^^^^^

This test case produces a 0.6 |Msun| white dwarf with a carbon-oxygen
dominated core and a stratified atmosphere dominated by hydrogen at
its surface. The final model produced by this test case also serves as
the starting model for :ref:`wd_diffusion` and :ref:`wd_cool_0.6M`.


:ref:`make_env`
^^^^^^^^^^^^^^^

This test case checks the creation and stability of a pure iron neutron star envelope.


:ref:`make_he_wd`
^^^^^^^^^^^^^^^^^

This test case checks the creation and evolution of a 0.15 |Msun| helium white dwarf.

:ref:`make_metals`
^^^^^^^^^^^^^^^^^^

This test case demonstrates the creation and evolution of 3 |Msun| model whose initial metallicity is Z = 0.


:ref:`make_o_ne_wd`
^^^^^^^^^^^^^^^^^^^

This test case produces a 1.05 |Msun| oxygen-neon-magnesium white dwarf using stellar engineering.

:ref:`make_planets`
^^^^^^^^^^^^^^^^^^^

This test case shows an example of a 1 Mjup model with a 10 Mearth core that is irradiated and evolved for 10Gyr.

:ref:`make_sdb`
^^^^^^^^^^^^^^^

This test case shows an example of making a 0.4 |Msun|, Z=0.02 metallicity, helium model - a B-type subdwarf (sdB) star.

:ref:`make_zams`
^^^^^^^^^^^^^^^^

This test case shows an example of creating a 4 |Msun|, Z = 0.01 metallicity, pre-main sequence model and evolving it to the zero age main sequence.

:ref:`make_zams_low_mass`
^^^^^^^^^^^^^^^^^^^^^^^^^

This test case shows an example of creating a 0.085 |Msun|, Z = 0.014 metallicity, pre-main sequence model and evolving it to the zero age main sequence.

:ref:`make_zams_ultra_high_mass`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This test case shows an example of creating a 250 |Msun|, Z = 1e-4 metallicity, model close to the main sequence.


:ref:`ns_h`
^^^^^^^^^^^

This test case shows an example of steady hydrogen burning within a neutron star envelope.

:ref:`ns_he`
^^^^^^^^^^^^

This test case shows an example of a helium flash within a neutron star envelope.

:ref:`ns_c`
^^^^^^^^^^^

This test case shows an example of a carbon flash within a neutron star envelope.

:ref:`other_physics_hooks`
^^^^^^^^^^^^^^^^^^^^^^^^^^

This test case exercises several of the ``other_*`` physics hooks simultaneously in a 1 |Msun|, Z=0.02 metallicity, model.
It provides an example of how to include your own physics code into a MESA run.


:ref:`pisn`
^^^^^^^^^^^^^^^^^^^^^^

This test case evolves an initially 200 |Msun| star from ZAMS untill it undergoes a pair instability supernovae (PISN).


:ref:`ppisn`
^^^^^^^^^^^^

This test case shows an example of a star undergoing a pulsational
pair-instability supernova. The model starts from a massive helium
star, and includes switches from hydrostatic to hydrodynamic models,
as well as the removal of ejected layers.


:ref:`R_CrB_star`
^^^^^^^^^^^^^^^^^

This test case creates and evolves a simple model of an R Corona
Borealis star and provides an example of how to use AESOPUS opacity
tables in MESA.


:ref:`radiative_levitation`
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This test case exercises radiative levitation and the OP mono opacities in the outer layers of a 0.466 |Msun|, Z=0.02 metallicity, B-type subdwarf (sdB) model.


:ref:`relax_composition_j_entropy`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This test calls the routines that relax the composition, angular
momentum and energy of a model to given target values.

:ref:`rsp_BEP`
^^^^^^^^^^^^^^

This test case checks the non-linear pulsation evolution of a 0.26 |Msun|, Teff = 6968 K, L = 33 Lsun, Z = 0.01 metallicity model - a binary evolution pulsator similar
the one shown in Smolec et al 2013, MNRAS.

:ref:`rsp_BLAP`
^^^^^^^^^^^^^^^

This test case checks the non-linear pulsation evolution of a 0.36 |Msun|, Teff = 26,000 K, L = 320 Lsun, Z = 0.05 metallicity model -
a blue large-amplitude pulsator model originally contributed by Alfred Gautschy.

:ref:`rsp_Cepheid`
^^^^^^^^^^^^^^^^^^

This test case checks the non-linear pulsation evolution of a 4.165 |Msun|, Teff = 6050 K, L = 1438.8 Lsun, Z = 0.007 metallicity model -
a classical Cepheid variable similar to CEP-227 shown in |Pilecki2013|.


:ref:`rsp_Delta_Scuti`
^^^^^^^^^^^^^^^^^^^^^^

This test case checks the non-linear pulsation evolution of a 2 |Msun|, Teff = 6900 K, L = 30 Lsun, Z = 0.02 metallicity -
a double-mode delta Scuti variable leaving the main-sequence phase originally contributed by Alfred Gautschy.

:ref:`rsp_RR_Lyrae`
^^^^^^^^^^^^^^^^^^^

This test case checks the non-linear pulsation evolution of a 0.65 |Msun|, Teff = 6500 K, L = 60 Lsun, Z = 0.004 metallicity -
a long-period RR Lyrae model contributed by Radek Smolec.

:ref:`rsp_Type_II_Cepheid`
^^^^^^^^^^^^^^^^^^^^^^^^^^

This test case checks the non-linear pulsation evolution of a 0.55 |Msun|, Teff = 6410 K, L = 136 Lsun, Z = 0.0001 metallicity model -
type-II Cepheid of BL Her type based on |Smolec14|.

:ref:`rsp_check_2nd_crossing`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This test case exercises the RSP model building and linear nonadiabatic stability analysis
to find the instability strip edges, and effective temperatures offset from the blue edge of the instability strip.

:ref:`rsp_save_and_load_file`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This test case checks that RSP models can be saved and loaded to produce the same results as test case :ref:`rsp_Cepheid`.

:ref:`semiconvection`
^^^^^^^^^^^^^^^^^^^^^

This test case checks placement of the convective and semiconvective boundaries when using the Ledoux criterion and predictive mixing,
see |MESA V|.The test vehicle is with a 1.5 |Msun|, Z=0.02 metallicity, model.

:ref:`simplex_solar_calibration`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This test case exercises the simplex framework with a check of the chi^2 value for 1.0 |Msun|, Z=0.02 metallicity, solar model.

:ref:`split_burn_big_net`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This test case tests MESA's ability to perfom a split-burn calculation in a 25 |Msun| star during silicon burning.

:ref:`starspots`
^^^^^^^^^^^^^^^^

This test case implements modifications to the surface structure of a 1 solar mass star based on a star spots formalism.

:ref:`T_tau_gradr`
^^^^^^^^^^^^^^^^^^

This test checks the implementation of the control
``use_T_tau_gradr_factor``, which modifies the radiative gradient so
that regions of low optical depth have a temperature that follows the
:math:`T(\tau)` relation specified by ``atm_T_tau_relation``.


:ref:`test_memory`
^^^^^^^^^^^^^^^^^^

This test case program checks MESA's memory management.
It is designed primarily to be run inside the valgrind leak-checking tool,
and is based on code provided originally by Warrick Ball.

:ref:`timing`
^^^^^^^^^^^^^

This test checks the counter and timing routines with a 1.5 |Msun|, Z=0.02 metallicity model.

:ref:`twin_studies`
^^^^^^^^^^^^^^^^^^^

This test case exercise the capability to simultaneously evolve two model stars.
The test vehicle is a pair of 15 |Msun|, Z=0.02 metallicity, models one with overshooting and one without overshooting.


:ref:`tzo`
^^^^^^^^^^

This test case creates a Thorne-Zytkow object (TZO) and evolves until the NS has accreted a small amount of material.

:ref:`wd_acc_small_dm`
^^^^^^^^^^^^^^^^^^^^^^

This test case models an accreting CO white dwarf (WD) and checks that
the composition of the accreted material is being correctly tracked.

:ref:`wd_aic`
^^^^^^^^^^^^^

This test case shows an accreting ONeMg white dwarf (WD) evolving
towards accretion induced collapse (AIC).  It also illustrates use of
the special weak rate implementation described in Section 8 of |MESA
III|.
   

:ref:`wd_c_core_ignition`
^^^^^^^^^^^^^^^^^^^^^^^^^

This test case the checks the onset of a thermonuclear runaway in an accreting Chandrasekhar mass carobon-oxygen white dwarf.

:ref:`wd_cool_0.6M`
^^^^^^^^^^^^^^^^^^^

This test case the checks the evolution of a cooling, element diffusing 0.6 |Msun| white dwarf.

:ref:`wd_diffusion`
^^^^^^^^^^^^^^^^^^^

This test case the checks element diffusion in a 0.6 |Msun| carbon-oxygen white dwarf.

:ref:`wd_he_shell_ignition`
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This test case the ignition of a helium layer in an accreting in a 0.96 |Msun| carbon-oxygen white dwarf model.

:ref:`wd_nova_burst`
^^^^^^^^^^^^^^^^^^^^

This test case checks the evolution of a nova outburst for one cycle.

:ref:`wd_stable_h_burn`
^^^^^^^^^^^^^^^^^^^^^^^

This test case checks the evolution stable hydrogen burning on a white dwarf.



Binary tests
------------

double_bh
^^^^^^^^^

Creates a binary black hole from two stars in a very close orbit
through the chemically-homogeneous evolution (CHE) mechanism.  Stars
evolve through overcontact phases, so they test the overcontact
prescription.

evolve_both_stars
^^^^^^^^^^^^^^^^^

Tests MESA evolving two stars simultaneously including mass transfer.

jdot_ml_check
^^^^^^^^^^^^^

Using pre-specified efficiency options, verifies that the evolution
follows the analytical result from `Tauris & van den Heuvel (2006)
<https://ui.adsabs.harvard.edu/abs/2006csxs.book..623T>`_. Shuts off
all other ``jdot`` sources.

jdot_gr_check
^^^^^^^^^^^^^

With all other ``jdot`` sources turned off, this verifies that the
orbital evolution due to GW emission follows the analytical result of
`Peters (1964) <https://ui.adsabs.harvard.edu/abs/1964PhRv..136.1224P>`_.

jdot_ls_check
^^^^^^^^^^^^^

Verifies that models with tidal evolution conserve angular momentum.

star_plus_point_mass
^^^^^^^^^^^^^^^^^^^^

Tests MESA evolving one star plus a point mass, including mass
transfer to the point mass.

star_plus_point_mass_explicit_mdot
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Same as above, but run using an explicit calculation for the mass
transfer rate.

wind_fed_hmxb
^^^^^^^^^^^^^

Model for a high mass X-ray binary, including both Roche lobe overflow
and wind mass transfer. Verifies the Eddington limit is working, and
that the accretion luminosity is computed correctly.

Astero tests
------------

astero_adipls and astero_gyre
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Demonstrate how to use the ``astero`` module to call ADIPLS or GYRE, respectively.

Both tests use the same main-sequence evolution of a 1.2 |Msun| star,
so the evolutionary outputs (e.g. ``final.mod``, ``LOGS/history.data``)
should be identical for both tests.  Tables of
mode frequencies are displayed in the terminal every 50 models and
should be the same to within about 0.1 μHz.  These can be compared by
``diff``\ ing the terminal output.

The tests compare one mode frequency (currently for (*ℓ*, *n*)=(0,4))
with a target value set by ``x_ctrl(1)`` with a tolerance of 3%.  Over time, cumulative
changes to microphysics might mean the targets are missed, causing
failure.  In this case, both targets can be updated to the same value.
Failure should be investigated if, e.g.,

* something is changed that shouldn't affect the main-sequence evolution of a 1.2 |Msun| star or
* one of the tests passes and the other fails.

Note that
GYRE can also be called directly, without using the ``astero`` module.
See the ``gyre_in_mesa_*`` test cases in ``star``'s test suite.

:ref:`example_astero`
^^^^^^^^^^^^^^^^^^^^^

An example optimisation run of the ``astero`` module, based on the CoRoT
target HD 49385.  This is the usual starting point if you want to
optimise model parameters using the ``astero`` module.

:ref:`fast_from_file`
^^^^^^^^^^^^^^^^^^^^^
:ref:`fast_newuoa`
^^^^^^^^^^^^^^^^^^
:ref:`fast_scan_grid`
^^^^^^^^^^^^^^^^^^^^^
:ref:`fast_simplex`
^^^^^^^^^^^^^^^^^^^

Each of these test cases runs a handful of iterations of a crude
optimisation, principally to increase test coverage across the
``astero`` module.

:ref:`surface_effects`
^^^^^^^^^^^^^^^^^^^^^^

Tests the implementation of the various surface effect corrections
available in MESA.


.. |Smolec14| replace:: `Smolec and Moskalik (2014) <https://ui.adsabs.harvard.edu/abs/2014MNRAS.441..101S/abstract>`__
.. |Pilecki2013| replace:: `Pilecki et al (2013) <https://ui.adsabs.harvard.edu/abs/2013MNRAS.436..953P/abstract>`__


Test Index
----------

This index only includes tests that are documented via a ``README.rst``.

.. toctree::
   :maxdepth: 1
   :glob:

   test_suite/*
