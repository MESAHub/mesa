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
creating a pre-main-sequence model of 100 Msun with Z=0.02, and then it
"relaxes" Z down to 1e-5 and the mass up to 110 Msun before starting the
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
a 1.3 Msun, metal-rich Z=0.04 model from the pre-main sequence to core hydrogen depletion.

:ref:`1.4M_ms_op_mono`
^^^^^^^^^^^^^^^^^^^^^^

The test checks the functionality of OP mono opacities. 
The test vehicle is a 1.4 Msun solar metallicity model.

:ref:`1.5M_with_diffusion`
^^^^^^^^^^^^^^^^^^^^^^^^^^

The test checks the functionality of element diffusion.
The test vehicle is a 1.5 Msun solar metallicity model.

:ref:`15M_dynamo`
^^^^^^^^^^^^^^^^^

The test checks the functionality of element rotation in a 15 Msun solar metallicity model.

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

This test case checks the evolution of a 1 Msun, Z=0.02 metallicity from the pre-main sequence to a white dwarf.

:ref:`1M_thermohaline`
^^^^^^^^^^^^^^^^^^^^^^

The test checks thermohaline mixing in a rotating, 1 Msun, Z=0.02 metallicity model.


:ref:`20M_z2m2_high_rotation`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This test case checks the evolution of a strongly rotating,
Omega/Omega_crit = 0.75, 20 Msun, Z=0.02 metallicity model from the
pre-main sequence to the end of core helium burning.


:ref:`5M_cepheid_blue_loop`
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This test case checks that the evolution of a 5 Msun, metal-poor Z = 0.008, helium-enriched Y=0.256 model
executes a blue-loop in the HR diagram and crosses the classical Cepheid instability strip boundaries three times.

:ref:`7M_prems_to_AGB`
^^^^^^^^^^^^^^^^^^^^^^

This test case checks that the evolution of a 7 Msun, metal-poor Z = 0.001, model reaches the AGB.


:ref:`accreted_material_j`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This test suite example checks the accretion of material and angular momentum onto a 20 Msun model.

:ref:`adjust_net`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This test suite example checks the functionality of the adaptive nuclear reaction network.


:ref:`carbon_acc`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This test suite example checks the acceration of carbon onto a pure carbon white dwarf.


:ref:`c13_pocket`
^^^^^^^^^^^^^^^^^

This test evolves a 2.0 |Msun| star through one thermal pulse on the
asymptotic giant branch (AGB) and illustrates third dredge up and the
formation of a :math:`^{13}{\rm C}` pocket.

      
:ref:`conductive_flame`
^^^^^^^^^^^^^^^^^^^^^^^

This test case models a conductively-propagated deflagration wave
("flame") in a high-density, degenerate carbon-oxygen mixture.  It
also provides an example for use of the ``other_build_initial_model``
and ``other_surface_PT`` hooks.


:ref:`hb_2M`
^^^^^^^^^^^^

This test case shows a 2 |Msun| stellar model evolving
on the horizontal branch (HB) through core helium burning.

:ref:`make_co_wd`
^^^^^^^^^^^^^^^^^

This test case produces a 0.6 |Msun| white dwarf with a carbon-oxygen
dominated core and a stratified atmosphere dominated by hydrogen at
its surface. The final model produced by this test case also serves as
the starting model for :ref:`wd_diffusion` and :ref:`wd_cool_0.6M`.

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

:ref:`T_tau_gradr`
^^^^^^^^^^^^^^^^^^

This test checks the implementation of the control
``use_T_tau_gradr_factor``, which modifies the radiative gradient so
that regions of low optical depth have a temperature that follows the
:math:`T(\tau)` relation specified by ``atm_T_tau_relation``.

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

surface_effects
^^^^^^^^^^^^^^^

Tests the implementation of the various surface effect corrections
available in MESA.


Test Index
----------

This index only includes tests that are documented via a ``README.rst``.

.. toctree::
   :maxdepth: 1
   :glob:

   test_suite/*
