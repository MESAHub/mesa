.. _make_co_wd:

**********
make_co_wd
**********

This test builds a carbon-oxygen white dwarf starting from a specified
pre main-sequence mass. It does this by running through a series of
inlists for different stages of evolution. By default, the first two
steps (pre main-sequence and main-sequence through core helium
burning) are skipped so that the test will run in a reasonable amount
of time, and the test produces a 0.6 |Msun| white dwarf from a roughly
3 |Msun| progenitor.

Starting from the beginning to produce a white dwarf descended from a
progenitor of different mass or metallicity can be accomplished by
setting the environment variable ``MESA_RUN_OPTIONAL=t`` and editing
the values assigned to ``initial_mass`` or ``initial_z`` in
``inlist_common``.

Larger values for ``initial_mass`` can be used to produce more massive
carbon-oxygen white dwarfs, but note that this test produces only an
approximate initial-to-final-mass relation (IFMR). The test is
designed to truncate the AGB evolution at the beginning of the
thermally-pulsing AGB (TP-AGB) stage, which can be difficult and/or
resource-intensive to model.

The five steps in this procedure are:

1. ``inlist_zams`` (optional)
-----------------------------

This step starts from the pre main-sequence to produce a ZAMS model
with the initial mass and metallicity specified in ``inlist_common``.


2. ``inlist_to_end_he_core_burn`` (optional)
--------------------------------------------

This step evolves from ZAMS through the end of core helium burning.


3. ``inlist_co_core``
---------------------

This step evolves to near the end of shell helium burning, producing
the final C/O core profile in the interior. This step terminates when
a thermal pulse occurs in the helium shell.


4. ``inlist_remove_env``
------------------------

This step removes the AGB envelope using the ``star_relax_to_star_cut``
method in ``src/run_star_extras.f90``. The location of the cut is
specified using ``x_ctrl(1)`` to leave :math:`10^{-3}` |Msun| of the
hydrogen envelope. After the relaxation ends, residual burning will
reduce the final hydrogen envelope mass to around
:math:`10^{-4}` |Msun| (for the default 0.6 |Msun| white dwarf).
      

5. ``inlist_settle``
--------------------

This step turns on diffusion in the young proto-WD model to allow the
model to settle into a stratified envelope structure with a pure
hydrogen atmosphere. This final step ends when the white dwarf has
cooled down to reach a luminosity of 1 |Lsun|.

