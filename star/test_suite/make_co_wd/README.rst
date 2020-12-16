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
carbon-oxygen white dwarfs, but note that this test is not designed to
produce a realistic initial-to-final-mass relation (IFMR). Instead,
this test is designed to robustly produce usable white dwarf models as
soon as evolution has yielded the desired core structure, and so this
procedure truncates the AGB evolution just before reaching the
thermally-pulsing AGB (TP-AGB) stage, which can be difficult and/or
resource-intensive to model. Skipping the TP-AGB stage means that the
procedure in this test case will generally require a somewhat larger
progenitor mass to produce a given white dwarf mass than what would be
expected from a realistic IFMR.

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
the He layer mass falls below 0.04 |Msun|, avoiding the TP-AGB phase
that would begin when the He layer becomes smaller. After the AGB
envelope is removed in the next step, residual burning in the He layer
will reduce its final mass to around 0.01 |Msun|.


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

