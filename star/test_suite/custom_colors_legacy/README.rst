.. _custom_colors:

**************************************
Legacy Colors (precomputed photometry)
**************************************

This test suite exercises the **legacy colors path** of the MESA ``colors``
module: precomputed photometry / bolometric-correction tables loaded through
``photo_precompute = .true.``.

This is **not** the standard runtime SED-convolution workflow. In this suite,
magnitudes come from tabulated bolometric corrections, not from interpolated
stellar atmosphere spectra convolved with filter transmission curves.

What this test demonstrates
===========================

This case is meant to show that the legacy colors implementation still works
and still produces photometric history columns during a normal stellar
evolution run.

Specifically, it demonstrates that:

* the ``&colors`` namelist is configured for legacy mode via
  ``photo_precompute = .true.``;
* legacy BC tables are loaded from the filenames listed in
  ``color_file_names(:)``;
* the photometric history columns are populated from those BC tables; and
* the run bypasses the standard SED / filter-convolution path.

How this suite proves it is using legacy colors
===============================================

The legacy branch in ``colors_history.f90`` has a distinctive output
signature:

* ``photo_precompute = .true.`` in ``inlist_colors``;
* ``Mag_bol`` is computed from the stellar luminosity;
* ``Flux_bol`` is written as ``-1``;
* ``Interp_rad`` is written as ``-1``; and
* the filter magnitudes are computed as ``Mag_bol - BC`` using the loaded
  legacy bolometric-correction tables.

By contrast, the standard colors path:

* loads atmosphere grids and filters;
* computes a full SED;
* computes a physical ``Flux_bol``;
* computes a physical ``Interp_rad``; and
* evaluates synthetic magnitudes by convolving the SED with filter response
  curves.

So, for this suite, the simplest post-run proof of legacy mode is:

* ``photo_precompute = .true.`` in the inlist,
* finite photometric magnitudes in ``LOGS/history.data``, and
* ``Flux_bol = -1`` and ``Interp_rad = -1`` throughout the run.

A validation helper script is provided below for exactly that purpose.

Running the test suite
======================

This is a standard MESA work directory.

.. code-block:: bash

   ./clean
   ./mk
   ./rn

The run evolves a 7 Msun model from a pre-main-sequence start until the
configured central-helium termination criterion is reached.

After a successful run you should have at least:

* ``LOGS/history.data``
* ``final.mod``

If you add the validation helper to ``python_helpers/``, you can then verify
that the run really used legacy colors with:

.. code-block:: bash

   cd python_helpers
   python check_legacy_mode.py

Key controls used by this suite
===============================

The important settings in ``inlist_colors`` are:

.. code-block:: fortran

   &colors
      use_colors = .true.
      photo_precompute = .true.

      color_num_files = 2
      color_file_names(1) = 'lcb98cor.dat'
      color_num_colors(1) = 11
      color_file_names(2) = 'blackbody_johnson.dat'
      color_num_colors(2) = 5
   /

These settings select the legacy bolometric-correction-table path.

Where the legacy tables are loaded from
=======================================

The legacy loader first tries the filename exactly as supplied in
``color_file_names(:)``. If that fails, it falls back to:

.. code-block:: text

   $MESA_DIR/data/colors_data/photo_precompute/<filename>

So, for the default configuration in this suite, MESA looks for:

* ``lcb98cor.dat``
* ``blackbody_johnson.dat``

in the current working directory first, then in MESA's
``data/colors_data/photo_precompute/`` directory.

Important differences from standard colors
==========================================

This suite is intentionally testing the legacy path, so several standard
``colors`` controls are not the main story here.

In legacy mode as currently implemented:

* ``instrument`` is not used;
* ``stellar_atm`` is not used;
* ``vega_sed`` is not used;
* ``make_csv`` / ``sed_per_model`` are not part of the legacy photometry path;
* ``mag_system`` is not used by the legacy branch in ``colors_history.f90``;
* ``Flux_bol`` is not computed and is written as ``-1``;
* ``Interp_rad`` is not computed and is written as ``-1``; and
* ``distance`` is not applied in the legacy branch.

That last point matters: the current legacy implementation computes
``Mag_bol`` from the luminosity scaling alone,

.. code-block:: text

   Mag_bol = Mbol_sun - 2.5 log10(L / L_sun)

with no distance term. So this suite should be documented as a test of legacy
precomputed photometry, not as a test of apparent-magnitude generation.

Recommended validation script
=============================

A good minimal regression check for this suite is:

1. confirm ``photo_precompute = .true.`` in ``inlist_colors``;
2. confirm that ``Mag_bol`` exists and is finite;
3. confirm that ``Flux_bol`` exists and is always ``-1``;
4. confirm that ``Interp_rad`` exists and is always ``-1``; and
5. confirm that at least one photometric column after ``Interp_rad`` contains
   finite values.

This verifies the characteristic legacy output signature without depending on
any one specific filter set.

Notes and caveats
=================

The current ``inlist_colors`` file contains a distance comment that should be
corrected:

.. code-block:: fortran

   distance = 3.0857d19

This value is **10 parsecs in cm**, not 1 parsec.

Also note that the file currently enables:

.. code-block:: fortran

   pgstar_flag = .true.

That may be useful for interactive local runs, but for automated CI-style test
execution you may prefer ``.false.``.
