.. _custom_colors:

*************
custom_colors
*************

This test case demonstrates the calculation of synthetic photometry during stellar evolution using the MESA colors module. It evolves a stellar model while computing bolometric magnitudes and synthetic magnitudes in multiple photometric filters, adding these as extra history columns.

This test case has one part. Click to see a larger view of a plot.

* Part 1 (``inlist``) evolves a 1 |Msun|, Z=0.02 metallicity model from the pre-main sequence, calculating synthetic photometry at each timestep using interpolated stellar atmosphere models. The test adds the following columns to the history file:

.. code-block:: console

   Mag_bol              ! Bolometric magnitude
   Flux_bol             ! Bolometric flux
   Gbp                  ! Gaia blue photometer magnitude
   G                    ! Gaia magnitude
   Grp                  ! Gaia red photometer magnitude

At the end of the run, the test reports synthetic magnitudes that should be within expected ranges for the final stellar parameters.

Physical Setup
==============

The colors module interpolates between pre-computed stellar atmosphere models (Kurucz 2003) to construct spectral energy distributions (SEDs) matching the evolving stellar parameters (Teff, log g, metallicity). These SEDs are then convolved with filter transmission curves to produce synthetic magnitudes.

The test uses:

* **Stellar atmosphere models**: Kurucz 2003 grid covering Teff = 3500-50000 K, log g = 0.0-5.0, [M/H] = -5.0 to +1.0
* **Filter system**: Gaia DR3 photometric bands (Gbp, G, Grp)
* **Reference spectrum**: Vega SED for magnitude zero-points
* **Distance**: 10 parsecs (for absolute magnitudes)

Configuration
=============

The test requires a ``&colors`` namelist in addition to standard MESA controls:

.. literalinclude:: ../../../star/test_suite/custom_colors/inlist_colors
   :start-after: &colors
   :end-before: ! end of colors namelist

Data Download
=============

The test automatically downloads required data files (~35MB) on first run:

.. code-block:: console

   ./mk      # Downloads atmosphere models, filters, and Vega spectrum
   ./rn      # Runs the test

The build script (``mk``) creates the necessary directory structure and downloads:

* Stellar atmosphere model grid and lookup table
* Filter transmission curves
* Vega reference spectrum

Verification Tools
==================

Python helper scripts are provided for monitoring and verification:

.. code-block:: console

   cd python_helpers
   python HISTORY_check.py      # Real-time monitoring of photometric evolution
   python static_HISTORY_check.py   # Static analysis of complete tracks
   python SED_check.py          # SED visualization and validation

Expected Outputs
================

The test produces standard MESA output plus additional photometric columns in the history file. Successful completion should show:

* Smooth evolution of synthetic magnitudes consistent with stellar parameter changes
* Bolometric magnitudes within expected ranges for stellar mass and evolutionary phase
* Filter magnitudes that track effective temperature variations

The ``run_star_extras.f90`` module handles the integration between MESA's stellar evolution and the colors module calculations.

Last-Updated: 17Jul2025 by Niall Miller
