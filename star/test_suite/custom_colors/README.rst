.. _custom_colors:

*************
custom_colors
*************

This test suite case demonstrates the calculation of synthetic photometry during stellar evolution using the MESA colors module. It evolves a 1 |Msun|, Z=0.02 metallicity stellar model from the pre-main sequence while computing bolometric magnitudes and synthetic magnitudes in multiple photometric filters, adding these as extra history columns.

The colors module, developed by Niall Miller and Meridith Joyce, enables synthetic photometry by interpolating stellar atmosphere models to construct spectral energy distributions (SEDs) matching the evolving stellar parameters (effective temperature \( T_{\eff} \), surface gravity \( \log g \), metallicity [M/H]). These SEDs are convolved with filter transmission curves to produce absolute magnitudes.

This module is introduced in MESA r25.10.1-rc1 (release candidate available at https://zenodo.org/records/17426065, with full release expected in November/December 2025). For details, see the module documentation at https://docs.mesastar.org/en/latest/colors/overview.html.

This test case has one part.

* Part 1 (``inlist``) creates a 1 |Msun|, Z=0.02 metallicity pre-main sequence model and evolves it, calculating synthetic photometry at each timestep using interpolated stellar atmosphere models. The test adds the following columns to the history file:

.. code-block:: console

   Mag_bol              ! Bolometric magnitude
   Flux_bol             ! Bolometric flux
   Gbp                  ! Gaia blue photometer magnitude
   G                    ! Gaia magnitude
   Grp                  ! Gaia red photometer magnitude

At the end of the run, the test reports synthetic magnitudes within expected ranges for the final stellar parameters.

Physical Setup
==============

The colors module interpolates pre-computed stellar atmosphere models to generate SEDs. These are integrated for bolometric quantities and convolved with filters for synthetic magnitudes.

- **Stellar Atmosphere Models**: Kurucz 2003 grid covering \( T_{\eff} = 3500-50000 \) K, \( \log g = 0.0-5.0 \), [M/H] = -5.0 to +1.0.
- **Filter System**: Gaia DR3 photometric bands (Gbp, G, Grp) by default.
- **Reference Spectrum**: Vega SED for magnitude zero-points.
- **Distance**: 10 parsecs for absolute magnitudes.

Data Download
=============

The test automatically downloads required data files (~35 MB) on first run:

.. code-block:: console

   ./mk      # Downloads atmosphere models, filters, and Vega spectrum
   ./rn      # Runs the test

The build script (``mk``) creates the necessary directory structure and downloads:

- Stellar atmosphere model grid and lookup table
- Filter transmission curves
- Vega reference spectrum

For additional atmospheres and filters from sources like SVO, MSG, MAST, use the SED_Tools repository at https://github.com/nialljmiller/SED_Tools, which provides a CLI for downloading and preparing them for use in the custom colors module.

Verification Tools
==================

Python helper scripts are provided for monitoring and verification:

.. code-block:: console

- python HISTORY_check.py      # Real-time monitoring of photometric evolution
- python SED_check.py          # Real-time SED visualization

Expected Outputs
================

The test produces standard MESA output plus additional photometric columns in the history file. 
- The history file will have "Mag_bol", "Flux_bol".
- A column for every filter in the user-selected instrument.
- If make_csv is .true., it will also output an SED for each filter band as a CSV file in colors_results/{filter}_SED.csv.


Last-Updated: 05Nov2025 by Niall Miller