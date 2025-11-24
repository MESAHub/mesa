.. _custom_colors:

*************
custom_colors
*************

This test suite case demonstrates the functionality of the MESA ``colors`` module, a framework introduced in MESA r25.10.1 for calculating synthetic photometry and bolometric quantities during stellar evolution.

1. What is MESA colors?
=======================

MESA colors is a post-processing and runtime module that allows users to generate "observer-ready" data directly from stellar evolution models. Instead of limiting output to theoretical quantities like Luminosity (:math:`L`) and Surface Temperature (:math:`T_{\rm eff}`), the colors module computes:

* **Bolometric Magnitude** (:math:`M_{\rm bol}`)
* **Bolometric Flux** (:math:`F_{\rm bol}`)
* **Synthetic Magnitudes** in specific photometric filters (e.g., Johnson V, Gaia G, 2MASS J).

This bridges the gap between theoretical evolutionary tracks and observational color-magnitude diagrams (CMDs).

2. How does the MESA colors module work?
========================================

The module operates by coupling the stellar structure model with pre-computed grids of stellar atmospheres.

1.  **Interpolation**: At each timestep, the module takes the star's current surface parameters—Effective Temperature (:math:`T_{\rm eff}`), Surface Gravity (:math:`\log g`), and Metallicity ([M/H])—and queries a user-specified library of stellar atmospheres (defined in ``stellar_atm``). It interpolates within this grid to construct a specific Spectral Energy Distribution (SED) for the stars current features.
    
2.  **Convolution**: This specific SED is then convolved with filter transmission curves (defined in ``instrument``) to calculate the flux passing through each filter.

3.  **Integration**: The fluxes are converted into magnitudes using the user-selected magnitude system (AB, ST, or Vega).

3. Running the Test Suite
=========================

The module automatically appends new columns to the ``history.data`` file. You **do not** need to manually add these columns to your ``history_columns.list``. The module handles this dynamically based on the filters found in your specified instrument directory.

The standard output includes:

* ``Mag_bol``: The absolute bolometric magnitude.
* ``Flux_bol``: The bolometric flux (in cgs).
* ``[Filter_Name]``: A magnitude column for every filter file found in the directory specified by ``instrument``.

For example, if your instrument folder contains ``V.dat`` and ``B.dat``, your history file will automatically gain columns named ``V`` and ``B``.

4. Inlist Options & Parameters
==============================

The colors module is controlled via the ``&colors`` namelist. Below is a detailed guide to the key parameters.

instrument
----------
**Type:** `string`
      instrument = 

**Default:** `'/colors/data/filters/Generic/Johnson'`

This points to the directory containing the filter transmission curves you wish to use. The path must be structured as ``facility/instrument``.

* The directory must contain a file named after the instrument (e.g., ``Johnson``) which acts as an index.
* The module will read every ``.dat`` file listed in that directory and create a corresponding history column for it.

**Example:**

.. code-block:: fortran

   instrument = '/colors/data/filters/GAIA/GAIA'


stellar_atm
-----------
**Type:** `string`
**Default:** `'/colors/data/stellar_models/Kurucz2003all/'`

Specifies the path to the directory containing the grid of stellar atmosphere models. This directory must contain:

1.  **lookup_table.csv**: A map linking filenames to physical parameters (:math:`T_{\rm eff}`, :math:`\log g`, [M/H]).
2.  **SED files**: The actual spectra (text or binary format).
3.  **flux_cube.bin**: (Optional but recommended) A binary cube for rapid interpolation.

The module queries this grid using the star's current parameters. If the star evolves outside the grid boundaries, the module may clamp to the nearest edge or extrapolate, depending on internal settings.


distance
--------
**Type:** `float` (double)
**Default:** `3.0857d19` (10 parsecs in cm)

The distance to the star in centimeters. 

* This value is used to convert surface flux to observed flux.
* **Default Behavior:** It defaults to 10 parsecs (:math:`3.0857 \times 10^{19}` cm), resulting in **Absolute Magnitudes**.
* **Custom Usage:** You can set this to a specific source distance (e.g., distance to Betelgeuse) to calculate Apparent Magnitudes.


make_csv
--------
**Type:** `logical`
**Default:** `.false.`

If set to ``.true.``, the module exports the full calculated SED at every profile interval. 

* **Destination:** Files are saved to the directory defined by ``colors_results_directory``.
* **Format:** CSV files containing Wavelength vs. Flux.
* **Use Case:** useful for debugging or plotting the full spectrum of the star at a specific age.


colors_results_directory
------------------------
**Type:** `string`
**Default:** `'SED'`

The folder where csv files (if ``make_csv = .true.``) and other debug outputs are saved.


mag_system
----------
**Type:** `string`
**Default:** `'ST'`

Defines the zero-point system for magnitude calculations. Options are:

* ``'AB'``: Based on a flat spectral flux density of 3631 Jy.
* ``'ST'``: Based on a flat spectral flux density per unit wavelength.
* ``'Vega'``: Calibrated such that the star Vega has magnitude 0 in all bands.


vega_sed
--------
**Type:** `string`
**Default:** `'/colors/data/stellar_models/vega_flam.csv'`

Required only if ``mag_system = 'Vega'``. This points to the reference SED file for Vega. The default path points to a file provided with the MESA data distribution.


5. Data Preparation (SED_Tools)
===============================

The ``colors`` module requires input data (atmospheres and filters) to be formatted in a specific way. To assist with this, we provide the **SED_Tools** repository.

**Repository:** `SED_Tools <https://github.com/nialljmiller/SED_Tools>`_

This tool allows you to:

1.  Download spectra from repositories like **SVO**, **MSG**, and **MAST**.
2.  Download filter profiles from the **SVO Filter Profile Service**.
3.  **Rebuild** raw data into the binary formats (``.bin``, ``.h5``) and directory structures required by MESA.

**Workflow Summary:**
Download Data (via SED_Tools) :math:`\rightarrow` Rebuild for MESA :math:`\rightarrow` Point ``inlist`` to new directories.

6. Defaults Reference
=====================

Below are the default values for the colors module parameters as defined in ``colors.defaults``. These are used if you do not override them in your inlist.

.. code-block:: fortran

      use_colors = .false.
      instrument = '/colors/data/filters/Generic/Johnson'
      vega_sed = '/colors/data/stellar_models/vega_flam.csv'
      stellar_atm = '/colors/data/stellar_models/Kurucz2003all/'
      distance = 3.0857d19  ! 10 parsecs in cm (Absolute Magnitude)
      make_csv = .false.
      colors_results_directory = 'SED'
      mag_system = 'ST'

Visual Summary of Data Flow
===========================

.. code-block:: text

   +----------------+       +------------------+       +-------------------+
   |  Stellar Model |       |  Stellar Atm     |       |   Filter Curves   |
   | (Teff, logg, Z)| ----> |     Library      | ----> |    (Instrument)   |
   +----------------+       +------------------+       +-------------------+
           |                         |                           |
           v                         v                           v
   +-----------------------------------------------------------------------+
   |                        MESA COLORS MODULE                             |
   | 1. Query Atm Grid -> Interpolate SED                                  |
   | 2. Apply Distance -> Flux_bol                                         |
   | 3. Convolve SED + Filters -> Band Flux                                |
   | 4. Apply Zero Point (Vega/AB/ST) -> Magnitude                         |
   +-----------------------------------------------------------------------+
           |
           v
   +----------------------+
   |    history.data      |
   | age, Teff, ...       |
   | Mag_bol, Flux_bol    |
   | V, B, I, ...         |
   +----------------------+




================================================
================================================

7. Python Helper Scripts
========================

A collection of Python scripts is provided in the ``python_helpers/`` directory to assist with real-time monitoring, visualization, and analysis of the colors module output.

Dependencies:
  * ``matplotlib``
  * ``numpy``
  * ``mesa_reader`` (ensure this is installed/accessible)
  * ``ffmpeg`` (optional, for movie generation)

HISTORY_check.py
----------------
**Usage:** ``python python_helpers/HISTORY_check.py``

A real-time dashboard that monitors your ``history.data`` file as MESA runs. It automatically refreshes when new data is written.

**Plots:**
  1. **HR Diagram (Color vs. Magnitude):** Points colored by evolutionary phase.
  2. **Theoretical HR (Teff vs. Log L):** Standard theoretical track.
  3. **Color Evolution:** Color index vs. Age.
  4. **Light Curves:** Absolute magnitude vs. Age for all filters.

**Requirements:** Ensure ``phase_of_evolution`` is present in your ``history_columns.list`` for phase-based coloring.

SED_check.py
------------
**Usage:** ``python python_helpers/SED_check.py``

Monitors the ``colors_results_directory`` (default: ``SED/``) for new CSV output.

**Features:**
  * Plots the full high-resolution stellar spectrum (black line).
  * Plots the filter-convolved fluxes (colored lines corresponding to filters).
  * Displays stellar parameters (Mass, Z, Distance) parsed from your inlist.
  * Updates automatically if ``make_csv = .true.`` in your inlist and MESA overwrites the files.

interactive_cmd_3d.py
---------------------
**Usage:** ``python python_helpers/interactive_cmd_3d.py``

Generates an interactive 3D Color-Magnitude Diagram.

* **X-axis:** Color Index (e.g., B-V or Gbp-Grp).
* **Y-axis:** Magnitude (e.g., V or G).
* **Z-axis:** User-selectable column from the history file (e.g., ``Interp_rad``, ``star_age``, ``log_R``). The script will prompt you to choose a column at runtime.

Movie Makers
------------
Scripts to generate MP4 animations of your run. Requires ``ffmpeg``.

* **make_history_movie.py**: Creates ``history.mp4``, an animated version of the ``HISTORY_check.py`` dashboard showing the evolution over time.
* **make_CMD_InterpRad_movie.py**: Creates ``cmd_interp_rad_rotation.mp4``, a rotating 3D view of the Color-Magnitude Diagram with the Interpolation Radius on the Z-axis. Useful for visualizing grid coverage and interpolation quality.    