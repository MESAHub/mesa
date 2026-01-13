.. _custom_colors:

*************
custom_colors
*************

This test suite was tested against SDK 25.12.1

This test suite case demonstrates the functionality of the MESA ``colors`` module, a framework introduced in MESA r25.10.1 for calculating synthetic photometry and bolometric quantities during stellar evolution.

What is MESA colors?
=======================

MESA colors is a post-processing and runtime module that allows users to generate "observer-ready" data directly from stellar evolution models. Instead of limiting output to theoretical quantities like Luminosity (:math:`L`) and Surface Temperature (:math:`T_{\rm eff}`), the colors module computes:

* **Bolometric Magnitude** (:math:`M_{\rm bol}`)
* **Bolometric Flux** (:math:`F_{\rm bol}`)
* **Synthetic Magnitudes** in specific photometric filters (e.g., Johnson V, Gaia G, 2MASS J).

This bridges the gap between theoretical evolutionary tracks and observational color-magnitude diagrams (CMDs).

How does the MESA colors module work?
========================================

The module operates by coupling the stellar structure model with pre-computed grids of stellar atmospheres.

1.  **Interpolation**: At each timestep, the module takes the star's current surface parameters—Effective Temperature (:math:`T_{\rm eff}`), Surface Gravity (:math:`\log g`), and Metallicity ([M/H])—and queries a user-specified library of stellar atmospheres (defined in ``stellar_atm``). It interpolates within this grid to construct a specific Spectral Energy Distribution (SED) for the stars current features.

2.  **Convolution**: This specific SED is then convolved with filter transmission curves (defined in ``instrument``) to calculate the flux passing through each filter.

3.  **Integration**: The fluxes are converted into magnitudes using the user-selected magnitude system (AB, ST, or Vega).

The Test Suite
=========================

This test suite evolves a complete stellar model (from the pre–main sequence onward) while the ``colors`` module runs *continuously* in the background.
At every timestep, MESA computes synthetic photometry by interpolating a stellar atmosphere grid and convolving the resulting SED with the filters you specify in the inlist.

During the run, the module automatically appends new photometric columns to the ``history.data`` file.
You **do not** need to list these in ``history_columns.list``—the module detects the available filters by inspecting the directory defined in the ``instrument`` parameter.

What the Test Suite Produces
----------------------------

The standard output of the test suite includes:

* ``Mag_bol``
  The bolometric magnitude computed from the star's instantaneous bolometric flux.

* ``Flux_bol``
  The bolometric flux (cgs units) after distance dilution is applied.
  The default distance is 10 pc, producing **absolute magnitudes** unless changed.

* ``[Filter_Name]``
  A synthetic magnitude column for **every** filter in the instrument directory.

For example, if your filter directory contains:

.. code-block:: text

   filters/Generic/Johnson/
       B.dat
       V.dat
       R.dat
       Johnson

then your history file will include:

``B``, ``V``, ``R``

as new magnitude columns generated automatically at runtime.

What the Test Suite Actually Does
---------------------------------

The provided ``inlist_colors`` configures a full evolution run with the colors module activated:

* Starts from a **pre–main-sequence model**
* Evolves the model through multiple phases while computing synthetic photometry
* Uses the Johnson filter set
* Uses the Kurucz2003all atmosphere grid
* Outputs bolometric and filter-specific magnitudes to ``history.data`` every step

Because the test suite's inlist also defines a set of PGSTAR panels, you automatically
get real-time plots of:

* HR diagram (log L vs. log Teff)
* A light curve based on any synthetic magnitude (the test suite uses ``V``)

Real-Time Visualization (Enabled by Default)
--------------------------------------------

The test suite's ``pgstar`` block is configured so that, as the star evolves:

* Panel 1: HR diagram (theoretical)
* Panel 2: Light curve in the Johnson V band

These update automatically as the model runs.

Purpose of the Test Suite
-------------------------

This test problem is designed to demonstrate:

1. That MESA can compute synthetic photometry **at runtime**, without external tools.
2. How atmosphere grids and filters affect magnitude evolution.
3. How to configure your own inlists for scientific use.
4. How the ``colors`` module integrates with PGSTAR visualization.
5. How synthetic magnitudes appear in ``history.data`` and how to use them for CMDs, light curves, and population modeling.



Inlist Options & Parameters
==============================

The colors module is controlled via the ``&colors`` namelist. Below is a detailed guide to the key parameters.

instrument
----------
**Default:** `'/data/colors_data/filters/Generic/Johnson'`

This points to the directory containing the filter transmission curves you wish to use. The path must be structured as ``facility/instrument``.

* The directory must contain a file named after the instrument (e.g., ``Johnson``) which acts as an index.
* The module will read every ``.dat`` file listed in that directory and create a corresponding history column for it.

**Example:**

.. code-block:: fortran

   instrument = '/data/colors_data/filters/GAIA/GAIA'


stellar_atm
-----------
**Default:** `'/data/colors_data/stellar_models/Kurucz2003all/'`

Specifies the path to the directory containing the grid of stellar atmosphere models. This directory must contain:

1.  **lookup_table.csv**: A map linking filenames to physical parameters (:math:`T_{\rm eff}`, :math:`\log g`, [M/H]).
2.  **SED files**: The actual spectra (text or binary format).
3.  **flux_cube.bin**: (Optional but recommended) A binary cube for rapid interpolation.

The module queries this grid using the star's current parameters. If the star evolves outside the grid boundaries, the module may clamp to the nearest edge or extrapolate, depending on internal settings.

**Example:**

.. code-block:: fortran

   stellar_atm = '/data/colors_data/stellar_models/sg-SPHINX/'


distance
--------
**Default:** `3.0857d19` (10 parsecs in cm)

The distance to the star in centimeters.

* This value is used to convert surface flux to observed flux.
* **Default Behavior:** It defaults to 10 parsecs (:math:`3.0857 \times 10^{19}` cm), resulting in **Absolute Magnitudes**.
* **Custom Usage:** You can set this to a specific source distance (e.g., distance to Betelgeuse) to calculate Apparent Magnitudes.

**Example:**

.. code-block:: fortran

    distance = 5.1839d20

make_csv
--------
**Default:** `.false.`

If set to ``.true.``, the module exports the full calculated SED at every profile interval.

* **Destination:** Files are saved to the directory defined by ``colors_results_directory``.
* **Format:** CSV files containing Wavelength vs. Flux.
* **Use Case:** useful for debugging or plotting the full spectrum of the star at a specific age.

**Example:**

.. code-block:: fortran

      make_csv = .true.


sed_per_model
--------

**Default:** `.false.`

If set to ``.true.`` AND ``make_csv`` set to ``.true.``, the module exports the full calculated SED at every profile interval WITH the model number suffix.

!!!WARNING: Enabling this feature will cause the colors_results_directory to drastically increase in size. DO NOT enable this without first ensureing you have appropriate storage. !!!

* **Destination:** Files are saved to the directory defined by ``colors_results_directory``.
* **Format:** CSV files containing Wavelength vs. Flux with model number suffic on file name. 
* **Use Case:** useful for seeing SED evolution.

**Example:**

.. code-block:: fortran

      sed_per_model = .true.      


colors_results_directory
------------------------
**Default:** `'SED'`

The folder where csv files (if ``make_csv = .true.``) and other debug outputs are saved.

**Example:**

.. code-block:: fortran

      colors_results_directory = 'sed'


mag_system
----------
**Default:** `'Vega'`

Defines the zero-point system for magnitude calculations. Options are:

* ``'AB'``: Based on a flat spectral flux density of 3631 Jy.
* ``'ST'``: Based on a flat spectral flux density per unit wavelength.
* ``'Vega'``: Calibrated such that the star Vega has magnitude 0 in all bands.

**Example:**

.. code-block:: fortran

      mag_system = 'AB'


vega_sed
--------
**Default:** `'/data/colors_data/stellar_models/vega_flam.csv'`

Required only if ``mag_system = 'Vega'``. This points to the reference SED file for Vega. The default path points to a file provided with the MESA data distribution.

**Example:**

.. code-block:: fortran

      vega_sed = '/another/file/for/vega_SED.csv'


Data Preparation (SED_Tools)
===============================

The ``colors`` module requires pre-processed stellar atmospheres and filter
profiles organized in a very specific directory structure. To automate this
entire workflow, we provide the dedicated repository:

**Repository:** `SED_Tools <https://github.com/nialljmiller/SED_Tools>`_

SED_Tools downloads, validates, and converts raw spectral atmosphere grids and
filter transmission curves from the following public archives:

* `SVO Filter Profile Service <http://svo2.cab.inta-csic.es/theory/fps/>`_
* `MAST BOSZ Stellar Atmosphere Library <https://archive.stsci.edu/prepds/bosz/>`_
* `MSG / Townsend Atmosphere Grids <https://www.astro.wisc.edu/~townsend/msg/>`_

These sources provide heterogeneous formats and file organizations. SED_Tools
standardizes them into the exact structure required by MESA:

* ``lookup_table.csv``
* Raw SED files (text or/and HDF5)
* ``flux_cube.bin`` (binary cube for fast interpolation)
* Filter index files and ``*.dat`` transmission curves


SED_Tools produces:

.. code-block:: text

    data/
    ├── stellar_models/<ModelName>/
    │   ├── flux_cube.bin
    │   ├── lookup_table.csv
    │   ├── *.txt / *.h5
    └── filters/<Facility>/<Instrument>/
        ├── *.dat
        └── <Instrument>

These directories can be copied or symlinked directly into your MESA
installation.

A browsable online mirror of the processed SED_Tools output is also available:

`SED Tools Web Interface (mirror) <https://nillmill.ddns.net/sed_tools/>`_

This server provides a live view of:

* All downloaded stellar atmosphere grids
* All available filter facilities and instruments
* File counts, disk usage, and metadata
* Direct links to the directory structure used by MESA


Defaults Reference
=====================

Below are the default values for the colors module parameters as defined in ``colors.defaults``. These are used if you do not override them in your inlist.

.. code-block:: fortran

      use_colors = .false.
      instrument = '/data/colors_data/filters/Generic/Johnson'
      vega_sed = '/data/colors_data/stellar_models/vega_flam.csv'
      stellar_atm = '/data/colors_data/stellar_models/Kurucz2003all/'
      distance = 3.0857d19  ! 10 parsecs in cm (Absolute Magnitude)
      make_csv = .false.
      colors_results_directory = 'SED'
      mag_system = 'Vega'
      vega_sed = '/data/colors_data/stellar_models/vega_flam.csv'
Visual Summary of Data Flow
===========================

.. code-block:: text

   +----------------+
   |   MESA Model   |
   | (Teff, logg, Z)|
   +----------------+
           |
           v
   +-------------------------------------------------------------------------+
   |                        MESA COLORS MODULE                               |
   | 1. Query Stellar Atmosphere Grid with input model                       |
   | 2. Interpolate grid to construct specific SED                           |
   | 3. Convolve SED with filters to generate band flux                      |
   | 2. Apply distance flux dilution to generate bolometric flux -> Flux_bol |
   | 4. Apply zero point (Vega/AB/ST) to  generate magnitudes                |
   |                                    (Both bolometric and per filter)     |
   +-------------------------------------------------------------------------+
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

Python Helper Scripts
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
