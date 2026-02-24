.. _custom_colors:

******
Colors
******

This test suite case demonstrates the functionality of the MESA ``colors`` module, a framework introduced in MESA r25.10.1 for calculating synthetic photometry and bolometric quantities during stellar evolution.

What is MESA colors?
====================

MESA colors is a post-processing and runtime module that allows users to generate "observer-ready" data directly from stellar evolution models. Instead of limiting output to theoretical quantities like Luminosity (:math:`L`) and Surface Temperature (:math:`T_{\rm eff}`), the colors module computes:

* **Bolometric Magnitude** (:math:`M_{\rm bol}`)
* **Bolometric Flux** (:math:`F_{\rm bol}`)
* **Synthetic Magnitudes** in specific photometric filters (e.g., Johnson V, Gaia G, 2MASS J).

This bridges the gap between theoretical evolutionary tracks and observational color-magnitude diagrams (CMDs).

How does the MESA colors module work?
=====================================

The module operates by coupling the stellar structure model with pre-computed grids of stellar atmospheres.

1.  **Interpolation**: At each timestep, the module takes the star's current surface parameters—Effective Temperature (:math:`T_{\rm eff}`), Surface Gravity (:math:`\log g`), and Metallicity ([M/H])—and queries a user-specified library of stellar atmospheres (defined in ``stellar_atm``). It interpolates within this grid to construct a specific Spectral Energy Distribution (SED) for the star's current parameters.

2.  **Convolution**: This specific SED is then convolved with filter transmission curves (defined in ``instrument``) to calculate the flux passing through each filter.

3.  **Integration**: The fluxes are converted into magnitudes using the user-selected magnitude system (AB, ST, or Vega).

Inlist Options & Parameters
===========================

The colors module is controlled via the ``&colors`` namelist. Below is a detailed guide to the key parameters.

use_colors
----------

**Default:** ``.false.``

Master switch for the module. Must be set to ``.true.`` to enable any photometric output.

**Example:**

.. code-block:: fortran

   use_colors = .true.


instrument
----------

**Default:** ``'data/colors_data/filters/Generic/Johnson'``

Path to the directory containing the filter transmission curves to use. The path must be structured as ``facility/instrument``.

* Paths may be relative to ``$MESA_DIR``, relative to the working directory, or absolute.
* The directory must contain a file named after the instrument (e.g., ``Johnson``) which acts as an index listing the filters to load.
* The module will read every ``.dat`` file listed in that index and create a corresponding history column for it.

**Example:**

.. code-block:: fortran

   instrument = 'data/colors_data/filters/GAIA/GAIA'


stellar_atm
-----------

**Default:** ``'data/colors_data/stellar_models/Kurucz2003all/'``

Path to the directory containing the grid of stellar atmosphere models. Paths may be relative to ``$MESA_DIR``, relative to the working directory, or absolute. This directory must contain:

1.  **lookup_table.csv**: A map linking filenames to physical parameters (:math:`T_{\rm eff}`, :math:`\log g`, [M/H]).
2.  **SED files**: The actual spectra (text or binary format).
3.  **flux_cube.bin**: (Optional but recommended) A binary cube for rapid interpolation.

The module queries this grid using the star's current parameters. If the star evolves outside the grid boundaries, the module will clamp to the nearest edge.

**Example:**

.. code-block:: fortran

   stellar_atm = 'data/colors_data/stellar_models/sg-SPHINX/'


distance
--------

**Default:** ``3.0857d19`` (10 parsecs in cm)

The distance to the star in centimetres, used to convert surface flux to observed flux.

* **Default Behaviour:** At 10 parsecs (:math:`3.0857 \times 10^{19}` cm) the output is **Absolute Magnitudes**.
* **Custom Usage:** Set this to a specific source distance to calculate Apparent Magnitudes.

**Example:**

.. code-block:: fortran

   distance = 5.1839d20


make_csv
--------

**Default:** ``.false.``

If set to ``.true.``, the module exports the full calculated SED at every profile interval.

* **Destination:** Files are saved to the directory defined by ``colors_results_directory``.
* **Format:** CSV files containing Wavelength vs. Flux.
* **Use Case:** Useful for debugging or plotting the full spectrum of the star at a specific evolutionary age.

**Example:**

.. code-block:: fortran

   make_csv = .true.


sed_per_model
-------------

**Default:** ``.false.``

Requires ``make_csv = .true.``. If set to ``.true.``, each exported SED file is stamped with the model number, preserving one SED file per model rather than overwriting a single file.

.. warning::

   Enabling this feature will cause the ``colors_results_directory`` to grow very rapidly. Do not enable it without first ensuring you have sufficient storage.

* **Destination:** Files are saved to the directory defined by ``colors_results_directory``.
* **Format:** CSV files containing Wavelength vs. Flux, with the model number as a filename suffix.
* **Use Case:** Useful for tracking the full SED evolution of the star over time.

**Example:**

.. code-block:: fortran

   sed_per_model = .true.


colors_results_directory
------------------------

**Default:** ``'SED'``

The folder where CSV files (if ``make_csv = .true.``) and other outputs are saved.

**Example:**

.. code-block:: fortran

   colors_results_directory = 'sed'


mag_system
----------

**Default:** ``'Vega'``

Defines the zero-point system for magnitude calculations. Options are:

* ``'AB'``: Based on a flat spectral flux density of 3631 Jy.
* ``'ST'``: Based on a flat spectral flux density per unit wavelength.
* ``'Vega'``: Calibrated such that Vega has magnitude 0 in all bands.

**Example:**

.. code-block:: fortran

   mag_system = 'AB'


vega_sed
--------

**Default:** ``'data/colors_data/stellar_models/vega_flam.csv'``

Required only if ``mag_system = 'Vega'``. Points to the reference SED file for Vega, used to compute photometric zero-points. Paths may be relative to ``$MESA_DIR``, relative to the working directory, or absolute.

**Example:**

.. code-block:: fortran

   vega_sed = '/path/to/my/vega_SED.csv'


Data Preparation (SED_Tools)
============================

The ``colors`` module requires pre-processed stellar atmospheres and filter
profiles organised in a specific directory structure. To automate this
entire workflow, we provide the dedicated repository:

**Repository:** `SED_Tools <https://github.com/nialljmiller/SED_Tools>`_

SED_Tools downloads, validates, and converts raw spectral atmosphere grids and
filter transmission curves from the following public archives:

* `SVO Filter Profile Service <http://svo2.cab.inta-csic.es/theory/fps/>`_
* `MAST BOSZ Stellar Atmosphere Library <https://archive.stsci.edu/prepds/bosz/>`_
* `MSG / Townsend Atmosphere Grids <https://www.astro.wisc.edu/~townsend/msg/>`_

These sources provide heterogeneous formats and file organisations. SED_Tools
standardises them into the exact structure required by MESA:

* ``lookup_table.csv``
* Raw SED files (text and/or HDF5)
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
==================

Below are the default values for all user-facing ``colors`` module parameters as defined in ``colors.defaults``.

.. code-block:: fortran

      use_colors = .false.
      instrument = 'data/colors_data/filters/Generic/Johnson'
      stellar_atm = 'data/colors_data/stellar_models/Kurucz2003all/'
      vega_sed = 'data/colors_data/stellar_models/vega_flam.csv'
      distance = 3.0857d19  ! 10 parsecs in cm (Absolute Magnitude)
      make_csv = .false.
      sed_per_model = .false.
      colors_results_directory = 'SED'
      mag_system = 'Vega'

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
   | 4. Apply distance flux dilution to generate bolometric flux -> Flux_bol |
   | 5. Apply zero point (Vega/AB/ST) to generate magnitudes                 |
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