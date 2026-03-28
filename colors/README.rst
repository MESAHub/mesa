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

1.  **Interpolation**: At each timestep, the module takes the star's current surface parameters, Effective Temperature (:math:`T_{\rm eff}`), Surface Gravity (:math:`\log g`), and Metallicity ([Fe/H]), and queries a user-specified library of stellar atmospheres (defined in ``stellar_atm``). It interpolates within this grid to construct a specific Spectral Energy Distribution (SED) for the stars current features.

2.  **Convolution**: This specific SED is then convolved with filter transmission curves (defined in ``instrument``) to calculate the flux passing through each filter.

3.  **Integration**: The fluxes are converted into magnitudes using the user-selected magnitude system (AB, ST, or Vega).

Inlist Options & Parameters
===========================

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

1.  **lookup_table.csv**: A map linking filenames to physical parameters (:math:`T_{\rm eff}`, :math:`\log g`, [Fe/H]).
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


z_over_x_ref
------------

**Default:** ``0.0230057173972`` (the GS98 solar :math:`Z/X` ratio)

Reference metal-to-hydrogen ratio used to convert the model photospheric
composition into the metallicity coordinate passed to the atmosphere grid:

.. math::

   [Fe/H] = \log_{10}\left(\frac{(Z/X)_\mathrm{phot}}{(Z/X)_\mathrm{ref}}\right)

The default is chosen to match the default ``stellar_atm =
'/data/colors_data/stellar_models/Kurucz2003all/'`` grid, which is based on
ATLAS9 models using the Grevesse & Sauval (1998) solar abundance scale. The
numeric value ``0.0230057173972`` is computed from the MESA GS98 constants
``Z_solar = 0.0169`` and ``Y_solar = 0.2485``, so
``(Z/X)_ref = 0.0169 / (1 - 0.0169 - 0.2485)``.

If you switch to a different atmosphere library, you must set ``z_over_x_ref``
so that its reference composition matches the metallicity labels in that
library's ``lookup_table.csv``. For custom abundance patterns, making this
choice is the responsibility of the user.

If the photospheric metallicity, expressed in terms of ``Z/X``, falls outside
the metallicity range available in the tabulated atmosphere lookup table, then
MESA uses the nearest metallicity in the table for the atmosphere
interpolation. If the photospheric hydrogen mass fraction or metal mass
fraction is not positive, then MESA cannot form
``[Fe/H] = log10((Z/X)/(Z/X)_ref)``, and so it uses the lowest metallicity in
the table instead. This is a fallback safeguard and not a dedicated treatment
for H-free atmospheres.

For the shipped ``Kurucz2003all`` table in this checkout, the tabulated
metallicity range is ``[Fe/H] = -2.5`` to ``[Fe/H] = 4.0``. With the default
``z_over_x_ref = 2.30057173972d-2``, the lower edge corresponds to an adopted
``Z/X`` of about ``7.275d-5``.

**Example:**

.. code-block:: fortran

      z_over_x_ref = 2.30057173972d-2

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
============================

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
==================

Below are the default values for the colors module parameters as defined in ``colors.defaults``. These are used if you do not override them in your inlist.

.. code-block:: fortran

      use_colors = .false.
      instrument = '/data/colors_data/filters/Generic/Johnson'
      vega_sed = '/data/colors_data/stellar_models/vega_flam.csv'
      stellar_atm = '/data/colors_data/stellar_models/Kurucz2003all/'
      distance = 3.0857d19  ! 10 parsecs in cm (Absolute Magnitude)
      z_over_x_ref = 2.30057173972d-2  ! GS98 solar Z/X for Kurucz2003all
      make_csv = .false.
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
