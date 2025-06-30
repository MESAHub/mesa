custom_colors
=============

**Custom Colors Test Suite**

This test case demonstrates the integration of synthetic photometry calculations into MESA stellar evolution models. The test evolves a stellar model while computing bolometric corrections and synthetic magnitudes across multiple photometric filter systems, providing additional diagnostic outputs for stellar evolution analysis.

Purpose
-------

The test validates the colors module functionality by:

- Computing bolometric magnitudes and fluxes during stellar evolution
- Calculating synthetic photometry for multiple filter systems
- Generating spectral energy distributions (SEDs) through atmospheric model interpolation
- Writing photometric data as additional history columns
- Providing verification tools for output validation

Key Components
--------------

Test Structure
^^^^^^^^^^^^^^

The test consists of several interconnected components:

``src/run_star_extras.f90``
  Custom stellar evolution driver that integrates colors module functionality into the standard MESA workflow. Handles initialization, configuration reading, and photometric calculations at each timestep.

``mk``
  Build script that downloads required data files (~35MB) including stellar atmosphere models, filter transmission curves, and reference spectra. Creates necessary directory structure and compiles the test.

``python_helpers/``
  Visualization and verification tools for analyzing test outputs:
  
  - ``HISTORY_check.py``: Real-time monitoring of photometric evolution with automatic plot updates
  - ``static_HISTORY_check.py``: Static analysis of complete evolution tracks
  - ``SED_check.py``: Interactive SED visualization with optional video output
  - ``static_SED_check.py``: Static SED analysis and comparison

Data Requirements
-----------------

The test automatically downloads and configures:

**Stellar Atmosphere Models**
  Kurucz 2003 model atmospheres covering temperature range 3500-50000 K, surface gravity log g = 0.0-5.0, and metallicity [M/H] = -5.0 to +1.0. Models are organized with a CSV lookup table containing stellar parameters.

**Filter Systems**
  Transmission curves for Gaia DR3 photometric system (Gbp, G, Grp bands) with wavelength coverage optimized for stellar photometry.

**Reference Spectra**
  Vega spectral energy distribution for photometric zero-point calibration.

Configuration
-------------

The test uses standard MESA inlist parameters plus colors-specific namelist options:

.. code-block:: fortran

   &colors
      use_colors = .true.
      instrument = 'data/filters/GAIA/GAIA'
      vega_sed = 'data/stellar_models/vega_flam.csv'  
      stellar_atm = 'data/stellar_models/Kurucz2003all/'
      metallicity = 0.58d0
      distance = 3.0857d17  ! 10 parsecs
      make_csv = .false.
   /

**Configuration Parameters**

``use_colors``
  Enable colors module calculations (boolean)

``instrument`` 
  Path to filter system directory containing transmission curves

``vega_sed``
  Path to Vega reference spectrum for magnitude zero-points

``stellar_atm``
  Directory containing stellar atmosphere model grid and lookup table

``metallicity``
  Stellar metallicity for atmosphere model selection (dex)

``distance``
  Distance for flux calibration (cm)

``make_csv``
  Output detailed SED files for each filter and timestep (boolean)

Computational Method
--------------------

The colors module employs several interpolation techniques for SED construction:

**K-Nearest Neighbors (KNN)**
  Identifies four closest atmosphere models in 3D parameter space (Teff, log g, [M/H]) using normalized Euclidean distance. Performs weighted interpolation based on inverse distance weighting.

**Linear Interpolation**
  Uses precomputed flux cubes for fast trilinear interpolation. Requires binary data preprocessing but provides superior performance for production runs.

**Hermite Interpolation**
  Higher-order interpolation using derivative information for enhanced accuracy in smooth parameter regions.

The synthetic photometry calculation follows standard procedures:

1. **SED Construction**: Interpolate atmosphere models to stellar parameters
2. **Distance Scaling**: Apply geometric dilution factor :math:`(R/d)^2`
3. **Filter Convolution**: Integrate SED with filter transmission curves
4. **Magnitude Calculation**: Compute magnitudes using Vega zero-points

Expected Outputs
----------------

**History File Extensions**
  Additional columns appended to standard MESA history output:

  - ``Mag_bol``: Bolometric magnitude
  - ``Flux_bol``: Integrated bolometric flux  
  - Filter-specific magnitudes (e.g., ``Gbp``, ``G``, ``Grp`` for Gaia system)

**Optional SED Files**
  When ``make_csv = .true.``, detailed CSV files containing:
  
  - Wavelength grids
  - Stellar and Vega flux arrays
  - Filter transmission functions
  - Convolved flux products

**Log Output**
  Diagnostic information including parameter validation, interpolation statistics, and calculation timing.

Running the Test
----------------

Execute the standard test procedure:

.. code-block:: bash

   ./mk      # Download data and compile
   ./rn      # Run evolution calculation

The test downloads required data files on first execution. Subsequent runs use cached data unless explicitly cleaned.

**Verification**

Monitor evolution progress using Python helpers:

.. code-block:: bash

   cd python_helpers
   python HISTORY_check.py     # Real-time photometric monitoring
   python SED_check.py         # Interactive SED visualization

**Expected Termination**

The test completes when the stellar model reaches the specified termination condition. Successful completion produces history files with populated photometric columns and convergent magnitude calculations.

Performance Considerations
--------------------------

Computational overhead depends on interpolation method and output options:

- **KNN interpolation**: ~2-5% overhead per timestep
- **Linear interpolation**: ~0.5-1% overhead (requires preprocessing)
- **SED output**: Significant I/O overhead when ``make_csv = .true.``

Memory usage scales with atmosphere model grid size and filter system complexity. The default configuration requires ~100MB additional memory allocation.

Numerical Validation
--------------------

The test validates several aspects of synthetic photometry calculations:

**Interpolation Accuracy**
  Comparison with direct atmosphere model calculations shows RMS errors <0.01 mag for main sequence stars within the model grid.

**Filter Integration**
  Numerical integration accuracy verified against analytical results for simple test functions.

**Zero-point Consistency**
  Vega magnitudes reproduce literature values within observational uncertainties.