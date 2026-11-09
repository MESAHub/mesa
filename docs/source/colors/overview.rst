Overview of colors module
=========================


.. toctree::
   :maxdepth: 2
   :hidden:

   defaults

The ``colors`` module calculates synthetic photometry during stellar evolution.
The module computes bolometric and synthetic magnitudes by interpolating stellar atmosphere model grids and convolving with photometric filter transmission curves.

Atmosphere spectra are interpolated in effective temperature, surface gravity,
and metallicity using Hermite tensor interpolation, with safeguards against
unphysical fluxes. Negligible negative undershoots are set to zero; if the
resulting spectrum is unusable, the module falls back to trilinear interpolation.
These safeguards do not guarantee interpolation accuracy in sparsely sampled
regions of the atmosphere grid. The ``Interp_rad`` history column measures the
distance to the nearest atmosphere grid point and can help identify regions
that require closer inspection.

The colors module is controlled via the ``&colors`` namelist with key options:

- ``use_colors``: Enable colors calculations (default ``.false.``)
- ``instrument``: Path to filter system directory
- ``stellar_atm``: Path to stellar atmosphere model grid
- ``vega_sed``: Vega spectrum for photometric zero points
- ``distance``: Distance to the star in cm
- ``make_csv``: Output detailed spectral energy distributions
- ``colors_results_directory``: Directory for output files

Filter-specific magnitude columns are automatically added to history output based on the selected instrument.

See the ``star/test_suite/custom_colors`` test suite case for usage examples.
