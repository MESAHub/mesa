Overview of colors module
=========================

The ``colors`` module calculates synthetic photometry during stellar evolution.
The module computes bolometric and synthetic magnitudes by interpolating stellar atmosphere model grids and convolving with photometric filter transmission curves.

The colors module is controlled via the ``&colors`` namelist with key options:

- ``use_colors``: Enable colors calculations (default ``.false.``)
- ``instrument``: Path to filter system directory
- ``stellar_atm``: Path to stellar atmosphere model grid
- ``vega_sed``: Vega spectrum for photometric zero points
- ``metallicity``: Metallicity of the star
- ``distance``: Distance to the star in cm
- ``make_csv``: Output detailed spectral energy distributions
- ``colors_results_directory``: Directory for output files

Filter-specific magnitude columns are automatically added to history output based on the selected instrument.

See the ``star/test_suite/custom_colors`` test suite case for usage examples.
