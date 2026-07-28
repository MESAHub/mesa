Overview of colors module
=========================


.. toctree::
   :maxdepth: 2
   :hidden:

   defaults

The ``colors`` module calculates synthetic photometry during stellar evolution.
The module computes bolometric and synthetic magnitudes by interpolating stellar atmosphere model grids and convolving with photometric filter transmission curves.

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


.. _colors_interpolation_limitations:

Interpolation limitations
-------------------------

The accuracy of the synthetic photometry depends on the atmosphere grid
resolving the variation of the spectral energy distribution (SED) with
:math:`T_{\rm eff}`, :math:`\log g`, and metallicity.  In particular, if the
SED changes substantially over a parameter interval narrower than the grid
spacing, values between the tabulated atmosphere models are underdetermined.
No choice of interpolation method can reconstruct a transition that is absent
from the grid.

Hermite interpolation can overshoot near a steep transition.  The bounded
Hermite implementation detects non-positive or sufficiently inconsistent
results and falls back to multilinear interpolation.  This prevents some
numerical failures, but does not establish that the fallback SED is physically
accurate.  A multilinear result can be positive and remain within the
pointwise range of the cell-corner SEDs while still mixing spectra from
different physical regimes.  The resulting continuum, spectral features,
bolometric flux, and filter magnitudes may therefore be inaccurate without
producing an interpolation error.

``Interp_rad`` is the normalized distance to the nearest atmosphere-grid
point.  It does not measure the local SED gradient, compare interpolation
methods, or detect an unresolved transition.  A small value of ``Interp_rad``
is consequently not a guarantee of interpolation accuracy.

Users should validate results in regions with rapid spectral changes against
the original atmosphere models and diagnostic SED output.  The atmosphere grid
must sample such transitions finely enough for the intended application.  If
it does not, switching between Hermite, multilinear, or nearest-neighbor
methods is not a reliable remedy; a more finely sampled atmosphere grid or an
application-specific physical treatment is required.

Filter transmission is treated as zero outside the wavelength interval
tabulated in each filter file.  Filter tables should cover the complete
intended passband.
