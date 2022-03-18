.. _surface_effects:

***************
surface_effects
***************

Tests the implementation of the various surface effect corrections
available in MESA using the output model from
``simplex_solar_calibration``.  First, it adds a surface effect
correction to the model frequencies and checks that the input
parameters are recovered precisely.  Second, it fits the surface
effect correction to observed BiSON frequencies and checks that the
fit parameters satisfy some very loose constraints.
The observed frequencies are from
`Broomhall et al. (2009) <https://ui.adsabs.harvard.edu/abs/2009MNRAS.396L.100B>`__
and
`Davies et al. (2014) <https://ui.adsabs.harvard.edu/abs/2014MNRAS.439.2025D/abstract>`__
for frequencies above and below 2000 μHz, respectively,
and can be retrieved in plain text from the
`BiSON Open Data Portal <http://bison.ph.bham.ac.uk/portal/frequencies>`__.

If the first test fails, it is almost certainly because of a bug
in the implementation of the failing surface term.

If the second test fails, it might be because the input model is
somehow no longer close enough to the Sun.  Try rerunning the solar
calibration and copying the best-fit model to ``s1.mod``.  You can
also update the target parameter values using the Python script
``get_targets.py``, whose output can be written to
``src/targets.inc``. i.e., run ::

  python3 get_targets.py > src/targets.inc

Last-Update: 2022-03-11 (mesa 4ee4cc0) by Warrick Ball
