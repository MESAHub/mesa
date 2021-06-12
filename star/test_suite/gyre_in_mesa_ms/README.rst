.. _gyre_in_mesa_ms:

***************
gyre_in_mesa_ms
***************

This test case checks the implementation of GYRE in MESA for a 1 Msun, Z=0.02 metallicity, model evolving from the zero-age main sequence to core hydrogen depletion.

This test case has 2 parts.

* Part 1 (``inlist_zams``) builds a 1 Msun, Z=0.02 metallicity, pre-main sequence model and evolves it to the main sequence.

* Part 2 (``inlist_gyre_in_mesa_ms``) continues the evolution until the central hydrogen mass fraction drops below 1e-3. During the evolution the ``run_star_extras.f90`` calls GYRE, processes the GYRE output, and searches for a p-mode frequency of 2e-4 Hz (200 microHz). Close matches to this target frequency are reported in the terminal:

.. code-block:: console

   model   order           freq (Hz)             P (sec)             P (day)        growth (day)              growth    cycles to double
      96       1          0.2049E-03           4880.2668              0.0565       56532479.2062              0.0000     1000848182.8715
 matched target frequency
      96       2          0.3296E-03           3033.8614              0.0351       25751771.8650              0.0000      733373357.5769

 ...

   model   order           freq (Hz)             P (sec)             P (day)        growth (day)              growth    cycles to double
     110       1          0.1963E-03           5093.3666              0.0590       52889722.8415              0.0000      897181064.5062
 matched target frequency
     110       2          0.3172E-03           3152.2376              0.0365       23156734.0926              0.0000      634705276.0032


Last-Updated: 11Jun2021 (MESA 5be9e57) by fxt.


