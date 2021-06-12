.. _gyre_in_mesa_spb:

****************
gyre_in_mesa_spb
****************

This test case checks the implementation of GYRE in MESA for a 5 Msun, Z=0.02 metallicity, model evolving from the zero-age main sequence to core hydrogen depletion; a slowly pulsating B-type star (SPB) stellar model.

This test case has 2 parts.

* Part 1 (``inlist_zams``) builds a 5 Msun, Z=0.02 metallicity, pre-main sequence model and evolves it to the main sequence.

* Part 2 (``inlist_gyre_in_mesa_ms``) continues the evolution until the central hydrogen mass fraction drops below 1e-2. During the evolution the ``run_star_extras.f90`` calls GYRE, processes the GYRE output, and searches for a non-radial g-mode frequency of 7e-6 Hz (7 microHz). Close matches to this target frequency are reported in the terminal:

.. code-block:: console

   model   order           freq (Hz)             P (sec)             P (day)        growth (day)              growth    cycles to double
      36     -19          0.6493E-05              0.0000              0.0000              stable
      36     -18          0.6823E-05         146555.6154              1.6962          50761.7451              0.0000          29925.9415
 matched target frequency
      36     -17          0.7220E-05         138502.8461              1.6030          25462.3183              0.0001          15883.7480

 ...

   model   order           freq (Hz)             P (sec)             P (day)        growth (day)              growth    cycles to double
      54     -19          0.7048E-05         141890.2396              1.6422         141050.7804              0.0000          85888.8354
 matched target frequency
      54     -18          0.7348E-05         136093.3780              1.5752          78902.0494              0.0000          50091.6148
      54     -17          0.7880E-05         126907.1215              1.4688         171001.2272              0.0000         116419.8340


Last-Updated: 11Jun2021 (MESA 5be9e57) by fxt.


