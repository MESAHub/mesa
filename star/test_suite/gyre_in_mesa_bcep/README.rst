.. _gyre_in_mesa_bcep:

*****************
gyre_in_mesa_bcep
*****************

This test case checks the implementation of GYRE in MESA for a beta Cephei model.

This test case has 2 parts.

* Part 1 (``inlist_zams``) builds a 12.0 Msun, Z=0.02 metallicity, pre-main sequence model and evolves it to the main sequence.

* Part 2 (``inlist_gyre_in_mesa_bcep``) continues the evolution until the central hydrogen mass fraction drops below 1e-3. During the evolution the ``run_star_extras.f90`` calls GYRE, processes the GYRE output, and searches for a frequency of 6e-5 Hz (60 microHz). Close matches to this target frequency are reported in the terminal:

.. code-block:: console

   model   order           freq (Hz)             P (sec)             P (day)        growth (day)              growth    cycles to double
      58       1          0.5846E-04          17105.7306              0.1980           7337.0594              0.0000          37059.0387
 matched target frequency
      58       2          0.7597E-04          13163.0272              0.1523            468.0420              0.0003           3072.1526

 ...

  model   order           freq (Hz)             P (sec)             P (day)        growth (day)              growth    cycles to double
      70       1          0.4636E-04          21572.1908              0.2497           2620.7283              0.0001          10496.4268
      70       2          0.5998E-04          16672.5212              0.1930            257.6324              0.0007           1335.0973
 matched target frequency

 ...

   model   order           freq (Hz)             P (sec)             P (day)        growth (day)              growth    cycles to double
      94       1          0.4663E-04          21443.2542              0.2482           1682.5284              0.0001           6779.3093
      94       2          0.6026E-04          16593.4740              0.1921            129.0619              0.0015            672.0079
 matched target frequency


Last-Updated: 10Jun2021 (MESA 5be9e57) by fxt.


