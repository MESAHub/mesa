.. _gyre_in_mesa_envelope:

*********************
gyre_in_mesa_envelope
*********************

This test case checks the implementation of GYRE in MESA for the envelope of beta Cephei model.

This test case has 4 parts.

* Part 1 (``inlist_zams``) builds a 12.0 Msun, Z=0.02 metallicity, pre-main sequence model and evolves it to the zero age main sequence.

* Part 2 (``inlist_tams``) continues the evolution until the central hydrogen mass fraction is less than 1e-4.

* Part 3 (``inlist_remove_center``) removes cells with a temperature greater than 2e6 K and setting new inner boundary conditions for mass, radius, velocity, and luminosity. This step excises about 11.96 Msun from the model and then run terminates (i.e., no evolution is done).

* Part 4 (``inlist_pulse``) continues the evolution for 5 timesteps. During this evolution the ``run_star_extras.f90`` calls GYRE, processes the GYRE output, and searches for a frequency of 6.2e-5 Hz (62 microHz). Close matches to this target frequency are reported in the terminal:

.. code-block:: console

   model   order           freq (Hz)             P (sec)             P (day)        growth (day)              growth    cycles to double
     317       1          0.6149E-04          16263.9143              0.1882            104.0585              0.0018            552.7976
 matched target frequency
     317       2          0.7806E-04              0.0000              0.0000              stable

 ...

   model   order           freq (Hz)             P (sec)             P (day)        growth (day)              growth    cycles to double
     321       1          0.6137E-04          16294.6111              0.1886            106.1335              0.0018            562.7586
 matched target frequency
     321       2          0.7790E-04              0.0000              0.0000              stable


Last-Updated: 11Jun2021 (MESA 5be9e57) by fxt.
