.. _gyre_in_mesa_wd:

***************
gyre_in_mesa_wd
***************

This test case checks the implementation of GYRE in MESA for a cooling 0.85 Msun white dwarf model.

This test case has 1 part.

* Part 1 (``inlist_pulse_wd``) loads a pre-built 0.85 Msun white dwarf model and evolves the model until the effective temperature drops below log10(Teff/K) = 4. During the evolution the ``run_star_extras.f90`` calls GYRE, processes the GYRE output, and searches for a non-radial g-mode frequency of 2.2e-3 Hz (period of 454 s). Close matches to this target frequency are reported in the terminal:

.. code-block:: console

   model   order           freq (Hz)             P (sec)             P (day)
    1292      -8          0.2204E-02              0.0000              0.0000
 matched target frequency
    1292      -7          0.2374E-02              0.0000              0.0000

 ...

   model   order           freq (Hz)             P (sec)             P (day)
    1296      -8          0.2163E-02              0.0000              0.0000
 matched target frequency
    1296      -7          0.2348E-02              0.0000              0.0000

Last-Updated: 11Jun2021 (MESA 5be9e57) by fxt.


