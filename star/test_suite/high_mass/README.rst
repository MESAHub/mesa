.. _high_mass:

*********
high_mass
*********

This test case checks the evolution of a 300 |Msun|, Z = 1e-4 metallicity, model through core hydrogen depletion.

This test case has 1 part. Click to see a larger version of a plot.

* Part 1 (``inlist_high_mass``) loads a 300 |Msun|, Z=1e-4 metallicity, zero-age main sequence model and evolves it until the central hydrogen mass fraction drops belw 0.05:

.. image:: ../../../star/test_suite/high_mass/docs/abund_000488.svg
   :width: 100%

.. image:: ../../../star/test_suite/high_mass/docs/hr000488.svg
   :width: 100%


pgstar commands used for the plots above:

.. literalinclude:: ../../../star/test_suite/high_mass/inlist_pgstar
   :language: console


Last-Updated: 12Jun2021 (MESA 5be9e57) by fxt.

Last-Run: 22Oct2024 (MESA 9b2017ca) by pmocz on C916PXT6XW in 84 seconds using 8 threads.
