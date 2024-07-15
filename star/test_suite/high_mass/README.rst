.. _high_mass:

*********
high_mass
*********

This test case checks the evolution of a 300 Msun, Z = 1e-4 metallicity, model through core hydrogen depletion.

This test case has 1 part. Click to see a larger version of a plot.

* Part 1 (``inlist_high_mass``) loads a 300 Msun, Z=1e-4 metallicity, zero-age main sequence model and evolves it until the central hydrogen mass fraction drops belw 0.05:

.. image:: ../../../star/test_suite/high_mass/docs/abund_000488.png
   :width: 100%

.. image:: ../../../star/test_suite/high_mass/docs/hr000488.png
   :width: 100%



pgstar commands used for the plots above:


.. literalinclude:: ../../../star/test_suite/high_mass/inlist_pgstar
  :language: console


Last-Updated: 15Jul2024 (MESA 6796e09f) by pmocz


