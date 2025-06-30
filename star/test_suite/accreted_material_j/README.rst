.. _accreted_material_j:

*******************
accreted_material_j
*******************

This test suite example checks the accretion of material and angular momentum onto a 20 |Msun| model.


This test case has two parts. Click to see a larger view of a plot.

* Part 1 (``inlist_zams``) creates a 20 |Msun|, Z=0.02 metallicity, main-sequence model.

* Part 2 (``inlist_accreted_material_j``) continues the evolution by first applying uniform rotation Omega/Omega_crit = 0.1, and then accreting material with the same composition as the outermost cell at a rate of 0.001 |Msun|/year with an angular momentum of 0.1 Keplerian. The model terminates when the mass reaches 25 |Msun|:

.. image:: ../../../star/test_suite/accreted_material_j/docs/track1_000452.svg
   :width: 100%


pgstar commands used for the plots above:

.. literalinclude:: ../../../star/test_suite/accreted_material_j/inlist_pgstar
  :language: console
  :start-at: &pgstar
  :end-at: ! end of pgstar namelist


Last-Updated: 30May2021 (MESA e37f76f) by fxt

Last-Run: 22Oct2024 (MESA 9b2017ca) by pmocz on C916PXT6XW in 178 seconds using 8 threads.
