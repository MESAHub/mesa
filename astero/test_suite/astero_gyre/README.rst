.. _astero_gyre:

***********
astero_gyre
***********

.. tags:: astero, asteroseismology, gyre, frequency-matching, seismic-constraints, low-mass, main-sequence

This test demonstrates using the ``astero`` module to call GYRE during
the evolution of a 1.2 |Msun|, Z=0.02 model based on the ``HD49385``
setup.

The test sets ``oscillation_code = 'gyre'`` and uses
``run_star_extras`` to sample mode frequencies during the run.  At the
end of the evolution it compares a selected mode frequency with the
target value stored in ``x_ctrl(1)``.

How this test passes
====================

The active test-suite entry expects the output string ``got ok match
for expected frequency``.  If GYRE is not enabled, the extras routine
intentionally emits a matching skip-style message.
