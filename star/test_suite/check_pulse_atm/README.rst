.. _check_pulse_atm:

***************
check_pulse_atm
***************

.. |Ttau| replace:: :math:`T(\tau)`
.. |tau| replace:: :math:`\tau`

This test checks that the atmosphere structure written to the
pulsation output matches the |Ttau| relation specified by
``atm_T_tau_relation``.  It is based on the ``T_tau_gradr`` test case
but checks reconstructed atmosphere in the pulsation data file, rather
than the structure of the interior model, which in this case doesn't
include the atmosphere.

The test changes ``atm_T_tau_relation`` every 5 steps and changes ``atm_T_tau_opacity`` every 20 steps to cycle through
all the available options.  Every 5 steps, it

1. saves the pulsation output,
2. reads it back in,
3. computes the root-mean-squared differences (rms, stored as ``T_rms``) between the pulsation data's temperature profile and the target |Ttau| relation in layers with |tau| < 0.1, and
4. compares this to a target value set by ``x_ctrl(1)``.
   
If new ``atm_T_tau_relation`` or ``atm_T_tau_opacity`` options are
added to MESA, they must be added to this test case by hand.  That is,
the implementation does not automatically track all the available
options.

If the test fails because ``T_rms`` is slightly larger than the tolerance,
there are two possibly benign explanations.

1. The interpolation error in ``T_face`` contributes too much to ``chi2``.
   The tolerance can be increased or the resolution of the atmosphere integration increased.
2. There are layers included in the sum that are becoming convective,
   in which case the temperature gradient won't follow the (radiative)
   |Ttau| relation.  The sum can be restricted to smaller optical
   depths.
   
If the test fails because ``T_rms`` is much larger (orders of magnitude
larger) than the tolerance, then there might be a bug
in the implementation of the integration of the atmosphere and how it's
written to the pulsation output.

Most options are deliberately left at their default values because
they shouldn't influence the test's result.

The test currently only checks the FGONG output.
It should be extended to cover the other pulsation data formats, too.

Last-Updated: 2022-03-09 (commit 97a1c26) by Warrick Ball
