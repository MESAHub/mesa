.. |Ttau| replace:: :math:`T(\tau)`
.. |tau| replace:: :math:`\tau`

This test checks the implementation of the control
``use_T_tau_gradr_factor``, which modifies the radiative gradient so
that regions of low optical depth have a temperature that follows the
|Ttau| relation specified by ``atm_T_tau_relation``.  This is useful
if you'd like to include regions of small optical depth as if they're
part of the interior model.

The test changes ``atm_T_tau_relation`` every 10 steps to cycle through
all the available options.  At each step, it computes the sum of
squared differences (stored as ``chi2``) between the stellar model's
temperature profile and the target |Ttau| relation in layers with |tau|
< 0.1.  If new ``atm_T_tau_relation`` options are added to MESA, they
must be added to this test case by hand.  That is, the implementation
does not automatically track all the available options.

The target |Ttau| relation is stored in the extra profile column
``T_check``, which can be compared to ``T_face`` and *not* ``T``, because
the optical depth ``tau`` is evaluated at the cell faces.

If the test fails because ``chi2`` is slightly larger than the tolerance,
there are two possibly benign explanations.

1. The interpolation error in ``T_face`` contributes too much to ``chi2``.
   The tolerance can be increased or the mesh resolution increased.
2. There are layers included in the sum that are becoming convective,
   in which case the temperature gradient won't follow the (radiative)
   |Ttau| relation.  The sum can be restricted to smaller optical
   depths.
   
If the test fails because ``chi2`` is much larger (orders of magnitude
larger) than the tolerance, then there might be a bug
in the implementation of the ``T_tau_gradr_factor``.

Most options are deliberately left at their default values because
they shouldn't influence the test's result.

Last-Updated: 2020-10-27 (mesa r14711) by Warrick Ball
