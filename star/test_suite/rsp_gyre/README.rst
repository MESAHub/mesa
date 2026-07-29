.. _rsp_gyre:

********
rsp_gyre
********

.. tags:: star, rsp, gyre, asteroseismology, radial-pulsation, linear-analysis, period-search, colors

This test exercises the interaction between RSP model construction,
linear pulsation analysis, and GYRE.

The case creates an RSP model, enables ``use_other_RSP_linear_analysis``,
and compares the selected pulsation period against the target stored in
``x_ctrl(1)``.  It also enables the colors module for synthetic
photometry output.

How this test passes
====================

The active test-suite entry expects the output string ``good match for
period``.  The case has ``.ignore_checksum``, so checksum reporting is
intentionally suppressed.
