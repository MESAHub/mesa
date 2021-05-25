.. _16M_conv_premix:

***************
16M_conv_premix
***************

This test suite example re-creates the 16-solar mass main-sequence
evolution with the inclusion of convective premixing (using the Ledoux
criterion), as detailed in Section 5.3 of the MESA V instrument paper
(Paxton et al 2019).

This test case has two parts.

* Part 1 (``inlist_start``) creates a 16 Msun pre-main-sequence model and evolves it for 10 time steps.

* Part 2 (``inlist_16M_conv_premix``) continues the evolution until core hydrogen depletion (mass fraction h1_center < 1e-6).

A Kippenhahn diagram shows the evolution of a retreating convective core on the main sequence, the blue region between model numbers 240 and 385.
This can be compared to Figure 41 in the MESA V instrument paper, and the predictive mixing Kippenhahn diagram in :ref:`16M_predictive_mix`.

.. image:: ../../../star/test_suite/16M_conv_premix/docs/kipp_00000388.svg


Last-Updated: 25May2021 (MESA ebecc10) by fxt

