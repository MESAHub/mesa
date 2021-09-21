.. _16M_conv_premix:

***************
16M_conv_premix
***************

This test suite example re-creates the 16-solar mass main-sequence
evolution with the inclusion of convective premixing (using the Ledoux
criterion), as detailed in Section 5.3 of the MESA V instrument paper
(Paxton et al 2019).

This test case has two parts. Click to see a larger view of a plot.

* Part 1 (``inlist_start``) creates a 16 Msun pre-main-sequence model and evolves it for 10 time steps.

* Part 2 (``inlist_16M_conv_premix``) continues the evolution until core hydrogen depletion (mass fraction h1_center < 1e-6).

A Kippenhahn diagram shows the evolution of a retreating convective core on the main sequence, the blue region between model numbers 240 and 385.
This can be compared to Figure 41 in the MESA V instrument paper, and the predictive mixing Kippenhahn diagram in :ref:`16M_predictive_mix`.

.. image:: ../../../star/test_suite/16M_conv_premix/docs/kipp_00000388.svg
   :scale: 100%


pgstar commands used for the plots above:

.. code-block:: console

 &pgstar

   file_white_on_black_flag = .true. ! white_on_black flags -- true means white foreground color on black background
   file_device = 'png'            ! png
   file_extension = 'png'

   !file_device = 'vcps'          ! postscript
   !file_extension = 'ps'

    pgstar_interval = 10
    file_digits = 8


   Kipp_win_flag = .true.

   ! window properties
   Kipp_win_width = 12
   Kipp_win_aspect_ratio = 0.75
   Kipp_txt_scale = 0.9
   Kipp_title = ''      

   ! y axis limits
   Kipp_mass_max = 16.0
   Kipp_mass_min = 0 
   Kipp_show_mass_boundaries = .true.

   ! x axis limits
   Kipp_xaxis_name = 'model_number'
   Kipp_xmax = -101              ! maximum step number.  negative means use default.
   Kipp_xmin = 0         ! minimum step number.  negative means use default.

   Kipp_show_mixing = .true.
   Kipp_show_burn = .true.
   Kipp_show_luminosities = .true.

   ! file output
   Kipp_file_flag = .true.
   Kipp_file_dir = 'kipp_png'
   Kipp_file_prefix = 'kipp_'
   Kipp_file_interval = 10     ! output when mod(model_number,file_interval)==0
   Kipp_file_width = 12        ! (inches) negative means use same value as for window
   Kipp_file_aspect_ratio = -1 ! negative means use same value as for window

 / ! end of pgstar namelist




Last-Updated: 27May2021 (MESA ebecc10) by fxt

