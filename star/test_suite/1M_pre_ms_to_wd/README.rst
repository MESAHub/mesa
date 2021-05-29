.. _1M_pre_ms_to_wd:

***************
1M_pre_ms_to_wd
***************

This test case checks the evolution of a 1 Msun, Z=0.02 metallicity from the pre-main sequence to a white dwarf. 

This test case has six parts. Click to see a larger view of a plot.

* Part 1 (``inlist_start_header``) builds a 1 Msun, Z=0.02 metallicity pre-main-sequence model.

* Part 2 (``inlist_to_end_core_h_burn``) continues the evolution until core hydrogen depletion (mass fraction h1_center < 1e-4).

.. image:: ../../../star/test_suite/1M_pre_ms_to_wd/docs/grid6_00000295.svg
   :scale: 100%

* Part 3 (``inlist_to_start_he_core_flash_header``) continues the evolution until the onset of core helium ignition (power from he burning > 10).

.. image:: ../../../star/test_suite/1M_pre_ms_to_wd/docs/grid6_00010474.svg
   :scale: 100%

* Part 4 (``inlist_to_end_core_he_burn``) continues the evolution until core helium depletion (mass fraction he4_center < 1e-4).

.. image:: ../../../star/test_suite/1M_pre_ms_to_wd/docs/grid6_00012472.svg
   :scale: 100%

* Part 5 (``inlist_to_end_agb``) continues the evolution through the thermal pulses until the end of the AGB phase of evolution (hydrogen-rich envelope mass < 1e-2 Msun).

.. image:: ../../../star/test_suite/1M_pre_ms_to_wd/docs/grid6_00013000.svg
   :scale: 100%


* Part 6 (``inlist_to_wd``) continues the evolution until the luminosity of the cooling white dwarf reaches L < 0.1 Lsun.

.. image:: ../../../star/test_suite/1M_pre_ms_to_wd/docs/grid6_00013356.svg
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

      Grid6_win_flag = .true.
      Grid6_win_width = 15
         
      Grid6_file_flag = .true.
      Grid6_file_dir = 'png'
      Grid6_file_prefix = 'grid6_'
      Grid6_file_interval = 1000 ! output when mod(model_number,Grid6_file_interval)==0
      Grid6_file_width = 15 ! (inches) negative means use same value as for window

      Abundance_log_mass_frac_min = -4.0 
      Abundance_legend_max_cnt = 0
      num_abundance_line_labels = 4
      Abundance_title = ''

      HR_title = ''

      TRho_title = '' 

      Summary_Burn_title = '' 
      Summary_Burn_xaxis_name = 'mass' 
      Summary_Burn_xaxis_reversed = .false.
      Summary_Burn_xmin = 0.00 ! -101d0 ! only used if /= -101d0
      Summary_Burn_xmax = 1.0  ! only used if /= -101d0


 / ! end of pgstar namelist




Last-Updated: 28May2021 (MESA ebecc10) by fxt

