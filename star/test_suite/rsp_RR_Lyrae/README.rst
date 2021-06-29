.. _rsp_RR_Lyrae:

************
rsp_RR_Lyrae
************

This test case checks the non-linear pulsation evolution of a 0.65 Msun, Teff = 6500 K, L = 60 Lsun, Z = 0.004 metallicity - 
a long-period RR Lyrae model contributed by Radek Smolec.

This test case has 1 part. Click to see a larger version of a plot.

* Part 1 (``inlist_rsp_rsp_RR_Lyrae``) creates the initial 0.65 Msun, Teff = 6500 K, L = 60 Lsun, Z = 0.004 metallicity model, and writes the results of conducting a linear nonadiabatic stability analysis to the LOGS directory (see Section 2.2 of |MESA V| for details). The evolution with RSP then begins. After 10 periods, the ``run_star_extras.f90`` checks if the energy conservation is less than 1e-5 and if fundamental period is within 1% of the expected 0.71262 day period. If these values are within bounds, then a message is written to the terminal and the run terminates:

.. code-block:: console

 rel_run_E_err   4.3566232972447423E-011
 good match for period  0.71063341179217354       0.71262000000000003


.. image:: ../../../star/test_suite/rsp_RR_Lyrae/docs/grid_0007220.svg
   :width: 100%

pgstar commands, in addition to those in ``inlist_rsp_pgstar_default``, used for the plot above:

.. code-block:: console

 &pgstar

  file_white_on_black_flag = .true. ! white_on_black flags -- true means white foreground color on black background
  !file_device = 'png'            ! png
  !file_extension = 'png'

  file_device = 'vcps'          ! postscript
  file_extension = 'ps'

  pgstar_interval = 100

      pgstar_age_scale = 0.8
      pgstar_age_lw = 3
      pgstar_age_disp = 3.9
      pgstar_age_coord = -0.11
      pgstar_age_fjust = 0.0

      pgstar_model_disp = 3.9

      History_Panels2_txt_scale = 0.7
      Profile_Panels2_txt_scale = 0.6
      logL_R_txt_scale = 0.7
      logL_v_txt_scale = 0.7
      logL_Teff_txt_scale = 0.7

       Grid2_win_flag = .true.
       Grid2_win_width = 12
       Grid2_title = 'rsp_RR_Lyrae'
       Grid2_txt_scale_factor(:) = 1.0

        Grid2_file_flag = .true.
        Grid2_file_dir = 'pgstar_out'
        Grid2_file_prefix = 'grid_'
        Grid2_file_interval = 10000
        Grid2_file_width = -1
        Grid2_file_aspect_ratio = -1

 / ! end of pgstar namelist

Last-Updated: 27Jun2021 (MESA e2acbc2) by fxt.
