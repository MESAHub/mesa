.. _rsp_Delta_Scuti:

***************
rsp_Delta_Scuti
***************

This test case checks the non-linear pulsation evolution of a 2 Msun, Teff = 6900 K, L = 30 Lsun, Z = 0.02 metallicity - 
a double-mode delta Scuti variable leaving the main-sequence phase originally contributed by Alfred Gautschy.

This test case has 1 part. Click to see a larger version of a plot.

* Part 1 (``inlist_rsp_Delta_Scuti``) creates the initial 2 Msun, Teff = 6900 K, L = 30 Lsun, Z = 0.02 metallicity model, and writes the results of conducting a linear nonadiabatic stability analysis to the LOGS directory (see Section 2.2 of |MESA V| for details). The evolution with RSP then begins. After 10 periods, the ``run_star_extras.f90`` checks if the energy conservation is less than 1e-5 and if fundamental period is within 1% of the expected 0.14715 day period. If these values are within bounds, then a message is written to the terminal and the run terminates:

.. code-block:: console

 rel_run_E_err  -7.2889211517560327E-011
 good match for period  0.14674260955972285       0.14715000000000000     

This test case also exercises the ability of the RSP module to write a file suitable for post-processing 
with gnuplot's pm3d. See the ``RSP_map_*`` inlist controls, with the output file being put in the ``LOGS`` directory.

.. image:: ../../../star/test_suite/rsp_Delta_Scuti/docs/grid_0007973.svg
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
       Grid2_title = 'rsp_Delta_Scuti'
       Grid2_txt_scale_factor(:) = 1.0

        Grid2_file_flag = .true.
        Grid2_file_dir = 'pgstar_out'
        Grid2_file_prefix = 'grid_'
        Grid2_file_interval = 10000
        Grid2_file_width = -1
        Grid2_file_aspect_ratio = -1

 / ! end of pgstar namelist

Last-Updated: 27Jun2021 (MESA e2acbc2) by fxt.
