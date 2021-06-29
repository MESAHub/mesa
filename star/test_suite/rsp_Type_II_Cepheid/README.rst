.. _rsp_Type_II_Cepheid:

*******************
rsp_Type_II_Cepheid
*******************

This test case checks the non-linear pulsation evolution of a 0.55 Msun, Teff = 6410 K, L = 136 Lsun, Z = 0.0001 metallicity model - 
type-II Cepheid of BL Her type based on |Smolec14|.

This test case has 1 part. Click to see a larger version of a plot.

* Part 1 (``inlist_rsp_Type_II_Cepheid``) creates the initial 0.65 Msun, 0.55 Msun, Teff = 6410 K, L = 136 Lsun, Z = 0.0001 metallicity model, and writes the results of conducting a linear nonadiabatic stability analysis to the LOGS directory (see Section 2.2 of |MESA V| for details). The evolution with RSP then begins. The convective parameters are such that the model should show deterministic chaos - lots of struggling between low and high amplitudes. The growth rates are large, so the attractor is approached only after ~20-40 cycles. After 10 periods, the ``run_star_extras.f90`` checks if the energy conservation is less than 1e-5 and if fundamental period is within 1% of the expected 0.71262 day period. If these values are within bounds, then a message is written to the terminal and the run terminates:

.. code-block:: console

 rel_run_E_err  -1.0724594936404145E-007
 good match for period   1.7372426555333260        1.7372399999999999


.. image:: ../../../star/test_suite/rsp_Type_II_Cepheid/docs/grid_0014680.svg
   :width: 100%

pgstar commands, in addition to those in ``inlist_rsp_pgstar_default`` which currently must be copied from one of the other ``rsp_*`` test suite cases, used for the plot above:

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
       Grid2_title = 'rsp_Type_II_Cepheid'
       Grid2_txt_scale_factor(:) = 1.0

        Grid2_file_flag = .true.
        Grid2_file_dir = 'pgstar_out'
        Grid2_file_prefix = 'grid_'
        Grid2_file_interval = 10000
        Grid2_file_width = -1
        Grid2_file_aspect_ratio = -1

 / ! end of pgstar namelist


.. |Smolec14| replace:: `Smolec and Moskalik (2014) <https://ui.adsabs.harvard.edu/abs/2014MNRAS.441..101S/abstract>`__


Last-Updated: 27Jun2021 (MESA e2acbc2) by fxt.
