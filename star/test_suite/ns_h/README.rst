.. _ns_h:

****
ns_h
****

This test case shows an example of steady hydrogen burning within a neutron star envelope.

This test case has 2 parts. Click to see a larger version of a plot.

* Part 1 (``inlist_add_he_layer``) first loads a prebuilt neutron star envelope model ``ns_env.mod``, see the :ref:`make_env` test case. Helium is then accreted at a rate of 1e-9 Msun/year for 0.5 seconds to create a thin helium layer:

.. image:: ../../../star/test_suite/ns_h/docs/grid_000108.png
   :width: 100%

* Part 2 (``inlist_to_steady_h_burn``) continues the evolution by accreting a hydrogen-rich solar mixture at a rate of 1e-11 Msun/year. Hydrogen burning eventually reaches a steady location, ``logxm`` = log10((Mstar - m) / Msun) = -12.8, within the neutron star envelope after 2e5 seconds of evolution:

.. image:: ../../../star/test_suite/ns_h/docs/grid_000381.png
   :width: 100%


pgstar commands used for the plot above:

.. code-block:: console

 &pgstar

  file_white_on_black_flag = .true. ! white_on_black flags -- true means white foreground color on black background
  !file_device = 'png'            ! png
  !file_extension = 'png'

  file_device = 'vcps'          ! postscript
  file_extension = 'ps'

  pgstar_interval = 10

  pgstar_grid_title_disp = 1.8

  Abundance_xaxis_name = 'logxm' 
  Abundance_xaxis_reversed = .true.
  Abundance_xmin = -16.0
  Abundance_xmax = -12.0

  Power_xaxis_name = 'logxm'
  Power_xaxis_reversed = .true
  Power_xmin = -16.0
  Power_xmax = -12.0

  Profile_Panels_title(1) = ''
  Profile_Panels_num_panels(1) = 2

  Profile_Panels_xaxis_name(1) = 'logxm'
  Profile_Panels_xaxis_reversed(1) = .true.
  Profile_Panels_xmin(1) = -101d0
  Profile_Panels_xmax(1) = -101d0

  Profile_Panels_yaxis_name(1, 1) = 'logRho'
  Profile_Panels_yaxis_name(1, 2) = 'logT'
  Profile_Panels_ymin(1, 1) = -101

  Profile_Panels_other_yaxis_name(1, 1) = 'logP'
  Profile_Panels_other_yaxis_name(1, 2) = 'entropy'
  Profile_Panels_other_ymin(1, 1) = -101

  Profile_Panels_title(2) = ''
  Profile_Panels_num_panels(2) = 2

  Profile_Panels_xaxis_name(2) = 'logxm'
  Profile_Panels_xaxis_reversed(2) = .true.
  Profile_Panels_xmin(2) = -101d0
  Profile_Panels_xmax(2) = -101d0

  Profile_Panels_yaxis_name(2, 1) = 'luminosity'
  Profile_Panels_yaxis_name(2, 2) = 'net_nuclear_energy'
  Profile_Panels_ymin(2, 1) = -101
  Profile_Panels_ymin(2, 2) = -101

  Profile_Panels_other_yaxis_name(2, 1) = 'opacity'
  Profile_Panels_other_yaxis_name(2, 2) = 'eps_nuc_neu_total'
  Profile_Panels_other_ymin(2, 1) = -101
  Profile_Panels_other_ymin(2, 2) = -101

  Text_Summary_txt_scale(1) = 5.5

  Text_Summary_num_rows(1) = 5
  Text_Summary_num_cols(1) = 3
  Text_Summary1_name(1, 1) = 'model_number'
  Text_Summary1_name(2, 1) = 'star_age_sec'
  Text_Summary1_name(3, 1) = 'time_step_sec'
  Text_Summary1_name(4, 1) = 'log_rel_run_E_err'
  Text_Summary1_name(5, 1) = 'total_energy'
  Text_Summary1_name(1, 2) = 'envelope_mass'
  Text_Summary1_name(2, 2) = 'log_abs_mdot'
  Text_Summary1_name(3, 2) = 'log_xmstar'
  Text_Summary1_name(4, 2) = 'm_center'
  Text_Summary1_name(5, 2) = 'r_center_km'
  Text_Summary1_name(1, 3) = 'num_zones'
  Text_Summary1_name(2, 3) = 'num_iters'
  Text_Summary1_name(3, 3) = 'num_retries'
  Text_Summary1_name(4, 3) = ' '
  Text_Summary1_name(5, 3) = ' '

  Grid_title(2) = 'ns_h'
  Grid_plot_name(2, 1) = 'Profile_Panels1'
  Grid_plot_name(2, 2) = 'Text_Summary1'
  Grid_plot_name(2, 3) = 'Abundance'
  Grid_plot_name(2, 4) = 'Power'
  Grid_plot_name(2, 5) = 'Profile_Panels2'
  Grid_plot_row(2, 1) = 1
  Grid_plot_row(2, 2) = 7
  Grid_plot_row(2, 3) = 1
  Grid_plot_row(2, 4) = 5
  Grid_plot_row(2, 5) = 1
  Grid_plot_rowspan(2, 1) = 6
  Grid_plot_rowspan(2, 2) = 2
  Grid_plot_rowspan(2, 3) = 4
  Grid_plot_rowspan(2, 4) = 4
  Grid_plot_rowspan(2, 5) = 6
  Grid_plot_col(2, 1) = 1
  Grid_plot_col(2, 2) = 1
  Grid_plot_col(2, 3) = 5
  Grid_plot_col(2, 4) = 5
  Grid_plot_col(2, 5) = 3
  Grid_plot_colspan(2, 1) = 2
  Grid_plot_colspan(2, 2) = 4
  Grid_plot_colspan(2, 3) = 3
  Grid_plot_colspan(2, 4) = 3
  Grid_plot_colspan(2, 5) = 2
  Grid_plot_pad_left(2, 1) = -0.02
  Grid_plot_pad_left(2, 2) = -0.08
  Grid_plot_pad_left(2, 3) = 0.14
  Grid_plot_pad_left(2, 4) = 0.14
  Grid_plot_pad_left(2, 5) = 0.06
  Grid_plot_pad_right(2, 1) = 0.07
  Grid_plot_pad_right(2, 2) = -0.12
  Grid_plot_pad_right(2, 3) = 0
  Grid_plot_pad_right(2, 4) = 0
  Grid_plot_pad_right(2, 5) = -0.01
  Grid_plot_pad_top(2, 1) = 0
  Grid_plot_pad_top(2, 2) = 0.08
  Grid_plot_pad_top(2, 3) = 0
  Grid_plot_pad_top(2, 4) = 0.06
  Grid_plot_pad_top(2, 5) = 0
  Grid_plot_pad_bot(2, 1) = 0
  Grid_plot_pad_bot(2, 2) = -0.04
  Grid_plot_pad_bot(2, 3) = 0.09
  Grid_plot_pad_bot(2, 4) = 0.03
  Grid_plot_pad_bot(2, 5) = 0
  Grid_txt_scale_factor(2, 1) = 0.65
  Grid_txt_scale_factor(2, 2) = 0.19
  Grid_txt_scale_factor(2, 3) = 0.65
  Grid_txt_scale_factor(2, 4) = 0.65
  Grid_txt_scale_factor(2, 5) = 0.65

  Grid_num_cols(2) = 7
  Grid_num_rows(2) = 8
  Grid_num_plots(2) = 5
  
  Grid_win_flag(2) = .true.
  Grid_win_width(2) = 16
  Grid_win_aspect_ratio(2) = 0.6
  
  Grid_file_flag(2) = .true.
  Grid_file_dir(2) = 'pgstar_out'
  Grid_file_prefix(2) = 'grid_'
  Grid_file_interval(2) = 10000
  Grid_file_width(2) = 20
  Grid_file_aspect_ratio(2) = -1

 / ! end of pgstar namelist


Last-Updated: 21Jun2021 (MESA e2acbc2) by fxt.
