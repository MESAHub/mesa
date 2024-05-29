.. _conserve_angular_momentum:

*************************
conserve_angular_momentum
*************************

This test suite example checks angular momentum conservation from the zero age main-sequence to the formation of a helium core in 1.0 Msun, Z=0.02 metallicity, model.

This test case has 2 parts. Click to see a larger view of a plot.

* Part 1 (``inlist_zams``) builds a 1.0 Msun, Z=0.02 metallicity, pre-main sequence model and evolves until the zero age main-sequence.

* Part 2 (``inlist_conserve_J``) imposes a uniform 10 km/s rotation profile, and continues the evolution until the helium core reaches a mass of 0.22 Msun. Throughout this evolution a relative change in the angular momentum (J_init - J)/J_init of order 3e-13 is reported in the terminal:

.. image:: ../../../star/test_suite/conserve_angular_momentum/docs/grid_000455.svg
   :scale: 100%


pgstar commands used for the first 7 plots:

.. code-block:: console

 &pgstar

   file_white_on_black_flag = .true. ! white_on_black flags -- true means white foreground color on black background
   file_device = 'png'            ! png
   file_extension = 'png'

   !file_device = 'vcps'          ! postscript
   !file_extension = 'ps'

    pgstar_interval = 10


 Abundance_win_flag = .false.
 Abundance_win_width = 15
 Abundance_win_aspect_ratio = 0.75
         
 Abundance_title = ''
 Abundance_num_isos_to_show = 6

 Abundance_which_isos_to_show(1)  = 'h1'
 Abundance_which_isos_to_show(2)  = 'he3'
 Abundance_which_isos_to_show(3)  = 'he4'
 Abundance_which_isos_to_show(4)  = 'c12'
 Abundance_which_isos_to_show(5)  = 'n14'
 Abundance_which_isos_to_show(6)  = 'o16'

 num_abundance_line_labels = 4
 Abundance_line_txt_scale_factor = 1.0
 Abundance_legend_max_cnt = 0
 Abundance_legend_txt_scale_factor = 0.6

 Abundance_xaxis_name = 'mass' 
 Abundance_xaxis_reversed = .false.
 Abundance_xmin = 0.0 
 Abundance_xmax = 1.0

 Abundance_log_mass_frac_min = -3.0 
 Abundance_log_mass_frac_max =  0.3

 Abundance_file_flag = .false.
 Abundance_file_dir = 'png'
 Abundance_file_prefix = 'abund_'
 Abundance_file_interval = 10     ! output when mod(model_number,file_interval)==0
 Abundance_file_width = -1        ! (inches) negative means use same value as for window
 Abundance_file_aspect_ratio = -1 ! negative means use same value as for window

 Profile_Panels_win_flag(2) = .false.
 Profile_Panels_num_panels(2) = 1
 Profile_Panels_show_grid(2) = .true.
 Profile_Panels_show_mix_regions_on_xaxis(2) = .true.

 Profile_Panels_win_width(2) = 15
 Profile_Panels_win_aspect_ratio(2) = 1.0
 Profile_Panels_txt_scale(2) = 1.0
 Profile_Panels_title(2) = ''

 Profile_Panels_xaxis_name(2) = 'mass'
 Profile_Panels_xmin(2) = 0.0
 Profile_Panels_xmax(2) = 1.0
 Profile_Panels_xaxis_reversed(2) = .false.

 Profile_Panels_yaxis_name(2, 1) = 'omega'
 Profile_Panels_ymin(2, 1) = 0.0
 Profile_Panels_ymax(2, 1) = 8.0e-4 ! only used if /= -101d0

 Profile_Panels_other_yaxis_name(2, 1) = 'log_j_rot'
 Profile_Panels_other_ymin(2, 1) = 10.0 ! only used if /= -101d0
 Profile_Panels_other_ymax(2, 1) = 17.0 ! only used if /= -101d0

 Profile_Panels_file_flag(2) = .false.
 Profile_Panels_file_dir(2) = 'png'
 Profile_Panels_file_prefix(2) = 'j_omega_'
 Profile_Panels_file_interval(2) = 10     ! output when mod(model_number,file_interval)==0
 Profile_Panels_file_width(2) = -1        ! (inches) negative means use same value as for window
 Profile_Panels_file_aspect_ratio(2) = -1 ! negative means use same value as for window



 Text_Summary_win_flag(1) = .false.
 Text_Summary_win_width(1) = 10
 Text_Summary_win_aspect_ratio(1) = 0.15

 Text_Summary_xleft(1) = 0.01
 Text_Summary_xright(1) = 0.99
 Text_Summary_ybot(1) = 0.0
 Text_Summary_ytop(1) = 1.0
 Text_Summary_txt_scale(1) = 1.0

 Text_Summary_num_rows(1) = 1 ! <= 20
 Text_Summary_num_cols(1) = 3 ! <= 20
 Text_Summary_name(1, :, :) = ''

 Text_Summary_name(1, 1, 1) = 'num_zones'
 Text_Summary_name(1, 1, 2) = 'total_angular_momentum'
 Text_Summary_name(1, 1, 3) = 'surf_avg_v_rot'

 Text_Summary_file_flag(1) = .false.
 Text_Summary_file_dir(1) = 'png'
 Text_Summary_file_prefix(1) = 'text_'
 Text_Summary_file_interval(1) = 10     ! output when mod(model_number,file_interval)==0
 Text_Summary_file_width(1) = -1        ! (inches) negative means use same value as for window
 Text_Summary_file_aspect_ratio(1) = -1 ! negative means use same value as for window


 Grid_win_flag(1) = .true.
 Grid_win_width(1) = 15
 Grid_win_aspect_ratio(1) = 0.6

 Grid_plot_name(1, :) = ''
 Grid_plot_row(1, :) = 1           ! number from 1 at top
 Grid_plot_rowspan(1, :) = 1       ! plot spans this number of rows
 Grid_plot_col(1, :) =  1          ! number from 1 at left
 Grid_plot_colspan(1, :) = 1       ! plot spans this number of columns
 Grid_plot_pad_left(1, :) = 0.0    ! fraction of full window width for padding on left
 Grid_plot_pad_right(1, :) = 0.0   ! fraction of full window width for padding on right
 Grid_plot_pad_top(1, :) = 0.0     ! fraction of full window height for padding at top
 Grid_plot_pad_bot(1, :) = 0.0     ! fraction of full window height for padding at bottom
 Grid_txt_scale_factor(1, :) = 0.7 ! multiply txt_scale for subplot by this


 Grid_num_cols(1) = 6 ! divide plotting region into this many equal width cols
 Grid_num_rows(1) = 2 ! divide plotting region into this many equal height rows
 Grid_num_plots(1) = 10 ! <= 10

 Grid_title(1) = 'inlist_conserve_J'

    pgstar_show_model_number = .true.
    pgstar_model_scale = 1.0
    pgstar_model_lw = 3
    pgstar_model_disp = 2.0
    pgstar_model_coord = 0.99
    pgstar_model_fjust = 1.0

    pgstar_show_age = .true.
    pgstar_age_scale = 1.0
    pgstar_age_lw = 3
    pgstar_age_disp = 2.0
    pgstar_age_coord = -0.10
    pgstar_age_fjust = 0.0



 Grid_plot_name(1, 1) = 'Text_Summary1'
 Grid_plot_row(1, 1) = 1           ! number from 1 at top
 Grid_plot_rowspan(1, 1) = 1       ! plot spans this number of rows
 Grid_plot_col(1, 1) =  1          ! number from 1 at left
 Grid_plot_colspan(1, 1) = 6       ! plot spans this number of columns

 Grid_plot_pad_left(1, 1) = -0.06    ! fraction of full window width for padding on left
 Grid_plot_pad_right(1, 1) = 0.05   ! fraction of full window width for padding on right
 Grid_plot_pad_top(1, 1) = -0.02     ! fraction of full window height for padding at top
 Grid_plot_pad_bot(1, 1) = 0.39     ! fraction of full window height for padding at bottom
 Grid_txt_scale_factor(1, 1) = 1.2 ! multiply txt_scale for subplot by this


 Grid_plot_name(1, 2) = 'Abundance'
 Grid_plot_row(1, 2) = 1           ! number from 1 at top
 Grid_plot_rowspan(1, 2) = 2       ! plot spans this number of rows
 Grid_plot_col(1, 2) =  1          ! number from 1 at left
 Grid_plot_colspan(1, 2) = 3       ! plot spans this number of columns

 Grid_plot_pad_left(1, 2) = -0.05    ! fraction of full window width for padding on left
 Grid_plot_pad_right(1, 2) = 0.10   ! fraction of full window width for padding on right
 Grid_plot_pad_top(1, 2) = 0.03     ! fraction of full window height for padding at top
 Grid_plot_pad_bot(1, 2) = 0.03     ! fraction of full window height for padding at bottom
 Grid_txt_scale_factor(1, 2) = 0.7 ! multiply txt_scale for subplot by this


 Grid_plot_name(1, 3) = 'Profile_Panels2'
 Grid_plot_row(1, 3) = 1          ! number from 1 at top
 Grid_plot_rowspan(1, 3) = 2       ! plot spans this number of rows
 Grid_plot_col(1, 3) =  5          ! Number from 1 at left
 Grid_plot_colspan(1, 3) = 3       ! plot spans this number of columns

 Grid_plot_pad_left(1, 3) = -0.15    ! fraction of full window width for padding on left
 Grid_plot_pad_right(1, 3) = 0.20   ! fraction of full window width for padding on right
 Grid_plot_pad_top(1, 3) = 0.03     ! fraction of full window height for padding at top
 Grid_plot_pad_bot(1, 3) = 0.03     ! fraction of full window height for padding at bottom
 Grid_txt_scale_factor(1, 3) = 0.7 ! multiply txt_scale for subplot by this


 Grid_file_flag(1) = .true.
 Grid_file_dir(1) = 'png'
 Grid_file_prefix(1) = 'grid_'
 Grid_file_interval(1) = 10     ! output when mod(model_number,Grid1_file_interval)==0
 Grid_file_width(1) = -1       ! (inches) negative means use same value as for window
 Grid_file_aspect_ratio(1) = -1 ! negative means use same value as for window




 / ! end of pgstar namelist


Last-Updated: 04Jun2021 (MESA 5be9e57) by fxt

