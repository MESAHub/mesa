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

 Profile_Panels2_win_flag = .false.
 Profile_Panels2_num_panels = 1
 Profile_Panels2_show_grid = .true.
 Profile_Panels2_show_mix_regions_on_xaxis = .true.

 Profile_Panels2_win_width = 15
 Profile_Panels2_win_aspect_ratio = 1.0
 Profile_Panels2_txt_scale = 1.0
 Profile_Panels2_title = ''      

 Profile_Panels2_xaxis_name = 'mass'
 Profile_Panels2_xmin = 0.0 
 Profile_Panels2_xmax = 1.0
 Profile_Panels2_xaxis_reversed = .false.

 Profile_Panels2_yaxis_name(1) = 'omega'   
 Profile_Panels2_ymin(1) = 0.0
 Profile_Panels2_ymax(1) = 8.0e-4 ! only used if /= -101d0

 Profile_Panels2_other_yaxis_name(1) = 'log_j_rot'   
 Profile_Panels2_other_ymin(1) = 10.0 ! only used if /= -101d0
 Profile_Panels2_other_ymax(1) = 17.0 ! only used if /= -101d0

 Profile_Panels2_file_flag = .false.
 Profile_Panels2_file_dir = 'png'
 Profile_Panels2_file_prefix = 'j_omega_'
 Profile_Panels2_file_interval = 10     ! output when mod(model_number,file_interval)==0
 Profile_Panels2_file_width = -1        ! (inches) negative means use same value as for window
 Profile_Panels2_file_aspect_ratio = -1 ! negative means use same value as for window



 Text_Summary1_win_flag = .false.
 Text_Summary1_win_width = 10
 Text_Summary1_win_aspect_ratio = 0.15

 Text_Summary1_xleft = 0.01
 Text_Summary1_xright = 0.99
 Text_Summary1_ybot = 0.0
 Text_Summary1_ytop = 1.0
 Text_Summary1_txt_scale = 1.0

 Text_Summary1_num_rows = 1 ! <= 20
 Text_Summary1_num_cols = 3 ! <= 20
 Text_Summary1_name(:,:) = ''

 Text_Summary1_name(1,1) = 'num_zones'
 Text_Summary1_name(1,2) = 'total_angular_momentum'
 Text_Summary1_name(1,3) = 'surf_avg_v_rot'

 Text_Summary1_file_flag = .false.
 Text_Summary1_file_dir = 'png'
 Text_Summary1_file_prefix = 'text_'
 Text_Summary1_file_interval = 10     ! output when mod(model_number,file_interval)==0
 Text_Summary1_file_width = -1        ! (inches) negative means use same value as for window
 Text_Summary1_file_aspect_ratio = -1 ! negative means use same value as for window


 Grid1_win_flag = .true.
 Grid1_win_width = 15
 Grid1_win_aspect_ratio = 0.6

 Grid1_plot_name(:) = ''
 Grid1_plot_row(:) = 1           ! number from 1 at top
 Grid1_plot_rowspan(:) = 1       ! plot spans this number of rows
 Grid1_plot_col(:) =  1          ! number from 1 at left
 Grid1_plot_colspan(:) = 1       ! plot spans this number of columns
 Grid1_plot_pad_left(:) = 0.0    ! fraction of full window width for padding on left
 Grid1_plot_pad_right(:) = 0.0   ! fraction of full window width for padding on right
 Grid1_plot_pad_top(:) = 0.0     ! fraction of full window height for padding at top
 Grid1_plot_pad_bot(:) = 0.0     ! fraction of full window height for padding at bottom
 Grid1_txt_scale_factor(:) = 0.7 ! multiply txt_scale for subplot by this


 Grid1_num_cols = 6 ! divide plotting region into this many equal width cols
 Grid1_num_rows = 2 ! divide plotting region into this many equal height rows
 Grid1_num_plots = 10 ! <= 10

 Grid1_title = 'inlist_conserve_J'

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



 Grid1_plot_name(1) = 'Text_Summary1'
 Grid1_plot_row(1) = 1           ! number from 1 at top
 Grid1_plot_rowspan(1) = 1       ! plot spans this number of rows
 Grid1_plot_col(1) =  1          ! number from 1 at left
 Grid1_plot_colspan(1) = 6       ! plot spans this number of columns

 Grid1_plot_pad_left(1) = -0.06    ! fraction of full window width for padding on left
 Grid1_plot_pad_right(1) = 0.05   ! fraction of full window width for padding on right
 Grid1_plot_pad_top(1) = -0.02     ! fraction of full window height for padding at top
 Grid1_plot_pad_bot(1) = 0.39     ! fraction of full window height for padding at bottom
 Grid1_txt_scale_factor(1) = 1.2 ! multiply txt_scale for subplot by this


 Grid1_plot_name(2) = 'Abundance'
 Grid1_plot_row(2) = 1           ! number from 1 at top
 Grid1_plot_rowspan(2) = 2       ! plot spans this number of rows
 Grid1_plot_col(2) =  1          ! number from 1 at left
 Grid1_plot_colspan(2) = 3       ! plot spans this number of columns

 Grid1_plot_pad_left(2) = -0.05    ! fraction of full window width for padding on left
 Grid1_plot_pad_right(2) = 0.10   ! fraction of full window width for padding on right
 Grid1_plot_pad_top(2) = 0.03     ! fraction of full window height for padding at top
 Grid1_plot_pad_bot(2) = 0.03     ! fraction of full window height for padding at bottom
 Grid1_txt_scale_factor(2) = 0.7 ! multiply txt_scale for subplot by this


 Grid1_plot_name(3) = 'Profile_Panels2'
 Grid1_plot_row(3) = 1          ! number from 1 at top
 Grid1_plot_rowspan(3) = 2       ! plot spans this number of rows
 Grid1_plot_col(3) =  5          ! Number from 1 at left
 Grid1_plot_colspan(3) = 3       ! plot spans this number of columns

 Grid1_plot_pad_left(3) = -0.15    ! fraction of full window width for padding on left
 Grid1_plot_pad_right(3) = 0.20   ! fraction of full window width for padding on right
 Grid1_plot_pad_top(3) = 0.03     ! fraction of full window height for padding at top
 Grid1_plot_pad_bot(3) = 0.03     ! fraction of full window height for padding at bottom
 Grid1_txt_scale_factor(3) = 0.7 ! multiply txt_scale for subplot by this


 Grid1_file_flag = .true.
 Grid1_file_dir = 'png'
 Grid1_file_prefix = 'grid_'
 Grid1_file_interval = 10     ! output when mod(model_number,Grid1_file_interval)==0
 Grid1_file_width = -1       ! (inches) negative means use same value as for window
 Grid1_file_aspect_ratio = -1 ! negative means use same value as for window




 / ! end of pgstar namelist


Last-Updated: 04Jun2021 (MESA 5be9e57) by fxt

