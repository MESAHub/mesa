!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! pgbinary controls

type (pgbinary_win_file_data) :: pgbinary_win_file_ptr(num_pgbinary_plots)

integer :: pgbinary_interval
character (len = 32) :: file_device
character (len = 32) :: file_extension
integer :: file_digits

logical :: file_white_on_black_flag, win_white_on_black_flag, &
      pgbinary_show_age, &
      pgbinary_show_age_in_seconds, &
      pgbinary_show_age_in_minutes, &
      pgbinary_show_age_in_hours, &
      pgbinary_show_age_in_days, &
      pgbinary_show_age_in_years, &
      pgbinary_show_log_age_in_years, &
      pgbinary_show_model_number, &
      pgbinary_grid_show_title, &
      pgbinary_show_title, &
      pgbinary_report_writing_files

integer :: pgbinary_lw, pgbinary_box_lw, &
      pgbinary_grid_title_lw, &
      pgbinary_title_lw, &
      pgbinary_age_lw, &
      pgbinary_model_lw, &
      pgbinary_xaxis_label_lw, &
      pgbinary_left_yaxis_label_lw, &
      pgbinary_right_yaxis_label_lw

real :: pgbinary_grid_title_scale, &
      pgbinary_grid_title_disp, &
      pgbinary_grid_title_coord, &
      pgbinary_grid_title_fjust, &
      pgbinary_title_scale, &
      pgbinary_title_disp, &
      pgbinary_title_coord, &
      pgbinary_title_fjust, &
      pgbinary_age_scale, &
      pgbinary_age_disp, &
      pgbinary_age_coord, &
      pgbinary_age_fjust, &
      pgbinary_model_scale, &
      pgbinary_model_disp, &
      pgbinary_model_coord, &
      pgbinary_model_fjust, &
      pgbinary_xaxis_label_disp, &
      pgbinary_xaxis_label_scale, &
      pgbinary_left_yaxis_label_scale, &
      pgbinary_right_yaxis_label_scale, &
      pgbinary_left_yaxis_label_disp, &
      pgbinary_right_yaxis_label_disp, &
      pgbinary_num_scale

logical :: pause
integer :: pause_interval
real :: pgbinary_sleep
logical :: clear_history

logical, dimension(pgbinary_array_length) :: Text_Summary_win_flag, Text_Summary_file_flag
logical, dimension(pgbinary_array_length)  :: do_Text_Summary_win, do_Text_Summary_file
integer, dimension(pgbinary_array_length)  :: id_Text_Summary_win, id_Text_Summary_file, Text_Summary_file_interval
character (len = strlen), dimension(pgbinary_array_length)  :: Text_Summary_title, Text_Summary_file_dir, &
      Text_Summary_file_prefix
integer, dimension(pgbinary_array_length)  :: Text_Summary_num_cols, Text_Summary_num_rows
character (len = 64), dimension(pgbinary_array_length, max_num_rows_Text_Summary, max_num_cols_Text_Summary)  :: &
      Text_Summary_name
real, dimension(pgbinary_array_length)  :: Text_Summary_xleft, Text_Summary_xright, &
      Text_Summary_ybot, Text_Summary_ytop, Text_Summary_txt_scale, &
      Text_Summary_win_width, Text_Summary_win_aspect_ratio, &
      prev_Text_Summary_win_width, prev_Text_Summary_win_ratio, &
      Text_Summary_file_width, Text_Summary_file_aspect_ratio, &
      prev_Text_Summary_file_width, prev_Text_Summary_file_ratio

logical, dimension(pgbinary_array_length) :: History_Track_win_flag, History_Track_file_flag
logical, dimension(pgbinary_array_length) :: do_History_Track_win, do_History_Track_file
logical, dimension(pgbinary_array_length) :: show_History_Track_annotation1, &
      show_History_Track_annotation2, show_History_Track_annotation3
logical, dimension(pgbinary_array_length) :: show_History_Track_target_box, &
      History_Track_reverse_xaxis, History_Track_reverse_yaxis, &
      History_Track_log_xaxis, History_Track_log_yaxis
integer, dimension(pgbinary_array_length) :: id_History_Track_win, id_History_Track_file, History_Track_file_interval
integer, dimension(pgbinary_array_length) :: History_Track_step_min, History_Track_step_max, History_Track_n_sigma
character (len = strlen), dimension(pgbinary_array_length) :: &
      History_Track_file_dir, History_Track_file_prefix, History_Track_title, &
      History_Track_xname, History_Track_xaxis_label, &
      History_Track_yname, History_Track_yaxis_label, History_Track_fname
real, dimension(pgbinary_array_length) :: History_Track_xmin, History_Track_xmax, &
      History_Track_ymin, History_Track_ymax, &
      History_Track_xmargin, History_Track_ymargin, &
      History_Track_dxmin, History_Track_dymin, &
      History_Track_xtarget, History_Track_xsigma, &
      History_Track_ytarget, History_Track_ysigma, &
      History_Track_win_width, History_Track_win_aspect_ratio, &
      History_Track_xleft, History_Track_xright, &
      History_Track_ybot, History_Track_ytop, History_Track_txt_scale, &
      prev_History_Track_win_width, prev_History_Track_win_ratio, &
      History_Track_file_width, History_Track_file_aspect_ratio, &
      prev_History_Track_file_width, prev_History_Track_file_ratio

logical, dimension(pgbinary_array_length) :: History_Panels_win_flag, History_Panels_file_flag
logical, dimension(pgbinary_array_length) :: do_History_Panels_win, do_History_Panels_file
integer, dimension(pgbinary_array_length) :: id_History_Panels_win, id_History_Panels_file, History_Panels_file_interval
character (len = strlen), dimension(pgbinary_array_length) :: History_Panels_file_dir, History_Panels_file_prefix
logical, dimension(pgbinary_array_length) :: History_Panels_xaxis_reversed, History_Panels_xaxis_log
character (len = strlen), dimension(pgbinary_array_length) :: History_Panels_xaxis_name, History_Panels_title
real, dimension(pgbinary_array_length) :: History_Panels_xmin, History_Panels_xmax, History_Panels_dxmin, &
      History_Panels_max_width, History_Panels_xmargin, &
      History_Panels_win_width, History_Panels_win_aspect_ratio, &
      History_Panels_xleft, History_Panels_xright, &
      History_Panels_ybot, History_Panels_ytop, History_Panels_txt_scale, &
      prev_History_Panels_win_width, prev_History_Panels_win_ratio, &
      History_Panels_file_width, History_Panels_file_aspect_ratio, &
      prev_History_Panels_file_width, prev_History_Panels_file_ratio
integer, dimension(pgbinary_array_length) :: History_Panels_num_panels
character (len = strlen), dimension(pgbinary_array_length, max_num_History_Panels) :: &
      History_Panels_yaxis_name, History_Panels_other_yaxis_name, &
      History_Panels_points_name
logical, dimension(pgbinary_array_length, max_num_History_Panels) :: &
      History_Panels_yaxis_reversed, History_Panels_other_yaxis_reversed, &
      History_Panels_yaxis_log, History_Panels_other_yaxis_log
real, dimension(pgbinary_array_length, max_num_History_Panels) :: &
      History_Panels_ymin, History_Panels_other_ymin, &
      History_Panels_ymax, History_Panels_other_ymax, &
      History_Panels_ymargin, History_Panels_other_ymargin, &
      History_Panels_dymin, History_Panels_other_dymin

logical :: History_Panel_points_error_bars
integer :: History_Panel_points_interval, History_Panel_points_marker, &
      History_Panel_points_ci, History_Panel_points_lw
real :: History_Panel_points_ch

logical, dimension(pgbinary_array_length) :: Star_History_track_win_flag, Star_History_track_file_flag
logical, dimension(pgbinary_array_length) :: do_Star_History_track_win, do_Star_History_track_file
logical, dimension(pgbinary_array_length) :: show_Star_History_track_annotation1, &
      show_Star_History_track_annotation2, show_Star_History_track_annotation3
logical, dimension(pgbinary_array_length) :: Star_History_track_reverse_xaxis, Star_History_track_reverse_yaxis, &
      Star_History_track_log_xaxis, Star_History_track_log_yaxis
integer, dimension(pgbinary_array_length) :: id_Star_History_track_win, id_Star_History_track_file, Star_History_track_file_interval
integer, dimension(pgbinary_array_length) :: Star_History_track_step_min, Star_History_track_step_max
character (len = strlen), dimension(pgbinary_array_length) :: &
      Star_History_track_file_dir, Star_History_track_file_prefix, Star_History_track_title, &
      Star_History_track_xname, Star_History_track_xaxis_label, &
      Star_History_track_yname, Star_History_track_yaxis_label, Star_History_track_fname
real, dimension(pgbinary_array_length) :: Star_History_track_xmin, Star_History_track_xmax, &
      Star_History_track_ymin, Star_History_track_ymax, &
      Star_History_track_xmargin, Star_History_track_ymargin, &
      Star_History_track_dxmin, Star_History_track_dymin, &
      Star_History_track_win_width, Star_History_track_win_aspect_ratio, &
      Star_History_track_xleft, Star_History_track_xright, &
      Star_History_track_ybot, Star_History_track_ytop, Star_History_track_txt_scale, &
      prev_Star_History_track_win_width, prev_Star_History_track_win_ratio, &
      Star_History_track_file_width, Star_History_track_file_aspect_ratio, &
      prev_Star_History_track_file_width, prev_Star_History_track_file_ratio

logical :: Summary_History_win_flag, Summary_History_file_flag
logical :: do_Summary_History_win, do_Summary_History_file, &
      Summary_History_scaled_value(max_num_Summary_History_Lines)
integer :: id_Summary_History_win, id_Summary_History_file, &
      Summary_History_num_lines, Summary_History_file_interval, &
      Summary_History_xmin, Summary_History_xmax, Summary_History_max_width
character (len = strlen) :: Summary_History_file_dir, Summary_History_file_prefix, &
      Summary_History_title, &
      Summary_History_name(max_num_Summary_History_Lines), &
      Summary_History_legend(max_num_Summary_History_Lines)
real :: Summary_History_win_width, Summary_History_win_aspect_ratio, &
      Summary_History_xleft, Summary_History_xright, &
      Summary_History_ybot, Summary_History_ytop, Summary_History_txt_scale, &
      prev_Summary_History_win_width, prev_Summary_History_win_ratio, &
      Summary_History_file_width, Summary_History_file_aspect_ratio, &
      prev_Summary_History_file_width, prev_Summary_History_file_ratio

logical, dimension(pgbinary_array_length) :: Grid_win_flag, Grid_file_flag
logical, dimension(pgbinary_array_length) :: do_Grid_win, do_Grid_file
integer, dimension(pgbinary_array_length) :: Grid_file_interval, Grid_step_max
integer, dimension(pgbinary_array_length) :: id_Grid_win, id_Grid_file
character (len = strlen), dimension(pgbinary_array_length) :: Grid_file_dir, Grid_file_prefix
real, dimension(pgbinary_array_length) :: Grid_win_width, Grid_win_aspect_ratio, &
      Grid_xleft, Grid_xright, &
      Grid_ybot, Grid_ytop, &
      prev_Grid_win_width, prev_Grid_win_ratio, &
      Grid_file_width, Grid_file_aspect_ratio, &
      prev_Grid_file_width, prev_Grid_file_ratio
character (len = strlen), dimension(pgbinary_array_length) :: Grid_title
character (len = strlen), dimension(pgbinary_array_length, max_num_pgbinary_grid_plots) :: Grid_plot_name
real, dimension(pgbinary_array_length, max_num_pgbinary_grid_plots) :: &
      Grid_txt_scale_factor, &
      Grid_plot_pad_left, Grid_plot_pad_right, &
      Grid_plot_pad_top, Grid_plot_pad_bot
integer, dimension(pgbinary_array_length):: Grid_num_cols, Grid_num_rows, Grid_num_plots
integer, dimension(pgbinary_array_length, max_num_pgbinary_grid_plots) :: &
      Grid_plot_row, Grid_plot_rowspan, &
      Grid_plot_col, Grid_plot_colspan

logical :: Star1_win_flag, Star1_file_flag
logical :: do_Star1_win, do_Star1_file, do_star1_box
integer :: Star1_file_interval, Star1_step_max
integer :: id_Star1_win, id_Star1_file
character (len = strlen) :: Star1_file_dir, Star1_file_prefix
real :: Star1_win_width, Star1_win_aspect_ratio, &
      Star1_xleft, Star1_xright, &
      Star1_ybot, Star1_ytop, &
      prev_Star1_win_width, prev_Star1_win_ratio, &
      Star1_file_width, Star1_file_aspect_ratio, &
      prev_Star1_file_width, prev_Star1_file_ratio, Star1_txt_scale_factor, &
      star1_box_pad_left, star1_box_pad_right, &
      star1_box_pad_bot, star1_box_pad_top
character (len = strlen) :: Star1_title, Star1_plot_name

logical :: Star2_win_flag, Star2_file_flag
logical :: do_Star2_win, do_Star2_file, do_star2_box
integer :: Star2_file_interval, Star2_step_max
integer :: id_Star2_win, id_Star2_file
character (len = strlen) :: Star2_file_dir, Star2_file_prefix
real :: Star2_win_width, Star2_win_aspect_ratio, &
      Star2_xleft, Star2_xright, &
      Star2_ybot, Star2_ytop, &
      prev_Star2_win_width, prev_Star2_win_ratio, &
      Star2_file_width, Star2_file_aspect_ratio, &
      prev_Star2_file_width, prev_Star2_file_ratio, Star2_txt_scale_factor, &
      star2_box_pad_left, star2_box_pad_right, &
      star2_box_pad_bot, star2_box_pad_top
character (len = strlen) :: Star2_title, Star2_plot_name

logical :: show_mtrans_status

logical :: Orbit_win_flag, Orbit_file_flag, do_Orbit_win, do_Orbit_file
logical :: Orbit_show_RL, Orbit_show_stars
integer :: Orbit_file_interval, id_Orbit_win, id_Orbit_file
character (len = strlen) :: Orbit_file_dir, Orbit_file_prefix, Orbit_title
real :: Orbit_win_width, Orbit_win_aspect_ratio, &
      Orbit_xleft, Orbit_xright, &
      Orbit_ybot, Orbit_ytop, &
      prev_Orbit_win_width, prev_Orbit_win_ratio, &
      Orbit_file_width, Orbit_file_aspect_ratio, &
      prev_Orbit_file_width, prev_Orbit_file_ratio, Orbit_txt_scale

integer :: annotation1_ci, annotation1_lw, annotation1_cf
real :: annotation1_ch
character (len = strlen) :: annotation1_text
character :: annotation1_side
real :: annotation1_disp, annotation1_coord, annotation1_fjust

integer :: annotation2_ci, annotation2_lw, annotation2_cf
real :: annotation2_ch
character (len = strlen) :: annotation2_text
character :: annotation2_side
real :: annotation2_disp, annotation2_coord, annotation2_fjust

integer :: annotation3_ci, annotation3_lw, annotation3_cf
real :: annotation3_ch
character (len = 356) :: annotation3_text
character :: annotation3_side
real :: annotation3_disp, annotation3_coord, annotation3_fjust

logical, dimension(max_extra_inlists) :: read_extra_pgbinary_inlist
character (len=strlen), dimension(max_extra_inlists) :: extra_pgbinary_inlist_name

logical, dimension(pgbinary_array_length) :: History_Panels_use_decorator, &
      History_Track_use_decorator, &
      Star_History_track_use_decorator
logical :: summary_history_use_decorator

