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

logical :: Text_Summary1_win_flag, Text_Summary1_file_flag
logical :: do_Text_Summary1_win, do_Text_Summary1_file
integer :: id_Text_Summary1_win, id_Text_Summary1_file, Text_Summary1_file_interval
character (len = strlen) :: Text_Summary1_title, Text_Summary1_file_dir, &
      Text_Summary1_file_prefix
integer :: Text_Summary1_num_cols, Text_Summary1_num_rows
character (len = 64) :: &
      Text_Summary1_name(max_num_rows_Text_Summary, max_num_cols_Text_Summary)
real :: Text_Summary1_xleft, Text_Summary1_xright, &
      Text_Summary1_ybot, Text_Summary1_ytop, Text_Summary1_txt_scale, &
      Text_Summary1_win_width, Text_Summary1_win_aspect_ratio, &
      prev_Text_Summary1_win_width, prev_Text_Summary1_win_ratio, &
      Text_Summary1_file_width, Text_Summary1_file_aspect_ratio, &
      prev_Text_Summary1_file_width, prev_Text_Summary1_file_ratio

logical :: Text_Summary2_win_flag, Text_Summary2_file_flag
logical :: do_Text_Summary2_win, do_Text_Summary2_file
integer :: id_Text_Summary2_win, id_Text_Summary2_file, Text_Summary2_file_interval
character (len = strlen) :: Text_Summary2_title, Text_Summary2_file_dir, &
      Text_Summary2_file_prefix
integer :: Text_Summary2_num_cols, Text_Summary2_num_rows
character (len = 64) :: &
      Text_Summary2_name(max_num_rows_Text_Summary, max_num_cols_Text_Summary)
real :: Text_Summary2_xleft, Text_Summary2_xright, &
      Text_Summary2_ybot, Text_Summary2_ytop, Text_Summary2_txt_scale, &
      Text_Summary2_win_width, Text_Summary2_win_aspect_ratio, &
      prev_Text_Summary2_win_width, prev_Text_Summary2_win_ratio, &
      Text_Summary2_file_width, Text_Summary2_file_aspect_ratio, &
      prev_Text_Summary2_file_width, prev_Text_Summary2_file_ratio

logical :: Text_Summary3_win_flag, Text_Summary3_file_flag
logical :: do_Text_Summary3_win, do_Text_Summary3_file
integer :: id_Text_Summary3_win, id_Text_Summary3_file, Text_Summary3_file_interval
character (len = strlen) :: Text_Summary3_title, Text_Summary3_file_dir, &
      Text_Summary3_file_prefix
integer :: Text_Summary3_num_cols, Text_Summary3_num_rows
character (len = 64) :: &
      Text_Summary3_name(max_num_rows_Text_Summary, max_num_cols_Text_Summary)
real :: Text_Summary3_xleft, Text_Summary3_xright, &
      Text_Summary3_ybot, Text_Summary3_ytop, Text_Summary3_txt_scale, &
      Text_Summary3_win_width, Text_Summary3_win_aspect_ratio, &
      prev_Text_Summary3_win_width, prev_Text_Summary3_win_ratio, &
      Text_Summary3_file_width, Text_Summary3_file_aspect_ratio, &
      prev_Text_Summary3_file_width, prev_Text_Summary3_file_ratio

logical :: Text_Summary4_win_flag, Text_Summary4_file_flag
logical :: do_Text_Summary4_win, do_Text_Summary4_file
integer :: id_Text_Summary4_win, id_Text_Summary4_file, Text_Summary4_file_interval
character (len = strlen) :: Text_Summary4_title, Text_Summary4_file_dir, &
      Text_Summary4_file_prefix
integer :: Text_Summary4_num_cols, Text_Summary4_num_rows
character (len = 64) :: &
      Text_Summary4_name(max_num_rows_Text_Summary, max_num_cols_Text_Summary)
real :: Text_Summary4_xleft, Text_Summary4_xright, &
      Text_Summary4_ybot, Text_Summary4_ytop, Text_Summary4_txt_scale, &
      Text_Summary4_win_width, Text_Summary4_win_aspect_ratio, &
      prev_Text_Summary4_win_width, prev_Text_Summary4_win_ratio, &
      Text_Summary4_file_width, Text_Summary4_file_aspect_ratio, &
      prev_Text_Summary4_file_width, prev_Text_Summary4_file_ratio

logical :: Text_Summary5_win_flag, Text_Summary5_file_flag
logical :: do_Text_Summary5_win, do_Text_Summary5_file
integer :: id_Text_Summary5_win, id_Text_Summary5_file, Text_Summary5_file_interval
character (len = strlen) :: Text_Summary5_title, Text_Summary5_file_dir, &
      Text_Summary5_file_prefix
integer :: Text_Summary5_num_cols, Text_Summary5_num_rows
character (len = 64) :: &
      Text_Summary5_name(max_num_rows_Text_Summary, max_num_cols_Text_Summary)
real :: Text_Summary5_xleft, Text_Summary5_xright, &
      Text_Summary5_ybot, Text_Summary5_ytop, Text_Summary5_txt_scale, &
      Text_Summary5_win_width, Text_Summary5_win_aspect_ratio, &
      prev_Text_Summary5_win_width, prev_Text_Summary5_win_ratio, &
      Text_Summary5_file_width, Text_Summary5_file_aspect_ratio, &
      prev_Text_Summary5_file_width, prev_Text_Summary5_file_ratio

logical :: Text_Summary6_win_flag, Text_Summary6_file_flag
logical :: do_Text_Summary6_win, do_Text_Summary6_file
integer :: id_Text_Summary6_win, id_Text_Summary6_file, Text_Summary6_file_interval
character (len = strlen) :: Text_Summary6_title, Text_Summary6_file_dir, &
      Text_Summary6_file_prefix
integer :: Text_Summary6_num_cols, Text_Summary6_num_rows
character (len = 64) :: &
      Text_Summary6_name(max_num_rows_Text_Summary, max_num_cols_Text_Summary)
real :: Text_Summary6_xleft, Text_Summary6_xright, &
      Text_Summary6_ybot, Text_Summary6_ytop, Text_Summary6_txt_scale, &
      Text_Summary6_win_width, Text_Summary6_win_aspect_ratio, &
      prev_Text_Summary6_win_width, prev_Text_Summary6_win_ratio, &
      Text_Summary6_file_width, Text_Summary6_file_aspect_ratio, &
      prev_Text_Summary6_file_width, prev_Text_Summary6_file_ratio

logical :: Text_Summary7_win_flag, Text_Summary7_file_flag
logical :: do_Text_Summary7_win, do_Text_Summary7_file
integer :: id_Text_Summary7_win, id_Text_Summary7_file, Text_Summary7_file_interval
character (len = strlen) :: Text_Summary7_title, Text_Summary7_file_dir, &
      Text_Summary7_file_prefix
integer :: Text_Summary7_num_cols, Text_Summary7_num_rows
character (len = 64) :: &
      Text_Summary7_name(max_num_rows_Text_Summary, max_num_cols_Text_Summary)
real :: Text_Summary7_xleft, Text_Summary7_xright, &
      Text_Summary7_ybot, Text_Summary7_ytop, Text_Summary7_txt_scale, &
      Text_Summary7_win_width, Text_Summary7_win_aspect_ratio, &
      prev_Text_Summary7_win_width, prev_Text_Summary7_win_ratio, &
      Text_Summary7_file_width, Text_Summary7_file_aspect_ratio, &
      prev_Text_Summary7_file_width, prev_Text_Summary7_file_ratio

logical :: Text_Summary8_win_flag, Text_Summary8_file_flag
logical :: do_Text_Summary8_win, do_Text_Summary8_file
integer :: id_Text_Summary8_win, id_Text_Summary8_file, Text_Summary8_file_interval
character (len = strlen) :: Text_Summary8_title, Text_Summary8_file_dir, &
      Text_Summary8_file_prefix
integer :: Text_Summary8_num_cols, Text_Summary8_num_rows
character (len = 64) :: &
      Text_Summary8_name(max_num_rows_Text_Summary, max_num_cols_Text_Summary)
real :: Text_Summary8_xleft, Text_Summary8_xright, &
      Text_Summary8_ybot, Text_Summary8_ytop, Text_Summary8_txt_scale, &
      Text_Summary8_win_width, Text_Summary8_win_aspect_ratio, &
      prev_Text_Summary8_win_width, prev_Text_Summary8_win_ratio, &
      Text_Summary8_file_width, Text_Summary8_file_aspect_ratio, &
      prev_Text_Summary8_file_width, prev_Text_Summary8_file_ratio

logical :: Text_Summary9_win_flag, Text_Summary9_file_flag
logical :: do_Text_Summary9_win, do_Text_Summary9_file
integer :: id_Text_Summary9_win, id_Text_Summary9_file, Text_Summary9_file_interval
character (len = strlen) :: Text_Summary9_title, Text_Summary9_file_dir, &
      Text_Summary9_file_prefix
integer :: Text_Summary9_num_cols, Text_Summary9_num_rows
character (len = 64) :: &
      Text_Summary9_name(max_num_rows_Text_Summary, max_num_cols_Text_Summary)
real :: Text_Summary9_xleft, Text_Summary9_xright, &
      Text_Summary9_ybot, Text_Summary9_ytop, Text_Summary9_txt_scale, &
      Text_Summary9_win_width, Text_Summary9_win_aspect_ratio, &
      prev_Text_Summary9_win_width, prev_Text_Summary9_win_ratio, &
      Text_Summary9_file_width, Text_Summary9_file_aspect_ratio, &
      prev_Text_Summary9_file_width, prev_Text_Summary9_file_ratio

logical :: History_Track1_win_flag, History_Track1_file_flag
logical :: do_History_Track1_win, do_History_Track1_file
logical :: show_History_Track1_annotation1, &
      show_History_Track1_annotation2, show_History_Track1_annotation3
logical :: show_History_Track1_target_box, &
      History_Track1_reverse_xaxis, History_Track1_reverse_yaxis, &
      History_Track1_log_xaxis, History_Track1_log_yaxis
integer :: id_History_Track1_win, id_History_Track1_file, History_Track1_file_interval
integer :: History_Track1_step_min, History_Track1_step_max, History_Track1_n_sigma
character (len = strlen) :: &
      History_Track1_file_dir, History_Track1_file_prefix, History_Track1_title, &
      History_Track1_xname, History_Track1_xaxis_label, &
      History_Track1_yname, History_Track1_yaxis_label, History_Track1_fname
real :: History_Track1_xmin, History_Track1_xmax, &
      History_Track1_ymin, History_Track1_ymax, &
      History_Track1_xmargin, History_Track1_ymargin, &
      History_Track1_dxmin, History_Track1_dymin, &
      History_Track1_xtarget, History_Track1_xsigma, &
      History_Track1_ytarget, History_Track1_ysigma, &
      History_Track1_win_width, History_Track1_win_aspect_ratio, &
      History_Track1_xleft, History_Track1_xright, &
      History_Track1_ybot, History_Track1_ytop, History_Track1_txt_scale, &
      prev_History_Track1_win_width, prev_History_Track1_win_ratio, &
      History_Track1_file_width, History_Track1_file_aspect_ratio, &
      prev_History_Track1_file_width, prev_History_Track1_file_ratio

logical :: History_Track2_win_flag, History_Track2_file_flag
logical :: do_History_Track2_win, do_History_Track2_file
logical :: show_History_Track2_annotation1, &
      show_History_Track2_annotation2, show_History_Track2_annotation3
logical :: show_History_Track2_target_box, &
      History_Track2_reverse_xaxis, History_Track2_reverse_yaxis, &
      History_Track2_log_xaxis, History_Track2_log_yaxis
integer :: id_History_Track2_win, id_History_Track2_file, History_Track2_file_interval
integer :: History_Track2_step_min, History_Track2_step_max, History_Track2_n_sigma
character (len = strlen) :: &
      History_Track2_file_dir, History_Track2_file_prefix, History_Track2_title, &
      History_Track2_xname, History_Track2_xaxis_label, &
      History_Track2_yname, History_Track2_yaxis_label, History_Track2_fname
real :: History_Track2_xmin, History_Track2_xmax, &
      History_Track2_ymin, History_Track2_ymax, &
      History_Track2_xmargin, History_Track2_ymargin, &
      History_Track2_dxmin, History_Track2_dymin, &
      History_Track2_xtarget, History_Track2_xsigma, &
      History_Track2_ytarget, History_Track2_ysigma, &
      History_Track2_win_width, History_Track2_win_aspect_ratio, &
      History_Track2_xleft, History_Track2_xright, &
      History_Track2_ybot, History_Track2_ytop, History_Track2_txt_scale, &
      prev_History_Track2_win_width, prev_History_Track2_win_ratio, &
      History_Track2_file_width, History_Track2_file_aspect_ratio, &
      prev_History_Track2_file_width, prev_History_Track2_file_ratio

logical :: History_Track3_win_flag, History_Track3_file_flag
logical :: do_History_Track3_win, do_History_Track3_file
logical :: show_History_Track3_annotation1, &
      show_History_Track3_annotation2, show_History_Track3_annotation3
logical :: show_History_Track3_target_box, &
      History_Track3_reverse_xaxis, History_Track3_reverse_yaxis, &
      History_Track3_log_xaxis, History_Track3_log_yaxis
integer :: id_History_Track3_win, id_History_Track3_file, History_Track3_file_interval
integer :: History_Track3_step_min, History_Track3_step_max, History_Track3_n_sigma
character (len = strlen) :: &
      History_Track3_file_dir, History_Track3_file_prefix, History_Track3_title, &
      History_Track3_xname, History_Track3_xaxis_label, &
      History_Track3_yname, History_Track3_yaxis_label, History_Track3_fname
real :: History_Track3_xmin, History_Track3_xmax, &
      History_Track3_ymin, History_Track3_ymax, &
      History_Track3_xmargin, History_Track3_ymargin, &
      History_Track3_dxmin, History_Track3_dymin, &
      History_Track3_xtarget, History_Track3_xsigma, &
      History_Track3_ytarget, History_Track3_ysigma, &
      History_Track3_win_width, History_Track3_win_aspect_ratio, &
      History_Track3_xleft, History_Track3_xright, &
      History_Track3_ybot, History_Track3_ytop, History_Track3_txt_scale, &
      prev_History_Track3_win_width, prev_History_Track3_win_ratio, &
      History_Track3_file_width, History_Track3_file_aspect_ratio, &
      prev_History_Track3_file_width, prev_History_Track3_file_ratio

logical :: History_Track4_win_flag, History_Track4_file_flag
logical :: do_History_Track4_win, do_History_Track4_file
logical :: show_History_Track4_annotation1, &
      show_History_Track4_annotation2, show_History_Track4_annotation3
logical :: show_History_Track4_target_box, &
      History_Track4_reverse_xaxis, History_Track4_reverse_yaxis, &
      History_Track4_log_xaxis, History_Track4_log_yaxis
integer :: id_History_Track4_win, id_History_Track4_file, History_Track4_file_interval
integer :: History_Track4_step_min, History_Track4_step_max, History_Track4_n_sigma
character (len = strlen) :: &
      History_Track4_file_dir, History_Track4_file_prefix, History_Track4_title, &
      History_Track4_xname, History_Track4_xaxis_label, &
      History_Track4_yname, History_Track4_yaxis_label, History_Track4_fname
real :: History_Track4_xmin, History_Track4_xmax, &
      History_Track4_ymin, History_Track4_ymax, &
      History_Track4_xmargin, History_Track4_ymargin, &
      History_Track4_dxmin, History_Track4_dymin, &
      History_Track4_xtarget, History_Track4_xsigma, &
      History_Track4_ytarget, History_Track4_ysigma, &
      History_Track4_win_width, History_Track4_win_aspect_ratio, &
      History_Track4_xleft, History_Track4_xright, &
      History_Track4_ybot, History_Track4_ytop, History_Track4_txt_scale, &
      prev_History_Track4_win_width, prev_History_Track4_win_ratio, &
      History_Track4_file_width, History_Track4_file_aspect_ratio, &
      prev_History_Track4_file_width, prev_History_Track4_file_ratio

logical :: History_Track5_win_flag, History_Track5_file_flag
logical :: do_History_Track5_win, do_History_Track5_file
logical :: show_History_Track5_annotation1, &
      show_History_Track5_annotation2, show_History_Track5_annotation3
logical :: show_History_Track5_target_box, &
      History_Track5_reverse_xaxis, History_Track5_reverse_yaxis, &
      History_Track5_log_xaxis, History_Track5_log_yaxis
integer :: id_History_Track5_win, id_History_Track5_file, History_Track5_file_interval
integer :: History_Track5_step_min, History_Track5_step_max, History_Track5_n_sigma
character (len = strlen) :: &
      History_Track5_file_dir, History_Track5_file_prefix, History_Track5_title, &
      History_Track5_xname, History_Track5_xaxis_label, &
      History_Track5_yname, History_Track5_yaxis_label, History_Track5_fname
real :: History_Track5_xmin, History_Track5_xmax, &
      History_Track5_ymin, History_Track5_ymax, &
      History_Track5_xmargin, History_Track5_ymargin, &
      History_Track5_dxmin, History_Track5_dymin, &
      History_Track5_xtarget, History_Track5_xsigma, &
      History_Track5_ytarget, History_Track5_ysigma, &
      History_Track5_win_width, History_Track5_win_aspect_ratio, &
      History_Track5_xleft, History_Track5_xright, &
      History_Track5_ybot, History_Track5_ytop, History_Track5_txt_scale, &
      prev_History_Track5_win_width, prev_History_Track5_win_ratio, &
      History_Track5_file_width, History_Track5_file_aspect_ratio, &
      prev_History_Track5_file_width, prev_History_Track5_file_ratio

logical :: History_Track6_win_flag, History_Track6_file_flag
logical :: do_History_Track6_win, do_History_Track6_file
logical :: show_History_Track6_annotation1, &
      show_History_Track6_annotation2, show_History_Track6_annotation3
logical :: show_History_Track6_target_box, &
      History_Track6_reverse_xaxis, History_Track6_reverse_yaxis, &
      History_Track6_log_xaxis, History_Track6_log_yaxis
integer :: id_History_Track6_win, id_History_Track6_file, History_Track6_file_interval
integer :: History_Track6_step_min, History_Track6_step_max, History_Track6_n_sigma
character (len = strlen) :: &
      History_Track6_file_dir, History_Track6_file_prefix, History_Track6_title, &
      History_Track6_xname, History_Track6_xaxis_label, &
      History_Track6_yname, History_Track6_yaxis_label, History_Track6_fname
real :: History_Track6_xmin, History_Track6_xmax, &
      History_Track6_ymin, History_Track6_ymax, &
      History_Track6_xmargin, History_Track6_ymargin, &
      History_Track6_dxmin, History_Track6_dymin, &
      History_Track6_xtarget, History_Track6_xsigma, &
      History_Track6_ytarget, History_Track6_ysigma, &
      History_Track6_win_width, History_Track6_win_aspect_ratio, &
      History_Track6_xleft, History_Track6_xright, &
      History_Track6_ybot, History_Track6_ytop, History_Track6_txt_scale, &
      prev_History_Track6_win_width, prev_History_Track6_win_ratio, &
      History_Track6_file_width, History_Track6_file_aspect_ratio, &
      prev_History_Track6_file_width, prev_History_Track6_file_ratio

logical :: History_Track7_win_flag, History_Track7_file_flag
logical :: do_History_Track7_win, do_History_Track7_file
logical :: show_History_Track7_annotation1, &
      show_History_Track7_annotation2, show_History_Track7_annotation3
logical :: show_History_Track7_target_box, &
      History_Track7_reverse_xaxis, History_Track7_reverse_yaxis, &
      History_Track7_log_xaxis, History_Track7_log_yaxis
integer :: id_History_Track7_win, id_History_Track7_file, History_Track7_file_interval
integer :: History_Track7_step_min, History_Track7_step_max, History_Track7_n_sigma
character (len = strlen) :: &
      History_Track7_file_dir, History_Track7_file_prefix, History_Track7_title, &
      History_Track7_xname, History_Track7_xaxis_label, &
      History_Track7_yname, History_Track7_yaxis_label, History_Track7_fname
real :: History_Track7_xmin, History_Track7_xmax, &
      History_Track7_ymin, History_Track7_ymax, &
      History_Track7_xmargin, History_Track7_ymargin, &
      History_Track7_dxmin, History_Track7_dymin, &
      History_Track7_xtarget, History_Track7_xsigma, &
      History_Track7_ytarget, History_Track7_ysigma, &
      History_Track7_win_width, History_Track7_win_aspect_ratio, &
      History_Track7_xleft, History_Track7_xright, &
      History_Track7_ybot, History_Track7_ytop, History_Track7_txt_scale, &
      prev_History_Track7_win_width, prev_History_Track7_win_ratio, &
      History_Track7_file_width, History_Track7_file_aspect_ratio, &
      prev_History_Track7_file_width, prev_History_Track7_file_ratio

logical :: History_Track8_win_flag, History_Track8_file_flag
logical :: do_History_Track8_win, do_History_Track8_file
logical :: show_History_Track8_annotation1, &
      show_History_Track8_annotation2, show_History_Track8_annotation3
logical :: show_History_Track8_target_box, &
      History_Track8_reverse_xaxis, History_Track8_reverse_yaxis, &
      History_Track8_log_xaxis, History_Track8_log_yaxis
integer :: id_History_Track8_win, id_History_Track8_file, History_Track8_file_interval
integer :: History_Track8_step_min, History_Track8_step_max, History_Track8_n_sigma
character (len = strlen) :: &
      History_Track8_file_dir, History_Track8_file_prefix, History_Track8_title, &
      History_Track8_xname, History_Track8_xaxis_label, &
      History_Track8_yname, History_Track8_yaxis_label, History_Track8_fname
real :: History_Track8_xmin, History_Track8_xmax, &
      History_Track8_ymin, History_Track8_ymax, &
      History_Track8_xmargin, History_Track8_ymargin, &
      History_Track8_dxmin, History_Track8_dymin, &
      History_Track8_xtarget, History_Track8_xsigma, &
      History_Track8_ytarget, History_Track8_ysigma, &
      History_Track8_win_width, History_Track8_win_aspect_ratio, &
      History_Track8_xleft, History_Track8_xright, &
      History_Track8_ybot, History_Track8_ytop, History_Track8_txt_scale, &
      prev_History_Track8_win_width, prev_History_Track8_win_ratio, &
      History_Track8_file_width, History_Track8_file_aspect_ratio, &
      prev_History_Track8_file_width, prev_History_Track8_file_ratio

logical :: History_Track9_win_flag, History_Track9_file_flag
logical :: do_History_Track9_win, do_History_Track9_file
logical :: show_History_Track9_annotation1, &
      show_History_Track9_annotation2, show_History_Track9_annotation3
logical :: show_History_Track9_target_box, &
      History_Track9_reverse_xaxis, History_Track9_reverse_yaxis, &
      History_Track9_log_xaxis, History_Track9_log_yaxis
integer :: id_History_Track9_win, id_History_Track9_file, History_Track9_file_interval
integer :: History_Track9_step_min, History_Track9_step_max, History_Track9_n_sigma
character (len = strlen) :: &
      History_Track9_file_dir, History_Track9_file_prefix, History_Track9_title, &
      History_Track9_xname, History_Track9_xaxis_label, &
      History_Track9_yname, History_Track9_yaxis_label, History_Track9_fname
real :: History_Track9_xmin, History_Track9_xmax, &
      History_Track9_ymin, History_Track9_ymax, &
      History_Track9_xmargin, History_Track9_ymargin, &
      History_Track9_dxmin, History_Track9_dymin, &
      History_Track9_xtarget, History_Track9_xsigma, &
      History_Track9_ytarget, History_Track9_ysigma, &
      History_Track9_win_width, History_Track9_win_aspect_ratio, &
      History_Track9_xleft, History_Track9_xright, &
      History_Track9_ybot, History_Track9_ytop, History_Track9_txt_scale, &
      prev_History_Track9_win_width, prev_History_Track9_win_ratio, &
      History_Track9_file_width, History_Track9_file_aspect_ratio, &
      prev_History_Track9_file_width, prev_History_Track9_file_ratio

logical :: History_Panels1_win_flag, History_Panels1_file_flag
logical :: do_History_Panels1_win, do_History_Panels1_file
integer :: id_History_Panels1_win, id_History_Panels1_file, History_Panels1_file_interval
character (len = strlen) :: History_Panels1_file_dir, History_Panels1_file_prefix
logical :: History_Panels1_xaxis_reversed, History_Panels1_xaxis_log
character (len = strlen) :: History_Panels1_xaxis_name, History_Panels1_title
real :: History_Panels1_xmin, History_Panels1_xmax, History_Panels1_dxmin, &
      History_Panels1_max_width, History_Panels1_xmargin, &
      History_Panels1_win_width, History_Panels1_win_aspect_ratio, &
      History_Panels1_xleft, History_Panels1_xright, &
      History_Panels1_ybot, History_Panels1_ytop, History_Panels1_txt_scale, &
      prev_History_Panels1_win_width, prev_History_Panels1_win_ratio, &
      History_Panels1_file_width, History_Panels1_file_aspect_ratio, &
      prev_History_Panels1_file_width, prev_History_Panels1_file_ratio
integer :: History_Panels1_num_panels
character (len = strlen), dimension(max_num_History_Panels) :: &
      History_Panels1_yaxis_name, History_Panels1_other_yaxis_name, &
      History_Panels1_points_name
logical, dimension(max_num_History_Panels) :: &
      History_Panels1_yaxis_reversed, History_Panels1_other_yaxis_reversed, &
      History_Panels1_yaxis_log, History_Panels1_other_yaxis_log
real, dimension(max_num_History_Panels) :: &
      History_Panels1_ymin, History_Panels1_other_ymin, &
      History_Panels1_ymax, History_Panels1_other_ymax, &
      History_Panels1_ymargin, History_Panels1_other_ymargin, &
      History_Panels1_dymin, History_Panels1_other_dymin

logical :: History_Panels2_win_flag, History_Panels2_file_flag
logical :: do_History_Panels2_win, do_History_Panels2_file
integer :: id_History_Panels2_win, id_History_Panels2_file, History_Panels2_file_interval
character (len = strlen) :: History_Panels2_file_dir, History_Panels2_file_prefix
logical :: History_Panels2_xaxis_reversed, History_Panels2_xaxis_log
character (len = strlen) :: History_Panels2_xaxis_name, History_Panels2_title
real :: History_Panels2_xmin, History_Panels2_xmax, History_Panels2_dxmin, &
      History_Panels2_max_width, History_Panels2_xmargin, &
      History_Panels2_win_width, History_Panels2_win_aspect_ratio, &
      History_Panels2_xleft, History_Panels2_xright, &
      History_Panels2_ybot, History_Panels2_ytop, History_Panels2_txt_scale, &
      prev_History_Panels2_win_width, prev_History_Panels2_win_ratio, &
      History_Panels2_file_width, History_Panels2_file_aspect_ratio, &
      prev_History_Panels2_file_width, prev_History_Panels2_file_ratio
integer :: History_Panels2_num_panels
character (len = strlen), dimension(max_num_History_Panels) :: &
      History_Panels2_yaxis_name, History_Panels2_other_yaxis_name, &
      History_Panels2_points_name
logical, dimension(max_num_History_Panels) :: &
      History_Panels2_yaxis_reversed, History_Panels2_other_yaxis_reversed, &
      History_Panels2_yaxis_log, History_Panels2_other_yaxis_log
real, dimension(max_num_History_Panels) :: &
      History_Panels2_ymin, History_Panels2_other_ymin, &
      History_Panels2_ymax, History_Panels2_other_ymax, &
      History_Panels2_ymargin, History_Panels2_other_ymargin, &
      History_Panels2_dymin, History_Panels2_other_dymin

logical :: History_Panels3_win_flag, History_Panels3_file_flag
logical :: do_History_Panels3_win, do_History_Panels3_file
integer :: id_History_Panels3_win, id_History_Panels3_file, History_Panels3_file_interval
character (len = strlen) :: History_Panels3_file_dir, History_Panels3_file_prefix
logical :: History_Panels3_xaxis_reversed, History_Panels3_xaxis_log
character (len = strlen) :: History_Panels3_xaxis_name, History_Panels3_title
real :: History_Panels3_xmin, History_Panels3_xmax, History_Panels3_dxmin, &
      History_Panels3_max_width, History_Panels3_xmargin, &
      History_Panels3_win_width, History_Panels3_win_aspect_ratio, &
      History_Panels3_xleft, History_Panels3_xright, &
      History_Panels3_ybot, History_Panels3_ytop, History_Panels3_txt_scale, &
      prev_History_Panels3_win_width, prev_History_Panels3_win_ratio, &
      History_Panels3_file_width, History_Panels3_file_aspect_ratio, &
      prev_History_Panels3_file_width, prev_History_Panels3_file_ratio
integer :: History_Panels3_num_panels
character (len = strlen), dimension(max_num_History_Panels) :: &
      History_Panels3_yaxis_name, History_Panels3_other_yaxis_name, &
      History_Panels3_points_name
logical, dimension(max_num_History_Panels) :: &
      History_Panels3_yaxis_reversed, History_Panels3_other_yaxis_reversed, &
      History_Panels3_yaxis_log, History_Panels3_other_yaxis_log
real, dimension(max_num_History_Panels) :: &
      History_Panels3_ymin, History_Panels3_other_ymin, &
      History_Panels3_ymax, History_Panels3_other_ymax, &
      History_Panels3_ymargin, History_Panels3_other_ymargin, &
      History_Panels3_dymin, History_Panels3_other_dymin

logical :: History_Panels4_win_flag, History_Panels4_file_flag
logical :: do_History_Panels4_win, do_History_Panels4_file
integer :: id_History_Panels4_win, id_History_Panels4_file, History_Panels4_file_interval
character (len = strlen) :: History_Panels4_file_dir, History_Panels4_file_prefix
logical :: History_Panels4_xaxis_reversed, History_Panels4_xaxis_log
character (len = strlen) :: History_Panels4_xaxis_name, History_Panels4_title
real :: History_Panels4_xmin, History_Panels4_xmax, History_Panels4_dxmin, &
      History_Panels4_max_width, History_Panels4_xmargin, &
      History_Panels4_win_width, History_Panels4_win_aspect_ratio, &
      History_Panels4_xleft, History_Panels4_xright, &
      History_Panels4_ybot, History_Panels4_ytop, History_Panels4_txt_scale, &
      prev_History_Panels4_win_width, prev_History_Panels4_win_ratio, &
      History_Panels4_file_width, History_Panels4_file_aspect_ratio, &
      prev_History_Panels4_file_width, prev_History_Panels4_file_ratio
integer :: History_Panels4_num_panels
character (len = strlen), dimension(max_num_History_Panels) :: &
      History_Panels4_yaxis_name, History_Panels4_other_yaxis_name, &
      History_Panels4_points_name
logical, dimension(max_num_History_Panels) :: &
      History_Panels4_yaxis_reversed, History_Panels4_other_yaxis_reversed, &
      History_Panels4_yaxis_log, History_Panels4_other_yaxis_log
real, dimension(max_num_History_Panels) :: &
      History_Panels4_ymin, History_Panels4_other_ymin, &
      History_Panels4_ymax, History_Panels4_other_ymax, &
      History_Panels4_ymargin, History_Panels4_other_ymargin, &
      History_Panels4_dymin, History_Panels4_other_dymin

logical :: History_Panels5_win_flag, History_Panels5_file_flag
logical :: do_History_Panels5_win, do_History_Panels5_file
integer :: id_History_Panels5_win, id_History_Panels5_file, History_Panels5_file_interval
character (len = strlen) :: History_Panels5_file_dir, History_Panels5_file_prefix
logical :: History_Panels5_xaxis_reversed, History_Panels5_xaxis_log
character (len = strlen) :: History_Panels5_xaxis_name, History_Panels5_title
real :: History_Panels5_xmin, History_Panels5_xmax, History_Panels5_dxmin, &
      History_Panels5_max_width, History_Panels5_xmargin, &
      History_Panels5_win_width, History_Panels5_win_aspect_ratio, &
      History_Panels5_xleft, History_Panels5_xright, &
      History_Panels5_ybot, History_Panels5_ytop, History_Panels5_txt_scale, &
      prev_History_Panels5_win_width, prev_History_Panels5_win_ratio, &
      History_Panels5_file_width, History_Panels5_file_aspect_ratio, &
      prev_History_Panels5_file_width, prev_History_Panels5_file_ratio
integer :: History_Panels5_num_panels
character (len = strlen), dimension(max_num_History_Panels) :: &
      History_Panels5_yaxis_name, History_Panels5_other_yaxis_name, &
      History_Panels5_points_name
logical, dimension(max_num_History_Panels) :: &
      History_Panels5_yaxis_reversed, History_Panels5_other_yaxis_reversed, &
      History_Panels5_yaxis_log, History_Panels5_other_yaxis_log
real, dimension(max_num_History_Panels) :: &
      History_Panels5_ymin, History_Panels5_other_ymin, &
      History_Panels5_ymax, History_Panels5_other_ymax, &
      History_Panels5_ymargin, History_Panels5_other_ymargin, &
      History_Panels5_dymin, History_Panels5_other_dymin

logical :: History_Panels6_win_flag, History_Panels6_file_flag
logical :: do_History_Panels6_win, do_History_Panels6_file
integer :: id_History_Panels6_win, id_History_Panels6_file, History_Panels6_file_interval
character (len = strlen) :: History_Panels6_file_dir, History_Panels6_file_prefix
logical :: History_Panels6_xaxis_reversed, History_Panels6_xaxis_log
character (len = strlen) :: History_Panels6_xaxis_name, History_Panels6_title
real :: History_Panels6_xmin, History_Panels6_xmax, History_Panels6_dxmin, &
      History_Panels6_max_width, History_Panels6_xmargin, &
      History_Panels6_win_width, History_Panels6_win_aspect_ratio, &
      History_Panels6_xleft, History_Panels6_xright, &
      History_Panels6_ybot, History_Panels6_ytop, History_Panels6_txt_scale, &
      prev_History_Panels6_win_width, prev_History_Panels6_win_ratio, &
      History_Panels6_file_width, History_Panels6_file_aspect_ratio, &
      prev_History_Panels6_file_width, prev_History_Panels6_file_ratio
integer :: History_Panels6_num_panels
character (len = strlen), dimension(max_num_History_Panels) :: &
      History_Panels6_yaxis_name, History_Panels6_other_yaxis_name, &
      History_Panels6_points_name
logical, dimension(max_num_History_Panels) :: &
      History_Panels6_yaxis_reversed, History_Panels6_other_yaxis_reversed, &
      History_Panels6_yaxis_log, History_Panels6_other_yaxis_log
real, dimension(max_num_History_Panels) :: &
      History_Panels6_ymin, History_Panels6_other_ymin, &
      History_Panels6_ymax, History_Panels6_other_ymax, &
      History_Panels6_ymargin, History_Panels6_other_ymargin, &
      History_Panels6_dymin, History_Panels6_other_dymin

logical :: History_Panels7_win_flag, History_Panels7_file_flag
logical :: do_History_Panels7_win, do_History_Panels7_file
integer :: id_History_Panels7_win, id_History_Panels7_file, History_Panels7_file_interval
character (len = strlen) :: History_Panels7_file_dir, History_Panels7_file_prefix
logical :: History_Panels7_xaxis_reversed, History_Panels7_xaxis_log
character (len = strlen) :: History_Panels7_xaxis_name, History_Panels7_title
real :: History_Panels7_xmin, History_Panels7_xmax, History_Panels7_dxmin, &
      History_Panels7_max_width, History_Panels7_xmargin, &
      History_Panels7_win_width, History_Panels7_win_aspect_ratio, &
      History_Panels7_xleft, History_Panels7_xright, &
      History_Panels7_ybot, History_Panels7_ytop, History_Panels7_txt_scale, &
      prev_History_Panels7_win_width, prev_History_Panels7_win_ratio, &
      History_Panels7_file_width, History_Panels7_file_aspect_ratio, &
      prev_History_Panels7_file_width, prev_History_Panels7_file_ratio
integer :: History_Panels7_num_panels
character (len = strlen), dimension(max_num_History_Panels) :: &
      History_Panels7_yaxis_name, History_Panels7_other_yaxis_name, &
      History_Panels7_points_name
logical, dimension(max_num_History_Panels) :: &
      History_Panels7_yaxis_reversed, History_Panels7_other_yaxis_reversed, &
      History_Panels7_yaxis_log, History_Panels7_other_yaxis_log
real, dimension(max_num_History_Panels) :: &
      History_Panels7_ymin, History_Panels7_other_ymin, &
      History_Panels7_ymax, History_Panels7_other_ymax, &
      History_Panels7_ymargin, History_Panels7_other_ymargin, &
      History_Panels7_dymin, History_Panels7_other_dymin

logical :: History_Panels8_win_flag, History_Panels8_file_flag
logical :: do_History_Panels8_win, do_History_Panels8_file
integer :: id_History_Panels8_win, id_History_Panels8_file, History_Panels8_file_interval
character (len = strlen) :: History_Panels8_file_dir, History_Panels8_file_prefix
logical :: History_Panels8_xaxis_reversed, History_Panels8_xaxis_log
character (len = strlen) :: History_Panels8_xaxis_name, History_Panels8_title
real :: History_Panels8_xmin, History_Panels8_xmax, History_Panels8_dxmin, &
      History_Panels8_max_width, History_Panels8_xmargin, &
      History_Panels8_win_width, History_Panels8_win_aspect_ratio, &
      History_Panels8_xleft, History_Panels8_xright, &
      History_Panels8_ybot, History_Panels8_ytop, History_Panels8_txt_scale, &
      prev_History_Panels8_win_width, prev_History_Panels8_win_ratio, &
      History_Panels8_file_width, History_Panels8_file_aspect_ratio, &
      prev_History_Panels8_file_width, prev_History_Panels8_file_ratio
integer :: History_Panels8_num_panels
character (len = strlen), dimension(max_num_History_Panels) :: &
      History_Panels8_yaxis_name, History_Panels8_other_yaxis_name, &
      History_Panels8_points_name
logical, dimension(max_num_History_Panels) :: &
      History_Panels8_yaxis_reversed, History_Panels8_other_yaxis_reversed, &
      History_Panels8_yaxis_log, History_Panels8_other_yaxis_log
real, dimension(max_num_History_Panels) :: &
      History_Panels8_ymin, History_Panels8_other_ymin, &
      History_Panels8_ymax, History_Panels8_other_ymax, &
      History_Panels8_ymargin, History_Panels8_other_ymargin, &
      History_Panels8_dymin, History_Panels8_other_dymin

logical :: History_Panels9_win_flag, History_Panels9_file_flag
logical :: do_History_Panels9_win, do_History_Panels9_file
integer :: id_History_Panels9_win, id_History_Panels9_file, History_Panels9_file_interval
character (len = strlen) :: History_Panels9_file_dir, History_Panels9_file_prefix
logical :: History_Panels9_xaxis_reversed, History_Panels9_xaxis_log
character (len = strlen) :: History_Panels9_xaxis_name, History_Panels9_title
real :: History_Panels9_xmin, History_Panels9_xmax, History_Panels9_dxmin, &
      History_Panels9_max_width, History_Panels9_xmargin, &
      History_Panels9_win_width, History_Panels9_win_aspect_ratio, &
      History_Panels9_xleft, History_Panels9_xright, &
      History_Panels9_ybot, History_Panels9_ytop, History_Panels9_txt_scale, &
      prev_History_Panels9_win_width, prev_History_Panels9_win_ratio, &
      History_Panels9_file_width, History_Panels9_file_aspect_ratio, &
      prev_History_Panels9_file_width, prev_History_Panels9_file_ratio
integer :: History_Panels9_num_panels
character (len = strlen), dimension(max_num_History_Panels) :: &
      History_Panels9_yaxis_name, History_Panels9_other_yaxis_name, &
      History_Panels9_points_name
logical, dimension(max_num_History_Panels) :: &
      History_Panels9_yaxis_reversed, History_Panels9_other_yaxis_reversed, &
      History_Panels9_yaxis_log, History_Panels9_other_yaxis_log
real, dimension(max_num_History_Panels) :: &
      History_Panels9_ymin, History_Panels9_other_ymin, &
      History_Panels9_ymax, History_Panels9_other_ymax, &
      History_Panels9_ymargin, History_Panels9_other_ymargin, &
      History_Panels9_dymin, History_Panels9_other_dymin
logical :: History_Panel_points_error_bars
integer :: History_Panel_points_interval, History_Panel_points_marker, &
      History_Panel_points_ci, History_Panel_points_lw
real :: History_Panel_points_ch

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

logical :: Grid1_win_flag, Grid1_file_flag
logical :: do_Grid1_win, do_Grid1_file
integer :: Grid1_file_interval, Grid1_step_max
integer :: id_Grid1_win, id_Grid1_file
character (len = strlen) :: Grid1_file_dir, Grid1_file_prefix
real :: Grid1_win_width, Grid1_win_aspect_ratio, &
      Grid1_xleft, Grid1_xright, &
      Grid1_ybot, Grid1_ytop, &
      prev_Grid1_win_width, prev_Grid1_win_ratio, &
      Grid1_file_width, Grid1_file_aspect_ratio, &
      prev_Grid1_file_width, prev_Grid1_file_ratio
character (len = strlen) :: Grid1_title, &
      Grid1_plot_name(max_num_pgbinary_grid_plots)
real, dimension(max_num_pgbinary_grid_plots) :: &
      Grid1_txt_scale_factor, &
      Grid1_plot_pad_left, Grid1_plot_pad_right, &
      Grid1_plot_pad_top, Grid1_plot_pad_bot
Integer :: Grid1_num_cols, Grid1_num_rows, Grid1_num_plots
integer, dimension(max_num_pgbinary_grid_plots) :: &
      Grid1_plot_row, Grid1_plot_rowspan, &
      Grid1_plot_col, Grid1_plot_colspan

logical :: Grid2_win_flag, Grid2_file_flag
logical :: do_Grid2_win, do_Grid2_file
integer :: Grid2_file_interval, Grid2_step_max
integer :: id_Grid2_win, id_Grid2_file
character (len = strlen) :: Grid2_file_dir, Grid2_file_prefix
real :: Grid2_win_width, Grid2_win_aspect_ratio, &
      Grid2_xleft, Grid2_xright, &
      Grid2_ybot, Grid2_ytop, &
      prev_Grid2_win_width, prev_Grid2_win_ratio, &
      Grid2_file_width, Grid2_file_aspect_ratio, &
      prev_Grid2_file_width, prev_Grid2_file_ratio
character (len = strlen) :: Grid2_title, &
      Grid2_plot_name(max_num_pgbinary_grid_plots)
real, dimension(max_num_pgbinary_grid_plots) :: &
      Grid2_txt_scale_factor, &
      Grid2_plot_pad_left, Grid2_plot_pad_right, &
      Grid2_plot_pad_top, Grid2_plot_pad_bot
Integer :: Grid2_num_cols, Grid2_num_rows, Grid2_num_plots
integer, dimension(max_num_pgbinary_grid_plots) :: &
      Grid2_plot_row, Grid2_plot_rowspan, &
      Grid2_plot_col, Grid2_plot_colspan

logical :: Grid3_win_flag, Grid3_file_flag
logical :: do_Grid3_win, do_Grid3_file
integer :: Grid3_file_interval, Grid3_step_max
integer :: id_Grid3_win, id_Grid3_file
character (len = strlen) :: Grid3_file_dir, Grid3_file_prefix
real :: Grid3_win_width, Grid3_win_aspect_ratio, &
      Grid3_xleft, Grid3_xright, &
      Grid3_ybot, Grid3_ytop, &
      prev_Grid3_win_width, prev_Grid3_win_ratio, &
      Grid3_file_width, Grid3_file_aspect_ratio, &
      prev_Grid3_file_width, prev_Grid3_file_ratio
character (len = strlen) :: Grid3_title, &
      Grid3_plot_name(max_num_pgbinary_grid_plots)
real, dimension(max_num_pgbinary_grid_plots) :: &
      Grid3_txt_scale_factor, &
      Grid3_plot_pad_left, Grid3_plot_pad_right, &
      Grid3_plot_pad_top, Grid3_plot_pad_bot
Integer :: Grid3_num_cols, Grid3_num_rows, Grid3_num_plots
integer, dimension(max_num_pgbinary_grid_plots) :: &
      Grid3_plot_row, Grid3_plot_rowspan, &
      Grid3_plot_col, Grid3_plot_colspan

logical :: Grid4_win_flag, Grid4_file_flag
logical :: do_Grid4_win, do_Grid4_file
integer :: Grid4_file_interval, Grid4_step_max
integer :: id_Grid4_win, id_Grid4_file
character (len = strlen) :: Grid4_file_dir, Grid4_file_prefix
real :: Grid4_win_width, Grid4_win_aspect_ratio, &
      Grid4_xleft, Grid4_xright, &
      Grid4_ybot, Grid4_ytop, &
      prev_Grid4_win_width, prev_Grid4_win_ratio, &
      Grid4_file_width, Grid4_file_aspect_ratio, &
      prev_Grid4_file_width, prev_Grid4_file_ratio
character (len = strlen) :: Grid4_title, &
      Grid4_plot_name(max_num_pgbinary_grid_plots)
real, dimension(max_num_pgbinary_grid_plots) :: &
      Grid4_txt_scale_factor, &
      Grid4_plot_pad_left, Grid4_plot_pad_right, &
      Grid4_plot_pad_top, Grid4_plot_pad_bot
Integer :: Grid4_num_cols, Grid4_num_rows, Grid4_num_plots
integer, dimension(max_num_pgbinary_grid_plots) :: &
      Grid4_plot_row, Grid4_plot_rowspan, &
      Grid4_plot_col, Grid4_plot_colspan

logical :: Grid5_win_flag, Grid5_file_flag
logical :: do_Grid5_win, do_Grid5_file
integer :: Grid5_file_interval, Grid5_step_max
integer :: id_Grid5_win, id_Grid5_file
character (len = strlen) :: Grid5_file_dir, Grid5_file_prefix
real :: Grid5_win_width, Grid5_win_aspect_ratio, &
      Grid5_xleft, Grid5_xright, &
      Grid5_ybot, Grid5_ytop, &
      prev_Grid5_win_width, prev_Grid5_win_ratio, &
      Grid5_file_width, Grid5_file_aspect_ratio, &
      prev_Grid5_file_width, prev_Grid5_file_ratio
character (len = strlen) :: Grid5_title, &
      Grid5_plot_name(max_num_pgbinary_grid_plots)
real, dimension(max_num_pgbinary_grid_plots) :: &
      Grid5_txt_scale_factor, &
      Grid5_plot_pad_left, Grid5_plot_pad_right, &
      Grid5_plot_pad_top, Grid5_plot_pad_bot
Integer :: Grid5_num_cols, Grid5_num_rows, Grid5_num_plots
integer, dimension(max_num_pgbinary_grid_plots) :: &
      Grid5_plot_row, Grid5_plot_rowspan, &
      Grid5_plot_col, Grid5_plot_colspan

logical :: Grid6_win_flag, Grid6_file_flag
logical :: do_Grid6_win, do_Grid6_file
integer :: Grid6_file_interval, Grid6_step_max
integer :: id_Grid6_win, id_Grid6_file
character (len = strlen) :: Grid6_file_dir, Grid6_file_prefix
real :: Grid6_win_width, Grid6_win_aspect_ratio, &
      Grid6_xleft, Grid6_xright, &
      Grid6_ybot, Grid6_ytop, &
      prev_Grid6_win_width, prev_Grid6_win_ratio, &
      Grid6_file_width, Grid6_file_aspect_ratio, &
      prev_Grid6_file_width, prev_Grid6_file_ratio
character (len = strlen) :: Grid6_title, &
      Grid6_plot_name(max_num_pgbinary_grid_plots)
real, dimension(max_num_pgbinary_grid_plots) :: &
      Grid6_txt_scale_factor, &
      Grid6_plot_pad_left, Grid6_plot_pad_right, &
      Grid6_plot_pad_top, Grid6_plot_pad_bot
Integer :: Grid6_num_cols, Grid6_num_rows, Grid6_num_plots
integer, dimension(max_num_pgbinary_grid_plots) :: &
      Grid6_plot_row, Grid6_plot_rowspan, &
      Grid6_plot_col, Grid6_plot_colspan

logical :: Grid7_win_flag, Grid7_file_flag
logical :: do_Grid7_win, do_Grid7_file
integer :: Grid7_file_interval, Grid7_step_max
integer :: id_Grid7_win, id_Grid7_file
character (len = strlen) :: Grid7_file_dir, Grid7_file_prefix
real :: Grid7_win_width, Grid7_win_aspect_ratio, &
      Grid7_xleft, Grid7_xright, &
      Grid7_ybot, Grid7_ytop, &
      prev_Grid7_win_width, prev_Grid7_win_ratio, &
      Grid7_file_width, Grid7_file_aspect_ratio, &
      prev_Grid7_file_width, prev_Grid7_file_ratio
character (len = strlen) :: Grid7_title, &
      Grid7_plot_name(max_num_pgbinary_grid_plots)
real, dimension(max_num_pgbinary_grid_plots) :: &
      Grid7_txt_scale_factor, &
      Grid7_plot_pad_left, Grid7_plot_pad_right, &
      Grid7_plot_pad_top, Grid7_plot_pad_bot
Integer :: Grid7_num_cols, Grid7_num_rows, Grid7_num_plots
integer, dimension(max_num_pgbinary_grid_plots) :: &
      Grid7_plot_row, Grid7_plot_rowspan, &
      Grid7_plot_col, Grid7_plot_colspan

logical :: Grid8_win_flag, Grid8_file_flag
logical :: do_Grid8_win, do_Grid8_file
integer :: Grid8_file_interval, Grid8_step_max
integer :: id_Grid8_win, id_Grid8_file
character (len = strlen) :: Grid8_file_dir, Grid8_file_prefix
real :: Grid8_win_width, Grid8_win_aspect_ratio, &
      Grid8_xleft, Grid8_xright, &
      Grid8_ybot, Grid8_ytop, &
      prev_Grid8_win_width, prev_Grid8_win_ratio, &
      Grid8_file_width, Grid8_file_aspect_ratio, &
      prev_Grid8_file_width, prev_Grid8_file_ratio
character (len = strlen) :: Grid8_title, &
      Grid8_plot_name(max_num_pgbinary_grid_plots)
real, dimension(max_num_pgbinary_grid_plots) :: &
      Grid8_txt_scale_factor, &
      Grid8_plot_pad_left, Grid8_plot_pad_right, &
      Grid8_plot_pad_top, Grid8_plot_pad_bot
Integer :: Grid8_num_cols, Grid8_num_rows, Grid8_num_plots
integer, dimension(max_num_pgbinary_grid_plots) :: &
      Grid8_plot_row, Grid8_plot_rowspan, &
      Grid8_plot_col, Grid8_plot_colspan

logical :: Grid9_win_flag, Grid9_file_flag
logical :: do_Grid9_win, do_Grid9_file
integer :: Grid9_file_interval, Grid9_step_max
integer :: id_Grid9_win, id_Grid9_file
character (len = strlen) :: Grid9_file_dir, Grid9_file_prefix
real :: Grid9_win_width, Grid9_win_aspect_ratio, &
      Grid9_xleft, Grid9_xright, &
      Grid9_ybot, Grid9_ytop, &
      prev_Grid9_win_width, prev_Grid9_win_ratio, &
      Grid9_file_width, Grid9_file_aspect_ratio, &
      prev_Grid9_file_width, prev_Grid9_file_ratio
character (len = strlen) :: Grid9_title, &
      Grid9_plot_name(max_num_pgbinary_grid_plots)
real, dimension(max_num_pgbinary_grid_plots) :: &
      Grid9_txt_scale_factor, &
      Grid9_plot_pad_left, Grid9_plot_pad_right, &
      Grid9_plot_pad_top, Grid9_plot_pad_bot
Integer :: Grid9_num_cols, Grid9_num_rows, Grid9_num_plots
integer, dimension(max_num_pgbinary_grid_plots) :: &
      Grid9_plot_row, Grid9_plot_rowspan, &
      Grid9_plot_col, Grid9_plot_colspan

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

logical :: History_Panels1_use_decorator, &
      History_Panels2_use_decorator, &
      History_Panels3_use_decorator, &
      History_Panels4_use_decorator, &
      History_Panels5_use_decorator, &
      History_Panels6_use_decorator, &
      History_Panels7_use_decorator, &
      History_Panels8_use_decorator, &
      History_Panels9_use_decorator, &
      History_Track1_use_decorator, &
      History_Track2_use_decorator, &
      History_Track3_use_decorator, &
      History_Track4_use_decorator, &
      History_Track5_use_decorator, &
      History_Track6_use_decorator, &
      History_Track7_use_decorator, &
      History_Track8_use_decorator, &
      History_Track9_use_decorator, &
      summary_history_use_decorator

