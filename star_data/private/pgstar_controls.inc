!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! PGStar state

type (pgstar_win_file_data) :: pgstar_win_file_ptr(num_pgstar_plots)

type (pgstar_hist_node), pointer :: pgstar_hist => null()


! PGstar controls
integer :: pgstar_interval
character (len = 32) :: file_device
character (len = 32) :: file_extension
integer :: file_digits

real :: delta_HR_limit_for_file_output

logical :: file_white_on_black_flag, win_white_on_black_flag, &
      pgstar_show_age, &
      pgstar_show_age_in_seconds, &
      pgstar_show_age_in_minutes, &
      pgstar_show_age_in_hours, &
      pgstar_show_age_in_days, &
      pgstar_show_age_in_years, &
      pgstar_show_log_age_in_years, &
      pgstar_show_model_number, &
      pgstar_grid_show_title, &
      pgstar_show_title, &
      pgstar_report_writing_files

integer :: pgstar_lw, pgstar_box_lw, &
      pgstar_grid_title_lw, &
      pgstar_title_lw, &
      pgstar_age_lw, &
      pgstar_model_lw, &
      pgstar_xaxis_label_lw, &
      pgstar_left_yaxis_label_lw, &
      pgstar_right_yaxis_label_lw, &
      pgstar_profile_line_style, &
      pgstar_history_line_style

real :: &
      pgstar_grid_title_scale, &
      pgstar_grid_title_disp, &
      pgstar_grid_title_coord, &
      pgstar_grid_title_fjust, &
      pgstar_title_scale, &
      pgstar_title_disp, &
      pgstar_title_coord, &
      pgstar_title_fjust, &
      pgstar_age_scale, &
      pgstar_age_disp, &
      pgstar_age_coord, &
      pgstar_age_fjust, &
      pgstar_model_scale, &
      pgstar_model_disp, &
      pgstar_model_coord, &
      pgstar_model_fjust, &
      pgstar_xaxis_label_disp, &
      pgstar_xaxis_label_scale, &
      pgstar_left_yaxis_label_scale, &
      pgstar_right_yaxis_label_scale, &
      pgstar_left_yaxis_label_disp, &
      pgstar_right_yaxis_label_disp, &
      pgstar_num_scale

logical :: pause
integer :: pause_interval
real :: pgstar_sleep
logical :: clear_history

logical :: Text_Summary1_win_flag, Text_Summary1_file_flag
logical :: do_Text_Summary1_win, do_Text_Summary1_file
integer :: id_Text_Summary1_win, id_Text_Summary1_file, Text_Summary1_file_interval
character (len = strlen) :: Text_Summary1_title, Text_Summary1_file_dir, &
      Text_Summary1_file_prefix
integer :: Text_Summary1_num_cols, Text_Summary1_num_rows
character (len = 64) :: &
      Text_Summary1_name(max_num_rows_Text_Summary, max_num_cols_Text_Summary)
real :: &
      Text_Summary1_xleft, Text_Summary1_xright, &
      Text_Summary1_ybot, Text_Summary1_ytop, Text_Summary1_txt_scale, Text_Summary1_dxval, &
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
real :: &
      Text_Summary2_xleft, Text_Summary2_xright, &
      Text_Summary2_ybot, Text_Summary2_ytop, Text_Summary2_txt_scale, Text_Summary2_dxval, &
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
real :: &
      Text_Summary3_xleft, Text_Summary3_xright, &
      Text_Summary3_ybot, Text_Summary3_ytop, Text_Summary3_txt_scale, Text_Summary3_dxval, &
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
real :: &
      Text_Summary4_xleft, Text_Summary4_xright, &
      Text_Summary4_ybot, Text_Summary4_ytop, Text_Summary4_txt_scale, Text_Summary4_dxval, &
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
real :: &
      Text_Summary5_xleft, Text_Summary5_xright, &
      Text_Summary5_ybot, Text_Summary5_ytop, Text_Summary5_txt_scale, Text_Summary5_dxval, &
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
real :: &
      Text_Summary6_xleft, Text_Summary6_xright, &
      Text_Summary6_ybot, Text_Summary6_ytop, Text_Summary6_txt_scale, Text_Summary6_dxval, &
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
real :: &
      Text_Summary7_xleft, Text_Summary7_xright, &
      Text_Summary7_ybot, Text_Summary7_ytop, Text_Summary7_txt_scale, Text_Summary7_dxval, &
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
real :: &
      Text_Summary8_xleft, Text_Summary8_xright, &
      Text_Summary8_ybot, Text_Summary8_ytop, Text_Summary8_txt_scale, Text_Summary8_dxval, &
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
real :: &
      Text_Summary9_xleft, Text_Summary9_xright, &
      Text_Summary9_ybot, Text_Summary9_ytop, Text_Summary9_txt_scale, Text_Summary9_dxval, &
      Text_Summary9_win_width, Text_Summary9_win_aspect_ratio, &
      prev_Text_Summary9_win_width, prev_Text_Summary9_win_ratio, &
      Text_Summary9_file_width, Text_Summary9_file_aspect_ratio, &
      prev_Text_Summary9_file_width, prev_Text_Summary9_file_ratio

logical :: dPg_dnu_win_flag, dPg_dnu_file_flag
logical :: do_dPg_dnu_win, do_dPg_dnu_file
logical :: show_dPg_dnu_annotation1, &
      show_dPg_dnu_annotation2, &
      show_dPg_dnu_annotation3
logical :: show_dPg_dnu_target_box
integer :: id_dPg_dnu_win, id_dPg_dnu_file, dPg_dnu_file_interval
integer :: dPg_dnu_step_min, dPg_dnu_step_max
character (len = strlen) :: &
      dPg_dnu_file_dir, dPg_dnu_file_prefix, dPg_dnu_title, dPg_dnu_fname
real :: &
      dPg_dnu_xleft, dPg_dnu_xright, &
      dPg_dnu_d_delta_nu_min, dPg_dnu_d_delta_Pg_min, &
      dPg_dnu_ybot, dPg_dnu_ytop, dPg_dnu_txt_scale, &
      dPg_dnu_delta_nu_min, dPg_dnu_delta_nu_max, dPg_dnu_delta_nu_margin, &
      dPg_dnu_delta_Pg_min, dPg_dnu_delta_Pg_max, dPg_dnu_delta_Pg_margin, &
      dPg_dnu_win_width, dPg_dnu_win_aspect_ratio, &
      prev_dPg_dnu_win_width, prev_dPg_dnu_win_ratio, &
      dPg_dnu_file_width, dPg_dnu_file_aspect_ratio, &
      prev_dPg_dnu_file_width, prev_dPg_dnu_file_ratio

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
real :: &
      History_Track1_xmin, History_Track1_xmax, &
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
real :: &
      History_Track2_xmin, History_Track2_xmax, &
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
real :: &
      History_Track3_xmin, History_Track3_xmax, &
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
real :: &
      History_Track4_xmin, History_Track4_xmax, &
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
real :: &
      History_Track5_xmin, History_Track5_xmax, &
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
real :: &
      History_Track6_xmin, History_Track6_xmax, &
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
real :: &
      History_Track7_xmin, History_Track7_xmax, &
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
real :: &
      History_Track8_xmin, History_Track8_xmax, &
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
real :: &
      History_Track9_xmin, History_Track9_xmax, &
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

logical :: HR_win_flag, HR_file_flag
logical :: do_HR_win, do_HR_file
logical :: show_HR_annotation1, show_HR_annotation2, show_HR_annotation3
logical :: show_HR_target_box, show_HR_classical_instability_strip, &
      show_HR_Mira_instability_region, show_HR_WD_instabilities
integer :: id_HR_win, id_HR_file, HR_file_interval
integer :: HR_step_min, HR_step_max, HR_target_n_sigma
character (len = strlen) :: HR_file_dir, HR_file_prefix, HR_title, HR_fname
real :: &
      HR_logT_min, HR_logT_max, HR_logT_margin, &
      HR_logL_min, HR_logL_max, HR_logL_margin, &
      HR_dlogT_min, HR_dlogL_min, &
      HR_target_logL, HR_target_logL_sigma, &
      HR_target_logT, HR_target_logT_sigma, &
      HR_win_width, HR_win_aspect_ratio, &
      HR_xleft, HR_xright, HR_ybot, HR_ytop, HR_txt_scale, &
      prev_HR_win_width, prev_HR_win_ratio, &
      HR_file_width, HR_file_aspect_ratio, &
      prev_HR_file_width, prev_HR_file_ratio

logical :: TRho_win_flag, TRho_file_flag
logical :: do_TRho_win, do_TRho_file
logical :: show_TRho_annotation1, show_TRho_annotation2, show_TRho_annotation3, &
      show_TRho_degeneracy_line
integer :: TRho_file_interval, TRho_step_max
integer :: id_TRho_win, id_TRho_file, TRho_step_min
character (len = strlen) :: TRho_file_dir, TRho_file_prefix, TRho_title, TRho_fname
real :: &
      TRho_logT_min, TRho_logT_max, TRho_logT_margin, &
      TRho_logRho_min, TRho_logRho_max, TRho_logRho_margin, &
      TRho_logRho_dlogRho_min, TRho_logT_dlogT_min, &
      TRho_win_width, TRho_win_aspect_ratio, &
      TRho_xleft, TRho_xright, &
      TRho_ybot, TRho_ytop, TRho_txt_scale, &
      prev_TRho_win_width, prev_TRho_win_ratio, &
      TRho_file_width, TRho_file_aspect_ratio, &
      prev_TRho_file_width, prev_TRho_file_ratio

logical :: TmaxRho_win_flag, TmaxRho_file_flag
logical :: do_TmaxRho_win, do_TmaxRho_file
logical :: show_TmaxRho_annotation1, show_TmaxRho_annotation2, show_TmaxRho_annotation3, &
      show_TmaxRho_degeneracy_line
integer :: TmaxRho_file_interval, TmaxRho_step_max
integer :: id_TmaxRho_win, id_TmaxRho_file, TmaxRho_step_min
character (len = strlen) :: TmaxRho_file_dir, TmaxRho_file_prefix, TmaxRho_title, TmaxRho_fname
real :: &
      TmaxRho_logT_min, TmaxRho_logT_max, TmaxRho_logT_margin, &
      TmaxRho_logRho_min, TmaxRho_logRho_max, TmaxRho_logRho_margin, &
      TmaxRho_logRho_dlogRho_min, TmaxRho_logT_dlogT_min, &
      TmaxRho_win_width, TmaxRho_win_aspect_ratio, &
      TmaxRho_xleft, TmaxRho_xright, &
      TmaxRho_ybot, TmaxRho_ytop, TmaxRho_txt_scale, &
      prev_TmaxRho_win_width, prev_TmaxRho_win_ratio, &
      TmaxRho_file_width, TmaxRho_file_aspect_ratio, &
      prev_TmaxRho_file_width, prev_TmaxRho_file_ratio

logical :: logL_R_win_flag, logL_R_file_flag
logical :: do_logL_R_win, do_logL_R_file
logical :: show_logL_R_target_box, show_logL_photosphere_r
logical :: show_logL_R_annotation1, show_logL_R_annotation2, show_logL_R_annotation3
integer :: logL_R_file_interval, logL_R_step_max, logL_R_target_n_sigma
integer :: id_logL_R_win, id_logL_R_file, logL_R_step_min
character (len = strlen) :: &
      logL_R_file_dir, logL_R_file_prefix, logL_R_title, logL_R_fname
real :: &
      logL_R_logL_min, logL_R_logL_max, logL_R_logL_margin, &
      logL_R_R_min, logL_R_R_max, logL_R_R_margin, &
      logL_R_win_width, logL_R_win_aspect_ratio, &
      logL_R_xleft, logL_R_xright, &
      logL_R_dR_min, logL_R_dlogL_min, &
      logL_R_ybot, logL_R_ytop, logL_R_txt_scale, &
      prev_logL_R_win_width, prev_logL_R_win_ratio, &
      logL_R_file_width, logL_R_file_aspect_ratio, &
      prev_logL_R_file_width, prev_logL_R_file_ratio, &
      logL_R_target_logL, logL_R_target_logL_sigma, &
      logL_R_target_R, logL_R_target_R_sigma

logical :: logL_Teff_win_flag, logL_Teff_file_flag
logical :: do_logL_Teff_win, do_logL_Teff_file
logical :: show_logL_Teff_target_box
logical :: show_logL_Teff_annotation1, show_logL_Teff_annotation2, show_logL_Teff_annotation3
integer :: logL_Teff_file_interval, logL_Teff_step_max, logL_Teff_target_n_sigma
integer :: id_logL_Teff_win, id_logL_Teff_file, logL_Teff_step_min
character (len = strlen) :: &
      logL_Teff_file_dir, logL_Teff_file_prefix, logL_Teff_title, logL_Teff_fname
real :: &
      logL_Teff_logL_min, logL_Teff_logL_max, logL_Teff_logL_margin, &
      logL_Teff_Teff_min, logL_Teff_Teff_max, logL_Teff_Teff_margin, &
      logL_Teff_win_width, logL_Teff_win_aspect_ratio, &
      logL_Teff_xleft, logL_Teff_xright, &
      logL_Teff_dTeff_min, logL_Teff_dlogL_min, &
      logL_Teff_ybot, logL_Teff_ytop, logL_Teff_txt_scale, &
      prev_logL_Teff_win_width, prev_logL_Teff_win_ratio, &
      logL_Teff_file_width, logL_Teff_file_aspect_ratio, &
      prev_logL_Teff_file_width, prev_logL_Teff_file_ratio, &
      logL_Teff_target_logL, logL_Teff_target_logL_sigma, &
      logL_Teff_target_Teff, logL_Teff_target_Teff_sigma

logical :: logL_v_win_flag, logL_v_file_flag
logical :: do_logL_v_win, do_logL_v_file
logical :: show_logL_v_target_box, show_logL_photosphere_v
logical :: show_logL_v_annotation1, show_logL_v_annotation2, show_logL_v_annotation3
integer :: logL_v_file_interval, logL_v_step_max, logL_v_target_n_sigma
integer :: id_logL_v_win, id_logL_v_file, logL_v_step_min
character (len = strlen) :: &
      logL_v_file_dir, logL_v_file_prefix, logL_v_title, logL_v_fname
real :: &
      logL_v_logL_min, logL_v_logL_max, logL_v_logL_margin, &
      logL_v_v_min, logL_v_v_max, logL_v_v_margin, &
      logL_v_win_width, logL_v_win_aspect_ratio, &
      logL_v_xleft, logL_v_xright, &
      logL_v_dv_min, logL_v_dlogL_min, &
      logL_v_ybot, logL_v_ytop, logL_v_txt_scale, &
      prev_logL_v_win_width, prev_logL_v_win_ratio, &
      logL_v_file_width, logL_v_file_aspect_ratio, &
      prev_logL_v_file_width, prev_logL_v_file_ratio, &
      logL_v_target_logL, logL_v_target_logL_sigma, &
      logL_v_target_v, logL_v_target_v_sigma

logical :: L_Teff_win_flag, L_Teff_file_flag
logical :: do_L_Teff_win, do_L_Teff_file
logical :: show_L_Teff_target_box
logical :: show_L_Teff_annotation1, show_L_Teff_annotation2, show_L_Teff_annotation3
integer :: L_Teff_file_interval, L_Teff_step_max, L_Teff_target_n_sigma
integer :: id_L_Teff_win, id_L_Teff_file, L_Teff_step_min
character (len = strlen) :: &
      L_Teff_file_dir, L_Teff_file_prefix, L_Teff_title, L_Teff_fname
real :: &
      L_Teff_L_min, L_Teff_L_max, L_Teff_L_margin, &
      L_Teff_Teff_min, L_Teff_Teff_max, L_Teff_Teff_margin, &
      L_Teff_win_width, L_Teff_win_aspect_ratio, &
      L_Teff_xleft, L_Teff_xright, &
      L_Teff_dTeff_min, L_Teff_dL_min, &
      L_Teff_ybot, L_Teff_ytop, L_Teff_txt_scale, &
      prev_L_Teff_win_width, prev_L_Teff_win_ratio, &
      L_Teff_file_width, L_Teff_file_aspect_ratio, &
      prev_L_Teff_file_width, prev_L_Teff_file_ratio, &
      L_Teff_target_L, L_Teff_target_L_sigma, &
      L_Teff_target_Teff, L_Teff_target_Teff_sigma

logical :: L_v_win_flag, L_v_file_flag
logical :: do_L_v_win, do_L_v_file
logical :: show_L_v_target_box
logical :: show_L_v_annotation1, show_L_v_annotation2, show_L_v_annotation3
integer :: L_v_file_interval, L_v_step_max, L_v_target_n_sigma
integer :: id_L_v_win, id_L_v_file, L_v_step_min
character (len = strlen) :: &
      L_v_file_dir, L_v_file_prefix, L_v_title, L_v_fname
real :: &
      L_v_L_min, L_v_L_max, L_v_L_margin, &
      L_v_v_min, L_v_v_max, L_v_v_margin, &
      L_v_win_width, L_v_win_aspect_ratio, &
      L_v_xleft, L_v_xright, &
      L_v_dv_min, L_v_dL_min, &
      L_v_ybot, L_v_ytop, L_v_txt_scale, &
      prev_L_v_win_width, prev_L_v_win_ratio, &
      L_v_file_width, L_v_file_aspect_ratio, &
      prev_L_v_file_width, prev_L_v_file_ratio, &
      L_v_target_L, L_v_target_L_sigma, &
      L_v_target_v, L_v_target_v_sigma

logical :: L_R_win_flag, L_R_file_flag
logical :: do_L_R_win, do_L_R_file
logical :: show_L_R_target_box
logical :: show_L_R_annotation1, show_L_R_annotation2, show_L_R_annotation3
integer :: L_R_file_interval, L_R_step_max, L_R_target_n_sigma
integer :: id_L_R_win, id_L_R_file, L_R_step_min
character (len = strlen) :: &
      L_R_file_dir, L_R_file_prefix, L_R_title, L_R_fname
real :: &
      L_R_L_min, L_R_L_max, L_R_L_margin, &
      L_R_R_min, L_R_R_max, L_R_R_margin, &
      L_R_win_width, L_R_win_aspect_ratio, &
      L_R_xleft, L_R_xright, &
      L_R_dR_min, L_R_dL_min, &
      L_R_ybot, L_R_ytop, L_R_txt_scale, &
      preR_L_R_win_width, preR_L_R_win_ratio, &
      L_R_file_width, L_R_file_aspect_ratio, &
      prev_L_R_file_width, prev_L_R_file_ratio, &
      L_R_target_L, L_R_target_L_sigma, &
      L_R_target_R, L_R_target_R_sigma

logical :: R_Teff_win_flag, R_Teff_file_flag
logical :: do_R_Teff_win, do_R_Teff_file
logical :: show_R_Teff_target_box
logical :: show_R_Teff_annotation1, show_R_Teff_annotation2, show_R_Teff_annotation3
integer :: R_Teff_file_interval, R_Teff_step_max, R_Teff_target_n_sigma
integer :: id_R_Teff_win, id_R_Teff_file, R_Teff_step_min
character (len = strlen) :: &
      R_Teff_file_dir, R_Teff_file_prefix, R_Teff_title, R_Teff_fname
real :: &
      R_Teff_R_min, R_Teff_R_max, R_Teff_R_margin, &
      R_Teff_Teff_min, R_Teff_Teff_max, R_Teff_Teff_margin, &
      R_Teff_win_width, R_Teff_win_aspect_ratio, &
      R_Teff_xleft, R_Teff_xright, &
      R_Teff_dTeff_min, R_Teff_dR_min, &
      R_Teff_ybot, R_Teff_ytop, R_Teff_txt_scale, &
      prev_R_Teff_win_width, prev_R_Teff_win_ratio, &
      R_Teff_file_width, R_Teff_file_aspect_ratio, &
      prev_R_Teff_file_width, prev_R_Teff_file_ratio, &
      R_Teff_target_R, R_Teff_target_R_sigma, &
      R_Teff_target_Teff, R_Teff_target_Teff_sigma

logical :: R_L_win_flag, R_L_file_flag
logical :: do_R_L_win, do_R_L_file
logical :: show_R_L_target_box
logical :: show_R_L_annotation1, show_R_L_annotation2, show_R_L_annotation3
integer :: R_L_file_interval, R_L_step_max, R_L_target_n_sigma
integer :: id_R_L_win, id_R_L_file, R_L_step_min
character (len = strlen) :: &
      R_L_file_dir, R_L_file_prefix, R_L_title, R_L_fname
real :: &
      R_L_R_min, R_L_R_max, R_L_R_margin, &
      R_L_L_min, R_L_L_max, R_L_L_margin, &
      R_L_win_width, R_L_win_aspect_ratio, &
      R_L_xleft, R_L_xright, &
      R_L_dL_min, R_L_dR_min, &
      R_L_ybot, R_L_ytop, R_L_txt_scale, &
      prev_R_L_win_width, prev_R_L_win_ratio, &
      R_L_file_width, R_L_file_aspect_ratio, &
      prev_R_L_file_width, prev_R_L_file_ratio, &
      R_L_target_R, R_L_target_R_sigma, &
      R_L_target_L, R_L_target_L_sigma

logical :: logg_Teff_win_flag, logg_Teff_file_flag
logical :: do_logg_Teff_win, do_logg_Teff_file
logical :: show_logg_Teff_target_box
logical :: show_logg_Teff_annotation1, show_logg_Teff_annotation2, show_logg_Teff_annotation3
integer :: logg_Teff_file_interval, logg_Teff_step_max, logg_Teff_target_n_sigma
integer :: id_logg_Teff_win, id_logg_Teff_file, logg_Teff_step_min
character (len = strlen) :: &
      logg_Teff_file_dir, logg_Teff_file_prefix, logg_Teff_title, logg_Teff_fname
real :: &
      logg_Teff_logg_min, logg_Teff_logg_max, logg_Teff_logg_margin, &
      logg_Teff_Teff_min, logg_Teff_Teff_max, logg_Teff_Teff_margin, &
      logg_Teff_win_width, logg_Teff_win_aspect_ratio, &
      logg_Teff_xleft, logg_Teff_xright, &
      logg_Teff_dTeff_min, logg_Teff_dlogg_min, &
      logg_Teff_ybot, logg_Teff_ytop, logg_Teff_txt_scale, &
      prev_logg_Teff_win_width, prev_logg_Teff_win_ratio, &
      logg_Teff_file_width, logg_Teff_file_aspect_ratio, &
      prev_logg_Teff_file_width, prev_logg_Teff_file_ratio, &
      logg_Teff_target_logg, logg_Teff_target_logg_sigma, &
      logg_Teff_target_Teff, logg_Teff_target_Teff_sigma

logical :: logg_logT_win_flag, logg_logT_file_flag
logical :: do_logg_logT_win, do_logg_logT_file
logical :: show_logg_logT_target_box
logical :: show_logg_logT_annotation1, show_logg_logT_annotation2, show_logg_logT_annotation3
integer :: logg_logT_file_interval, logg_logT_step_max, logg_logT_target_n_sigma
integer :: id_logg_logT_win, id_logg_logT_file, logg_logT_step_min
character (len = strlen) :: &
      logg_logT_file_dir, logg_logT_file_prefix, logg_logT_title, logg_logT_fname
real :: &
      logg_logT_logg_min, logg_logT_logg_max, logg_logT_logg_margin, &
      logg_logT_logT_min, logg_logT_logT_max, logg_logT_logT_margin, &
      logg_logT_win_width, logg_logT_win_aspect_ratio, &
      logg_logT_xleft, logg_logT_xright, &
      logg_logT_dlogT_min, logg_logT_dlogg_min, &
      logg_logT_ybot, logg_logT_ytop, logg_logT_txt_scale, &
      prev_logg_logT_win_width, prev_logg_logT_win_ratio, &
      logg_logT_file_width, logg_logT_file_aspect_ratio, &
      prev_logg_logT_file_width, prev_logg_logT_file_ratio, &
      logg_logT_target_logg, logg_logT_target_logg_sigma, &
      logg_logT_target_logT, logg_logT_target_logT_sigma

logical :: rti_win_flag, rti_file_flag
logical :: rti_xaxis_log, rti_xaxis_reversed, &
      rti_xaxis_in_Myr, rti_xaxis_in_seconds, rti_xaxis_time_from_present
logical :: do_rti_win, do_rti_file
integer :: id_rti_win, id_rti_file, rti_file_interval
integer :: rti_step_xmin, rti_step_xmax
integer :: rti_max_width, rti_interval, rti_line_weight
character (len = strlen) :: rti_file_dir, rti_file_prefix, rti_title, rti_xaxis_name
logical :: show_rti_annotation1, show_rti_annotation2, show_rti_annotation3
real :: &
      rti_mass_max, rti_mass_min, rti_mass_margin, &
      rti_win_width, rti_win_aspect_ratio, &
      rti_xleft, rti_xright, &
      rti_ybot, rti_ytop, rti_txt_scale, &
      prev_rti_win_width, prev_rti_win_ratio, &
      rti_file_width, rti_file_aspect_ratio, &
      prev_rti_file_width, prev_rti_file_ratio, &
      rti_xmax, rti_xmin, rti_xmargin

logical :: Kipp_win_flag, Kipp_file_flag
logical :: Kipp_xaxis_log, Kipp_xaxis_reversed, &
      Kipp_xaxis_in_Myr, Kipp_xaxis_in_seconds, Kipp_xaxis_time_from_present
logical :: do_Kipp_win, do_Kipp_file
integer :: id_Kipp_win, id_Kipp_file, Kipp_file_interval
integer :: Kipp_step_xmin, Kipp_step_xmax
integer :: Kipp_max_width, Kipp_mix_interval, &
      Kipp_burn_line_weight, Kipp_mix_line_weight, &
      Kipp_luminosities_line_weight, Kipp_masses_line_weight
character (len = strlen) :: Kipp_file_dir, Kipp_file_prefix, Kipp_title, Kipp_xaxis_name
logical :: &
      Kipp_show_burn, Kipp_show_mixing, &
      Kipp_show_luminosities, Kipp_show_mass_boundaries
logical :: show_Kipp_annotation1, show_Kipp_annotation2, show_Kipp_annotation3
real :: &
      Kipp_mass_max, Kipp_mass_min, Kipp_mass_margin, &
      Kipp_lgL_max, Kipp_lgL_min, Kipp_lgL_margin, &
      Kipp_win_width, Kipp_win_aspect_ratio, &
      Kipp_xleft, Kipp_xright, &
      Kipp_ybot, Kipp_ytop, Kipp_txt_scale, &
      prev_Kipp_win_width, prev_Kipp_win_ratio, &
      Kipp_file_width, Kipp_file_aspect_ratio, &
      prev_Kipp_file_width, prev_Kipp_file_ratio, &
      Kipp_xmax, Kipp_xmin, Kipp_xmargin, &
      Kipp_burn_type_cutoff

logical :: TRho_Profile_win_flag, TRho_Profile_file_flag, TRho_switch_to_Column_Depth
logical :: TRho_switch_to_mass
logical :: do_TRho_Profile_win, do_TRho_Profile_file
integer :: id_TRho_Profile_win, id_TRho_Profile_file, TRho_Profile_file_interval
character (len = strlen) :: &
      TRho_Profile_file_dir, TRho_Profile_file_prefix, TRho_Profile_title, TRho_Profile_fname
logical :: show_TRho_Profile_text_info, show_TRho_Profile_legend, show_TRho_Profile_mass_locs
logical :: show_TRho_Profile_burn_labels, show_TRho_Profile_gamma1_4_3rd
logical :: show_TRho_Profile_kap_regions, show_TRho_Profile_eos_regions
logical :: show_TRho_Profile_degeneracy_line, show_TRho_Profile_Pgas_Prad_line
logical :: show_TRho_Profile_burn_lines
logical :: show_TRho_Profile_annotation1, show_TRho_Profile_annotation2, show_TRho_Profile_annotation3
real :: &
      TRho_Profile_text_info_xfac, TRho_Profile_text_info_dxfac, &
      TRho_Profile_text_info_yfac, TRho_Profile_text_info_dyfac, &
      TRho_Profile_xmin, TRho_Profile_xmax, TRho_Profile_ymin, TRho_Profile_ymax, &
      TRho_Profile_legend_txt_scale, TRho_Profile_legend_coord, TRho_Profile_legend_fjust, &
      TRho_Profile_legend_disp1, TRho_Profile_legend_del_disp, &
      TRho_Profile_win_width, TRho_Profile_win_aspect_ratio, &
      prev_TRho_Profile_win_width, prev_TRho_Profile_win_ratio, &
      TRho_Profile_file_width, TRho_Profile_file_aspect_ratio, &
      TRho_Profile_xleft, TRho_Profile_xright, &
      TRho_Profile_ybot, TRho_Profile_ytop, TRho_Profile_txt_scale, &
      prev_TRho_Profile_file_width, prev_TRho_Profile_file_ratio
logical :: show_TRho_accretion_mesh_borders

integer :: num_profile_mass_points
real :: profile_mass_point_q(max_num_profile_mass_points)
integer :: profile_mass_point_color_index(max_num_profile_mass_points)
integer :: profile_mass_point_symbol(max_num_profile_mass_points)
real :: profile_mass_point_symbol_scale(max_num_profile_mass_points)
character (len = 50) :: profile_mass_point_str(max_num_profile_mass_points)
integer :: profile_mass_point_str_clr(max_num_profile_mass_points)
real :: profile_mass_point_str_scale(max_num_profile_mass_points)

logical :: Profile_Panels_show_Mach_1_location, &
      Profile_Panels_show_photosphere_location
real :: Profile_Panels_xwidth_left_div_shock_value, &
      Profile_Panels_xwidth_right_div_shock_value, &
      Profile_Panels_xwidth_left_of_shock, &
      Profile_Panels_xwidth_right_of_shock

logical :: Profile_Panels1_win_flag, Profile_Panels1_file_flag
logical :: do_Profile_Panels1_win, do_Profile_Panels1_file
integer :: id_Profile_Panels1_win, id_Profile_Panels1_file, Profile_Panels1_file_interval
character (len = strlen) :: Profile_Panels1_file_dir, Profile_Panels1_file_prefix
logical :: Profile_Panels1_xaxis_reversed, Profile_Panels1_show_mix_regions_on_xaxis
character (len = strlen) :: Profile_Panels1_xaxis_name, Profile_Panels1_title
real :: &
      Profile_Panels1_xmin, Profile_Panels1_xmax, Profile_Panels1_xmargin, &
      Profile_Panels1_win_width, Profile_Panels1_win_aspect_ratio, &
      Profile_Panels1_xleft, Profile_Panels1_xright, &
      Profile_Panels1_ybot, Profile_Panels1_ytop, Profile_Panels1_txt_scale, &
      prev_Profile_Panels1_win_width, prev_Profile_Panels1_win_ratio, &
      Profile_Panels1_file_width, Profile_Panels1_file_aspect_ratio, &
      prev_Profile_Panels1_file_width, prev_Profile_Panels1_file_ratio
integer :: Profile_Panels1_num_panels
character (len = strlen), dimension(max_num_Profile_Panels) :: &
      Profile_Panels1_yaxis_name, Profile_Panels1_other_yaxis_name
logical, dimension(max_num_Profile_Panels) :: &
      Profile_Panels1_yaxis_reversed, Profile_Panels1_other_yaxis_reversed, &
      Profile_Panels1_yaxis_log, Profile_Panels1_other_yaxis_log, &
      Profile_Panels1_same_yaxis_range
real, dimension(max_num_Profile_Panels) :: &
      Profile_Panels1_ymin, Profile_Panels1_other_ymin, &
      Profile_Panels1_ymax, Profile_Panels1_other_ymax, &
      Profile_Panels1_ycenter, Profile_Panels1_other_ycenter, &
      Profile_Panels1_ymargin, Profile_Panels1_other_ymargin, &
      Profile_Panels1_dymin, Profile_Panels1_other_dymin
logical :: Profile_Panels1_show_grid

logical :: Profile_Panels2_win_flag, Profile_Panels2_file_flag
logical :: do_Profile_Panels2_win, do_Profile_Panels2_file
integer :: id_Profile_Panels2_win, id_Profile_Panels2_file, Profile_Panels2_file_interval
character (len = strlen) :: Profile_Panels2_file_dir, Profile_Panels2_file_prefix
logical :: Profile_Panels2_xaxis_reversed, Profile_Panels2_show_mix_regions_on_xaxis
character (len = strlen) :: Profile_Panels2_xaxis_name, Profile_Panels2_title
real :: &
      Profile_Panels2_xmin, Profile_Panels2_xmax, Profile_Panels2_xmargin, &
      Profile_Panels2_win_width, Profile_Panels2_win_aspect_ratio, &
      Profile_Panels2_xleft, Profile_Panels2_xright, &
      Profile_Panels2_ybot, Profile_Panels2_ytop, Profile_Panels2_txt_scale, &
      prev_Profile_Panels2_win_width, prev_Profile_Panels2_win_ratio, &
      Profile_Panels2_file_width, Profile_Panels2_file_aspect_ratio, &
      prev_Profile_Panels2_file_width, prev_Profile_Panels2_file_ratio
integer :: Profile_Panels2_num_panels
character (len = strlen), dimension(max_num_Profile_Panels) :: &
      Profile_Panels2_yaxis_name, Profile_Panels2_other_yaxis_name
logical, dimension(max_num_Profile_Panels) :: &
      Profile_Panels2_yaxis_reversed, Profile_Panels2_other_yaxis_reversed, &
      Profile_Panels2_yaxis_log, Profile_Panels2_other_yaxis_log, &
      Profile_Panels2_same_yaxis_range
real, dimension(max_num_Profile_Panels) :: &
      Profile_Panels2_ymin, Profile_Panels2_other_ymin, &
      Profile_Panels2_ymax, Profile_Panels2_other_ymax, &
      Profile_Panels2_ycenter, Profile_Panels2_other_ycenter, &
      Profile_Panels2_ymargin, Profile_Panels2_other_ymargin, &
      Profile_Panels2_dymin, Profile_Panels2_other_dymin
logical :: Profile_Panels2_show_grid

logical :: Profile_Panels3_win_flag, Profile_Panels3_file_flag
logical :: do_Profile_Panels3_win, do_Profile_Panels3_file
integer :: id_Profile_Panels3_win, id_Profile_Panels3_file, Profile_Panels3_file_interval
character (len = strlen) :: Profile_Panels3_file_dir, Profile_Panels3_file_prefix
logical :: Profile_Panels3_xaxis_reversed, Profile_Panels3_show_mix_regions_on_xaxis
character (len = strlen) :: Profile_Panels3_xaxis_name, Profile_Panels3_title
real :: &
      Profile_Panels3_xmin, Profile_Panels3_xmax, Profile_Panels3_xmargin, &
      Profile_Panels3_win_width, Profile_Panels3_win_aspect_ratio, &
      Profile_Panels3_xleft, Profile_Panels3_xright, &
      Profile_Panels3_ybot, Profile_Panels3_ytop, Profile_Panels3_txt_scale, &
      prev_Profile_Panels3_win_width, prev_Profile_Panels3_win_ratio, &
      Profile_Panels3_file_width, Profile_Panels3_file_aspect_ratio, &
      prev_Profile_Panels3_file_width, prev_Profile_Panels3_file_ratio
integer :: Profile_Panels3_num_panels
character (len = strlen), dimension(max_num_Profile_Panels) :: &
      Profile_Panels3_yaxis_name, Profile_Panels3_other_yaxis_name
logical, dimension(max_num_Profile_Panels) :: &
      Profile_Panels3_yaxis_reversed, Profile_Panels3_other_yaxis_reversed, &
      Profile_Panels3_yaxis_log, Profile_Panels3_other_yaxis_log, &
      Profile_Panels3_same_yaxis_range
real, dimension(max_num_Profile_Panels) :: &
      Profile_Panels3_ymin, Profile_Panels3_other_ymin, &
      Profile_Panels3_ymax, Profile_Panels3_other_ymax, &
      Profile_Panels3_ycenter, Profile_Panels3_other_ycenter, &
      Profile_Panels3_ymargin, Profile_Panels3_other_ymargin, &
      Profile_Panels3_dymin, Profile_Panels3_other_dymin
logical :: Profile_Panels3_show_grid

logical :: Profile_Panels4_win_flag, Profile_Panels4_file_flag
logical :: do_Profile_Panels4_win, do_Profile_Panels4_file
integer :: id_Profile_Panels4_win, id_Profile_Panels4_file, Profile_Panels4_file_interval
character (len = strlen) :: Profile_Panels4_file_dir, Profile_Panels4_file_prefix
logical :: Profile_Panels4_xaxis_reversed, Profile_Panels4_show_mix_regions_on_xaxis
character (len = strlen) :: Profile_Panels4_xaxis_name, Profile_Panels4_title
real :: &
      Profile_Panels4_xmin, Profile_Panels4_xmax, Profile_Panels4_xmargin, &
      Profile_Panels4_win_width, Profile_Panels4_win_aspect_ratio, &
      Profile_Panels4_xleft, Profile_Panels4_xright, &
      Profile_Panels4_ybot, Profile_Panels4_ytop, Profile_Panels4_txt_scale, &
      prev_Profile_Panels4_win_width, prev_Profile_Panels4_win_ratio, &
      Profile_Panels4_file_width, Profile_Panels4_file_aspect_ratio, &
      prev_Profile_Panels4_file_width, prev_Profile_Panels4_file_ratio
integer :: Profile_Panels4_num_panels
character (len = strlen), dimension(max_num_Profile_Panels) :: &
      Profile_Panels4_yaxis_name, Profile_Panels4_other_yaxis_name
logical, dimension(max_num_Profile_Panels) :: &
      Profile_Panels4_yaxis_reversed, Profile_Panels4_other_yaxis_reversed, &
      Profile_Panels4_yaxis_log, Profile_Panels4_other_yaxis_log, &
      Profile_Panels4_same_yaxis_range
real, dimension(max_num_Profile_Panels) :: &
      Profile_Panels4_ymin, Profile_Panels4_other_ymin, &
      Profile_Panels4_ymax, Profile_Panels4_other_ymax, &
      Profile_Panels4_ycenter, Profile_Panels4_other_ycenter, &
      Profile_Panels4_ymargin, Profile_Panels4_other_ymargin, &
      Profile_Panels4_dymin, Profile_Panels4_other_dymin
logical :: Profile_Panels4_show_grid

logical :: Profile_Panels5_win_flag, Profile_Panels5_file_flag
logical :: do_Profile_Panels5_win, do_Profile_Panels5_file
integer :: id_Profile_Panels5_win, id_Profile_Panels5_file, Profile_Panels5_file_interval
character (len = strlen) :: Profile_Panels5_file_dir, Profile_Panels5_file_prefix
logical :: Profile_Panels5_xaxis_reversed, Profile_Panels5_show_mix_regions_on_xaxis
character (len = strlen) :: Profile_Panels5_xaxis_name, Profile_Panels5_title
real :: &
      Profile_Panels5_xmin, Profile_Panels5_xmax, Profile_Panels5_xmargin, &
      Profile_Panels5_win_width, Profile_Panels5_win_aspect_ratio, &
      Profile_Panels5_xleft, Profile_Panels5_xright, &
      Profile_Panels5_ybot, Profile_Panels5_ytop, Profile_Panels5_txt_scale, &
      prev_Profile_Panels5_win_width, prev_Profile_Panels5_win_ratio, &
      Profile_Panels5_file_width, Profile_Panels5_file_aspect_ratio, &
      prev_Profile_Panels5_file_width, prev_Profile_Panels5_file_ratio
integer :: Profile_Panels5_num_panels
character (len = strlen), dimension(max_num_Profile_Panels) :: &
      Profile_Panels5_yaxis_name, Profile_Panels5_other_yaxis_name
logical, dimension(max_num_Profile_Panels) :: &
      Profile_Panels5_yaxis_reversed, Profile_Panels5_other_yaxis_reversed, &
      Profile_Panels5_yaxis_log, Profile_Panels5_other_yaxis_log, &
      Profile_Panels5_same_yaxis_range
real, dimension(max_num_Profile_Panels) :: &
      Profile_Panels5_ymin, Profile_Panels5_other_ymin, &
      Profile_Panels5_ymax, Profile_Panels5_other_ymax, &
      Profile_Panels5_ycenter, Profile_Panels5_other_ycenter, &
      Profile_Panels5_ymargin, Profile_Panels5_other_ymargin, &
      Profile_Panels5_dymin, Profile_Panels5_other_dymin
logical :: Profile_Panels5_show_grid

logical :: Profile_Panels6_win_flag, Profile_Panels6_file_flag
logical :: do_Profile_Panels6_win, do_Profile_Panels6_file
integer :: id_Profile_Panels6_win, id_Profile_Panels6_file, Profile_Panels6_file_interval
character (len = strlen) :: Profile_Panels6_file_dir, Profile_Panels6_file_prefix
logical :: Profile_Panels6_xaxis_reversed, Profile_Panels6_show_mix_regions_on_xaxis
character (len = strlen) :: Profile_Panels6_xaxis_name, Profile_Panels6_title
real :: &
      Profile_Panels6_xmin, Profile_Panels6_xmax, Profile_Panels6_xmargin, &
      Profile_Panels6_win_width, Profile_Panels6_win_aspect_ratio, &
      Profile_Panels6_xleft, Profile_Panels6_xright, &
      Profile_Panels6_ybot, Profile_Panels6_ytop, Profile_Panels6_txt_scale, &
      prev_Profile_Panels6_win_width, prev_Profile_Panels6_win_ratio, &
      Profile_Panels6_file_width, Profile_Panels6_file_aspect_ratio, &
      prev_Profile_Panels6_file_width, prev_Profile_Panels6_file_ratio
integer :: Profile_Panels6_num_panels
character (len = strlen), dimension(max_num_Profile_Panels) :: &
      Profile_Panels6_yaxis_name, Profile_Panels6_other_yaxis_name
logical, dimension(max_num_Profile_Panels) :: &
      Profile_Panels6_yaxis_reversed, Profile_Panels6_other_yaxis_reversed, &
      Profile_Panels6_yaxis_log, Profile_Panels6_other_yaxis_log, &
      Profile_Panels6_same_yaxis_range
real, dimension(max_num_Profile_Panels) :: &
      Profile_Panels6_ymin, Profile_Panels6_other_ymin, &
      Profile_Panels6_ymax, Profile_Panels6_other_ymax, &
      Profile_Panels6_ycenter, Profile_Panels6_other_ycenter, &
      Profile_Panels6_ymargin, Profile_Panels6_other_ymargin, &
      Profile_Panels6_dymin, Profile_Panels6_other_dymin
logical :: Profile_Panels6_show_grid

logical :: Profile_Panels7_win_flag, Profile_Panels7_file_flag
logical :: do_Profile_Panels7_win, do_Profile_Panels7_file
integer :: id_Profile_Panels7_win, id_Profile_Panels7_file, Profile_Panels7_file_interval
character (len = strlen) :: Profile_Panels7_file_dir, Profile_Panels7_file_prefix
logical :: Profile_Panels7_xaxis_reversed, Profile_Panels7_show_mix_regions_on_xaxis
character (len = strlen) :: Profile_Panels7_xaxis_name, Profile_Panels7_title
real :: &
      Profile_Panels7_xmin, Profile_Panels7_xmax, Profile_Panels7_xmargin, &
      Profile_Panels7_win_width, Profile_Panels7_win_aspect_ratio, &
      Profile_Panels7_xleft, Profile_Panels7_xright, &
      Profile_Panels7_ybot, Profile_Panels7_ytop, Profile_Panels7_txt_scale, &
      prev_Profile_Panels7_win_width, prev_Profile_Panels7_win_ratio, &
      Profile_Panels7_file_width, Profile_Panels7_file_aspect_ratio, &
      prev_Profile_Panels7_file_width, prev_Profile_Panels7_file_ratio
integer :: Profile_Panels7_num_panels
character (len = strlen), dimension(max_num_Profile_Panels) :: &
      Profile_Panels7_yaxis_name, Profile_Panels7_other_yaxis_name
logical, dimension(max_num_Profile_Panels) :: &
      Profile_Panels7_yaxis_reversed, Profile_Panels7_other_yaxis_reversed, &
      Profile_Panels7_yaxis_log, Profile_Panels7_other_yaxis_log, &
      Profile_Panels7_same_yaxis_range
real, dimension(max_num_Profile_Panels) :: &
      Profile_Panels7_ymin, Profile_Panels7_other_ymin, &
      Profile_Panels7_ymax, Profile_Panels7_other_ymax, &
      Profile_Panels7_ycenter, Profile_Panels7_other_ycenter, &
      Profile_Panels7_ymargin, Profile_Panels7_other_ymargin, &
      Profile_Panels7_dymin, Profile_Panels7_other_dymin
logical :: Profile_Panels7_show_grid

logical :: Profile_Panels8_win_flag, Profile_Panels8_file_flag
logical :: do_Profile_Panels8_win, do_Profile_Panels8_file
integer :: id_Profile_Panels8_win, id_Profile_Panels8_file, Profile_Panels8_file_interval
character (len = strlen) :: Profile_Panels8_file_dir, Profile_Panels8_file_prefix
logical :: Profile_Panels8_xaxis_reversed, Profile_Panels8_show_mix_regions_on_xaxis
character (len = strlen) :: Profile_Panels8_xaxis_name, Profile_Panels8_title
real :: &
      Profile_Panels8_xmin, Profile_Panels8_xmax, Profile_Panels8_xmargin, &
      Profile_Panels8_win_width, Profile_Panels8_win_aspect_ratio, &
      Profile_Panels8_xleft, Profile_Panels8_xright, &
      Profile_Panels8_ybot, Profile_Panels8_ytop, Profile_Panels8_txt_scale, &
      prev_Profile_Panels8_win_width, prev_Profile_Panels8_win_ratio, &
      Profile_Panels8_file_width, Profile_Panels8_file_aspect_ratio, &
      prev_Profile_Panels8_file_width, prev_Profile_Panels8_file_ratio
integer :: Profile_Panels8_num_panels
character (len = strlen), dimension(max_num_Profile_Panels) :: &
      Profile_Panels8_yaxis_name, Profile_Panels8_other_yaxis_name
logical, dimension(max_num_Profile_Panels) :: &
      Profile_Panels8_yaxis_reversed, Profile_Panels8_other_yaxis_reversed, &
      Profile_Panels8_yaxis_log, Profile_Panels8_other_yaxis_log, &
      Profile_Panels8_same_yaxis_range
real, dimension(max_num_Profile_Panels) :: &
      Profile_Panels8_ymin, Profile_Panels8_other_ymin, &
      Profile_Panels8_ymax, Profile_Panels8_other_ymax, &
      Profile_Panels8_ycenter, Profile_Panels8_other_ycenter, &
      Profile_Panels8_ymargin, Profile_Panels8_other_ymargin, &
      Profile_Panels8_dymin, Profile_Panels8_other_dymin
logical :: Profile_Panels8_show_grid

logical :: Profile_Panels9_win_flag, Profile_Panels9_file_flag
logical :: do_Profile_Panels9_win, do_Profile_Panels9_file
integer :: id_Profile_Panels9_win, id_Profile_Panels9_file, Profile_Panels9_file_interval
character (len = strlen) :: Profile_Panels9_file_dir, Profile_Panels9_file_prefix
logical :: Profile_Panels9_xaxis_reversed, Profile_Panels9_show_mix_regions_on_xaxis
character (len = strlen) :: Profile_Panels9_xaxis_name, Profile_Panels9_title
real :: &
      Profile_Panels9_xmin, Profile_Panels9_xmax, Profile_Panels9_xmargin, &
      Profile_Panels9_win_width, Profile_Panels9_win_aspect_ratio, &
      Profile_Panels9_xleft, Profile_Panels9_xright, &
      Profile_Panels9_ybot, Profile_Panels9_ytop, Profile_Panels9_txt_scale, &
      prev_Profile_Panels9_win_width, prev_Profile_Panels9_win_ratio, &
      Profile_Panels9_file_width, Profile_Panels9_file_aspect_ratio, &
      prev_Profile_Panels9_file_width, prev_Profile_Panels9_file_ratio
integer :: Profile_Panels9_num_panels
character (len = strlen), dimension(max_num_Profile_Panels) :: &
      Profile_Panels9_yaxis_name, Profile_Panels9_other_yaxis_name
logical, dimension(max_num_Profile_Panels) :: &
      Profile_Panels9_yaxis_reversed, Profile_Panels9_other_yaxis_reversed, &
      Profile_Panels9_yaxis_log, Profile_Panels9_other_yaxis_log, &
      Profile_Panels9_same_yaxis_range
real, dimension(max_num_Profile_Panels) :: &
      Profile_Panels9_ymin, Profile_Panels9_other_ymin, &
      Profile_Panels9_ymax, Profile_Panels9_other_ymax, &
      Profile_Panels9_ycenter, Profile_Panels9_other_ycenter, &
      Profile_Panels9_ymargin, Profile_Panels9_other_ymargin, &
      Profile_Panels9_dymin, Profile_Panels9_other_dymin
logical :: Profile_Panels9_show_grid

logical :: Dynamo_win_flag, Dynamo_file_flag
logical :: do_Dynamo_win, do_Dynamo_file
integer :: id_Dynamo_win, id_Dynamo_file, Dynamo_file_interval
character (len = strlen) :: Dynamo_file_dir, Dynamo_file_prefix
logical :: show_Dynamo_annotation1, show_Dynamo_annotation2, show_Dynamo_annotation3, &
      Dynamo_xaxis_reversed
character (len = strlen) :: Dynamo_xaxis_name, Dynamo_title
real :: &
      Dynamo_xmin, Dynamo_xmax, Dynamo_legend_txt_scale_factor, &
      Dynamo_ymin_left, Dynamo_ymax_left, Dynamo_dymin_left, &
      Dynamo_ymin_right, Dynamo_ymax_right, Dynamo_dymin_right, &
      Dynamo_win_width, Dynamo_win_aspect_ratio, &
      Dynamo_xleft, Dynamo_xright, &
      Dynamo_ybot, Dynamo_ytop, Dynamo_txt_scale, &
      prev_Dynamo_win_width, prev_Dynamo_win_ratio, &
      Dynamo_file_width, Dynamo_file_aspect_ratio, &
      prev_Dynamo_file_width, prev_Dynamo_file_ratio

logical :: Mixing_win_flag, Mixing_file_flag
logical :: do_Mixing_win, do_Mixing_file, &
      Mixing_xaxis_reversed, Mixing_show_rotation_details
integer :: id_Mixing_win, id_Mixing_file, Mixing_file_interval
character (len = strlen) :: Mixing_file_dir, Mixing_file_prefix
logical :: show_Mixing_annotation1, show_Mixing_annotation2, show_Mixing_annotation3
character (len = strlen) :: Mixing_xaxis_name, Mixing_title
real :: &
      Mixing_xmin, Mixing_xmax, Mixing_dymin, &
      Mixing_ymin, Mixing_ymax, Mixing_legend_txt_scale_factor, &
      Mixing_win_width, Mixing_win_aspect_ratio, &
      Mixing_xleft, Mixing_xright, &
      Mixing_ybot, Mixing_ytop, Mixing_txt_scale, &
      prev_Mixing_win_width, prev_Mixing_win_ratio, &
      Mixing_file_width, Mixing_file_aspect_ratio, &
      prev_Mixing_file_width, prev_Mixing_file_ratio

logical :: History_Panels1_win_flag, History_Panels1_file_flag
logical :: do_History_Panels1_win, do_History_Panels1_file
integer :: id_History_Panels1_win, id_History_Panels1_file, History_Panels1_file_interval
character (len = strlen) :: History_Panels1_file_dir, History_Panels1_file_prefix
logical :: History_Panels1_xaxis_reversed, History_Panels1_automatic_star_age_units, History_Panels1_xaxis_log
character (len = strlen) :: History_Panels1_xaxis_name, History_Panels1_title
real :: &
      History_Panels1_xmin, History_Panels1_xmax, History_Panels1_dxmin, &
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
      History_Panels1_yaxis_log, History_Panels1_other_yaxis_log, &
      History_Panels1_same_yaxis_range
real, dimension(max_num_History_Panels) :: &
      History_Panels1_ymin, History_Panels1_other_ymin, &
      History_Panels1_ymax, History_Panels1_other_ymax, &
      History_Panels1_ymargin, History_Panels1_other_ymargin, &
      History_Panels1_dymin, History_Panels1_other_dymin

logical :: History_Panels2_win_flag, History_Panels2_file_flag
logical :: do_History_Panels2_win, do_History_Panels2_file
integer :: id_History_Panels2_win, id_History_Panels2_file, History_Panels2_file_interval
character (len = strlen) :: History_Panels2_file_dir, History_Panels2_file_prefix
logical :: History_Panels2_xaxis_reversed, History_Panels2_automatic_star_age_units, History_Panels2_xaxis_log
character (len = strlen) :: History_Panels2_xaxis_name, History_Panels2_title
real :: &
      History_Panels2_xmin, History_Panels2_xmax, History_Panels2_dxmin, &
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
      History_Panels2_yaxis_log, History_Panels2_other_yaxis_log, &
      History_Panels2_same_yaxis_range
real, dimension(max_num_History_Panels) :: &
      History_Panels2_ymin, History_Panels2_other_ymin, &
      History_Panels2_ymax, History_Panels2_other_ymax, &
      History_Panels2_ymargin, History_Panels2_other_ymargin, &
      History_Panels2_dymin, History_Panels2_other_dymin

logical :: History_Panels3_win_flag, History_Panels3_file_flag
logical :: do_History_Panels3_win, do_History_Panels3_file
integer :: id_History_Panels3_win, id_History_Panels3_file, History_Panels3_file_interval
character (len = strlen) :: History_Panels3_file_dir, History_Panels3_file_prefix
logical :: History_Panels3_xaxis_reversed, History_Panels3_automatic_star_age_units, History_Panels3_xaxis_log
character (len = strlen) :: History_Panels3_xaxis_name, History_Panels3_title
real :: &
      History_Panels3_xmin, History_Panels3_xmax, History_Panels3_dxmin, &
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
      History_Panels3_yaxis_log, History_Panels3_other_yaxis_log, &
      History_Panels3_same_yaxis_range
real, dimension(max_num_History_Panels) :: &
      History_Panels3_ymin, History_Panels3_other_ymin, &
      History_Panels3_ymax, History_Panels3_other_ymax, &
      History_Panels3_ymargin, History_Panels3_other_ymargin, &
      History_Panels3_dymin, History_Panels3_other_dymin

logical :: History_Panels4_win_flag, History_Panels4_file_flag
logical :: do_History_Panels4_win, do_History_Panels4_file
integer :: id_History_Panels4_win, id_History_Panels4_file, History_Panels4_file_interval
character (len = strlen) :: History_Panels4_file_dir, History_Panels4_file_prefix
logical :: History_Panels4_xaxis_reversed, History_Panels4_automatic_star_age_units, History_Panels4_xaxis_log
character (len = strlen) :: History_Panels4_xaxis_name, History_Panels4_title
real :: &
      History_Panels4_xmin, History_Panels4_xmax, History_Panels4_dxmin, &
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
      History_Panels4_yaxis_log, History_Panels4_other_yaxis_log, &
      History_Panels4_same_yaxis_range
real, dimension(max_num_History_Panels) :: &
      History_Panels4_ymin, History_Panels4_other_ymin, &
      History_Panels4_ymax, History_Panels4_other_ymax, &
      History_Panels4_ymargin, History_Panels4_other_ymargin, &
      History_Panels4_dymin, History_Panels4_other_dymin

logical :: History_Panels5_win_flag, History_Panels5_file_flag
logical :: do_History_Panels5_win, do_History_Panels5_file
integer :: id_History_Panels5_win, id_History_Panels5_file, History_Panels5_file_interval
character (len = strlen) :: History_Panels5_file_dir, History_Panels5_file_prefix
logical :: History_Panels5_xaxis_reversed, History_Panels5_automatic_star_age_units, History_Panels5_xaxis_log
character (len = strlen) :: History_Panels5_xaxis_name, History_Panels5_title
real :: &
      History_Panels5_xmin, History_Panels5_xmax, History_Panels5_dxmin, &
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
      History_Panels5_yaxis_log, History_Panels5_other_yaxis_log, &
      History_Panels5_same_yaxis_range
real, dimension(max_num_History_Panels) :: &
      History_Panels5_ymin, History_Panels5_other_ymin, &
      History_Panels5_ymax, History_Panels5_other_ymax, &
      History_Panels5_ymargin, History_Panels5_other_ymargin, &
      History_Panels5_dymin, History_Panels5_other_dymin

logical :: History_Panels6_win_flag, History_Panels6_file_flag
logical :: do_History_Panels6_win, do_History_Panels6_file
integer :: id_History_Panels6_win, id_History_Panels6_file, History_Panels6_file_interval
character (len = strlen) :: History_Panels6_file_dir, History_Panels6_file_prefix
logical :: History_Panels6_xaxis_reversed, History_Panels6_automatic_star_age_units, History_Panels6_xaxis_log
character (len = strlen) :: History_Panels6_xaxis_name, History_Panels6_title
real :: &
      History_Panels6_xmin, History_Panels6_xmax, History_Panels6_dxmin, &
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
      History_Panels6_yaxis_log, History_Panels6_other_yaxis_log, &
      History_Panels6_same_yaxis_range
real, dimension(max_num_History_Panels) :: &
      History_Panels6_ymin, History_Panels6_other_ymin, &
      History_Panels6_ymax, History_Panels6_other_ymax, &
      History_Panels6_ymargin, History_Panels6_other_ymargin, &
      History_Panels6_dymin, History_Panels6_other_dymin

logical :: History_Panels7_win_flag, History_Panels7_file_flag
logical :: do_History_Panels7_win, do_History_Panels7_file
integer :: id_History_Panels7_win, id_History_Panels7_file, History_Panels7_file_interval
character (len = strlen) :: History_Panels7_file_dir, History_Panels7_file_prefix
logical :: History_Panels7_xaxis_reversed, History_Panels7_automatic_star_age_units, History_Panels7_xaxis_log
character (len = strlen) :: History_Panels7_xaxis_name, History_Panels7_title
real :: &
      History_Panels7_xmin, History_Panels7_xmax, History_Panels7_dxmin, &
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
      History_Panels7_yaxis_log, History_Panels7_other_yaxis_log, &
      History_Panels7_same_yaxis_range
real, dimension(max_num_History_Panels) :: &
      History_Panels7_ymin, History_Panels7_other_ymin, &
      History_Panels7_ymax, History_Panels7_other_ymax, &
      History_Panels7_ymargin, History_Panels7_other_ymargin, &
      History_Panels7_dymin, History_Panels7_other_dymin

logical :: History_Panels8_win_flag, History_Panels8_file_flag
logical :: do_History_Panels8_win, do_History_Panels8_file
integer :: id_History_Panels8_win, id_History_Panels8_file, History_Panels8_file_interval
character (len = strlen) :: History_Panels8_file_dir, History_Panels8_file_prefix
logical :: History_Panels8_xaxis_reversed, History_Panels8_automatic_star_age_units, History_Panels8_xaxis_log
character (len = strlen) :: History_Panels8_xaxis_name, History_Panels8_title
real :: &
      History_Panels8_xmin, History_Panels8_xmax, History_Panels8_dxmin, &
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
      History_Panels8_yaxis_log, History_Panels8_other_yaxis_log, &
      History_Panels8_same_yaxis_range
real, dimension(max_num_History_Panels) :: &
      History_Panels8_ymin, History_Panels8_other_ymin, &
      History_Panels8_ymax, History_Panels8_other_ymax, &
      History_Panels8_ymargin, History_Panels8_other_ymargin, &
      History_Panels8_dymin, History_Panels8_other_dymin

logical :: History_Panels9_win_flag, History_Panels9_file_flag
logical :: do_History_Panels9_win, do_History_Panels9_file
integer :: id_History_Panels9_win, id_History_Panels9_file, History_Panels9_file_interval
character (len = strlen) :: History_Panels9_file_dir, History_Panels9_file_prefix
logical :: History_Panels9_xaxis_reversed, History_Panels9_automatic_star_age_units, History_Panels9_xaxis_log
character (len = strlen) :: History_Panels9_xaxis_name, History_Panels9_title
real :: &
      History_Panels9_xmin, History_Panels9_xmax, History_Panels9_dxmin, &
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
      History_Panels9_yaxis_log, History_Panels9_other_yaxis_log, &
      History_Panels9_same_yaxis_range
real, dimension(max_num_History_Panels) :: &
      History_Panels9_ymin, History_Panels9_other_ymin, &
      History_Panels9_ymax, History_Panels9_other_ymax, &
      History_Panels9_ymargin, History_Panels9_other_ymargin, &
      History_Panels9_dymin, History_Panels9_other_dymin
logical :: History_Panel_points_error_bars
integer :: History_Panel_points_interval, History_Panel_points_marker, &
      History_Panel_points_ci, History_Panel_points_lw
real :: History_Panel_points_ch

logical :: Color_Magnitude1_win_flag, Color_Magnitude1_file_flag
logical :: do_Color_Magnitude1_win, do_Color_Magnitude1_file
integer :: id_Color_Magnitude1_win, id_Color_Magnitude1_file, Color_Magnitude1_file_interval
character (len = strlen) :: Color_Magnitude1_file_dir, Color_Magnitude1_file_prefix
logical :: Color_Magnitude1_xaxis_reversed, Color_Magnitude1_xaxis_log
character (len = strlen) :: Color_Magnitude1_xaxis1_name, Color_Magnitude1_xaxis2_name, &
      Color_Magnitude1_title
real :: &
      Color_Magnitude1_xmin, Color_Magnitude1_xmax, Color_Magnitude1_dxmin, &
      Color_Magnitude1_max_width, Color_Magnitude1_xmargin, &
      Color_Magnitude1_win_width, Color_Magnitude1_win_aspect_ratio, &
      Color_Magnitude1_xleft, Color_Magnitude1_xright, &
      Color_Magnitude1_ybot, Color_Magnitude1_ytop, Color_Magnitude1_txt_scale, &
      prev_Color_Magnitude1_win_width, prev_Color_Magnitude1_win_ratio, &
      Color_Magnitude1_file_width, Color_Magnitude1_file_aspect_ratio, &
      prev_Color_Magnitude1_file_width, prev_Color_Magnitude1_file_ratio
integer :: Color_Magnitude1_num_panels
character (len = strlen), dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude1_yaxis1_name, Color_Magnitude1_other_yaxis1_name, &
      Color_Magnitude1_yaxis2_name, Color_Magnitude1_other_yaxis2_name
logical, dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude1_yaxis_reversed, Color_Magnitude1_other_yaxis_reversed, &
      Color_Magnitude1_yaxis_log, Color_Magnitude1_other_yaxis_log
real, dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude1_ymin, Color_Magnitude1_other_ymin, &
      Color_Magnitude1_ymax, Color_Magnitude1_other_ymax, &
      Color_Magnitude1_ymargin, Color_Magnitude1_other_ymargin, &
      Color_Magnitude1_dymin, Color_Magnitude1_other_dymin

logical :: Color_Magnitude2_win_flag, Color_Magnitude2_file_flag
logical :: do_Color_Magnitude2_win, do_Color_Magnitude2_file
integer :: id_Color_Magnitude2_win, id_Color_Magnitude2_file, Color_Magnitude2_file_interval
character (len = strlen) :: Color_Magnitude2_file_dir, Color_Magnitude2_file_prefix
logical :: Color_Magnitude2_xaxis_reversed, Color_Magnitude2_xaxis_log
character (len = strlen) :: Color_Magnitude2_xaxis1_name, Color_Magnitude2_xaxis2_name, &
      Color_Magnitude2_title
real :: &
      Color_Magnitude2_xmin, Color_Magnitude2_xmax, Color_Magnitude2_dxmin, &
      Color_Magnitude2_max_width, Color_Magnitude2_xmargin, &
      Color_Magnitude2_win_width, Color_Magnitude2_win_aspect_ratio, &
      Color_Magnitude2_xleft, Color_Magnitude2_xright, &
      Color_Magnitude2_ybot, Color_Magnitude2_ytop, Color_Magnitude2_txt_scale, &
      prev_Color_Magnitude2_win_width, prev_Color_Magnitude2_win_ratio, &
      Color_Magnitude2_file_width, Color_Magnitude2_file_aspect_ratio, &
      prev_Color_Magnitude2_file_width, prev_Color_Magnitude2_file_ratio
integer :: Color_Magnitude2_num_panels
character (len = strlen), dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude2_yaxis1_name, Color_Magnitude2_other_yaxis1_name, &
      Color_Magnitude2_yaxis2_name, Color_Magnitude2_other_yaxis2_name
logical, dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude2_yaxis_reversed, Color_Magnitude2_other_yaxis_reversed, &
      Color_Magnitude2_yaxis_log, Color_Magnitude2_other_yaxis_log
real, dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude2_ymin, Color_Magnitude2_other_ymin, &
      Color_Magnitude2_ymax, Color_Magnitude2_other_ymax, &
      Color_Magnitude2_ymargin, Color_Magnitude2_other_ymargin, &
      Color_Magnitude2_dymin, Color_Magnitude2_other_dymin

logical :: Color_Magnitude3_win_flag, Color_Magnitude3_file_flag
logical :: do_Color_Magnitude3_win, do_Color_Magnitude3_file
integer :: id_Color_Magnitude3_win, id_Color_Magnitude3_file, Color_Magnitude3_file_interval
character (len = strlen) :: Color_Magnitude3_file_dir, Color_Magnitude3_file_prefix
logical :: Color_Magnitude3_xaxis_reversed, Color_Magnitude3_xaxis_log
character (len = strlen) :: Color_Magnitude3_xaxis1_name, Color_Magnitude3_xaxis2_name, &
      Color_Magnitude3_title
real :: &
      Color_Magnitude3_xmin, Color_Magnitude3_xmax, Color_Magnitude3_dxmin, &
      Color_Magnitude3_max_width, Color_Magnitude3_xmargin, &
      Color_Magnitude3_win_width, Color_Magnitude3_win_aspect_ratio, &
      Color_Magnitude3_xleft, Color_Magnitude3_xright, &
      Color_Magnitude3_ybot, Color_Magnitude3_ytop, Color_Magnitude3_txt_scale, &
      prev_Color_Magnitude3_win_width, prev_Color_Magnitude3_win_ratio, &
      Color_Magnitude3_file_width, Color_Magnitude3_file_aspect_ratio, &
      prev_Color_Magnitude3_file_width, prev_Color_Magnitude3_file_ratio
integer :: Color_Magnitude3_num_panels
character (len = strlen), dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude3_yaxis1_name, Color_Magnitude3_other_yaxis1_name, &
      Color_Magnitude3_yaxis2_name, Color_Magnitude3_other_yaxis2_name
logical, dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude3_yaxis_reversed, Color_Magnitude3_other_yaxis_reversed, &
      Color_Magnitude3_yaxis_log, Color_Magnitude3_other_yaxis_log
real, dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude3_ymin, Color_Magnitude3_other_ymin, &
      Color_Magnitude3_ymax, Color_Magnitude3_other_ymax, &
      Color_Magnitude3_ymargin, Color_Magnitude3_other_ymargin, &
      Color_Magnitude3_dymin, Color_Magnitude3_other_dymin

logical :: Color_Magnitude4_win_flag, Color_Magnitude4_file_flag
logical :: do_Color_Magnitude4_win, do_Color_Magnitude4_file
integer :: id_Color_Magnitude4_win, id_Color_Magnitude4_file, Color_Magnitude4_file_interval
character (len = strlen) :: Color_Magnitude4_file_dir, Color_Magnitude4_file_prefix
logical :: Color_Magnitude4_xaxis_reversed, Color_Magnitude4_xaxis_log
character (len = strlen) :: Color_Magnitude4_xaxis1_name, Color_Magnitude4_xaxis2_name, &
      Color_Magnitude4_title
real :: &
      Color_Magnitude4_xmin, Color_Magnitude4_xmax, Color_Magnitude4_dxmin, &
      Color_Magnitude4_max_width, Color_Magnitude4_xmargin, &
      Color_Magnitude4_win_width, Color_Magnitude4_win_aspect_ratio, &
      Color_Magnitude4_xleft, Color_Magnitude4_xright, &
      Color_Magnitude4_ybot, Color_Magnitude4_ytop, Color_Magnitude4_txt_scale, &
      prev_Color_Magnitude4_win_width, prev_Color_Magnitude4_win_ratio, &
      Color_Magnitude4_file_width, Color_Magnitude4_file_aspect_ratio, &
      prev_Color_Magnitude4_file_width, prev_Color_Magnitude4_file_ratio
integer :: Color_Magnitude4_num_panels
character (len = strlen), dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude4_yaxis1_name, Color_Magnitude4_other_yaxis1_name, &
      Color_Magnitude4_yaxis2_name, Color_Magnitude4_other_yaxis2_name
logical, dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude4_yaxis_reversed, Color_Magnitude4_other_yaxis_reversed, &
      Color_Magnitude4_yaxis_log, Color_Magnitude4_other_yaxis_log
real, dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude4_ymin, Color_Magnitude4_other_ymin, &
      Color_Magnitude4_ymax, Color_Magnitude4_other_ymax, &
      Color_Magnitude4_ymargin, Color_Magnitude4_other_ymargin, &
      Color_Magnitude4_dymin, Color_Magnitude4_other_dymin

logical :: Color_Magnitude5_win_flag, Color_Magnitude5_file_flag
logical :: do_Color_Magnitude5_win, do_Color_Magnitude5_file
integer :: id_Color_Magnitude5_win, id_Color_Magnitude5_file, Color_Magnitude5_file_interval
character (len = strlen) :: Color_Magnitude5_file_dir, Color_Magnitude5_file_prefix
logical :: Color_Magnitude5_xaxis_reversed, Color_Magnitude5_xaxis_log
character (len = strlen) :: Color_Magnitude5_xaxis1_name, Color_Magnitude5_xaxis2_name, &
      Color_Magnitude5_title
real :: &
      Color_Magnitude5_xmin, Color_Magnitude5_xmax, Color_Magnitude5_dxmin, &
      Color_Magnitude5_max_width, Color_Magnitude5_xmargin, &
      Color_Magnitude5_win_width, Color_Magnitude5_win_aspect_ratio, &
      Color_Magnitude5_xleft, Color_Magnitude5_xright, &
      Color_Magnitude5_ybot, Color_Magnitude5_ytop, Color_Magnitude5_txt_scale, &
      prev_Color_Magnitude5_win_width, prev_Color_Magnitude5_win_ratio, &
      Color_Magnitude5_file_width, Color_Magnitude5_file_aspect_ratio, &
      prev_Color_Magnitude5_file_width, prev_Color_Magnitude5_file_ratio
integer :: Color_Magnitude5_num_panels
character (len = strlen), dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude5_yaxis1_name, Color_Magnitude5_other_yaxis1_name, &
      Color_Magnitude5_yaxis2_name, Color_Magnitude5_other_yaxis2_name
logical, dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude5_yaxis_reversed, Color_Magnitude5_other_yaxis_reversed, &
      Color_Magnitude5_yaxis_log, Color_Magnitude5_other_yaxis_log
real, dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude5_ymin, Color_Magnitude5_other_ymin, &
      Color_Magnitude5_ymax, Color_Magnitude5_other_ymax, &
      Color_Magnitude5_ymargin, Color_Magnitude5_other_ymargin, &
      Color_Magnitude5_dymin, Color_Magnitude5_other_dymin

logical :: Color_Magnitude6_win_flag, Color_Magnitude6_file_flag
logical :: do_Color_Magnitude6_win, do_Color_Magnitude6_file
integer :: id_Color_Magnitude6_win, id_Color_Magnitude6_file, Color_Magnitude6_file_interval
character (len = strlen) :: Color_Magnitude6_file_dir, Color_Magnitude6_file_prefix
logical :: Color_Magnitude6_xaxis_reversed, Color_Magnitude6_xaxis_log
character (len = strlen) :: Color_Magnitude6_xaxis1_name, Color_Magnitude6_xaxis2_name, &
      Color_Magnitude6_title
real :: &
      Color_Magnitude6_xmin, Color_Magnitude6_xmax, Color_Magnitude6_dxmin, &
      Color_Magnitude6_max_width, Color_Magnitude6_xmargin, &
      Color_Magnitude6_win_width, Color_Magnitude6_win_aspect_ratio, &
      Color_Magnitude6_xleft, Color_Magnitude6_xright, &
      Color_Magnitude6_ybot, Color_Magnitude6_ytop, Color_Magnitude6_txt_scale, &
      prev_Color_Magnitude6_win_width, prev_Color_Magnitude6_win_ratio, &
      Color_Magnitude6_file_width, Color_Magnitude6_file_aspect_ratio, &
      prev_Color_Magnitude6_file_width, prev_Color_Magnitude6_file_ratio
integer :: Color_Magnitude6_num_panels
character (len = strlen), dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude6_yaxis1_name, Color_Magnitude6_other_yaxis1_name, &
      Color_Magnitude6_yaxis2_name, Color_Magnitude6_other_yaxis2_name
logical, dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude6_yaxis_reversed, Color_Magnitude6_other_yaxis_reversed, &
      Color_Magnitude6_yaxis_log, Color_Magnitude6_other_yaxis_log
real, dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude6_ymin, Color_Magnitude6_other_ymin, &
      Color_Magnitude6_ymax, Color_Magnitude6_other_ymax, &
      Color_Magnitude6_ymargin, Color_Magnitude6_other_ymargin, &
      Color_Magnitude6_dymin, Color_Magnitude6_other_dymin

logical :: Color_Magnitude7_win_flag, Color_Magnitude7_file_flag
logical :: do_Color_Magnitude7_win, do_Color_Magnitude7_file
integer :: id_Color_Magnitude7_win, id_Color_Magnitude7_file, Color_Magnitude7_file_interval
character (len = strlen) :: Color_Magnitude7_file_dir, Color_Magnitude7_file_prefix
logical :: Color_Magnitude7_xaxis_reversed, Color_Magnitude7_xaxis_log
character (len = strlen) :: Color_Magnitude7_xaxis1_name, Color_Magnitude7_xaxis2_name, &
      Color_Magnitude7_title
real :: &
      Color_Magnitude7_xmin, Color_Magnitude7_xmax, Color_Magnitude7_dxmin, &
      Color_Magnitude7_max_width, Color_Magnitude7_xmargin, &
      Color_Magnitude7_win_width, Color_Magnitude7_win_aspect_ratio, &
      Color_Magnitude7_xleft, Color_Magnitude7_xright, &
      Color_Magnitude7_ybot, Color_Magnitude7_ytop, Color_Magnitude7_txt_scale, &
      prev_Color_Magnitude7_win_width, prev_Color_Magnitude7_win_ratio, &
      Color_Magnitude7_file_width, Color_Magnitude7_file_aspect_ratio, &
      prev_Color_Magnitude7_file_width, prev_Color_Magnitude7_file_ratio
integer :: Color_Magnitude7_num_panels
character (len = strlen), dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude7_yaxis1_name, Color_Magnitude7_other_yaxis1_name, &
      Color_Magnitude7_yaxis2_name, Color_Magnitude7_other_yaxis2_name
logical, dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude7_yaxis_reversed, Color_Magnitude7_other_yaxis_reversed, &
      Color_Magnitude7_yaxis_log, Color_Magnitude7_other_yaxis_log
real, dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude7_ymin, Color_Magnitude7_other_ymin, &
      Color_Magnitude7_ymax, Color_Magnitude7_other_ymax, &
      Color_Magnitude7_ymargin, Color_Magnitude7_other_ymargin, &
      Color_Magnitude7_dymin, Color_Magnitude7_other_dymin

logical :: Color_Magnitude8_win_flag, Color_Magnitude8_file_flag
logical :: do_Color_Magnitude8_win, do_Color_Magnitude8_file
integer :: id_Color_Magnitude8_win, id_Color_Magnitude8_file, Color_Magnitude8_file_interval
character (len = strlen) :: Color_Magnitude8_file_dir, Color_Magnitude8_file_prefix
logical :: Color_Magnitude8_xaxis_reversed, Color_Magnitude8_xaxis_log
character (len = strlen) :: Color_Magnitude8_xaxis1_name, Color_Magnitude8_xaxis2_name, &
      Color_Magnitude8_title
real :: &
      Color_Magnitude8_xmin, Color_Magnitude8_xmax, Color_Magnitude8_dxmin, &
      Color_Magnitude8_max_width, Color_Magnitude8_xmargin, &
      Color_Magnitude8_win_width, Color_Magnitude8_win_aspect_ratio, &
      Color_Magnitude8_xleft, Color_Magnitude8_xright, &
      Color_Magnitude8_ybot, Color_Magnitude8_ytop, Color_Magnitude8_txt_scale, &
      prev_Color_Magnitude8_win_width, prev_Color_Magnitude8_win_ratio, &
      Color_Magnitude8_file_width, Color_Magnitude8_file_aspect_ratio, &
      prev_Color_Magnitude8_file_width, prev_Color_Magnitude8_file_ratio
integer :: Color_Magnitude8_num_panels
character (len = strlen), dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude8_yaxis1_name, Color_Magnitude8_other_yaxis1_name, &
      Color_Magnitude8_yaxis2_name, Color_Magnitude8_other_yaxis2_name
logical, dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude8_yaxis_reversed, Color_Magnitude8_other_yaxis_reversed, &
      Color_Magnitude8_yaxis_log, Color_Magnitude8_other_yaxis_log
real, dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude8_ymin, Color_Magnitude8_other_ymin, &
      Color_Magnitude8_ymax, Color_Magnitude8_other_ymax, &
      Color_Magnitude8_ymargin, Color_Magnitude8_other_ymargin, &
      Color_Magnitude8_dymin, Color_Magnitude8_other_dymin

logical :: Color_Magnitude9_win_flag, Color_Magnitude9_file_flag
logical :: do_Color_Magnitude9_win, do_Color_Magnitude9_file
integer :: id_Color_Magnitude9_win, id_Color_Magnitude9_file, Color_Magnitude9_file_interval
character (len = strlen) :: Color_Magnitude9_file_dir, Color_Magnitude9_file_prefix
logical :: Color_Magnitude9_xaxis_reversed, Color_Magnitude9_xaxis_log
character (len = strlen) :: Color_Magnitude9_xaxis1_name, Color_Magnitude9_xaxis2_name, &
      Color_Magnitude9_title
real :: &
      Color_Magnitude9_xmin, Color_Magnitude9_xmax, Color_Magnitude9_dxmin, &
      Color_Magnitude9_max_width, Color_Magnitude9_xmargin, &
      Color_Magnitude9_win_width, Color_Magnitude9_win_aspect_ratio, &
      Color_Magnitude9_xleft, Color_Magnitude9_xright, &
      Color_Magnitude9_ybot, Color_Magnitude9_ytop, Color_Magnitude9_txt_scale, &
      prev_Color_Magnitude9_win_width, prev_Color_Magnitude9_win_ratio, &
      Color_Magnitude9_file_width, Color_Magnitude9_file_aspect_ratio, &
      prev_Color_Magnitude9_file_width, prev_Color_Magnitude9_file_ratio
integer :: Color_Magnitude9_num_panels
character (len = strlen), dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude9_yaxis1_name, Color_Magnitude9_other_yaxis1_name, &
      Color_Magnitude9_yaxis2_name, Color_Magnitude9_other_yaxis2_name
logical, dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude9_yaxis_reversed, Color_Magnitude9_other_yaxis_reversed, &
      Color_Magnitude9_yaxis_log, Color_Magnitude9_other_yaxis_log
real, dimension(max_num_Color_Magnitude) :: &
      Color_Magnitude9_ymin, Color_Magnitude9_other_ymin, &
      Color_Magnitude9_ymax, Color_Magnitude9_other_ymax, &
      Color_Magnitude9_ymargin, Color_Magnitude9_other_ymargin, &
      Color_Magnitude9_dymin, Color_Magnitude9_other_dymin

logical :: Mode_Prop_win_flag, Mode_Prop_file_flag
logical :: do_Mode_Prop_win, do_Mode_Prop_file, Mode_Prop_xaxis_reversed
integer :: id_Mode_Prop_win, id_Mode_Prop_file, Mode_Prop_file_interval
character (len = strlen) :: Mode_Prop_file_dir, Mode_Prop_file_prefix, &
      Mode_Prop_xaxis_name, Mode_Prop_title
real :: &
      Mode_Prop_nu_max_obs, Mode_Prop_xmin, Mode_Prop_xmax, Mode_Prop_ymin, Mode_Prop_ymax, &
      Mode_Prop_win_width, Mode_Prop_win_aspect_ratio, &
      Mode_Prop_xleft, Mode_Prop_xright, &
      Mode_Prop_ybot, Mode_Prop_ytop, Mode_Prop_txt_scale, &
      prev_Mode_Prop_win_width, prev_Mode_Prop_win_ratio, &
      Mode_Prop_file_width, Mode_Prop_file_aspect_ratio, &
      prev_Mode_Prop_file_width, prev_Mode_Prop_file_ratio

logical :: Summary_Burn_win_flag, Summary_Burn_file_flag
logical :: do_Summary_Burn_win, do_Summary_Burn_file, &
      Summary_Burn_xaxis_reversed
integer :: id_Summary_Burn_win, id_Summary_Burn_file, &
      Summary_Burn_num_lines, Summary_Burn_file_interval
character (len = strlen) :: Summary_Burn_file_dir, Summary_Burn_file_prefix, &
      Summary_Burn_xaxis_name, Summary_Burn_title
real :: &
      Summary_Burn_xmin, Summary_Burn_xmax, Summary_Burn_title_shift, &
      Summary_Burn_win_width, Summary_Burn_win_aspect_ratio, &
      Summary_Burn_xleft, Summary_Burn_xright, &
      Summary_Burn_ybot, Summary_Burn_ytop, Summary_Burn_txt_scale, &
      prev_Summary_Burn_win_width, prev_Summary_Burn_win_ratio, &
      Summary_Burn_file_width, Summary_Burn_file_aspect_ratio, &
      prev_Summary_Burn_file_width, prev_Summary_Burn_file_ratio

logical :: Summary_Profile_win_flag, Summary_Profile_file_flag
logical :: do_Summary_Profile_win, do_Summary_Profile_file, &
      Summary_Profile_xaxis_reversed, &
      Summary_Profile_scaled_value(max_num_Summary_Profile_Lines)
integer :: id_Summary_Profile_win, id_Summary_Profile_file, &
      Summary_Profile_num_lines, Summary_Profile_file_interval
character (len = strlen) :: Summary_Profile_file_dir, Summary_Profile_file_prefix, &
      Summary_Profile_xaxis_name, Summary_Profile_title, &
      Summary_Profile_name(max_num_Summary_Profile_Lines), &
      Summary_Profile_legend(max_num_Summary_Profile_Lines)
real :: &
      Summary_Profile_xmin, Summary_Profile_xmax, &
      Summary_Profile_win_width, Summary_Profile_win_aspect_ratio, &
      Summary_Profile_xleft, Summary_Profile_xright, &
      Summary_Profile_ybot, Summary_Profile_ytop, Summary_Profile_txt_scale, &
      prev_Summary_Profile_win_width, prev_Summary_Profile_win_ratio, &
      Summary_Profile_file_width, Summary_Profile_file_aspect_ratio, &
      prev_Summary_Profile_file_width, prev_Summary_Profile_file_ratio

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
real :: &
      Summary_History_win_width, Summary_History_win_aspect_ratio, &
      Summary_History_xleft, Summary_History_xright, &
      Summary_History_ybot, Summary_History_ytop, Summary_History_txt_scale, &
      prev_Summary_History_win_width, prev_Summary_History_win_ratio, &
      Summary_History_file_width, Summary_History_file_aspect_ratio, &
      prev_Summary_History_file_width, prev_Summary_History_file_ratio

logical :: Power_win_flag, Power_file_flag
logical :: do_Power_win, do_Power_file, Power_xaxis_reversed
integer :: id_Power_win, id_Power_file, Power_file_interval, &
      Power_legend_max_cnt
character (len = strlen) :: Power_file_dir, Power_file_prefix, &
      Power_xaxis_name, Power_title
real :: &
      Power_xmin, Power_xmax, Power_ymin, Power_ymax, &
      Power_win_width, Power_win_aspect_ratio, &
      Power_xleft, Power_xright, Power_legend_txt_scale_factor, &
      Power_ybot, Power_ytop, Power_txt_scale, &
      prev_Power_win_width, prev_Power_win_ratio, &
      Power_file_width, Power_file_aspect_ratio, &
      prev_Power_file_width, prev_Power_file_ratio

logical :: Abundance_win_flag, Abundance_file_flag
logical :: do_Abundance_win, do_Abundance_file, Abundance_xaxis_reversed, &
      Abundance_show_photosphere_location
integer :: id_Abundance_win, id_Abundance_file, Abundance_file_interval
integer :: Abundance_num_isos_to_show, num_abundance_line_labels, &
      Abundance_legend_max_cnt
character(len = iso_name_length) :: &
      Abundance_which_isos_to_show(max_Abundance_num_isos_to_show)
character (len = strlen) :: Abundance_file_dir, Abundance_file_prefix, &
      Abundance_xaxis_name, Abundance_title
real :: &
      Abundance_xmin, Abundance_xmax, Abundance_legend_txt_scale_factor, &
      Abundance_line_txt_scale_factor, &
      Abundance_log_mass_frac_min, Abundance_log_mass_frac_max, &
      Abundance_win_width, Abundance_win_aspect_ratio, &
      Abundance_xleft, Abundance_xright, &
      Abundance_ybot, Abundance_ytop, Abundance_txt_scale, &
      prev_Abundance_win_width, prev_Abundance_win_ratio, &
      Abundance_file_width, Abundance_file_aspect_ratio, &
      prev_Abundance_file_width, prev_Abundance_file_ratio

logical :: Network_win_flag, Network_file_flag
logical :: do_Network_win, do_Network_file, Network_show_colorbar
logical :: Network_show_element_names, Network_show_mass_fraction
integer :: id_Network_win, id_Network_file, Network_file_interval
character (len = strlen) :: Network_file_dir, Network_file_prefix, &
      Network_title
real :: &
      Network_zmin, Network_zmax, &
      Network_nmin, Network_nmax, &
      Network_log_mass_frac_min, Network_log_mass_frac_max, &
      Network_win_width, Network_win_aspect_ratio, &
      Network_xleft, Network_xright, &
      Network_ybot, Network_ytop, Network_txt_scale, &
      prev_Network_win_width, prev_Network_win_ratio, &
      Network_file_width, Network_file_aspect_ratio, &
      prev_Network_file_width, prev_Network_file_ratio

logical :: Production_win_flag, Production_file_flag
logical :: do_Production_win, do_Production_file
logical :: Production_show_element_names
integer :: id_Production_win, id_Production_file, Production_file_interval
character (len = strlen) :: Production_file_dir, Production_file_prefix, &
      Production_title
real :: &
      Production_amin, Production_amax, &
      Production_ymin, Production_ymax, &
      Production_min_mass, Production_max_mass, &
      Production_min_mass_frac, &
      Production_win_width, Production_win_aspect_ratio, &
      Production_xleft, Production_xright, &
      Production_ybot, Production_ytop, Production_txt_scale, &
      prev_Production_win_width, prev_Production_win_ratio, &
      Production_file_width, Production_file_aspect_ratio, &
      prev_Production_file_width, prev_Production_file_ratio

logical :: Grid1_win_flag, Grid1_file_flag
logical :: do_Grid1_win, do_Grid1_file
integer :: Grid1_file_interval, Grid1_step_max
integer :: id_Grid1_win, id_Grid1_file
character (len = strlen) :: Grid1_file_dir, Grid1_file_prefix
real :: &
      Grid1_win_width, Grid1_win_aspect_ratio, &
      Grid1_xleft, Grid1_xright, &
      Grid1_ybot, Grid1_ytop, &
      prev_Grid1_win_width, prev_Grid1_win_ratio, &
      Grid1_file_width, Grid1_file_aspect_ratio, &
      prev_Grid1_file_width, prev_Grid1_file_ratio
character (len = strlen) :: Grid1_title, &
      Grid1_plot_name(max_num_pgstar_grid_plots)
real, dimension(max_num_pgstar_grid_plots) :: &
      Grid1_txt_scale_factor, &
      Grid1_plot_pad_left, Grid1_plot_pad_right, &
      Grid1_plot_pad_top, Grid1_plot_pad_bot
integer :: &
      Grid1_num_cols, Grid1_num_rows, Grid1_num_plots
integer, dimension(max_num_pgstar_grid_plots) :: &
      Grid1_plot_row, Grid1_plot_rowspan, &
      Grid1_plot_col, Grid1_plot_colspan

logical :: Grid2_win_flag, Grid2_file_flag
logical :: do_Grid2_win, do_Grid2_file
integer :: Grid2_file_interval, Grid2_step_max
integer :: id_Grid2_win, id_Grid2_file
character (len = strlen) :: Grid2_file_dir, Grid2_file_prefix
real :: &
      Grid2_win_width, Grid2_win_aspect_ratio, &
      Grid2_xleft, Grid2_xright, &
      Grid2_ybot, Grid2_ytop, &
      prev_Grid2_win_width, prev_Grid2_win_ratio, &
      Grid2_file_width, Grid2_file_aspect_ratio, &
      prev_Grid2_file_width, prev_Grid2_file_ratio
character (len = strlen) :: Grid2_title, &
      Grid2_plot_name(max_num_pgstar_grid_plots)
real, dimension(max_num_pgstar_grid_plots) :: &
      Grid2_txt_scale_factor, &
      Grid2_plot_pad_left, Grid2_plot_pad_right, &
      Grid2_plot_pad_top, Grid2_plot_pad_bot
integer :: &
      Grid2_num_cols, Grid2_num_rows, Grid2_num_plots
integer, dimension(max_num_pgstar_grid_plots) :: &
      Grid2_plot_row, Grid2_plot_rowspan, &
      Grid2_plot_col, Grid2_plot_colspan

logical :: Grid3_win_flag, Grid3_file_flag
logical :: do_Grid3_win, do_Grid3_file
integer :: Grid3_file_interval, Grid3_step_max
integer :: id_Grid3_win, id_Grid3_file
character (len = strlen) :: Grid3_file_dir, Grid3_file_prefix
real :: &
      Grid3_win_width, Grid3_win_aspect_ratio, &
      Grid3_xleft, Grid3_xright, &
      Grid3_ybot, Grid3_ytop, &
      prev_Grid3_win_width, prev_Grid3_win_ratio, &
      Grid3_file_width, Grid3_file_aspect_ratio, &
      prev_Grid3_file_width, prev_Grid3_file_ratio
character (len = strlen) :: Grid3_title, &
      Grid3_plot_name(max_num_pgstar_grid_plots)
real, dimension(max_num_pgstar_grid_plots) :: &
      Grid3_txt_scale_factor, &
      Grid3_plot_pad_left, Grid3_plot_pad_right, &
      Grid3_plot_pad_top, Grid3_plot_pad_bot
integer :: &
      Grid3_num_cols, Grid3_num_rows, Grid3_num_plots
integer, dimension(max_num_pgstar_grid_plots) :: &
      Grid3_plot_row, Grid3_plot_rowspan, &
      Grid3_plot_col, Grid3_plot_colspan

logical :: Grid4_win_flag, Grid4_file_flag
logical :: do_Grid4_win, do_Grid4_file
integer :: Grid4_file_interval, Grid4_step_max
integer :: id_Grid4_win, id_Grid4_file
character (len = strlen) :: Grid4_file_dir, Grid4_file_prefix
real :: &
      Grid4_win_width, Grid4_win_aspect_ratio, &
      Grid4_xleft, Grid4_xright, &
      Grid4_ybot, Grid4_ytop, &
      prev_Grid4_win_width, prev_Grid4_win_ratio, &
      Grid4_file_width, Grid4_file_aspect_ratio, &
      prev_Grid4_file_width, prev_Grid4_file_ratio
character (len = strlen) :: Grid4_title, &
      Grid4_plot_name(max_num_pgstar_grid_plots)
real, dimension(max_num_pgstar_grid_plots) :: &
      Grid4_txt_scale_factor, &
      Grid4_plot_pad_left, Grid4_plot_pad_right, &
      Grid4_plot_pad_top, Grid4_plot_pad_bot
integer :: &
      Grid4_num_cols, Grid4_num_rows, Grid4_num_plots
integer, dimension(max_num_pgstar_grid_plots) :: &
      Grid4_plot_row, Grid4_plot_rowspan, &
      Grid4_plot_col, Grid4_plot_colspan

logical :: Grid5_win_flag, Grid5_file_flag
logical :: do_Grid5_win, do_Grid5_file
integer :: Grid5_file_interval, Grid5_step_max
integer :: id_Grid5_win, id_Grid5_file
character (len = strlen) :: Grid5_file_dir, Grid5_file_prefix
real :: &
      Grid5_win_width, Grid5_win_aspect_ratio, &
      Grid5_xleft, Grid5_xright, &
      Grid5_ybot, Grid5_ytop, &
      prev_Grid5_win_width, prev_Grid5_win_ratio, &
      Grid5_file_width, Grid5_file_aspect_ratio, &
      prev_Grid5_file_width, prev_Grid5_file_ratio
character (len = strlen) :: Grid5_title, &
      Grid5_plot_name(max_num_pgstar_grid_plots)
real, dimension(max_num_pgstar_grid_plots) :: &
      Grid5_txt_scale_factor, &
      Grid5_plot_pad_left, Grid5_plot_pad_right, &
      Grid5_plot_pad_top, Grid5_plot_pad_bot
integer :: &
      Grid5_num_cols, Grid5_num_rows, Grid5_num_plots
integer, dimension(max_num_pgstar_grid_plots) :: &
      Grid5_plot_row, Grid5_plot_rowspan, &
      Grid5_plot_col, Grid5_plot_colspan

logical :: Grid6_win_flag, Grid6_file_flag
logical :: do_Grid6_win, do_Grid6_file
integer :: Grid6_file_interval, Grid6_step_max
integer :: id_Grid6_win, id_Grid6_file
character (len = strlen) :: Grid6_file_dir, Grid6_file_prefix
real :: &
      Grid6_win_width, Grid6_win_aspect_ratio, &
      Grid6_xleft, Grid6_xright, &
      Grid6_ybot, Grid6_ytop, &
      prev_Grid6_win_width, prev_Grid6_win_ratio, &
      Grid6_file_width, Grid6_file_aspect_ratio, &
      prev_Grid6_file_width, prev_Grid6_file_ratio
character (len = strlen) :: Grid6_title, &
      Grid6_plot_name(max_num_pgstar_grid_plots)
real, dimension(max_num_pgstar_grid_plots) :: &
      Grid6_txt_scale_factor, &
      Grid6_plot_pad_left, Grid6_plot_pad_right, &
      Grid6_plot_pad_top, Grid6_plot_pad_bot
integer :: &
      Grid6_num_cols, Grid6_num_rows, Grid6_num_plots
integer, dimension(max_num_pgstar_grid_plots) :: &
      Grid6_plot_row, Grid6_plot_rowspan, &
      Grid6_plot_col, Grid6_plot_colspan

logical :: Grid7_win_flag, Grid7_file_flag
logical :: do_Grid7_win, do_Grid7_file
integer :: Grid7_file_interval, Grid7_step_max
integer :: id_Grid7_win, id_Grid7_file
character (len = strlen) :: Grid7_file_dir, Grid7_file_prefix
real :: &
      Grid7_win_width, Grid7_win_aspect_ratio, &
      Grid7_xleft, Grid7_xright, &
      Grid7_ybot, Grid7_ytop, &
      prev_Grid7_win_width, prev_Grid7_win_ratio, &
      Grid7_file_width, Grid7_file_aspect_ratio, &
      prev_Grid7_file_width, prev_Grid7_file_ratio
character (len = strlen) :: Grid7_title, &
      Grid7_plot_name(max_num_pgstar_grid_plots)
real, dimension(max_num_pgstar_grid_plots) :: &
      Grid7_txt_scale_factor, &
      Grid7_plot_pad_left, Grid7_plot_pad_right, &
      Grid7_plot_pad_top, Grid7_plot_pad_bot
integer :: &
      Grid7_num_cols, Grid7_num_rows, Grid7_num_plots
integer, dimension(max_num_pgstar_grid_plots) :: &
      Grid7_plot_row, Grid7_plot_rowspan, &
      Grid7_plot_col, Grid7_plot_colspan

logical :: Grid8_win_flag, Grid8_file_flag
logical :: do_Grid8_win, do_Grid8_file
integer :: Grid8_file_interval, Grid8_step_max
integer :: id_Grid8_win, id_Grid8_file
character (len = strlen) :: Grid8_file_dir, Grid8_file_prefix
real :: &
      Grid8_win_width, Grid8_win_aspect_ratio, &
      Grid8_xleft, Grid8_xright, &
      Grid8_ybot, Grid8_ytop, &
      prev_Grid8_win_width, prev_Grid8_win_ratio, &
      Grid8_file_width, Grid8_file_aspect_ratio, &
      prev_Grid8_file_width, prev_Grid8_file_ratio
character (len = strlen) :: Grid8_title, &
      Grid8_plot_name(max_num_pgstar_grid_plots)
real, dimension(max_num_pgstar_grid_plots) :: &
      Grid8_txt_scale_factor, &
      Grid8_plot_pad_left, Grid8_plot_pad_right, &
      Grid8_plot_pad_top, Grid8_plot_pad_bot
integer :: &
      Grid8_num_cols, Grid8_num_rows, Grid8_num_plots
integer, dimension(max_num_pgstar_grid_plots) :: &
      Grid8_plot_row, Grid8_plot_rowspan, &
      Grid8_plot_col, Grid8_plot_colspan

logical :: Grid9_win_flag, Grid9_file_flag
logical :: do_Grid9_win, do_Grid9_file
integer :: Grid9_file_interval, Grid9_step_max
integer :: id_Grid9_win, id_Grid9_file
character (len = strlen) :: Grid9_file_dir, Grid9_file_prefix
real :: &
      Grid9_win_width, Grid9_win_aspect_ratio, &
      Grid9_xleft, Grid9_xright, &
      Grid9_ybot, Grid9_ytop, &
      prev_Grid9_win_width, prev_Grid9_win_ratio, &
      Grid9_file_width, Grid9_file_aspect_ratio, &
      prev_Grid9_file_width, prev_Grid9_file_ratio
character (len = strlen) :: Grid9_title, &
      Grid9_plot_name(max_num_pgstar_grid_plots)
real, dimension(max_num_pgstar_grid_plots) :: &
      Grid9_txt_scale_factor, &
      Grid9_plot_pad_left, Grid9_plot_pad_right, &
      Grid9_plot_pad_top, Grid9_plot_pad_bot
integer :: &
      Grid9_num_cols, Grid9_num_rows, Grid9_num_plots
integer, dimension(max_num_pgstar_grid_plots) :: &
      Grid9_plot_row, Grid9_plot_rowspan, &
      Grid9_plot_col, Grid9_plot_colspan

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

logical :: read_extra_pgstar_inlist1
character (len = strlen) :: extra_pgstar_inlist1_name

logical :: read_extra_pgstar_inlist2
character (len = strlen) :: extra_pgstar_inlist2_name

logical :: read_extra_pgstar_inlist3
character (len = strlen) :: extra_pgstar_inlist3_name

logical :: read_extra_pgstar_inlist4
character (len = strlen) :: extra_pgstar_inlist4_name

logical :: read_extra_pgstar_inlist5
character (len = strlen) :: extra_pgstar_inlist5_name

logical :: Abundance_use_decorator, &
      Color_Magnitude1_use_decorator, &
      Color_Magnitude2_use_decorator, &
      Color_Magnitude3_use_decorator, &
      Color_Magnitude4_use_decorator, &
      Color_Magnitude5_use_decorator, &
      Color_Magnitude6_use_decorator, &
      Color_Magnitude7_use_decorator, &
      Color_Magnitude8_use_decorator, &
      Color_Magnitude9_use_decorator, &
      dPg_dnu_use_decorator, &
      Dynamo_use_decorator, &
      History_Panels1_use_decorator, &
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
      HR_use_decorator, &
      Kipp_use_decorator, &
      logg_logt_use_decorator, &
      logg_teff_use_decorator, &
      logl_r_use_decorator, &
      logl_teff_use_decorator, &
      logl_v_use_decorator, &
      l_r_use_decorator, &
      l_teff_use_decorator, &
      l_v_use_decorator, &
      mixing_use_decorator, &
      mode_prop_use_decorator, &
      network_use_decorator, &
      power_use_decorator, &
      production_use_decorator, &
      Profile_Panels1_use_decorator, &
      Profile_Panels2_use_decorator, &
      Profile_Panels3_use_decorator, &
      Profile_Panels4_use_decorator, &
      Profile_Panels5_use_decorator, &
      Profile_Panels6_use_decorator, &
      Profile_Panels7_use_decorator, &
      Profile_Panels8_use_decorator, &
      Profile_Panels9_use_decorator, &
      R_L_use_decorator, &
      R_Teff_use_decorator, &
      rti_use_decorator, &
      summary_burn_use_decorator, &
      summary_profile_use_decorator, &
      summary_history_use_decorator, &
      trho_use_decorator, &
      tmaxrho_use_decorator, &
      trho_profile_use_decorator
