! ***********************************************************************
!
!   Copyright (C) 2010-2022  Bill Paxton, Matthias Fabry
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

module pgbinary_ctrls_io

   use const_def
   use binary_private_def

   implicit none

   include "pgbinary_controls.inc"

   namelist /pgbinary/ &

      file_device, &
      file_extension, &
      file_digits, &
      pgbinary_interval, &
      pause, &
      pause_interval, &
      pgbinary_sleep, &
      clear_history, &

      file_white_on_black_flag, &
      win_white_on_black_flag, &
      pgbinary_show_model_number, &
      pgbinary_show_age, &
      pgbinary_show_age_in_seconds, &
      pgbinary_show_age_in_minutes, &
      pgbinary_show_age_in_hours, &
      pgbinary_show_age_in_days, &
      pgbinary_show_age_in_years, &
      pgbinary_show_log_age_in_years, &

      pgbinary_report_writing_files, &

      pgbinary_show_title, &
      pgbinary_title_scale, &
      pgbinary_title_disp, &
      pgbinary_title_coord, &
      pgbinary_title_fjust, &
      pgbinary_title_lw, &

      pgbinary_grid_show_title, &
      pgbinary_grid_title_scale, &
      pgbinary_grid_title_disp, &
      pgbinary_grid_title_coord, &
      pgbinary_grid_title_fjust, &
      pgbinary_grid_title_lw, &

      pgbinary_age_scale, &
      pgbinary_age_disp, &
      pgbinary_age_coord, &
      pgbinary_age_fjust, &
      pgbinary_age_lw, &

      pgbinary_model_scale, &
      pgbinary_model_disp, &
      pgbinary_model_coord, &
      pgbinary_model_fjust, &
      pgbinary_model_lw, &

      pgbinary_xaxis_label_scale, &
      pgbinary_left_yaxis_label_scale, &
      pgbinary_right_yaxis_label_scale, &
      pgbinary_xaxis_label_lw, &
      pgbinary_left_yaxis_label_lw, &
      pgbinary_right_yaxis_label_lw, &
      pgbinary_xaxis_label_disp, &
      pgbinary_left_yaxis_label_disp, &
      pgbinary_right_yaxis_label_disp, &
      pgbinary_num_scale, &
      pgbinary_lw, &
      pgbinary_model_lw, &
      pgbinary_box_lw, &


      Text_Summary1_win_flag, &
      Text_Summary1_file_flag, &
      Text_Summary1_file_interval, &
      Text_Summary1_file_dir, &
      Text_Summary1_file_prefix, &
      Text_Summary1_num_cols, Text_Summary1_num_rows, Text_Summary1_name, &
      Text_Summary1_win_width, &
      Text_Summary1_win_aspect_ratio, &
      Text_Summary1_file_width, &
      Text_Summary1_file_aspect_ratio, &
      Text_Summary1_title, Text_Summary1_xleft, Text_Summary1_xright, &
      Text_Summary1_ybot, Text_Summary1_ytop, Text_Summary1_txt_scale, &

      Text_Summary2_win_flag, &
      Text_Summary2_file_flag, &
      Text_Summary2_file_interval, &
      Text_Summary2_file_dir, &
      Text_Summary2_file_prefix, &
      Text_Summary2_num_cols, Text_Summary2_num_rows, Text_Summary2_name, &
      Text_Summary2_win_width, &
      Text_Summary2_win_aspect_ratio, &
      Text_Summary2_file_width, &
      Text_Summary2_file_aspect_ratio, &
      Text_Summary2_title, Text_Summary2_xleft, Text_Summary2_xright, &
      Text_Summary2_ybot, Text_Summary2_ytop, Text_Summary2_txt_scale, &

      Text_Summary3_win_flag, &
      Text_Summary3_file_flag, &
      Text_Summary3_file_interval, &
      Text_Summary3_file_dir, &
      Text_Summary3_file_prefix, &
      Text_Summary3_num_cols, Text_Summary3_num_rows, Text_Summary3_name, &
      Text_Summary3_win_width, &
      Text_Summary3_win_aspect_ratio, &
      Text_Summary3_file_width, &
      Text_Summary3_file_aspect_ratio, &
      Text_Summary3_title, Text_Summary3_xleft, Text_Summary3_xright, &
      Text_Summary3_ybot, Text_Summary3_ytop, Text_Summary3_txt_scale, &

      Text_Summary4_win_flag, &
      Text_Summary4_file_flag, &
      Text_Summary4_file_interval, &
      Text_Summary4_file_dir, &
      Text_Summary4_file_prefix, &
      Text_Summary4_num_cols, Text_Summary4_num_rows, Text_Summary4_name, &
      Text_Summary4_win_width, &
      Text_Summary4_win_aspect_ratio, &
      Text_Summary4_file_width, &
      Text_Summary4_file_aspect_ratio, &
      Text_Summary4_title, Text_Summary4_xleft, Text_Summary4_xright, &
      Text_Summary4_ybot, Text_Summary4_ytop, Text_Summary4_txt_scale, &

      Text_Summary5_win_flag, &
      Text_Summary5_file_flag, &
      Text_Summary5_file_interval, &
      Text_Summary5_file_dir, &
      Text_Summary5_file_prefix, &
      Text_Summary5_num_cols, Text_Summary5_num_rows, Text_Summary5_name, &
      Text_Summary5_win_width, &
      Text_Summary5_win_aspect_ratio, &
      Text_Summary5_file_width, &
      Text_Summary5_file_aspect_ratio, &
      Text_Summary5_title, Text_Summary5_xleft, Text_Summary5_xright, &
      Text_Summary5_ybot, Text_Summary5_ytop, Text_Summary5_txt_scale, &

      Text_Summary6_win_flag, &
      Text_Summary6_file_flag, &
      Text_Summary6_file_interval, &
      Text_Summary6_file_dir, &
      Text_Summary6_file_prefix, &
      Text_Summary6_num_cols, Text_Summary6_num_rows, Text_Summary6_name, &
      Text_Summary6_win_width, &
      Text_Summary6_win_aspect_ratio, &
      Text_Summary6_file_width, &
      Text_Summary6_file_aspect_ratio, &
      Text_Summary6_title, Text_Summary6_xleft, Text_Summary6_xright, &
      Text_Summary6_ybot, Text_Summary6_ytop, Text_Summary6_txt_scale, &

      Text_Summary7_win_flag, &
      Text_Summary7_file_flag, &
      Text_Summary7_file_interval, &
      Text_Summary7_file_dir, &
      Text_Summary7_file_prefix, &
      Text_Summary7_num_cols, Text_Summary7_num_rows, Text_Summary7_name, &
      Text_Summary7_win_width, &
      Text_Summary7_win_aspect_ratio, &
      Text_Summary7_file_width, &
      Text_Summary7_file_aspect_ratio, &
      Text_Summary7_title, Text_Summary7_xleft, Text_Summary7_xright, &
      Text_Summary7_ybot, Text_Summary7_ytop, Text_Summary7_txt_scale, &

      Text_Summary8_win_flag, &
      Text_Summary8_file_flag, &
      Text_Summary8_file_interval, &
      Text_Summary8_file_dir, &
      Text_Summary8_file_prefix, &
      Text_Summary8_num_cols, Text_Summary8_num_rows, Text_Summary8_name, &
      Text_Summary8_win_width, &
      Text_Summary8_win_aspect_ratio, &
      Text_Summary8_file_width, &
      Text_Summary8_file_aspect_ratio, &
      Text_Summary8_title, Text_Summary8_xleft, Text_Summary8_xright, &
      Text_Summary8_ybot, Text_Summary8_ytop, Text_Summary8_txt_scale, &

      Text_Summary9_win_flag, &
      Text_Summary9_file_flag, &
      Text_Summary9_file_interval, &
      Text_Summary9_file_dir, &
      Text_Summary9_file_prefix, &
      Text_Summary9_num_cols, Text_Summary9_num_rows, Text_Summary9_name, &
      Text_Summary9_win_width, &
      Text_Summary9_win_aspect_ratio, &
      Text_Summary9_file_width, &
      Text_Summary9_file_aspect_ratio, &
      Text_Summary9_title, Text_Summary9_xleft, Text_Summary9_xright, &
      Text_Summary9_ybot, Text_Summary9_ytop, Text_Summary9_txt_scale, &

      History_Track1_win_flag, &
      History_Track1_file_flag, &
      History_Track1_file_interval, &
      History_Track1_step_min, &
      History_Track1_step_max, &
      show_History_Track1_target_box, &
      History_Track1_n_sigma, &
      History_Track1_xtarget, &
      History_Track1_xsigma, &
      History_Track1_ytarget, &
      History_Track1_ysigma, &
      History_Track1_xname, &
      History_Track1_xaxis_label, &
      History_Track1_yname, &
      History_Track1_yaxis_label, &
      History_Track1_file_dir, &
      History_Track1_file_prefix, &
      show_History_Track1_annotation1, &
      show_History_Track1_annotation2, &
      show_History_Track1_annotation3, &
      History_Track1_fname, &
      History_Track1_reverse_xaxis, &
      History_Track1_reverse_yaxis, &
      History_Track1_log_xaxis, &
      History_Track1_log_yaxis, &
      History_Track1_xmin, &
      History_Track1_xmax, &
      History_Track1_ymin, &
      History_Track1_ymax, &
      History_Track1_xmargin, &
      History_Track1_ymargin, &
      History_Track1_dxmin, &
      History_Track1_dymin, &
      History_Track1_win_width, &
      History_Track1_win_aspect_ratio, &
      History_Track1_file_width, &
      History_Track1_file_aspect_ratio, &
      History_Track1_xleft, &
      History_Track1_xright, &
      History_Track1_ybot, &
      History_Track1_ytop, &
      History_Track1_txt_scale, &
      History_Track1_title, &
      History_Track1_use_decorator, &

      History_Track2_win_flag, &
      History_Track2_file_flag, &
      History_Track2_file_interval, &
      History_Track2_step_min, &
      History_Track2_step_max, &
      show_History_Track2_target_box, &
      History_Track2_n_sigma, &
      History_Track2_xtarget, &
      History_Track2_xsigma, &
      History_Track2_ytarget, &
      History_Track2_ysigma, &
      History_Track2_xname, &
      History_Track2_xaxis_label, &
      History_Track2_yname, &
      History_Track2_yaxis_label, &
      History_Track2_file_dir, &
      History_Track2_file_prefix, &
      show_History_Track2_annotation1, &
      show_History_Track2_annotation2, &
      show_History_Track2_annotation3, &
      History_Track2_fname, &
      History_Track2_reverse_xaxis, &
      History_Track2_reverse_yaxis, &
      History_Track2_log_xaxis, &
      History_Track2_log_yaxis, &
      History_Track2_xmin, &
      History_Track2_xmax, &
      History_Track2_ymin, &
      History_Track2_ymax, &
      History_Track2_xmargin, &
      History_Track2_ymargin, &
      History_Track2_dxmin, &
      History_Track2_dymin, &
      History_Track2_win_width, &
      History_Track2_win_aspect_ratio, &
      History_Track2_file_width, &
      History_Track2_file_aspect_ratio, &
      History_Track2_xleft, &
      History_Track2_xright, &
      History_Track2_ybot, &
      History_Track2_ytop, &
      History_Track2_txt_scale, &
      History_Track2_title, &
      History_Track2_use_decorator, &

      History_Track3_win_flag, &
      History_Track3_file_flag, &
      History_Track3_file_interval, &
      History_Track3_step_min, &
      History_Track3_step_max, &
      show_History_Track3_target_box, &
      History_Track3_n_sigma, &
      History_Track3_xtarget, &
      History_Track3_xsigma, &
      History_Track3_ytarget, &
      History_Track3_ysigma, &
      History_Track3_xname, &
      History_Track3_xaxis_label, &
      History_Track3_yname, &
      History_Track3_yaxis_label, &
      History_Track3_file_dir, &
      History_Track3_file_prefix, &
      show_History_Track3_annotation1, &
      show_History_Track3_annotation2, &
      show_History_Track3_annotation3, &
      History_Track3_fname, &
      History_Track3_reverse_xaxis, &
      History_Track3_reverse_yaxis, &
      History_Track3_log_xaxis, &
      History_Track3_log_yaxis, &
      History_Track3_xmin, &
      History_Track3_xmax, &
      History_Track3_ymin, &
      History_Track3_ymax, &
      History_Track3_xmargin, &
      History_Track3_ymargin, &
      History_Track3_dxmin, &
      History_Track3_dymin, &
      History_Track3_win_width, &
      History_Track3_win_aspect_ratio, &
      History_Track3_file_width, &
      History_Track3_file_aspect_ratio, &
      History_Track3_xleft, &
      History_Track3_xright, &
      History_Track3_ybot, &
      History_Track3_ytop, &
      History_Track3_txt_scale, &
      History_Track3_title, &
      History_Track3_use_decorator, &

      History_Track4_win_flag, &
      History_Track4_file_flag, &
      History_Track4_file_interval, &
      History_Track4_step_min, &
      History_Track4_step_max, &
      show_History_Track4_target_box, &
      History_Track4_n_sigma, &
      History_Track4_xtarget, &
      History_Track4_xsigma, &
      History_Track4_ytarget, &
      History_Track4_ysigma, &
      History_Track4_xname, &
      History_Track4_xaxis_label, &
      History_Track4_yname, &
      History_Track4_yaxis_label, &
      History_Track4_file_dir, &
      History_Track4_file_prefix, &
      show_History_Track4_annotation1, &
      show_History_Track4_annotation2, &
      show_History_Track4_annotation3, &
      History_Track4_fname, &
      History_Track4_reverse_xaxis, &
      History_Track4_reverse_yaxis, &
      History_Track4_log_xaxis, &
      History_Track4_log_yaxis, &
      History_Track4_xmin, &
      History_Track4_xmax, &
      History_Track4_ymin, &
      History_Track4_ymax, &
      History_Track4_xmargin, &
      History_Track4_ymargin, &
      History_Track4_dxmin, &
      History_Track4_dymin, &
      History_Track4_win_width, &
      History_Track4_win_aspect_ratio, &
      History_Track4_file_width, &
      History_Track4_file_aspect_ratio, &
      History_Track4_xleft, &
      History_Track4_xright, &
      History_Track4_ybot, &
      History_Track4_ytop, &
      History_Track4_txt_scale, &
      History_Track4_title, &
      History_Track4_use_decorator, &

      History_Track5_win_flag, &
      History_Track5_file_flag, &
      History_Track5_file_interval, &
      History_Track5_step_min, &
      History_Track5_step_max, &
      show_History_Track5_target_box, &
      History_Track5_n_sigma, &
      History_Track5_xtarget, &
      History_Track5_xsigma, &
      History_Track5_ytarget, &
      History_Track5_ysigma, &
      History_Track5_xname, &
      History_Track5_xaxis_label, &
      History_Track5_yname, &
      History_Track5_yaxis_label, &
      History_Track5_file_dir, &
      History_Track5_file_prefix, &
      show_History_Track5_annotation1, &
      show_History_Track5_annotation2, &
      show_History_Track5_annotation3, &
      History_Track5_fname, &
      History_Track5_reverse_xaxis, &
      History_Track5_reverse_yaxis, &
      History_Track5_log_xaxis, &
      History_Track5_log_yaxis, &
      History_Track5_xmin, &
      History_Track5_xmax, &
      History_Track5_ymin, &
      History_Track5_ymax, &
      History_Track5_xmargin, &
      History_Track5_ymargin, &
      History_Track5_dxmin, &
      History_Track5_dymin, &
      History_Track5_win_width, &
      History_Track5_win_aspect_ratio, &
      History_Track5_file_width, &
      History_Track5_file_aspect_ratio, &
      History_Track5_xleft, &
      History_Track5_xright, &
      History_Track5_ybot, &
      History_Track5_ytop, &
      History_Track5_txt_scale, &
      History_Track5_title, &
      History_Track5_use_decorator, &

      History_Track6_win_flag, &
      History_Track6_file_flag, &
      History_Track6_file_interval, &
      History_Track6_step_min, &
      History_Track6_step_max, &
      show_History_Track6_target_box, &
      History_Track6_n_sigma, &
      History_Track6_xtarget, &
      History_Track6_xsigma, &
      History_Track6_ytarget, &
      History_Track6_ysigma, &
      History_Track6_xname, &
      History_Track6_xaxis_label, &
      History_Track6_yname, &
      History_Track6_yaxis_label, &
      History_Track6_file_dir, &
      History_Track6_file_prefix, &
      show_History_Track6_annotation1, &
      show_History_Track6_annotation2, &
      show_History_Track6_annotation3, &
      History_Track6_fname, &
      History_Track6_reverse_xaxis, &
      History_Track6_reverse_yaxis, &
      History_Track6_log_xaxis, &
      History_Track6_log_yaxis, &
      History_Track6_xmin, &
      History_Track6_xmax, &
      History_Track6_ymin, &
      History_Track6_ymax, &
      History_Track6_xmargin, &
      History_Track6_ymargin, &
      History_Track6_dxmin, &
      History_Track6_dymin, &
      History_Track6_win_width, &
      History_Track6_win_aspect_ratio, &
      History_Track6_file_width, &
      History_Track6_file_aspect_ratio, &
      History_Track6_xleft, &
      History_Track6_xright, &
      History_Track6_ybot, &
      History_Track6_ytop, &
      History_Track6_txt_scale, &
      History_Track6_title, &
      History_Track6_use_decorator, &

      History_Track7_win_flag, &
      History_Track7_file_flag, &
      History_Track7_file_interval, &
      History_Track7_step_min, &
      History_Track7_step_max, &
      show_History_Track7_target_box, &
      History_Track7_n_sigma, &
      History_Track7_xtarget, &
      History_Track7_xsigma, &
      History_Track7_ytarget, &
      History_Track7_ysigma, &
      History_Track7_xname, &
      History_Track7_xaxis_label, &
      History_Track7_yname, &
      History_Track7_yaxis_label, &
      History_Track7_file_dir, &
      History_Track7_file_prefix, &
      show_History_Track7_annotation1, &
      show_History_Track7_annotation2, &
      show_History_Track7_annotation3, &
      History_Track7_fname, &
      History_Track7_reverse_xaxis, &
      History_Track7_reverse_yaxis, &
      History_Track7_log_xaxis, &
      History_Track7_log_yaxis, &
      History_Track7_xmin, &
      History_Track7_xmax, &
      History_Track7_ymin, &
      History_Track7_ymax, &
      History_Track7_xmargin, &
      History_Track7_ymargin, &
      History_Track7_dxmin, &
      History_Track7_dymin, &
      History_Track7_win_width, &
      History_Track7_win_aspect_ratio, &
      History_Track7_file_width, &
      History_Track7_file_aspect_ratio, &
      History_Track7_xleft, &
      History_Track7_xright, &
      History_Track7_ybot, &
      History_Track7_ytop, &
      History_Track7_txt_scale, &
      History_Track7_title, &
      History_Track7_use_decorator, &

      History_Track8_win_flag, &
      History_Track8_file_flag, &
      History_Track8_file_interval, &
      History_Track8_step_min, &
      History_Track8_step_max, &
      show_History_Track8_target_box, &
      History_Track8_n_sigma, &
      History_Track8_xtarget, &
      History_Track8_xsigma, &
      History_Track8_ytarget, &
      History_Track8_ysigma, &
      History_Track8_xname, &
      History_Track8_xaxis_label, &
      History_Track8_yname, &
      History_Track8_yaxis_label, &
      History_Track8_file_dir, &
      History_Track8_file_prefix, &
      show_History_Track8_annotation1, &
      show_History_Track8_annotation2, &
      show_History_Track8_annotation3, &
      History_Track8_fname, &
      History_Track8_reverse_xaxis, &
      History_Track8_reverse_yaxis, &
      History_Track8_log_xaxis, &
      History_Track8_log_yaxis, &
      History_Track8_xmin, &
      History_Track8_xmax, &
      History_Track8_ymin, &
      History_Track8_ymax, &
      History_Track8_xmargin, &
      History_Track8_ymargin, &
      History_Track8_dxmin, &
      History_Track8_dymin, &
      History_Track8_win_width, &
      History_Track8_win_aspect_ratio, &
      History_Track8_file_width, &
      History_Track8_file_aspect_ratio, &
      History_Track8_xleft, &
      History_Track8_xright, &
      History_Track8_ybot, &
      History_Track8_ytop, &
      History_Track8_txt_scale, &
      History_Track8_title, &
      History_Track8_use_decorator, &

      History_Track9_win_flag, &
      History_Track9_file_flag, &
      History_Track9_file_interval, &
      History_Track9_step_min, &
      History_Track9_step_max, &
      show_History_Track9_target_box, &
      History_Track9_n_sigma, &
      History_Track9_xtarget, &
      History_Track9_xsigma, &
      History_Track9_ytarget, &
      History_Track9_ysigma, &
      History_Track9_xname, &
      History_Track9_xaxis_label, &
      History_Track9_yname, &
      History_Track9_yaxis_label, &
      History_Track9_file_dir, &
      History_Track9_file_prefix, &
      show_History_Track9_annotation1, &
      show_History_Track9_annotation2, &
      show_History_Track9_annotation3, &
      History_Track9_fname, &
      History_Track9_reverse_xaxis, &
      History_Track9_reverse_yaxis, &
      History_Track9_log_xaxis, &
      History_Track9_log_yaxis, &
      History_Track9_xmin, &
      History_Track9_xmax, &
      History_Track9_ymin, &
      History_Track9_ymax, &
      History_Track9_xmargin, &
      History_Track9_ymargin, &
      History_Track9_dxmin, &
      History_Track9_dymin, &
      History_Track9_win_width, &
      History_Track9_win_aspect_ratio, &
      History_Track9_file_width, &
      History_Track9_file_aspect_ratio, &
      History_Track9_xleft, &
      History_Track9_xright, &
      History_Track9_ybot, &
      History_Track9_ytop, &
      History_Track9_txt_scale, &
      History_Track9_title, &
      History_Track9_use_decorator, &

      History_Panels1_win_flag, &
      History_Panels1_win_width, &
      History_Panels1_win_aspect_ratio, &
      History_Panels1_xleft, &
      History_Panels1_xright, &
      History_Panels1_ybot, &
      History_Panels1_ytop, &
      History_Panels1_txt_scale, &
      History_Panels1_title, &
      History_Panels1_xmax, &
      History_Panels1_xmin, &
      History_Panels1_dxmin, &
      History_Panels1_max_width, &
      History_Panels1_num_panels, &
      History_Panels1_xaxis_name, &
      History_Panels1_yaxis_name, &
      History_Panels1_xaxis_reversed, &
      History_Panels1_yaxis_reversed, &
      History_Panels1_yaxis_log, &
      History_Panels1_ymin, &
      History_Panels1_ymax, &
      History_Panels1_dymin, &
      History_Panels1_other_yaxis_name, &
      History_Panels1_other_yaxis_reversed, &
      History_Panels1_xaxis_log, &
      History_Panels1_other_yaxis_log, &
      History_Panels1_other_ymin, &
      History_Panels1_other_ymax, &
      History_Panels1_other_dymin, &
      History_Panels1_points_name, &
      History_Panels1_file_flag, &
      History_Panels1_file_dir, &
      History_Panels1_file_prefix, &
      History_Panels1_file_interval, &
      History_Panels1_file_width, &
      History_Panels1_file_aspect_ratio, &
      History_Panels1_xmargin, &
      History_Panels1_ymargin, &
      History_Panels1_other_ymargin, &
      History_Panels1_use_decorator, &

      History_Panels2_win_flag, &
      History_Panels2_win_width, &
      History_Panels2_win_aspect_ratio, &
      History_Panels2_xleft, &
      History_Panels2_xright, &
      History_Panels2_ybot, &
      History_Panels2_ytop, &
      History_Panels2_txt_scale, &
      History_Panels2_title, &
      History_Panels2_xmax, &
      History_Panels2_xmin, &
      History_Panels2_dxmin, &
      History_Panels2_max_width, &
      History_Panels2_num_panels, &
      History_Panels2_xaxis_name, &
      History_Panels2_yaxis_name, &
      History_Panels2_xaxis_reversed, &
      History_Panels2_yaxis_reversed, &
      History_Panels2_yaxis_log, &
      History_Panels2_ymin, &
      History_Panels2_ymax, &
      History_Panels2_dymin, &
      History_Panels2_other_yaxis_name, &
      History_Panels2_other_yaxis_reversed, &
      History_Panels2_xaxis_log, &
      History_Panels2_other_yaxis_log, &
      History_Panels2_other_ymin, &
      History_Panels2_other_ymax, &
      History_Panels2_other_dymin, &
      History_Panels2_points_name, &
      History_Panels2_file_flag, &
      History_Panels2_file_dir, &
      History_Panels2_file_prefix, &
      History_Panels2_file_interval, &
      History_Panels2_file_width, &
      History_Panels2_file_aspect_ratio, &
      History_Panels2_xmargin, &
      History_Panels2_ymargin, &
      History_Panels2_other_ymargin, &
      History_Panels2_use_decorator, &

      History_Panels3_win_flag, &
      History_Panels3_win_width, &
      History_Panels3_win_aspect_ratio, &
      History_Panels3_xleft, &
      History_Panels3_xright, &
      History_Panels3_ybot, &
      History_Panels3_ytop, &
      History_Panels3_txt_scale, &
      History_Panels3_title, &
      History_Panels3_xmax, &
      History_Panels3_xmin, &
      History_Panels3_dxmin, &
      History_Panels3_max_width, &
      History_Panels3_num_panels, &
      History_Panels3_xaxis_name, &
      History_Panels3_yaxis_name, &
      History_Panels3_xaxis_reversed, &
      History_Panels3_yaxis_reversed, &
      History_Panels3_yaxis_log, &
      History_Panels3_ymin, &
      History_Panels3_ymax, &
      History_Panels3_dymin, &
      History_Panels3_other_yaxis_name, &
      History_Panels3_other_yaxis_reversed, &
      History_Panels3_xaxis_log, &
      History_Panels3_other_yaxis_log, &
      History_Panels3_other_ymin, &
      History_Panels3_other_ymax, &
      History_Panels3_other_dymin, &
      History_Panels3_points_name, &
      History_Panels3_file_flag, &
      History_Panels3_file_dir, &
      History_Panels3_file_prefix, &
      History_Panels3_file_interval, &
      History_Panels3_file_width, &
      History_Panels3_file_aspect_ratio, &
      History_Panels3_xmargin, &
      History_Panels3_ymargin, &
      History_Panels3_other_ymargin, &
      History_Panels3_use_decorator, &

      History_Panels4_win_flag, &
      History_Panels4_win_width, &
      History_Panels4_win_aspect_ratio, &
      History_Panels4_xleft, &
      History_Panels4_xright, &
      History_Panels4_ybot, &
      History_Panels4_ytop, &
      History_Panels4_txt_scale, &
      History_Panels4_title, &
      History_Panels4_xmax, &
      History_Panels4_xmin, &
      History_Panels4_dxmin, &
      History_Panels4_max_width, &
      History_Panels4_num_panels, &
      History_Panels4_xaxis_name, &
      History_Panels4_yaxis_name, &
      History_Panels4_xaxis_reversed, &
      History_Panels4_yaxis_reversed, &
      History_Panels4_yaxis_log, &
      History_Panels4_ymin, &
      History_Panels4_ymax, &
      History_Panels4_dymin, &
      History_Panels4_other_yaxis_name, &
      History_Panels4_other_yaxis_reversed, &
      History_Panels4_xaxis_log, &
      History_Panels4_other_yaxis_log, &
      History_Panels4_other_ymin, &
      History_Panels4_other_ymax, &
      History_Panels4_other_dymin, &
      History_Panels4_points_name, &
      History_Panels4_file_flag, &
      History_Panels4_file_dir, &
      History_Panels4_file_prefix, &
      History_Panels4_file_interval, &
      History_Panels4_file_width, &
      History_Panels4_file_aspect_ratio, &
      History_Panels4_xmargin, &
      History_Panels4_ymargin, &
      History_Panels4_other_ymargin, &
      History_Panels4_use_decorator, &

      History_Panels5_win_flag, &
      History_Panels5_win_width, &
      History_Panels5_win_aspect_ratio, &
      History_Panels5_xleft, &
      History_Panels5_xright, &
      History_Panels5_ybot, &
      History_Panels5_ytop, &
      History_Panels5_txt_scale, &
      History_Panels5_title, &
      History_Panels5_xmax, &
      History_Panels5_xmin, &
      History_Panels5_dxmin, &
      History_Panels5_max_width, &
      History_Panels5_num_panels, &
      History_Panels5_xaxis_name, &
      History_Panels5_yaxis_name, &
      History_Panels5_xaxis_reversed, &
      History_Panels5_yaxis_reversed, &
      History_Panels5_yaxis_log, &
      History_Panels5_ymin, &
      History_Panels5_ymax, &
      History_Panels5_dymin, &
      History_Panels5_other_yaxis_name, &
      History_Panels5_other_yaxis_reversed, &
      History_Panels5_xaxis_log, &
      History_Panels5_other_yaxis_log, &
      History_Panels5_other_ymin, &
      History_Panels5_other_ymax, &
      History_Panels5_other_dymin, &
      History_Panels5_points_name, &
      History_Panels5_file_flag, &
      History_Panels5_file_dir, &
      History_Panels5_file_prefix, &
      History_Panels5_file_interval, &
      History_Panels5_file_width, &
      History_Panels5_file_aspect_ratio, &
      History_Panels5_xmargin, &
      History_Panels5_ymargin, &
      History_Panels5_other_ymargin, &
      History_Panels5_use_decorator, &

      History_Panels6_win_flag, &
      History_Panels6_win_width, &
      History_Panels6_win_aspect_ratio, &
      History_Panels6_xleft, &
      History_Panels6_xright, &
      History_Panels6_ybot, &
      History_Panels6_ytop, &
      History_Panels6_txt_scale, &
      History_Panels6_title, &
      History_Panels6_xmax, &
      History_Panels6_xmin, &
      History_Panels6_dxmin, &
      History_Panels6_max_width, &
      History_Panels6_num_panels, &
      History_Panels6_xaxis_name, &
      History_Panels6_yaxis_name, &
      History_Panels6_xaxis_reversed, &
      History_Panels6_yaxis_reversed, &
      History_Panels6_yaxis_log, &
      History_Panels6_ymin, &
      History_Panels6_ymax, &
      History_Panels6_dymin, &
      History_Panels6_other_yaxis_name, &
      History_Panels6_other_yaxis_reversed, &
      History_Panels6_xaxis_log, &
      History_Panels6_other_yaxis_log, &
      History_Panels6_other_ymin, &
      History_Panels6_other_ymax, &
      History_Panels6_other_dymin, &
      History_Panels6_points_name, &
      History_Panels6_file_flag, &
      History_Panels6_file_dir, &
      History_Panels6_file_prefix, &
      History_Panels6_file_interval, &
      History_Panels6_file_width, &
      History_Panels6_file_aspect_ratio, &
      History_Panels6_xmargin, &
      History_Panels6_ymargin, &
      History_Panels6_other_ymargin, &
      History_Panels6_use_decorator, &

      History_Panels7_win_flag, &
      History_Panels7_win_width, &
      History_Panels7_win_aspect_ratio, &
      History_Panels7_xleft, &
      History_Panels7_xright, &
      History_Panels7_ybot, &
      History_Panels7_ytop, &
      History_Panels7_txt_scale, &
      History_Panels7_title, &
      History_Panels7_xmax, &
      History_Panels7_xmin, &
      History_Panels7_dxmin, &
      History_Panels7_max_width, &
      History_Panels7_num_panels, &
      History_Panels7_xaxis_name, &
      History_Panels7_yaxis_name, &
      History_Panels7_xaxis_reversed, &
      History_Panels7_yaxis_reversed, &
      History_Panels7_yaxis_log, &
      History_Panels7_ymin, &
      History_Panels7_ymax, &
      History_Panels7_dymin, &
      History_Panels7_other_yaxis_name, &
      History_Panels7_other_yaxis_reversed, &
      History_Panels7_xaxis_log, &
      History_Panels7_other_yaxis_log, &
      History_Panels7_other_ymin, &
      History_Panels7_other_ymax, &
      History_Panels7_other_dymin, &
      History_Panels7_points_name, &
      History_Panels7_file_flag, &
      History_Panels7_file_dir, &
      History_Panels7_file_prefix, &
      History_Panels7_file_interval, &
      History_Panels7_file_width, &
      History_Panels7_file_aspect_ratio, &
      History_Panels7_xmargin, &
      History_Panels7_ymargin, &
      History_Panels7_other_ymargin, &
      History_Panels7_use_decorator, &

      History_Panels8_win_flag, &
      History_Panels8_win_width, &
      History_Panels8_win_aspect_ratio, &
      History_Panels8_xleft, &
      History_Panels8_xright, &
      History_Panels8_ybot, &
      History_Panels8_ytop, &
      History_Panels8_txt_scale, &
      History_Panels8_title, &
      History_Panels8_xmax, &
      History_Panels8_xmin, &
      History_Panels8_dxmin, &
      History_Panels8_max_width, &
      History_Panels8_num_panels, &
      History_Panels8_xaxis_name, &
      History_Panels8_yaxis_name, &
      History_Panels8_xaxis_reversed, &
      History_Panels8_yaxis_reversed, &
      History_Panels8_yaxis_log, &
      History_Panels8_ymin, &
      History_Panels8_ymax, &
      History_Panels8_dymin, &
      History_Panels8_other_yaxis_name, &
      History_Panels8_other_yaxis_reversed, &
      History_Panels8_xaxis_log, &
      History_Panels8_other_yaxis_log, &
      History_Panels8_other_ymin, &
      History_Panels8_other_ymax, &
      History_Panels8_other_dymin, &
      History_Panels8_points_name, &
      History_Panels8_file_flag, &
      History_Panels8_file_dir, &
      History_Panels8_file_prefix, &
      History_Panels8_file_interval, &
      History_Panels8_file_width, &
      History_Panels8_file_aspect_ratio, &
      History_Panels8_xmargin, &
      History_Panels8_ymargin, &
      History_Panels8_other_ymargin, &
      History_Panels8_use_decorator, &

      History_Panels9_win_flag, &
      History_Panels9_win_width, &
      History_Panels9_win_aspect_ratio, &
      History_Panels9_xleft, &
      History_Panels9_xright, &
      History_Panels9_ybot, &
      History_Panels9_ytop, &
      History_Panels9_txt_scale, &
      History_Panels9_title, &
      History_Panels9_xmax, &
      History_Panels9_xmin, &
      History_Panels9_dxmin, &
      History_Panels9_max_width, &
      History_Panels9_num_panels, &
      History_Panels9_xaxis_name, &
      History_Panels9_yaxis_name, &
      History_Panels9_xaxis_reversed, &
      History_Panels9_yaxis_reversed, &
      History_Panels9_yaxis_log, &
      History_Panels9_ymin, &
      History_Panels9_ymax, &
      History_Panels9_dymin, &
      History_Panels9_other_yaxis_name, &
      History_Panels9_other_yaxis_reversed, &
      History_Panels9_xaxis_log, &
      History_Panels9_other_yaxis_log, &
      History_Panels9_other_ymin, &
      History_Panels9_other_ymax, &
      History_Panels9_other_dymin, &
      History_Panels9_points_name, &
      History_Panels9_file_flag, &
      History_Panels9_file_dir, &
      History_Panels9_file_prefix, &
      History_Panels9_file_interval, &
      History_Panels9_file_width, &
      History_Panels9_file_aspect_ratio, &
      History_Panels9_xmargin, &
      History_Panels9_ymargin, &
      History_Panels9_other_ymargin, &
      History_Panels9_use_decorator, &

      History_Panel_points_error_bars, &
      History_Panel_points_interval, &
      History_Panel_points_marker, &
      History_Panel_points_ci, &
      History_Panel_points_lw, &
      History_Panel_points_ch, &

      Summary_History_win_flag, &
      Summary_History_file_flag, &
      Summary_History_file_interval, &
      Summary_History_file_dir, &
      Summary_History_file_prefix, &
      Summary_History_scaled_value, &
      Summary_History_xmin, &
      Summary_History_xmax, &
      Summary_History_max_width, &
      Summary_History_win_width, &
      Summary_History_win_aspect_ratio, &
      Summary_History_file_width, &
      Summary_History_file_aspect_ratio, &
      Summary_History_xleft, &
      Summary_History_xright, &
      Summary_History_ybot, &
      Summary_History_ytop, &
      Summary_History_txt_scale, &
      Summary_History_title, &
      Summary_History_name, &
      Summary_History_legend, &
      Summary_History_num_lines, &
      Summary_History_use_decorator, &

      Grid1_win_flag, &
      Grid1_win_width, &
      Grid1_win_aspect_ratio, &
      Grid1_xleft, &
      Grid1_xright, &
      Grid1_ybot, &
      Grid1_ytop, &
      Grid1_title, &
      Grid1_txt_scale_factor, &
      Grid1_num_cols, &
      Grid1_num_rows, &
      Grid1_num_plots, &
      Grid1_plot_name, &
      Grid1_plot_row, &
      Grid1_plot_rowspan, &
      Grid1_plot_col, &
      Grid1_plot_colspan, &
      Grid1_plot_pad_left, &
      Grid1_plot_pad_right, &
      Grid1_plot_pad_top, &
      Grid1_plot_pad_bot, &
      Grid1_file_flag, &
      Grid1_file_dir, &
      Grid1_file_prefix, &
      Grid1_file_interval, &
      Grid1_file_width, &
      Grid1_file_aspect_ratio, &

      Grid2_win_flag, &
      Grid2_win_width, &
      Grid2_win_aspect_ratio, &
      Grid2_xleft, &
      Grid2_xright, &
      Grid2_ybot, &
      Grid2_ytop, &
      Grid2_title, &
      Grid2_txt_scale_factor, &
      Grid2_num_cols, &
      Grid2_num_rows, &
      Grid2_num_plots, &
      Grid2_plot_name, &
      Grid2_plot_row, &
      Grid2_plot_rowspan, &
      Grid2_plot_col, &
      Grid2_plot_colspan, &
      Grid2_plot_pad_left, &
      Grid2_plot_pad_right, &
      Grid2_plot_pad_top, &
      Grid2_plot_pad_bot, &
      Grid2_file_flag, &
      Grid2_file_dir, &
      Grid2_file_prefix, &
      Grid2_file_interval, &
      Grid2_file_width, &
      Grid2_file_aspect_ratio, &

      Grid3_win_flag, &
      Grid3_win_width, &
      Grid3_win_aspect_ratio, &
      Grid3_xleft, &
      Grid3_xright, &
      Grid3_ybot, &
      Grid3_ytop, &
      Grid3_title, &
      Grid3_txt_scale_factor, &
      Grid3_num_cols, &
      Grid3_num_rows, &
      Grid3_num_plots, &
      Grid3_plot_name, &
      Grid3_plot_row, &
      Grid3_plot_rowspan, &
      Grid3_plot_col, &
      Grid3_plot_colspan, &
      Grid3_plot_pad_left, &
      Grid3_plot_pad_right, &
      Grid3_plot_pad_top, &
      Grid3_plot_pad_bot, &
      Grid3_file_flag, &
      Grid3_file_dir, &
      Grid3_file_prefix, &
      Grid3_file_interval, &
      Grid3_file_width, &
      Grid3_file_aspect_ratio, &

      Grid4_win_flag, &
      Grid4_win_width, &
      Grid4_win_aspect_ratio, &
      Grid4_xleft, &
      Grid4_xright, &
      Grid4_ybot, &
      Grid4_ytop, &
      Grid4_title, &
      Grid4_txt_scale_factor, &
      Grid4_num_cols, &
      Grid4_num_rows, &
      Grid4_num_plots, &
      Grid4_plot_name, &
      Grid4_plot_row, &
      Grid4_plot_rowspan, &
      Grid4_plot_col, &
      Grid4_plot_colspan, &
      Grid4_plot_pad_left, &
      Grid4_plot_pad_right, &
      Grid4_plot_pad_top, &
      Grid4_plot_pad_bot, &
      Grid4_file_flag, &
      Grid4_file_dir, &
      Grid4_file_prefix, &
      Grid4_file_interval, &
      Grid4_file_width, &
      Grid4_file_aspect_ratio, &

      Grid5_win_flag, &
      Grid5_win_width, &
      Grid5_win_aspect_ratio, &
      Grid5_xleft, &
      Grid5_xright, &
      Grid5_ybot, &
      Grid5_ytop, &
      Grid5_title, &
      Grid5_txt_scale_factor, &
      Grid5_num_cols, &
      Grid5_num_rows, &
      Grid5_num_plots, &
      Grid5_plot_name, &
      Grid5_plot_row, &
      Grid5_plot_rowspan, &
      Grid5_plot_col, &
      Grid5_plot_colspan, &
      Grid5_plot_pad_left, &
      Grid5_plot_pad_right, &
      Grid5_plot_pad_top, &
      Grid5_plot_pad_bot, &
      Grid5_file_flag, &
      Grid5_file_dir, &
      Grid5_file_prefix, &
      Grid5_file_interval, &
      Grid5_file_width, &
      Grid5_file_aspect_ratio, &

      Grid6_win_flag, &
      Grid6_win_width, &
      Grid6_win_aspect_ratio, &
      Grid6_xleft, &
      Grid6_xright, &
      Grid6_ybot, &
      Grid6_ytop, &
      Grid6_title, &
      Grid6_txt_scale_factor, &
      Grid6_num_cols, &
      Grid6_num_rows, &
      Grid6_num_plots, &
      Grid6_plot_name, &
      Grid6_plot_row, &
      Grid6_plot_rowspan, &
      Grid6_plot_col, &
      Grid6_plot_colspan, &
      Grid6_plot_pad_left, &
      Grid6_plot_pad_right, &
      Grid6_plot_pad_top, &
      Grid6_plot_pad_bot, &
      Grid6_file_flag, &
      Grid6_file_dir, &
      Grid6_file_prefix, &
      Grid6_file_interval, &
      Grid6_file_width, &
      Grid6_file_aspect_ratio, &

      Grid7_win_flag, &
      Grid7_win_width, &
      Grid7_win_aspect_ratio, &
      Grid7_xleft, &
      Grid7_xright, &
      Grid7_ybot, &
      Grid7_ytop, &
      Grid7_title, &
      Grid7_txt_scale_factor, &
      Grid7_num_cols, &
      Grid7_num_rows, &
      Grid7_num_plots, &
      Grid7_plot_name, &
      Grid7_plot_row, &
      Grid7_plot_rowspan, &
      Grid7_plot_col, &
      Grid7_plot_colspan, &
      Grid7_plot_pad_left, &
      Grid7_plot_pad_right, &
      Grid7_plot_pad_top, &
      Grid7_plot_pad_bot, &
      Grid7_file_flag, &
      Grid7_file_dir, &
      Grid7_file_prefix, &
      Grid7_file_interval, &
      Grid7_file_width, &
      Grid7_file_aspect_ratio, &

      Grid8_win_flag, &
      Grid8_win_width, &
      Grid8_win_aspect_ratio, &
      Grid8_xleft, &
      Grid8_xright, &
      Grid8_ybot, &
      Grid8_ytop, &
      Grid8_title, &
      Grid8_txt_scale_factor, &
      Grid8_num_cols, &
      Grid8_num_rows, &
      Grid8_num_plots, &
      Grid8_plot_name, &
      Grid8_plot_row, &
      Grid8_plot_rowspan, &
      Grid8_plot_col, &
      Grid8_plot_colspan, &
      Grid8_plot_pad_left, &
      Grid8_plot_pad_right, &
      Grid8_plot_pad_top, &
      Grid8_plot_pad_bot, &
      Grid8_file_flag, &
      Grid8_file_dir, &
      Grid8_file_prefix, &
      Grid8_file_interval, &
      Grid8_file_width, &
      Grid8_file_aspect_ratio, &

      Grid9_win_flag, &
      Grid9_win_width, &
      Grid9_win_aspect_ratio, &
      Grid9_xleft, &
      Grid9_xright, &
      Grid9_ybot, &
      Grid9_ytop, &
      Grid9_title, &
      Grid9_txt_scale_factor, &
      Grid9_num_cols, &
      Grid9_num_rows, &
      Grid9_num_plots, &
      Grid9_plot_name, &
      Grid9_plot_row, &
      Grid9_plot_rowspan, &
      Grid9_plot_col, &
      Grid9_plot_colspan, &
      Grid9_plot_pad_left, &
      Grid9_plot_pad_right, &
      Grid9_plot_pad_top, &
      Grid9_plot_pad_bot, &
      Grid9_file_flag, &
      Grid9_file_dir, &
      Grid9_file_prefix, &
      Grid9_file_interval, &
      Grid9_file_width, &
      Grid9_file_aspect_ratio, &

      Star1_win_flag, &
      Star1_file_flag, &
      Star1_file_interval, &
      Star1_file_dir, &
      Star1_file_prefix, &
      Star1_win_width, &
      Star1_win_aspect_ratio, &
      Star1_xleft, &
      Star1_xright, &
      Star1_ybot, &
      Star1_ytop, &
      Star1_file_width, &
      Star1_file_aspect_ratio, &
      Star1_txt_scale_factor, &
      Star1_title, &
      Star1_plot_name, &
      do_star1_box, &
      star1_box_pad_left, &
      star1_box_pad_right, &
      star1_box_pad_bot, &
      star1_box_pad_top, &

      Star2_win_flag, &
      Star2_file_flag, &
      Star2_file_interval, &
      Star2_file_dir, &
      Star2_file_prefix, &
      Star2_win_width, &
      Star2_win_aspect_ratio, &
      Star2_xleft, &
      Star2_xright, &
      Star2_ybot, &
      Star2_ytop, &
      Star2_file_width, &
      Star2_file_aspect_ratio, &
      Star2_txt_scale_factor, &
      Star2_title, &
      Star2_plot_name, &
      do_star2_box, &
      star2_box_pad_left, &
      star2_box_pad_right, &
      star2_box_pad_bot, &
      star2_box_pad_top, &

      annotation1_ci, &
      annotation1_ch, &
      annotation1_lw, &
      annotation1_cf, &
      annotation1_text, &
      annotation1_side, &
      annotation1_disp, &
      annotation1_coord, &
      annotation1_fjust, &

      annotation2_ci, &
      annotation2_ch, &
      annotation2_lw, &
      annotation2_cf, &
      annotation2_text, &
      annotation2_side, &
      annotation2_disp, &
      annotation2_coord, &
      annotation2_fjust, &

      annotation3_ci, &
      annotation3_ch, &
      annotation3_lw, &
      annotation3_cf, &
      annotation3_text, &
      annotation3_side, &
      annotation3_disp, &
      annotation3_coord, &
      annotation3_fjust, &

      read_extra_pgbinary_inlist1, &
      extra_pgbinary_inlist1_name, &

      read_extra_pgbinary_inlist2, &
      extra_pgbinary_inlist2_name, &

      read_extra_pgbinary_inlist3, &
      extra_pgbinary_inlist3_name, &

      read_extra_pgbinary_inlist4, &
      extra_pgbinary_inlist4_name, &

      read_extra_pgbinary_inlist5, &
      extra_pgbinary_inlist5_name


contains


   subroutine read_pgbinary(b, filename, ierr)
      use binary_private_def
      use utils_lib
      type (binary_info), pointer :: b
      character(*), intent(in) :: filename
      integer, intent(out) :: ierr
!      character (len = strlen) :: pgbinary_namelist_name
!      pgbinary_namelist_name = ''
      ierr = 0
      call set_default_pgbinary_controls
      call read_pgbinary_file(b, filename, 1, ierr)
   end subroutine read_pgbinary


   recursive subroutine read_pgbinary_file(b, filename, level, ierr)
      use binary_private_def
      use utils_lib
      character(*), intent(in) :: filename
      type (binary_info), pointer :: b
      integer, intent(in) :: level
      integer, intent(out) :: ierr
      logical :: read_extra1, read_extra2, read_extra3, read_extra4, read_extra5
      character (len = strlen) :: message, extra1, extra2, extra3, extra4, extra5
      integer :: unit

      ierr = 0

      if (level >= 10) then
         write(*, *) 'ERROR: too many levels of nested extra pgbinary inlist files'
         ierr = -1
         return
      end if
      if (len_trim(filename) > 0) then
         open(newunit = unit, file = trim(filename), action = 'read', delim = 'quote', status = 'old', iostat = ierr)
         if (ierr /= 0) then
            write(*, *) 'Failed to open pgbinary namelist file ', trim(filename)
            return
         end if
         read(unit, nml = pgbinary, iostat = ierr)
         close(unit)
         if (ierr /= 0) then
            write(*, *)
            write(*, *)
            write(*, '(a)') &
               'Failed while trying to read pgbinary namelist file: ' // trim(filename)
            write(*, '(a)') &
               'Perhaps the following runtime error message will help you find the problem.'
            write(*, *)
            open(newunit = unit, file = trim(filename), action = 'read', delim = 'quote', status = 'old', iostat = ierr)
            read(unit, nml = pgbinary)
            close(unit)
            return
         end if
      end if

      call store_pgbinary_controls(b, ierr)

      ! recursive calls to read other inlists

      read_extra1 = read_extra_pgbinary_inlist1
      read_extra_pgbinary_inlist1 = .false.
      extra1 = extra_pgbinary_inlist1_name
      extra_pgbinary_inlist1_name = 'undefined'

      read_extra2 = read_extra_pgbinary_inlist2
      read_extra_pgbinary_inlist2 = .false.
      extra2 = extra_pgbinary_inlist2_name
      extra_pgbinary_inlist2_name = 'undefined'

      read_extra3 = read_extra_pgbinary_inlist3
      read_extra_pgbinary_inlist3 = .false.
      extra3 = extra_pgbinary_inlist3_name
      extra_pgbinary_inlist3_name = 'undefined'

      read_extra4 = read_extra_pgbinary_inlist4
      read_extra_pgbinary_inlist4 = .false.
      extra4 = extra_pgbinary_inlist4_name
      extra_pgbinary_inlist4_name = 'undefined'

      read_extra5 = read_extra_pgbinary_inlist5
      read_extra_pgbinary_inlist5 = .false.
      extra5 = extra_pgbinary_inlist5_name
      extra_pgbinary_inlist5_name = 'undefined'

      if (read_extra1) then
         call read_pgbinary_file(b, extra1, level + 1, ierr)
         if (ierr /= 0) return
      end if

      if (read_extra2) then
         call read_pgbinary_file(b, extra2, level + 1, ierr)
         if (ierr /= 0) return
      end if

      if (read_extra3) then
         call read_pgbinary_file(b, extra3, level + 1, ierr)
         if (ierr /= 0) return
      end if

      if (read_extra4) then
         call read_pgbinary_file(b, extra4, level + 1, ierr)
         if (ierr /= 0) return
      end if

      if (read_extra5) then
         call read_pgbinary_file(b, extra5, level + 1, ierr)
         if (ierr /= 0) return
      end if

   end subroutine read_pgbinary_file


   subroutine store_pgbinary_controls(b, ierr)
      use binary_private_def
      type (binary_info), pointer :: b
      integer, intent(out) :: ierr

      ierr = 0

      b% file_device = file_device
      b% file_extension = file_extension
      b% file_digits = file_digits
      b% pgbinary_interval = pgbinary_interval
      b% pause = pause
      b% pause_interval = pause_interval
      b% pgbinary_sleep = pgbinary_sleep
      b% clear_history = clear_history

      b% file_white_on_black_flag = file_white_on_black_flag
      b% win_white_on_black_flag = win_white_on_black_flag
      b% pgbinary_show_model_number = pgbinary_show_model_number
      b% pgbinary_show_age = pgbinary_show_age
      b% pgbinary_show_age_in_seconds = pgbinary_show_age_in_seconds
      b% pgbinary_show_age_in_minutes = pgbinary_show_age_in_minutes
      b% pgbinary_show_age_in_hours = pgbinary_show_age_in_hours
      b% pgbinary_show_age_in_days = pgbinary_show_age_in_days
      b% pgbinary_show_age_in_years = pgbinary_show_age_in_years
      b% pgbinary_show_log_age_in_years = pgbinary_show_log_age_in_years

      b% pgbinary_report_writing_files = pgbinary_report_writing_files

      b% pgbinary_show_title = pgbinary_show_title
      b% pgbinary_title_scale = pgbinary_title_scale
      b% pgbinary_title_disp = pgbinary_title_disp
      b% pgbinary_title_coord = pgbinary_title_coord
      b% pgbinary_title_fjust = pgbinary_title_fjust
      b% pgbinary_title_lw = pgbinary_title_lw

      b% pgbinary_grid_show_title = pgbinary_grid_show_title
      b% pgbinary_grid_title_scale = pgbinary_grid_title_scale
      b% pgbinary_grid_title_disp = pgbinary_grid_title_disp
      b% pgbinary_grid_title_coord = pgbinary_grid_title_coord
      b% pgbinary_grid_title_fjust = pgbinary_grid_title_fjust
      b% pgbinary_grid_title_lw = pgbinary_grid_title_lw

      b% pgbinary_age_scale = pgbinary_age_scale
      b% pgbinary_age_disp = pgbinary_age_disp
      b% pgbinary_age_coord = pgbinary_age_coord
      b% pgbinary_age_fjust = pgbinary_age_fjust
      b% pgbinary_age_lw = pgbinary_age_lw

      b% pgbinary_model_scale = pgbinary_model_scale
      b% pgbinary_model_disp = pgbinary_model_disp
      b% pgbinary_model_coord = pgbinary_model_coord
      b% pgbinary_model_fjust = pgbinary_model_fjust
      b% pgbinary_model_lw = pgbinary_model_lw

      b% pgbinary_xaxis_label_scale = pgbinary_xaxis_label_scale
      b% pgbinary_left_yaxis_label_scale = pgbinary_left_yaxis_label_scale
      b% pgbinary_right_yaxis_label_scale = pgbinary_right_yaxis_label_scale
      b% pgbinary_xaxis_label_lw = pgbinary_xaxis_label_lw
      b% pgbinary_left_yaxis_label_lw = pgbinary_left_yaxis_label_lw
      b% pgbinary_right_yaxis_label_lw = pgbinary_right_yaxis_label_lw
      b% pgbinary_xaxis_label_disp = pgbinary_xaxis_label_disp
      b% pgbinary_left_yaxis_label_disp = pgbinary_left_yaxis_label_disp
      b% pgbinary_right_yaxis_label_disp = pgbinary_right_yaxis_label_disp
      b% pgbinary_num_scale = pgbinary_num_scale
      b% pgbinary_lw = pgbinary_lw
      b% pgbinary_model_lw = pgbinary_model_lw
      b% pgbinary_box_lw = pgbinary_box_lw

      b% Text_Summary1_win_flag = Text_Summary1_win_flag
      b% Text_Summary1_file_flag = Text_Summary1_file_flag
      b% Text_Summary1_file_interval = Text_Summary1_file_interval
      b% Text_Summary1_file_dir = Text_Summary1_file_dir
      b% Text_Summary1_file_prefix = Text_Summary1_file_prefix
      b% Text_Summary1_num_cols = Text_Summary1_num_cols
      b% Text_Summary1_num_rows = Text_Summary1_num_rows
      b% Text_Summary1_name = Text_Summary1_name
      b% Text_Summary1_win_width = Text_Summary1_win_width
      b% Text_Summary1_win_aspect_ratio = Text_Summary1_win_aspect_ratio
      b% Text_Summary1_file_width = Text_Summary1_file_width
      b% Text_Summary1_file_aspect_ratio = Text_Summary1_file_aspect_ratio
      b% Text_Summary1_title = Text_Summary1_title
      b% Text_Summary1_xleft = Text_Summary1_xleft
      b% Text_Summary1_xright = Text_Summary1_xright
      b% Text_Summary1_ybot = Text_Summary1_ybot
      b% Text_Summary1_ytop = Text_Summary1_ytop
      b% Text_Summary1_txt_scale = Text_Summary1_txt_scale

      b% Text_Summary2_win_flag = Text_Summary2_win_flag
      b% Text_Summary2_file_flag = Text_Summary2_file_flag
      b% Text_Summary2_file_interval = Text_Summary2_file_interval
      b% Text_Summary2_file_dir = Text_Summary2_file_dir
      b% Text_Summary2_file_prefix = Text_Summary2_file_prefix
      b% Text_Summary2_num_cols = Text_Summary2_num_cols
      b% Text_Summary2_num_rows = Text_Summary2_num_rows
      b% Text_Summary2_name = Text_Summary2_name
      b% Text_Summary2_win_width = Text_Summary2_win_width
      b% Text_Summary2_win_aspect_ratio = Text_Summary2_win_aspect_ratio
      b% Text_Summary2_file_width = Text_Summary2_file_width
      b% Text_Summary2_file_aspect_ratio = Text_Summary2_file_aspect_ratio
      b% Text_Summary2_title = Text_Summary2_title
      b% Text_Summary2_xleft = Text_Summary2_xleft
      b% Text_Summary2_xright = Text_Summary2_xright
      b% Text_Summary2_ybot = Text_Summary2_ybot
      b% Text_Summary2_ytop = Text_Summary2_ytop
      b% Text_Summary2_txt_scale = Text_Summary2_txt_scale

      b% Text_Summary3_win_flag = Text_Summary3_win_flag
      b% Text_Summary3_file_flag = Text_Summary3_file_flag
      b% Text_Summary3_file_interval = Text_Summary3_file_interval
      b% Text_Summary3_file_dir = Text_Summary3_file_dir
      b% Text_Summary3_file_prefix = Text_Summary3_file_prefix
      b% Text_Summary3_num_cols = Text_Summary3_num_cols
      b% Text_Summary3_num_rows = Text_Summary3_num_rows
      b% Text_Summary3_name = Text_Summary3_name
      b% Text_Summary3_win_width = Text_Summary3_win_width
      b% Text_Summary3_win_aspect_ratio = Text_Summary3_win_aspect_ratio
      b% Text_Summary3_file_width = Text_Summary3_file_width
      b% Text_Summary3_file_aspect_ratio = Text_Summary3_file_aspect_ratio
      b% Text_Summary3_title = Text_Summary3_title
      b% Text_Summary3_xleft = Text_Summary3_xleft
      b% Text_Summary3_xright = Text_Summary3_xright
      b% Text_Summary3_ybot = Text_Summary3_ybot
      b% Text_Summary3_ytop = Text_Summary3_ytop
      b% Text_Summary3_txt_scale = Text_Summary3_txt_scale

      b% Text_Summary4_win_flag = Text_Summary4_win_flag
      b% Text_Summary4_file_flag = Text_Summary4_file_flag
      b% Text_Summary4_file_interval = Text_Summary4_file_interval
      b% Text_Summary4_file_dir = Text_Summary4_file_dir
      b% Text_Summary4_file_prefix = Text_Summary4_file_prefix
      b% Text_Summary4_num_cols = Text_Summary4_num_cols
      b% Text_Summary4_num_rows = Text_Summary4_num_rows
      b% Text_Summary4_name = Text_Summary4_name
      b% Text_Summary4_win_width = Text_Summary4_win_width
      b% Text_Summary4_win_aspect_ratio = Text_Summary4_win_aspect_ratio
      b% Text_Summary4_file_width = Text_Summary4_file_width
      b% Text_Summary4_file_aspect_ratio = Text_Summary4_file_aspect_ratio
      b% Text_Summary4_title = Text_Summary4_title
      b% Text_Summary4_xleft = Text_Summary4_xleft
      b% Text_Summary4_xright = Text_Summary4_xright
      b% Text_Summary4_ybot = Text_Summary4_ybot
      b% Text_Summary4_ytop = Text_Summary4_ytop
      b% Text_Summary4_txt_scale = Text_Summary4_txt_scale

      b% Text_Summary5_win_flag = Text_Summary5_win_flag
      b% Text_Summary5_file_flag = Text_Summary5_file_flag
      b% Text_Summary5_file_interval = Text_Summary5_file_interval
      b% Text_Summary5_file_dir = Text_Summary5_file_dir
      b% Text_Summary5_file_prefix = Text_Summary5_file_prefix
      b% Text_Summary5_num_cols = Text_Summary5_num_cols
      b% Text_Summary5_num_rows = Text_Summary5_num_rows
      b% Text_Summary5_name = Text_Summary5_name
      b% Text_Summary5_win_width = Text_Summary5_win_width
      b% Text_Summary5_win_aspect_ratio = Text_Summary5_win_aspect_ratio
      b% Text_Summary5_file_width = Text_Summary5_file_width
      b% Text_Summary5_file_aspect_ratio = Text_Summary5_file_aspect_ratio
      b% Text_Summary5_title = Text_Summary5_title
      b% Text_Summary5_xleft = Text_Summary5_xleft
      b% Text_Summary5_xright = Text_Summary5_xright
      b% Text_Summary5_ybot = Text_Summary5_ybot
      b% Text_Summary5_ytop = Text_Summary5_ytop
      b% Text_Summary5_txt_scale = Text_Summary5_txt_scale

      b% Text_Summary6_win_flag = Text_Summary6_win_flag
      b% Text_Summary6_file_flag = Text_Summary6_file_flag
      b% Text_Summary6_file_interval = Text_Summary6_file_interval
      b% Text_Summary6_file_dir = Text_Summary6_file_dir
      b% Text_Summary6_file_prefix = Text_Summary6_file_prefix
      b% Text_Summary6_num_cols = Text_Summary6_num_cols
      b% Text_Summary6_num_rows = Text_Summary6_num_rows
      b% Text_Summary6_name = Text_Summary6_name
      b% Text_Summary6_win_width = Text_Summary6_win_width
      b% Text_Summary6_win_aspect_ratio = Text_Summary6_win_aspect_ratio
      b% Text_Summary6_file_width = Text_Summary6_file_width
      b% Text_Summary6_file_aspect_ratio = Text_Summary6_file_aspect_ratio
      b% Text_Summary6_title = Text_Summary6_title
      b% Text_Summary6_xleft = Text_Summary6_xleft
      b% Text_Summary6_xright = Text_Summary6_xright
      b% Text_Summary6_ybot = Text_Summary6_ybot
      b% Text_Summary6_ytop = Text_Summary6_ytop
      b% Text_Summary6_txt_scale = Text_Summary6_txt_scale

      b% Text_Summary7_win_flag = Text_Summary7_win_flag
      b% Text_Summary7_file_flag = Text_Summary7_file_flag
      b% Text_Summary7_file_interval = Text_Summary7_file_interval
      b% Text_Summary7_file_dir = Text_Summary7_file_dir
      b% Text_Summary7_file_prefix = Text_Summary7_file_prefix
      b% Text_Summary7_num_cols = Text_Summary7_num_cols
      b% Text_Summary7_num_rows = Text_Summary7_num_rows
      b% Text_Summary7_name = Text_Summary7_name
      b% Text_Summary7_win_width = Text_Summary7_win_width
      b% Text_Summary7_win_aspect_ratio = Text_Summary7_win_aspect_ratio
      b% Text_Summary7_file_width = Text_Summary7_file_width
      b% Text_Summary7_file_aspect_ratio = Text_Summary7_file_aspect_ratio
      b% Text_Summary7_title = Text_Summary7_title
      b% Text_Summary7_xleft = Text_Summary7_xleft
      b% Text_Summary7_xright = Text_Summary7_xright
      b% Text_Summary7_ybot = Text_Summary7_ybot
      b% Text_Summary7_ytop = Text_Summary7_ytop
      b% Text_Summary7_txt_scale = Text_Summary7_txt_scale

      b% Text_Summary8_win_flag = Text_Summary8_win_flag
      b% Text_Summary8_file_flag = Text_Summary8_file_flag
      b% Text_Summary8_file_interval = Text_Summary8_file_interval
      b% Text_Summary8_file_dir = Text_Summary8_file_dir
      b% Text_Summary8_file_prefix = Text_Summary8_file_prefix
      b% Text_Summary8_num_cols = Text_Summary8_num_cols
      b% Text_Summary8_num_rows = Text_Summary8_num_rows
      b% Text_Summary8_name = Text_Summary8_name
      b% Text_Summary8_win_width = Text_Summary8_win_width
      b% Text_Summary8_win_aspect_ratio = Text_Summary8_win_aspect_ratio
      b% Text_Summary8_file_width = Text_Summary8_file_width
      b% Text_Summary8_file_aspect_ratio = Text_Summary8_file_aspect_ratio
      b% Text_Summary8_title = Text_Summary8_title
      b% Text_Summary8_xleft = Text_Summary8_xleft
      b% Text_Summary8_xright = Text_Summary8_xright
      b% Text_Summary8_ybot = Text_Summary8_ybot
      b% Text_Summary8_ytop = Text_Summary8_ytop
      b% Text_Summary8_txt_scale = Text_Summary8_txt_scale

      b% Text_Summary9_win_flag = Text_Summary9_win_flag
      b% Text_Summary9_file_flag = Text_Summary9_file_flag
      b% Text_Summary9_file_interval = Text_Summary9_file_interval
      b% Text_Summary9_file_dir = Text_Summary9_file_dir
      b% Text_Summary9_file_prefix = Text_Summary9_file_prefix
      b% Text_Summary9_num_cols = Text_Summary9_num_cols
      b% Text_Summary9_num_rows = Text_Summary9_num_rows
      b% Text_Summary9_name = Text_Summary9_name
      b% Text_Summary9_win_width = Text_Summary9_win_width
      b% Text_Summary9_win_aspect_ratio = Text_Summary9_win_aspect_ratio
      b% Text_Summary9_file_width = Text_Summary9_file_width
      b% Text_Summary9_file_aspect_ratio = Text_Summary9_file_aspect_ratio
      b% Text_Summary9_title = Text_Summary9_title
      b% Text_Summary9_xleft = Text_Summary9_xleft
      b% Text_Summary9_xright = Text_Summary9_xright
      b% Text_Summary9_ybot = Text_Summary9_ybot
      b% Text_Summary9_ytop = Text_Summary9_ytop
      b% Text_Summary9_txt_scale = Text_Summary9_txt_scale

      b% History_Track1_win_flag = History_Track1_win_flag
      b% History_Track1_file_flag = History_Track1_file_flag
      b% History_Track1_file_interval = History_Track1_file_interval
      b% History_Track1_step_min = History_Track1_step_min
      b% History_Track1_step_max = History_Track1_step_max
      b% show_History_Track1_target_box = show_History_Track1_target_box
      b% History_Track1_n_sigma = History_Track1_n_sigma
      b% History_Track1_xtarget = History_Track1_xtarget
      b% History_Track1_xsigma = History_Track1_xsigma
      b% History_Track1_ytarget = History_Track1_ytarget
      b% History_Track1_ysigma = History_Track1_ysigma
      b% History_Track1_file_dir = History_Track1_file_dir
      b% History_Track1_file_prefix = History_Track1_file_prefix
      b% show_History_Track1_annotation1 = show_History_Track1_annotation1
      b% show_History_Track1_annotation2 = show_History_Track1_annotation2
      b% show_History_Track1_annotation3 = show_History_Track1_annotation3
      b% History_Track1_fname = History_Track1_fname
      b% History_Track1_xname = History_Track1_xname
      b% History_Track1_xaxis_label = History_Track1_xaxis_label
      b% History_Track1_yname = History_Track1_yname
      b% History_Track1_yaxis_label = History_Track1_yaxis_label
      b% History_Track1_reverse_xaxis = History_Track1_reverse_xaxis
      b% History_Track1_reverse_yaxis = History_Track1_reverse_yaxis
      b% History_Track1_log_xaxis = History_Track1_log_xaxis
      b% History_Track1_log_yaxis = History_Track1_log_yaxis
      b% History_Track1_xmin = History_Track1_xmin
      b% History_Track1_xmax = History_Track1_xmax
      b% History_Track1_ymin = History_Track1_ymin
      b% History_Track1_ymax = History_Track1_ymax
      b% History_Track1_xmargin = History_Track1_xmargin
      b% History_Track1_ymargin = History_Track1_ymargin
      b% History_Track1_dxmin = History_Track1_dxmin
      b% History_Track1_dymin = History_Track1_dymin
      b% History_Track1_win_width = History_Track1_win_width
      b% History_Track1_win_aspect_ratio = History_Track1_win_aspect_ratio
      b% History_Track1_file_width = History_Track1_file_width
      b% History_Track1_file_aspect_ratio = History_Track1_file_aspect_ratio
      b% History_Track1_xleft = History_Track1_xleft
      b% History_Track1_xright = History_Track1_xright
      b% History_Track1_ybot = History_Track1_ybot
      b% History_Track1_ytop = History_Track1_ytop
      b% History_Track1_txt_scale = History_Track1_txt_scale
      b% History_Track1_title = History_Track1_title
      b% History_Track1_use_decorator = History_Track1_use_decorator

      b% History_Track2_win_flag = History_Track2_win_flag
      b% History_Track2_file_flag = History_Track2_file_flag
      b% History_Track2_file_interval = History_Track2_file_interval
      b% History_Track2_step_min = History_Track2_step_min
      b% History_Track2_step_max = History_Track2_step_max
      b% show_History_Track2_target_box = show_History_Track2_target_box
      b% History_Track2_n_sigma = History_Track2_n_sigma
      b% History_Track2_xtarget = History_Track2_xtarget
      b% History_Track2_xsigma = History_Track2_xsigma
      b% History_Track2_ytarget = History_Track2_ytarget
      b% History_Track2_ysigma = History_Track2_ysigma
      b% History_Track2_xname = History_Track2_xname
      b% History_Track2_xaxis_label = History_Track2_xaxis_label
      b% History_Track2_yname = History_Track2_yname
      b% History_Track2_yaxis_label = History_Track2_yaxis_label
      b% History_Track2_file_dir = History_Track2_file_dir
      b% History_Track2_file_prefix = History_Track2_file_prefix
      b% show_History_Track2_annotation1 = show_History_Track2_annotation1
      b% show_History_Track2_annotation2 = show_History_Track2_annotation2
      b% show_History_Track2_annotation3 = show_History_Track2_annotation3
      b% History_Track2_fname = History_Track2_fname
      b% History_Track2_reverse_xaxis = History_Track2_reverse_xaxis
      b% History_Track2_reverse_yaxis = History_Track2_reverse_yaxis
      b% History_Track2_log_xaxis = History_Track2_log_xaxis
      b% History_Track2_log_yaxis = History_Track2_log_yaxis
      b% History_Track2_xmin = History_Track2_xmin
      b% History_Track2_xmax = History_Track2_xmax
      b% History_Track2_ymin = History_Track2_ymin
      b% History_Track2_ymax = History_Track2_ymax
      b% History_Track2_xmargin = History_Track2_xmargin
      b% History_Track2_ymargin = History_Track2_ymargin
      b% History_Track2_dxmin = History_Track2_dxmin
      b% History_Track2_dymin = History_Track2_dymin
      b% History_Track2_win_width = History_Track2_win_width
      b% History_Track2_win_aspect_ratio = History_Track2_win_aspect_ratio
      b% History_Track2_file_width = History_Track2_file_width
      b% History_Track2_file_aspect_ratio = History_Track2_file_aspect_ratio
      b% History_Track2_xleft = History_Track2_xleft
      b% History_Track2_xright = History_Track2_xright
      b% History_Track2_ybot = History_Track2_ybot
      b% History_Track2_ytop = History_Track2_ytop
      b% History_Track2_txt_scale = History_Track2_txt_scale
      b% History_Track2_title = History_Track2_title
      b% History_Track2_use_decorator = History_Track2_use_decorator

      b% History_Track3_win_flag = History_Track3_win_flag
      b% History_Track3_file_flag = History_Track3_file_flag
      b% History_Track3_file_interval = History_Track3_file_interval
      b% History_Track3_step_min = History_Track3_step_min
      b% History_Track3_step_max = History_Track3_step_max
      b% show_History_Track3_target_box = show_History_Track3_target_box
      b% History_Track3_n_sigma = History_Track3_n_sigma
      b% History_Track3_xtarget = History_Track3_xtarget
      b% History_Track3_xsigma = History_Track3_xsigma
      b% History_Track3_ytarget = History_Track3_ytarget
      b% History_Track3_ysigma = History_Track3_ysigma
      b% History_Track3_xname = History_Track3_xname
      b% History_Track3_xaxis_label = History_Track3_xaxis_label
      b% History_Track3_yname = History_Track3_yname
      b% History_Track3_yaxis_label = History_Track3_yaxis_label
      b% History_Track3_file_dir = History_Track3_file_dir
      b% History_Track3_file_prefix = History_Track3_file_prefix
      b% show_History_Track3_annotation1 = show_History_Track3_annotation1
      b% show_History_Track3_annotation2 = show_History_Track3_annotation2
      b% show_History_Track3_annotation3 = show_History_Track3_annotation3
      b% History_Track3_fname = History_Track3_fname
      b% History_Track3_reverse_xaxis = History_Track3_reverse_xaxis
      b% History_Track3_reverse_yaxis = History_Track3_reverse_yaxis
      b% History_Track3_log_xaxis = History_Track3_log_xaxis
      b% History_Track3_log_yaxis = History_Track3_log_yaxis
      b% History_Track3_xmin = History_Track3_xmin
      b% History_Track3_xmax = History_Track3_xmax
      b% History_Track3_ymin = History_Track3_ymin
      b% History_Track3_ymax = History_Track3_ymax
      b% History_Track3_xmargin = History_Track3_xmargin
      b% History_Track3_ymargin = History_Track3_ymargin
      b% History_Track3_dxmin = History_Track3_dxmin
      b% History_Track3_dymin = History_Track3_dymin
      b% History_Track3_win_width = History_Track3_win_width
      b% History_Track3_win_aspect_ratio = History_Track3_win_aspect_ratio
      b% History_Track3_file_width = History_Track3_file_width
      b% History_Track3_file_aspect_ratio = History_Track3_file_aspect_ratio
      b% History_Track3_xleft = History_Track3_xleft
      b% History_Track3_xright = History_Track3_xright
      b% History_Track3_ybot = History_Track3_ybot
      b% History_Track3_ytop = History_Track3_ytop
      b% History_Track3_txt_scale = History_Track3_txt_scale
      b% History_Track3_title = History_Track3_title
      b% History_Track3_use_decorator = History_Track3_use_decorator

      b% History_Track4_win_flag = History_Track4_win_flag
      b% History_Track4_file_flag = History_Track4_file_flag
      b% History_Track4_file_interval = History_Track4_file_interval
      b% History_Track4_step_min = History_Track4_step_min
      b% History_Track4_step_max = History_Track4_step_max
      b% show_History_Track4_target_box = show_History_Track4_target_box
      b% History_Track4_n_sigma = History_Track4_n_sigma
      b% History_Track4_xtarget = History_Track4_xtarget
      b% History_Track4_xsigma = History_Track4_xsigma
      b% History_Track4_ytarget = History_Track4_ytarget
      b% History_Track4_ysigma = History_Track4_ysigma
      b% History_Track4_xname = History_Track4_xname
      b% History_Track4_xaxis_label = History_Track4_xaxis_label
      b% History_Track4_yname = History_Track4_yname
      b% History_Track4_yaxis_label = History_Track4_yaxis_label
      b% History_Track4_file_dir = History_Track4_file_dir
      b% History_Track4_file_prefix = History_Track4_file_prefix
      b% show_History_Track4_annotation1 = show_History_Track4_annotation1
      b% show_History_Track4_annotation2 = show_History_Track4_annotation2
      b% show_History_Track4_annotation3 = show_History_Track4_annotation3
      b% History_Track4_fname = History_Track4_fname
      b% History_Track4_reverse_xaxis = History_Track4_reverse_xaxis
      b% History_Track4_reverse_yaxis = History_Track4_reverse_yaxis
      b% History_Track4_log_xaxis = History_Track4_log_xaxis
      b% History_Track4_log_yaxis = History_Track4_log_yaxis
      b% History_Track4_xmin = History_Track4_xmin
      b% History_Track4_xmax = History_Track4_xmax
      b% History_Track4_ymin = History_Track4_ymin
      b% History_Track4_ymax = History_Track4_ymax
      b% History_Track4_xmargin = History_Track4_xmargin
      b% History_Track4_ymargin = History_Track4_ymargin
      b% History_Track4_dxmin = History_Track4_dxmin
      b% History_Track4_dymin = History_Track4_dymin
      b% History_Track4_win_width = History_Track4_win_width
      b% History_Track4_win_aspect_ratio = History_Track4_win_aspect_ratio
      b% History_Track4_file_width = History_Track4_file_width
      b% History_Track4_file_aspect_ratio = History_Track4_file_aspect_ratio
      b% History_Track4_xleft = History_Track4_xleft
      b% History_Track4_xright = History_Track4_xright
      b% History_Track4_ybot = History_Track4_ybot
      b% History_Track4_ytop = History_Track4_ytop
      b% History_Track4_txt_scale = History_Track4_txt_scale
      b% History_Track4_title = History_Track4_title
      b% History_Track4_use_decorator = History_Track4_use_decorator

      b% History_Track5_win_flag = History_Track5_win_flag
      b% History_Track5_file_flag = History_Track5_file_flag
      b% History_Track5_file_interval = History_Track5_file_interval
      b% History_Track5_step_min = History_Track5_step_min
      b% History_Track5_step_max = History_Track5_step_max
      b% show_History_Track5_target_box = show_History_Track5_target_box
      b% History_Track5_n_sigma = History_Track5_n_sigma
      b% History_Track5_xtarget = History_Track5_xtarget
      b% History_Track5_xsigma = History_Track5_xsigma
      b% History_Track5_ytarget = History_Track5_ytarget
      b% History_Track5_ysigma = History_Track5_ysigma
      b% History_Track5_xname = History_Track5_xname
      b% History_Track5_xaxis_label = History_Track5_xaxis_label
      b% History_Track5_yname = History_Track5_yname
      b% History_Track5_yaxis_label = History_Track5_yaxis_label
      b% History_Track5_file_dir = History_Track5_file_dir
      b% History_Track5_file_prefix = History_Track5_file_prefix
      b% show_History_Track5_annotation1 = show_History_Track5_annotation1
      b% show_History_Track5_annotation2 = show_History_Track5_annotation2
      b% show_History_Track5_annotation3 = show_History_Track5_annotation3
      b% History_Track5_fname = History_Track5_fname
      b% History_Track5_reverse_xaxis = History_Track5_reverse_xaxis
      b% History_Track5_reverse_yaxis = History_Track5_reverse_yaxis
      b% History_Track5_log_xaxis = History_Track5_log_xaxis
      b% History_Track5_log_yaxis = History_Track5_log_yaxis
      b% History_Track5_xmin = History_Track5_xmin
      b% History_Track5_xmax = History_Track5_xmax
      b% History_Track5_ymin = History_Track5_ymin
      b% History_Track5_ymax = History_Track5_ymax
      b% History_Track5_xmargin = History_Track5_xmargin
      b% History_Track5_ymargin = History_Track5_ymargin
      b% History_Track5_dxmin = History_Track5_dxmin
      b% History_Track5_dymin = History_Track5_dymin
      b% History_Track5_win_width = History_Track5_win_width
      b% History_Track5_win_aspect_ratio = History_Track5_win_aspect_ratio
      b% History_Track5_file_width = History_Track5_file_width
      b% History_Track5_file_aspect_ratio = History_Track5_file_aspect_ratio
      b% History_Track5_xleft = History_Track5_xleft
      b% History_Track5_xright = History_Track5_xright
      b% History_Track5_ybot = History_Track5_ybot
      b% History_Track5_ytop = History_Track5_ytop
      b% History_Track5_txt_scale = History_Track5_txt_scale
      b% History_Track5_title = History_Track5_title
      b% History_Track5_use_decorator = History_Track5_use_decorator

      b% History_Track6_win_flag = History_Track6_win_flag
      b% History_Track6_file_flag = History_Track6_file_flag
      b% History_Track6_file_interval = History_Track6_file_interval
      b% History_Track6_step_min = History_Track6_step_min
      b% History_Track6_step_max = History_Track6_step_max
      b% show_History_Track6_target_box = show_History_Track6_target_box
      b% History_Track6_n_sigma = History_Track6_n_sigma
      b% History_Track6_xtarget = History_Track6_xtarget
      b% History_Track6_xsigma = History_Track6_xsigma
      b% History_Track6_ytarget = History_Track6_ytarget
      b% History_Track6_ysigma = History_Track6_ysigma
      b% History_Track6_xname = History_Track6_xname
      b% History_Track6_xaxis_label = History_Track6_xaxis_label
      b% History_Track6_yname = History_Track6_yname
      b% History_Track6_yaxis_label = History_Track6_yaxis_label
      b% History_Track6_file_dir = History_Track6_file_dir
      b% History_Track6_file_prefix = History_Track6_file_prefix
      b% show_History_Track6_annotation1 = show_History_Track6_annotation1
      b% show_History_Track6_annotation2 = show_History_Track6_annotation2
      b% show_History_Track6_annotation3 = show_History_Track6_annotation3
      b% History_Track6_fname = History_Track6_fname
      b% History_Track6_reverse_xaxis = History_Track6_reverse_xaxis
      b% History_Track6_reverse_yaxis = History_Track6_reverse_yaxis
      b% History_Track6_log_xaxis = History_Track6_log_xaxis
      b% History_Track6_log_yaxis = History_Track6_log_yaxis
      b% History_Track6_xmin = History_Track6_xmin
      b% History_Track6_xmax = History_Track6_xmax
      b% History_Track6_ymin = History_Track6_ymin
      b% History_Track6_ymax = History_Track6_ymax
      b% History_Track6_xmargin = History_Track6_xmargin
      b% History_Track6_ymargin = History_Track6_ymargin
      b% History_Track6_dxmin = History_Track6_dxmin
      b% History_Track6_dymin = History_Track6_dymin
      b% History_Track6_win_width = History_Track6_win_width
      b% History_Track6_win_aspect_ratio = History_Track6_win_aspect_ratio
      b% History_Track6_file_width = History_Track6_file_width
      b% History_Track6_file_aspect_ratio = History_Track6_file_aspect_ratio
      b% History_Track6_xleft = History_Track6_xleft
      b% History_Track6_xright = History_Track6_xright
      b% History_Track6_ybot = History_Track6_ybot
      b% History_Track6_ytop = History_Track6_ytop
      b% History_Track6_txt_scale = History_Track6_txt_scale
      b% History_Track6_title = History_Track6_title
      b% History_Track6_use_decorator = History_Track6_use_decorator

      b% History_Track7_win_flag = History_Track7_win_flag
      b% History_Track7_file_flag = History_Track7_file_flag
      b% History_Track7_file_interval = History_Track7_file_interval
      b% History_Track7_step_min = History_Track7_step_min
      b% History_Track7_step_max = History_Track7_step_max
      b% show_History_Track7_target_box = show_History_Track7_target_box
      b% History_Track7_n_sigma = History_Track7_n_sigma
      b% History_Track7_xtarget = History_Track7_xtarget
      b% History_Track7_xsigma = History_Track7_xsigma
      b% History_Track7_ytarget = History_Track7_ytarget
      b% History_Track7_ysigma = History_Track7_ysigma
      b% History_Track7_xname = History_Track7_xname
      b% History_Track7_xaxis_label = History_Track7_xaxis_label
      b% History_Track7_yname = History_Track7_yname
      b% History_Track7_yaxis_label = History_Track7_yaxis_label
      b% History_Track7_file_dir = History_Track7_file_dir
      b% History_Track7_file_prefix = History_Track7_file_prefix
      b% show_History_Track7_annotation1 = show_History_Track7_annotation1
      b% show_History_Track7_annotation2 = show_History_Track7_annotation2
      b% show_History_Track7_annotation3 = show_History_Track7_annotation3
      b% History_Track7_fname = History_Track7_fname
      b% History_Track7_reverse_xaxis = History_Track7_reverse_xaxis
      b% History_Track7_reverse_yaxis = History_Track7_reverse_yaxis
      b% History_Track7_log_xaxis = History_Track7_log_xaxis
      b% History_Track7_log_yaxis = History_Track7_log_yaxis
      b% History_Track7_xmin = History_Track7_xmin
      b% History_Track7_xmax = History_Track7_xmax
      b% History_Track7_ymin = History_Track7_ymin
      b% History_Track7_ymax = History_Track7_ymax
      b% History_Track7_xmargin = History_Track7_xmargin
      b% History_Track7_ymargin = History_Track7_ymargin
      b% History_Track7_dxmin = History_Track7_dxmin
      b% History_Track7_dymin = History_Track7_dymin
      b% History_Track7_win_width = History_Track7_win_width
      b% History_Track7_win_aspect_ratio = History_Track7_win_aspect_ratio
      b% History_Track7_file_width = History_Track7_file_width
      b% History_Track7_file_aspect_ratio = History_Track7_file_aspect_ratio
      b% History_Track7_xleft = History_Track7_xleft
      b% History_Track7_xright = History_Track7_xright
      b% History_Track7_ybot = History_Track7_ybot
      b% History_Track7_ytop = History_Track7_ytop
      b% History_Track7_txt_scale = History_Track7_txt_scale
      b% History_Track7_title = History_Track7_title
      b% History_Track7_use_decorator = History_Track7_use_decorator

      b% History_Track8_win_flag = History_Track8_win_flag
      b% History_Track8_file_flag = History_Track8_file_flag
      b% History_Track8_file_interval = History_Track8_file_interval
      b% History_Track8_step_min = History_Track8_step_min
      b% History_Track8_step_max = History_Track8_step_max
      b% show_History_Track8_target_box = show_History_Track8_target_box
      b% History_Track8_n_sigma = History_Track8_n_sigma
      b% History_Track8_xtarget = History_Track8_xtarget
      b% History_Track8_xsigma = History_Track8_xsigma
      b% History_Track8_ytarget = History_Track8_ytarget
      b% History_Track8_ysigma = History_Track8_ysigma
      b% History_Track8_xname = History_Track8_xname
      b% History_Track8_xaxis_label = History_Track8_xaxis_label
      b% History_Track8_yname = History_Track8_yname
      b% History_Track8_yaxis_label = History_Track8_yaxis_label
      b% History_Track8_file_dir = History_Track8_file_dir
      b% History_Track8_file_prefix = History_Track8_file_prefix
      b% show_History_Track8_annotation1 = show_History_Track8_annotation1
      b% show_History_Track8_annotation2 = show_History_Track8_annotation2
      b% show_History_Track8_annotation3 = show_History_Track8_annotation3
      b% History_Track8_fname = History_Track8_fname
      b% History_Track8_reverse_xaxis = History_Track8_reverse_xaxis
      b% History_Track8_reverse_yaxis = History_Track8_reverse_yaxis
      b% History_Track8_log_xaxis = History_Track8_log_xaxis
      b% History_Track8_log_yaxis = History_Track8_log_yaxis
      b% History_Track8_xmin = History_Track8_xmin
      b% History_Track8_xmax = History_Track8_xmax
      b% History_Track8_ymin = History_Track8_ymin
      b% History_Track8_ymax = History_Track8_ymax
      b% History_Track8_xmargin = History_Track8_xmargin
      b% History_Track8_ymargin = History_Track8_ymargin
      b% History_Track8_dxmin = History_Track8_dxmin
      b% History_Track8_dymin = History_Track8_dymin
      b% History_Track8_win_width = History_Track8_win_width
      b% History_Track8_win_aspect_ratio = History_Track8_win_aspect_ratio
      b% History_Track8_file_width = History_Track8_file_width
      b% History_Track8_file_aspect_ratio = History_Track8_file_aspect_ratio
      b% History_Track8_xleft = History_Track8_xleft
      b% History_Track8_xright = History_Track8_xright
      b% History_Track8_ybot = History_Track8_ybot
      b% History_Track8_ytop = History_Track8_ytop
      b% History_Track8_txt_scale = History_Track8_txt_scale
      b% History_Track8_title = History_Track8_title
      b% History_Track8_use_decorator = History_Track8_use_decorator

      b% History_Track9_win_flag = History_Track9_win_flag
      b% History_Track9_file_flag = History_Track9_file_flag
      b% History_Track9_file_interval = History_Track9_file_interval
      b% History_Track9_step_min = History_Track9_step_min
      b% History_Track9_step_max = History_Track9_step_max
      b% show_History_Track9_target_box = show_History_Track9_target_box
      b% History_Track9_n_sigma = History_Track9_n_sigma
      b% History_Track9_xtarget = History_Track9_xtarget
      b% History_Track9_xsigma = History_Track9_xsigma
      b% History_Track9_ytarget = History_Track9_ytarget
      b% History_Track9_ysigma = History_Track9_ysigma
      b% History_Track9_xname = History_Track9_xname
      b% History_Track9_xaxis_label = History_Track9_xaxis_label
      b% History_Track9_yname = History_Track9_yname
      b% History_Track9_yaxis_label = History_Track9_yaxis_label
      b% History_Track9_file_dir = History_Track9_file_dir
      b% History_Track9_file_prefix = History_Track9_file_prefix
      b% show_History_Track9_annotation1 = show_History_Track9_annotation1
      b% show_History_Track9_annotation2 = show_History_Track9_annotation2
      b% show_History_Track9_annotation3 = show_History_Track9_annotation3
      b% History_Track9_fname = History_Track9_fname
      b% History_Track9_reverse_xaxis = History_Track9_reverse_xaxis
      b% History_Track9_reverse_yaxis = History_Track9_reverse_yaxis
      b% History_Track9_log_xaxis = History_Track9_log_xaxis
      b% History_Track9_log_yaxis = History_Track9_log_yaxis
      b% History_Track9_xmin = History_Track9_xmin
      b% History_Track9_xmax = History_Track9_xmax
      b% History_Track9_ymin = History_Track9_ymin
      b% History_Track9_ymax = History_Track9_ymax
      b% History_Track9_xmargin = History_Track9_xmargin
      b% History_Track9_ymargin = History_Track9_ymargin
      b% History_Track9_dxmin = History_Track9_dxmin
      b% History_Track9_dymin = History_Track9_dymin
      b% History_Track9_win_width = History_Track9_win_width
      b% History_Track9_win_aspect_ratio = History_Track9_win_aspect_ratio
      b% History_Track9_file_width = History_Track9_file_width
      b% History_Track9_file_aspect_ratio = History_Track9_file_aspect_ratio
      b% History_Track9_xleft = History_Track9_xleft
      b% History_Track9_xright = History_Track9_xright
      b% History_Track9_ybot = History_Track9_ybot
      b% History_Track9_ytop = History_Track9_ytop
      b% History_Track9_txt_scale = History_Track9_txt_scale
      b% History_Track9_title = History_Track9_title
      b% History_Track9_use_decorator = History_Track9_use_decorator

      b% History_Panels1_win_flag = History_Panels1_win_flag
      b% History_Panels1_win_width = History_Panels1_win_width
      b% History_Panels1_win_aspect_ratio = History_Panels1_win_aspect_ratio
      b% History_Panels1_xleft = History_Panels1_xleft
      b% History_Panels1_xright = History_Panels1_xright
      b% History_Panels1_ybot = History_Panels1_ybot
      b% History_Panels1_ytop = History_Panels1_ytop
      b% History_Panels1_txt_scale = History_Panels1_txt_scale
      b% History_Panels1_title = History_Panels1_title
      b% History_Panels1_xmax = History_Panels1_xmax
      b% History_Panels1_xmin = History_Panels1_xmin
      b% History_Panels1_dxmin = History_Panels1_dxmin
      b% History_Panels1_max_width = History_Panels1_max_width
      b% History_Panels1_num_panels = History_Panels1_num_panels
      b% History_Panels1_xaxis_name = History_Panels1_xaxis_name
      b% History_Panels1_yaxis_name = History_Panels1_yaxis_name
      b% History_Panels1_xaxis_reversed = History_Panels1_xaxis_reversed
      b% History_Panels1_yaxis_reversed = History_Panels1_yaxis_reversed
      b% History_Panels1_xaxis_log = History_Panels1_xaxis_log
      b% History_Panels1_yaxis_log = History_Panels1_yaxis_log
      b% History_Panels1_ymin = History_Panels1_ymin
      b% History_Panels1_ymax = History_Panels1_ymax
      b% History_Panels1_dymin = History_Panels1_dymin
      b% History_Panels1_other_yaxis_name = History_Panels1_other_yaxis_name
      b% History_Panels1_other_yaxis_reversed = History_Panels1_other_yaxis_reversed
      b% History_Panels1_other_yaxis_log = History_Panels1_other_yaxis_log
      b% History_Panels1_other_ymin = History_Panels1_other_ymin
      b% History_Panels1_other_ymax = History_Panels1_other_ymax
      b% History_Panels1_other_dymin = History_Panels1_other_dymin
      b% History_Panels1_file_flag = History_Panels1_file_flag
      b% History_Panels1_points_name = History_Panels1_points_name
      b% History_Panels1_file_dir = History_Panels1_file_dir
      b% History_Panels1_file_prefix = History_Panels1_file_prefix
      b% History_Panels1_file_interval = History_Panels1_file_interval
      b% History_Panels1_file_width = History_Panels1_file_width
      b% History_Panels1_file_aspect_ratio = History_Panels1_file_aspect_ratio
      b% History_Panels1_xmargin = History_Panels1_xmargin
      b% History_Panels1_ymargin = History_Panels1_ymargin
      b% History_Panels1_other_ymargin = History_Panels1_other_ymargin
      b% History_Panels1_use_decorator = History_Panels1_use_decorator

      b% History_Panels2_win_flag = History_Panels2_win_flag
      b% History_Panels2_win_width = History_Panels2_win_width
      b% History_Panels2_win_aspect_ratio = History_Panels2_win_aspect_ratio
      b% History_Panels2_xleft = History_Panels2_xleft
      b% History_Panels2_xright = History_Panels2_xright
      b% History_Panels2_ybot = History_Panels2_ybot
      b% History_Panels2_ytop = History_Panels2_ytop
      b% History_Panels2_txt_scale = History_Panels2_txt_scale
      b% History_Panels2_title = History_Panels2_title
      b% History_Panels2_xmax = History_Panels2_xmax
      b% History_Panels2_xmin = History_Panels2_xmin
      b% History_Panels2_dxmin = History_Panels2_dxmin
      b% History_Panels2_max_width = History_Panels2_max_width
      b% History_Panels2_num_panels = History_Panels2_num_panels
      b% History_Panels2_xaxis_name = History_Panels2_xaxis_name
      b% History_Panels2_yaxis_name = History_Panels2_yaxis_name
      b% History_Panels2_xaxis_reversed = History_Panels2_xaxis_reversed
      b% History_Panels2_yaxis_reversed = History_Panels2_yaxis_reversed
      b% History_Panels2_xaxis_log = History_Panels2_xaxis_log
      b% History_Panels2_yaxis_log = History_Panels2_yaxis_log
      b% History_Panels2_ymin = History_Panels2_ymin
      b% History_Panels2_ymax = History_Panels2_ymax
      b% History_Panels2_dymin = History_Panels2_dymin
      b% History_Panels2_other_yaxis_name = History_Panels2_other_yaxis_name
      b% History_Panels2_other_yaxis_reversed = History_Panels2_other_yaxis_reversed
      b% History_Panels2_other_yaxis_log = History_Panels2_other_yaxis_log
      b% History_Panels2_other_ymin = History_Panels2_other_ymin
      b% History_Panels2_other_ymax = History_Panels2_other_ymax
      b% History_Panels2_other_dymin = History_Panels2_other_dymin
      b% History_Panels2_file_flag = History_Panels2_file_flag
      b% History_Panels2_points_name = History_Panels2_points_name
      b% History_Panels2_file_dir = History_Panels2_file_dir
      b% History_Panels2_file_prefix = History_Panels2_file_prefix
      b% History_Panels2_file_interval = History_Panels2_file_interval
      b% History_Panels2_file_width = History_Panels2_file_width
      b% History_Panels2_file_aspect_ratio = History_Panels2_file_aspect_ratio
      b% History_Panels2_xmargin = History_Panels2_xmargin
      b% History_Panels2_ymargin = History_Panels2_ymargin
      b% History_Panels2_other_ymargin = History_Panels2_other_ymargin
      b% History_Panels2_use_decorator = History_Panels2_use_decorator

      b% History_Panels3_win_flag = History_Panels3_win_flag
      b% History_Panels3_win_width = History_Panels3_win_width
      b% History_Panels3_win_aspect_ratio = History_Panels3_win_aspect_ratio
      b% History_Panels3_xleft = History_Panels3_xleft
      b% History_Panels3_xright = History_Panels3_xright
      b% History_Panels3_ybot = History_Panels3_ybot
      b% History_Panels3_ytop = History_Panels3_ytop
      b% History_Panels3_txt_scale = History_Panels3_txt_scale
      b% History_Panels3_title = History_Panels3_title
      b% History_Panels3_xmax = History_Panels3_xmax
      b% History_Panels3_xmin = History_Panels3_xmin
      b% History_Panels3_dxmin = History_Panels3_dxmin
      b% History_Panels3_max_width = History_Panels3_max_width
      b% History_Panels3_num_panels = History_Panels3_num_panels
      b% History_Panels3_xaxis_name = History_Panels3_xaxis_name
      b% History_Panels3_yaxis_name = History_Panels3_yaxis_name
      b% History_Panels3_xaxis_reversed = History_Panels3_xaxis_reversed
      b% History_Panels3_yaxis_reversed = History_Panels3_yaxis_reversed
      b% History_Panels3_xaxis_log = History_Panels3_xaxis_log
      b% History_Panels3_yaxis_log = History_Panels3_yaxis_log
      b% History_Panels3_ymin = History_Panels3_ymin
      b% History_Panels3_ymax = History_Panels3_ymax
      b% History_Panels3_dymin = History_Panels3_dymin
      b% History_Panels3_other_yaxis_name = History_Panels3_other_yaxis_name
      b% History_Panels3_other_yaxis_reversed = History_Panels3_other_yaxis_reversed
      b% History_Panels3_other_yaxis_log = History_Panels3_other_yaxis_log
      b% History_Panels3_other_ymin = History_Panels3_other_ymin
      b% History_Panels3_other_ymax = History_Panels3_other_ymax
      b% History_Panels3_other_dymin = History_Panels3_other_dymin
      b% History_Panels3_file_flag = History_Panels3_file_flag
      b% History_Panels3_points_name = History_Panels3_points_name
      b% History_Panels3_file_dir = History_Panels3_file_dir
      b% History_Panels3_file_prefix = History_Panels3_file_prefix
      b% History_Panels3_file_interval = History_Panels3_file_interval
      b% History_Panels3_file_width = History_Panels3_file_width
      b% History_Panels3_file_aspect_ratio = History_Panels3_file_aspect_ratio
      b% History_Panels3_xmargin = History_Panels3_xmargin
      b% History_Panels3_ymargin = History_Panels3_ymargin
      b% History_Panels3_other_ymargin = History_Panels3_other_ymargin
      b% History_Panels3_use_decorator = History_Panels3_use_decorator

      b% History_Panels4_win_flag = History_Panels4_win_flag
      b% History_Panels4_win_width = History_Panels4_win_width
      b% History_Panels4_win_aspect_ratio = History_Panels4_win_aspect_ratio
      b% History_Panels4_xleft = History_Panels4_xleft
      b% History_Panels4_xright = History_Panels4_xright
      b% History_Panels4_ybot = History_Panels4_ybot
      b% History_Panels4_ytop = History_Panels4_ytop
      b% History_Panels4_txt_scale = History_Panels4_txt_scale
      b% History_Panels4_title = History_Panels4_title
      b% History_Panels4_xmax = History_Panels4_xmax
      b% History_Panels4_xmin = History_Panels4_xmin
      b% History_Panels4_dxmin = History_Panels4_dxmin
      b% History_Panels4_max_width = History_Panels4_max_width
      b% History_Panels4_num_panels = History_Panels4_num_panels
      b% History_Panels4_xaxis_name = History_Panels4_xaxis_name
      b% History_Panels4_yaxis_name = History_Panels4_yaxis_name
      b% History_Panels4_xaxis_reversed = History_Panels4_xaxis_reversed
      b% History_Panels4_yaxis_reversed = History_Panels4_yaxis_reversed
      b% History_Panels4_xaxis_log = History_Panels4_xaxis_log
      b% History_Panels4_yaxis_log = History_Panels4_yaxis_log
      b% History_Panels4_ymin = History_Panels4_ymin
      b% History_Panels4_ymax = History_Panels4_ymax
      b% History_Panels4_dymin = History_Panels4_dymin
      b% History_Panels4_other_yaxis_name = History_Panels4_other_yaxis_name
      b% History_Panels4_other_yaxis_reversed = History_Panels4_other_yaxis_reversed
      b% History_Panels4_other_yaxis_log = History_Panels4_other_yaxis_log
      b% History_Panels4_other_ymin = History_Panels4_other_ymin
      b% History_Panels4_other_ymax = History_Panels4_other_ymax
      b% History_Panels4_other_dymin = History_Panels4_other_dymin
      b% History_Panels4_file_flag = History_Panels4_file_flag
      b% History_Panels4_points_name = History_Panels4_points_name
      b% History_Panels4_file_dir = History_Panels4_file_dir
      b% History_Panels4_file_prefix = History_Panels4_file_prefix
      b% History_Panels4_file_interval = History_Panels4_file_interval
      b% History_Panels4_file_width = History_Panels4_file_width
      b% History_Panels4_file_aspect_ratio = History_Panels4_file_aspect_ratio
      b% History_Panels4_xmargin = History_Panels4_xmargin
      b% History_Panels4_ymargin = History_Panels4_ymargin
      b% History_Panels4_other_ymargin = History_Panels4_other_ymargin
      b% History_Panels4_use_decorator = History_Panels4_use_decorator

      b% History_Panels5_win_flag = History_Panels5_win_flag
      b% History_Panels5_win_width = History_Panels5_win_width
      b% History_Panels5_win_aspect_ratio = History_Panels5_win_aspect_ratio
      b% History_Panels5_xleft = History_Panels5_xleft
      b% History_Panels5_xright = History_Panels5_xright
      b% History_Panels5_ybot = History_Panels5_ybot
      b% History_Panels5_ytop = History_Panels5_ytop
      b% History_Panels5_txt_scale = History_Panels5_txt_scale
      b% History_Panels5_title = History_Panels5_title
      b% History_Panels5_xmax = History_Panels5_xmax
      b% History_Panels5_xmin = History_Panels5_xmin
      b% History_Panels5_dxmin = History_Panels5_dxmin
      b% History_Panels5_max_width = History_Panels5_max_width
      b% History_Panels5_num_panels = History_Panels5_num_panels
      b% History_Panels5_xaxis_name = History_Panels5_xaxis_name
      b% History_Panels5_yaxis_name = History_Panels5_yaxis_name
      b% History_Panels5_xaxis_reversed = History_Panels5_xaxis_reversed
      b% History_Panels5_yaxis_reversed = History_Panels5_yaxis_reversed
      b% History_Panels5_xaxis_log = History_Panels5_xaxis_log
      b% History_Panels5_yaxis_log = History_Panels5_yaxis_log
      b% History_Panels5_ymin = History_Panels5_ymin
      b% History_Panels5_ymax = History_Panels5_ymax
      b% History_Panels5_dymin = History_Panels5_dymin
      b% History_Panels5_other_yaxis_name = History_Panels5_other_yaxis_name
      b% History_Panels5_other_yaxis_reversed = History_Panels5_other_yaxis_reversed
      b% History_Panels5_other_yaxis_log = History_Panels5_other_yaxis_log
      b% History_Panels5_other_ymin = History_Panels5_other_ymin
      b% History_Panels5_other_ymax = History_Panels5_other_ymax
      b% History_Panels5_other_dymin = History_Panels5_other_dymin
      b% History_Panels5_file_flag = History_Panels5_file_flag
      b% History_Panels5_points_name = History_Panels5_points_name
      b% History_Panels5_file_dir = History_Panels5_file_dir
      b% History_Panels5_file_prefix = History_Panels5_file_prefix
      b% History_Panels5_file_interval = History_Panels5_file_interval
      b% History_Panels5_file_width = History_Panels5_file_width
      b% History_Panels5_file_aspect_ratio = History_Panels5_file_aspect_ratio
      b% History_Panels5_xmargin = History_Panels5_xmargin
      b% History_Panels5_ymargin = History_Panels5_ymargin
      b% History_Panels5_other_ymargin = History_Panels5_other_ymargin
      b% History_Panels5_use_decorator = History_Panels5_use_decorator

      b% History_Panels6_win_flag = History_Panels6_win_flag
      b% History_Panels6_win_width = History_Panels6_win_width
      b% History_Panels6_win_aspect_ratio = History_Panels6_win_aspect_ratio
      b% History_Panels6_xleft = History_Panels6_xleft
      b% History_Panels6_xright = History_Panels6_xright
      b% History_Panels6_ybot = History_Panels6_ybot
      b% History_Panels6_ytop = History_Panels6_ytop
      b% History_Panels6_txt_scale = History_Panels6_txt_scale
      b% History_Panels6_title = History_Panels6_title
      b% History_Panels6_xmax = History_Panels6_xmax
      b% History_Panels6_xmin = History_Panels6_xmin
      b% History_Panels6_dxmin = History_Panels6_dxmin
      b% History_Panels6_max_width = History_Panels6_max_width
      b% History_Panels6_num_panels = History_Panels6_num_panels
      b% History_Panels6_xaxis_name = History_Panels6_xaxis_name
      b% History_Panels6_yaxis_name = History_Panels6_yaxis_name
      b% History_Panels6_xaxis_reversed = History_Panels6_xaxis_reversed
      b% History_Panels6_yaxis_reversed = History_Panels6_yaxis_reversed
      b% History_Panels6_xaxis_log = History_Panels6_xaxis_log
      b% History_Panels6_yaxis_log = History_Panels6_yaxis_log
      b% History_Panels6_ymin = History_Panels6_ymin
      b% History_Panels6_ymax = History_Panels6_ymax
      b% History_Panels6_dymin = History_Panels6_dymin
      b% History_Panels6_other_yaxis_name = History_Panels6_other_yaxis_name
      b% History_Panels6_other_yaxis_reversed = History_Panels6_other_yaxis_reversed
      b% History_Panels6_other_yaxis_log = History_Panels6_other_yaxis_log
      b% History_Panels6_other_ymin = History_Panels6_other_ymin
      b% History_Panels6_other_ymax = History_Panels6_other_ymax
      b% History_Panels6_other_dymin = History_Panels6_other_dymin
      b% History_Panels6_file_flag = History_Panels6_file_flag
      b% History_Panels6_points_name = History_Panels6_points_name
      b% History_Panels6_file_dir = History_Panels6_file_dir
      b% History_Panels6_file_prefix = History_Panels6_file_prefix
      b% History_Panels6_file_interval = History_Panels6_file_interval
      b% History_Panels6_file_width = History_Panels6_file_width
      b% History_Panels6_file_aspect_ratio = History_Panels6_file_aspect_ratio
      b% History_Panels6_xmargin = History_Panels6_xmargin
      b% History_Panels6_ymargin = History_Panels6_ymargin
      b% History_Panels6_other_ymargin = History_Panels6_other_ymargin
      b% History_Panels6_use_decorator = History_Panels6_use_decorator

      b% History_Panels7_win_flag = History_Panels7_win_flag
      b% History_Panels7_win_width = History_Panels7_win_width
      b% History_Panels7_win_aspect_ratio = History_Panels7_win_aspect_ratio
      b% History_Panels7_xleft = History_Panels7_xleft
      b% History_Panels7_xright = History_Panels7_xright
      b% History_Panels7_ybot = History_Panels7_ybot
      b% History_Panels7_ytop = History_Panels7_ytop
      b% History_Panels7_txt_scale = History_Panels7_txt_scale
      b% History_Panels7_title = History_Panels7_title
      b% History_Panels7_xmax = History_Panels7_xmax
      b% History_Panels7_xmin = History_Panels7_xmin
      b% History_Panels7_dxmin = History_Panels7_dxmin
      b% History_Panels7_max_width = History_Panels7_max_width
      b% History_Panels7_num_panels = History_Panels7_num_panels
      b% History_Panels7_xaxis_name = History_Panels7_xaxis_name
      b% History_Panels7_yaxis_name = History_Panels7_yaxis_name
      b% History_Panels7_xaxis_reversed = History_Panels7_xaxis_reversed
      b% History_Panels7_yaxis_reversed = History_Panels7_yaxis_reversed
      b% History_Panels7_xaxis_log = History_Panels7_xaxis_log
      b% History_Panels7_yaxis_log = History_Panels7_yaxis_log
      b% History_Panels7_ymin = History_Panels7_ymin
      b% History_Panels7_ymax = History_Panels7_ymax
      b% History_Panels7_dymin = History_Panels7_dymin
      b% History_Panels7_other_yaxis_name = History_Panels7_other_yaxis_name
      b% History_Panels7_other_yaxis_reversed = History_Panels7_other_yaxis_reversed
      b% History_Panels7_other_yaxis_log = History_Panels7_other_yaxis_log
      b% History_Panels7_other_ymin = History_Panels7_other_ymin
      b% History_Panels7_other_ymax = History_Panels7_other_ymax
      b% History_Panels7_other_dymin = History_Panels7_other_dymin
      b% History_Panels7_file_flag = History_Panels7_file_flag
      b% History_Panels7_points_name = History_Panels7_points_name
      b% History_Panels7_file_dir = History_Panels7_file_dir
      b% History_Panels7_file_prefix = History_Panels7_file_prefix
      b% History_Panels7_file_interval = History_Panels7_file_interval
      b% History_Panels7_file_width = History_Panels7_file_width
      b% History_Panels7_file_aspect_ratio = History_Panels7_file_aspect_ratio
      b% History_Panels7_xmargin = History_Panels7_xmargin
      b% History_Panels7_ymargin = History_Panels7_ymargin
      b% History_Panels7_other_ymargin = History_Panels7_other_ymargin
      b% History_Panels7_use_decorator = History_Panels7_use_decorator

      b% History_Panels8_win_flag = History_Panels8_win_flag
      b% History_Panels8_win_width = History_Panels8_win_width
      b% History_Panels8_win_aspect_ratio = History_Panels8_win_aspect_ratio
      b% History_Panels8_xleft = History_Panels8_xleft
      b% History_Panels8_xright = History_Panels8_xright
      b% History_Panels8_ybot = History_Panels8_ybot
      b% History_Panels8_ytop = History_Panels8_ytop
      b% History_Panels8_txt_scale = History_Panels8_txt_scale
      b% History_Panels8_title = History_Panels8_title
      b% History_Panels8_xmax = History_Panels8_xmax
      b% History_Panels8_xmin = History_Panels8_xmin
      b% History_Panels8_dxmin = History_Panels8_dxmin
      b% History_Panels8_max_width = History_Panels8_max_width
      b% History_Panels8_num_panels = History_Panels8_num_panels
      b% History_Panels8_xaxis_name = History_Panels8_xaxis_name
      b% History_Panels8_yaxis_name = History_Panels8_yaxis_name
      b% History_Panels8_xaxis_reversed = History_Panels8_xaxis_reversed
      b% History_Panels8_yaxis_reversed = History_Panels8_yaxis_reversed
      b% History_Panels8_xaxis_log = History_Panels8_xaxis_log
      b% History_Panels8_yaxis_log = History_Panels8_yaxis_log
      b% History_Panels8_ymin = History_Panels8_ymin
      b% History_Panels8_ymax = History_Panels8_ymax
      b% History_Panels8_dymin = History_Panels8_dymin
      b% History_Panels8_other_yaxis_name = History_Panels8_other_yaxis_name
      b% History_Panels8_other_yaxis_reversed = History_Panels8_other_yaxis_reversed
      b% History_Panels8_other_yaxis_log = History_Panels8_other_yaxis_log
      b% History_Panels8_other_ymin = History_Panels8_other_ymin
      b% History_Panels8_other_ymax = History_Panels8_other_ymax
      b% History_Panels8_other_dymin = History_Panels8_other_dymin
      b% History_Panels8_file_flag = History_Panels8_file_flag
      b% History_Panels8_points_name = History_Panels8_points_name
      b% History_Panels8_file_dir = History_Panels8_file_dir
      b% History_Panels8_file_prefix = History_Panels8_file_prefix
      b% History_Panels8_file_interval = History_Panels8_file_interval
      b% History_Panels8_file_width = History_Panels8_file_width
      b% History_Panels8_file_aspect_ratio = History_Panels8_file_aspect_ratio
      b% History_Panels8_xmargin = History_Panels8_xmargin
      b% History_Panels8_ymargin = History_Panels8_ymargin
      b% History_Panels8_other_ymargin = History_Panels8_other_ymargin
      b% History_Panels8_use_decorator = History_Panels8_use_decorator

      b% History_Panels9_win_flag = History_Panels9_win_flag
      b% History_Panels9_win_width = History_Panels9_win_width
      b% History_Panels9_win_aspect_ratio = History_Panels9_win_aspect_ratio
      b% History_Panels9_xleft = History_Panels9_xleft
      b% History_Panels9_xright = History_Panels9_xright
      b% History_Panels9_ybot = History_Panels9_ybot
      b% History_Panels9_ytop = History_Panels9_ytop
      b% History_Panels9_txt_scale = History_Panels9_txt_scale
      b% History_Panels9_title = History_Panels9_title
      b% History_Panels9_xmax = History_Panels9_xmax
      b% History_Panels9_xmin = History_Panels9_xmin
      b% History_Panels9_dxmin = History_Panels9_dxmin
      b% History_Panels9_max_width = History_Panels9_max_width
      b% History_Panels9_num_panels = History_Panels9_num_panels
      b% History_Panels9_xaxis_name = History_Panels9_xaxis_name
      b% History_Panels9_yaxis_name = History_Panels9_yaxis_name
      b% History_Panels9_xaxis_reversed = History_Panels9_xaxis_reversed
      b% History_Panels9_yaxis_reversed = History_Panels9_yaxis_reversed
      b% History_Panels9_xaxis_log = History_Panels9_xaxis_log
      b% History_Panels9_yaxis_log = History_Panels9_yaxis_log
      b% History_Panels9_ymin = History_Panels9_ymin
      b% History_Panels9_ymax = History_Panels9_ymax
      b% History_Panels9_dymin = History_Panels9_dymin
      b% History_Panels9_other_yaxis_name = History_Panels9_other_yaxis_name
      b% History_Panels9_other_yaxis_reversed = History_Panels9_other_yaxis_reversed
      b% History_Panels9_other_yaxis_log = History_Panels9_other_yaxis_log
      b% History_Panels9_other_ymin = History_Panels9_other_ymin
      b% History_Panels9_other_ymax = History_Panels9_other_ymax
      b% History_Panels9_other_dymin = History_Panels9_other_dymin
      b% History_Panels9_file_flag = History_Panels9_file_flag
      b% History_Panels9_points_name = History_Panels9_points_name
      b% History_Panels9_file_flag = History_Panels9_file_flag
      b% History_Panels9_file_dir = History_Panels9_file_dir
      b% History_Panels9_file_prefix = History_Panels9_file_prefix
      b% History_Panels9_file_interval = History_Panels9_file_interval
      b% History_Panels9_file_width = History_Panels9_file_width
      b% History_Panels9_file_aspect_ratio = History_Panels9_file_aspect_ratio
      b% History_Panels9_xmargin = History_Panels9_xmargin
      b% History_Panels9_ymargin = History_Panels9_ymargin
      b% History_Panels9_other_ymargin = History_Panels9_other_ymargin
      b% History_Panels9_use_decorator = History_Panels9_use_decorator

      b% History_Panel_points_error_bars = History_Panel_points_error_bars
      b% History_Panel_points_interval = History_Panel_points_interval
      b% History_Panel_points_marker = History_Panel_points_marker
      b% History_Panel_points_ci = History_Panel_points_ci
      b% History_Panel_points_lw = History_Panel_points_lw
      b% History_Panel_points_ch = History_Panel_points_ch

      b% Summary_History_win_flag = Summary_History_win_flag
      b% Summary_History_file_flag = Summary_History_file_flag
      b% Summary_History_file_interval = Summary_History_file_interval
      b% Summary_History_file_dir = Summary_History_file_dir
      b% Summary_History_file_prefix = Summary_History_file_prefix
      b% Summary_History_scaled_value = Summary_History_scaled_value
      b% Summary_History_xmin = Summary_History_xmin
      b% Summary_History_xmax = Summary_History_xmax
      b% Summary_History_max_width = Summary_History_max_width
      b% Summary_History_win_width = Summary_History_win_width
      b% Summary_History_win_aspect_ratio = Summary_History_win_aspect_ratio
      b% Summary_History_file_width = Summary_History_file_width
      b% Summary_History_file_aspect_ratio = Summary_History_file_aspect_ratio
      b% Summary_History_xleft = Summary_History_xleft
      b% Summary_History_xright = Summary_History_xright
      b% Summary_History_ybot = Summary_History_ybot
      b% Summary_History_ytop = Summary_History_ytop
      b% Summary_History_txt_scale = Summary_History_txt_scale
      b% Summary_History_title = Summary_History_title
      b% Summary_History_name = Summary_History_name
      b% Summary_History_legend = Summary_History_legend
      b% Summary_History_num_lines = Summary_History_num_lines
      b% Summary_History_use_decorator = Summary_History_use_decorator

      b% Grid1_win_flag = Grid1_win_flag
      b% Grid1_win_width = Grid1_win_width
      b% Grid1_win_aspect_ratio = Grid1_win_aspect_ratio
      b% Grid1_xleft = Grid1_xleft
      b% Grid1_xright = Grid1_xright
      b% Grid1_ybot = Grid1_ybot
      b% Grid1_ytop = Grid1_ytop
      b% Grid1_title = Grid1_title
      b% Grid1_txt_scale_factor = Grid1_txt_scale_factor
      b% Grid1_num_cols = Grid1_num_cols
      b% Grid1_num_rows = Grid1_num_rows
      b% Grid1_num_plots = Grid1_num_plots
      b% Grid1_plot_name = Grid1_plot_name
      b% Grid1_plot_row = Grid1_plot_row
      b% Grid1_plot_rowspan = Grid1_plot_rowspan
      b% Grid1_plot_col = Grid1_plot_col
      b% Grid1_plot_colspan = Grid1_plot_colspan
      b% Grid1_plot_pad_left = Grid1_plot_pad_left
      b% Grid1_plot_pad_right = Grid1_plot_pad_right
      b% Grid1_plot_pad_top = Grid1_plot_pad_top
      b% Grid1_plot_pad_bot = Grid1_plot_pad_bot
      b% Grid1_file_flag = Grid1_file_flag
      b% Grid1_file_dir = Grid1_file_dir
      b% Grid1_file_prefix = Grid1_file_prefix
      b% Grid1_file_interval = Grid1_file_interval
      b% Grid1_file_width = Grid1_file_width
      b% Grid1_file_aspect_ratio = Grid1_file_aspect_ratio

      b% Grid2_win_flag = Grid2_win_flag
      b% Grid2_win_width = Grid2_win_width
      b% Grid2_win_aspect_ratio = Grid2_win_aspect_ratio
      b% Grid2_xleft = Grid2_xleft
      b% Grid2_xright = Grid2_xright
      b% Grid2_ybot = Grid2_ybot
      b% Grid2_ytop = Grid2_ytop
      b% Grid2_title = Grid2_title
      b% Grid2_txt_scale_factor = Grid2_txt_scale_factor
      b% Grid2_num_cols = Grid2_num_cols
      b% Grid2_num_rows = Grid2_num_rows
      b% Grid2_num_plots = Grid2_num_plots
      b% Grid2_plot_name = Grid2_plot_name
      b% Grid2_plot_row = Grid2_plot_row
      b% Grid2_plot_rowspan = Grid2_plot_rowspan
      b% Grid2_plot_col = Grid2_plot_col
      b% Grid2_plot_colspan = Grid2_plot_colspan
      b% Grid2_plot_pad_left = Grid2_plot_pad_left
      b% Grid2_plot_pad_right = Grid2_plot_pad_right
      b% Grid2_plot_pad_top = Grid2_plot_pad_top
      b% Grid2_plot_pad_bot = Grid2_plot_pad_bot
      b% Grid2_file_flag = Grid2_file_flag
      b% Grid2_file_dir = Grid2_file_dir
      b% Grid2_file_prefix = Grid2_file_prefix
      b% Grid2_file_interval = Grid2_file_interval
      b% Grid2_file_width = Grid2_file_width
      b% Grid2_file_aspect_ratio = Grid2_file_aspect_ratio

      b% Grid3_win_flag = Grid3_win_flag
      b% Grid3_win_width = Grid3_win_width
      b% Grid3_win_aspect_ratio = Grid3_win_aspect_ratio
      b% Grid3_xleft = Grid3_xleft
      b% Grid3_xright = Grid3_xright
      b% Grid3_ybot = Grid3_ybot
      b% Grid3_ytop = Grid3_ytop
      b% Grid3_title = Grid3_title
      b% Grid3_txt_scale_factor = Grid3_txt_scale_factor
      b% Grid3_num_cols = Grid3_num_cols
      b% Grid3_num_rows = Grid3_num_rows
      b% Grid3_num_plots = Grid3_num_plots
      b% Grid3_plot_name = Grid3_plot_name
      b% Grid3_plot_row = Grid3_plot_row
      b% Grid3_plot_rowspan = Grid3_plot_rowspan
      b% Grid3_plot_col = Grid3_plot_col
      b% Grid3_plot_colspan = Grid3_plot_colspan
      b% Grid3_plot_pad_left = Grid3_plot_pad_left
      b% Grid3_plot_pad_right = Grid3_plot_pad_right
      b% Grid3_plot_pad_top = Grid3_plot_pad_top
      b% Grid3_plot_pad_bot = Grid3_plot_pad_bot
      b% Grid3_file_flag = Grid3_file_flag
      b% Grid3_file_dir = Grid3_file_dir
      b% Grid3_file_prefix = Grid3_file_prefix
      b% Grid3_file_interval = Grid3_file_interval
      b% Grid3_file_width = Grid3_file_width
      b% Grid3_file_aspect_ratio = Grid3_file_aspect_ratio

      b% Grid4_win_flag = Grid4_win_flag
      b% Grid4_win_width = Grid4_win_width
      b% Grid4_win_aspect_ratio = Grid4_win_aspect_ratio
      b% Grid4_xleft = Grid4_xleft
      b% Grid4_xright = Grid4_xright
      b% Grid4_ybot = Grid4_ybot
      b% Grid4_ytop = Grid4_ytop
      b% Grid4_title = Grid4_title
      b% Grid4_txt_scale_factor = Grid4_txt_scale_factor
      b% Grid4_num_cols = Grid4_num_cols
      b% Grid4_num_rows = Grid4_num_rows
      b% Grid4_num_plots = Grid4_num_plots
      b% Grid4_plot_name = Grid4_plot_name
      b% Grid4_plot_row = Grid4_plot_row
      b% Grid4_plot_rowspan = Grid4_plot_rowspan
      b% Grid4_plot_col = Grid4_plot_col
      b% Grid4_plot_colspan = Grid4_plot_colspan
      b% Grid4_plot_pad_left = Grid4_plot_pad_left
      b% Grid4_plot_pad_right = Grid4_plot_pad_right
      b% Grid4_plot_pad_top = Grid4_plot_pad_top
      b% Grid4_plot_pad_bot = Grid4_plot_pad_bot
      b% Grid4_file_flag = Grid4_file_flag
      b% Grid4_file_dir = Grid4_file_dir
      b% Grid4_file_prefix = Grid4_file_prefix
      b% Grid4_file_interval = Grid4_file_interval
      b% Grid4_file_width = Grid4_file_width
      b% Grid4_file_aspect_ratio = Grid4_file_aspect_ratio

      b% Grid5_win_flag = Grid5_win_flag
      b% Grid5_win_width = Grid5_win_width
      b% Grid5_win_aspect_ratio = Grid5_win_aspect_ratio
      b% Grid5_xleft = Grid5_xleft
      b% Grid5_xright = Grid5_xright
      b% Grid5_ybot = Grid5_ybot
      b% Grid5_ytop = Grid5_ytop
      b% Grid5_title = Grid5_title
      b% Grid5_txt_scale_factor = Grid5_txt_scale_factor
      b% Grid5_num_cols = Grid5_num_cols
      b% Grid5_num_rows = Grid5_num_rows
      b% Grid5_num_plots = Grid5_num_plots
      b% Grid5_plot_name = Grid5_plot_name
      b% Grid5_plot_row = Grid5_plot_row
      b% Grid5_plot_rowspan = Grid5_plot_rowspan
      b% Grid5_plot_col = Grid5_plot_col
      b% Grid5_plot_colspan = Grid5_plot_colspan
      b% Grid5_plot_pad_left = Grid5_plot_pad_left
      b% Grid5_plot_pad_right = Grid5_plot_pad_right
      b% Grid5_plot_pad_top = Grid5_plot_pad_top
      b% Grid5_plot_pad_bot = Grid5_plot_pad_bot
      b% Grid5_file_flag = Grid5_file_flag
      b% Grid5_file_dir = Grid5_file_dir
      b% Grid5_file_prefix = Grid5_file_prefix
      b% Grid5_file_interval = Grid5_file_interval
      b% Grid5_file_width = Grid5_file_width
      b% Grid5_file_aspect_ratio = Grid5_file_aspect_ratio

      b% Grid6_win_flag = Grid6_win_flag
      b% Grid6_win_width = Grid6_win_width
      b% Grid6_win_aspect_ratio = Grid6_win_aspect_ratio
      b% Grid6_xleft = Grid6_xleft
      b% Grid6_xright = Grid6_xright
      b% Grid6_ybot = Grid6_ybot
      b% Grid6_ytop = Grid6_ytop
      b% Grid6_title = Grid6_title
      b% Grid6_txt_scale_factor = Grid6_txt_scale_factor
      b% Grid6_num_cols = Grid6_num_cols
      b% Grid6_num_rows = Grid6_num_rows
      b% Grid6_num_plots = Grid6_num_plots
      b% Grid6_plot_name = Grid6_plot_name
      b% Grid6_plot_row = Grid6_plot_row
      b% Grid6_plot_rowspan = Grid6_plot_rowspan
      b% Grid6_plot_col = Grid6_plot_col
      b% Grid6_plot_colspan = Grid6_plot_colspan
      b% Grid6_plot_pad_left = Grid6_plot_pad_left
      b% Grid6_plot_pad_right = Grid6_plot_pad_right
      b% Grid6_plot_pad_top = Grid6_plot_pad_top
      b% Grid6_plot_pad_bot = Grid6_plot_pad_bot
      b% Grid6_file_flag = Grid6_file_flag
      b% Grid6_file_dir = Grid6_file_dir
      b% Grid6_file_prefix = Grid6_file_prefix
      b% Grid6_file_interval = Grid6_file_interval
      b% Grid6_file_width = Grid6_file_width
      b% Grid6_file_aspect_ratio = Grid6_file_aspect_ratio

      b% Grid7_win_flag = Grid7_win_flag
      b% Grid7_win_width = Grid7_win_width
      b% Grid7_win_aspect_ratio = Grid7_win_aspect_ratio
      b% Grid7_xleft = Grid7_xleft
      b% Grid7_xright = Grid7_xright
      b% Grid7_ybot = Grid7_ybot
      b% Grid7_ytop = Grid7_ytop
      b% Grid7_title = Grid7_title
      b% Grid7_txt_scale_factor = Grid7_txt_scale_factor
      b% Grid7_num_cols = Grid7_num_cols
      b% Grid7_num_rows = Grid7_num_rows
      b% Grid7_num_plots = Grid7_num_plots
      b% Grid7_plot_name = Grid7_plot_name
      b% Grid7_plot_row = Grid7_plot_row
      b% Grid7_plot_rowspan = Grid7_plot_rowspan
      b% Grid7_plot_col = Grid7_plot_col
      b% Grid7_plot_colspan = Grid7_plot_colspan
      b% Grid7_plot_pad_left = Grid7_plot_pad_left
      b% Grid7_plot_pad_right = Grid7_plot_pad_right
      b% Grid7_plot_pad_top = Grid7_plot_pad_top
      b% Grid7_plot_pad_bot = Grid7_plot_pad_bot
      b% Grid7_file_flag = Grid7_file_flag
      b% Grid7_file_dir = Grid7_file_dir
      b% Grid7_file_prefix = Grid7_file_prefix
      b% Grid7_file_interval = Grid7_file_interval
      b% Grid7_file_width = Grid7_file_width
      b% Grid7_file_aspect_ratio = Grid7_file_aspect_ratio

      b% Grid8_win_flag = Grid8_win_flag
      b% Grid8_win_width = Grid8_win_width
      b% Grid8_win_aspect_ratio = Grid8_win_aspect_ratio
      b% Grid8_xleft = Grid8_xleft
      b% Grid8_xright = Grid8_xright
      b% Grid8_ybot = Grid8_ybot
      b% Grid8_ytop = Grid8_ytop
      b% Grid8_title = Grid8_title
      b% Grid8_txt_scale_factor = Grid8_txt_scale_factor
      b% Grid8_num_cols = Grid8_num_cols
      b% Grid8_num_rows = Grid8_num_rows
      b% Grid8_num_plots = Grid8_num_plots
      b% Grid8_plot_name = Grid8_plot_name
      b% Grid8_plot_row = Grid8_plot_row
      b% Grid8_plot_rowspan = Grid8_plot_rowspan
      b% Grid8_plot_col = Grid8_plot_col
      b% Grid8_plot_colspan = Grid8_plot_colspan
      b% Grid8_plot_pad_left = Grid8_plot_pad_left
      b% Grid8_plot_pad_right = Grid8_plot_pad_right
      b% Grid8_plot_pad_top = Grid8_plot_pad_top
      b% Grid8_plot_pad_bot = Grid8_plot_pad_bot
      b% Grid8_file_flag = Grid8_file_flag
      b% Grid8_file_dir = Grid8_file_dir
      b% Grid8_file_prefix = Grid8_file_prefix
      b% Grid8_file_interval = Grid8_file_interval
      b% Grid8_file_width = Grid8_file_width
      b% Grid8_file_aspect_ratio = Grid8_file_aspect_ratio

      b% Grid9_win_flag = Grid9_win_flag
      b% Grid9_win_width = Grid9_win_width
      b% Grid9_win_aspect_ratio = Grid9_win_aspect_ratio
      b% Grid9_xleft = Grid9_xleft
      b% Grid9_xright = Grid9_xright
      b% Grid9_ybot = Grid9_ybot
      b% Grid9_ytop = Grid9_ytop
      b% Grid9_title = Grid9_title
      b% Grid9_txt_scale_factor = Grid9_txt_scale_factor
      b% Grid9_num_cols = Grid9_num_cols
      b% Grid9_num_rows = Grid9_num_rows
      b% Grid9_num_plots = Grid9_num_plots
      b% Grid9_plot_name = Grid9_plot_name
      b% Grid9_plot_row = Grid9_plot_row
      b% Grid9_plot_rowspan = Grid9_plot_rowspan
      b% Grid9_plot_col = Grid9_plot_col
      b% Grid9_plot_colspan = Grid9_plot_colspan
      b% Grid9_plot_pad_left = Grid9_plot_pad_left
      b% Grid9_plot_pad_right = Grid9_plot_pad_right
      b% Grid9_plot_pad_top = Grid9_plot_pad_top
      b% Grid9_plot_pad_bot = Grid9_plot_pad_bot
      b% Grid9_file_flag = Grid9_file_flag
      b% Grid9_file_dir = Grid9_file_dir
      b% Grid9_file_prefix = Grid9_file_prefix
      b% Grid9_file_interval = Grid9_file_interval
      b% Grid9_file_width = Grid9_file_width
      b% Grid9_file_aspect_ratio = Grid9_file_aspect_ratio
      
      b% Star1_win_flag = Star1_win_flag
      b% Star1_file_flag = Star1_file_flag
      b% Star1_file_interval = Star1_file_interval
      b% Star1_file_dir = Star1_file_dir
      b% Star1_file_prefix = Star1_file_prefix
      b% Star1_win_width = Star1_win_width
      b% Star1_win_aspect_ratio = Star1_win_aspect_ratio
      b% Star1_xleft = Star1_xleft
      b% Star1_xright = Star1_xright
      b% Star1_ybot = Star1_ybot
      b% Star1_ytop = Star1_ytop
      b% Star1_file_width = Star1_file_width
      b% Star1_file_aspect_ratio = Star1_file_aspect_ratio
      b% Star1_txt_scale_factor = Star1_txt_scale_factor
      b% Star1_title = Star1_title
      b% Star1_plot_name = Star1_plot_name
      b% do_star1_box = do_star1_box
      b% star1_box_pad_left = star1_box_pad_left
      b% star1_box_pad_right = star1_box_pad_right
      b% star1_box_pad_bot = star1_box_pad_bot
      b% star1_box_pad_top = star1_box_pad_top

      b% Star2_win_flag = Star2_win_flag
      b% Star2_file_flag = Star2_file_flag
      b% Star2_file_interval = Star2_file_interval
      b% Star2_file_dir = Star2_file_dir
      b% Star2_file_prefix = Star2_file_prefix
      b% Star2_win_width = Star2_win_width
      b% Star2_win_aspect_ratio = Star2_win_aspect_ratio
      b% Star2_xleft = Star2_xleft
      b% Star2_xright = Star2_xright
      b% Star2_ybot = Star2_ybot
      b% Star2_ytop = Star2_ytop
      b% Star2_file_width = Star2_file_width
      b% Star2_file_aspect_ratio = Star2_file_aspect_ratio
      b% Star2_txt_scale_factor = Star2_txt_scale_factor
      b% Star2_title = Star2_title
      b% Star2_plot_name = Star2_plot_name
      b% do_star2_box = do_star2_box
      b% star2_box_pad_left = star2_box_pad_left
      b% star2_box_pad_right = star2_box_pad_right
      b% star2_box_pad_bot = star2_box_pad_bot
      b% star2_box_pad_top = star2_box_pad_top

      b% annotation1_ci = annotation1_ci
      b% annotation1_ch = annotation1_ch
      b% annotation1_lw = annotation1_lw
      b% annotation1_cf = annotation1_cf
      b% annotation1_text = annotation1_text
      b% annotation1_side = annotation1_side
      b% annotation1_disp = annotation1_disp
      b% annotation1_coord = annotation1_coord
      b% annotation1_fjust = annotation1_fjust

      b% annotation2_ci = annotation2_ci
      b% annotation2_ch = annotation2_ch
      b% annotation2_lw = annotation2_lw
      b% annotation2_cf = annotation2_cf
      b% annotation2_text = annotation2_text
      b% annotation2_side = annotation2_side
      b% annotation2_disp = annotation2_disp
      b% annotation2_coord = annotation2_coord
      b% annotation2_fjust = annotation2_fjust

      b% annotation3_ci = annotation3_ci
      b% annotation3_ch = annotation3_ch
      b% annotation3_lw = annotation3_lw
      b% annotation3_cf = annotation3_cf
      b% annotation3_text = annotation3_text
      b% annotation3_side = annotation3_side
      b% annotation3_disp = annotation3_disp
      b% annotation3_coord = annotation3_coord
      b% annotation3_fjust = annotation3_fjust

      b% read_extra_pgbinary_inlist1 = read_extra_pgbinary_inlist1
      b% extra_pgbinary_inlist1_name = extra_pgbinary_inlist1_name

      b% read_extra_pgbinary_inlist2 = read_extra_pgbinary_inlist2
      b% extra_pgbinary_inlist2_name = extra_pgbinary_inlist2_name

      b% read_extra_pgbinary_inlist3 = read_extra_pgbinary_inlist3
      b% extra_pgbinary_inlist3_name = extra_pgbinary_inlist3_name

      b% read_extra_pgbinary_inlist4 = read_extra_pgbinary_inlist4
      b% extra_pgbinary_inlist4_name = extra_pgbinary_inlist4_name

      b% read_extra_pgbinary_inlist5 = read_extra_pgbinary_inlist5
      b% extra_pgbinary_inlist5_name = extra_pgbinary_inlist5_name

   end subroutine store_pgbinary_controls


   subroutine set_default_pgbinary_controls

      include 'pgbinary.defaults'

   end subroutine set_default_pgbinary_controls


end module pgbinary_ctrls_io

