! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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

      module pgstar_ctrls_io


      use const_def
      use star_private_def
      use star_pgstar

      implicit none

      include "pgstar_controls.inc"

      namelist /pgstar/ &

            file_device, &
            file_extension, &
            file_digits, &
            pgstar_interval, &
            pause, &
            pause_interval, &
            pgstar_sleep, &
            clear_history, &

            delta_HR_limit_for_file_output, &
            file_white_on_black_flag, &
            win_white_on_black_flag, &
            pgstar_show_model_number, &
            pgstar_show_age, &
            pgstar_show_age_in_seconds, &
            pgstar_show_age_in_minutes, &
            pgstar_show_age_in_hours, &
            pgstar_show_age_in_days, &
            pgstar_show_age_in_years, &
            pgstar_show_log_age_in_years, &

            pgstar_report_writing_files, &

            pgstar_show_title, &
            pgstar_title_scale, &
            pgstar_title_disp, &
            pgstar_title_coord, &
            pgstar_title_fjust, &
            pgstar_title_lw, &

            pgstar_grid_show_title, &
            pgstar_grid_title_scale, &
            pgstar_grid_title_disp, &
            pgstar_grid_title_coord, &
            pgstar_grid_title_fjust, &
            pgstar_grid_title_lw, &

            pgstar_age_scale, &
            pgstar_age_disp, &
            pgstar_age_coord, &
            pgstar_age_fjust, &
            pgstar_age_lw, &

            pgstar_model_scale, &
            pgstar_model_disp, &
            pgstar_model_coord, &
            pgstar_model_fjust, &
            pgstar_model_lw, &

            pgstar_xaxis_label_scale, &
            pgstar_left_yaxis_label_scale, &
            pgstar_right_yaxis_label_scale, &
            pgstar_xaxis_label_lw, &
            pgstar_left_yaxis_label_lw, &
            pgstar_right_yaxis_label_lw, &
            pgstar_xaxis_label_disp, &
            pgstar_left_yaxis_label_disp, &
            pgstar_right_yaxis_label_disp, &
            pgstar_num_scale, &
            pgstar_lw, &
            pgstar_profile_line_style, &
            pgstar_history_line_style, &
            pgstar_model_lw, &
            pgstar_box_lw, &

            Profile_Panels_show_Mach_1_location, &
            Profile_Panels_show_photosphere_location, &
            Profile_Panels_xwidth_left_div_shock_value, &
            Profile_Panels_xwidth_right_div_shock_value, &
            Profile_Panels_xwidth_left_of_shock, &
            Profile_Panels_xwidth_right_of_shock, &

            Profile_Panels_win_flag, &
            Profile_Panels_file_flag, &
            Profile_Panels_file_interval, &
            Profile_Panels_file_dir, &
            Profile_Panels_file_prefix, &
            Profile_Panels_xaxis_reversed, &
            Profile_Panels_xaxis_name, &
            Profile_Panels_title, &
            Profile_Panels_xmin, &
            Profile_Panels_xmax, &
            Profile_Panels_xmargin, &
            Profile_Panels_show_mix_regions_on_xaxis, &
            Profile_Panels_win_width, &
            Profile_Panels_win_aspect_ratio, &
            Profile_Panels_xleft, &
            Profile_Panels_xright, &
            Profile_Panels_ybot, &
            Profile_Panels_ytop, &
            Profile_Panels_txt_scale, &
            prev_Profile_Panels_win_width, &
            prev_Profile_Panels_win_ratio, &
            Profile_Panels_file_width, &
            Profile_Panels_file_aspect_ratio, &
            prev_Profile_Panels_file_width, &
            prev_Profile_Panels_file_ratio, &
            Profile_Panels_num_panels, &
            Profile_Panels_yaxis_name, &
            Profile_Panels_other_yaxis_name, &
            Profile_Panels_yaxis_reversed, &
            Profile_Panels_other_yaxis_reversed, &
            Profile_Panels_yaxis_log, &
            Profile_Panels_other_yaxis_log, &
            Profile_Panels_same_yaxis_range, &
            Profile_Panels_ymin, &
            Profile_Panels_other_ymin, &
            Profile_Panels_ymax, &
            Profile_Panels_other_ymax, &
            Profile_Panels_ycenter, &
            Profile_Panels_other_ycenter, &
            Profile_Panels_ymargin, &
            Profile_Panels_other_ymargin, &
            Profile_Panels_dymin, &
            Profile_Panels_other_dymin, &
            Profile_Panels_show_grid, &
            Profile_Panels_use_decorator, &


            Text_Summary_win_flag, &
            Text_Summary_file_flag, &
            Text_Summary_file_interval, &
            Text_Summary_file_dir, &
            Text_Summary_file_prefix, &
            Text_Summary_num_cols, &
            Text_Summary_num_rows, &
            Text_Summary_name, &
            Text_Summary_win_width, &
            Text_Summary_win_aspect_ratio, &
            Text_Summary_file_width, &
            Text_Summary_file_aspect_ratio, &
            Text_Summary_title, &
            Text_Summary_xleft, &
            Text_Summary_xright, &
            Text_Summary_ybot, &
            Text_Summary_ytop, &
            Text_Summary_txt_scale, &
            Text_Summary_dxval, &

            logg_Teff_win_flag, &
            logg_Teff_file_flag, &
            show_logg_Teff_target_box, &
            logg_Teff_target_n_sigma, &
            logg_Teff_target_logg, &
            logg_Teff_target_logg_sigma, &
            logg_Teff_target_Teff, &
            logg_Teff_target_Teff_sigma, &
            logg_Teff_file_interval, &
            logg_Teff_step_min, &
            logg_Teff_step_max, &
            logg_Teff_file_dir, &
            logg_Teff_file_prefix, &
            show_logg_Teff_annotation1, &
            show_logg_Teff_annotation2, &
            show_logg_Teff_annotation3, &
            logg_Teff_fname, &
            logg_Teff_title, &
            logg_Teff_logg_min, &
            logg_Teff_logg_max, &
            logg_Teff_Teff_min, &
            logg_Teff_Teff_max, &
            logg_Teff_Teff_margin, &
            logg_Teff_logg_margin, &
            logg_Teff_dTeff_min, &
            logg_Teff_dlogg_min, &
            logg_Teff_win_width, &
            logg_Teff_xleft, &
            logg_Teff_xright, &
            logg_Teff_ybot, &
            logg_Teff_ytop, &
            logg_Teff_txt_scale, &
            logg_Teff_win_aspect_ratio, &
            logg_Teff_file_width, &
            logg_Teff_file_aspect_ratio, &
            logg_Teff_use_decorator, &
            
            logL_Teff_win_flag, &
            logL_Teff_file_flag, &
            show_logL_Teff_target_box, &
            logL_Teff_target_n_sigma, &
            logL_Teff_target_logL, &
            logL_Teff_target_logL_sigma, &
            logL_Teff_target_Teff, &
            logL_Teff_target_Teff_sigma, &
            logL_Teff_file_interval, &
            logL_Teff_step_min, &
            logL_Teff_step_max, &
            logL_Teff_file_dir, &
            logL_Teff_file_prefix, &
            show_logL_Teff_annotation1, &
            show_logL_Teff_annotation2, &
            show_logL_Teff_annotation3, &
            logL_Teff_fname, &
            logL_Teff_title, &
            logL_Teff_logL_min, &
            logL_Teff_logL_max, &
            logL_Teff_Teff_min, &
            logL_Teff_Teff_max, &
            logL_Teff_Teff_margin, &
            logL_Teff_logL_margin, &
            logL_Teff_dTeff_min, &
            logL_Teff_dlogL_min, &
            logL_Teff_win_width, &
            logL_Teff_xleft, &
            logL_Teff_xright, &
            logL_Teff_ybot, &
            logL_Teff_ytop, &
            logL_Teff_txt_scale, &
            logL_Teff_win_aspect_ratio, &
            logL_Teff_file_width, &
            logL_Teff_file_aspect_ratio, &
            logL_Teff_use_decorator, &
            
            logL_R_win_flag, &
            logL_R_file_flag, &
            show_logL_R_target_box, &
            logL_R_target_n_sigma, &
            logL_R_target_logL, &
            logL_R_target_logL_sigma, &
            logL_R_target_R, &
            logL_R_target_R_sigma, &
            logL_R_file_interval, &
            logL_R_step_min, &
            logL_R_step_max, &
            logL_R_file_dir, &
            logL_R_file_prefix, &
            show_logL_R_annotation1, &
            show_logL_R_annotation2, &
            show_logL_R_annotation3, &
            logL_R_fname, &
            logL_R_title, &
            show_logL_photosphere_r, &
            logL_R_logL_min, &
            logL_R_logL_max, &
            logL_R_R_min, &
            logL_R_R_max, &
            logL_R_R_margin, &
            logL_R_logL_margin, &
            logL_R_dR_min, &
            logL_R_dlogL_min, &
            logL_R_win_width, &
            logL_R_xleft, &
            logL_R_xright, &
            logL_R_ybot, &
            logL_R_ytop, &
            logL_R_txt_scale, &
            logL_R_win_aspect_ratio, &
            logL_R_file_width, &
            logL_R_file_aspect_ratio, &
            logL_R_use_decorator, &

            logL_v_win_flag, &
            logL_v_file_flag, &
            show_logL_v_target_box, &
            logL_v_target_n_sigma, &
            logL_v_target_logL, &
            logL_v_target_logL_sigma, &
            logL_v_target_v, &
            logL_v_target_v_sigma, &
            logL_v_file_interval, &
            logL_v_step_min, &
            logL_v_step_max, &
            logL_v_file_dir, &
            logL_v_file_prefix, &
            show_logL_v_annotation1, &
            show_logL_v_annotation2, &
            show_logL_v_annotation3, &
            logL_v_fname, &
            logL_v_title, &
            show_logL_photosphere_v, &
            logL_v_logL_min, &
            logL_v_logL_max, &
            logL_v_v_min, &
            logL_v_v_max, &
            logL_v_v_margin, &
            logL_v_logL_margin, &
            logL_v_dv_min, &
            logL_v_dlogL_min, &
            logL_v_win_width, &
            logL_v_xleft, &
            logL_v_xright, &
            logL_v_ybot, &
            logL_v_ytop, &
            logL_v_txt_scale, &
            logL_v_win_aspect_ratio, &
            logL_v_file_width, &
            logL_v_file_aspect_ratio, &
            logL_v_use_decorator, &

            L_Teff_win_flag, &
            L_Teff_file_flag, &
            show_L_Teff_target_box, &
            L_Teff_target_n_sigma, &
            L_Teff_target_L, &
            L_Teff_target_L_sigma, &
            L_Teff_target_Teff, &
            L_Teff_target_Teff_sigma, &
            L_Teff_file_interval, &
            L_Teff_step_min, &
            L_Teff_step_max, &
            L_Teff_file_dir, &
            L_Teff_file_prefix, &
            show_L_Teff_annotation1, &
            show_L_Teff_annotation2, &
            show_L_Teff_annotation3, &
            L_Teff_fname, &
            L_Teff_title, &
            L_Teff_L_min, &
            L_Teff_L_max, &
            L_Teff_Teff_min, &
            L_Teff_Teff_max, &
            L_Teff_Teff_margin, &
            L_Teff_L_margin, &
            L_Teff_dTeff_min, &
            L_Teff_dL_min, &
            L_Teff_win_width, &
            L_Teff_xleft, &
            L_Teff_xright, &
            L_Teff_ybot, &
            L_Teff_ytop, &
            L_Teff_txt_scale, &
            L_Teff_win_aspect_ratio, &
            L_Teff_file_width, &
            L_Teff_file_aspect_ratio, &
            L_Teff_use_decorator, &

            L_v_win_flag, &
            L_v_file_flag, &
            show_L_v_target_box, &
            L_v_target_n_sigma, &
            L_v_target_L, &
            L_v_target_L_sigma, &
            L_v_target_v, &
            L_v_target_v_sigma, &
            L_v_file_interval, &
            L_v_step_min, &
            L_v_step_max, &
            L_v_file_dir, &
            L_v_file_prefix, &
            show_L_v_annotation1, &
            show_L_v_annotation2, &
            show_L_v_annotation3, &
            L_v_fname, &
            L_v_title, &
            L_v_L_min, &
            L_v_L_max, &
            L_v_v_min, &
            L_v_v_max, &
            L_v_v_margin, &
            L_v_L_margin, &
            L_v_dv_min, &
            L_v_dL_min, &
            L_v_win_width, &
            L_v_xleft, &
            L_v_xright, &
            L_v_ybot, &
            L_v_ytop, &
            L_v_txt_scale, &
            L_v_win_aspect_ratio, &
            L_v_file_width, &
            L_v_file_aspect_ratio, &
            L_V_use_decorator, &

            L_R_win_flag, &
            L_R_file_flag, &
            show_L_R_target_box, &
            L_R_target_n_sigma, &
            L_R_target_L, &
            L_R_target_L_sigma, &
            L_R_target_R, &
            L_R_target_R_sigma, &
            L_R_file_interval, &
            L_R_step_min, &
            L_R_step_max, &
            L_R_file_dir, &
            L_R_file_prefix, &
            show_L_R_annotation1, &
            show_L_R_annotation2, &
            show_L_R_annotation3, &
            L_R_fname, &
            L_R_title, &
            L_R_L_min, &
            L_R_L_max, &
            L_R_R_min, &
            L_R_R_max, &
            L_R_R_margin, &
            L_R_L_margin, &
            L_R_dR_min, &
            L_R_dL_min, &
            L_R_win_width, &
            L_R_xleft, &
            L_R_xright, &
            L_R_ybot, &
            L_R_ytop, &
            L_R_txt_scale, &
            L_R_win_aspect_ratio, &
            L_R_file_width, &
            L_R_file_aspect_ratio, &
            L_R_use_decorator, &

            R_Teff_win_flag, &
            R_Teff_file_flag, &
            show_R_Teff_target_box, &
            R_Teff_target_n_sigma, &
            R_Teff_target_R, &
            R_Teff_target_R_sigma, &
            R_Teff_target_Teff, &
            R_Teff_target_Teff_sigma, &
            R_Teff_file_interval, &
            R_Teff_step_min, &
            R_Teff_step_max, &
            R_Teff_file_dir, &
            R_Teff_file_prefix, &
            show_R_Teff_annotation1, &
            show_R_Teff_annotation2, &
            show_R_Teff_annotation3, &
            R_Teff_fname, &
            R_Teff_title, &
            R_Teff_R_min, &
            R_Teff_R_max, &
            R_Teff_Teff_min, &
            R_Teff_Teff_max, &
            R_Teff_Teff_margin, &
            R_Teff_R_margin, &
            R_Teff_dTeff_min, &
            R_Teff_dR_min, &
            R_Teff_win_width, &
            R_Teff_xleft, &
            R_Teff_xright, &
            R_Teff_ybot, &
            R_Teff_ytop, &
            R_Teff_txt_scale, &
            R_Teff_win_aspect_ratio, &
            R_Teff_file_width, &
            R_Teff_file_aspect_ratio, &
            R_Teff_use_decorator, &

            R_L_win_flag, &
            R_L_file_flag, &
            show_R_L_target_box, &
            R_L_target_n_sigma, &
            R_L_target_R, &
            R_L_target_R_sigma, &
            R_L_target_L, &
            R_L_target_L_sigma, &
            R_L_file_interval, &
            R_L_step_min, &
            R_L_step_max, &
            R_L_file_dir, &
            R_L_file_prefix, &
            show_R_L_annotation1, &
            show_R_L_annotation2, &
            show_R_L_annotation3, &
            R_L_fname, &
            R_L_title, &
            R_L_R_min, &
            R_L_R_max, &
            R_L_L_min, &
            R_L_L_max, &
            R_L_L_margin, &
            R_L_R_margin, &
            R_L_dL_min, &
            R_L_dR_min, &
            R_L_win_width, &
            R_L_xleft, &
            R_L_xright, &
            R_L_ybot, &
            R_L_ytop, &
            R_L_txt_scale, &
            R_L_win_aspect_ratio, &
            R_L_file_width, &
            R_L_file_aspect_ratio, &
            R_L_use_decorator, &

            logg_logT_win_flag, &
            logg_logT_file_flag, &
            show_logg_logT_target_box, &
            logg_logT_target_n_sigma, &
            logg_logT_target_logg, &
            logg_logT_target_logg_sigma, &
            logg_logT_target_logT, &
            logg_logT_target_logT_sigma, &
            logg_logT_file_interval, &
            logg_logT_step_min, &
            logg_logT_step_max, &
            logg_logT_file_dir, &
            logg_logT_file_prefix, &
            show_logg_logT_annotation1, &
            show_logg_logT_annotation2, &
            show_logg_logT_annotation3, &
            logg_logT_fname, &
            logg_logT_logg_min, &
            logg_logT_logg_max, &
            logg_logT_logT_min, &
            logg_logT_logT_max, &
            logg_logT_logg_margin, &
            logg_logT_logT_margin, &
            logg_logT_dlogT_min, &
            logg_logT_dlogg_min, &
            logg_logT_win_width, &
            logg_logT_win_aspect_ratio, &
            logg_logT_file_width, &
            logg_logT_file_aspect_ratio, &
            logg_logT_xleft, &
            logg_logT_xright, &
            logg_logT_ybot, &
            logg_logT_ytop, &
            logg_logT_txt_scale, &
            logg_logT_title, &
            logg_logT_use_decorator, &

            HR_win_flag, &
            HR_file_flag, &
            HR_file_interval, &
            HR_step_min, &
            HR_step_max, &
            show_HR_classical_instability_strip, &
            show_HR_Mira_instability_region, &
            show_HR_WD_instabilities, &
            show_HR_target_box, &
            HR_target_n_sigma, &
            HR_target_logL, &
            HR_target_logL_sigma, &
            HR_target_logT, &
            HR_target_logT_sigma, &
            HR_file_dir, &
            HR_file_prefix, &
            show_HR_annotation1, &
            show_HR_annotation2, &
            show_HR_annotation3, &
            HR_fname, &
            HR_logT_min, &
            HR_logT_max, &
            HR_logL_min, &
            HR_logL_max, &
            HR_logL_margin, &
            HR_logT_margin, &
            HR_dlogT_min, &
            HR_dlogL_min, &
            HR_win_width, &
            HR_win_aspect_ratio, &
            HR_file_width, &
            HR_file_aspect_ratio, &
            HR_xleft, &
            HR_xright, &
            HR_ybot, &
            HR_ytop, &
            HR_txt_scale, &
            HR_title, &
            HR_use_decorator, &

            TRho_win_flag, &
            TRho_file_flag, &
            TRho_file_interval, &
            TRho_step_max, &
            TRho_step_min, &
            TRho_file_dir, &
            TRho_file_prefix, &
            show_TRho_annotation1, &
            show_TRho_annotation2, &
            show_TRho_annotation3, &
            show_TRho_degeneracy_line, &
            TRho_fname, &
            TRho_logT_min, &
            TRho_logT_max, &
            TRho_logRho_min, &
            TRho_logRho_max, &
            TRho_logT_margin, &
            TRho_logRho_margin, &
            TRho_logRho_dlogRho_min, &
            TRho_logT_dlogT_min, &
            TRho_win_width, &
            TRho_win_aspect_ratio, &
            TRho_file_width, &
            TRho_file_aspect_ratio, &
            TRho_xleft, &
            TRho_xright, &
            TRho_ybot, &
            TRho_ytop, &
            TRho_txt_scale, &
            TRho_title, &
            TRho_use_decorator, &

            TmaxRho_win_flag, &
            TmaxRho_file_flag, &
            TmaxRho_file_interval, &
            TmaxRho_step_max, &
            TmaxRho_step_min, &
            TmaxRho_file_dir, &
            TmaxRho_file_prefix, &
            show_TmaxRho_annotation1, &
            show_TmaxRho_annotation2, &
            show_TmaxRho_annotation3, &
            show_TmaxRho_degeneracy_line, &
            TmaxRho_fname, &
            TmaxRho_logT_min, &
            TmaxRho_logT_max, &
            TmaxRho_logRho_min, &
            TmaxRho_logRho_max, &
            TmaxRho_logT_margin, &
            TmaxRho_logRho_margin, &
            TmaxRho_logRho_dlogRho_min, &
            TmaxRho_logT_dlogT_min, &
            TmaxRho_win_width, &
            TmaxRho_win_aspect_ratio, &
            TmaxRho_file_width, &
            TmaxRho_file_aspect_ratio, &
            TmaxRho_xleft, &
            TmaxRho_xright, &
            TmaxRho_ybot, &
            TmaxRho_ytop, &
            TmaxRho_txt_scale, &
            TmaxRho_title, &
            TmaxRho_use_decorator, &

            History_Track_win_flag, &
            History_Track_file_flag, &
            History_Track_file_interval, &
            History_Track_step_min, &
            History_Track_step_max, &
            show_History_Track_target_box, &
            History_Track_n_sigma, &
            History_Track_xtarget, &
            History_Track_xsigma, &
            History_Track_ytarget, &
            History_Track_ysigma, &
            History_Track_xname, &
            History_Track_xaxis_label, &
            History_Track_yname, &
            History_Track_yaxis_label, &
            History_Track_file_dir, &
            History_Track_file_prefix, &
            show_History_Track_annotation1, &
            show_History_Track_annotation2, &
            show_History_Track_annotation3, &
            History_Track_fname, &
            History_Track_reverse_xaxis, &
            History_Track_reverse_yaxis, &
            History_Track_log_xaxis, &
            History_Track_log_yaxis, &
            History_Track_xmin, &
            History_Track_xmax, &
            History_Track_ymin, &
            History_Track_ymax, &
            History_Track_xmargin, &
            History_Track_ymargin, &
            History_Track_dxmin, &
            History_Track_dymin, &
            History_Track_win_width, &
            History_Track_win_aspect_ratio, &
            History_Track_file_width, &
            History_Track_file_aspect_ratio, &
            History_Track_xleft, &
            History_Track_xright, &
            History_Track_ybot, &
            History_Track_ytop, &
            History_Track_txt_scale, &
            History_Track_title, &
            History_Track_use_decorator, &

            Kipp_win_flag, &
            Kipp_file_flag, &
            Kipp_file_interval, &
            Kipp_xmax, &
            Kipp_xmin, &
            Kipp_max_width, &
            Kipp_step_xmax, &
            Kipp_step_xmin, &
            Kipp_xaxis_name, &
            Kipp_xaxis_log, &
            Kipp_xmargin, &
            Kipp_xaxis_reversed, &
            Kipp_xaxis_in_seconds, &
            Kipp_xaxis_in_Myr, &
            Kipp_xaxis_time_from_present, &
            Kipp_file_dir, &
            Kipp_file_prefix, &
            show_Kipp_annotation1, &
            show_Kipp_annotation2, &
            show_Kipp_annotation3, &
            Kipp_show_burn, &
            Kipp_show_mixing, &
            Kipp_show_luminosities, &
            Kipp_show_mass_boundaries, &
            Kipp_mix_line_weight, &
            Kipp_mix_interval, &
            Kipp_burn_line_weight, &
            Kipp_burn_type_cutoff, &
            Kipp_luminosities_line_weight, &
            Kipp_masses_line_weight, &
            Kipp_mass_max, &
            Kipp_mass_min, &
            Kipp_lgL_max, &
            Kipp_lgL_min, &
            Kipp_mass_margin, &
            Kipp_lgL_margin, &
            Kipp_win_width, &
            Kipp_win_aspect_ratio, &
            Kipp_file_width, &
            Kipp_file_aspect_ratio, &
            Kipp_xleft, &
            Kipp_xright, &
            Kipp_ybot, &
            Kipp_ytop, &
            Kipp_txt_scale, &
            Kipp_title, &
            kipp_use_decorator, &

            rti_win_flag, &
            rti_file_flag, &
            rti_file_interval, &
            rti_xmax, &
            rti_xmin, &
            rti_max_width, &
            rti_step_xmax, &
            rti_step_xmin, &
            rti_xaxis_name, &
            rti_xaxis_log, &
            rti_xmargin, &
            rti_xaxis_reversed, &
            rti_xaxis_in_seconds, &
            rti_xaxis_in_Myr, &
            rti_xaxis_time_from_present, &
            rti_mass_max, &
            rti_mass_min, &
            rti_mass_margin, &
            rti_file_dir, &
            rti_file_prefix, &
            show_rti_annotation1, &
            show_rti_annotation2, &
            show_rti_annotation3, &
            rti_line_weight, &
            rti_interval, &
            rti_win_width, &
            rti_win_aspect_ratio, &
            rti_file_width, &
            rti_file_aspect_ratio, &
            rti_xleft, &
            rti_xright, &
            rti_ybot, &
            rti_ytop, &
            rti_txt_scale, &
            rti_title, &
            rti_use_decorator, &

            TRho_Profile_win_flag, &
            TRho_switch_to_Column_Depth, &
            TRho_switch_to_mass, &
            TRho_Profile_file_flag, &
            TRho_Profile_file_interval, &
            TRho_Profile_file_dir, &
            TRho_Profile_file_prefix, &
            show_TRho_Profile_text_info, &
            show_TRho_Profile_legend, &
            show_TRho_Profile_mass_locs, &
            show_TRho_Profile_burn_labels, &
            show_TRho_Profile_kap_regions, &
            show_TRho_Profile_eos_regions, &
            show_TRho_Profile_gamma1_4_3rd, &
            show_TRho_Profile_burn_lines, &
            show_TRho_Profile_degeneracy_line, &
            show_TRho_Profile_Pgas_Prad_line, &
            show_TRho_Profile_annotation1, &
            show_TRho_Profile_annotation2, &
            show_TRho_Profile_annotation3, &
            TRho_Profile_fname, &
            show_TRho_accretion_mesh_borders, &
            TRho_Profile_text_info_xfac, &
            TRho_Profile_text_info_dxfac, &
            TRho_Profile_text_info_yfac, &
            TRho_Profile_text_info_dyfac, &
            TRho_Profile_xmin, &
            TRho_Profile_xmax, &
            TRho_Profile_ymin, &
            TRho_Profile_ymax, &
            TRho_Profile_legend_coord, &
            TRho_Profile_legend_fjust, &
            TRho_Profile_legend_disp1, &
            TRho_Profile_legend_del_disp, &
            TRho_Profile_legend_txt_scale, &
            TRho_Profile_win_width, &
            TRho_Profile_win_aspect_ratio, &
            TRho_Profile_xleft, &
            TRho_Profile_xright, &
            TRho_Profile_ybot, &
            TRho_Profile_ytop, &
            TRho_Profile_txt_scale, &
            TRho_Profile_title, &
            TRho_Profile_file_width, &
            TRho_Profile_file_aspect_ratio, &
            num_profile_mass_points, &
            profile_mass_point_q, &
            profile_mass_point_color_index, &
            profile_mass_point_symbol, &
            profile_mass_point_symbol_scale, &
            profile_mass_point_str, &
            profile_mass_point_str_clr, &
            profile_mass_point_str_scale, &
            TRho_Profile_use_decorator, &

            Dynamo_win_flag, &
            Dynamo_file_flag, &
            Dynamo_file_interval, &
            Dynamo_file_dir, &
            Dynamo_file_prefix, &
            show_Dynamo_annotation1, &
            show_Dynamo_annotation2, &
            show_Dynamo_annotation3, &
            Dynamo_xaxis_name, &
            Dynamo_xaxis_reversed, &
            Dynamo_xmin, &
            Dynamo_xmax, &
            Dynamo_ymin_left, &
            Dynamo_ymax_left, &
            Dynamo_dymin_left, &
            Dynamo_ymin_right, &
            Dynamo_ymax_right, &
            Dynamo_dymin_right, &
            Dynamo_win_width, &
            Dynamo_win_aspect_ratio, &
            Dynamo_file_width, &
            Dynamo_file_aspect_ratio, &
            Dynamo_xleft, &
            Dynamo_xright, &
            Dynamo_ybot, &
            Dynamo_ytop, &
            Dynamo_txt_scale, &
            Dynamo_title, &
            Dynamo_legend_txt_scale_factor, &
            Dynamo_use_decorator, &

            Mixing_win_flag, &
            Mixing_file_flag, &
            Mixing_file_interval, &
            Mixing_file_dir, &
            Mixing_file_prefix, &
            show_Mixing_annotation1, &
            show_Mixing_annotation2, &
            show_Mixing_annotation3, &
            Mixing_xaxis_name, &
            Mixing_xaxis_reversed, &
            Mixing_legend_txt_scale_factor, &
            Mixing_show_rotation_details, &
            Mixing_xmin, &
            Mixing_xmax, &
            Mixing_ymin, &
            Mixing_ymax, &
            Mixing_dymin, &
            Mixing_win_width, &
            Mixing_win_aspect_ratio, &
            Mixing_file_width, &
            Mixing_file_aspect_ratio, &
            Mixing_xleft, &
            Mixing_xright, &
            Mixing_ybot, &
            Mixing_ytop, &
            Mixing_txt_scale, &
            Mixing_title, &
            Mixing_use_decorator, &

            History_Panels_win_flag, &
            History_Panels_win_width, &
            History_Panels_win_aspect_ratio, &
            History_Panels_xleft, &
            History_Panels_xright, &
            History_Panels_ybot, &
            History_Panels_ytop, &
            History_Panels_txt_scale, &
            History_Panels_title, &
            History_Panels_xmax, &
            History_Panels_xmin, &
            History_Panels_dxmin, &
            History_Panels_max_width, &
            History_Panels_num_panels, &
            History_Panels_xaxis_name, &
            History_Panels_automatic_star_age_units, &
            History_Panels_yaxis_name, &
            History_Panels_xaxis_reversed, &
            History_Panels_yaxis_reversed, &
            History_Panels_yaxis_log, &
            History_Panels_ymin, &
            History_Panels_ymax, &
            History_Panels_dymin, &
            History_Panels_other_yaxis_name, &
            History_Panels_other_yaxis_reversed, &
            History_Panels_xaxis_log, &
            History_Panels_other_yaxis_log, &
            History_Panels_same_yaxis_range, &
            History_Panels_other_ymin, &
            History_Panels_other_ymax, &
            History_Panels_other_dymin, &
            History_Panels_points_name, &
            History_Panels_file_flag, &
            History_Panels_file_dir, &
            History_Panels_file_prefix, &
            History_Panels_file_interval, &
            History_Panels_file_width, &
            History_Panels_file_aspect_ratio, &
            History_Panels_xmargin, &
            History_Panels_ymargin, &
            History_Panels_other_ymargin, &
            History_Panels_use_decorator, &
            
            History_Panel_points_error_bars, &
            History_Panel_points_interval, &
            History_Panel_points_marker, &
            History_Panel_points_ci, &
            History_Panel_points_lw, &
            History_Panel_points_ch, &

            Color_Magnitude_win_flag, &
            Color_Magnitude_win_width, &
            Color_Magnitude_win_aspect_ratio, &
            Color_Magnitude_xleft, &
            Color_Magnitude_xright, &
            Color_Magnitude_ybot, &
            Color_Magnitude_ytop, &
            Color_Magnitude_txt_scale, &
            Color_Magnitude_title, &
            Color_Magnitude_xmax, &
            Color_Magnitude_xmin, &
            Color_Magnitude_dxmin, &
            Color_Magnitude_max_width, &
            Color_Magnitude_num_panels, &
            Color_Magnitude_xaxis1_name, &
            Color_Magnitude_xaxis2_name, &
            Color_Magnitude_yaxis1_name, &
            Color_Magnitude_yaxis2_name, &
            Color_Magnitude_xaxis_reversed, &
            Color_Magnitude_yaxis_reversed, &
            Color_Magnitude_yaxis_log, &
            Color_Magnitude_ymin, &
            Color_Magnitude_ymax, &
            Color_Magnitude_dymin, &
            Color_Magnitude_other_yaxis1_name, &
            Color_Magnitude_other_yaxis2_name, &
            Color_Magnitude_other_yaxis_reversed, &
            Color_Magnitude_xaxis_log, &
            Color_Magnitude_other_yaxis_log, &
            Color_Magnitude_other_ymin, &
            Color_Magnitude_other_ymax, &
            Color_Magnitude_other_dymin, &
            Color_Magnitude_file_flag, &
            Color_Magnitude_file_dir, &
            Color_Magnitude_file_prefix, &
            Color_Magnitude_file_interval, &
            Color_Magnitude_file_width, &
            Color_Magnitude_file_aspect_ratio, &
            Color_Magnitude_xmargin, &
            Color_Magnitude_ymargin, &
            Color_Magnitude_other_ymargin, &
            Color_Magnitude_use_decorator, &

            Mode_Prop_win_flag, &
            Mode_Prop_file_flag, &
            Mode_Prop_file_interval, &
            Mode_Prop_file_dir, &
            Mode_Prop_file_prefix, &
            Mode_Prop_xaxis_name, &
            Mode_Prop_xaxis_reversed, &
            Mode_Prop_xmin, &
            Mode_Prop_xmax, &
            Mode_Prop_ymin, &
            Mode_Prop_ymax, &
            Mode_Prop_win_width, &
            Mode_Prop_win_aspect_ratio, &
            Mode_Prop_file_width, &
            Mode_Prop_file_aspect_ratio, &
            Mode_Prop_nu_max_obs, &
            Mode_Prop_xleft, &
            Mode_Prop_xright, &
            Mode_Prop_ybot, &
            Mode_Prop_ytop, &
            Mode_Prop_txt_scale, &
            Mode_Prop_title, &
            Mode_Prop_use_decorator, &

            Network_win_flag, &
            Network_file_flag, &
            Network_file_interval, &
            Network_file_dir, &
            Network_file_prefix, &
            Network_nmin, &
            Network_nmax, &
            Network_zmin, &
            Network_zmax, &
            Network_show_mass_fraction, &
            Network_log_mass_frac_min, &
            Network_log_mass_frac_max, &
            Network_show_element_names, &
            Network_win_width, &
            Network_win_aspect_ratio, &
            Network_file_width, &
            Network_file_aspect_ratio, &
            Network_xleft, &
            Network_xright, &
            Network_ybot, &
            Network_ytop, &
            Network_txt_scale, &
            Network_title, &
            Network_use_decorator, &
            Network_show_colorbar, &

            Production_win_flag, &
            Production_file_flag, &
            Production_file_interval, &
            Production_file_dir, &
            Production_file_prefix, &
            Production_amin, &
            Production_amax, &
            Production_ymin, &
            Production_ymax, &
            Production_min_mass, &
            Production_max_mass, &
            Production_min_mass_frac, &
            Production_show_element_names, &
            Production_win_width, &
            Production_win_aspect_ratio, &
            Production_file_width, &
            Production_file_aspect_ratio, &
            Production_xleft, &
            Production_xright, &
            Production_ybot, &
            Production_ytop, &
            Production_txt_scale, &
            Production_title, &
            Production_use_decorator, &

            Summary_Burn_win_flag, &
            Summary_Burn_file_flag, &
            Summary_Burn_file_interval, &
            Summary_Burn_file_dir, &
            Summary_Burn_file_prefix, &
            Summary_Burn_xaxis_name, &
            Summary_Burn_xaxis_reversed, &
            Summary_Burn_xmin, &
            Summary_Burn_xmax, &
            Summary_Burn_win_width, &
            Summary_Burn_win_aspect_ratio, &
            Summary_Burn_file_width, &
            Summary_Burn_file_aspect_ratio, &
            Summary_Burn_xleft, &
            Summary_Burn_xright, &
            Summary_Burn_ybot, &
            Summary_Burn_ytop, &
            Summary_Burn_txt_scale, &
            Summary_Burn_title, &
            Summary_Burn_title_shift, &
            Summary_Burn_use_decorator, &
            
            Summary_Profile_win_flag, &
            Summary_Profile_file_flag, &
            Summary_Profile_file_interval, &
            Summary_Profile_file_dir, &
            Summary_Profile_file_prefix, &
            Summary_Profile_xaxis_name, &
            Summary_Profile_xaxis_reversed, &
            Summary_Profile_scaled_value, &
            Summary_Profile_xmin, &
            Summary_Profile_xmax, &
            Summary_Profile_win_width, &
            Summary_Profile_win_aspect_ratio, &
            Summary_Profile_file_width, &
            Summary_Profile_file_aspect_ratio, &
            Summary_Profile_xleft, &
            Summary_Profile_xright, &
            Summary_Profile_ybot, &
            Summary_Profile_ytop, &
            Summary_Profile_txt_scale, &
            Summary_Profile_title, &
            Summary_Profile_name, &
            Summary_Profile_legend, &
            Summary_Profile_num_lines, &
            Summary_Profile_use_decorator, &
            
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

            Power_win_flag, &
            Power_file_flag, &
            Power_file_interval, &
            Power_file_dir, &
            Power_file_prefix, &
            Power_xaxis_name, &
            Power_xaxis_reversed, &
            Power_legend_max_cnt, &
            Power_legend_txt_scale_factor, &
            Power_xmin, &
            Power_xmax, &
            Power_ymin, &
            Power_ymax, &
            Power_win_width, &
            Power_win_aspect_ratio, &
            Power_file_width, &
            Power_file_aspect_ratio, &
            Power_xleft, &
            Power_xright, &
            Power_ybot, &
            Power_ytop, &
            Power_txt_scale, &
            Power_title, &
            Power_use_decorator, &

            Abundance_win_flag, &
            Abundance_file_flag, &
            Abundance_file_interval, &
            Abundance_file_dir, &
            Abundance_file_prefix, &
            Abundance_num_isos_to_show, &
            num_abundance_line_labels, &
            Abundance_legend_txt_scale_factor, &
            Abundance_line_txt_scale_factor, &
            Abundance_legend_max_cnt, &
            Abundance_which_isos_to_show, &
            Abundance_xaxis_name, &
            Abundance_xaxis_reversed, &
            Abundance_show_photosphere_location, &
            Abundance_xmin, &
            Abundance_xmax, &
            Abundance_log_mass_frac_min, &
            Abundance_log_mass_frac_max, &
            Abundance_use_decorator, &
            Abundance_win_width, &
            Abundance_win_aspect_ratio, &
            Abundance_file_width, &
            Abundance_file_aspect_ratio, &
            Abundance_xleft, &
            Abundance_xright, &
            Abundance_ybot, &
            Abundance_ytop, &
            Abundance_txt_scale, &
            Abundance_title, &
            Abundance_use_decorator, &


            dPg_dnu_win_flag, &
            dPg_dnu_file_flag, &
            dPg_dnu_file_interval, &
            dPg_dnu_step_min, &
            dPg_dnu_step_max, &
            dPg_dnu_file_dir, &
            dPg_dnu_file_prefix, &
            show_dPg_dnu_annotation1, &
            show_dPg_dnu_annotation2, &
            show_dPg_dnu_annotation3, &
            dPg_dnu_fname, &
            dPg_dnu_delta_nu_min, &
            dPg_dnu_delta_nu_max, &
            dPg_dnu_delta_Pg_min, &
            dPg_dnu_delta_Pg_max, &
            dPg_dnu_delta_nu_margin, &
            dPg_dnu_delta_Pg_margin, &
            dPg_dnu_d_delta_nu_min, &
            dPg_dnu_d_delta_Pg_min, &
            dPg_dnu_win_width, &
            dPg_dnu_win_aspect_ratio, &
            dPg_dnu_xleft, &
            dPg_dnu_xright, &
            dPg_dnu_ybot, &
            dPg_dnu_ytop, &
            dPg_dnu_txt_scale, &
            dPg_dnu_file_width, &
            dPg_dnu_file_aspect_ratio, &
            dPg_dnu_title, &
            dPg_dnu_use_decorator, &

            Grid_win_flag, &
            Grid_win_width, &
            Grid_win_aspect_ratio, &
            Grid_xleft, &
            Grid_xright, &
            Grid_ybot, &
            Grid_ytop, &
            Grid_title, &
            Grid_txt_scale_factor, &
            Grid_num_cols, &
            Grid_num_rows, &
            Grid_num_plots, &
            Grid_plot_name, &
            Grid_plot_row, &
            Grid_plot_rowspan, &
            Grid_plot_col, &
            Grid_plot_colspan, &
            Grid_plot_pad_left, &
            Grid_plot_pad_right, &
            Grid_plot_pad_top, &
            Grid_plot_pad_bot, &
            Grid_file_flag, &
            Grid_file_dir, &
            Grid_file_prefix, &
            Grid_file_interval, &
            Grid_file_width, &
            Grid_file_aspect_ratio, &

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

            read_extra_pgstar_inlist, &
            extra_pgstar_inlist_name

      contains


      subroutine read_pgstar(s, filename, ierr)
         use star_private_def
         use utils_lib
         type (star_info), pointer :: s
         character(*), intent(in) :: filename
         integer, intent(out) :: ierr
         character (len=strlen) :: pgstar_namelist_name
         pgstar_namelist_name = ''
         ierr = 0
         call set_default_pgstar_controls
         call read_pgstar_file(s, filename, 1, ierr)
      end subroutine read_pgstar


      recursive subroutine read_pgstar_file(s, filename, level, ierr)
         use star_private_def
         use utils_lib
         character(*), intent(in) :: filename
         type (star_info), pointer :: s
         integer, intent(in) :: level
         integer, intent(out) :: ierr
         logical, dimension(max_extra_inlists) :: read_extra
         character (len=strlen) :: message
         character (len=strlen), dimension(max_extra_inlists) :: extra
         integer :: unit, i

         ierr = 0

         if (level >= 10) then
            write(*,*) 'ERROR: too many levels of nested extra pgstar inlist files'
            ierr = -1
            return
         end if

         if (len_trim(filename) > 0) then
            open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
            if (ierr /= 0) then
               write(*, *) 'Failed to open pgstar namelist file ', trim(filename)
               return
            end if
            read(unit, nml=pgstar, iostat=ierr)
            close(unit)
            if (ierr /= 0) then
               write(*, *)
               write(*, *)
               write(*, *)
               write(*, *)
               write(*, '(a)') &
                  'Failed while trying to read pgstar namelist file: ' // trim(filename)
               write(*, '(a)') &
                  'Perhaps the following runtime error message will help you find the problem.'
               write(*, *)
               open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
               read(unit, nml=pgstar)
               close(unit)
               return
            end if
         end if

         call store_pgstar_controls(s, ierr)

         ! recursive calls to read other inlists
         do i=1, max_extra_inlists
            read_extra(i) = read_extra_pgstar_inlist(i)
            read_extra_pgstar_inlist(i) = .false.
            extra(i) = extra_pgstar_inlist_name(i)
            extra_pgstar_inlist_name(i) = 'undefined'
   
            if (read_extra(i)) then
               call read_pgstar_file(s, extra(i), level+1, ierr)
               if (ierr /= 0) return
            end if
         end do

      end subroutine read_pgstar_file


      subroutine store_pgstar_controls(s, ierr)
         use star_private_def
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         ierr = 0

         s% pg% file_device = file_device
         s% pg% file_extension = file_extension
         s% pg% file_digits = file_digits
         s% pg% pgstar_interval = pgstar_interval
         s% pg% pause = pause
         s% pg% pause_interval = pause_interval
         s% pg% pgstar_sleep = pgstar_sleep
         s% pg% clear_history = clear_history

         s% pg% file_white_on_black_flag = file_white_on_black_flag
         s% pg% delta_HR_limit_for_file_output = delta_HR_limit_for_file_output
         s% pg% win_white_on_black_flag = win_white_on_black_flag
         s% pg% pgstar_show_model_number = pgstar_show_model_number
         s% pg% pgstar_show_age = pgstar_show_age
         s% pg% pgstar_show_age_in_seconds = pgstar_show_age_in_seconds
         s% pg% pgstar_show_age_in_minutes = pgstar_show_age_in_minutes
         s% pg% pgstar_show_age_in_hours = pgstar_show_age_in_hours
         s% pg% pgstar_show_age_in_days = pgstar_show_age_in_days
         s% pg% pgstar_show_age_in_years = pgstar_show_age_in_years
         s% pg% pgstar_show_log_age_in_years = pgstar_show_log_age_in_years

         s% pg% pgstar_report_writing_files = pgstar_report_writing_files

         s% pg% pgstar_show_title = pgstar_show_title
         s% pg% pgstar_title_scale = pgstar_title_scale
         s% pg% pgstar_title_disp = pgstar_title_disp
         s% pg% pgstar_title_coord = pgstar_title_coord
         s% pg% pgstar_title_fjust = pgstar_title_fjust
         s% pg% pgstar_title_lw = pgstar_title_lw

         s% pg% pgstar_grid_show_title = pgstar_grid_show_title
         s% pg% pgstar_grid_title_scale = pgstar_grid_title_scale
         s% pg% pgstar_grid_title_disp = pgstar_grid_title_disp
         s% pg% pgstar_grid_title_coord = pgstar_grid_title_coord
         s% pg% pgstar_grid_title_fjust = pgstar_grid_title_fjust
         s% pg% pgstar_grid_title_lw = pgstar_grid_title_lw

         s% pg% pgstar_age_scale = pgstar_age_scale
         s% pg% pgstar_age_disp = pgstar_age_disp
         s% pg% pgstar_age_coord = pgstar_age_coord
         s% pg% pgstar_age_fjust = pgstar_age_fjust
         s% pg% pgstar_age_lw = pgstar_age_lw

         s% pg% pgstar_model_scale = pgstar_model_scale
         s% pg% pgstar_model_disp = pgstar_model_disp
         s% pg% pgstar_model_coord = pgstar_model_coord
         s% pg% pgstar_model_fjust = pgstar_model_fjust
         s% pg% pgstar_model_lw = pgstar_model_lw

         s% pg% pgstar_xaxis_label_scale = pgstar_xaxis_label_scale
         s% pg% pgstar_left_yaxis_label_scale = pgstar_left_yaxis_label_scale
         s% pg% pgstar_right_yaxis_label_scale = pgstar_right_yaxis_label_scale
         s% pg% pgstar_xaxis_label_lw = pgstar_xaxis_label_lw
         s% pg% pgstar_left_yaxis_label_lw = pgstar_left_yaxis_label_lw
         s% pg% pgstar_right_yaxis_label_lw = pgstar_right_yaxis_label_lw
         s% pg% pgstar_xaxis_label_disp = pgstar_xaxis_label_disp
         s% pg% pgstar_left_yaxis_label_disp = pgstar_left_yaxis_label_disp
         s% pg% pgstar_right_yaxis_label_disp = pgstar_right_yaxis_label_disp
         s% pg% pgstar_num_scale = pgstar_num_scale
         s% pg% pgstar_lw = pgstar_lw
         s% pg% pgstar_profile_line_style = pgstar_profile_line_style
         s% pg% pgstar_history_line_style = pgstar_history_line_style
         s% pg% pgstar_model_lw = pgstar_model_lw
         s% pg% pgstar_box_lw = pgstar_box_lw

         s% pg% Profile_Panels_show_Mach_1_location = &
            Profile_Panels_show_Mach_1_location
         s% pg% Profile_Panels_show_photosphere_location = &
            Profile_Panels_show_photosphere_location
         s% pg% Profile_Panels_xwidth_left_div_shock_value = &
            Profile_Panels_xwidth_left_div_shock_value
         s% pg% Profile_Panels_xwidth_right_div_shock_value = &
            Profile_Panels_xwidth_right_div_shock_value
         s% pg% Profile_Panels_xwidth_left_of_shock = &
            Profile_Panels_xwidth_left_of_shock
         s% pg% Profile_Panels_xwidth_right_of_shock = &
            Profile_Panels_xwidth_right_of_shock

         s% pg% Profile_Panels_win_flag = Profile_Panels_win_flag
         s% pg% Profile_Panels_file_flag = Profile_Panels_file_flag
         s% pg% Profile_Panels_file_interval = Profile_Panels_file_interval
         s% pg% Profile_Panels_file_dir = Profile_Panels_file_dir
         s% pg% Profile_Panels_file_prefix = Profile_Panels_file_prefix
         s% pg% Profile_Panels_xaxis_reversed = Profile_Panels_xaxis_reversed
         s% pg% Profile_Panels_xaxis_name = Profile_Panels_xaxis_name
         s% pg% Profile_Panels_title = Profile_Panels_title
         s% pg% Profile_Panels_xmin = Profile_Panels_xmin
         s% pg% Profile_Panels_xmax = Profile_Panels_xmax
         s% pg% Profile_Panels_xmargin = Profile_Panels_xmargin
         s% pg% Profile_Panels_show_mix_regions_on_xaxis = &
            Profile_Panels_show_mix_regions_on_xaxis
         s% pg% Profile_Panels_win_width = Profile_Panels_win_width
         s% pg% Profile_Panels_win_aspect_ratio = Profile_Panels_win_aspect_ratio
         s% pg% Profile_Panels_xleft = Profile_Panels_xleft
         s% pg% Profile_Panels_xright = Profile_Panels_xright
         s% pg% Profile_Panels_ybot = Profile_Panels_ybot
         s% pg% Profile_Panels_ytop = Profile_Panels_ytop
         s% pg% Profile_Panels_txt_scale = Profile_Panels_txt_scale
         s% pg% prev_Profile_Panels_win_width = prev_Profile_Panels_win_width
         s% pg% prev_Profile_Panels_win_ratio = prev_Profile_Panels_win_ratio
         s% pg% Profile_Panels_file_width = Profile_Panels_file_width
         s% pg% Profile_Panels_file_aspect_ratio = Profile_Panels_file_aspect_ratio
         s% pg% prev_Profile_Panels_file_width = prev_Profile_Panels_file_width
         s% pg% prev_Profile_Panels_file_ratio = prev_Profile_Panels_file_ratio
         s% pg% Profile_Panels_num_panels = Profile_Panels_num_panels
         s% pg% Profile_Panels_yaxis_name = Profile_Panels_yaxis_name
         s% pg% Profile_Panels_other_yaxis_name = Profile_Panels_other_yaxis_name
         s% pg% Profile_Panels_yaxis_reversed = Profile_Panels_yaxis_reversed
         s% pg% Profile_Panels_other_yaxis_reversed = Profile_Panels_other_yaxis_reversed
         s% pg% Profile_Panels_yaxis_log = Profile_Panels_yaxis_log
         s% pg% Profile_Panels_other_yaxis_log = Profile_Panels_other_yaxis_log
         s% pg% Profile_Panels_same_yaxis_range = Profile_Panels_same_yaxis_range
         s% pg% Profile_Panels_ymin = Profile_Panels_ymin
         s% pg% Profile_Panels_other_ymin = Profile_Panels_other_ymin
         s% pg% Profile_Panels_ymax = Profile_Panels_ymax
         s% pg% Profile_Panels_other_ymax = Profile_Panels_other_ymax
         s% pg% Profile_Panels_ycenter = Profile_Panels_ycenter
         s% pg% Profile_Panels_other_ycenter = Profile_Panels_other_ycenter
         s% pg% Profile_Panels_ymargin = Profile_Panels_ymargin
         s% pg% Profile_Panels_other_ymargin = Profile_Panels_other_ymargin
         s% pg% Profile_Panels_dymin = Profile_Panels_dymin
         s% pg% Profile_Panels_other_dymin = Profile_Panels_other_dymin
         s% pg% Profile_Panels_show_grid = Profile_Panels_show_grid
         s% pg% Profile_Panels_use_decorator = Profile_Panels_use_decorator

         s% pg% Text_Summary_win_flag = Text_Summary_win_flag
         s% pg% Text_Summary_file_flag = Text_Summary_file_flag
         s% pg% Text_Summary_file_interval = Text_Summary_file_interval
         s% pg% Text_Summary_file_dir = Text_Summary_file_dir
         s% pg% Text_Summary_file_prefix = Text_Summary_file_prefix
         s% pg% Text_Summary_num_cols = Text_Summary_num_cols
         s% pg% Text_Summary_num_rows = Text_Summary_num_rows
         s% pg% Text_Summary_name = Text_Summary_name
         s% pg% Text_Summary_win_width = Text_Summary_win_width
         s% pg% Text_Summary_win_aspect_ratio = Text_Summary_win_aspect_ratio
         s% pg% Text_Summary_file_width = Text_Summary_file_width
         s% pg% Text_Summary_file_aspect_ratio = Text_Summary_file_aspect_ratio
         s% pg% Text_Summary_title = Text_Summary_title
         s% pg% Text_Summary_xleft = Text_Summary_xleft
         s% pg% Text_Summary_xright = Text_Summary_xright
         s% pg% Text_Summary_ybot = Text_Summary_ybot
         s% pg% Text_Summary_ytop = Text_Summary_ytop
         s% pg% Text_Summary_txt_scale = Text_Summary_txt_scale
         s% pg% Text_Summary_dxval = Text_Summary_dxval

         s% pg% logg_Teff_win_flag = logg_Teff_win_flag
         s% pg% logg_Teff_file_flag = logg_Teff_file_flag
         s% pg% show_logg_Teff_target_box = show_logg_Teff_target_box
         s% pg% logg_Teff_target_n_sigma = logg_Teff_target_n_sigma
         s% pg% logg_Teff_target_logg = logg_Teff_target_logg
         s% pg% logg_Teff_target_logg_sigma = logg_Teff_target_logg_sigma
         s% pg% logg_Teff_target_Teff = logg_Teff_target_Teff
         s% pg% logg_Teff_target_Teff_sigma = logg_Teff_target_Teff_sigma
         s% pg% logg_Teff_file_interval = logg_Teff_file_interval
         s% pg% logg_Teff_step_min = logg_Teff_step_min
         s% pg% logg_Teff_step_max = logg_Teff_step_max
         s% pg% logg_Teff_file_dir = logg_Teff_file_dir
         s% pg% logg_Teff_file_prefix = logg_Teff_file_prefix
         s% pg% show_logg_Teff_annotation1 = show_logg_Teff_annotation1
         s% pg% show_logg_Teff_annotation2 = show_logg_Teff_annotation2
         s% pg% show_logg_Teff_annotation3 = show_logg_Teff_annotation3
         s% pg% logg_Teff_fname = logg_Teff_fname
         s% pg% logg_Teff_title = logg_Teff_title
         s% pg% logg_Teff_logg_min = logg_Teff_logg_min
         s% pg% logg_Teff_logg_max = logg_Teff_logg_max
         s% pg% logg_Teff_Teff_min = logg_Teff_Teff_min
         s% pg% logg_Teff_Teff_max = logg_Teff_Teff_max
         s% pg% logg_Teff_Teff_margin = logg_Teff_Teff_margin
         s% pg% logg_Teff_logg_margin = logg_Teff_logg_margin
         s% pg% logg_Teff_dTeff_min = logg_Teff_dTeff_min
         s% pg% logg_Teff_dlogg_min = logg_Teff_dlogg_min
         s% pg% logg_Teff_win_width = logg_Teff_win_width
         s% pg% logg_Teff_win_aspect_ratio = logg_Teff_win_aspect_ratio
         s% pg% logg_Teff_xleft = logg_Teff_xleft
         s% pg% logg_Teff_xright = logg_Teff_xright
         s% pg% logg_Teff_ybot = logg_Teff_ybot
         s% pg% logg_Teff_ytop = logg_Teff_ytop
         s% pg% logg_Teff_txt_scale = logg_Teff_txt_scale
         s% pg% logg_Teff_file_width = logg_Teff_file_width
         s% pg% logg_Teff_file_aspect_ratio = logg_Teff_file_aspect_ratio
         s% pg% logg_Teff_use_decorator = logg_Teff_use_decorator

         s% pg% logL_Teff_win_flag = logL_Teff_win_flag
         s% pg% logL_Teff_file_flag = logL_Teff_file_flag
         s% pg% show_logL_Teff_target_box = show_logL_Teff_target_box
         s% pg% logL_Teff_target_n_sigma = logL_Teff_target_n_sigma
         s% pg% logL_Teff_target_logL = logL_Teff_target_logL
         s% pg% logL_Teff_target_logL_sigma = logL_Teff_target_logL_sigma
         s% pg% logL_Teff_target_Teff = logL_Teff_target_Teff
         s% pg% logL_Teff_target_Teff_sigma = logL_Teff_target_Teff_sigma
         s% pg% logL_Teff_file_interval = logL_Teff_file_interval
         s% pg% logL_Teff_step_min = logL_Teff_step_min
         s% pg% logL_Teff_step_max = logL_Teff_step_max
         s% pg% logL_Teff_file_dir = logL_Teff_file_dir
         s% pg% logL_Teff_file_prefix = logL_Teff_file_prefix
         s% pg% show_logL_Teff_annotation1 = show_logL_Teff_annotation1
         s% pg% show_logL_Teff_annotation2 = show_logL_Teff_annotation2
         s% pg% show_logL_Teff_annotation3 = show_logL_Teff_annotation3
         s% pg% logL_Teff_fname = logL_Teff_fname
         s% pg% logL_Teff_title = logL_Teff_title
         s% pg% logL_Teff_logL_min = logL_Teff_logL_min
         s% pg% logL_Teff_logL_max = logL_Teff_logL_max
         s% pg% logL_Teff_Teff_min = logL_Teff_Teff_min
         s% pg% logL_Teff_Teff_max = logL_Teff_Teff_max
         s% pg% logL_Teff_Teff_margin = logL_Teff_Teff_margin
         s% pg% logL_Teff_logL_margin = logL_Teff_logL_margin
         s% pg% logL_Teff_dTeff_min = logL_Teff_dTeff_min
         s% pg% logL_Teff_dlogL_min = logL_Teff_dlogL_min
         s% pg% logL_Teff_win_width = logL_Teff_win_width
         s% pg% logL_Teff_win_aspect_ratio = logL_Teff_win_aspect_ratio
         s% pg% logL_Teff_xleft = logL_Teff_xleft
         s% pg% logL_Teff_xright = logL_Teff_xright
         s% pg% logL_Teff_ybot = logL_Teff_ybot
         s% pg% logL_Teff_ytop = logL_Teff_ytop
         s% pg% logL_Teff_txt_scale = logL_Teff_txt_scale
         s% pg% logL_Teff_file_width = logL_Teff_file_width
         s% pg% logL_Teff_file_aspect_ratio = logL_Teff_file_aspect_ratio
         s% pg% logL_Teff_use_decorator = logL_Teff_use_decorator

         s% pg% logL_R_win_flag = logL_R_win_flag
         s% pg% logL_R_file_flag = logL_R_file_flag
         s% pg% show_logL_R_target_box = show_logL_R_target_box
         s% pg% logL_R_target_n_sigma = logL_R_target_n_sigma
         s% pg% logL_R_target_logL = logL_R_target_logL
         s% pg% logL_R_target_logL_sigma = logL_R_target_logL_sigma
         s% pg% logL_R_target_R = logL_R_target_R
         s% pg% logL_R_target_R_sigma = logL_R_target_R_sigma
         s% pg% logL_R_file_interval = logL_R_file_interval
         s% pg% logL_R_step_min = logL_R_step_min
         s% pg% logL_R_step_max = logL_R_step_max
         s% pg% logL_R_file_dir = logL_R_file_dir
         s% pg% logL_R_file_prefix = logL_R_file_prefix
         s% pg% show_logL_R_annotation1 = show_logL_R_annotation1
         s% pg% show_logL_R_annotation2 = show_logL_R_annotation2
         s% pg% show_logL_R_annotation3 = show_logL_R_annotation3
         s% pg% logL_R_fname = logL_R_fname
         s% pg% logL_R_title = logL_R_title
         s% pg% show_logL_photosphere_r = show_logL_photosphere_r
         s% pg% logL_R_logL_min = logL_R_logL_min
         s% pg% logL_R_logL_max = logL_R_logL_max
         s% pg% logL_R_R_min = logL_R_R_min
         s% pg% logL_R_R_max = logL_R_R_max
         s% pg% logL_R_R_margin = logL_R_R_margin
         s% pg% logL_R_logL_margin = logL_R_logL_margin
         s% pg% logL_R_dR_min = logL_R_dR_min
         s% pg% logL_R_dlogL_min = logL_R_dlogL_min
         s% pg% logL_R_win_width = logL_R_win_width
         s% pg% logL_R_win_aspect_ratio = logL_R_win_aspect_ratio
         s% pg% logL_R_xleft = logL_R_xleft
         s% pg% logL_R_xright = logL_R_xright
         s% pg% logL_R_ybot = logL_R_ybot
         s% pg% logL_R_ytop = logL_R_ytop
         s% pg% logL_R_txt_scale = logL_R_txt_scale
         s% pg% logL_R_file_width = logL_R_file_width
         s% pg% logL_R_file_aspect_ratio = logL_R_file_aspect_ratio
         s% pg% logL_R_use_decorator = logL_R_use_decorator

         s% pg% logL_v_win_flag = logL_v_win_flag
         s% pg% logL_v_file_flag = logL_v_file_flag
         s% pg% show_logL_v_target_box = show_logL_v_target_box
         s% pg% logL_v_target_n_sigma = logL_v_target_n_sigma
         s% pg% logL_v_target_logL = logL_v_target_logL
         s% pg% logL_v_target_logL_sigma = logL_v_target_logL_sigma
         s% pg% logL_v_target_v = logL_v_target_v
         s% pg% logL_v_target_v_sigma = logL_v_target_v_sigma
         s% pg% logL_v_file_interval = logL_v_file_interval
         s% pg% logL_v_step_min = logL_v_step_min
         s% pg% logL_v_step_max = logL_v_step_max
         s% pg% logL_v_file_dir = logL_v_file_dir
         s% pg% logL_v_file_prefix = logL_v_file_prefix
         s% pg% show_logL_v_annotation1 = show_logL_v_annotation1
         s% pg% show_logL_v_annotation2 = show_logL_v_annotation2
         s% pg% show_logL_v_annotation3 = show_logL_v_annotation3
         s% pg% logL_v_fname = logL_v_fname
         s% pg% logL_v_title = logL_v_title
         s% pg% show_logL_photosphere_v = show_logL_photosphere_v
         s% pg% logL_v_logL_min = logL_v_logL_min
         s% pg% logL_v_logL_max = logL_v_logL_max
         s% pg% logL_v_v_min = logL_v_v_min
         s% pg% logL_v_v_max = logL_v_v_max
         s% pg% logL_v_v_margin = logL_v_v_margin
         s% pg% logL_v_logL_margin = logL_v_logL_margin
         s% pg% logL_v_dv_min = logL_v_dv_min
         s% pg% logL_v_dlogL_min = logL_v_dlogL_min
         s% pg% logL_v_win_width = logL_v_win_width
         s% pg% logL_v_win_aspect_ratio = logL_v_win_aspect_ratio
         s% pg% logL_v_xleft = logL_v_xleft
         s% pg% logL_v_xright = logL_v_xright
         s% pg% logL_v_ybot = logL_v_ybot
         s% pg% logL_v_ytop = logL_v_ytop
         s% pg% logL_v_txt_scale = logL_v_txt_scale
         s% pg% logL_v_file_width = logL_v_file_width
         s% pg% logL_v_file_aspect_ratio = logL_v_file_aspect_ratio
         s% pg% logL_v_use_decorator = logL_v_use_decorator

         s% pg% L_Teff_win_flag = L_Teff_win_flag
         s% pg% L_Teff_file_flag = L_Teff_file_flag
         s% pg% show_L_Teff_target_box = show_L_Teff_target_box
         s% pg% L_Teff_target_n_sigma = L_Teff_target_n_sigma
         s% pg% L_Teff_target_L = L_Teff_target_L
         s% pg% L_Teff_target_L_sigma = L_Teff_target_L_sigma
         s% pg% L_Teff_target_Teff = L_Teff_target_Teff
         s% pg% L_Teff_target_Teff_sigma = L_Teff_target_Teff_sigma
         s% pg% L_Teff_file_interval = L_Teff_file_interval
         s% pg% L_Teff_step_min = L_Teff_step_min
         s% pg% L_Teff_step_max = L_Teff_step_max
         s% pg% L_Teff_file_dir = L_Teff_file_dir
         s% pg% L_Teff_file_prefix = L_Teff_file_prefix
         s% pg% show_L_Teff_annotation1 = show_L_Teff_annotation1
         s% pg% show_L_Teff_annotation2 = show_L_Teff_annotation2
         s% pg% show_L_Teff_annotation3 = show_L_Teff_annotation3
         s% pg% L_Teff_fname = L_Teff_fname
         s% pg% L_Teff_title = L_Teff_title
         s% pg% L_Teff_L_min = L_Teff_L_min
         s% pg% L_Teff_L_max = L_Teff_L_max
         s% pg% L_Teff_Teff_min = L_Teff_Teff_min
         s% pg% L_Teff_Teff_max = L_Teff_Teff_max
         s% pg% L_Teff_Teff_margin = L_Teff_Teff_margin
         s% pg% L_Teff_L_margin = L_Teff_L_margin
         s% pg% L_Teff_dTeff_min = L_Teff_dTeff_min
         s% pg% L_Teff_dL_min = L_Teff_dL_min
         s% pg% L_Teff_win_width = L_Teff_win_width
         s% pg% L_Teff_win_aspect_ratio = L_Teff_win_aspect_ratio
         s% pg% L_Teff_xleft = L_Teff_xleft
         s% pg% L_Teff_xright = L_Teff_xright
         s% pg% L_Teff_ybot = L_Teff_ybot
         s% pg% L_Teff_ytop = L_Teff_ytop
         s% pg% L_Teff_txt_scale = L_Teff_txt_scale
         s% pg% L_Teff_file_width = L_Teff_file_width
         s% pg% L_Teff_file_aspect_ratio = L_Teff_file_aspect_ratio
         s% pg% L_Teff_use_decorator = L_Teff_use_decorator

         s% pg% L_v_win_flag = L_v_win_flag
         s% pg% L_v_file_flag = L_v_file_flag
         s% pg% show_L_v_target_box = show_L_v_target_box
         s% pg% L_v_target_n_sigma = L_v_target_n_sigma
         s% pg% L_v_target_L = L_v_target_L
         s% pg% L_v_target_L_sigma = L_v_target_L_sigma
         s% pg% L_v_target_v = L_v_target_v
         s% pg% L_v_target_v_sigma = L_v_target_v_sigma
         s% pg% L_v_file_interval = L_v_file_interval
         s% pg% L_v_step_min = L_v_step_min
         s% pg% L_v_step_max = L_v_step_max
         s% pg% L_v_file_dir = L_v_file_dir
         s% pg% L_v_file_prefix = L_v_file_prefix
         s% pg% show_L_v_annotation1 = show_L_v_annotation1
         s% pg% show_L_v_annotation2 = show_L_v_annotation2
         s% pg% show_L_v_annotation3 = show_L_v_annotation3
         s% pg% L_v_fname = L_v_fname
         s% pg% L_v_title = L_v_title
         s% pg% L_v_L_min = L_v_L_min
         s% pg% L_v_L_max = L_v_L_max
         s% pg% L_v_v_min = L_v_v_min
         s% pg% L_v_v_max = L_v_v_max
         s% pg% L_v_v_margin = L_v_v_margin
         s% pg% L_v_L_margin = L_v_L_margin
         s% pg% L_v_dv_min = L_v_dv_min
         s% pg% L_v_dL_min = L_v_dL_min
         s% pg% L_v_win_width = L_v_win_width
         s% pg% L_v_win_aspect_ratio = L_v_win_aspect_ratio
         s% pg% L_v_xleft = L_v_xleft
         s% pg% L_v_xright = L_v_xright
         s% pg% L_v_ybot = L_v_ybot
         s% pg% L_v_ytop = L_v_ytop
         s% pg% L_v_txt_scale = L_v_txt_scale
         s% pg% L_v_file_width = L_v_file_width
         s% pg% L_v_file_aspect_ratio = L_v_file_aspect_ratio
         s% pg% L_v_use_decorator = L_v_use_decorator

         s% pg% L_R_win_flag = L_R_win_flag
         s% pg% L_R_file_flag = L_R_file_flag
         s% pg% show_L_R_target_box = show_L_R_target_box
         s% pg% L_R_target_n_sigma = L_R_target_n_sigma
         s% pg% L_R_target_L = L_R_target_L
         s% pg% L_R_target_L_sigma = L_R_target_L_sigma
         s% pg% L_R_target_R = L_R_target_R
         s% pg% L_R_target_R_sigma = L_R_target_R_sigma
         s% pg% L_R_file_interval = L_R_file_interval
         s% pg% L_R_step_min = L_R_step_min
         s% pg% L_R_step_max = L_R_step_max
         s% pg% L_R_file_dir = L_R_file_dir
         s% pg% L_R_file_prefix = L_R_file_prefix
         s% pg% show_L_R_annotation1 = show_L_R_annotation1
         s% pg% show_L_R_annotation2 = show_L_R_annotation2
         s% pg% show_L_R_annotation3 = show_L_R_annotation3
         s% pg% L_R_fname = L_R_fname
         s% pg% L_R_title = L_R_title
         s% pg% L_R_L_min = L_R_L_min
         s% pg% L_R_L_max = L_R_L_max
         s% pg% L_R_R_min = L_R_R_min
         s% pg% L_R_R_max = L_R_R_max
         s% pg% L_R_R_margin = L_R_R_margin
         s% pg% L_R_L_margin = L_R_L_margin
         s% pg% L_R_dR_min = L_R_dR_min
         s% pg% L_R_dL_min = L_R_dL_min
         s% pg% L_R_win_width = L_R_win_width
         s% pg% L_R_win_aspect_ratio = L_R_win_aspect_ratio
         s% pg% L_R_xleft = L_R_xleft
         s% pg% L_R_xright = L_R_xright
         s% pg% L_R_ybot = L_R_ybot
         s% pg% L_R_ytop = L_R_ytop
         s% pg% L_R_txt_scale = L_R_txt_scale
         s% pg% L_R_file_width = L_R_file_width
         s% pg% L_R_file_aspect_ratio = L_R_file_aspect_ratio
         s% pg% L_R_use_decorator = L_R_use_decorator

         s% pg% R_Teff_win_flag = R_Teff_win_flag
         s% pg% R_Teff_file_flag = R_Teff_file_flag
         s% pg% show_R_Teff_target_box = show_R_Teff_target_box
         s% pg% R_Teff_target_n_sigma = R_Teff_target_n_sigma
         s% pg% R_Teff_target_R = R_Teff_target_R
         s% pg% R_Teff_target_R_sigma = R_Teff_target_R_sigma
         s% pg% R_Teff_target_Teff = R_Teff_target_Teff
         s% pg% R_Teff_target_Teff_sigma = R_Teff_target_Teff_sigma
         s% pg% R_Teff_file_interval = R_Teff_file_interval
         s% pg% R_Teff_step_min = R_Teff_step_min
         s% pg% R_Teff_step_max = R_Teff_step_max
         s% pg% R_Teff_file_dir = R_Teff_file_dir
         s% pg% R_Teff_file_prefix = R_Teff_file_prefix
         s% pg% show_R_Teff_annotation1 = show_R_Teff_annotation1
         s% pg% show_R_Teff_annotation2 = show_R_Teff_annotation2
         s% pg% show_R_Teff_annotation3 = show_R_Teff_annotation3
         s% pg% R_Teff_fname = R_Teff_fname
         s% pg% R_Teff_title = R_Teff_title
         s% pg% R_Teff_R_min = R_Teff_R_min
         s% pg% R_Teff_R_max = R_Teff_R_max
         s% pg% R_Teff_Teff_min = R_Teff_Teff_min
         s% pg% R_Teff_Teff_max = R_Teff_Teff_max
         s% pg% R_Teff_Teff_margin = R_Teff_Teff_margin
         s% pg% R_Teff_R_margin = R_Teff_R_margin
         s% pg% R_Teff_dTeff_min = R_Teff_dTeff_min
         s% pg% R_Teff_dR_min = R_Teff_dR_min
         s% pg% R_Teff_win_width = R_Teff_win_width
         s% pg% R_Teff_win_aspect_ratio = R_Teff_win_aspect_ratio
         s% pg% R_Teff_xleft = R_Teff_xleft
         s% pg% R_Teff_xright = R_Teff_xright
         s% pg% R_Teff_ybot = R_Teff_ybot
         s% pg% R_Teff_ytop = R_Teff_ytop
         s% pg% R_Teff_txt_scale = R_Teff_txt_scale
         s% pg% R_Teff_file_width = R_Teff_file_width
         s% pg% R_Teff_file_aspect_ratio = R_Teff_file_aspect_ratio
         s% pg% R_Teff_use_decorator = R_Teff_use_decorator

         s% pg% R_L_win_flag = R_L_win_flag
         s% pg% R_L_file_flag = R_L_file_flag
         s% pg% show_R_L_target_box = show_R_L_target_box
         s% pg% R_L_target_n_sigma = R_L_target_n_sigma
         s% pg% R_L_target_R = R_L_target_R
         s% pg% R_L_target_R_sigma = R_L_target_R_sigma
         s% pg% R_L_target_L = R_L_target_L
         s% pg% R_L_target_L_sigma = R_L_target_L_sigma
         s% pg% R_L_file_interval = R_L_file_interval
         s% pg% R_L_step_min = R_L_step_min
         s% pg% R_L_step_max = R_L_step_max
         s% pg% R_L_file_dir = R_L_file_dir
         s% pg% R_L_file_prefix = R_L_file_prefix
         s% pg% show_R_L_annotation1 = show_R_L_annotation1
         s% pg% show_R_L_annotation2 = show_R_L_annotation2
         s% pg% show_R_L_annotation3 = show_R_L_annotation3
         s% pg% R_L_fname = R_L_fname
         s% pg% R_L_title = R_L_title
         s% pg% R_L_R_min = R_L_R_min
         s% pg% R_L_R_max = R_L_R_max
         s% pg% R_L_L_min = R_L_L_min
         s% pg% R_L_L_max = R_L_L_max
         s% pg% R_L_L_margin = R_L_L_margin
         s% pg% R_L_R_margin = R_L_R_margin
         s% pg% R_L_dL_min = R_L_dL_min
         s% pg% R_L_dR_min = R_L_dR_min
         s% pg% R_L_win_width = R_L_win_width
         s% pg% R_L_win_aspect_ratio = R_L_win_aspect_ratio
         s% pg% R_L_xleft = R_L_xleft
         s% pg% R_L_xright = R_L_xright
         s% pg% R_L_ybot = R_L_ybot
         s% pg% R_L_ytop = R_L_ytop
         s% pg% R_L_txt_scale = R_L_txt_scale
         s% pg% R_L_file_width = R_L_file_width
         s% pg% R_L_file_aspect_ratio = R_L_file_aspect_ratio
         s% pg% R_L_use_decorator = R_L_use_decorator

         s% pg% logg_logT_win_flag = logg_logT_win_flag
         s% pg% logg_logT_file_flag = logg_logT_file_flag
         s% pg% show_logg_logT_target_box = show_logg_logT_target_box
         s% pg% logg_logT_target_n_sigma = logg_logT_target_n_sigma
         s% pg% logg_logT_target_logg = logg_logT_target_logg
         s% pg% logg_logT_target_logg_sigma = logg_logT_target_logg_sigma
         s% pg% logg_logT_target_logT = logg_logT_target_logT
         s% pg% logg_logT_target_logT_sigma = logg_logT_target_logT_sigma
         s% pg% logg_logT_file_interval = logg_logT_file_interval
         s% pg% logg_logT_step_min = logg_logT_step_min
         s% pg% logg_logT_step_max = logg_logT_step_max
         s% pg% logg_logT_file_dir = logg_logT_file_dir
         s% pg% logg_logT_file_prefix = logg_logT_file_prefix
         s% pg% show_logg_logT_annotation1 = show_logg_logT_annotation1
         s% pg% show_logg_logT_annotation2 = show_logg_logT_annotation2
         s% pg% show_logg_logT_annotation3 = show_logg_logT_annotation3
         s% pg% logg_logT_fname = logg_logT_fname
         s% pg% logg_logT_logg_min = logg_logT_logg_min
         s% pg% logg_logT_logg_max = logg_logT_logg_max
         s% pg% logg_logT_logT_min = logg_logT_logT_min
         s% pg% logg_logT_logT_max = logg_logT_logT_max
         s% pg% logg_logT_logg_margin = logg_logT_logg_margin
         s% pg% logg_logT_logT_margin = logg_logT_logT_margin
         s% pg% logg_logT_dlogT_min = logg_logT_dlogT_min
         s% pg% logg_logT_dlogg_min = logg_logT_dlogg_min
         s% pg% logg_logT_win_width = logg_logT_win_width
         s% pg% logg_logT_win_aspect_ratio = logg_logT_win_aspect_ratio
         s% pg% logg_logT_file_width = logg_logT_file_width
         s% pg% logg_logT_file_aspect_ratio = logg_logT_file_aspect_ratio
         s% pg% logg_logT_xleft = logg_logT_xleft
         s% pg% logg_logT_xright = logg_logT_xright
         s% pg% logg_logT_ybot = logg_logT_ybot
         s% pg% logg_logT_ytop = logg_logT_ytop
         s% pg% logg_logT_txt_scale = logg_logT_txt_scale
         s% pg% logg_logT_title = logg_logT_title
         s% pg% logg_logT_use_decorator = logg_LogT_use_decorator

         s% pg% HR_win_flag = HR_win_flag
         s% pg% HR_file_flag = HR_file_flag
         s% pg% HR_file_interval = HR_file_interval
         s% pg% HR_step_min = HR_step_min
         s% pg% HR_step_max = HR_step_max
         s% pg% show_HR_classical_instability_strip = show_HR_classical_instability_strip
         s% pg% show_HR_Mira_instability_region = show_HR_Mira_instability_region
         s% pg% show_HR_WD_instabilities = show_HR_WD_instabilities
         s% pg% show_HR_target_box = show_HR_target_box
         s% pg% HR_target_n_sigma = HR_target_n_sigma
         s% pg% HR_target_logL = HR_target_logL
         s% pg% HR_target_logL_sigma = HR_target_logL_sigma
         s% pg% HR_target_logT = HR_target_logT
         s% pg% HR_target_logT_sigma = HR_target_logT_sigma
         s% pg% HR_file_dir = HR_file_dir
         s% pg% HR_file_prefix = HR_file_prefix
         s% pg% show_HR_annotation1 = show_HR_annotation1
         s% pg% show_HR_annotation2 = show_HR_annotation2
         s% pg% show_HR_annotation3 = show_HR_annotation3
         s% pg% HR_fname = HR_fname
         s% pg% HR_logT_min = HR_logT_min
         s% pg% HR_logT_max = HR_logT_max
         s% pg% HR_logL_min = HR_logL_min
         s% pg% HR_logL_max = HR_logL_max
         s% pg% HR_logL_margin = HR_logL_margin
         s% pg% HR_logT_margin = HR_logT_margin
         s% pg% HR_dlogT_min = HR_dlogT_min
         s% pg% HR_dlogL_min = HR_dlogL_min
         s% pg% HR_win_width = HR_win_width
         s% pg% HR_win_aspect_ratio = HR_win_aspect_ratio
         s% pg% HR_file_width = HR_file_width
         s% pg% HR_file_aspect_ratio = HR_file_aspect_ratio
         s% pg% HR_xleft = HR_xleft
         s% pg% HR_xright = HR_xright
         s% pg% HR_ybot = HR_ybot
         s% pg% HR_ytop = HR_ytop
         s% pg% HR_txt_scale = HR_txt_scale
         s% pg% HR_title = HR_title
         s% pg% HR_use_decorator = HR_use_decorator

         s% pg% TRho_win_flag = TRho_win_flag
         s% pg% TRho_file_flag = TRho_file_flag
         s% pg% TRho_file_interval = TRho_file_interval
         s% pg% TRho_step_max = TRho_step_max
         s% pg% TRho_step_min = TRho_step_min
         s% pg% TRho_file_dir = TRho_file_dir
         s% pg% TRho_file_prefix = TRho_file_prefix
         s% pg% show_TRho_annotation1 = show_TRho_annotation1
         s% pg% show_TRho_annotation2 = show_TRho_annotation2
         s% pg% show_TRho_annotation3 = show_TRho_annotation3
         s% pg% show_TRho_degeneracy_line = show_TRho_degeneracy_line
         s% pg% TRho_fname = TRho_fname
         s% pg% TRho_logT_min = TRho_logT_min
         s% pg% TRho_logT_max = TRho_logT_max
         s% pg% TRho_logRho_min = TRho_logRho_min
         s% pg% TRho_logRho_max = TRho_logRho_max
         s% pg% TRho_logT_margin = TRho_logT_margin
         s% pg% TRho_logRho_margin = TRho_logRho_margin
         s% pg% TRho_logRho_dlogRho_min = TRho_logRho_dlogRho_min
         s% pg% TRho_logT_dlogT_min = TRho_logT_dlogT_min
         s% pg% TRho_win_width = TRho_win_width
         s% pg% TRho_win_aspect_ratio = TRho_win_aspect_ratio
         s% pg% TRho_file_width = TRho_file_width
         s% pg% TRho_file_aspect_ratio = TRho_file_aspect_ratio
         s% pg% TRho_xleft = TRho_xleft
         s% pg% TRho_xright = TRho_xright
         s% pg% TRho_ybot = TRho_ybot
         s% pg% TRho_ytop = TRho_ytop
         s% pg% TRho_txt_scale = TRho_txt_scale
         s% pg% TRho_title = TRho_title
         s% pg% TRho_use_decorator = TRho_use_decorator

         s% pg% TmaxRho_win_flag = TmaxRho_win_flag
         s% pg% TmaxRho_file_flag = TmaxRho_file_flag
         s% pg% TmaxRho_file_interval = TmaxRho_file_interval
         s% pg% TmaxRho_step_max = TmaxRho_step_max
         s% pg% TmaxRho_step_min = TmaxRho_step_min
         s% pg% TmaxRho_file_dir = TmaxRho_file_dir
         s% pg% TmaxRho_file_prefix = TmaxRho_file_prefix
         s% pg% show_TmaxRho_annotation1 = show_TmaxRho_annotation1
         s% pg% show_TmaxRho_annotation2 = show_TmaxRho_annotation2
         s% pg% show_TmaxRho_annotation3 = show_TmaxRho_annotation3
         s% pg% show_TmaxRho_degeneracy_line = show_TmaxRho_degeneracy_line
         s% pg% TmaxRho_fname = TmaxRho_fname
         s% pg% TmaxRho_logT_min = TmaxRho_logT_min
         s% pg% TmaxRho_logT_max = TmaxRho_logT_max
         s% pg% TmaxRho_logRho_min = TmaxRho_logRho_min
         s% pg% TmaxRho_logRho_max = TmaxRho_logRho_max
         s% pg% TmaxRho_logT_margin = TmaxRho_logT_margin
         s% pg% TmaxRho_logRho_margin = TmaxRho_logRho_margin
         s% pg% TmaxRho_logRho_dlogRho_min = TmaxRho_logRho_dlogRho_min
         s% pg% TmaxRho_logT_dlogT_min = TmaxRho_logT_dlogT_min
         s% pg% TmaxRho_win_width = TmaxRho_win_width
         s% pg% TmaxRho_win_aspect_ratio = TmaxRho_win_aspect_ratio
         s% pg% TmaxRho_file_width = TmaxRho_file_width
         s% pg% TmaxRho_file_aspect_ratio = TmaxRho_file_aspect_ratio
         s% pg% TmaxRho_xleft = TmaxRho_xleft
         s% pg% TmaxRho_xright = TmaxRho_xright
         s% pg% TmaxRho_ybot = TmaxRho_ybot
         s% pg% TmaxRho_ytop = TmaxRho_ytop
         s% pg% TmaxRho_txt_scale = TmaxRho_txt_scale
         s% pg% TmaxRho_title = TmaxRho_title
         s% pg% TmaxRho_use_decorator = TmaxRho_use_decorator

         s% pg% Dynamo_win_flag = Dynamo_win_flag
         s% pg% Dynamo_file_flag = Dynamo_file_flag
         s% pg% Dynamo_file_interval = Dynamo_file_interval
         s% pg% Dynamo_file_dir = Dynamo_file_dir
         s% pg% Dynamo_file_prefix = Dynamo_file_prefix
         s% pg% show_Dynamo_annotation1 = show_Dynamo_annotation1
         s% pg% show_Dynamo_annotation2 = show_Dynamo_annotation2
         s% pg% show_Dynamo_annotation3 = show_Dynamo_annotation3
         s% pg% Dynamo_xaxis_name = Dynamo_xaxis_name
         s% pg% Dynamo_xaxis_reversed = Dynamo_xaxis_reversed
         s% pg% Dynamo_xmin = Dynamo_xmin
         s% pg% Dynamo_xmax = Dynamo_xmax
         s% pg% Dynamo_ymin_left = Dynamo_ymin_left
         s% pg% Dynamo_ymax_left = Dynamo_ymax_left
         s% pg% Dynamo_dymin_left = Dynamo_dymin_left
         s% pg% Dynamo_ymin_right = Dynamo_ymin_right
         s% pg% Dynamo_ymax_right = Dynamo_ymax_right
         s% pg% Dynamo_dymin_right = Dynamo_dymin_right
         s% pg% Dynamo_win_width = Dynamo_win_width
         s% pg% Dynamo_win_aspect_ratio = Dynamo_win_aspect_ratio
         s% pg% Dynamo_file_width = Dynamo_file_width
         s% pg% Dynamo_file_aspect_ratio = Dynamo_file_aspect_ratio
         s% pg% Dynamo_xleft = Dynamo_xleft
         s% pg% Dynamo_xright = Dynamo_xright
         s% pg% Dynamo_ybot = Dynamo_ybot
         s% pg% Dynamo_ytop = Dynamo_ytop
         s% pg% Dynamo_txt_scale = Dynamo_txt_scale
         s% pg% Dynamo_title = Dynamo_title
         s% pg% Dynamo_legend_txt_scale_factor = Dynamo_legend_txt_scale_factor
         s% pg% Dynamo_use_decorator = Dynamo_use_decorator

         s% pg% Mixing_win_flag = Mixing_win_flag
         s% pg% Mixing_file_flag = Mixing_file_flag
         s% pg% Mixing_file_interval = Mixing_file_interval
         s% pg% Mixing_file_dir = Mixing_file_dir
         s% pg% Mixing_file_prefix = Mixing_file_prefix
         s% pg% show_Mixing_annotation1 = show_Mixing_annotation1
         s% pg% show_Mixing_annotation2 = show_Mixing_annotation2
         s% pg% show_Mixing_annotation3 = show_Mixing_annotation3
         s% pg% Mixing_xaxis_name = Mixing_xaxis_name
         s% pg% Mixing_xaxis_reversed = Mixing_xaxis_reversed
         s% pg% Mixing_legend_txt_scale_factor = Mixing_legend_txt_scale_factor
         s% pg% Mixing_show_rotation_details = Mixing_show_rotation_details
         s% pg% Mixing_xmin = Mixing_xmin
         s% pg% Mixing_xmax = Mixing_xmax
         s% pg% Mixing_ymin = Mixing_ymin
         s% pg% Mixing_ymax = Mixing_ymax
         s% pg% Mixing_dymin = Mixing_dymin
         s% pg% Mixing_win_width = Mixing_win_width
         s% pg% Mixing_win_aspect_ratio = Mixing_win_aspect_ratio
         s% pg% Mixing_file_width = Mixing_file_width
         s% pg% Mixing_file_aspect_ratio = Mixing_file_aspect_ratio
         s% pg% Mixing_xleft = Mixing_xleft
         s% pg% Mixing_xright = Mixing_xright
         s% pg% Mixing_ybot = Mixing_ybot
         s% pg% Mixing_ytop = Mixing_ytop
         s% pg% Mixing_txt_scale = Mixing_txt_scale
         s% pg% Mixing_title = Mixing_title
         s% pg% Mixing_use_decorator = Mixing_use_decorator

         s% pg% Network_win_flag = Network_win_flag
         s% pg% Network_file_flag = Network_file_flag
         s% pg% Network_file_interval = Network_file_interval
         s% pg% Network_file_dir = Network_file_dir
         s% pg% Network_file_prefix = Network_file_prefix
         s% pg% Network_nmin = Network_nmin
         s% pg% Network_nmax = Network_nmax
         s% pg% Network_zmin = Network_zmin
         s% pg% Network_zmax = Network_zmax
         s% pg% Network_show_mass_fraction = Network_show_mass_fraction
         s% pg% Network_log_mass_frac_min = Network_log_mass_frac_min
         s% pg% Network_log_mass_frac_max = Network_log_mass_frac_max
         s% pg% Network_show_element_names = Network_show_element_names
         s% pg% Network_win_width = Network_win_width
         s% pg% Network_win_aspect_ratio = Network_win_aspect_ratio
         s% pg% Network_file_width = Network_file_width
         s% pg% Network_file_aspect_ratio = Network_file_aspect_ratio
         s% pg% Network_xleft = Network_xleft
         s% pg% Network_xright = Network_xright
         s% pg% Network_ybot = Network_ybot
         s% pg% Network_ytop = Network_ytop
         s% pg% Network_txt_scale = Network_txt_scale
         s% pg% Network_title = Network_title
         s% pg% Network_use_decorator = Network_use_decorator
         s% pg% Network_show_colorbar = Network_show_colorbar

         s% pg% Production_win_flag = Production_win_flag
         s% pg% Production_file_flag = Production_file_flag
         s% pg% Production_file_interval = Production_file_interval
         s% pg% Production_file_dir = Production_file_dir
         s% pg% Production_file_prefix = Production_file_prefix
         s% pg% Production_amin = Production_amin
         s% pg% Production_amax = Production_amax
         s% pg% Production_ymin = Production_ymin
         s% pg% Production_ymax = Production_ymax
         s% pg% Production_min_mass = Production_min_mass
         s% pg% Production_max_mass = Production_max_mass
         s% pg% Production_min_mass_frac = Production_min_mass_frac
         s% pg% Production_show_element_names = Production_show_element_names
         s% pg% Production_win_width = Production_win_width
         s% pg% Production_win_aspect_ratio = Production_win_aspect_ratio
         s% pg% Production_file_width = Production_file_width
         s% pg% Production_file_aspect_ratio = Production_file_aspect_ratio
         s% pg% Production_xleft = Production_xleft
         s% pg% Production_xright = Production_xright
         s% pg% Production_ybot = Production_ybot
         s% pg% Production_ytop = Production_ytop
         s% pg% Production_txt_scale = Production_txt_scale
         s% pg% Production_title = Production_title
         s% pg% Production_use_decorator = Production_use_decorator


         s% pg% History_Track_win_flag = History_Track_win_flag
         s% pg% History_Track_file_flag = History_Track_file_flag
         s% pg% History_Track_file_interval = History_Track_file_interval
         s% pg% History_Track_step_min = History_Track_step_min
         s% pg% History_Track_step_max = History_Track_step_max
         s% pg% show_History_Track_target_box = show_History_Track_target_box
         s% pg% History_Track_n_sigma = History_Track_n_sigma
         s% pg% History_Track_xtarget = History_Track_xtarget
         s% pg% History_Track_xsigma = History_Track_xsigma
         s% pg% History_Track_ytarget = History_Track_ytarget
         s% pg% History_Track_ysigma = History_Track_ysigma
         s% pg% History_Track_file_dir = History_Track_file_dir
         s% pg% History_Track_file_prefix = History_Track_file_prefix
         s% pg% show_History_Track_annotation1 = show_History_Track_annotation1
         s% pg% show_History_Track_annotation2 = show_History_Track_annotation2
         s% pg% show_History_Track_annotation3 = show_History_Track_annotation3
         s% pg% History_Track_fname = History_Track_fname
         s% pg% History_Track_xname = History_Track_xname
         s% pg% History_Track_xaxis_label = History_Track_xaxis_label
         s% pg% History_Track_yname = History_Track_yname
         s% pg% History_Track_yaxis_label = History_Track_yaxis_label
         s% pg% History_Track_reverse_xaxis = History_Track_reverse_xaxis
         s% pg% History_Track_reverse_yaxis = History_Track_reverse_yaxis
         s% pg% History_Track_log_xaxis = History_Track_log_xaxis
         s% pg% History_Track_log_yaxis = History_Track_log_yaxis
         s% pg% History_Track_xmin = History_Track_xmin
         s% pg% History_Track_xmax = History_Track_xmax
         s% pg% History_Track_ymin = History_Track_ymin
         s% pg% History_Track_ymax = History_Track_ymax
         s% pg% History_Track_xmargin = History_Track_xmargin
         s% pg% History_Track_ymargin = History_Track_ymargin
         s% pg% History_Track_dxmin = History_Track_dxmin
         s% pg% History_Track_dymin = History_Track_dymin
         s% pg% History_Track_win_width = History_Track_win_width
         s% pg% History_Track_win_aspect_ratio = History_Track_win_aspect_ratio
         s% pg% History_Track_file_width = History_Track_file_width
         s% pg% History_Track_file_aspect_ratio = History_Track_file_aspect_ratio
         s% pg% History_Track_xleft = History_Track_xleft
         s% pg% History_Track_xright = History_Track_xright
         s% pg% History_Track_ybot = History_Track_ybot
         s% pg% History_Track_ytop = History_Track_ytop
         s% pg% History_Track_txt_scale = History_Track_txt_scale
         s% pg% History_Track_title = History_Track_title
         s% pg% History_Track_use_decorator = History_Track_use_decorator

         s% pg% Kipp_win_flag = Kipp_win_flag
         s% pg% Kipp_file_flag = Kipp_file_flag
         s% pg% Kipp_file_interval = Kipp_file_interval
         s% pg% Kipp_xmax = Kipp_xmax
         s% pg% Kipp_xmin = Kipp_xmin
         s% pg% Kipp_max_width = Kipp_max_width
         s% pg% Kipp_step_xmax = Kipp_step_xmax
         s% pg% Kipp_step_xmin = Kipp_step_xmin
         s% pg% Kipp_xaxis_name = Kipp_xaxis_name
         s% pg% Kipp_xaxis_log = Kipp_xaxis_log
         s% pg% Kipp_xmargin = Kipp_xmargin
         s% pg% Kipp_xaxis_reversed = Kipp_xaxis_reversed
         s% pg% Kipp_xaxis_in_seconds = Kipp_xaxis_in_seconds
         s% pg% Kipp_xaxis_in_Myr = Kipp_xaxis_in_Myr
         s% pg% Kipp_xaxis_time_from_present = Kipp_xaxis_time_from_present
         s% pg% Kipp_file_dir = Kipp_file_dir
         s% pg% Kipp_file_prefix = Kipp_file_prefix
         s% pg% show_Kipp_annotation1 = show_Kipp_annotation1
         s% pg% show_Kipp_annotation2 = show_Kipp_annotation2
         s% pg% show_Kipp_annotation3 = show_Kipp_annotation3
         s% pg% Kipp_show_burn = Kipp_show_burn
         s% pg% Kipp_show_mixing = Kipp_show_mixing
         s% pg% Kipp_show_luminosities = Kipp_show_luminosities
         s% pg% Kipp_show_mass_boundaries = Kipp_show_mass_boundaries
         s% pg% Kipp_mix_line_weight = Kipp_mix_line_weight
         s% pg% Kipp_mix_interval = Kipp_mix_interval
         s% pg% Kipp_burn_line_weight = Kipp_burn_line_weight
         s% pg% Kipp_burn_type_cutoff = Kipp_burn_type_cutoff
         s% pg% Kipp_luminosities_line_weight = Kipp_luminosities_line_weight
         s% pg% Kipp_masses_line_weight = Kipp_masses_line_weight
         s% pg% Kipp_mass_max = Kipp_mass_max
         s% pg% Kipp_mass_min = Kipp_mass_min
         s% pg% Kipp_lgL_max = Kipp_lgL_max
         s% pg% Kipp_lgL_min = Kipp_lgL_min
         s% pg% Kipp_mass_margin = Kipp_mass_margin
         s% pg% Kipp_lgL_margin = Kipp_lgL_margin
         s% pg% Kipp_win_width = Kipp_win_width
         s% pg% Kipp_win_aspect_ratio = Kipp_win_aspect_ratio
         s% pg% Kipp_file_width = Kipp_file_width
         s% pg% Kipp_file_aspect_ratio = Kipp_file_aspect_ratio
         s% pg% Kipp_xleft = Kipp_xleft
         s% pg% Kipp_xright = Kipp_xright
         s% pg% Kipp_ybot = Kipp_ybot
         s% pg% Kipp_ytop = Kipp_ytop
         s% pg% Kipp_txt_scale = Kipp_txt_scale
         s% pg% Kipp_title = Kipp_title
         s% pg% kipp_use_decorator = kipp_use_decorator

         s% pg% rti_win_flag = rti_win_flag
         s% pg% rti_file_flag = rti_file_flag
         s% pg% rti_file_interval = rti_file_interval
         s% pg% rti_xmax = rti_xmax
         s% pg% rti_xmin = rti_xmin
         s% pg% rti_max_width = rti_max_width
         s% pg% rti_step_xmax = rti_step_xmax
         s% pg% rti_step_xmin = rti_step_xmin
         s% pg% rti_xaxis_name = rti_xaxis_name
         s% pg% rti_xaxis_log = rti_xaxis_log
         s% pg% rti_xmargin = rti_xmargin
         s% pg% rti_xaxis_reversed = rti_xaxis_reversed
         s% pg% rti_xaxis_in_seconds = rti_xaxis_in_seconds
         s% pg% rti_xaxis_in_Myr = rti_xaxis_in_Myr
         s% pg% rti_xaxis_time_from_present = rti_xaxis_time_from_present
         s% pg% rti_mass_max = rti_mass_max
         s% pg% rti_mass_min = rti_mass_min
         s% pg% rti_mass_margin = rti_mass_margin
         s% pg% rti_file_dir = rti_file_dir
         s% pg% rti_file_prefix = rti_file_prefix
         s% pg% show_rti_annotation1 = show_rti_annotation1
         s% pg% show_rti_annotation2 = show_rti_annotation2
         s% pg% show_rti_annotation3 = show_rti_annotation3
         s% pg% rti_line_weight = rti_line_weight
         s% pg% rti_interval = rti_interval
         s% pg% rti_win_width = rti_win_width
         s% pg% rti_win_aspect_ratio = rti_win_aspect_ratio
         s% pg% rti_file_width = rti_file_width
         s% pg% rti_file_aspect_ratio = rti_file_aspect_ratio
         s% pg% rti_xleft = rti_xleft
         s% pg% rti_xright = rti_xright
         s% pg% rti_ybot = rti_ybot
         s% pg% rti_ytop = rti_ytop
         s% pg% rti_txt_scale = rti_txt_scale
         s% pg% rti_title = rti_title
         s% pg% rti_use_decorator = rti_use_decorator

         s% pg% TRho_Profile_win_flag = TRho_Profile_win_flag
         s% pg% TRho_switch_to_Column_Depth = TRho_switch_to_Column_Depth
         s% pg% TRho_switch_to_mass = TRho_switch_to_mass
         s% pg% TRho_Profile_file_flag = TRho_Profile_file_flag
         s% pg% TRho_Profile_file_interval = TRho_Profile_file_interval
         s% pg% TRho_Profile_file_dir = TRho_Profile_file_dir
         s% pg% TRho_Profile_file_prefix = TRho_Profile_file_prefix
         s% pg% show_TRho_Profile_text_info = show_TRho_Profile_text_info
         s% pg% show_TRho_Profile_legend = show_TRho_Profile_legend
         s% pg% show_TRho_Profile_mass_locs = show_TRho_Profile_mass_locs
         s% pg% show_TRho_Profile_burn_labels = show_TRho_Profile_burn_labels
         s% pg% show_TRho_Profile_kap_regions = show_TRho_Profile_kap_regions
         s% pg% show_TRho_Profile_eos_regions = show_TRho_Profile_eos_regions
         s% pg% show_TRho_Profile_gamma1_4_3rd = show_TRho_Profile_gamma1_4_3rd
         s% pg% show_TRho_Profile_burn_lines = show_TRho_Profile_burn_lines
         s% pg% show_TRho_Profile_degeneracy_line = show_TRho_Profile_degeneracy_line
         s% pg% show_TRho_Profile_Pgas_Prad_line = show_TRho_Profile_Pgas_Prad_line
         s% pg% show_TRho_Profile_annotation1 = show_TRho_Profile_annotation1
         s% pg% show_TRho_Profile_annotation2 = show_TRho_Profile_annotation2
         s% pg% show_TRho_Profile_annotation3 = show_TRho_Profile_annotation3
         s% pg% TRho_Profile_fname = TRho_Profile_fname
         s% pg% show_TRho_accretion_mesh_borders = show_TRho_accretion_mesh_borders
         s% pg% TRho_Profile_text_info_xfac = TRho_Profile_text_info_xfac
         s% pg% TRho_Profile_text_info_dxfac = TRho_Profile_text_info_dxfac
         s% pg% TRho_Profile_text_info_yfac = TRho_Profile_text_info_yfac
         s% pg% TRho_Profile_text_info_dyfac = TRho_Profile_text_info_dyfac
         s% pg% TRho_Profile_xmin = TRho_Profile_xmin
         s% pg% TRho_Profile_xmax = TRho_Profile_xmax
         s% pg% TRho_Profile_ymin = TRho_Profile_ymin
         s% pg% TRho_Profile_ymax = TRho_Profile_ymax
         s% pg% TRho_Profile_legend_coord = TRho_Profile_legend_coord
         s% pg% TRho_Profile_legend_fjust = TRho_Profile_legend_fjust
         s% pg% TRho_Profile_legend_disp1 = TRho_Profile_legend_disp1
         s% pg% TRho_Profile_legend_del_disp = TRho_Profile_legend_del_disp
         s% pg% TRho_Profile_legend_txt_scale = TRho_Profile_legend_txt_scale
         s% pg% TRho_Profile_win_width = TRho_Profile_win_width
         s% pg% TRho_Profile_win_aspect_ratio = TRho_Profile_win_aspect_ratio
         s% pg% TRho_Profile_xleft = TRho_Profile_xleft
         s% pg% TRho_Profile_xright = TRho_Profile_xright
         s% pg% TRho_Profile_ybot = TRho_Profile_ybot
         s% pg% TRho_Profile_ytop = TRho_Profile_ytop
         s% pg% TRho_Profile_txt_scale = TRho_Profile_txt_scale
         s% pg% TRho_Profile_file_width = TRho_Profile_file_width
         s% pg% TRho_Profile_file_aspect_ratio = TRho_Profile_file_aspect_ratio
         s% pg% TRho_Profile_title = TRho_Profile_title
         s% pg% num_profile_mass_points = num_profile_mass_points
         s% pg% profile_mass_point_q = profile_mass_point_q
         s% pg% profile_mass_point_color_index = profile_mass_point_color_index
         s% pg% profile_mass_point_symbol = profile_mass_point_symbol
         s% pg% profile_mass_point_symbol_scale = profile_mass_point_symbol_scale
         s% pg% profile_mass_point_str = profile_mass_point_str
         s% pg% profile_mass_point_str_clr = profile_mass_point_str_clr
         s% pg% profile_mass_point_str_scale = profile_mass_point_str_scale
         s% pg% TRho_Profile_use_decorator = TRho_Profile_use_decorator

         s% pg% History_Panels_win_flag = History_Panels_win_flag
         s% pg% History_Panels_win_width = History_Panels_win_width
         s% pg% History_Panels_win_aspect_ratio = History_Panels_win_aspect_ratio
         s% pg% History_Panels_xleft = History_Panels_xleft
         s% pg% History_Panels_xright = History_Panels_xright
         s% pg% History_Panels_ybot = History_Panels_ybot
         s% pg% History_Panels_ytop = History_Panels_ytop
         s% pg% History_Panels_txt_scale = History_Panels_txt_scale
         s% pg% History_Panels_title = History_Panels_title
         s% pg% History_Panels_xmax = History_Panels_xmax
         s% pg% History_Panels_xmin = History_Panels_xmin
         s% pg% History_Panels_dxmin = History_Panels_dxmin
         s% pg% History_Panels_max_width = History_Panels_max_width
         s% pg% History_Panels_num_panels = History_Panels_num_panels
         s% pg% History_Panels_xaxis_name = History_Panels_xaxis_name
         s% pg% History_Panels_automatic_star_age_units = History_Panels_automatic_star_age_units
         s% pg% History_Panels_yaxis_name = History_Panels_yaxis_name
         s% pg% History_Panels_xaxis_reversed = History_Panels_xaxis_reversed
         s% pg% History_Panels_yaxis_reversed = History_Panels_yaxis_reversed
         s% pg% History_Panels_xaxis_log = History_Panels_xaxis_log
         s% pg% History_Panels_yaxis_log = History_Panels_yaxis_log
         s% pg% History_Panels_ymin = History_Panels_ymin
         s% pg% History_Panels_ymax = History_Panels_ymax
         s% pg% History_Panels_dymin = History_Panels_dymin
         s% pg% History_Panels_other_yaxis_name = History_Panels_other_yaxis_name
         s% pg% History_Panels_other_yaxis_reversed = History_Panels_other_yaxis_reversed
         s% pg% History_Panels_other_yaxis_log = History_Panels_other_yaxis_log
         s% pg% History_Panels_same_yaxis_range = History_Panels_same_yaxis_range
         s% pg% History_Panels_other_ymin = History_Panels_other_ymin
         s% pg% History_Panels_other_ymax = History_Panels_other_ymax
         s% pg% History_Panels_other_dymin = History_Panels_other_dymin
         s% pg% History_Panels_file_flag = History_Panels_file_flag
         s% pg% History_Panels_points_name = History_Panels_points_name
         s% pg% History_Panels_file_dir = History_Panels_file_dir
         s% pg% History_Panels_file_prefix = History_Panels_file_prefix
         s% pg% History_Panels_file_interval = History_Panels_file_interval
         s% pg% History_Panels_file_width = History_Panels_file_width
         s% pg% History_Panels_file_aspect_ratio = History_Panels_file_aspect_ratio
         s% pg% History_Panels_xmargin = History_Panels_xmargin
         s% pg% History_Panels_ymargin = History_Panels_ymargin
         s% pg% History_Panels_other_ymargin = History_Panels_other_ymargin
         s% pg% History_Panels_use_decorator = History_Panels_use_decorator

         s% pg% History_Panel_points_error_bars = History_Panel_points_error_bars
         s% pg% History_Panel_points_interval = History_Panel_points_interval
         s% pg% History_Panel_points_marker = History_Panel_points_marker
         s% pg% History_Panel_points_ci = History_Panel_points_ci
         s% pg% History_Panel_points_lw = History_Panel_points_lw
         s% pg% History_Panel_points_ch = History_Panel_points_ch

         s% pg% Color_Magnitude_win_flag = Color_Magnitude_win_flag
         s% pg% Color_Magnitude_win_width = Color_Magnitude_win_width
         s% pg% Color_Magnitude_win_aspect_ratio = Color_Magnitude_win_aspect_ratio
         s% pg% Color_Magnitude_xleft = Color_Magnitude_xleft
         s% pg% Color_Magnitude_xright = Color_Magnitude_xright
         s% pg% Color_Magnitude_ybot = Color_Magnitude_ybot
         s% pg% Color_Magnitude_ytop = Color_Magnitude_ytop
         s% pg% Color_Magnitude_txt_scale = Color_Magnitude_txt_scale
         s% pg% Color_Magnitude_title = Color_Magnitude_title
         s% pg% Color_Magnitude_xmax = Color_Magnitude_xmax
         s% pg% Color_Magnitude_xmin = Color_Magnitude_xmin
         s% pg% Color_Magnitude_dxmin = Color_Magnitude_dxmin
         s% pg% Color_Magnitude_max_width = Color_Magnitude_max_width
         s% pg% Color_Magnitude_num_panels = Color_Magnitude_num_panels
         s% pg% Color_Magnitude_xaxis1_name = Color_Magnitude_xaxis1_name
         s% pg% Color_Magnitude_xaxis2_name = Color_Magnitude_xaxis2_name
         s% pg% Color_Magnitude_yaxis1_name = Color_Magnitude_yaxis1_name
         s% pg% Color_Magnitude_yaxis2_name = Color_Magnitude_yaxis2_name
         s% pg% Color_Magnitude_xaxis_reversed = Color_Magnitude_xaxis_reversed
         s% pg% Color_Magnitude_yaxis_reversed = Color_Magnitude_yaxis_reversed
         s% pg% Color_Magnitude_xaxis_log = Color_Magnitude_xaxis_log
         s% pg% Color_Magnitude_yaxis_log = Color_Magnitude_yaxis_log
         s% pg% Color_Magnitude_ymin = Color_Magnitude_ymin
         s% pg% Color_Magnitude_ymax = Color_Magnitude_ymax
         s% pg% Color_Magnitude_dymin = Color_Magnitude_dymin
         s% pg% Color_Magnitude_other_yaxis1_name = Color_Magnitude_other_yaxis1_name
         s% pg% Color_Magnitude_other_yaxis2_name = Color_Magnitude_other_yaxis2_name
         s% pg% Color_Magnitude_other_yaxis_reversed = Color_Magnitude_other_yaxis_reversed
         s% pg% Color_Magnitude_other_yaxis_log = Color_Magnitude_other_yaxis_log
         s% pg% Color_Magnitude_other_ymin = Color_Magnitude_other_ymin
         s% pg% Color_Magnitude_other_ymax = Color_Magnitude_other_ymax
         s% pg% Color_Magnitude_other_dymin = Color_Magnitude_other_dymin
         s% pg% Color_Magnitude_file_flag = Color_Magnitude_file_flag
         s% pg% Color_Magnitude_file_dir = Color_Magnitude_file_dir
         s% pg% Color_Magnitude_file_prefix = Color_Magnitude_file_prefix
         s% pg% Color_Magnitude_file_interval = Color_Magnitude_file_interval
         s% pg% Color_Magnitude_file_width = Color_Magnitude_file_width
         s% pg% Color_Magnitude_file_aspect_ratio = Color_Magnitude_file_aspect_ratio
         s% pg% Color_Magnitude_xmargin = Color_Magnitude_xmargin
         s% pg% Color_Magnitude_ymargin = Color_Magnitude_ymargin
         s% pg% Color_Magnitude_other_ymargin = Color_Magnitude_other_ymargin
         s% pg% Color_Magnitude_use_decorator = Color_Magnitude_use_decorator

         s% pg% Mode_Prop_win_flag = Mode_Prop_win_flag
         s% pg% Mode_Prop_file_flag = Mode_Prop_file_flag
         s% pg% Mode_Prop_file_interval = Mode_Prop_file_interval
         s% pg% Mode_Prop_file_dir = Mode_Prop_file_dir
         s% pg% Mode_Prop_file_prefix = Mode_Prop_file_prefix
         s% pg% Mode_Prop_xaxis_name = Mode_Prop_xaxis_name
         s% pg% Mode_Prop_xaxis_reversed = Mode_Prop_xaxis_reversed
         s% pg% Mode_Prop_nu_max_obs = Mode_Prop_nu_max_obs
         s% pg% Mode_Prop_xmin = Mode_Prop_xmin
         s% pg% Mode_Prop_xmax = Mode_Prop_xmax
         s% pg% Mode_Prop_ymin = Mode_Prop_ymin
         s% pg% Mode_Prop_ymax = Mode_Prop_ymax
         s% pg% Mode_Prop_win_width = Mode_Prop_win_width
         s% pg% Mode_Prop_win_aspect_ratio = Mode_Prop_win_aspect_ratio
         s% pg% Mode_Prop_file_width = Mode_Prop_file_width
         s% pg% Mode_Prop_file_aspect_ratio = Mode_Prop_file_aspect_ratio
         s% pg% Mode_Prop_xleft = Mode_Prop_xleft
         s% pg% Mode_Prop_xright = Mode_Prop_xright
         s% pg% Mode_Prop_ybot = Mode_Prop_ybot
         s% pg% Mode_Prop_ytop = Mode_Prop_ytop
         s% pg% Mode_Prop_txt_scale = Mode_Prop_txt_scale
         s% pg% Mode_Prop_title = Mode_Prop_title
         s% pg% Mode_Prop_use_decorator = Mode_Prop_use_decorator

         s% pg% Summary_Burn_win_flag = Summary_Burn_win_flag
         s% pg% Summary_Burn_file_flag = Summary_Burn_file_flag
         s% pg% Summary_Burn_file_interval = Summary_Burn_file_interval
         s% pg% Summary_Burn_file_dir = Summary_Burn_file_dir
         s% pg% Summary_Burn_file_prefix = Summary_Burn_file_prefix
         s% pg% Summary_Burn_xaxis_name = Summary_Burn_xaxis_name
         s% pg% Summary_Burn_xaxis_reversed = Summary_Burn_xaxis_reversed
         s% pg% Summary_Burn_xmin = Summary_Burn_xmin
         s% pg% Summary_Burn_xmax = Summary_Burn_xmax
         s% pg% Summary_Burn_win_width = Summary_Burn_win_width
         s% pg% Summary_Burn_win_aspect_ratio = Summary_Burn_win_aspect_ratio
         s% pg% Summary_Burn_file_width = Summary_Burn_file_width
         s% pg% Summary_Burn_file_aspect_ratio = Summary_Burn_file_aspect_ratio
         s% pg% Summary_Burn_xleft = Summary_Burn_xleft
         s% pg% Summary_Burn_xright = Summary_Burn_xright
         s% pg% Summary_Burn_ybot = Summary_Burn_ybot
         s% pg% Summary_Burn_ytop = Summary_Burn_ytop
         s% pg% Summary_Burn_txt_scale = Summary_Burn_txt_scale
         s% pg% Summary_Burn_title = Summary_Burn_title
         s% pg% Summary_Burn_title_shift = Summary_Burn_title_shift
         s% pg% Summary_Burn_use_decorator = Summary_Burn_use_decorator
         
         s% pg% Summary_Profile_win_flag = Summary_Profile_win_flag
         s% pg% Summary_Profile_file_flag = Summary_Profile_file_flag
         s% pg% Summary_Profile_file_interval = Summary_Profile_file_interval
         s% pg% Summary_Profile_file_dir = Summary_Profile_file_dir
         s% pg% Summary_Profile_file_prefix = Summary_Profile_file_prefix
         s% pg% Summary_Profile_xaxis_name = Summary_Profile_xaxis_name
         s% pg% Summary_Profile_xaxis_reversed = Summary_Profile_xaxis_reversed
         s% pg% Summary_Profile_scaled_value = Summary_Profile_scaled_value
         s% pg% Summary_Profile_xmin = Summary_Profile_xmin
         s% pg% Summary_Profile_xmax = Summary_Profile_xmax
         s% pg% Summary_Profile_win_width = Summary_Profile_win_width
         s% pg% Summary_Profile_win_aspect_ratio = Summary_Profile_win_aspect_ratio
         s% pg% Summary_Profile_file_width = Summary_Profile_file_width
         s% pg% Summary_Profile_file_aspect_ratio = Summary_Profile_file_aspect_ratio
         s% pg% Summary_Profile_xleft = Summary_Profile_xleft
         s% pg% Summary_Profile_xright = Summary_Profile_xright
         s% pg% Summary_Profile_ybot = Summary_Profile_ybot
         s% pg% Summary_Profile_ytop = Summary_Profile_ytop
         s% pg% Summary_Profile_txt_scale = Summary_Profile_txt_scale
         s% pg% Summary_Profile_title = Summary_Profile_title
         s% pg% Summary_Profile_name = Summary_Profile_name
         s% pg% Summary_Profile_legend = Summary_Profile_legend
         s% pg% Summary_Profile_num_lines = Summary_Profile_num_lines
         s% pg% Summary_Profile_use_decorator = Summary_Profile_use_decorator

         s% pg% Summary_History_win_flag = Summary_History_win_flag
         s% pg% Summary_History_file_flag = Summary_History_file_flag
         s% pg% Summary_History_file_interval = Summary_History_file_interval
         s% pg% Summary_History_file_dir = Summary_History_file_dir
         s% pg% Summary_History_file_prefix = Summary_History_file_prefix
         s% pg% Summary_History_scaled_value = Summary_History_scaled_value
         s% pg% Summary_History_xmin = Summary_History_xmin
         s% pg% Summary_History_xmax = Summary_History_xmax
         s% pg% Summary_History_max_width = Summary_History_max_width
         s% pg% Summary_History_win_width = Summary_History_win_width
         s% pg% Summary_History_win_aspect_ratio = Summary_History_win_aspect_ratio
         s% pg% Summary_History_file_width = Summary_History_file_width
         s% pg% Summary_History_file_aspect_ratio = Summary_History_file_aspect_ratio
         s% pg% Summary_History_xleft = Summary_History_xleft
         s% pg% Summary_History_xright = Summary_History_xright
         s% pg% Summary_History_ybot = Summary_History_ybot
         s% pg% Summary_History_ytop = Summary_History_ytop
         s% pg% Summary_History_txt_scale = Summary_History_txt_scale
         s% pg% Summary_History_title = Summary_History_title
         s% pg% Summary_History_name = Summary_History_name
         s% pg% Summary_History_legend = Summary_History_legend
         s% pg% Summary_History_num_lines = Summary_History_num_lines
         s% pg% Summary_History_use_decorator = Summary_History_use_decorator

         s% pg% Power_win_flag = Power_win_flag
         s% pg% Power_file_flag = Power_file_flag
         s% pg% Power_file_interval = Power_file_interval
         s% pg% Power_file_dir = Power_file_dir
         s% pg% Power_file_prefix = Power_file_prefix
         s% pg% Power_xaxis_name = Power_xaxis_name
         s% pg% Power_xaxis_reversed = Power_xaxis_reversed
         s% pg% Power_legend_max_cnt = Power_legend_max_cnt
         s% pg% Power_legend_txt_scale_factor = Power_legend_txt_scale_factor
         s% pg% Power_xmin = Power_xmin
         s% pg% Power_xmax = Power_xmax
         s% pg% Power_ymin = Power_ymin
         s% pg% Power_ymax = Power_ymax
         s% pg% Power_win_width = Power_win_width
         s% pg% Power_win_aspect_ratio = Power_win_aspect_ratio
         s% pg% Power_file_width = Power_file_width
         s% pg% Power_file_aspect_ratio = Power_file_aspect_ratio
         s% pg% Power_xleft = Power_xleft
         s% pg% Power_xright = Power_xright
         s% pg% Power_ybot = Power_ybot
         s% pg% Power_ytop = Power_ytop
         s% pg% Power_txt_scale = Power_txt_scale
         s% pg% Power_title = Power_title
         s% pg% Power_use_decorator = Power_use_decorator

         s% pg% Abundance_win_flag = Abundance_win_flag
         s% pg% Abundance_file_flag = Abundance_file_flag
         s% pg% Abundance_file_interval = Abundance_file_interval
         s% pg% Abundance_file_dir = Abundance_file_dir
         s% pg% Abundance_file_prefix = Abundance_file_prefix
         s% pg% Abundance_num_isos_to_show = Abundance_num_isos_to_show
         s% pg% num_abundance_line_labels = num_abundance_line_labels
         s% pg% Abundance_legend_txt_scale_factor = Abundance_legend_txt_scale_factor
         s% pg% Abundance_line_txt_scale_factor = Abundance_line_txt_scale_factor
         s% pg% Abundance_line_txt_scale_factor = Abundance_line_txt_scale_factor
         s% pg% Abundance_legend_max_cnt = Abundance_legend_max_cnt
         s% pg% Abundance_which_isos_to_show = Abundance_which_isos_to_show
         s% pg% Abundance_xaxis_name = Abundance_xaxis_name
         s% pg% Abundance_xaxis_reversed = Abundance_xaxis_reversed
         s% pg% Abundance_show_photosphere_location = Abundance_show_photosphere_location
         s% pg% Abundance_xmin = Abundance_xmin
         s% pg% Abundance_xmax = Abundance_xmax
         s% pg% Abundance_log_mass_frac_min = Abundance_log_mass_frac_min
         s% pg% Abundance_log_mass_frac_max = Abundance_log_mass_frac_max
         s% pg% Abundance_win_width = Abundance_win_width
         s% pg% Abundance_win_aspect_ratio = Abundance_win_aspect_ratio
         s% pg% Abundance_file_width = Abundance_file_width
         s% pg% Abundance_file_aspect_ratio = Abundance_file_aspect_ratio
         s% pg% Abundance_xleft = Abundance_xleft
         s% pg% Abundance_xright = Abundance_xright
         s% pg% Abundance_ybot = Abundance_ybot
         s% pg% Abundance_ytop = Abundance_ytop
         s% pg% Abundance_txt_scale = Abundance_txt_scale
         s% pg% Abundance_title = Abundance_title
         s% pg% Abundance_use_decorator = Abundance_use_decorator
         
         s% pg% dPg_dnu_win_flag = dPg_dnu_win_flag
         s% pg% dPg_dnu_file_flag = dPg_dnu_file_flag
         s% pg% dPg_dnu_xleft = dPg_dnu_xleft
         s% pg% dPg_dnu_xright = dPg_dnu_xright
         s% pg% dPg_dnu_ybot = dPg_dnu_ybot
         s% pg% dPg_dnu_ytop = dPg_dnu_ytop
         s% pg% dPg_dnu_txt_scale = dPg_dnu_txt_scale
         s% pg% dPg_dnu_title = dPg_dnu_title
         s% pg% dPg_dnu_file_interval = dPg_dnu_file_interval
         s% pg% dPg_dnu_step_min = dPg_dnu_step_min
         s% pg% dPg_dnu_step_max = dPg_dnu_step_max
         s% pg% dPg_dnu_file_dir = dPg_dnu_file_dir
         s% pg% dPg_dnu_file_prefix = dPg_dnu_file_prefix
         s% pg% show_dPg_dnu_annotation1 = show_dPg_dnu_annotation1
         s% pg% show_dPg_dnu_annotation2 = show_dPg_dnu_annotation2
         s% pg% show_dPg_dnu_annotation3 = show_dPg_dnu_annotation3
         s% pg% dPg_dnu_fname = dPg_dnu_fname
         s% pg% dPg_dnu_delta_nu_min = dPg_dnu_delta_nu_min
         s% pg% dPg_dnu_delta_nu_max = dPg_dnu_delta_nu_max
         s% pg% dPg_dnu_delta_Pg_min = dPg_dnu_delta_Pg_min
         s% pg% dPg_dnu_delta_Pg_max = dPg_dnu_delta_Pg_max
         s% pg% dPg_dnu_delta_nu_margin = dPg_dnu_delta_nu_margin
         s% pg% dPg_dnu_delta_Pg_margin = dPg_dnu_delta_Pg_margin
         s% pg% dPg_dnu_d_delta_nu_min = dPg_dnu_d_delta_nu_min
         s% pg% dPg_dnu_d_delta_Pg_min = dPg_dnu_d_delta_Pg_min
         s% pg% dPg_dnu_win_width = dPg_dnu_win_width
         s% pg% dPg_dnu_win_aspect_ratio = dPg_dnu_win_aspect_ratio
         s% pg% dPg_dnu_file_width = dPg_dnu_file_width
         s% pg% dPg_dnu_file_aspect_ratio = dPg_dnu_file_aspect_ratio
         s% pg% dPg_dnu_use_decorator = dPg_dnu_use_decorator

         s% pg% Grid_win_flag = Grid_win_flag
         s% pg% Grid_win_width = Grid_win_width
         s% pg% Grid_win_aspect_ratio = Grid_win_aspect_ratio
         s% pg% Grid_xleft = Grid_xleft
         s% pg% Grid_xright = Grid_xright
         s% pg% Grid_ybot = Grid_ybot
         s% pg% Grid_ytop = Grid_ytop
         s% pg% Grid_title = Grid_title
         s% pg% Grid_txt_scale_factor = Grid_txt_scale_factor
         s% pg% Grid_num_cols = Grid_num_cols
         s% pg% Grid_num_rows = Grid_num_rows
         s% pg% Grid_num_plots = Grid_num_plots
         s% pg% Grid_plot_name = Grid_plot_name
         s% pg% Grid_plot_row = Grid_plot_row
         s% pg% Grid_plot_rowspan = Grid_plot_rowspan
         s% pg% Grid_plot_col = Grid_plot_col
         s% pg% Grid_plot_colspan = Grid_plot_colspan
         s% pg% Grid_plot_pad_left = Grid_plot_pad_left
         s% pg% Grid_plot_pad_right = Grid_plot_pad_right
         s% pg% Grid_plot_pad_top = Grid_plot_pad_top
         s% pg% Grid_plot_pad_bot = Grid_plot_pad_bot
         s% pg% Grid_file_flag = Grid_file_flag
         s% pg% Grid_file_dir = Grid_file_dir
         s% pg% Grid_file_prefix = Grid_file_prefix
         s% pg% Grid_file_interval = Grid_file_interval
         s% pg% Grid_file_width = Grid_file_width
         s% pg% Grid_file_aspect_ratio = Grid_file_aspect_ratio

         s% pg% annotation1_ci = annotation1_ci
         s% pg% annotation1_ch = annotation1_ch
         s% pg% annotation1_lw = annotation1_lw
         s% pg% annotation1_cf = annotation1_cf
         s% pg% annotation1_text = annotation1_text
         s% pg% annotation1_side = annotation1_side
         s% pg% annotation1_disp = annotation1_disp
         s% pg% annotation1_coord = annotation1_coord
         s% pg% annotation1_fjust = annotation1_fjust

         s% pg% annotation2_ci = annotation2_ci
         s% pg% annotation2_ch = annotation2_ch
         s% pg% annotation2_lw = annotation2_lw
         s% pg% annotation2_cf = annotation2_cf
         s% pg% annotation2_text = annotation2_text
         s% pg% annotation2_side = annotation2_side
         s% pg% annotation2_disp = annotation2_disp
         s% pg% annotation2_coord = annotation2_coord
         s% pg% annotation2_fjust = annotation2_fjust

         s% pg% annotation3_ci = annotation3_ci
         s% pg% annotation3_ch = annotation3_ch
         s% pg% annotation3_lw = annotation3_lw
         s% pg% annotation3_cf = annotation3_cf
         s% pg% annotation3_text = annotation3_text
         s% pg% annotation3_side = annotation3_side
         s% pg% annotation3_disp = annotation3_disp
         s% pg% annotation3_coord = annotation3_coord
         s% pg% annotation3_fjust = annotation3_fjust

         s% pg% read_extra_pgstar_inlist = read_extra_pgstar_inlist
         s% pg% extra_pgstar_inlist_name = extra_pgstar_inlist_name

      end subroutine store_pgstar_controls


      subroutine set_default_pgstar_controls

         include 'pgstar.defaults'

      end subroutine set_default_pgstar_controls

      end module pgstar_ctrls_io

