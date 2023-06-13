      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! PGStar state

      type (pgstar_win_file_data) :: pgstar_win_file_ptr(num_pgstar_plots)

      type (pgstar_hist_node), pointer :: pgstar_hist => null()

      ! PGstar controls
      integer :: pgstar_interval
      character (len=32) :: file_device
      character (len=32) :: file_extension
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

      logical, dimension(pgstar_array_length) :: Text_Summary_win_flag, Text_Summary_file_flag
      integer, dimension(pgstar_array_length) :: id_Text_Summary_win, id_Text_Summary_file, Text_Summary_file_interval
      character (len=strlen), dimension(pgstar_array_length) :: Text_Summary_title, Text_Summary_file_dir, &
         Text_Summary_file_prefix
      integer, dimension(pgstar_array_length) :: Text_Summary_num_cols, Text_Summary_num_rows
      character (len=64), dimension(pgstar_array_length, max_num_rows_Text_Summary, max_num_cols_Text_Summary) :: &
         Text_Summary_name
      real, dimension(pgstar_array_length) :: &
         Text_Summary_xleft, Text_Summary_xright, &
         Text_Summary_ybot, Text_Summary_ytop, Text_Summary_txt_scale, Text_Summary_dxval, &
         Text_Summary_win_width, Text_Summary_win_aspect_ratio, &
         prev_Text_Summary_win_width, prev_Text_Summary_win_ratio, &
         Text_Summary_file_width, Text_Summary_file_aspect_ratio, &
         prev_Text_Summary_file_width, prev_Text_Summary_file_ratio

      logical :: dPg_dnu_win_flag, dPg_dnu_file_flag
      logical :: show_dPg_dnu_annotation1, &
         show_dPg_dnu_annotation2, &
         show_dPg_dnu_annotation3
      logical :: show_dPg_dnu_target_box
      integer :: id_dPg_dnu_win, id_dPg_dnu_file, dPg_dnu_file_interval
      integer :: dPg_dnu_step_min, dPg_dnu_step_max
      character (len=strlen) :: &
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

      logical, dimension(pgstar_array_length) :: History_Track_win_flag, History_Track_file_flag
      logical, dimension(pgstar_array_length):: show_History_Track_annotation1, &
         show_History_Track_annotation2, show_History_Track_annotation3
      logical, dimension(pgstar_array_length) :: show_History_Track_target_box, &
         History_Track_reverse_xaxis, History_Track_reverse_yaxis, &
         History_Track_log_xaxis, History_Track_log_yaxis
      integer, dimension(pgstar_array_length) :: id_History_Track_win, id_History_Track_file, History_Track_file_interval
      integer, dimension(pgstar_array_length) :: History_Track_step_min, History_Track_step_max, History_Track_n_sigma
      character (len=strlen), dimension(pgstar_array_length) :: &
         History_Track_file_dir, History_Track_file_prefix, History_Track_title, &
         History_Track_xname, History_Track_xaxis_label, &
         History_Track_yname, History_Track_yaxis_label, History_Track_fname
      real, dimension(pgstar_array_length) :: &
         History_Track_xmin, History_Track_xmax, &
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

      logical :: HR_win_flag, HR_file_flag
      logical :: show_HR_annotation1, show_HR_annotation2, show_HR_annotation3
      logical :: show_HR_target_box, show_HR_classical_instability_strip, &
         show_HR_Mira_instability_region, show_HR_WD_instabilities
      integer :: id_HR_win, id_HR_file, HR_file_interval
      integer :: HR_step_min, HR_step_max, HR_target_n_sigma
      character (len=strlen) :: HR_file_dir, HR_file_prefix, HR_title, HR_fname
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
      logical :: show_TRho_annotation1, show_TRho_annotation2, show_TRho_annotation3, &
         show_TRho_degeneracy_line
      integer :: TRho_file_interval, TRho_step_max
      integer :: id_TRho_win, id_TRho_file, TRho_step_min
      character (len=strlen) :: TRho_file_dir, TRho_file_prefix, TRho_title, TRho_fname
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
      logical :: show_TmaxRho_annotation1, show_TmaxRho_annotation2, show_TmaxRho_annotation3, &
         show_TmaxRho_degeneracy_line
      integer :: TmaxRho_file_interval, TmaxRho_step_max
      integer :: id_TmaxRho_win, id_TmaxRho_file, TmaxRho_step_min
      character (len=strlen) :: TmaxRho_file_dir, TmaxRho_file_prefix, TmaxRho_title, TmaxRho_fname
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
      logical :: show_logL_R_target_box, show_logL_photosphere_r
      logical :: show_logL_R_annotation1, show_logL_R_annotation2, show_logL_R_annotation3
      integer :: logL_R_file_interval, logL_R_step_max, logL_R_target_n_sigma
      integer :: id_logL_R_win, id_logL_R_file, logL_R_step_min
      character (len=strlen) :: &
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
      logical :: show_logL_Teff_target_box
      logical :: show_logL_Teff_annotation1, show_logL_Teff_annotation2, show_logL_Teff_annotation3
      integer :: logL_Teff_file_interval, logL_Teff_step_max, logL_Teff_target_n_sigma
      integer :: id_logL_Teff_win, id_logL_Teff_file, logL_Teff_step_min
      character (len=strlen) :: &
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
      logical :: show_logL_v_target_box, show_logL_photosphere_v
      logical :: show_logL_v_annotation1, show_logL_v_annotation2, show_logL_v_annotation3
      integer :: logL_v_file_interval, logL_v_step_max, logL_v_target_n_sigma
      integer :: id_logL_v_win, id_logL_v_file, logL_v_step_min
      character (len=strlen) :: &
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
      logical :: show_L_Teff_target_box
      logical :: show_L_Teff_annotation1, show_L_Teff_annotation2, show_L_Teff_annotation3
      integer :: L_Teff_file_interval, L_Teff_step_max, L_Teff_target_n_sigma
      integer :: id_L_Teff_win, id_L_Teff_file, L_Teff_step_min
      character (len=strlen) :: &
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
      logical :: show_L_v_target_box
      logical :: show_L_v_annotation1, show_L_v_annotation2, show_L_v_annotation3
      integer :: L_v_file_interval, L_v_step_max, L_v_target_n_sigma
      integer :: id_L_v_win, id_L_v_file, L_v_step_min
      character (len=strlen) :: &
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
      logical :: show_L_R_target_box
      logical :: show_L_R_annotation1, show_L_R_annotation2, show_L_R_annotation3
      integer :: L_R_file_interval, L_R_step_max, L_R_target_n_sigma
      integer :: id_L_R_win, id_L_R_file, L_R_step_min
      character (len=strlen) :: &
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
      logical :: show_R_Teff_target_box
      logical :: show_R_Teff_annotation1, show_R_Teff_annotation2, show_R_Teff_annotation3
      integer :: R_Teff_file_interval, R_Teff_step_max, R_Teff_target_n_sigma
      integer :: id_R_Teff_win, id_R_Teff_file, R_Teff_step_min
      character (len=strlen) :: &
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
      logical :: show_R_L_target_box
      logical :: show_R_L_annotation1, show_R_L_annotation2, show_R_L_annotation3
      integer :: R_L_file_interval, R_L_step_max, R_L_target_n_sigma
      integer :: id_R_L_win, id_R_L_file, R_L_step_min
      character (len=strlen) :: &
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
      logical :: show_logg_Teff_target_box
      logical :: show_logg_Teff_annotation1, show_logg_Teff_annotation2, show_logg_Teff_annotation3
      integer :: logg_Teff_file_interval, logg_Teff_step_max, logg_Teff_target_n_sigma
      integer :: id_logg_Teff_win, id_logg_Teff_file, logg_Teff_step_min
      character (len=strlen) :: &
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
      character (len=strlen) :: &
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
      integer :: id_rti_win, id_rti_file, rti_file_interval
      integer :: rti_step_xmin,rti_step_xmax
      integer :: rti_max_width, rti_interval, rti_line_weight
      character (len=strlen) :: rti_file_dir, rti_file_prefix, rti_title, rti_xaxis_name
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
      integer :: id_Kipp_win, id_Kipp_file, Kipp_file_interval
      integer :: Kipp_step_xmin,Kipp_step_xmax
      integer :: Kipp_max_width, Kipp_mix_interval, &
         Kipp_burn_line_weight, Kipp_mix_line_weight, &
         Kipp_luminosities_line_weight, Kipp_masses_line_weight
      character (len=strlen) :: Kipp_file_dir, Kipp_file_prefix, Kipp_title, Kipp_xaxis_name
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
         Kipp_xmax, Kipp_xmin,Kipp_xmargin, &
         Kipp_burn_type_cutoff

      logical :: TRho_Profile_win_flag, TRho_Profile_file_flag, TRho_switch_to_Column_Depth
      logical :: TRho_switch_to_mass
      integer :: id_TRho_Profile_win, id_TRho_Profile_file, TRho_Profile_file_interval
      character (len=strlen) :: &
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
      character (len=50) :: profile_mass_point_str(max_num_profile_mass_points)
      integer :: profile_mass_point_str_clr(max_num_profile_mass_points)
      real :: profile_mass_point_str_scale(max_num_profile_mass_points)

      logical :: Profile_Panels_show_Mach_1_location, &
         Profile_Panels_show_photosphere_location
      real :: Profile_Panels_xwidth_left_div_shock_value, &
         Profile_Panels_xwidth_right_div_shock_value, &
         Profile_Panels_xwidth_left_of_shock, &
         Profile_Panels_xwidth_right_of_shock

      logical, dimension(pgstar_array_length) :: Profile_Panels_win_flag, Profile_Panels_file_flag
      integer, dimension(pgstar_array_length) :: id_Profile_Panels_win, id_Profile_Panels_file, Profile_Panels_file_interval
      character (len=strlen), dimension(pgstar_array_length) :: Profile_Panels_file_dir, Profile_Panels_file_prefix
      logical, dimension(pgstar_array_length) :: Profile_Panels_xaxis_reversed, Profile_Panels_show_mix_regions_on_xaxis
      character (len=strlen), dimension(pgstar_array_length) :: Profile_Panels_xaxis_name, Profile_Panels_title
      real, dimension(pgstar_array_length) :: &
         Profile_Panels_xmin, Profile_Panels_xmax, Profile_Panels_xmargin, &
         Profile_Panels_win_width, Profile_Panels_win_aspect_ratio, &
         Profile_Panels_xleft, Profile_Panels_xright, &
         Profile_Panels_ybot, Profile_Panels_ytop, Profile_Panels_txt_scale, &
         prev_Profile_Panels_win_width, prev_Profile_Panels_win_ratio, &
         Profile_Panels_file_width, Profile_Panels_file_aspect_ratio, &
         prev_Profile_Panels_file_width, prev_Profile_Panels_file_ratio
      integer, dimension(pgstar_array_length) :: Profile_Panels_num_panels
      character (len=strlen), dimension(pgstar_array_length, max_num_Profile_Panels) :: &
         Profile_Panels_yaxis_name, Profile_Panels_other_yaxis_name
      logical, dimension(pgstar_array_length, max_num_Profile_Panels) :: &
         Profile_Panels_yaxis_reversed, Profile_Panels_other_yaxis_reversed, &
         Profile_Panels_yaxis_log, Profile_Panels_other_yaxis_log, &
         Profile_Panels_same_yaxis_range
      real, dimension(pgstar_array_length, max_num_Profile_Panels) :: &
         Profile_Panels_ymin, Profile_Panels_other_ymin, &
         Profile_Panels_ymax, Profile_Panels_other_ymax, &
         Profile_Panels_ycenter, Profile_Panels_other_ycenter, &
         Profile_Panels_ymargin, Profile_Panels_other_ymargin, &
         Profile_Panels_dymin, Profile_Panels_other_dymin
      logical, dimension(pgstar_array_length) :: Profile_Panels_show_grid

      logical :: Dynamo_win_flag, Dynamo_file_flag
      integer :: id_Dynamo_win, id_Dynamo_file, Dynamo_file_interval
      character (len=strlen) :: Dynamo_file_dir, Dynamo_file_prefix
      logical :: show_Dynamo_annotation1, show_Dynamo_annotation2, show_Dynamo_annotation3, &
         Dynamo_xaxis_reversed
      character (len=strlen) :: Dynamo_xaxis_name, Dynamo_title
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
      logical :: Mixing_xaxis_reversed, Mixing_show_rotation_details
      integer :: id_Mixing_win, id_Mixing_file, Mixing_file_interval
      character (len=strlen) :: Mixing_file_dir, Mixing_file_prefix
      logical :: show_Mixing_annotation1, show_Mixing_annotation2, show_Mixing_annotation3
      character (len=strlen) :: Mixing_xaxis_name, Mixing_title
      real :: &
         Mixing_xmin, Mixing_xmax, Mixing_dymin, &
         Mixing_ymin, Mixing_ymax, Mixing_legend_txt_scale_factor, &
         Mixing_win_width, Mixing_win_aspect_ratio, &
         Mixing_xleft, Mixing_xright, &
         Mixing_ybot, Mixing_ytop, Mixing_txt_scale, &
         prev_Mixing_win_width, prev_Mixing_win_ratio, &
         Mixing_file_width, Mixing_file_aspect_ratio, &
         prev_Mixing_file_width, prev_Mixing_file_ratio

      logical, dimension(pgstar_array_length) :: History_Panels_win_flag, History_Panels_file_flag
      integer, dimension(pgstar_array_length) :: id_History_Panels_win, id_History_Panels_file, History_Panels_file_interval
      character (len=strlen), dimension(pgstar_array_length) :: History_Panels_file_dir, History_Panels_file_prefix
      logical, dimension(pgstar_array_length) :: History_Panels_xaxis_reversed, History_Panels_automatic_star_age_units, History_Panels_xaxis_log
      character (len=strlen), dimension(pgstar_array_length) :: History_Panels_xaxis_name, History_Panels_title
      real, dimension(pgstar_array_length) :: &
         History_Panels_xmin, History_Panels_xmax, History_Panels_dxmin,&
         History_Panels_max_width,History_Panels_xmargin, &
         History_Panels_win_width, History_Panels_win_aspect_ratio, &
         History_Panels_xleft, History_Panels_xright, &
         History_Panels_ybot, History_Panels_ytop, History_Panels_txt_scale, &
         prev_History_Panels_win_width, prev_History_Panels_win_ratio, &
         History_Panels_file_width, History_Panels_file_aspect_ratio, &
         prev_History_Panels_file_width, prev_History_Panels_file_ratio
      integer, dimension(pgstar_array_length) :: History_Panels_num_panels
      character (len=strlen), dimension(pgstar_array_length, max_num_History_Panels) :: &
         History_Panels_yaxis_name, History_Panels_other_yaxis_name, &
         History_Panels_points_name
      logical, dimension(pgstar_array_length, max_num_History_Panels) :: &
         History_Panels_yaxis_reversed, History_Panels_other_yaxis_reversed, &
         History_Panels_yaxis_log, History_Panels_other_yaxis_log, &
         History_Panels_same_yaxis_range
      real, dimension(pgstar_array_length, max_num_History_Panels) :: &
         History_Panels_ymin, History_Panels_other_ymin, &
         History_Panels_ymax, History_Panels_other_ymax, &
         History_Panels_ymargin, History_Panels_other_ymargin, &
         History_Panels_dymin, History_Panels_other_dymin
      logical :: History_Panel_points_error_bars
      integer :: History_Panel_points_interval, History_Panel_points_marker, &
         History_Panel_points_ci, History_Panel_points_lw
      real :: History_Panel_points_ch

      logical, dimension(pgstar_array_length) :: Color_Magnitude_win_flag, Color_Magnitude_file_flag
      integer, dimension(pgstar_array_length) :: id_Color_Magnitude_win, id_Color_Magnitude_file, Color_Magnitude_file_interval
      character (len=strlen), dimension(pgstar_array_length) :: Color_Magnitude_file_dir, Color_Magnitude_file_prefix
      logical, dimension(pgstar_array_length) :: Color_Magnitude_xaxis_reversed, Color_Magnitude_xaxis_log
      character (len=strlen), dimension(pgstar_array_length) :: Color_Magnitude_xaxis1_name, Color_Magnitude_xaxis2_name, &
         Color_Magnitude_title
      real, dimension(pgstar_array_length) :: &
         Color_Magnitude_xmin, Color_Magnitude_xmax, Color_Magnitude_dxmin, Color_Magnitude_xmargin, &
         Color_Magnitude_max_width, Color_Magnitude_win_width, Color_Magnitude_win_aspect_ratio, &
         Color_Magnitude_xleft, Color_Magnitude_xright, &
         Color_Magnitude_ybot, Color_Magnitude_ytop, Color_Magnitude_txt_scale, &
         prev_Color_Magnitude_win_width, prev_Color_Magnitude_win_ratio, &
         Color_Magnitude_file_width, Color_Magnitude_file_aspect_ratio, &
         prev_Color_Magnitude_file_width, prev_Color_Magnitude_file_ratio
      integer, dimension(pgstar_array_length) :: Color_Magnitude_num_panels
      character (len=strlen), dimension(pgstar_array_length, max_num_Color_Magnitude) :: &
         Color_Magnitude_yaxis1_name, Color_Magnitude_other_yaxis1_name, &
         Color_Magnitude_yaxis2_name, Color_Magnitude_other_yaxis2_name
      logical, dimension(pgstar_array_length, max_num_Color_Magnitude) :: &
         Color_Magnitude_yaxis_reversed, Color_Magnitude_other_yaxis_reversed, &
         Color_Magnitude_yaxis_log, Color_Magnitude_other_yaxis_log
      real, dimension(pgstar_array_length, max_num_Color_Magnitude) :: &
         Color_Magnitude_ymin, Color_Magnitude_other_ymin, &
         Color_Magnitude_ymax, Color_Magnitude_other_ymax, &
         Color_Magnitude_ymargin, Color_Magnitude_other_ymargin, &
         Color_Magnitude_dymin, Color_Magnitude_other_dymin

      logical :: Mode_Prop_win_flag, Mode_Prop_file_flag
      logical :: Mode_Prop_xaxis_reversed
      integer :: id_Mode_Prop_win, id_Mode_Prop_file, Mode_Prop_file_interval
      character (len=strlen) :: Mode_Prop_file_dir, Mode_Prop_file_prefix, &
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
      logical :: Summary_Burn_xaxis_reversed
      integer :: id_Summary_Burn_win, id_Summary_Burn_file, &
         Summary_Burn_num_lines, Summary_Burn_file_interval
      character (len=strlen) :: Summary_Burn_file_dir, Summary_Burn_file_prefix, &
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
      logical :: Summary_Profile_xaxis_reversed, &
         Summary_Profile_scaled_value(max_num_Summary_Profile_Lines)
      integer :: id_Summary_Profile_win, id_Summary_Profile_file, &
         Summary_Profile_num_lines, Summary_Profile_file_interval
      character (len=strlen) :: Summary_Profile_file_dir, Summary_Profile_file_prefix, &
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
      logical :: Summary_History_scaled_value(max_num_Summary_History_Lines)
      integer :: id_Summary_History_win, id_Summary_History_file, &
         Summary_History_num_lines, Summary_History_file_interval, &
         Summary_History_xmin, Summary_History_xmax, Summary_History_max_width
      character (len=strlen) :: Summary_History_file_dir, Summary_History_file_prefix, &
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
      logical :: Power_xaxis_reversed
      integer :: id_Power_win, id_Power_file, Power_file_interval, &
         Power_legend_max_cnt
      character (len=strlen) :: Power_file_dir, Power_file_prefix, &
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
      logical :: Abundance_xaxis_reversed, Abundance_show_photosphere_location
      integer :: id_Abundance_win, id_Abundance_file, Abundance_file_interval
      integer :: Abundance_num_isos_to_show, num_abundance_line_labels, &
         Abundance_legend_max_cnt
      character(len=iso_name_length) :: &
         Abundance_which_isos_to_show(max_Abundance_num_isos_to_show)
      character (len=strlen) :: Abundance_file_dir, Abundance_file_prefix, &
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
      logical :: Network_show_colorbar, Network_show_element_names, &
            Network_show_mass_fraction
      integer :: id_Network_win, id_Network_file, Network_file_interval
      character (len=strlen) :: Network_file_dir, Network_file_prefix, &
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
      logical :: Production_show_element_names
      integer :: id_Production_win, id_Production_file, Production_file_interval
      character (len=strlen) :: Production_file_dir, Production_file_prefix, &
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

      logical, dimension(pgstar_array_length) :: Grid_win_flag, Grid_file_flag
      integer, dimension(pgstar_array_length) :: Grid_file_interval, Grid_step_max
      integer, dimension(pgstar_array_length) :: id_Grid_win, id_Grid_file
      character (len=strlen), dimension(pgstar_array_length) :: Grid_file_dir, Grid_file_prefix, Grid_title
      real, dimension(pgstar_array_length) :: &
         Grid_win_width, Grid_win_aspect_ratio, &
         Grid_xleft, Grid_xright, &
         Grid_ybot, Grid_ytop, &
         prev_Grid_win_width, prev_Grid_win_ratio, &
         Grid_file_width, Grid_file_aspect_ratio, &
         prev_Grid_file_width, prev_Grid_file_ratio
      character (len=strlen), dimension(pgstar_array_length, max_num_pgstar_grid_plots) :: &
         Grid_plot_name
      real, dimension(pgstar_array_length, max_num_pgstar_grid_plots) :: &
         Grid_txt_scale_factor, &
         Grid_plot_pad_left, Grid_plot_pad_right, &
         Grid_plot_pad_top, Grid_plot_pad_bot
      integer, dimension(pgstar_array_length) :: &
         Grid_num_cols, Grid_num_rows, Grid_num_plots
      integer, dimension(pgstar_array_length, max_num_pgstar_grid_plots) :: &
         Grid_plot_row, Grid_plot_rowspan, &
         Grid_plot_col, Grid_plot_colspan

      integer :: annotation1_ci, annotation1_lw, annotation1_cf
      real :: annotation1_ch
      character (len=strlen) :: annotation1_text
      character :: annotation1_side
      real :: annotation1_disp, annotation1_coord, annotation1_fjust

      integer :: annotation2_ci, annotation2_lw, annotation2_cf
      real :: annotation2_ch
      character (len=strlen) :: annotation2_text
      character :: annotation2_side
      real :: annotation2_disp, annotation2_coord, annotation2_fjust

      integer :: annotation3_ci, annotation3_lw, annotation3_cf
      real :: annotation3_ch
      character (len=356) :: annotation3_text
      character :: annotation3_side
      real :: annotation3_disp, annotation3_coord, annotation3_fjust

      logical, dimension(max_extra_inlists) :: read_extra_pgstar_inlist
      character (len=strlen), dimension(max_extra_inlists) :: extra_pgstar_inlist_name

      logical, dimension(pgstar_array_length) :: &
            Color_Magnitude_use_decorator, &
            History_Panels_use_decorator, &
            History_Track_use_decorator, &
            Profile_Panels_use_decorator

      logical :: Abundance_use_decorator, &
            dPg_dnu_use_decorator, &
            Dynamo_use_decorator, &
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
            R_L_use_decorator, &
            R_Teff_use_decorator, &
            rti_use_decorator, &
            summary_burn_use_decorator, &
            summary_profile_use_decorator, &
            summary_history_use_decorator, &
            trho_use_decorator, &
            tmaxrho_use_decorator, &
            trho_profile_use_decorator
