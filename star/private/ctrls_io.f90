! ***********************************************************************
!
! Copyright (C) 2010 The Mesa Team
!
! MESA is free software; you can use it and/or modify
! it under the combined terms and restrictions of the MESA MANIFESTO
! and the GNU General Library Public License as published
! by the Free Software Foundation; either version 2 of the License,
! or (at your option) any later version.
!
! You should have received a copy of the MESA MANIFESTO along with
! this software; if not, it is available at the mesa website:
! http://mesa.sourceforge.net/
!
! MESA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
!
! You should have received a copy of the GNU Library General Public License
! along with this software; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

 module ctrls_io

 use const_def
 use star_private_def

 implicit none

 include 'star_controls.inc'
 include 'star_controls_dev.inc'

 logical :: read_extra_controls_inlist1
 character (len=strlen) :: extra_controls_inlist1_name

 logical :: read_extra_controls_inlist2
 character (len=strlen) :: extra_controls_inlist2_name

 logical :: read_extra_controls_inlist3
 character (len=strlen) :: extra_controls_inlist3_name

 logical :: read_extra_controls_inlist4
 character (len=strlen) :: extra_controls_inlist4_name

 logical :: read_extra_controls_inlist5
 character (len=strlen) :: extra_controls_inlist5_name

 logical :: save_controls_namelist
 character (len=strlen) :: controls_namelist_name

 namelist /controls/ &
 
    ! where to start
    initial_mass, initial_z, initial_y, initial_he3, &
    
    ! definition of core boundaries
    he_core_boundary_h1_fraction, co_core_boundary_he4_fraction, one_core_boundary_he4_c12_fraction, &
    fe_core_boundary_si28_fraction, neutron_rich_core_boundary_Ye_max, min_boundary_fraction, &
    
    ! when to stop
    max_model_number, relax_max_number_retries, max_number_retries, max_age, max_age_in_seconds, max_age_in_days, &
    num_adjusted_dt_steps_before_max_age, dt_years_for_steps_before_max_age, max_abs_rel_run_E_err, &
    reduction_factor_for_max_timestep, when_to_stop_rtol, when_to_stop_atol, &
    gamma_center_limit, eta_center_limit, log_center_temp_limit, log_max_temp_upper_limit, &
    log_max_temp_lower_limit , log_center_temp_lower_limit, log_center_density_limit, &
    log_center_density_lower_limit, center_entropy_limit, center_entropy_lower_limit, &
    max_entropy_limit, max_entropy_lower_limit, min_timestep_limit, non_fe_core_rebound_limit, &
    fe_core_infall_limit, center_Ye_lower_limit, center_R_lower_limit, non_fe_core_infall_limit, &
    v_div_csound_surf_limit, v_div_csound_max_limit, Lnuc_div_L_upper_limit, Lnuc_div_L_lower_limit,&
    v_surf_div_v_kh_upper_limit, v_surf_div_v_kh_lower_limit, v_surf_div_v_esc_limit, v_surf_kms_limit, &
    stop_near_zams, Lnuc_div_L_zams_limit, Pgas_div_P_limit, Pgas_div_P_limit_max_q, gamma1_limit, gamma1_limit_max_q, &
    stop_at_phase_PreMS, stop_at_phase_ZAMS, stop_at_phase_IAMS, stop_at_phase_TAMS, gamma1_limit_max_v_div_vesc, &
    stop_at_phase_He_Burn, stop_at_phase_ZACHeB, stop_at_phase_TACHeB, &
    stop_at_phase_TP_AGB, stop_at_phase_C_Burn, stop_at_phase_Ne_Burn, &
    stop_at_phase_O_Burn, stop_at_phase_Si_Burn, stop_at_phase_WDCS, &
    peak_burn_vconv_div_cs_limit, omega_div_omega_crit_limit, delta_nu_lower_limit, &
    delta_nu_upper_limit, delta_Pg_lower_limit, delta_Pg_upper_limit, shock_mass_upper_limit, &
    mach1_mass_upper_limit, stop_when_reach_this_cumulative_extra_heating, &
    xa_central_lower_limit_species, xa_central_lower_limit, xa_central_upper_limit_species, xa_central_upper_limit, &
    xa_surface_lower_limit_species, xa_surface_lower_limit, xa_surface_upper_limit_species, xa_surface_upper_limit, &
    xa_average_lower_limit_species, xa_average_lower_limit, xa_average_upper_limit_species, xa_average_upper_limit, &
    star_species_mass_min_limit, star_species_mass_min_limit_iso, star_species_mass_max_limit, star_species_mass_max_limit_iso, &
    xmstar_min_limit, xmstar_max_limit, envelope_mass_limit, envelope_fraction_left_limit, &
    he_core_mass_limit, co_core_mass_limit, one_core_mass_limit, &
    fe_core_mass_limit, neutron_rich_core_mass_limit, HB_limit, star_mass_min_limit, star_mass_max_limit, &
    he_layer_mass_lower_limit, abs_diff_lg_LH_lg_Ls_limit, Teff_upper_limit, Teff_lower_limit, &
    photosphere_m_upper_limit, photosphere_m_lower_limit, photosphere_m_sub_M_center_limit, &
    photosphere_r_upper_limit, photosphere_r_lower_limit, log_Teff_upper_limit, log_Teff_lower_limit, &
    log_Tsurf_upper_limit, log_Tsurf_lower_limit, log_Rsurf_upper_limit, log_Rsurf_lower_limit, &
    log_Psurf_upper_limit, log_Psurf_lower_limit, remnant_mass_min_limit, ejecta_mass_max_limit, &
    log_Dsurf_upper_limit, log_Dsurf_lower_limit, log_L_upper_limit, log_L_lower_limit, &
    log_g_upper_limit, log_g_lower_limit, power_nuc_burn_upper_limit, power_h_burn_upper_limit, &
    power_he_burn_upper_limit, power_z_burn_upper_limit, power_nuc_burn_lower_limit, &
    power_h_burn_lower_limit, power_he_burn_lower_limit, power_z_burn_lower_limit, &
    
    ! max timesteps
    max_timestep, max_years_for_timestep, &
    hi_T_max_years_for_timestep, max_timestep_hi_T_limit, &
    
    ! output of "snapshots" for restarts
    photo_interval, photo_digits, photo_directory, &
    
    ! output of logs and profiles
    do_history_file, history_interval, write_header_frequency, terminal_interval, &
    terminal_show_age_units, terminal_show_timestep_units, terminal_show_log_dt, terminal_show_log_age, &
    num_trace_history_values, trace_history_value_name, write_profiles_flag, profile_interval, &
    priority_profile_interval, log_directory, star_history_name, star_history_header_name, &
    star_history_dbl_format, star_history_int_format, star_history_txt_format, extra_terminal_output_file, &
    profiles_index_name, profile_data_prefix, profile_data_suffix, profile_data_header_suffix, &
    profile_int_format, profile_txt_format, profile_dbl_format, profile_header_include_sys_details, &
    profile_model, max_num_profile_models, max_num_profile_zones, &
    write_controls_info_with_profile, controls_data_prefix, controls_data_suffix, &
    write_pulse_data_with_profile, pulse_data_format, add_atmosphere_to_pulse_data, &
    add_center_point_to_pulse_data, keep_surface_point_for_pulse_data, add_double_points_to_pulse_data, &
    interpolate_rho_for_pulse_data, threshold_grad_mu_for_double_point, max_number_of_double_points,&
    fgong_header, fgong_ivers, &
    max_num_gyre_points, format_for_OSC_data, &
    fgong_zero_A_inside_r, use_other_export_pulse_data, use_other_get_pulse_data, use_other_edit_pulse_data, &
    write_model_with_profile, model_data_prefix, model_data_suffix, &
    mixing_D_limit_for_log, trace_mass_location, min_tau_for_max_abs_v_location, &
    min_q_for_inner_mach1_location, max_q_for_outer_mach1_location, &
    mass_depth_for_L_surf, conv_core_gap_dq_limit, &
    alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, &
    
    ! burn zone eps definitions for use in logs and profiles
    burn_min1, burn_min2, &
    max_conv_vel_div_csound_maxq, width_for_limit_conv_vel, max_q_for_limit_conv_vel, &
    max_mass_in_gm_for_limit_conv_vel, max_r_in_cm_for_limit_conv_vel, &
    
    ! for reported surface/center abundances
    surface_avg_abundance_dq, center_avg_value_dq, &
    
    ! mixing parameters 
    min_convective_gap, min_thermohaline_gap, min_semiconvection_gap, min_thermohaline_dropout, &
    max_dropout_gradL_sub_grada, remove_embedded_semiconvection, recalc_mix_info_after_evolve, remove_mixing_glitches, &
    okay_to_remove_mixing_singleton, prune_bad_cz_min_Hp_height, prune_bad_cz_min_log_eps_nuc, &
    redo_conv_for_dr_lt_mixing_length, alpha_semiconvection, okay_to_reduce_gradT_excess, &
    semiconvection_option, use_Ledoux_criterion, D_mix_zero_region_bottom_q, &
    num_cells_for_smooth_gradL_composition_term, threshold_for_smooth_gradL_composition_term, clip_D_limit, &
   gradT_excess_f1, gradT_excess_f2, gradT_excess_age_fraction, gradT_excess_max_change, gradT_excess_lambda1, &
   gradT_excess_beta1, gradT_excess_lambda2, gradT_excess_beta2, gradT_excess_dlambda, gradT_excess_dbeta, gradT_excess_max_center_h1, &
   gradT_excess_min_center_he4, gradT_excess_max_logT, gradT_excess_min_log_tau_full_on, gradT_excess_max_log_tau_full_off, &
    use_superad_reduction, superad_reduction_gamma_limit, superad_reduction_gamma_limit_scale, D_mix_zero_region_top_q, &
    superad_reduction_gamma_inv_scale, superad_reduction_diff_grads_limit, superad_reduction_limit, &
    fix_eps_grav_transition_to_grid, make_gradr_sticky_in_solver_iters, min_logT_for_make_gradr_sticky_in_solver_iters, &
    max_logT_for_mlt, thermohaline_coeff, thermohaline_option, mixing_length_alpha, remove_small_D_limit, &
    alt_scale_height_flag, Henyey_MLT_y_param, Henyey_MLT_nu_param, no_MLT_below_shock, mlt_make_surface_no_mixing, &
    MLT_option, mlt_use_rotation_correction, mlt_Pturb_factor, do_normalize_dqs_as_part_of_set_qs, &
    max_Y_for_burn_z_mix_region, max_X_for_burn_he_mix_region, &
    limit_overshoot_Hp_using_size_of_convection_zone, RSP_min_tau_for_turbulent_flux, &
    predictive_mix, predictive_superad_thresh, predictive_avoid_reversal, predictive_limit_ingestion,&
    predictive_ingestion_factor, predictive_zone_type, predictive_zone_loc, predictive_bdy_loc, &
    predictive_bdy_q_min, predictive_bdy_q_max, T_mix_limit, RSP_report_undercorrections, &
    do_conv_premix, conv_premix_avoid_increase, conv_premix_time_factor, &
    conv_premix_fix_pgas, conv_premix_dump_snapshots, do_premix_heating, &
    overshoot_f, overshoot_f0, overshoot_D0, RSP_Qvisc_linear, dq_D_mix_zero_at_H_He_crossover, &
    overshoot_Delta0, overshoot_mass_full_on, overshoot_mass_full_off, dq_D_mix_zero_at_H_C_crossover, &
    overshoot_scheme, overshoot_zone_type, overshoot_zone_loc, RSP_Qvisc_quadratic, &
    overshoot_bdy_loc, overshoot_D_min, overshoot_brunt_B_max, mlt_gradT_fraction, max_conv_vel_div_csound, &
    max_v_for_convection, max_q_for_convection_with_hydro_on, alpha_RTI_src_max_q, &
    max_v_div_cs_for_convection, max_abs_du_div_cs_for_convection, RSP_max_dt, RSP_relax_dm_tolerance, &
    calculate_Brunt_B, calculate_Brunt_N2, brunt_N2_coefficient, num_cells_for_smooth_brunt_B, &
    threshold_for_smooth_brunt_B, min_magnitude_brunt_B, RSP_max_dt_times_min_rad_diff_time, &
    min_overshoot_q, overshoot_alpha, RSP_target_steps_per_cycle, &
    RSP_max_num_periods, RSP_min_max_R_for_periods, RSP_min_deltaR_for_periods, &
    RSP_min_PERIOD_div_PERIODLIN, RSP_report_limit_dt, RSP_mode_for_setting_PERIODLIN, RSP_initial_dt_factor, &
    RSP_v_div_cs_threshold_for_dt_limit, RSP_max_dt_times_min_dr_div_cs, RSP_thetae, &
    RSP_alfa, RSP_thetaq, RSP_default_PERIODLIN, &
   RSP_theta, RSP_thetat, RSP_thetau, RSP_wtr, RSP_wtc, RSP_wtt, RSP_gam, RSP_max_retries_per_step, RSP_Qvisc_linear_static, &
   RSP_alfa, RSP_alfap, RSP_alfam, RSP_alfat, RSP_alfas, RSP_alfac, RSP_alfad, RSP_gammar, RSP_nz_div_IBOTOM, &
   RSP_efl0, RSP_cq, RSP_zsh, RSP_tol_max_corr, RSP_tol_max_resid, RSP_max_iters_per_try, &
    RTI_smooth_mass, RTI_smooth_iterations, RTI_smooth_fraction, RSP_dq_1_factor, &
    alpha_RTI_diffusion_factor, dudt_RTI_diffusion_factor, dedt_RTI_diffusion_factor, alpha_RTI_src_min_v_div_cs, &
    dlnddt_RTI_diffusion_factor, composition_RTI_diffusion_factor, max_M_RTI_factors_full_on, min_M_RTI_factors_full_off, &
    alpha_RTI_src_min_v_div_cs, alpha_RTI_src_max_q, alpha_RTI_src_min_q, RSP_kick_vsurf_km_per_sec, &
    RSP_nz_outer, RSP_nz, RSP_T_anchor, RSP_T_inner, RSP_Teff, RSP_mass, RSP_L, RSP_X, RSP_Z, RSP_T_inner_tolerance, &
    RSP_relax_max_tries, RSP_hydro_only, &
   RSP_tau_surf_for_atm_grey_with_kap, RSP_fixed_Psurf, RSP_Avel, RSP_Arnd, RSP_relax_adjust_inner_mass_distribution, &
   use_other_RSP_linear_analysis, RSP_use_atm_grey_with_kap_for_Psurf, &
   RSP_fraction_1st_overtone,RSP_fraction_2nd_overtone, RSP_testing, RSP_use_Prad_for_Psurf, RSP_map_zone_interval, &
   RSP_write_map, RSP_map_filename, RSP_map_history_filename, RSP_map_first_period, RSP_map_last_period, &
   use_other_RSP_build_model, RSP_Psurf, RSP_work_period, RSP_work_filename, RSP_nmodes, RSP_surface_tau, &
   set_RSP_Psurf_to_multiple_of_initial_P1, use_RSP_new_start_scheme, RSP_do_check_omega, RSP_report_check_omega_changes, &
    RSP_relax_initial_model, RSP_trace_RSP_build_model, &
   RSP_GREKM_avg_abs_limit, RSP_GREKM_avg_abs_frac_new, RSP_kap_density_factor, RSP_map_columns_filename, &
   RSP_relax_alfap_before_alfat, RSP_max_outer_dm_tries, RSP_max_inner_scale_tries, RSP_T_anchor_tolerance, &
    ! mass gain or loss  
    mass_change, mass_change_full_on_dt, mass_change_full_off_dt, trace_dt_control_mass_change, &
    min_wind, max_wind, use_accreted_material_j, accreted_material_j, D_omega_mixing_rate, &
    D_omega_mixing_across_convection_boundary, max_q_for_D_omega_zero_in_convection_region, nu_omega_mixing_rate, &
    nu_omega_mixing_across_convection_boundary, max_q_for_nu_omega_zero_in_convection_region, &
    mdot_omega_power, max_rotational_mdot_boost, max_mdot_jump_for_rotation, &
    lim_trace_rotational_mdot_boost, rotational_mdot_boost_fac, rotational_mdot_kh_fac, surf_avg_tau, surf_avg_tau_min, &
    max_tries_for_implicit_wind, iwind_tolerance, iwind_lambda, super_eddington_scaling_factor, &
    super_eddington_wind_Ledd_factor, wind_boost_full_off_L_div_Ledd, wind_boost_full_on_L_div_Ledd, &
    super_eddington_wind_max_boost, trace_super_eddington_wind_boost, &
    rlo_scaling_factor, rlo_wind_min_L, rlo_wind_max_Teff, rlo_wind_roche_lobe_radius, &
    roche_lobe_xfer_full_on, roche_lobe_xfer_full_off, rlo_wind_base_mdot, rlo_wind_scale_height, &
    hot_wind_scheme, cool_wind_RGB_scheme, cool_wind_AGB_scheme, RGB_to_AGB_wind_switch, &
    Reimers_scaling_factor, Blocker_scaling_factor, de_Jager_scaling_factor, van_Loon_scaling_factor, &
    Nieuwenhuijzen_scaling_factor, Vink_scaling_factor, &
    Dutch_scaling_factor, Dutch_wind_lowT_scheme, wind_He_layer_limit, &
    wind_H_envelope_limit, wind_H_He_envelope_limit, hot_wind_full_on_T, cool_wind_full_on_T, &
    
    ! composition of added mass
    accrete_same_as_surface, &
    accrete_given_mass_fractions, num_accretion_species, accretion_species_id, accretion_species_xa, &
    accretion_h1, accretion_h2, accretion_he3, accretion_he4, accretion_zfracs, accretion_dump_missing_metals_into_heaviest, &
    
    ! special list of z fractions
    z_fraction_li, z_fraction_be, z_fraction_b, z_fraction_c, z_fraction_n,&
    z_fraction_o, z_fraction_f, z_fraction_ne, z_fraction_na, z_fraction_mg, z_fraction_al, &
    z_fraction_si, z_fraction_p, z_fraction_s, z_fraction_cl, z_fraction_ar, z_fraction_k, &
    z_fraction_ca, z_fraction_sc, z_fraction_ti, z_fraction_v, z_fraction_cr, z_fraction_mn, &
    z_fraction_fe, z_fraction_co, z_fraction_ni, z_fraction_cu, z_fraction_zn, &
    lgT_lo_for_set_new_abundances, lgT_hi_for_set_new_abundances, &
    
    ! automatic stops for mass loss/gain
    max_star_mass_for_gain, min_star_mass_for_loss, max_T_center_for_any_mass_loss, max_T_center_for_full_mass_loss, &
    
    ! extra power source
    extra_power_source, &
    
    ! relaxation parameters
    relax_dlnZ, relax_dY, &
    
    ! mesh adjustment
    show_mesh_changes, okay_to_remesh, restore_mesh_on_retry, num_steps_to_hold_mesh_after_retry, &
    max_rel_delta_IE_for_mesh_total_energy_balance, &
    trace_mesh_adjust_error_in_conservation, max_allowed_nz, mesh_max_allowed_ratio, &
    remesh_max_allowed_logT, max_delta_x_for_merge, &
    mesh_ok_to_merge, mesh_max_k_old_for_split, mesh_min_k_old_for_split, &
    mesh_adjust_get_T_from_E, &
    max_dq, min_dq, min_dq_for_split, min_dq_for_xa, min_dq_for_xa_convective, &
    mesh_min_dlnR, merge_if_dlnR_too_small, min_dq_for_logT, &
    mesh_min_dr_div_dRstar, merge_if_dr_div_dRstar_too_small, &
    mesh_min_dr_div_cs, merge_if_dr_div_cs_too_small, &
    max_center_cell_dq, max_surface_cell_dq, max_num_subcells, max_num_merge_cells, &
    mesh_delta_coeff, mesh_delta_coeff_for_highT, &
    logT_max_for_standard_mesh_delta_coeff, logT_min_for_highT_mesh_delta_coeff, remesh_dt_limit, &
    mesh_Pgas_div_P_exponent, &
    E_function_weight, E_function_param, P_function_weight, &
    mesh_logX_species, &
    mesh_logX_min_for_extra, mesh_dlogX_dlogP_extra, mesh_dlogX_dlogP_full_on, mesh_dlogX_dlogP_full_off, &

    mesh_dlog_eps_min_for_extra, mesh_dlog_eps_dlogP_full_on, mesh_dlog_eps_dlogP_full_off, &
    mesh_dlog_pp_dlogP_extra, mesh_dlog_cno_dlogP_extra, mesh_dlog_3alf_dlogP_extra, &
    mesh_dlog_burn_c_dlogP_extra, mesh_dlog_burn_n_dlogP_extra, mesh_dlog_burn_o_dlogP_extra, &
    mesh_dlog_burn_ne_dlogP_extra, mesh_dlog_burn_na_dlogP_extra, mesh_dlog_burn_mg_dlogP_extra, &
    mesh_dlog_burn_si_dlogP_extra, mesh_dlog_burn_s_dlogP_extra, mesh_dlog_burn_ar_dlogP_extra, &
    mesh_dlog_burn_ca_dlogP_extra, mesh_dlog_burn_ti_dlogP_extra, mesh_dlog_burn_cr_dlogP_extra, &
    mesh_dlog_burn_fe_dlogP_extra, mesh_dlog_cc_dlogP_extra, mesh_dlog_co_dlogP_extra, mesh_dlog_oo_dlogP_extra, &
    mesh_dlog_pnhe4_dlogP_extra, mesh_dlog_photo_dlogP_extra, mesh_dlog_other_dlogP_extra, &
    T_function1_weight, T_function2_weight, T_function2_param, mesh_delta_coeff_factor_smooth_iters, &
    R_function_weight, R_function_param, &
    R_function2_weight, R_function2_param1, R_function2_param2, &
    R_function3_weight, M_function_weight, M_function_param, &
    gradT_function_weight, log_tau_function_weight, log_kap_function_weight, omega_function_weight, &
    gam_function_weight, gam_function_param1, gam_function_param2, &
    xa_function_species, xa_function_weight, xa_function_param, xa_mesh_delta_coeff, split_merge_amr_mesh_delta_coeff, &
    use_split_merge_amr, split_merge_amr_nz_baseline, split_merge_amr_log_zoning, split_merge_amr_hybrid_zoning, &
    split_merge_amr_flipped_hybrid_zoning, split_merge_amr_logtau_zoning, split_merge_amr_okay_to_split_nz, split_merge_amr_nz_r_core, &
    split_merge_amr_okay_to_split_1, merge_amr_inhibit_at_jumps, split_merge_amr_MaxLong, split_merge_amr_nz_r_core_fraction, &
    split_merge_amr_MaxShort, merge_amr_max_abs_du_div_cs, &
    merge_amr_ignore_surface_cells, merge_amr_k_for_ignore_surface_cells, &
    merge_amr_du_div_cs_limit_only_for_compression, split_merge_amr_avoid_repeated_remesh, split_merge_amr_r_core_cm, &
    split_merge_amr_dq_min, split_merge_amr_dq_max, split_merge_amr_max_iters, trace_split_merge_amr, equal_split_density_amr, &

    ! nuclear reaction parameters
    screening_mode, default_net_name, net_logTcut_lo, net_logTcut_lim, &
    eps_nuc_factor, op_split_burn_eps_nuc_infall_limit, eps_WD_sedimentation_factor, &
    max_abs_eps_nuc, dxdt_nuc_factor, max_abar_for_burning, mix_factor, &
    fe56ec_fake_factor, min_T_for_fe56ec_fake_factor, weak_rate_factor, &
    sig_term_limit, sig_min_factor_for_high_Tcenter, Tcenter_min_for_sig_min_factor_full_on,&
    Tcenter_max_for_sig_min_factor_full_off, max_delta_m_to_bdy_for_sig_min_factor, &
    delta_m_lower_for_sig_min_factor, delta_m_upper_for_sig_min_factor, &
    am_sig_term_limit, am_D_mix_factor, am_gradmu_factor, am_nu_factor, &
    D_visc_factor, D_DSI_factor, D_SH_factor, D_SSI_factor, D_ES_factor, D_GSF_factor, D_ST_factor, &
    am_nu_non_rotation_factor, skip_rotation_in_convection_zones, am_nu_DSI_factor, am_nu_SH_factor,&
    am_nu_SSI_factor, am_nu_ES_factor, am_nu_GSF_factor, am_nu_ST_factor, am_nu_visc_factor, smooth_am_nu_rot, &
    ST_angsml, ST_angsmt, am_nu_omega_rot_factor, am_nu_omega_non_rot_factor, am_nu_j_rot_factor, am_nu_j_non_rot_factor, &
    smooth_nu_ST, smooth_D_ST, smooth_D_SH, smooth_D_DSI, smooth_D_ES, smooth_D_SSI, smooth_D_GSF, smooth_D_omega, &
    do_adjust_J_lost, premix_omega, angular_momentum_error_warn, angular_momentum_error_retry, &
    simple_i_rot_flag, recalc_mixing_info_each_substep, adjust_J_fraction, &
    min_q_for_adjust_J_lost, min_J_div_delta_J, max_mdot_redo_cnt, mdot_revise_factor, &
    implicit_mdot_boost, min_years_dt_for_redo_mdot, surf_omega_div_omega_crit_limit, surf_omega_div_omega_crit_tol, &
    w_div_wcrit_max, w_div_wcrit_max2, &
    fp_min, ft_min, fp_error_limit, ft_error_limit, &
    D_mix_rotation_max_logT_full_on, D_mix_rotation_min_logT_full_off, &
    set_uniform_am_nu_non_rot, uniform_am_nu_non_rot, &
    set_min_am_nu_non_rot, min_am_nu_non_rot, min_center_Ye_for_min_am_nu_non_rot, &
    set_min_D_mix_below_Tmax, min_D_mix_below_Tmax, set_min_D_mix_in_H_He, min_D_mix_in_H_He, &
    set_min_D_mix, mass_lower_limit_for_min_D_mix, mass_upper_limit_for_min_D_mix, &
    min_D_mix, min_center_Ye_for_min_D_mix, &
    smooth_outer_xa_big, smooth_outer_xa_small, nonlocal_NiCo_kap_gamma, nonlocal_NiCo_decay_heat, &
    dtau_gamma_NiCo_decay_heat, max_logT_for_net, reaction_neuQs_factor, &
    
    ! element diffusion parameters
    diffusion_use_iben_macdonald, diffusion_use_paquette, diffusion_use_cgs_solver, &
    diffusion_use_full_net, do_WD_sedimentation_heating, min_xa_for_WD_sedimentation_heating, &
    do_diffusion_heating, do_element_diffusion, &
    cgs_thermal_diffusion_eta_full_on, cgs_thermal_diffusion_eta_full_off, diffusion_min_dq_at_surface, &
    diffusion_min_T_at_surface, diffusion_min_dq_ratio_at_surface, diffusion_dt_limit, &
    diffusion_min_X_hard_limit, diffusion_X_total_atol, diffusion_X_total_rtol, &
    diffusion_upwind_abs_v_limit, diffusion_dt_div_timescale, diffusion_min_num_substeps, &
    diffusion_max_iters_per_substep, diffusion_max_retries_per_substep, diffusion_v_max, &
    diffusion_gamma_full_off, diffusion_gamma_full_on, diffusion_T_full_on, diffusion_T_full_off, &
    D_mix_ignore_diffusion, diffusion_calculates_ionization, diffusion_nsmooth_typical_charge, &
    diffusion_tol_correction_max, diffusion_tol_correction_norm, &
    diffusion_AD_dm_full_on, diffusion_AD_dm_full_off, diffusion_AD_boost_factor, &
    diffusion_SIG_factor, diffusion_GT_factor, &
    diffusion_Vlimit_dm_full_on, diffusion_Vlimit_dm_full_off, diffusion_Vlimit, &
    diffusion_max_T_for_radaccel, diffusion_min_T_for_radaccel, diffusion_max_Z_for_radaccel, &
    diffusion_min_Z_for_radaccel, diffusion_screening_for_radaccel, &
    op_mono_data_path, op_mono_data_cache_filename, &
    show_diffusion_info, show_diffusion_substep_info, show_diffusion_timing, &
    diffusion_num_classes, diffusion_class_representative, diffusion_class_A_max, &
    diffusion_class_typical_charge, diffusion_class_factor, &
    diffusion_use_isolve, diffusion_rtol_for_isolve, diffusion_atol_for_isolve, &
    diffusion_maxsteps_for_isolve, diffusion_isolve_solver, &

    ! WD phase separation
    do_phase_separation, &
    do_phase_separation_heating, &
    phase_separation_mixing_use_brunt, &
    phase_separation_no_diffusion, &
    
    ! eos controls
    fix_d_eos_dxa_partials, &

    ! opacity controls
    use_simple_es_for_kap, use_starting_composition_for_kap, &
    min_kap_for_dPrad_dm_eqn, low_logT_op_mono_full_off, low_logT_op_mono_full_on, high_logT_op_mono_full_off, &
    high_logT_op_mono_full_on, op_mono_min_X_to_include, use_op_mono_alt_get_kap, &
    
    
    include_L_in_correction_limits, include_v_in_correction_limits, include_u_in_correction_limits, include_w_in_correction_limits, &
    
    ! asteroseismology controls
    get_delta_nu_from_scaled_solar, nu_max_sun, delta_nu_sun, Teff_sun, delta_Pg_mode_freq, &
    
    ! hydro parameters
    energy_eqn_option, &
    opacity_factor, opacity_max, min_logT_for_opacity_factor_off, min_logT_for_opacity_factor_on, &
    max_logT_for_opacity_factor_on, max_logT_for_opacity_factor_off, &
    non_nuc_neu_factor, &
    use_time_centered_eps_grav, &
    use_mass_corrections, use_gravity_rotation_correction, eps_grav_factor, eps_mdot_factor, &
    include_composition_in_eps_grav, no_dedt_form_during_relax, &
    max_abs_rel_change_surf_lnS, &
    max_num_surf_revisions, Gamma_lnS_eps_grav_full_off, Gamma_lnS_eps_grav_full_on, &
    use_dPrad_dm_form_of_T_gradient_eqn, use_gradT_actual_vs_gradT_MLT_for_T_gradient_eqn, dedt_eqn_r_scale, &
    RTI_A, RTI_B, RTI_C, RTI_D, RTI_max_alpha, RTI_C_X_factor, RTI_C_X0_frac, steps_before_use_velocity_time_centering, &
    RTI_dm_for_center_eta_nondecreasing, RTI_min_dm_behind_shock_for_full_on, RTI_energy_floor, &
    RTI_D_mix_floor, RTI_min_m_for_D_mix_floor, RTI_log_max_boost, RTI_m_full_boost, RTI_m_no_boost, &
    include_P_in_velocity_time_centering, include_L_in_velocity_time_centering, &
    P_theta_for_velocity_time_centering, L_theta_for_velocity_time_centering, &
    steps_before_use_TDC, use_P_d_1_div_rho_form_of_work_when_time_centering_velocity, compare_TDC_to_MLT, &
    velocity_logT_lower_bound, max_dt_yrs_for_velocity_logT_lower_bound, velocity_q_upper_bound, &
    retry_for_v_above_clight, &

    ! hydro solver
    use_gold2_tolerances, gold2_solver_iters_timestep_limit, steps_before_use_gold2_tolerances, &
    gold2_iter_for_resid_tol2, gold2_iter_for_resid_tol3, gold2_tol_residual_norm1, gold2_tol_max_residual1, &
    gold2_tol_residual_norm2, gold2_tol_max_residual2, gold2_tol_residual_norm3, gold2_tol_max_residual3, &
    tol_correction_norm, tol_max_correction, correction_xa_limit, &
    tol_correction_high_T_limit, tol_correction_norm_high_T, tol_max_correction_high_T, &
    tol_correction_extreme_T_limit, tol_correction_norm_extreme_T, tol_max_correction_extreme_T, &
    tol_bad_max_correction, bad_max_correction_series_limit, &
    tol_max_residual1, tol_residual_norm1, tol_max_residual2, &
    tol_residual_norm2, tol_max_residual3, tol_residual_norm3, &
    warning_limit_for_max_residual, trace_solver_damping, &
    relax_use_gold_tolerances, relax_tol_correction_norm, relax_tol_max_correction, &
    relax_solver_iters_timestep_limit, &
    relax_iter_for_resid_tol2, relax_tol_residual_norm1, relax_tol_max_residual1, &
    relax_iter_for_resid_tol3, relax_tol_residual_norm2, relax_tol_max_residual2, &
    relax_tol_residual_norm3, relax_tol_max_residual3, relax_maxT_for_gold_tolerances, &
    use_gold_tolerances, gold_solver_iters_timestep_limit, maxT_for_gold_tolerances, &
    gold_tol_residual_norm1, gold_tol_max_residual1, gold_iter_for_resid_tol2, &
    gold_tol_residual_norm2, gold_tol_max_residual2, gold_iter_for_resid_tol3, &
    gold_tol_residual_norm3, gold_tol_max_residual3 , steps_before_use_gold_tolerances, &
    include_rotation_in_total_energy, convergence_ignore_equL_residuals, convergence_ignore_alpha_RTI_residuals, &
    iter_for_resid_tol2, iter_for_resid_tol3, &
    solver_itermin, solver_itermin_until_reduce_min_corr_coeff, &
    solver_reduced_min_corr_coeff, do_solver_damping_for_neg_xa, &
    scale_max_correction_for_negative_surf_lum, max_frac_for_negative_surf_lum, &
    hydro_mtx_max_allowed_abs_dlogT, hydro_mtx_max_allowed_abs_dlogRho, &
    min_logT_for_hydro_mtx_max_allowed, hydro_mtx_max_allowed_logT, &
    hydro_mtx_max_allowed_logRho, report_min_rcond_from_DGESXV, &
    hydro_mtx_min_allowed_logT, hydro_mtx_min_allowed_logRho, use_DGESVX_in_bcyclic, use_equilibration_in_DGESVX, &
    op_split_burn, op_split_burn_min_T, op_split_burn_eps, op_split_burn_odescal, &
    op_split_burn_min_T_for_variable_T_solver, solver_test_partials_show_dx_var_name, &
    tiny_corr_coeff_limit, scale_correction_norm, corr_param_factor, num_times_solver_reuse_mtx, &
    scale_max_correction, ignore_min_corr_coeff_for_scale_max_correction, ignore_too_large_correction, ignore_species_in_max_correction, &
    corr_norm_jump_limit, max_corr_jump_limit, resid_norm_jump_limit, max_resid_jump_limit, RSP2_use_mass_interp_face_values, &
    corr_coeff_limit, tiny_corr_factor, solver_test_partials_call_number, solver_test_partials_iter_number, &
    max_tries1, solver_max_tries_before_reject, max_tries_for_retry, max_tries_after_5_retries, solver_test_partials_sink_name, &
    max_tries_after_10_retries, max_tries_after_20_retries, retry_limit, redo_limit, use_Pvsc_art_visc, Pvsc_cq, Pvsc_zsh, &
    min_xa_hard_limit, min_xa_hard_limit_for_highT, logT_max_for_min_xa_hard_limit, logT_min_for_min_xa_hard_limit_for_highT, &
    sum_xa_hard_limit, sum_xa_hard_limit_for_highT, logT_max_for_sum_xa_hard_limit, logT_min_for_sum_xa_hard_limit_for_highT, &
    xa_clip_limit, report_solver_progress, solver_test_partials_k_high, RSP2_use_L_eqn_at_surface, RSP2_use_RSP_eqn_for_Y_face, &
    solver_epsder_chem, solver_epsder_struct, solver_numerical_jacobian, energy_conservation_dump_model_number, &
    solver_jacobian_nzlo, solver_jacobian_nzhi, solver_check_everything, solver_inspect_soln_flag, RSP2_assume_HSE, &
    solver_test_partials_dx_0, solver_test_partials_k, solver_show_correction_info, eps_mdot_leak_frac_factor, &
    solver_test_partials_write_eos_call_info, solver_save_photo_call_number, RSP2_min_Lc_div_L_for_convective_mixing_type, &
    solver_test_partials_var_name, solver_test_partials_equ_name, RSP2_min_Lt_div_L_for_overshooting_mixing_type, &
    solver_test_eos_partials, solver_test_kap_partials, solver_test_net_partials, solver_test_atm_partials, &
    fill_arrays_with_NaNs, zero_when_allocate, warn_when_large_rel_run_E_err, solver_test_partials_k_low, &
    warn_when_large_virial_thm_rel_err, warn_when_get_a_bad_eos_result, warn_rates_for_high_temp, max_safe_logT_for_rates, &
    RSP2_alfap, RSP2_alfat, RSP2_alfam, RSP2_alfar, RSP2_Lsurf_factor, RSP2_use_Stellingwerf_Lr, RSP2_remesh_when_load, &
    RSP2_alfad, RSP2_num_outermost_cells_forced_nonturbulent, RSP2_num_innermost_cells_forced_nonturbulent, &
    RSP2_target_steps_per_cycle, RSP2_max_num_periods, RSP2_work_period, RSP2_map_first_period, RSP2_map_last_period, &
    RSP2_min_max_R_for_periods, RSP2_GREKM_avg_abs_frac_new, RSP2_GREKM_avg_abs_limit, RSP2_map_zone_interval, &
    RSP2_work_filename, RSP2_map_columns_filename, RSP2_map_filename, RSP2_map_history_filename, RSP2_write_map, &
    RSP2_T_anchor, RSP2_dq_1_factor, RSP2_nz, RSP2_nz_outer, RSP2_nz_div_IBOTOM, RSP2_report_adjust_w, &
    RSP2_w_min_for_damping, RSP2_source_seed, RSP2_w_fix_if_neg, max_X_for_conv_timescale, min_X_for_conv_timescale, &
    max_q_for_conv_timescale, min_q_for_conv_timescale, max_q_for_QHSE_timescale, min_q_for_QHSE_timescale, &
    
    
    ! timestep
    time_delta_coeff, min_timestep_factor, max_timestep_factor, timestep_factor_for_retries, retry_hold, &
    neg_mass_fraction_hold, timestep_dt_factor, use_dt_low_pass_controller, &
    force_timestep_min, force_timestep_min_years, force_timestep_min_factor, force_timestep, force_timestep_years, &
    varcontrol_target, min_allowed_varcontrol_target, varcontrol_dt_limit_ratio_hard_max, xa_scale, &
    solver_iters_timestep_limit, burn_steps_limit, burn_steps_hard_limit, &
    diffusion_steps_limit, diffusion_steps_hard_limit, diffusion_iters_limit, diffusion_iters_hard_limit, &
    dt_div_dt_cell_collapse_limit, dt_div_dt_cell_collapse_hard_limit, &
    dt_div_min_dr_div_cs_limit, dt_div_min_dr_div_cs_hard_limit, min_abs_u_div_cs_for_dt_div_min_dr_div_cs_limit, &
    min_abs_du_div_cs_for_dt_div_min_dr_div_cs_limit, min_k_for_dt_div_min_dr_div_cs_limit, &
    min_q_for_dt_div_min_dr_div_cs_limit, max_q_for_dt_div_min_dr_div_cs_limit, check_remnant_only_for_dt_div_min_dr_div_cs_limit, &
    dX_mix_dist_limit, dH_limit_min_H, dH_limit, dH_hard_limit, dH_div_H_limit_min_H, &
    dH_div_H_limit, dH_div_H_hard_limit, dH_decreases_only, max_timestep_factor_at_high_T, &
    dHe_limit_min_He, dHe_limit, dHe_hard_limit, dHe_div_He_limit_min_He, &
    dHe_div_He_limit, dHe_div_He_hard_limit, dHe_decreases_only, min_logT_for_max_timestep_factor_at_high_T, &
    dHe3_limit_min_He3, dHe3_limit, dHe3_hard_limit, dHe3_div_He3_limit_min_He3, &
    dHe3_div_He3_limit, dHe3_div_He3_hard_limit, dHe3_decreases_only, dX_div_X_at_high_T_limit_lgT_min, &
    dX_limit_min_X, dX_limit, dX_hard_limit, dX_div_X_limit_min_X, dX_div_X_at_high_T_hard_limit, &
    dX_div_X_limit, dX_div_X_hard_limit, dX_decreases_only, dX_div_X_at_high_T_limit, &
    dX_nuc_drop_min_X_limit, dX_nuc_drop_max_A_limit, dX_nuc_drop_limit, &
    dX_nuc_drop_limit_at_high_T, dX_nuc_drop_hard_limit, dX_nuc_drop_min_yrs_for_dt, &
    dL_div_L_limit_min_L, dL_div_L_limit, dL_div_L_hard_limit, &
    delta_lgP_limit, delta_lgP_hard_limit, delta_lgP_limit_min_lgP, &
    delta_lgRho_limit, delta_lgRho_hard_limit, delta_lgRho_limit_min_lgRho, &
    delta_lgT_limit, delta_lgT_hard_limit, delta_lgT_limit_min_lgT, &
    delta_lgE_limit, delta_lgE_hard_limit, delta_lgE_limit_min_lgE, &
    delta_lgR_limit, delta_lgR_hard_limit, delta_lgR_limit_min_lgR, &
    delta_Ye_highT_limit, delta_Ye_highT_hard_limit, minT_for_highT_Ye_limit, &
    delta_lgL_nuc_cat_limit, delta_lgL_nuc_cat_hard_limit, lgL_nuc_cat_burn_min, lgL_nuc_mix_dist_limit, &
    delta_lgT_max_limit_only_after_near_zams, &
    delta_lgL_H_limit, delta_lgL_H_hard_limit, lgL_H_burn_min, lgL_H_drop_factor, lgL_H_burn_relative_limit, &
    delta_lgL_He_limit, delta_lgL_He_hard_limit, lgL_He_burn_min, lgL_He_drop_factor, lgL_He_burn_relative_limit, &
    delta_lgL_z_limit, delta_lgL_z_hard_limit, lgL_z_burn_min, lgL_z_drop_factor, lgL_z_burn_relative_limit, &
    delta_lgL_power_photo_limit, delta_lgL_power_photo_hard_limit, lgL_power_photo_burn_min, lgL_power_photo_drop_factor, &
    delta_lgL_nuc_limit, delta_lgL_nuc_hard_limit, lgL_nuc_burn_min, lgL_nuc_drop_factor, min_lgT_for_lgL_power_photo_limit, &
    delta_lgL_nuc_at_high_T_limit, delta_lgL_nuc_at_high_T_hard_limit, delta_lgL_nuc_at_high_T_limit_lgT_min, &
    delta_lgRho_cntr_limit, delta_lgRho_cntr_hard_limit, delta_lgP_cntr_limit, delta_lgP_cntr_hard_limit, &
    delta_lgT_cntr_limit, delta_lgT_cntr_hard_limit, delta_lgT_cntr_limit_only_after_near_zams, &
    delta_lgT_max_limit, delta_lgT_max_hard_limit, max_lgT_for_lgL_nuc_limit, &
    delta_log_eps_nuc_limit, delta_log_eps_nuc_hard_limit, delta_lgT_max_limit_lgT_min, &
    delta_lgT_max_at_high_T_limit, delta_lgT_max_at_high_T_hard_limit, delta_lgT_max_at_high_T_limit_lgT_min, &
    delta_dX_div_X_cntr_min, delta_dX_div_X_cntr_max, delta_dX_div_X_cntr_limit, delta_dX_div_X_cntr_hard_limit, &
    delta_dX_div_X_drop_only, delta_lg_XH_drop_only, &
    delta_lg_XHe_drop_only, delta_lg_XC_drop_only, delta_lg_XNe_drop_only, delta_lg_XO_drop_only, delta_lg_XSi_drop_only, &
    delta_XH_drop_only, delta_XHe_drop_only, delta_XC_drop_only, delta_XNe_drop_only, delta_XO_drop_only, delta_XSi_drop_only, &    
    delta_lg_XH_cntr_min, delta_lg_XH_cntr_max, delta_lg_XH_cntr_limit, delta_lg_XH_cntr_hard_limit, &
    delta_lg_XHe_cntr_min, delta_lg_XHe_cntr_max, delta_lg_XHe_cntr_limit, delta_lg_XHe_cntr_hard_limit, &
    delta_lg_XC_cntr_min, delta_lg_XC_cntr_max, delta_lg_XC_cntr_limit, delta_lg_XC_cntr_hard_limit, &
    delta_lg_XNe_cntr_limit, delta_lg_XNe_cntr_hard_limit, delta_lg_XNe_cntr_min, delta_lg_XNe_cntr_max, &
    delta_lg_XO_cntr_limit, delta_lg_XO_cntr_hard_limit, delta_lg_XO_cntr_min, delta_lg_XO_cntr_max, &
    delta_lg_XSi_cntr_limit, delta_lg_XSi_cntr_hard_limit, delta_lg_XSi_cntr_min, delta_lg_XSi_cntr_max, &
    delta_XH_cntr_limit, delta_XH_cntr_hard_limit, delta_XHe_cntr_limit, delta_XHe_cntr_hard_limit, &
    delta_XC_cntr_limit, delta_XC_cntr_hard_limit, delta_XNe_cntr_limit, delta_XNe_cntr_hard_limit, &
    delta_XO_cntr_limit, delta_XO_cntr_hard_limit, delta_XSi_cntr_limit, delta_XSi_cntr_hard_limit, &
    delta_lgTeff_limit, delta_lgTeff_hard_limit, delta_lgL_limit, delta_lgL_limit_L_min, delta_lgL_hard_limit, &
    delta_HR_ds_L, delta_HR_ds_Teff, delta_HR_limit, delta_HR_hard_limit, &
    delta_lg_star_mass_limit, delta_lg_star_mass_hard_limit, &
    delta_mdot_atol, delta_mdot_rtol, delta_mdot_limit, delta_mdot_hard_limit, &
    adjust_J_q_limit, adjust_J_q_hard_limit, &
    never_skip_hard_limits, relax_hard_limits_after_retry, &
    report_dt_hard_limit_retries, report_min_dr_div_cs, report_solver_dt_info, &
    limit_for_rel_error_in_energy_conservation, hard_limit_for_rel_error_in_energy_conservation, &

    ! atmosphere -- surface boundary conditions

    atm_option, atm_off_table_option, Pextra_factor, &
    atm_fixed_Teff, atm_fixed_Psurf, atm_fixed_Tsurf, &

    atm_T_tau_relation, atm_T_tau_opacity, atm_T_tau_errtol, atm_T_tau_max_iters, &
    atm_T_tau_max_steps, &

    atm_table, &

    atm_irradiated_opacity, atm_irradiated_errtol, atm_irradiated_T_eq, &
    atm_irradiated_kap_v, atm_irradiated_kap_v_div_kap_th, atm_irradiated_P_surf, &
    atm_irradiated_max_iters, &

    use_compression_outer_BC, use_momentum_outer_BC, Tsurf_factor, use_zero_Pgas_outer_BC, &
    fixed_Psurf, use_fixed_Psurf_outer_BC, fixed_vsurf, use_fixed_vsurf_outer_BC, &
    
    atm_build_tau_outer, atm_build_dlogtau, atm_build_errtol, &

    use_T_tau_gradr_factor, &
    
    ! extra heat near surface to model irradiation
    irradiation_flux, column_depth_for_irradiation, &
    
    ! uniform extra heat
    inject_uniform_extra_heat, min_q_for_uniform_extra_heat, max_q_for_uniform_extra_heat, &
    inject_extra_ergs_sec, base_of_inject_extra_ergs_sec, total_mass_for_inject_extra_ergs_sec, &
    start_time_for_inject_extra_ergs_sec, duration_for_inject_extra_ergs_sec, &
    inject_until_reach_model_with_total_energy, &
    
    ! mass gain or loss
    no_wind_if_no_rotation, max_logT_for_k_below_const_q, &
    max_q_for_k_below_const_q, min_q_for_k_below_const_q, max_logT_for_k_const_mass, &
    min_q_for_k_const_mass, max_q_for_k_const_mass, &
    
    ! info for debugging
    stop_for_bad_nums, report_ierr, report_bad_negative_xa, diffusion_dump_call_number, &
     
    ! controls for the evolve routine
    trace_evolve, &

    ! misc
    min_chem_eqn_scale, zams_filename, set_rho_to_dm_div_dV, use_other_momentum_implicit, &
    use_other_surface_PT, use_other_mlt_results, use_other_kap, use_other_diffusion, use_other_diffusion_factor, &
    use_other_adjust_mdot, use_other_j_for_adjust_J_lost, use_other_alpha_mlt, use_other_remove_surface, &
    use_other_am_mixing, use_other_brunt, use_other_brunt_smoothing, use_other_solver_monitor, &
    use_other_build_initial_model, use_other_cgrav, use_other_energy_implicit, use_other_momentum, &
    use_other_energy, use_other_mesh_functions, use_other_eps_grav, use_other_gradr_factor, &
    use_other_D_mix, use_other_neu, use_other_net_get, use_other_opacity_factor, use_other_pressure, &
    use_other_diffusion_coefficients, use_other_pgstar_plots, use_other_eval_fp_ft, use_other_eval_i_rot, use_other_torque, &
    use_other_torque_implicit, use_other_wind, use_other_accreting_state, use_other_after_struct_burn_mix, use_other_mesh_delta_coeff_factor, &
    use_other_before_struct_burn_mix, use_other_astero_freq_corr, use_other_timestep_limit, use_other_set_pgstar_controls, &
    use_other_screening, &
    x_ctrl, x_integer_ctrl, x_logical_ctrl, x_character_ctrl, &
    
    ! extra files
    read_extra_controls_inlist1, extra_controls_inlist1_name, read_extra_controls_inlist2, &
    extra_controls_inlist2_name, read_extra_controls_inlist3, extra_controls_inlist3_name, &
    read_extra_controls_inlist4, extra_controls_inlist4_name, read_extra_controls_inlist5, extra_controls_inlist5_name, &
    save_controls_namelist, controls_namelist_name


 contains


 subroutine write_controls(s, fname, ierr)
    type (star_info), pointer :: s
 character (len=*), intent(in) :: fname
 integer, intent(out) :: ierr 
 integer :: iounit
 character (len=256) :: filename
 
 ierr = 0
 filename = fname

 if (len_trim(filename) == 0) filename = 'dump_controls.txt'

 open(newunit=iounit, file=trim(filename), &
    action='write', status='replace', iostat=ierr)
 if (ierr /= 0) then
    write(*,*) 'failed to open ' // trim(filename)
    return
 endif
 
 call set_controls_for_writing(s, ierr)
 if (ierr /= 0) then
    close(iounit)
    return
 end if

 write(iounit, nml=controls, iostat=ierr)
 
 write(*,*) 'write controls namelist values to "' // trim(filename)//'"'
    
 close(iounit)

 end subroutine write_controls


 subroutine do_one_setup(id, inlist, ierr)
    use utils_lib
    character (len=*), intent(in) :: inlist
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    include 'formats'
    call get_star_ptr(id, s, ierr)
    if (ierr /= 0) return
    call set_default_controls
    call read_controls(id, inlist, ierr)
    if (ierr /= 0) return

   if (save_controls_namelist) then
      call write_controls(s, controls_namelist_name, ierr)
      if (ierr /= 0) then
         write(*,*) 'ierr from write_controls ' // &
            trim(controls_namelist_name)
         return
      end if
   end if

 end subroutine do_one_setup


 subroutine read_controls(id, filename, ierr)
 use star_private_def
 use utils_lib
 character(*), intent(in) :: filename
 integer, intent(in) :: id
 integer, intent(out) :: ierr

 type (star_info), pointer :: s
 ierr = 0
 call get_star_ptr(id, s, ierr)
 if (ierr /= 0) return

 call read_controls_file(s, filename, 1, ierr)
 call check_controls(s, ierr)

 end subroutine read_controls


 subroutine check_controls(s, ierr)
    type (star_info), pointer :: s
    integer, intent(out) :: ierr

    ierr = 0

    if (.not. (trim(s% ctrl% energy_eqn_option) == 'dedt' .or. trim(s% ctrl% energy_eqn_option) == 'eps_grav')) then
       write(*,'(A)')
       write(*,*) "Invalid choice for energy_eqn_option"
       write(*,*) "Available options are 'dedt' or 'eps_grav'"
       write(*,'(A)')
       ierr = -1
       return
    end if

 end subroutine check_controls


 recursive subroutine read_controls_file(s, filename, level, ierr)
 use star_private_def
 use utils_lib
 character(*), intent(in) :: filename
 type (star_info), pointer :: s
 integer, intent(in) :: level
 integer, intent(out) :: ierr
 logical :: read_extra1, read_extra2, read_extra3, read_extra4, read_extra5
 character (len=strlen) :: message, extra1, extra2, extra3, extra4, extra5
 integer :: unit

 ierr = 0

 if (level >= 10) then
 write(*,*) 'ERROR: too many levels of nested extra controls inlist files'
 ierr = -1
 return
 end if

 if (len_trim(filename) > 0) then
    open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
    if (ierr /= 0) then
       write(*, *) 'Failed to open control namelist file ', trim(filename)
       return
    end if
    read(unit, nml=controls, iostat=ierr)
    close(unit)
    if (ierr /= 0) then
       write(*, *)
       write(*, *)
       write(*, *)
       write(*, *)
       write(*, '(a)') 'Failed while trying to read control namelist file: ' // trim(filename)
       write(*, '(a)') 'Perhaps the following runtime error message will help you find the problem.'
       write(*, *)
       open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
       read(unit, nml=controls)
       close(unit)
       return
    end if
 end if

 call store_controls(s, ierr)

 ! recursive calls to read other inlists

 read_extra1 = read_extra_controls_inlist1
 read_extra_controls_inlist1 = .false.
 extra1 = extra_controls_inlist1_name
 extra_controls_inlist1_name = 'undefined'

 read_extra2 = read_extra_controls_inlist2
 read_extra_controls_inlist2 = .false.
 extra2 = extra_controls_inlist2_name
 extra_controls_inlist2_name = 'undefined'

 read_extra3 = read_extra_controls_inlist3
 read_extra_controls_inlist3 = .false.
 extra3 = extra_controls_inlist3_name
 extra_controls_inlist3_name = 'undefined'

 read_extra4 = read_extra_controls_inlist4
 read_extra_controls_inlist4 = .false.
 extra4 = extra_controls_inlist4_name
 extra_controls_inlist4_name = 'undefined'

 read_extra5 = read_extra_controls_inlist5
 read_extra_controls_inlist5 = .false.
 extra5 = extra_controls_inlist5_name
 extra_controls_inlist5_name = 'undefined'

 if (read_extra1) then
 write(*,*) 'read ' // trim(extra1)
 call read_controls_file(s, extra1, level+1, ierr)
 if (ierr /= 0) return
 end if

 if (read_extra2) then
 write(*,*) 'read ' // trim(extra2)
 call read_controls_file(s, extra2, level+1, ierr)
 if (ierr /= 0) return
 end if

 if (read_extra3) then
 write(*,*) 'read ' // trim(extra3)
 call read_controls_file(s, extra3, level+1, ierr)
 if (ierr /= 0) return
 end if

 if (read_extra4) then
    write(*,*) 'read ' // trim(extra4)
    call read_controls_file(s, extra4, level+1, ierr)
    if (ierr /= 0) return
 end if

 if (read_extra5) then
    write(*,*) 'read ' // trim(extra5)
    call read_controls_file(s, extra5, level+1, ierr)
    if (ierr /= 0) return
 end if

 end subroutine read_controls_file
 

 subroutine set_default_controls

 xa_central_lower_limit_species(:) = ''
 xa_central_lower_limit(:) = 0
 xa_central_upper_limit_species(:) = ''
 xa_central_upper_limit(:) = 0
 xa_surface_lower_limit_species(:) = ''
 xa_surface_lower_limit(:) = 0
 xa_surface_upper_limit_species(:) = ''
 xa_surface_upper_limit(:) = 0
 xa_average_lower_limit_species(:) = ''
 xa_average_lower_limit(:) = 0
 xa_average_upper_limit_species(:) = ''
 xa_average_upper_limit(:) = 0

 trace_history_value_name(:) = ''
 accretion_species_id(:) = ''
 accretion_species_xa(:) = 0
 mesh_logX_species(:) = ''
 mesh_logX_min_for_extra(:) = 0
 mesh_dlogX_dlogP_extra(:) = 0
 mesh_dlogX_dlogP_full_on(:) = 0
 mesh_dlogX_dlogP_full_off(:) = 0

 xa_function_species(:) = ''
 xa_function_weight(:) = 0
 xa_function_param(:) = 0
 xa_mesh_delta_coeff(:) = 0

 diffusion_class_representative(:) = ''
 diffusion_class_A_max(:) = 0
 diffusion_class_typical_charge(:) = 0
 diffusion_class_factor(:) = 0

 x_character_ctrl(:) = ''

 predictive_mix(:) = .FALSE.
 predictive_zone_type(:) = ''
 predictive_zone_loc(:) = ''
 predictive_bdy_loc(:) = ''
 predictive_bdy_q_min(:) = 0d0
 predictive_bdy_q_max(:) = 1d0
 predictive_superad_thresh(:) = 0d0
 predictive_avoid_reversal(:) = ''
 predictive_limit_ingestion(:) = ''
 predictive_ingestion_factor(:) = 0d0

 include 'controls.defaults'
 include 'controls_dev.defaults'

 end subroutine set_default_controls


 subroutine store_controls(s, ierr)
 use star_private_def
 use chem_def ! categories
 use utils_lib, only: mkdir
 type (star_info), pointer :: s
 integer, intent(out) :: ierr

 ierr = 0

 ! where to start
 s% ctrl% initial_mass = initial_mass
 s% ctrl% initial_z = initial_z
 s% ctrl% initial_y = initial_y
 s% ctrl% initial_he3 = initial_he3

 ! definition of core boundaries
 s% ctrl% he_core_boundary_h1_fraction = he_core_boundary_h1_fraction
 s% ctrl% co_core_boundary_he4_fraction = co_core_boundary_he4_fraction
 s% ctrl% one_core_boundary_he4_c12_fraction = one_core_boundary_he4_c12_fraction
 s% ctrl% fe_core_boundary_si28_fraction = fe_core_boundary_si28_fraction
 s% ctrl% neutron_rich_core_boundary_Ye_max = neutron_rich_core_boundary_Ye_max
 s% ctrl% min_boundary_fraction = min_boundary_fraction

 ! when to stop
 s% ctrl% max_model_number = max_model_number
 s% ctrl% max_number_retries = max_number_retries
 s% ctrl% max_abs_rel_run_E_err = max_abs_rel_run_E_err
 s% ctrl% relax_max_number_retries = relax_max_number_retries
 s% ctrl% max_age = max_age
 s% ctrl% max_age_in_days = max_age_in_days
 s% ctrl% max_age_in_seconds = max_age_in_seconds
 s% ctrl% num_adjusted_dt_steps_before_max_age = num_adjusted_dt_steps_before_max_age
 s% ctrl% dt_years_for_steps_before_max_age = dt_years_for_steps_before_max_age
 s% ctrl% reduction_factor_for_max_timestep = reduction_factor_for_max_timestep
 s% ctrl% when_to_stop_rtol = when_to_stop_rtol
 s% ctrl% when_to_stop_atol = when_to_stop_atol
 s% ctrl% gamma_center_limit = gamma_center_limit
 s% ctrl% eta_center_limit = eta_center_limit
 s% ctrl% log_center_temp_limit = log_center_temp_limit
 s% ctrl% log_max_temp_upper_limit = log_max_temp_upper_limit
 s% ctrl% log_max_temp_lower_limit = log_max_temp_lower_limit
 s% ctrl% log_center_temp_lower_limit = log_center_temp_lower_limit
 s% ctrl% log_center_density_limit = log_center_density_limit
 s% ctrl% log_center_density_lower_limit = log_center_density_lower_limit
 s% ctrl% min_timestep_limit = min_timestep_limit

 s% ctrl% center_entropy_limit = center_entropy_limit
 s% ctrl% center_entropy_lower_limit = center_entropy_lower_limit
 s% ctrl% max_entropy_limit = max_entropy_limit
 s% ctrl% max_entropy_lower_limit = max_entropy_lower_limit

 s% ctrl% fe_core_infall_limit = fe_core_infall_limit
 s% ctrl% center_Ye_lower_limit = center_Ye_lower_limit
 s% ctrl% center_R_lower_limit = center_R_lower_limit
 s% ctrl% non_fe_core_infall_limit = non_fe_core_infall_limit
 s% ctrl% non_fe_core_rebound_limit = non_fe_core_rebound_limit
 s% ctrl% v_div_csound_surf_limit = v_div_csound_surf_limit
 s% ctrl% v_div_csound_max_limit = v_div_csound_max_limit
 s% ctrl% Lnuc_div_L_upper_limit = Lnuc_div_L_upper_limit
 s% ctrl% Lnuc_div_L_lower_limit = Lnuc_div_L_lower_limit
 s% ctrl% v_surf_div_v_kh_upper_limit = v_surf_div_v_kh_upper_limit
 s% ctrl% v_surf_div_v_kh_lower_limit = v_surf_div_v_kh_lower_limit
 s% ctrl% v_surf_div_v_esc_limit = v_surf_div_v_esc_limit
 s% ctrl% v_surf_kms_limit = v_surf_kms_limit

 s% ctrl% stop_near_zams = stop_near_zams
 s% ctrl% stop_at_phase_PreMS = stop_at_phase_PreMS
 s% ctrl% stop_at_phase_ZAMS = stop_at_phase_ZAMS
 s% ctrl% stop_at_phase_IAMS = stop_at_phase_IAMS
 s% ctrl% stop_at_phase_TAMS = stop_at_phase_TAMS
 s% ctrl% stop_at_phase_He_Burn = stop_at_phase_He_Burn
 s% ctrl% stop_at_phase_ZACHeB = stop_at_phase_ZACHeB
 s% ctrl% stop_at_phase_TACHeB = stop_at_phase_TACHeB
 s% ctrl% stop_at_phase_TP_AGB = stop_at_phase_TP_AGB
 s% ctrl% stop_at_phase_C_Burn = stop_at_phase_C_Burn
 s% ctrl% stop_at_phase_Ne_Burn = stop_at_phase_Ne_Burn
 s% ctrl% stop_at_phase_O_Burn = stop_at_phase_O_Burn
 s% ctrl% stop_at_phase_Si_Burn = stop_at_phase_Si_Burn
 s% ctrl% stop_at_phase_WDCS = stop_at_phase_WDCS
 s% ctrl% Lnuc_div_L_zams_limit = Lnuc_div_L_zams_limit
 s% ctrl% gamma1_limit = gamma1_limit
 s% ctrl% gamma1_limit_max_q = gamma1_limit_max_q
 s% ctrl% gamma1_limit_max_v_div_vesc = gamma1_limit_max_v_div_vesc
 s% ctrl% Pgas_div_P_limit = Pgas_div_P_limit
 s% ctrl% Pgas_div_P_limit_max_q = Pgas_div_P_limit_max_q
 s% ctrl% peak_burn_vconv_div_cs_limit = peak_burn_vconv_div_cs_limit
 s% ctrl% omega_div_omega_crit_limit = omega_div_omega_crit_limit
 s% ctrl% delta_nu_lower_limit = delta_nu_lower_limit
 s% ctrl% delta_nu_upper_limit = delta_nu_upper_limit
 s% ctrl% delta_Pg_lower_limit = delta_Pg_lower_limit
 s% ctrl% delta_Pg_upper_limit = delta_Pg_upper_limit
 s% ctrl% shock_mass_upper_limit = shock_mass_upper_limit
 s% ctrl% mach1_mass_upper_limit = mach1_mass_upper_limit
 s% ctrl% stop_when_reach_this_cumulative_extra_heating = stop_when_reach_this_cumulative_extra_heating

 s% ctrl% xa_central_lower_limit_species = xa_central_lower_limit_species
 s% ctrl% xa_central_lower_limit = xa_central_lower_limit

 s% ctrl% xa_central_upper_limit_species = xa_central_upper_limit_species
 s% ctrl% xa_central_upper_limit = xa_central_upper_limit

 s% ctrl% xa_surface_lower_limit_species = xa_surface_lower_limit_species
 s% ctrl% xa_surface_lower_limit = xa_surface_lower_limit

 s% ctrl% xa_surface_upper_limit_species = xa_surface_upper_limit_species
 s% ctrl% xa_surface_upper_limit = xa_surface_upper_limit

 s% ctrl% xa_average_lower_limit_species = xa_average_lower_limit_species
 s% ctrl% xa_average_lower_limit = xa_average_lower_limit

 s% ctrl% xa_average_upper_limit_species = xa_average_upper_limit_species
 s% ctrl% xa_average_upper_limit = xa_average_upper_limit

 s% ctrl% HB_limit = HB_limit

 s% ctrl% star_mass_max_limit = star_mass_max_limit
 s% ctrl% star_mass_min_limit = star_mass_min_limit
 s% ctrl% ejecta_mass_max_limit = ejecta_mass_max_limit
 s% ctrl% remnant_mass_min_limit = remnant_mass_min_limit
 
 s% ctrl% star_species_mass_min_limit = star_species_mass_min_limit
 s% ctrl% star_species_mass_min_limit_iso = star_species_mass_min_limit_iso
 s% ctrl% star_species_mass_max_limit = star_species_mass_max_limit
 s% ctrl% star_species_mass_max_limit_iso = star_species_mass_max_limit_iso
 
 s% ctrl% xmstar_min_limit = xmstar_min_limit
 s% ctrl% xmstar_max_limit = xmstar_max_limit
 s% ctrl% envelope_mass_limit = envelope_mass_limit
 s% ctrl% envelope_fraction_left_limit = envelope_fraction_left_limit

 s% ctrl% he_core_mass_limit = he_core_mass_limit
 s% ctrl% co_core_mass_limit = co_core_mass_limit
 s% ctrl% one_core_mass_limit = one_core_mass_limit
 s% ctrl% fe_core_mass_limit = fe_core_mass_limit
 s% ctrl% neutron_rich_core_mass_limit = neutron_rich_core_mass_limit

 s% ctrl% he_layer_mass_lower_limit = he_layer_mass_lower_limit
 s% ctrl% abs_diff_lg_LH_lg_Ls_limit = abs_diff_lg_LH_lg_Ls_limit
 s% ctrl% Teff_upper_limit = Teff_upper_limit
 s% ctrl% Teff_lower_limit = Teff_lower_limit
 s% ctrl% photosphere_m_upper_limit = photosphere_m_upper_limit
 s% ctrl% photosphere_m_lower_limit = photosphere_m_lower_limit
 s% ctrl% photosphere_m_sub_M_center_limit = photosphere_m_sub_M_center_limit
 s% ctrl% photosphere_r_upper_limit = photosphere_r_upper_limit
 s% ctrl% photosphere_r_lower_limit = photosphere_r_lower_limit
 s% ctrl% log_Teff_upper_limit = log_Teff_upper_limit
 s% ctrl% log_Teff_lower_limit = log_Teff_lower_limit
 s% ctrl% log_Tsurf_upper_limit = log_Tsurf_upper_limit
 s% ctrl% log_Tsurf_lower_limit = log_Tsurf_lower_limit
 s% ctrl% log_Rsurf_upper_limit = log_Rsurf_upper_limit
 s% ctrl% log_Rsurf_lower_limit = log_Rsurf_lower_limit
 s% ctrl% log_Psurf_upper_limit = log_Psurf_upper_limit
 s% ctrl% log_Psurf_lower_limit = log_Psurf_lower_limit
 s% ctrl% log_Dsurf_upper_limit = log_Dsurf_upper_limit
 s% ctrl% log_Dsurf_lower_limit = log_Dsurf_lower_limit
 s% ctrl% log_L_upper_limit = log_L_upper_limit
 s% ctrl% log_L_lower_limit = log_L_lower_limit
 s% ctrl% log_g_upper_limit = log_g_upper_limit
 s% ctrl% log_g_lower_limit = log_g_lower_limit

 s% ctrl% power_nuc_burn_upper_limit = power_nuc_burn_upper_limit
 s% ctrl% power_h_burn_upper_limit = power_h_burn_upper_limit
 s% ctrl% power_he_burn_upper_limit = power_he_burn_upper_limit
 s% ctrl% power_z_burn_upper_limit = power_z_burn_upper_limit
 s% ctrl% power_nuc_burn_lower_limit = power_nuc_burn_lower_limit
 s% ctrl% power_h_burn_lower_limit = power_h_burn_lower_limit
 s% ctrl% power_he_burn_lower_limit = power_he_burn_lower_limit
 s% ctrl% power_z_burn_lower_limit = power_z_burn_lower_limit

 ! output of "snapshots" for restarts
 s% ctrl% photo_interval = photo_interval
 s% ctrl% photo_digits = photo_digits
 s% ctrl% photo_directory = photo_directory
 ! output of history and profiles.
 s% ctrl% do_history_file = do_history_file
 s% ctrl% history_interval = history_interval

 s% ctrl% write_header_frequency = write_header_frequency
 s% ctrl% terminal_interval = terminal_interval
 s% ctrl% terminal_show_age_units = terminal_show_age_units
 s% ctrl% terminal_show_timestep_units = terminal_show_timestep_units
 s% ctrl% terminal_show_log_dt = terminal_show_log_dt
 s% ctrl% terminal_show_log_age = terminal_show_log_age
 s% ctrl% extra_terminal_output_file = extra_terminal_output_file
 s% ctrl% num_trace_history_values = num_trace_history_values
 s% ctrl% trace_history_value_name = trace_history_value_name

 s% ctrl% log_directory = log_directory
 s% ctrl% star_history_name = star_history_name
 s% ctrl% star_history_header_name = star_history_header_name
 s% ctrl% star_history_dbl_format = star_history_dbl_format
 s% ctrl% star_history_int_format = star_history_int_format
 s% ctrl% star_history_txt_format = star_history_txt_format

 s% ctrl% profiles_index_name = profiles_index_name
 s% ctrl% profile_data_prefix = profile_data_prefix
 s% ctrl% profile_data_suffix = profile_data_suffix
 s% ctrl% profile_data_header_suffix = profile_data_header_suffix
 s% ctrl% profile_int_format = profile_int_format
 s% ctrl% profile_txt_format = profile_txt_format
 s% ctrl% profile_dbl_format = profile_dbl_format
 s% ctrl% profile_header_include_sys_details = profile_header_include_sys_details
 s% ctrl% write_profiles_flag = write_profiles_flag
 s% ctrl% profile_interval = profile_interval
 s% ctrl% priority_profile_interval = priority_profile_interval
 s% ctrl% profile_model = profile_model
 s% ctrl% max_num_profile_models = max_num_profile_models
 s% ctrl% max_num_profile_zones = max_num_profile_zones

 s% ctrl% write_controls_info_with_profile = write_controls_info_with_profile
 s% ctrl% controls_data_prefix = controls_data_prefix
 s% ctrl% controls_data_suffix = controls_data_suffix

 s% ctrl% write_pulse_data_with_profile = write_pulse_data_with_profile
 s% ctrl% pulse_data_format = pulse_data_format
 s% ctrl% add_atmosphere_to_pulse_data = add_atmosphere_to_pulse_data
 s% ctrl% add_center_point_to_pulse_data = add_center_point_to_pulse_data
 s% ctrl% keep_surface_point_for_pulse_data = keep_surface_point_for_pulse_data
 s% ctrl% add_double_points_to_pulse_data = add_double_points_to_pulse_data
 s% ctrl% interpolate_rho_for_pulse_data = interpolate_rho_for_pulse_data
 s% ctrl% threshold_grad_mu_for_double_point = threshold_grad_mu_for_double_point
 s% ctrl% max_number_of_double_points = max_number_of_double_points

 s% ctrl% fgong_header = fgong_header
 s% ctrl% fgong_ivers = fgong_ivers

 s% ctrl% max_num_gyre_points = max_num_gyre_points
 s% ctrl% format_for_OSC_data = format_for_OSC_data
 s% ctrl% fgong_zero_A_inside_r = fgong_zero_A_inside_r
 s% ctrl% use_other_export_pulse_data = use_other_export_pulse_data
 s% ctrl% use_other_get_pulse_data = use_other_get_pulse_data
 s% ctrl% use_other_edit_pulse_data = use_other_edit_pulse_data

 s% ctrl% write_model_with_profile = write_model_with_profile
 s% ctrl% model_data_prefix = model_data_prefix
 s% ctrl% model_data_suffix = model_data_suffix

 s% ctrl% mixing_D_limit_for_log = mixing_D_limit_for_log
 s% ctrl% trace_mass_location = trace_mass_location
 s% ctrl% min_tau_for_max_abs_v_location = min_tau_for_max_abs_v_location
 s% ctrl% min_q_for_inner_mach1_location = min_q_for_inner_mach1_location
 s% ctrl% max_q_for_outer_mach1_location = max_q_for_outer_mach1_location
 
 s% ctrl% mass_depth_for_L_surf = mass_depth_for_L_surf
 s% ctrl% conv_core_gap_dq_limit = conv_core_gap_dq_limit

 ! burn zone eps definitions for use in logs and profiles
 s% ctrl% burn_min1 = burn_min1
 s% ctrl% burn_min2 = burn_min2

 s% ctrl% max_conv_vel_div_csound_maxq = max_conv_vel_div_csound_maxq
 s% ctrl% width_for_limit_conv_vel = width_for_limit_conv_vel
 s% ctrl% max_q_for_limit_conv_vel = max_q_for_limit_conv_vel
 s% ctrl% max_mass_in_gm_for_limit_conv_vel = max_mass_in_gm_for_limit_conv_vel
 s% ctrl% max_r_in_cm_for_limit_conv_vel = max_r_in_cm_for_limit_conv_vel

 ! for reported average values
 s% ctrl% surface_avg_abundance_dq = surface_avg_abundance_dq
 s% ctrl% center_avg_value_dq = center_avg_value_dq

 ! mixing parameters
 s% ctrl% min_convective_gap = min_convective_gap
 s% ctrl% min_thermohaline_gap = min_thermohaline_gap
 s% ctrl% min_semiconvection_gap = min_semiconvection_gap
 s% ctrl% min_thermohaline_dropout = min_thermohaline_dropout
 s% ctrl% max_dropout_gradL_sub_grada = max_dropout_gradL_sub_grada
 s% ctrl% remove_embedded_semiconvection = remove_embedded_semiconvection
 s% ctrl% recalc_mix_info_after_evolve = recalc_mix_info_after_evolve
 s% ctrl% remove_mixing_glitches = remove_mixing_glitches
 s% ctrl% okay_to_remove_mixing_singleton = okay_to_remove_mixing_singleton
 s% ctrl% prune_bad_cz_min_Hp_height = prune_bad_cz_min_Hp_height
 s% ctrl% prune_bad_cz_min_log_eps_nuc = prune_bad_cz_min_log_eps_nuc
 s% ctrl% redo_conv_for_dr_lt_mixing_length = redo_conv_for_dr_lt_mixing_length

 s% ctrl% alpha_semiconvection = alpha_semiconvection
 s% ctrl% semiconvection_option = semiconvection_option
 s% ctrl% use_Ledoux_criterion = use_Ledoux_criterion
 s% ctrl% num_cells_for_smooth_gradL_composition_term = num_cells_for_smooth_gradL_composition_term
 s% ctrl% threshold_for_smooth_gradL_composition_term = threshold_for_smooth_gradL_composition_term
 s% ctrl% clip_D_limit = clip_D_limit
 s% ctrl% fix_eps_grav_transition_to_grid = fix_eps_grav_transition_to_grid

s% ctrl% okay_to_reduce_gradT_excess = okay_to_reduce_gradT_excess
s% ctrl% gradT_excess_f1 = gradT_excess_f1
s% ctrl% gradT_excess_f2 = gradT_excess_f2
s% ctrl% gradT_excess_age_fraction = gradT_excess_age_fraction
s% ctrl% gradT_excess_max_change = gradT_excess_max_change
s% ctrl% gradT_excess_lambda1 = gradT_excess_lambda1
s% ctrl% gradT_excess_beta1 = gradT_excess_beta1
s% ctrl% gradT_excess_lambda2 = gradT_excess_lambda2
s% ctrl% gradT_excess_beta2 = gradT_excess_beta2
s% ctrl% gradT_excess_dlambda = gradT_excess_dlambda
s% ctrl% gradT_excess_dbeta = gradT_excess_dbeta
s% ctrl% gradT_excess_max_center_h1 = gradT_excess_max_center_h1
s% ctrl% gradT_excess_min_center_he4 = gradT_excess_min_center_he4
s% ctrl% gradT_excess_max_logT = gradT_excess_max_logT
s% ctrl% gradT_excess_min_log_tau_full_on = gradT_excess_min_log_tau_full_on
s% ctrl% gradT_excess_max_log_tau_full_off = gradT_excess_max_log_tau_full_off

 s% ctrl% D_mix_zero_region_bottom_q = D_mix_zero_region_bottom_q
 s% ctrl% D_mix_zero_region_top_q = D_mix_zero_region_top_q
 s% ctrl% dq_D_mix_zero_at_H_He_crossover = dq_D_mix_zero_at_H_He_crossover
 s% ctrl% dq_D_mix_zero_at_H_C_crossover = dq_D_mix_zero_at_H_C_crossover

 s% ctrl% use_superad_reduction = use_superad_reduction
 s% ctrl% superad_reduction_gamma_limit = superad_reduction_gamma_limit
 s% ctrl% superad_reduction_gamma_limit_scale = superad_reduction_gamma_limit_scale
 s% ctrl% superad_reduction_gamma_inv_scale = superad_reduction_gamma_inv_scale
 s% ctrl% superad_reduction_diff_grads_limit = superad_reduction_diff_grads_limit
 s% ctrl% superad_reduction_limit = superad_reduction_limit
 
 s% ctrl% max_logT_for_mlt = max_logT_for_mlt
 s% ctrl% mlt_make_surface_no_mixing = mlt_make_surface_no_mixing
 s% ctrl% do_normalize_dqs_as_part_of_set_qs = do_normalize_dqs_as_part_of_set_qs

 s% ctrl% thermohaline_coeff = thermohaline_coeff
 s% ctrl% thermohaline_option = thermohaline_option
 s% ctrl% mixing_length_alpha = mixing_length_alpha
 s% ctrl% remove_small_D_limit = remove_small_D_limit
 s% ctrl% alt_scale_height_flag = alt_scale_height_flag
 s% ctrl% Henyey_MLT_y_param = Henyey_MLT_y_param
 s% ctrl% Henyey_MLT_nu_param = Henyey_MLT_nu_param
 s% ctrl% make_gradr_sticky_in_solver_iters = make_gradr_sticky_in_solver_iters
 s% ctrl% min_logT_for_make_gradr_sticky_in_solver_iters = min_logT_for_make_gradr_sticky_in_solver_iters
 s% ctrl% no_MLT_below_shock = no_MLT_below_shock
 s% ctrl% MLT_option = MLT_option
 s% ctrl% steps_before_use_TDC = steps_before_use_TDC
 s% ctrl% mlt_use_rotation_correction = mlt_use_rotation_correction
 s% ctrl% mlt_Pturb_factor = mlt_Pturb_factor

 s% ctrl% burn_z_mix_region_logT = burn_z_mix_region_logT
 s% ctrl% burn_he_mix_region_logT = burn_he_mix_region_logT
 s% ctrl% burn_h_mix_region_logT = burn_h_mix_region_logT
 s% ctrl% max_Y_for_burn_z_mix_region = max_Y_for_burn_z_mix_region
 s% ctrl% max_X_for_burn_he_mix_region = max_X_for_burn_he_mix_region
 
 s% ctrl% limit_overshoot_Hp_using_size_of_convection_zone = limit_overshoot_Hp_using_size_of_convection_zone

 s% ctrl% predictive_mix = predictive_mix
 s% ctrl% predictive_superad_thresh = predictive_superad_thresh
 s% ctrl% predictive_avoid_reversal = predictive_avoid_reversal
 s% ctrl% predictive_limit_ingestion = predictive_limit_ingestion
 s% ctrl% predictive_ingestion_factor = predictive_ingestion_factor
 s% ctrl% predictive_zone_type = predictive_zone_type
 s% ctrl% predictive_zone_loc = predictive_zone_loc
 s% ctrl% predictive_bdy_loc = predictive_bdy_loc
 s% ctrl% predictive_bdy_q_min = predictive_bdy_q_min
 s% ctrl% predictive_bdy_q_max = predictive_bdy_q_max

 s% ctrl% do_conv_premix = do_conv_premix
 s% ctrl% conv_premix_avoid_increase = conv_premix_avoid_increase
 s% ctrl% conv_premix_time_factor = conv_premix_time_factor
 s% ctrl% conv_premix_fix_pgas = conv_premix_fix_pgas
 s% ctrl% conv_premix_dump_snapshots = conv_premix_dump_snapshots
 s% ctrl% do_premix_heating = do_premix_heating

 s% ctrl% overshoot_f = overshoot_f
 s% ctrl% overshoot_f0 = overshoot_f0
 s% ctrl% overshoot_D0 = overshoot_D0
 s% ctrl% overshoot_Delta0 = overshoot_Delta0
 s% ctrl% overshoot_mass_full_on = overshoot_mass_full_on
 s% ctrl% overshoot_mass_full_off = overshoot_mass_full_off
 s% ctrl% overshoot_scheme = overshoot_scheme
 s% ctrl% overshoot_zone_type = overshoot_zone_type
 s% ctrl% overshoot_zone_loc = overshoot_zone_loc
 s% ctrl% overshoot_bdy_loc = overshoot_bdy_loc
 s% ctrl% overshoot_D_min = overshoot_D_min
 s% ctrl% overshoot_brunt_B_max = overshoot_brunt_B_max

 s% ctrl% max_conv_vel_div_csound = max_conv_vel_div_csound
 s% ctrl% max_v_for_convection = max_v_for_convection
 s% ctrl% max_q_for_convection_with_hydro_on = max_q_for_convection_with_hydro_on
 s% ctrl% max_v_div_cs_for_convection = max_v_div_cs_for_convection
 s% ctrl% max_abs_du_div_cs_for_convection = max_abs_du_div_cs_for_convection

 s% ctrl% calculate_Brunt_B = calculate_Brunt_B
 s% ctrl% calculate_Brunt_N2 = calculate_Brunt_N2
 s% ctrl% brunt_N2_coefficient = brunt_N2_coefficient
 s% ctrl% num_cells_for_smooth_brunt_B = num_cells_for_smooth_brunt_B
 s% ctrl% threshold_for_smooth_brunt_B = threshold_for_smooth_brunt_B
 s% ctrl% min_magnitude_brunt_B = min_magnitude_brunt_B

 s% ctrl% min_overshoot_q = min_overshoot_q
 s% ctrl% overshoot_alpha = overshoot_alpha

   s% ctrl% RSP_max_num_periods = RSP_max_num_periods
   s% ctrl% RSP_target_steps_per_cycle = RSP_target_steps_per_cycle
   s% ctrl% RSP_min_max_R_for_periods = RSP_min_max_R_for_periods
   s% ctrl% RSP_min_deltaR_for_periods = RSP_min_deltaR_for_periods
   s% ctrl% RSP_default_PERIODLIN = RSP_default_PERIODLIN
   s% ctrl% RSP_min_PERIOD_div_PERIODLIN = RSP_min_PERIOD_div_PERIODLIN
   s% ctrl% RSP_GREKM_avg_abs_frac_new = RSP_GREKM_avg_abs_frac_new
   s% ctrl% RSP_GREKM_avg_abs_limit = RSP_GREKM_avg_abs_limit
   s% ctrl% RSP_theta = RSP_theta
   s% ctrl% RSP_thetat = RSP_thetat
   s% ctrl% RSP_thetau = RSP_thetau
   s% ctrl% RSP_thetae = RSP_thetae
   s% ctrl% RSP_thetaq = RSP_thetaq
   s% ctrl% RSP_wtr = RSP_wtr
   s% ctrl% RSP_wtc = RSP_wtc
   s% ctrl% RSP_wtt = RSP_wtt
   s% ctrl% RSP_gam = RSP_gam
   s% ctrl% RSP_alfa = RSP_alfa
   s% ctrl% RSP_alfap = RSP_alfap
   s% ctrl% RSP_alfam = RSP_alfam
   s% ctrl% RSP_alfat = RSP_alfat
   s% ctrl% RSP_alfas = RSP_alfas
   s% ctrl% RSP_alfac = RSP_alfac
   s% ctrl% RSP_alfad = RSP_alfad
   s% ctrl% RSP_gammar = RSP_gammar
   s% ctrl% RSP_efl0 = RSP_efl0
   s% ctrl% RSP_min_tau_for_turbulent_flux = RSP_min_tau_for_turbulent_flux
   s% ctrl% RSP_cq = RSP_cq
   s% ctrl% RSP_zsh = RSP_zsh
   s% ctrl% RSP_Qvisc_quadratic = RSP_Qvisc_quadratic
   s% ctrl% RSP_Qvisc_linear = RSP_Qvisc_linear
   s% ctrl% RSP_Qvisc_linear_static = RSP_Qvisc_linear_static
   s% ctrl% RSP_tol_max_corr = RSP_tol_max_corr
   s% ctrl% RSP_tol_max_resid = RSP_tol_max_resid
   s% ctrl% RSP_max_iters_per_try = RSP_max_iters_per_try
   s% ctrl% RSP_max_retries_per_step = RSP_max_retries_per_step
   s% ctrl% RSP_nz_div_IBOTOM = RSP_nz_div_IBOTOM
   s% ctrl% RSP_kick_vsurf_km_per_sec = RSP_kick_vsurf_km_per_sec
   s% ctrl% RSP_fraction_1st_overtone = RSP_fraction_1st_overtone
   s% ctrl% RSP_fraction_2nd_overtone = RSP_fraction_2nd_overtone
   s% ctrl% RSP_Avel = RSP_Avel
   s% ctrl% RSP_Arnd = RSP_Arnd
   s% ctrl% RSP_mode_for_setting_PERIODLIN = RSP_mode_for_setting_PERIODLIN
   s% ctrl% RSP_initial_dt_factor = RSP_initial_dt_factor
   s% ctrl% RSP_v_div_cs_threshold_for_dt_limit = RSP_v_div_cs_threshold_for_dt_limit
   s% ctrl% RSP_max_dt_times_min_dr_div_cs = RSP_max_dt_times_min_dr_div_cs
   s% ctrl% RSP_max_dt_times_min_rad_diff_time = RSP_max_dt_times_min_rad_diff_time
   s% ctrl% RSP_max_dt = RSP_max_dt
   s% ctrl% RSP_testing = RSP_testing
   s% ctrl% RSP_report_limit_dt = RSP_report_limit_dt
   s% ctrl% RSP_use_Prad_for_Psurf = RSP_use_Prad_for_Psurf
   s% ctrl% RSP_report_undercorrections = RSP_report_undercorrections
   s% ctrl% RSP_use_atm_grey_with_kap_for_Psurf = RSP_use_atm_grey_with_kap_for_Psurf
   s% ctrl% use_other_RSP_linear_analysis = use_other_RSP_linear_analysis
   s% ctrl% use_other_RSP_build_model = use_other_RSP_build_model
   s% ctrl% RSP_kap_density_factor = RSP_kap_density_factor
   s% ctrl% RSP_fixed_Psurf = RSP_fixed_Psurf
   s% ctrl% RSP_hydro_only = RSP_hydro_only
   s% ctrl% RSP_tau_surf_for_atm_grey_with_kap = RSP_tau_surf_for_atm_grey_with_kap
   s% ctrl% RSP_Psurf = RSP_Psurf
   s% ctrl% set_RSP_Psurf_to_multiple_of_initial_P1 = set_RSP_Psurf_to_multiple_of_initial_P1
   s% ctrl% RSP_surface_tau = RSP_surface_tau
   s% ctrl% RSP_write_map = RSP_write_map
   s% ctrl% RSP_trace_RSP_build_model = RSP_trace_RSP_build_model
   s% ctrl% RSP_map_filename = RSP_map_filename
   s% ctrl% RSP_map_columns_filename = RSP_map_columns_filename
   s% ctrl% RSP_map_history_filename = RSP_map_history_filename
   s% ctrl% RSP_map_first_period = RSP_map_first_period
   s% ctrl% RSP_map_last_period = RSP_map_last_period
   s% ctrl% RSP_map_zone_interval = RSP_map_zone_interval
   s% ctrl% RSP_nmodes = RSP_nmodes
   s% ctrl% RSP_work_period = RSP_work_period
   s% ctrl% RSP_work_filename = RSP_work_filename
   s% ctrl% RSP_nz_outer = RSP_nz_outer
   s% ctrl% RSP_max_outer_dm_tries = RSP_max_outer_dm_tries
   s% ctrl% RSP_max_inner_scale_tries = RSP_max_inner_scale_tries
   s% ctrl% RSP_relax_max_tries = RSP_relax_max_tries
   s% ctrl% RSP_T_anchor_tolerance = RSP_T_anchor_tolerance
   s% ctrl% RSP_T_inner_tolerance = RSP_T_inner_tolerance
   s% ctrl% RSP_relax_dm_tolerance = RSP_relax_dm_tolerance
   s% ctrl% RSP_dq_1_factor = RSP_dq_1_factor
   s% ctrl% use_RSP_new_start_scheme = use_RSP_new_start_scheme
   s% ctrl% RSP_do_check_omega = RSP_do_check_omega
   s% ctrl% RSP_report_check_omega_changes = RSP_report_check_omega_changes
   s% ctrl% RSP_nz = RSP_nz
   s% ctrl% RSP_T_anchor = RSP_T_anchor
   s% ctrl% RSP_T_inner = RSP_T_inner
   s% ctrl% RSP_relax_initial_model = RSP_relax_initial_model
   s% ctrl% RSP_relax_alfap_before_alfat = RSP_relax_alfap_before_alfat
   s% ctrl% RSP_relax_adjust_inner_mass_distribution = RSP_relax_adjust_inner_mass_distribution
   s% ctrl% RSP_Teff = RSP_Teff
   s% ctrl% RSP_mass = RSP_mass
   s% ctrl% RSP_L = RSP_L
   s% ctrl% RSP_X = RSP_X
   s% ctrl% RSP_Z = RSP_Z

 s% ctrl% RTI_smooth_mass = RTI_smooth_mass
 s% ctrl% RTI_smooth_iterations = RTI_smooth_iterations
 s% ctrl% RTI_smooth_fraction = RTI_smooth_fraction

 s% ctrl% alpha_RTI_diffusion_factor = alpha_RTI_diffusion_factor
 s% ctrl% dudt_RTI_diffusion_factor = dudt_RTI_diffusion_factor
 s% ctrl% dedt_RTI_diffusion_factor = dedt_RTI_diffusion_factor
 s% ctrl% dlnddt_RTI_diffusion_factor = dlnddt_RTI_diffusion_factor
 s% ctrl% composition_RTI_diffusion_factor = composition_RTI_diffusion_factor
 s% ctrl% max_M_RTI_factors_full_on = max_M_RTI_factors_full_on
 s% ctrl% min_M_RTI_factors_full_off = min_M_RTI_factors_full_off

 s% ctrl% alpha_RTI_src_min_v_div_cs = alpha_RTI_src_min_v_div_cs
 s% ctrl% alpha_RTI_src_max_q = alpha_RTI_src_max_q
 s% ctrl% alpha_RTI_src_min_q = alpha_RTI_src_min_q

 s% ctrl% T_mix_limit = T_mix_limit
 s% ctrl% mlt_gradT_fraction = mlt_gradT_fraction

 ! atmosphere -- surface boundary conditions
 s% ctrl% atm_option = atm_option
 s% ctrl% atm_off_table_option = atm_off_table_option
 s% ctrl% Pextra_factor = Pextra_factor
 s% ctrl% atm_fixed_Teff = atm_fixed_Teff
 s% ctrl% atm_fixed_Psurf = atm_fixed_Psurf
 s% ctrl% atm_fixed_Tsurf = atm_fixed_Tsurf

 s% ctrl% atm_T_tau_relation = atm_T_tau_relation
 s% ctrl% atm_T_tau_opacity = atm_T_tau_opacity
 s% ctrl% atm_T_tau_errtol = atm_T_tau_errtol
 s% ctrl% atm_T_tau_max_iters = atm_T_tau_max_iters
 s% ctrl% atm_T_tau_max_steps = atm_T_tau_max_steps

 s% ctrl% atm_table = atm_table

 s% ctrl% atm_irradiated_opacity = atm_irradiated_opacity
 s% ctrl% atm_irradiated_errtol = atm_irradiated_errtol
 s% ctrl% atm_irradiated_T_eq = atm_irradiated_T_eq
 s% ctrl% atm_irradiated_kap_v = atm_irradiated_kap_v
 s% ctrl% atm_irradiated_kap_v_div_kap_th = atm_irradiated_kap_v_div_kap_th
 s% ctrl% atm_irradiated_P_surf = atm_irradiated_P_surf
 s% ctrl% atm_irradiated_max_iters = atm_irradiated_max_iters

 s% ctrl% use_compression_outer_BC = use_compression_outer_BC
 s% ctrl% use_momentum_outer_BC = use_momentum_outer_BC
 s% ctrl% Tsurf_factor = Tsurf_factor
 s% ctrl% use_zero_Pgas_outer_BC = use_zero_Pgas_outer_BC
 s% ctrl% fixed_vsurf = fixed_vsurf
 s% ctrl% use_fixed_vsurf_outer_BC = use_fixed_vsurf_outer_BC
 s% ctrl% fixed_Psurf = fixed_Psurf
 s% ctrl% use_fixed_Psurf_outer_BC = use_fixed_Psurf_outer_BC

 s% ctrl% atm_build_tau_outer = atm_build_tau_outer
 s% ctrl% atm_build_dlogtau = atm_build_dlogtau
 s% ctrl% atm_build_errtol = atm_build_errtol

 s% ctrl% use_T_tau_gradr_factor = use_T_tau_gradr_factor

 ! extra heat near surface to model irradiation
 s% ctrl% irradiation_flux = irradiation_flux
 s% ctrl% column_depth_for_irradiation = column_depth_for_irradiation

 ! extra heat
 s% ctrl% inject_uniform_extra_heat = inject_uniform_extra_heat
 s% ctrl% min_q_for_uniform_extra_heat = min_q_for_uniform_extra_heat
 s% ctrl% max_q_for_uniform_extra_heat = max_q_for_uniform_extra_heat
 s% ctrl% inject_extra_ergs_sec = inject_extra_ergs_sec
 s% ctrl% base_of_inject_extra_ergs_sec = base_of_inject_extra_ergs_sec
 s% ctrl% total_mass_for_inject_extra_ergs_sec = total_mass_for_inject_extra_ergs_sec
 s% ctrl% start_time_for_inject_extra_ergs_sec = start_time_for_inject_extra_ergs_sec
 s% ctrl% duration_for_inject_extra_ergs_sec = duration_for_inject_extra_ergs_sec
 s% ctrl% inject_until_reach_model_with_total_energy = inject_until_reach_model_with_total_energy

 ! mass gain or loss
 s% ctrl% mass_change = mass_change
 s% ctrl% mass_change_full_off_dt = mass_change_full_off_dt
 s% ctrl% mass_change_full_on_dt = mass_change_full_on_dt
 s% ctrl% trace_dt_control_mass_change = trace_dt_control_mass_change
 s% ctrl% no_wind_if_no_rotation = no_wind_if_no_rotation

 s% ctrl% min_wind = min_wind
 s% ctrl% max_wind = max_wind
 s% ctrl% use_accreted_material_j = use_accreted_material_j
 s% ctrl% accreted_material_j = accreted_material_j
 s% ctrl% D_omega_mixing_rate = D_omega_mixing_rate
 s% ctrl% D_omega_mixing_across_convection_boundary = D_omega_mixing_across_convection_boundary
 s% ctrl% max_q_for_D_omega_zero_in_convection_region = max_q_for_D_omega_zero_in_convection_region
 s% ctrl% nu_omega_mixing_rate = nu_omega_mixing_rate
 s% ctrl% nu_omega_mixing_across_convection_boundary = nu_omega_mixing_across_convection_boundary
 s% ctrl% max_q_for_nu_omega_zero_in_convection_region = max_q_for_nu_omega_zero_in_convection_region

 s% ctrl% mdot_omega_power = mdot_omega_power
 s% ctrl% max_rotational_mdot_boost = max_rotational_mdot_boost
 s% ctrl% max_mdot_jump_for_rotation = max_mdot_jump_for_rotation
 s% ctrl% lim_trace_rotational_mdot_boost = lim_trace_rotational_mdot_boost
 s% ctrl% rotational_mdot_boost_fac = rotational_mdot_boost_fac
 s% ctrl% rotational_mdot_kh_fac = rotational_mdot_kh_fac
 s% ctrl% surf_avg_tau = surf_avg_tau
 s% ctrl% surf_avg_tau_min = surf_avg_tau_min

 s% ctrl% super_eddington_scaling_factor = super_eddington_scaling_factor
 s% ctrl% super_eddington_wind_Ledd_factor = super_eddington_wind_Ledd_factor
 s% ctrl% wind_boost_full_off_L_div_Ledd = wind_boost_full_off_L_div_Ledd
 s% ctrl% wind_boost_full_on_L_div_Ledd = wind_boost_full_on_L_div_Ledd
 s% ctrl% super_eddington_wind_max_boost = super_eddington_wind_max_boost
 s% ctrl% trace_super_eddington_wind_boost = trace_super_eddington_wind_boost
 
 s% ctrl% max_tries_for_implicit_wind = max_tries_for_implicit_wind
 s% ctrl% iwind_tolerance = iwind_tolerance
 s% ctrl% iwind_lambda = iwind_lambda

 s% ctrl% cool_wind_full_on_T = cool_wind_full_on_T
 s% ctrl% hot_wind_full_on_T = hot_wind_full_on_T

 s% ctrl% rlo_scaling_factor = rlo_scaling_factor
 s% ctrl% rlo_wind_min_L = rlo_wind_min_L
 s% ctrl% rlo_wind_max_Teff = rlo_wind_max_Teff
 s% ctrl% rlo_wind_roche_lobe_radius = rlo_wind_roche_lobe_radius
 s% ctrl% roche_lobe_xfer_full_on = roche_lobe_xfer_full_on
 s% ctrl% roche_lobe_xfer_full_off = roche_lobe_xfer_full_off
 s% ctrl% rlo_wind_base_mdot = rlo_wind_base_mdot
 s% ctrl% rlo_wind_scale_height = rlo_wind_scale_height

 s% ctrl% hot_wind_scheme = hot_wind_scheme
 s% ctrl% cool_wind_RGB_scheme = cool_wind_RGB_scheme
 s% ctrl% cool_wind_AGB_scheme = cool_wind_AGB_scheme
 s% ctrl% RGB_to_AGB_wind_switch = RGB_to_AGB_wind_switch
 s% ctrl% Reimers_scaling_factor = Reimers_scaling_factor
 s% ctrl% Blocker_scaling_factor = Blocker_scaling_factor
 s% ctrl% de_Jager_scaling_factor = de_Jager_scaling_factor
 s% ctrl% van_Loon_scaling_factor = van_Loon_scaling_factor
 s% ctrl% Nieuwenhuijzen_scaling_factor = Nieuwenhuijzen_scaling_factor
 s% ctrl% Vink_scaling_factor = Vink_scaling_factor
 s% ctrl% Dutch_scaling_factor = Dutch_scaling_factor
 s% ctrl% Dutch_wind_lowT_scheme = Dutch_wind_lowT_scheme

 s% ctrl% wind_H_envelope_limit = wind_H_envelope_limit
 s% ctrl% wind_H_He_envelope_limit = wind_H_He_envelope_limit
 s% ctrl% wind_He_layer_limit = wind_He_layer_limit

 s% ctrl% max_logT_for_k_below_const_q = max_logT_for_k_below_const_q
 s% ctrl% max_q_for_k_below_const_q = max_q_for_k_below_const_q
 s% ctrl% min_q_for_k_below_const_q = min_q_for_k_below_const_q
 s% ctrl% max_logT_for_k_const_mass = max_logT_for_k_const_mass
 s% ctrl% min_q_for_k_const_mass = min_q_for_k_const_mass
 s% ctrl% max_q_for_k_const_mass = max_q_for_k_const_mass

 ! composition of added mass
 s% ctrl% accrete_same_as_surface = accrete_same_as_surface

 s% ctrl% accrete_given_mass_fractions = accrete_given_mass_fractions
 s% ctrl% num_accretion_species = num_accretion_species
 s% ctrl% accretion_species_id = accretion_species_id
 s% ctrl% accretion_species_xa = accretion_species_xa

 s% ctrl% accretion_h1 = accretion_h1
 s% ctrl% accretion_h2 = accretion_h2
 s% ctrl% accretion_he3 = accretion_he3
 s% ctrl% accretion_he4 = accretion_he4
 s% ctrl% accretion_zfracs = accretion_zfracs
 s% ctrl% accretion_dump_missing_metals_into_heaviest = accretion_dump_missing_metals_into_heaviest

 ! special list of z fractions
 s% ctrl% z_fraction_li = z_fraction_li
 s% ctrl% z_fraction_be = z_fraction_be
 s% ctrl% z_fraction_b = z_fraction_b
 s% ctrl% z_fraction_c = z_fraction_c
 s% ctrl% z_fraction_n = z_fraction_n
 s% ctrl% z_fraction_o = z_fraction_o
 s% ctrl% z_fraction_f = z_fraction_f
 s% ctrl% z_fraction_ne = z_fraction_ne
 s% ctrl% z_fraction_na = z_fraction_na
 s% ctrl% z_fraction_mg = z_fraction_mg
 s% ctrl% z_fraction_al = z_fraction_al
 s% ctrl% z_fraction_si = z_fraction_si
 s% ctrl% z_fraction_p = z_fraction_p
 s% ctrl% z_fraction_s = z_fraction_s
 s% ctrl% z_fraction_cl = z_fraction_cl
 s% ctrl% z_fraction_ar = z_fraction_ar
 s% ctrl% z_fraction_k = z_fraction_k
 s% ctrl% z_fraction_ca = z_fraction_ca
 s% ctrl% z_fraction_sc = z_fraction_sc
 s% ctrl% z_fraction_ti = z_fraction_ti
 s% ctrl% z_fraction_v = z_fraction_v
 s% ctrl% z_fraction_cr = z_fraction_cr
 s% ctrl% z_fraction_mn = z_fraction_mn
 s% ctrl% z_fraction_fe = z_fraction_fe
 s% ctrl% z_fraction_co = z_fraction_co
 s% ctrl% z_fraction_ni = z_fraction_ni
 s% ctrl% z_fraction_cu = z_fraction_cu
 s% ctrl% z_fraction_zn = z_fraction_zn

 s% ctrl% lgT_lo_for_set_new_abundances = lgT_lo_for_set_new_abundances
 s% ctrl% lgT_hi_for_set_new_abundances = lgT_hi_for_set_new_abundances

 ! automatic stops for mass loss/gain
 s% ctrl% max_star_mass_for_gain = max_star_mass_for_gain
 s% ctrl% min_star_mass_for_loss = min_star_mass_for_loss
 s% ctrl% max_T_center_for_any_mass_loss = max_T_center_for_any_mass_loss
 s% ctrl% max_T_center_for_full_mass_loss = max_T_center_for_full_mass_loss

 ! extra power source
 s% ctrl% extra_power_source = extra_power_source

 ! relaxation parameters
 s% ctrl% relax_dlnZ = relax_dlnZ
 s% ctrl% relax_dY = relax_dY

 ! mesh adjustment
 s% ctrl% show_mesh_changes = show_mesh_changes
 s% ctrl% okay_to_remesh = okay_to_remesh
 s% ctrl% restore_mesh_on_retry = restore_mesh_on_retry
 s% ctrl% num_steps_to_hold_mesh_after_retry = num_steps_to_hold_mesh_after_retry
 s% ctrl% trace_mesh_adjust_error_in_conservation = trace_mesh_adjust_error_in_conservation
 s% ctrl% max_rel_delta_IE_for_mesh_total_energy_balance = max_rel_delta_IE_for_mesh_total_energy_balance
 s% ctrl% max_allowed_nz = max_allowed_nz
 s% ctrl% mesh_max_allowed_ratio = mesh_max_allowed_ratio
 s% ctrl% remesh_max_allowed_logT = remesh_max_allowed_logT
 s% ctrl% max_delta_x_for_merge = max_delta_x_for_merge

 s% ctrl% mesh_ok_to_merge = mesh_ok_to_merge
 s% ctrl% mesh_max_k_old_for_split = mesh_max_k_old_for_split
 s% ctrl% mesh_min_k_old_for_split = mesh_min_k_old_for_split
 s% ctrl% mesh_adjust_get_T_from_E = mesh_adjust_get_T_from_E

 s% ctrl% max_dq = max_dq
 s% ctrl% min_dq = min_dq
 s% ctrl% min_dq_for_split = min_dq_for_split
 s% ctrl% min_dq_for_xa = min_dq_for_xa
 s% ctrl% min_dq_for_xa_convective = min_dq_for_xa_convective
 s% ctrl% min_dq_for_logT = min_dq_for_logT

 s% ctrl% mesh_min_dlnR = mesh_min_dlnR
 s% ctrl% merge_if_dlnR_too_small = merge_if_dlnR_too_small

 s% ctrl% mesh_min_dr_div_dRstar = mesh_min_dr_div_dRstar
 s% ctrl% merge_if_dr_div_dRstar_too_small = merge_if_dr_div_dRstar_too_small

 s% ctrl% mesh_min_dr_div_cs = mesh_min_dr_div_cs
 s% ctrl% merge_if_dr_div_cs_too_small = merge_if_dr_div_cs_too_small

 s% ctrl% max_center_cell_dq = max_center_cell_dq
 s% ctrl% max_surface_cell_dq = max_surface_cell_dq
 s% ctrl% max_num_subcells = max_num_subcells
 s% ctrl% max_num_merge_cells = max_num_merge_cells

 s% ctrl% mesh_delta_coeff = mesh_delta_coeff
 s% ctrl% mesh_delta_coeff_for_highT = mesh_delta_coeff_for_highT
 s% ctrl% logT_max_for_standard_mesh_delta_coeff = logT_max_for_standard_mesh_delta_coeff
 s% ctrl% logT_min_for_highT_mesh_delta_coeff = logT_min_for_highT_mesh_delta_coeff
 s% ctrl% mesh_Pgas_div_P_exponent = mesh_Pgas_div_P_exponent

 s% ctrl% remesh_dt_limit = remesh_dt_limit

 s% ctrl% E_function_weight = E_function_weight
 s% ctrl% E_function_param = E_function_param
 s% ctrl% P_function_weight = P_function_weight

 s% ctrl% mesh_logX_species = mesh_logX_species
 s% ctrl% mesh_logX_min_for_extra = mesh_logX_min_for_extra
 s% ctrl% mesh_dlogX_dlogP_extra = mesh_dlogX_dlogP_extra
 s% ctrl% mesh_dlogX_dlogP_full_on = mesh_dlogX_dlogP_full_on
 s% ctrl% mesh_dlogX_dlogP_full_off = mesh_dlogX_dlogP_full_off

 s% ctrl% mesh_dlog_eps_min_for_extra = mesh_dlog_eps_min_for_extra
 s% ctrl% mesh_dlog_eps_dlogP_full_on = mesh_dlog_eps_dlogP_full_on
 s% ctrl% mesh_dlog_eps_dlogP_full_off = mesh_dlog_eps_dlogP_full_off

 s% ctrl% mesh_dlog_pp_dlogP_extra = mesh_dlog_pp_dlogP_extra
 s% ctrl% mesh_dlog_cno_dlogP_extra = mesh_dlog_cno_dlogP_extra
 s% ctrl% mesh_dlog_3alf_dlogP_extra = mesh_dlog_3alf_dlogP_extra

 s% ctrl% mesh_dlog_burn_c_dlogP_extra = mesh_dlog_burn_c_dlogP_extra
 s% ctrl% mesh_dlog_burn_n_dlogP_extra = mesh_dlog_burn_n_dlogP_extra
 s% ctrl% mesh_dlog_burn_o_dlogP_extra = mesh_dlog_burn_o_dlogP_extra
 s% ctrl% mesh_dlog_burn_ne_dlogP_extra = mesh_dlog_burn_ne_dlogP_extra
 s% ctrl% mesh_dlog_burn_na_dlogP_extra = mesh_dlog_burn_na_dlogP_extra
 s% ctrl% mesh_dlog_burn_mg_dlogP_extra = mesh_dlog_burn_mg_dlogP_extra
 s% ctrl% mesh_dlog_burn_si_dlogP_extra = mesh_dlog_burn_si_dlogP_extra
 s% ctrl% mesh_dlog_burn_s_dlogP_extra = mesh_dlog_burn_s_dlogP_extra
 s% ctrl% mesh_dlog_burn_ar_dlogP_extra = mesh_dlog_burn_ar_dlogP_extra
 s% ctrl% mesh_dlog_burn_ca_dlogP_extra = mesh_dlog_burn_ca_dlogP_extra
 s% ctrl% mesh_dlog_burn_ti_dlogP_extra = mesh_dlog_burn_ti_dlogP_extra
 s% ctrl% mesh_dlog_burn_cr_dlogP_extra = mesh_dlog_burn_cr_dlogP_extra
 s% ctrl% mesh_dlog_burn_fe_dlogP_extra = mesh_dlog_burn_fe_dlogP_extra

 s% ctrl% mesh_dlog_cc_dlogP_extra = mesh_dlog_cc_dlogP_extra
 s% ctrl% mesh_dlog_co_dlogP_extra = mesh_dlog_co_dlogP_extra
 s% ctrl% mesh_dlog_oo_dlogP_extra = mesh_dlog_oo_dlogP_extra

 s% ctrl% mesh_dlog_pnhe4_dlogP_extra = mesh_dlog_pnhe4_dlogP_extra
 s% ctrl% mesh_dlog_photo_dlogP_extra = mesh_dlog_photo_dlogP_extra
 s% ctrl% mesh_dlog_other_dlogP_extra = mesh_dlog_other_dlogP_extra
 
 s% ctrl% mesh_delta_coeff_factor_smooth_iters = mesh_delta_coeff_factor_smooth_iters

 s% ctrl% T_function1_weight = T_function1_weight
 s% ctrl% T_function2_weight = T_function2_weight
 s% ctrl% T_function2_param = T_function2_param

 s% ctrl% R_function_weight = R_function_weight
 s% ctrl% R_function_param = R_function_param

 s% ctrl% R_function2_weight = R_function2_weight
 s% ctrl% R_function2_param1 = R_function2_param1
 s% ctrl% R_function2_param2 = R_function2_param2

 s% ctrl% R_function3_weight = R_function3_weight

 s% ctrl% M_function_weight = M_function_weight
 s% ctrl% M_function_param = M_function_param

 s% ctrl% gradT_function_weight = gradT_function_weight
 s% ctrl% log_tau_function_weight = log_tau_function_weight
 s% ctrl% log_kap_function_weight = log_kap_function_weight
 s% ctrl% omega_function_weight = omega_function_weight

 s% ctrl% gam_function_weight = gam_function_weight
 s% ctrl% gam_function_param1 = gam_function_param1
 s% ctrl% gam_function_param2 = gam_function_param2

 s% ctrl% xa_function_species = xa_function_species
 s% ctrl% xa_function_weight = xa_function_weight
 s% ctrl% xa_function_param = xa_function_param
 s% ctrl% xa_mesh_delta_coeff = xa_mesh_delta_coeff
 
 s% ctrl% use_split_merge_amr = use_split_merge_amr
 s% ctrl% split_merge_amr_nz_baseline = split_merge_amr_nz_baseline
 s% ctrl% split_merge_amr_nz_r_core = split_merge_amr_nz_r_core
 s% ctrl% split_merge_amr_nz_r_core_fraction = split_merge_amr_nz_r_core_fraction
 s% ctrl% split_merge_amr_mesh_delta_coeff = split_merge_amr_mesh_delta_coeff
 s% ctrl% split_merge_amr_log_zoning = split_merge_amr_log_zoning
 s% ctrl% split_merge_amr_hybrid_zoning = split_merge_amr_hybrid_zoning
 s% ctrl% split_merge_amr_flipped_hybrid_zoning = split_merge_amr_flipped_hybrid_zoning
 s% ctrl% split_merge_amr_logtau_zoning = split_merge_amr_logtau_zoning
 s% ctrl% split_merge_amr_okay_to_split_nz = split_merge_amr_okay_to_split_nz
 s% ctrl% split_merge_amr_okay_to_split_1 = split_merge_amr_okay_to_split_1
 s% ctrl% merge_amr_inhibit_at_jumps = merge_amr_inhibit_at_jumps
 s% ctrl% split_merge_amr_MaxLong = split_merge_amr_MaxLong
 s% ctrl% split_merge_amr_MaxShort = split_merge_amr_MaxShort
 s% ctrl% merge_amr_max_abs_du_div_cs = merge_amr_max_abs_du_div_cs
 s% ctrl% merge_amr_ignore_surface_cells = merge_amr_ignore_surface_cells
 s% ctrl% merge_amr_du_div_cs_limit_only_for_compression = merge_amr_du_div_cs_limit_only_for_compression
 s% ctrl% split_merge_amr_avoid_repeated_remesh = split_merge_amr_avoid_repeated_remesh
 s% ctrl% merge_amr_k_for_ignore_surface_cells = merge_amr_k_for_ignore_surface_cells
 s% ctrl% split_merge_amr_dq_min = split_merge_amr_dq_min
 s% ctrl% split_merge_amr_dq_max = split_merge_amr_dq_max
 s% ctrl% split_merge_amr_r_core_cm = split_merge_amr_r_core_cm
 s% ctrl% split_merge_amr_max_iters = split_merge_amr_max_iters
 s% ctrl% trace_split_merge_amr = trace_split_merge_amr
 s% ctrl% equal_split_density_amr = equal_split_density_amr

 ! nuclear reaction parameters
 s% ctrl% screening_mode = screening_mode
 s% ctrl% default_net_name = default_net_name

 s% ctrl% net_logTcut_lo = net_logTcut_lo
 s% ctrl% net_logTcut_lim = net_logTcut_lim

 s% ctrl% eps_nuc_factor = eps_nuc_factor
 s% ctrl% op_split_burn_eps_nuc_infall_limit = op_split_burn_eps_nuc_infall_limit
 s% ctrl% eps_WD_sedimentation_factor = eps_WD_sedimentation_factor
 s% ctrl% max_abs_eps_nuc = max_abs_eps_nuc
 s% ctrl% dxdt_nuc_factor = dxdt_nuc_factor
 s% ctrl% max_abar_for_burning = max_abar_for_burning
 s% ctrl% fe56ec_fake_factor = fe56ec_fake_factor
 s% ctrl% min_T_for_fe56ec_fake_factor = min_T_for_fe56ec_fake_factor
 s% ctrl% weak_rate_factor = weak_rate_factor

 ! mixing
 s% ctrl% mix_factor = mix_factor

 s% ctrl% sig_term_limit = sig_term_limit

 s% ctrl% sig_min_factor_for_high_Tcenter = sig_min_factor_for_high_Tcenter
 s% ctrl% Tcenter_min_for_sig_min_factor_full_on = Tcenter_min_for_sig_min_factor_full_on
 s% ctrl% Tcenter_max_for_sig_min_factor_full_off = Tcenter_max_for_sig_min_factor_full_off
 s% ctrl% max_delta_m_to_bdy_for_sig_min_factor = max_delta_m_to_bdy_for_sig_min_factor
 s% ctrl% delta_m_lower_for_sig_min_factor = delta_m_lower_for_sig_min_factor
 s% ctrl% delta_m_upper_for_sig_min_factor = delta_m_upper_for_sig_min_factor

 s% ctrl% am_sig_term_limit = am_sig_term_limit
 s% ctrl% am_D_mix_factor = am_D_mix_factor
 s% ctrl% am_gradmu_factor = am_gradmu_factor
 s% ctrl% am_nu_factor = am_nu_factor

 s% ctrl% D_visc_factor = D_visc_factor
 s% ctrl% D_DSI_factor = D_DSI_factor
 s% ctrl% D_SH_factor = D_SH_factor
 s% ctrl% D_SSI_factor = D_SSI_factor
 s% ctrl% D_ES_factor = D_ES_factor
 s% ctrl% D_GSF_factor = D_GSF_factor
 s% ctrl% D_ST_factor = D_ST_factor

 s% ctrl% am_nu_non_rotation_factor = am_nu_non_rotation_factor
 s% ctrl% skip_rotation_in_convection_zones = skip_rotation_in_convection_zones
 s% ctrl% am_nu_DSI_factor = am_nu_DSI_factor
 s% ctrl% am_nu_SH_factor = am_nu_SH_factor
 s% ctrl% am_nu_SSI_factor = am_nu_SSI_factor
 s% ctrl% am_nu_ES_factor = am_nu_ES_factor
 s% ctrl% am_nu_GSF_factor = am_nu_GSF_factor
 s% ctrl% am_nu_ST_factor = am_nu_ST_factor
 s% ctrl% am_nu_visc_factor = am_nu_visc_factor

 s% ctrl% am_nu_omega_rot_factor = am_nu_omega_rot_factor
 s% ctrl% am_nu_omega_non_rot_factor = am_nu_omega_non_rot_factor
 s% ctrl% am_nu_j_rot_factor = am_nu_j_rot_factor
 s% ctrl% am_nu_j_non_rot_factor = am_nu_j_non_rot_factor

 s% ctrl% smooth_nu_ST = smooth_nu_ST
 s% ctrl% smooth_D_ST = smooth_D_ST
 s% ctrl% smooth_D_DSI = smooth_D_DSI
 s% ctrl% smooth_D_SSI = smooth_D_SSI
 s% ctrl% smooth_D_SH = smooth_D_SH
 s% ctrl% smooth_D_GSF = smooth_D_GSF
 s% ctrl% smooth_D_ES = smooth_D_ES
 s% ctrl% smooth_D_omega = smooth_D_omega
 s% ctrl% smooth_am_nu_rot = smooth_am_nu_rot
 s% ctrl% ST_angsmt = ST_angsmt
 s% ctrl% ST_angsml = ST_angsml

 s% ctrl% simple_i_rot_flag = simple_i_rot_flag
 s% ctrl% do_adjust_J_lost = do_adjust_J_lost
 s% ctrl% premix_omega = premix_omega
 s% ctrl% angular_momentum_error_warn = angular_momentum_error_warn
 s% ctrl% angular_momentum_error_retry = angular_momentum_error_retry
 s% ctrl% recalc_mixing_info_each_substep = recalc_mixing_info_each_substep
 s% ctrl% adjust_J_fraction = adjust_J_fraction
 s% ctrl% min_q_for_adjust_J_lost = min_q_for_adjust_J_lost
 s% ctrl% min_J_div_delta_J = min_J_div_delta_J
 s% ctrl% max_mdot_redo_cnt = max_mdot_redo_cnt
 s% ctrl% mdot_revise_factor = mdot_revise_factor
 s% ctrl% implicit_mdot_boost = implicit_mdot_boost
 s% ctrl% min_years_dt_for_redo_mdot = min_years_dt_for_redo_mdot
 s% ctrl% surf_omega_div_omega_crit_limit = surf_omega_div_omega_crit_limit
 s% ctrl% surf_omega_div_omega_crit_tol = surf_omega_div_omega_crit_tol
 s% ctrl% w_div_wcrit_max = w_div_wcrit_max
 s% ctrl% w_div_wcrit_max2 = w_div_wcrit_max2
 s% ctrl% fp_min = fp_min
 s% ctrl% ft_min = ft_min
 s% ctrl% fp_error_limit = fp_error_limit
 s% ctrl% ft_error_limit = ft_error_limit

 s% ctrl% D_mix_rotation_max_logT_full_on = D_mix_rotation_max_logT_full_on
 s% ctrl% D_mix_rotation_min_logT_full_off = D_mix_rotation_min_logT_full_off

 s% ctrl% set_uniform_am_nu_non_rot = set_uniform_am_nu_non_rot
 s% ctrl% uniform_am_nu_non_rot = uniform_am_nu_non_rot

 s% ctrl% set_min_am_nu_non_rot = set_min_am_nu_non_rot
 s% ctrl% min_am_nu_non_rot = min_am_nu_non_rot
 s% ctrl% min_center_Ye_for_min_am_nu_non_rot = min_center_Ye_for_min_am_nu_non_rot

 s% ctrl% set_min_D_mix = set_min_D_mix
 s% ctrl% mass_lower_limit_for_min_D_mix = mass_lower_limit_for_min_D_mix
 s% ctrl% mass_upper_limit_for_min_D_mix = mass_upper_limit_for_min_D_mix
 s% ctrl% min_D_mix = min_D_mix
 s% ctrl% set_min_D_mix_below_Tmax = set_min_D_mix_below_Tmax
 s% ctrl% min_D_mix_below_Tmax = min_D_mix_below_Tmax
 s% ctrl% set_min_D_mix_in_H_He = set_min_D_mix_in_H_He
 s% ctrl% min_D_mix_in_H_He = min_D_mix_in_H_He
 s% ctrl% min_center_Ye_for_min_D_mix = min_center_Ye_for_min_D_mix
 s% ctrl% reaction_neuQs_factor = reaction_neuQs_factor
 s% ctrl% nonlocal_NiCo_kap_gamma = nonlocal_NiCo_kap_gamma
 s% ctrl% nonlocal_NiCo_decay_heat = nonlocal_NiCo_decay_heat
 s% ctrl% dtau_gamma_NiCo_decay_heat = dtau_gamma_NiCo_decay_heat
 s% ctrl% max_logT_for_net = max_logT_for_net
 s% ctrl% smooth_outer_xa_big = smooth_outer_xa_big
 s% ctrl% smooth_outer_xa_small = smooth_outer_xa_small

 ! element diffusion parameters
 s% ctrl% diffusion_use_iben_macdonald = diffusion_use_iben_macdonald
 s% ctrl% diffusion_use_paquette = diffusion_use_paquette
 s% ctrl% diffusion_use_cgs_solver = diffusion_use_cgs_solver
 s% ctrl% diffusion_use_full_net = diffusion_use_full_net
 s% ctrl% do_WD_sedimentation_heating = do_WD_sedimentation_heating
 s% ctrl% min_xa_for_WD_sedimentation_heating = min_xa_for_WD_sedimentation_heating
 s% ctrl% do_diffusion_heating = do_diffusion_heating
 s% ctrl% do_element_diffusion = do_element_diffusion
 s% ctrl% cgs_thermal_diffusion_eta_full_on = cgs_thermal_diffusion_eta_full_on
 s% ctrl% cgs_thermal_diffusion_eta_full_off = cgs_thermal_diffusion_eta_full_off
 s% ctrl% diffusion_min_dq_at_surface = diffusion_min_dq_at_surface
 s% ctrl% diffusion_min_T_at_surface = diffusion_min_T_at_surface
 s% ctrl% diffusion_min_dq_ratio_at_surface = diffusion_min_dq_ratio_at_surface
 s% ctrl% diffusion_dt_limit = diffusion_dt_limit

 s% ctrl% diffusion_min_X_hard_limit = diffusion_min_X_hard_limit
 s% ctrl% diffusion_X_total_atol = diffusion_X_total_atol
 s% ctrl% diffusion_X_total_rtol = diffusion_X_total_rtol
 s% ctrl% diffusion_upwind_abs_v_limit = diffusion_upwind_abs_v_limit
 s% ctrl% diffusion_dt_div_timescale = diffusion_dt_div_timescale
 s% ctrl% diffusion_min_num_substeps = diffusion_min_num_substeps
 s% ctrl% diffusion_max_iters_per_substep = diffusion_max_iters_per_substep
 s% ctrl% diffusion_max_retries_per_substep = diffusion_max_retries_per_substep
 s% ctrl% diffusion_v_max = diffusion_v_max
 s% ctrl% diffusion_gamma_full_off = diffusion_gamma_full_off
 s% ctrl% diffusion_gamma_full_on = diffusion_gamma_full_on
 s% ctrl% diffusion_T_full_off = diffusion_T_full_off
 s% ctrl% D_mix_ignore_diffusion = D_mix_ignore_diffusion
 s% ctrl% diffusion_T_full_on = diffusion_T_full_on
 s% ctrl% diffusion_calculates_ionization = diffusion_calculates_ionization
 s% ctrl% diffusion_nsmooth_typical_charge = diffusion_nsmooth_typical_charge
 s% ctrl% diffusion_tol_correction_max = diffusion_tol_correction_max
 s% ctrl% diffusion_tol_correction_norm = diffusion_tol_correction_norm

 s% ctrl% diffusion_AD_dm_full_on = diffusion_AD_dm_full_on
 s% ctrl% diffusion_AD_dm_full_off = diffusion_AD_dm_full_off
 s% ctrl% diffusion_AD_boost_factor = diffusion_AD_boost_factor

 s% ctrl% diffusion_SIG_factor = diffusion_SIG_factor
 s% ctrl% diffusion_GT_factor = diffusion_GT_factor

 s% ctrl% diffusion_Vlimit_dm_full_on = diffusion_Vlimit_dm_full_on
 s% ctrl% diffusion_Vlimit_dm_full_off = diffusion_Vlimit_dm_full_off
 s% ctrl% diffusion_Vlimit = diffusion_Vlimit

 s% ctrl% diffusion_max_T_for_radaccel = diffusion_max_T_for_radaccel
 s% ctrl% diffusion_min_T_for_radaccel = diffusion_min_T_for_radaccel
 s% ctrl% diffusion_max_Z_for_radaccel = diffusion_max_Z_for_radaccel
 s% ctrl% diffusion_min_Z_for_radaccel = diffusion_min_Z_for_radaccel
 s% ctrl% diffusion_screening_for_radaccel = diffusion_screening_for_radaccel
 s% ctrl% op_mono_data_path = op_mono_data_path
 s% ctrl% op_mono_data_cache_filename = op_mono_data_cache_filename

 s% ctrl% show_diffusion_info = show_diffusion_info
 s% ctrl% show_diffusion_substep_info = show_diffusion_substep_info
 s% ctrl% show_diffusion_timing = show_diffusion_timing

 s% ctrl% diffusion_num_classes = diffusion_num_classes
 s% ctrl% diffusion_class_representative = diffusion_class_representative
 s% ctrl% diffusion_class_A_max = diffusion_class_A_max
 s% ctrl% diffusion_class_typical_charge = diffusion_class_typical_charge
 s% ctrl% diffusion_class_factor = diffusion_class_factor

 s% ctrl% diffusion_use_isolve = diffusion_use_isolve
 s% ctrl% diffusion_rtol_for_isolve = diffusion_rtol_for_isolve
 s% ctrl% diffusion_atol_for_isolve = diffusion_atol_for_isolve
 s% ctrl% diffusion_maxsteps_for_isolve = diffusion_maxsteps_for_isolve
 s% ctrl% diffusion_isolve_solver = diffusion_isolve_solver

 ! WD phase separation
 s% ctrl% do_phase_separation = do_phase_separation
 s% ctrl% do_phase_separation_heating = do_phase_separation_heating
 s% ctrl% phase_separation_mixing_use_brunt = phase_separation_mixing_use_brunt
 s% ctrl% phase_separation_no_diffusion = phase_separation_no_diffusion

 ! eos controls
 s% ctrl% fix_d_eos_dxa_partials = fix_d_eos_dxa_partials

 ! opacity controls
 s% ctrl% use_simple_es_for_kap = use_simple_es_for_kap
 s% ctrl% use_starting_composition_for_kap = use_starting_composition_for_kap

 s% ctrl% min_kap_for_dPrad_dm_eqn = min_kap_for_dPrad_dm_eqn
 s% ctrl% low_logT_op_mono_full_off = low_logT_op_mono_full_off
 s% ctrl% low_logT_op_mono_full_on = low_logT_op_mono_full_on
 s% ctrl% high_logT_op_mono_full_off = high_logT_op_mono_full_off
 s% ctrl% high_logT_op_mono_full_on = high_logT_op_mono_full_on
 s% ctrl% op_mono_min_X_to_include = op_mono_min_X_to_include
 s% ctrl% use_op_mono_alt_get_kap = use_op_mono_alt_get_kap
  
 s% ctrl% include_L_in_correction_limits = include_L_in_correction_limits
 s% ctrl% include_v_in_correction_limits = include_v_in_correction_limits
 s% ctrl% include_u_in_correction_limits = include_u_in_correction_limits
 s% ctrl% include_w_in_correction_limits = include_w_in_correction_limits

 ! asteroseismology controls

 s% ctrl% get_delta_nu_from_scaled_solar = get_delta_nu_from_scaled_solar
 s% ctrl% nu_max_sun = nu_max_sun
 s% ctrl% delta_nu_sun = delta_nu_sun
 s% ctrl% Teff_sun = Teff_sun
 s% ctrl% delta_Pg_mode_freq = delta_Pg_mode_freq


 ! hydro parameters
 s% ctrl% energy_eqn_option = energy_eqn_option
 s% ctrl% opacity_factor = opacity_factor
 s% ctrl% opacity_max = opacity_max
 s% ctrl% min_logT_for_opacity_factor_off = min_logT_for_opacity_factor_off
 s% ctrl% min_logT_for_opacity_factor_on = min_logT_for_opacity_factor_on
 s% ctrl% max_logT_for_opacity_factor_on = max_logT_for_opacity_factor_on
 s% ctrl% max_logT_for_opacity_factor_off = max_logT_for_opacity_factor_off

 s% ctrl% dxdt_nuc_factor = dxdt_nuc_factor
 s% ctrl% non_nuc_neu_factor = non_nuc_neu_factor
 s% ctrl% use_time_centered_eps_grav = use_time_centered_eps_grav
 s% ctrl% no_dedt_form_during_relax = no_dedt_form_during_relax
 s% ctrl% dedt_eqn_r_scale = dedt_eqn_r_scale
 s% ctrl% use_mass_corrections = use_mass_corrections
 s% ctrl% use_gravity_rotation_correction = use_gravity_rotation_correction
 s% ctrl% eps_grav_factor = eps_grav_factor
 s% ctrl% eps_mdot_factor = eps_mdot_factor
 s% ctrl% include_composition_in_eps_grav = include_composition_in_eps_grav
 s% ctrl% max_abs_rel_change_surf_lnS = max_abs_rel_change_surf_lnS
 s% ctrl% max_num_surf_revisions = max_num_surf_revisions
 s% ctrl% Gamma_lnS_eps_grav_full_off = Gamma_lnS_eps_grav_full_off
 s% ctrl% Gamma_lnS_eps_grav_full_on = Gamma_lnS_eps_grav_full_on

 s% ctrl% use_dPrad_dm_form_of_T_gradient_eqn = use_dPrad_dm_form_of_T_gradient_eqn
 s% ctrl% use_gradT_actual_vs_gradT_MLT_for_T_gradient_eqn = use_gradT_actual_vs_gradT_MLT_for_T_gradient_eqn
 s% ctrl% include_P_in_velocity_time_centering = include_P_in_velocity_time_centering
 s% ctrl% include_L_in_velocity_time_centering = include_L_in_velocity_time_centering
 s% ctrl% use_P_d_1_div_rho_form_of_work_when_time_centering_velocity = use_P_d_1_div_rho_form_of_work_when_time_centering_velocity
 s% ctrl% steps_before_use_velocity_time_centering = steps_before_use_velocity_time_centering
 s% ctrl% P_theta_for_velocity_time_centering = P_theta_for_velocity_time_centering
 s% ctrl% L_theta_for_velocity_time_centering = L_theta_for_velocity_time_centering

 s% ctrl% RTI_A = RTI_A
 s% ctrl% RTI_B = RTI_B
 s% ctrl% RTI_C = RTI_C
 s% ctrl% RTI_D = RTI_D
 s% ctrl% RTI_max_alpha = RTI_max_alpha
 s% ctrl% RTI_C_X_factor = RTI_C_X_factor
 s% ctrl% RTI_C_X0_frac = RTI_C_X0_frac
 s% ctrl% RTI_dm_for_center_eta_nondecreasing = RTI_dm_for_center_eta_nondecreasing
 s% ctrl% RTI_min_dm_behind_shock_for_full_on = RTI_min_dm_behind_shock_for_full_on
 s% ctrl% RTI_energy_floor = RTI_energy_floor
 s% ctrl% RTI_D_mix_floor = RTI_D_mix_floor
 s% ctrl% RTI_min_m_for_D_mix_floor = RTI_min_m_for_D_mix_floor
 s% ctrl% RTI_log_max_boost = RTI_log_max_boost 
 s% ctrl% RTI_m_full_boost = RTI_m_full_boost
 s% ctrl% RTI_m_no_boost = RTI_m_no_boost

 s% ctrl% velocity_logT_lower_bound = velocity_logT_lower_bound
 s% ctrl% max_dt_yrs_for_velocity_logT_lower_bound = max_dt_yrs_for_velocity_logT_lower_bound
 s% ctrl% velocity_q_upper_bound = velocity_q_upper_bound

 s% ctrl% retry_for_v_above_clight = retry_for_v_above_clight

 ! solvers

 s% ctrl% tol_correction_norm = tol_correction_norm
 s% ctrl% tol_max_correction = tol_max_correction
 s% ctrl% correction_xa_limit = correction_xa_limit

 s% ctrl% tol_correction_high_T_limit = tol_correction_high_T_limit
 s% ctrl% tol_correction_norm_high_T = tol_correction_norm_high_T
 s% ctrl% tol_max_correction_high_T = tol_max_correction_high_T

 s% ctrl% tol_correction_extreme_T_limit = tol_correction_extreme_T_limit
 s% ctrl% tol_correction_norm_extreme_T = tol_correction_norm_extreme_T
 s% ctrl% tol_max_correction_extreme_T = tol_max_correction_extreme_T
 
 s% ctrl% tol_bad_max_correction = tol_bad_max_correction
 s% ctrl% bad_max_correction_series_limit = bad_max_correction_series_limit

 s% ctrl% tol_residual_norm1 = tol_residual_norm1
 s% ctrl% tol_max_residual1 = tol_max_residual1
 s% ctrl% tol_residual_norm2 = tol_residual_norm2
 s% ctrl% tol_max_residual2 = tol_max_residual2
 s% ctrl% tol_residual_norm3 = tol_residual_norm3
 s% ctrl% tol_max_residual3 = tol_max_residual3
 s% ctrl% warning_limit_for_max_residual = warning_limit_for_max_residual
 s% ctrl% trace_solver_damping = trace_solver_damping
 
 s% ctrl% relax_use_gold_tolerances = relax_use_gold_tolerances
 s% ctrl% relax_tol_correction_norm = relax_tol_correction_norm
 s% ctrl% relax_tol_max_correction = relax_tol_max_correction
 s% ctrl% relax_solver_iters_timestep_limit = relax_solver_iters_timestep_limit
 s% ctrl% relax_iter_for_resid_tol2 = relax_iter_for_resid_tol2
 s% ctrl% relax_tol_residual_norm1 = relax_tol_residual_norm1
 s% ctrl% relax_tol_max_residual1 = relax_tol_max_residual1
 s% ctrl% relax_iter_for_resid_tol3 = relax_iter_for_resid_tol3
 s% ctrl% relax_tol_residual_norm2 = relax_tol_residual_norm2
 s% ctrl% relax_tol_max_residual2 = relax_tol_max_residual2
 s% ctrl% relax_tol_residual_norm3 = relax_tol_residual_norm3
 s% ctrl% relax_tol_max_residual3 = relax_tol_max_residual3
 s% ctrl% relax_maxT_for_gold_tolerances = relax_maxT_for_gold_tolerances
 
 s% ctrl% use_gold_tolerances = use_gold_tolerances
 s% ctrl% gold_solver_iters_timestep_limit = gold_solver_iters_timestep_limit 
 s% ctrl% maxT_for_gold_tolerances = maxT_for_gold_tolerances
 s% ctrl% gold_tol_residual_norm1 = gold_tol_residual_norm1
 s% ctrl% gold_tol_max_residual1 = gold_tol_max_residual1
 s% ctrl% gold_iter_for_resid_tol2 = gold_iter_for_resid_tol2
 s% ctrl% gold_tol_residual_norm2 = gold_tol_residual_norm2
 s% ctrl% gold_tol_max_residual2 = gold_tol_max_residual2
 s% ctrl% gold_iter_for_resid_tol3 = gold_iter_for_resid_tol3
 s% ctrl% gold_tol_residual_norm3 = gold_tol_residual_norm3
 s% ctrl% gold_tol_max_residual3 = gold_tol_max_residual3
 s% ctrl% steps_before_use_gold_tolerances = steps_before_use_gold_tolerances
 
 s% ctrl% use_gold2_tolerances = use_gold2_tolerances
 s% ctrl% gold2_solver_iters_timestep_limit = gold2_solver_iters_timestep_limit 
 s% ctrl% gold2_tol_residual_norm1 = gold2_tol_residual_norm1
 s% ctrl% gold2_tol_max_residual1 = gold2_tol_max_residual1
 s% ctrl% gold2_iter_for_resid_tol2 = gold2_iter_for_resid_tol2
 s% ctrl% gold2_tol_residual_norm2 = gold2_tol_residual_norm2
 s% ctrl% gold2_tol_max_residual2 = gold2_tol_max_residual2
 s% ctrl% gold2_iter_for_resid_tol3 = gold2_iter_for_resid_tol3
 s% ctrl% gold2_tol_residual_norm3 = gold2_tol_residual_norm3
 s% ctrl% gold2_tol_max_residual3 = gold2_tol_max_residual3
 s% ctrl% steps_before_use_gold2_tolerances = steps_before_use_gold2_tolerances
 
 s% ctrl% include_rotation_in_total_energy = include_rotation_in_total_energy

 s% ctrl% convergence_ignore_equL_residuals = convergence_ignore_equL_residuals
 s% ctrl% convergence_ignore_alpha_RTI_residuals = convergence_ignore_alpha_RTI_residuals

 s% ctrl% iter_for_resid_tol2 = iter_for_resid_tol2
 s% ctrl% iter_for_resid_tol3 = iter_for_resid_tol3

 s% ctrl% solver_itermin = solver_itermin
 s% ctrl% solver_itermin_until_reduce_min_corr_coeff = solver_itermin_until_reduce_min_corr_coeff
 s% ctrl% solver_reduced_min_corr_coeff = solver_reduced_min_corr_coeff
 s% ctrl% do_solver_damping_for_neg_xa = do_solver_damping_for_neg_xa
 s% ctrl% scale_max_correction_for_negative_surf_lum = scale_max_correction_for_negative_surf_lum
 s% ctrl% max_frac_for_negative_surf_lum = max_frac_for_negative_surf_lum
 s% ctrl% hydro_mtx_max_allowed_abs_dlogT = hydro_mtx_max_allowed_abs_dlogT
 s% ctrl% hydro_mtx_max_allowed_abs_dlogRho = hydro_mtx_max_allowed_abs_dlogRho
 s% ctrl% min_logT_for_hydro_mtx_max_allowed = min_logT_for_hydro_mtx_max_allowed
 s% ctrl% hydro_mtx_max_allowed_logT = hydro_mtx_max_allowed_logT
 s% ctrl% hydro_mtx_max_allowed_logRho = hydro_mtx_max_allowed_logRho
 s% ctrl% hydro_mtx_min_allowed_logT = hydro_mtx_min_allowed_logT
 s% ctrl% hydro_mtx_min_allowed_logRho = hydro_mtx_min_allowed_logRho
 
 s% ctrl% use_DGESVX_in_bcyclic = use_DGESVX_in_bcyclic
 s% ctrl% use_equilibration_in_DGESVX = use_equilibration_in_DGESVX
 s% ctrl% report_min_rcond_from_DGESXV = report_min_rcond_from_DGESXV
 
 s% ctrl% op_split_burn = op_split_burn
 s% ctrl% op_split_burn_min_T = op_split_burn_min_T
 s% ctrl% op_split_burn_eps = op_split_burn_eps
 s% ctrl% op_split_burn_odescal = op_split_burn_odescal
 s% ctrl% op_split_burn_min_T_for_variable_T_solver = op_split_burn_min_T_for_variable_T_solver

 s% ctrl% tiny_corr_coeff_limit = tiny_corr_coeff_limit
 s% ctrl% scale_correction_norm = scale_correction_norm
 s% ctrl% num_times_solver_reuse_mtx = num_times_solver_reuse_mtx
 s% ctrl% corr_param_factor = corr_param_factor
 s% ctrl% scale_max_correction = scale_max_correction
 s% ctrl% ignore_min_corr_coeff_for_scale_max_correction = ignore_min_corr_coeff_for_scale_max_correction
 s% ctrl% ignore_too_large_correction = ignore_too_large_correction
 s% ctrl% ignore_species_in_max_correction = ignore_species_in_max_correction

 s% ctrl% corr_norm_jump_limit = corr_norm_jump_limit
 s% ctrl% max_corr_jump_limit = max_corr_jump_limit
 s% ctrl% resid_norm_jump_limit = resid_norm_jump_limit
 s% ctrl% max_resid_jump_limit = max_resid_jump_limit

 s% ctrl% corr_coeff_limit = corr_coeff_limit
 s% ctrl% tiny_corr_factor = tiny_corr_factor

 s% ctrl% solver_max_tries_before_reject = solver_max_tries_before_reject
 s% ctrl% max_tries1 = max_tries1
 s% ctrl% max_tries_for_retry = max_tries_for_retry
 s% ctrl% max_tries_after_5_retries = max_tries_after_5_retries
 s% ctrl% max_tries_after_10_retries = max_tries_after_10_retries
 s% ctrl% max_tries_after_20_retries = max_tries_after_20_retries
 s% ctrl% retry_limit = retry_limit
 s% ctrl% redo_limit = redo_limit

 s% ctrl% use_Pvsc_art_visc = use_Pvsc_art_visc
 s% ctrl% Pvsc_cq = Pvsc_cq
 s% ctrl% Pvsc_zsh = Pvsc_zsh

 s% ctrl% min_xa_hard_limit = min_xa_hard_limit
 s% ctrl% min_xa_hard_limit_for_highT = min_xa_hard_limit_for_highT
 s% ctrl% logT_max_for_min_xa_hard_limit = logT_max_for_min_xa_hard_limit
 s% ctrl% logT_min_for_min_xa_hard_limit_for_highT = logT_min_for_min_xa_hard_limit_for_highT

 s% ctrl% sum_xa_hard_limit = sum_xa_hard_limit
 s% ctrl% sum_xa_hard_limit_for_highT = sum_xa_hard_limit_for_highT
 s% ctrl% logT_max_for_sum_xa_hard_limit = logT_max_for_sum_xa_hard_limit
 s% ctrl% logT_min_for_sum_xa_hard_limit_for_highT = logT_min_for_sum_xa_hard_limit_for_highT

 s% ctrl% xa_clip_limit = xa_clip_limit
 s% ctrl% report_solver_progress = report_solver_progress
 s% ctrl% solver_test_partials_call_number = solver_test_partials_call_number
 s% ctrl% solver_test_partials_iter_number = solver_test_partials_iter_number
 s% ctrl% solver_epsder_chem = solver_epsder_chem
 s% ctrl% solver_epsder_struct = solver_epsder_struct
 s% ctrl% solver_numerical_jacobian = solver_numerical_jacobian
 s% ctrl% solver_jacobian_nzlo = solver_jacobian_nzlo
 s% ctrl% solver_jacobian_nzhi = solver_jacobian_nzhi
 s% ctrl% solver_check_everything = solver_check_everything
 s% ctrl% energy_conservation_dump_model_number = energy_conservation_dump_model_number
 s% ctrl% solver_inspect_soln_flag = solver_inspect_soln_flag
 s% ctrl% solver_test_partials_dx_0 = solver_test_partials_dx_0
 s% ctrl% solver_test_partials_k = solver_test_partials_k
 s% ctrl% solver_test_partials_k_low = solver_test_partials_k_low
 s% ctrl% solver_test_partials_k_high = solver_test_partials_k_high
 s% ctrl% solver_show_correction_info = solver_show_correction_info
 s% ctrl% solver_test_partials_write_eos_call_info = solver_test_partials_write_eos_call_info
 s% ctrl% solver_test_eos_partials = solver_test_eos_partials
 s% ctrl% solver_test_kap_partials = solver_test_kap_partials
 s% ctrl% solver_test_net_partials = solver_test_net_partials
 s% ctrl% solver_test_atm_partials = solver_test_atm_partials
 s% ctrl% solver_test_partials_var_name = solver_test_partials_var_name
 s% ctrl% solver_test_partials_sink_name = solver_test_partials_sink_name
 s% ctrl% solver_test_partials_equ_name = solver_test_partials_equ_name
 s% ctrl% solver_test_partials_show_dx_var_name = solver_test_partials_show_dx_var_name
 s% ctrl% solver_save_photo_call_number = solver_save_photo_call_number
 s% ctrl% fill_arrays_with_NaNs = fill_arrays_with_NaNs
 s% ctrl% zero_when_allocate = zero_when_allocate
 s% ctrl% warn_when_large_rel_run_E_err = warn_when_large_rel_run_E_err
 s% ctrl% warn_when_large_virial_thm_rel_err = warn_when_large_virial_thm_rel_err
 s% ctrl% warn_when_get_a_bad_eos_result = warn_when_get_a_bad_eos_result
 s% ctrl% warn_rates_for_high_temp = warn_rates_for_high_temp
 s% ctrl% max_safe_logT_for_rates = max_safe_logT_for_rates
 s% ctrl% eps_mdot_leak_frac_factor = eps_mdot_leak_frac_factor

 s% ctrl% alpha_TDC_DAMP = alpha_TDC_DAMP
 s% ctrl% alpha_TDC_DAMPR = alpha_TDC_DAMPR
 s% ctrl% alpha_TDC_PtdVdt = alpha_TDC_PtdVdt
 s% ctrl% compare_TDC_to_MLT = compare_TDC_to_MLT

 s% ctrl% RSP2_alfap = RSP2_alfap
 s% ctrl% RSP2_alfad = RSP2_alfad
 s% ctrl% RSP2_alfat = RSP2_alfat 
 s% ctrl% RSP2_alfam = RSP2_alfam
 s% ctrl% RSP2_alfar = RSP2_alfar
 s% ctrl% RSP2_min_Lt_div_L_for_overshooting_mixing_type = RSP2_min_Lt_div_L_for_overshooting_mixing_type
 s% ctrl% RSP2_min_Lc_div_L_for_convective_mixing_type = RSP2_min_Lc_div_L_for_convective_mixing_type
 s% ctrl% RSP2_Lsurf_factor = RSP2_Lsurf_factor
 s% ctrl% RSP2_use_Stellingwerf_Lr = RSP2_use_Stellingwerf_Lr
 s% ctrl% RSP2_remesh_when_load = RSP2_remesh_when_load
 s% ctrl% RSP2_use_L_eqn_at_surface = RSP2_use_L_eqn_at_surface
 s% ctrl% RSP2_report_adjust_w = RSP2_report_adjust_w
 s% ctrl% RSP2_assume_HSE = RSP2_assume_HSE
 s% ctrl% RSP2_use_RSP_eqn_for_Y_face = RSP2_use_RSP_eqn_for_Y_face
 s% ctrl% RSP2_use_mass_interp_face_values = RSP2_use_mass_interp_face_values
 s% ctrl% RSP2_num_outermost_cells_forced_nonturbulent = RSP2_num_outermost_cells_forced_nonturbulent
 s% ctrl% RSP2_num_innermost_cells_forced_nonturbulent = RSP2_num_innermost_cells_forced_nonturbulent
 s% ctrl% RSP2_T_anchor = RSP2_T_anchor
 s% ctrl% RSP2_dq_1_factor = RSP2_dq_1_factor
 s% ctrl% RSP2_nz = RSP2_nz
 s% ctrl% RSP2_nz_outer = RSP2_nz_outer
 s% ctrl% RSP2_nz_div_IBOTOM = RSP2_nz_div_IBOTOM
 s% ctrl% RSP2_target_steps_per_cycle = RSP2_target_steps_per_cycle
 s% ctrl% RSP2_max_num_periods = RSP2_max_num_periods
 s% ctrl% RSP2_work_period = RSP2_work_period
 s% ctrl% RSP2_map_first_period = RSP2_map_first_period
 s% ctrl% RSP2_map_last_period = RSP2_map_last_period
 s% ctrl% RSP2_min_max_R_for_periods = RSP2_min_max_R_for_periods
 s% ctrl% RSP2_GREKM_avg_abs_frac_new = RSP2_GREKM_avg_abs_frac_new
 s% ctrl% RSP2_GREKM_avg_abs_limit = RSP2_GREKM_avg_abs_limit
 s% ctrl% RSP2_map_zone_interval = RSP2_map_zone_interval
 s% ctrl% RSP2_work_filename = RSP2_work_filename
 s% ctrl% RSP2_map_columns_filename = RSP2_map_columns_filename
 s% ctrl% RSP2_map_filename = RSP2_map_filename
 s% ctrl% RSP2_map_history_filename = RSP2_map_history_filename
 s% ctrl% RSP2_write_map = RSP2_write_map
 s% ctrl% RSP2_w_min_for_damping = RSP2_w_min_for_damping
 s% ctrl% RSP2_source_seed = RSP2_source_seed
 s% ctrl% RSP2_w_fix_if_neg = RSP2_w_fix_if_neg
 
 s% ctrl% max_X_for_conv_timescale = max_X_for_conv_timescale
 s% ctrl% min_X_for_conv_timescale = min_X_for_conv_timescale
 s% ctrl% max_q_for_conv_timescale = max_q_for_conv_timescale
 s% ctrl% min_q_for_conv_timescale = min_q_for_conv_timescale
 s% ctrl% max_q_for_QHSE_timescale = max_q_for_QHSE_timescale
 s% ctrl% min_q_for_QHSE_timescale = min_q_for_QHSE_timescale

 ! timestep
 s% ctrl% max_timestep = max_timestep
 s% ctrl% max_years_for_timestep = max_years_for_timestep

 s% ctrl% hi_T_max_years_for_timestep = hi_T_max_years_for_timestep
 s% ctrl% max_timestep_hi_T_limit = max_timestep_hi_T_limit

 s% ctrl% min_timestep_factor = min_timestep_factor
 s% ctrl% max_timestep_factor = max_timestep_factor
 s% ctrl% max_timestep_factor_at_high_T = max_timestep_factor_at_high_T
 s% ctrl% min_logT_for_max_timestep_factor_at_high_T = min_logT_for_max_timestep_factor_at_high_T
 s% ctrl% time_delta_coeff = time_delta_coeff
 s% ctrl% timestep_factor_for_retries = timestep_factor_for_retries
 s% ctrl% retry_hold = retry_hold
 s% ctrl% neg_mass_fraction_hold = neg_mass_fraction_hold
 s% ctrl% timestep_dt_factor = timestep_dt_factor
 s% ctrl% use_dt_low_pass_controller = use_dt_low_pass_controller
 
 s% ctrl% force_timestep_min = force_timestep_min
 s% ctrl% force_timestep_min_years = force_timestep_min_years
 s% ctrl% force_timestep_min_factor = force_timestep_min_factor
 s% ctrl% force_timestep = force_timestep
 s% ctrl% force_timestep_years = force_timestep_years

 s% ctrl% varcontrol_target = varcontrol_target
 s% ctrl% min_allowed_varcontrol_target = min_allowed_varcontrol_target
 s% ctrl% varcontrol_dt_limit_ratio_hard_max = varcontrol_dt_limit_ratio_hard_max
 s% ctrl% xa_scale = xa_scale

 s% ctrl% solver_iters_timestep_limit = solver_iters_timestep_limit

 s% ctrl% burn_steps_limit = burn_steps_limit
 s% ctrl% burn_steps_hard_limit = burn_steps_hard_limit

 s% ctrl% diffusion_steps_limit = diffusion_steps_limit
 s% ctrl% diffusion_steps_hard_limit = diffusion_steps_hard_limit
 s% ctrl% diffusion_iters_limit = diffusion_iters_limit
 s% ctrl% diffusion_iters_hard_limit = diffusion_iters_hard_limit

 s% ctrl% dt_div_dt_cell_collapse_limit = dt_div_dt_cell_collapse_limit
 s% ctrl% dt_div_dt_cell_collapse_hard_limit = dt_div_dt_cell_collapse_hard_limit
 s% ctrl% dt_div_min_dr_div_cs_limit = dt_div_min_dr_div_cs_limit
 s% ctrl% dt_div_min_dr_div_cs_hard_limit = dt_div_min_dr_div_cs_hard_limit
 
 s% ctrl% min_abs_du_div_cs_for_dt_div_min_dr_div_cs_limit = min_abs_du_div_cs_for_dt_div_min_dr_div_cs_limit
 s% ctrl% min_abs_u_div_cs_for_dt_div_min_dr_div_cs_limit = min_abs_u_div_cs_for_dt_div_min_dr_div_cs_limit
 s% ctrl% min_k_for_dt_div_min_dr_div_cs_limit = min_k_for_dt_div_min_dr_div_cs_limit
 s% ctrl% min_q_for_dt_div_min_dr_div_cs_limit = min_q_for_dt_div_min_dr_div_cs_limit
 s% ctrl% max_q_for_dt_div_min_dr_div_cs_limit = max_q_for_dt_div_min_dr_div_cs_limit
 s% ctrl% check_remnant_only_for_dt_div_min_dr_div_cs_limit = check_remnant_only_for_dt_div_min_dr_div_cs_limit

 s% ctrl% dX_mix_dist_limit = dX_mix_dist_limit

 s% ctrl% dH_limit_min_H = dH_limit_min_H
 s% ctrl% dH_limit = dH_limit
 s% ctrl% dH_hard_limit = dH_hard_limit
 s% ctrl% dH_div_H_limit_min_H = dH_div_H_limit_min_H
 s% ctrl% dH_div_H_limit = dH_div_H_limit
 s% ctrl% dH_div_H_hard_limit = dH_div_H_hard_limit
 s% ctrl% dH_decreases_only = dH_decreases_only

 s% ctrl% dHe_limit_min_He = dHe_limit_min_He
 s% ctrl% dHe_limit = dHe_limit
 s% ctrl% dHe_hard_limit = dHe_hard_limit
 s% ctrl% dHe_div_He_limit_min_He = dHe_div_He_limit_min_He
 s% ctrl% dHe_div_He_limit = dHe_div_He_limit
 s% ctrl% dHe_div_He_hard_limit = dHe_div_He_hard_limit
 s% ctrl% dHe_decreases_only = dHe_decreases_only

 s% ctrl% dHe3_limit_min_He3 = dHe3_limit_min_He3
 s% ctrl% dHe3_limit = dHe3_limit
 s% ctrl% dHe3_hard_limit = dHe3_hard_limit
 s% ctrl% dHe3_div_He3_limit_min_He3 = dHe3_div_He3_limit_min_He3
 s% ctrl% dHe3_div_He3_limit = dHe3_div_He3_limit
 s% ctrl% dHe3_div_He3_hard_limit = dHe3_div_He3_hard_limit
 s% ctrl% dHe3_decreases_only = dHe3_decreases_only

 s% ctrl% dX_limit_min_X = dX_limit_min_X
 s% ctrl% dX_limit = dX_limit
 s% ctrl% dX_hard_limit = dX_hard_limit
 s% ctrl% dX_div_X_limit_min_X = dX_div_X_limit_min_X
 s% ctrl% dX_div_X_limit = dX_div_X_limit
 s% ctrl% dX_div_X_hard_limit = dX_div_X_hard_limit
 s% ctrl% dX_div_X_at_high_T_limit = dX_div_X_at_high_T_limit
 s% ctrl% dX_div_X_at_high_T_hard_limit = dX_div_X_at_high_T_hard_limit
 s% ctrl% dX_div_X_at_high_T_limit_lgT_min = dX_div_X_at_high_T_limit_lgT_min
 
 s% ctrl% dX_decreases_only = dX_decreases_only

 s% ctrl% dX_nuc_drop_min_X_limit = dX_nuc_drop_min_X_limit
 s% ctrl% dX_nuc_drop_max_A_limit = dX_nuc_drop_max_A_limit
 s% ctrl% dX_nuc_drop_limit = dX_nuc_drop_limit
 s% ctrl% dX_nuc_drop_limit_at_high_T = dX_nuc_drop_limit_at_high_T
 s% ctrl% dX_nuc_drop_hard_limit = dX_nuc_drop_hard_limit
 s% ctrl% dX_nuc_drop_min_yrs_for_dt = dX_nuc_drop_min_yrs_for_dt

 s% ctrl% dL_div_L_limit_min_L = dL_div_L_limit_min_L
 s% ctrl% dL_div_L_limit = dL_div_L_limit
 s% ctrl% dL_div_L_hard_limit = dL_div_L_hard_limit

 s% ctrl% delta_lgP_limit = delta_lgP_limit
 s% ctrl% delta_lgP_hard_limit = delta_lgP_hard_limit
 s% ctrl% delta_lgP_limit_min_lgP = delta_lgP_limit_min_lgP

 s% ctrl% delta_lgRho_limit = delta_lgRho_limit
 s% ctrl% delta_lgRho_hard_limit = delta_lgRho_hard_limit
 s% ctrl% delta_lgRho_limit_min_lgRho = delta_lgRho_limit_min_lgRho

 s% ctrl% delta_lgT_limit = delta_lgT_limit
 s% ctrl% delta_lgT_hard_limit = delta_lgT_hard_limit
 s% ctrl% delta_lgT_limit_min_lgT = delta_lgT_limit_min_lgT

 s% ctrl% delta_lgE_limit = delta_lgE_limit
 s% ctrl% delta_lgE_hard_limit = delta_lgE_hard_limit
 s% ctrl% delta_lgE_limit_min_lgE = delta_lgE_limit_min_lgE

 s% ctrl% delta_lgR_limit = delta_lgR_limit
 s% ctrl% delta_lgR_hard_limit = delta_lgR_hard_limit
 s% ctrl% delta_lgR_limit_min_lgR = delta_lgR_limit_min_lgR

 s% ctrl% delta_Ye_highT_limit = delta_Ye_highT_limit
 s% ctrl% delta_Ye_highT_hard_limit = delta_Ye_highT_hard_limit
 s% ctrl% minT_for_highT_Ye_limit = minT_for_highT_Ye_limit

 s% ctrl% delta_lgL_nuc_cat_limit = delta_lgL_nuc_cat_limit
 s% ctrl% delta_lgL_nuc_cat_hard_limit = delta_lgL_nuc_cat_hard_limit
 s% ctrl% lgL_nuc_cat_burn_min = lgL_nuc_cat_burn_min
 s% ctrl% lgL_nuc_mix_dist_limit = lgL_nuc_mix_dist_limit

 s% ctrl% delta_lgL_H_limit = delta_lgL_H_limit
 s% ctrl% delta_lgL_H_hard_limit = delta_lgL_H_hard_limit
 s% ctrl% lgL_H_burn_min = lgL_H_burn_min
 s% ctrl% lgL_H_drop_factor = lgL_H_drop_factor
 s% ctrl% lgL_H_burn_relative_limit = lgL_H_burn_relative_limit

 s% ctrl% delta_lgL_He_limit = delta_lgL_He_limit
 s% ctrl% delta_lgL_He_hard_limit = delta_lgL_He_hard_limit
 s% ctrl% lgL_He_burn_min = lgL_He_burn_min
 s% ctrl% lgL_He_drop_factor = lgL_He_drop_factor
 s% ctrl% lgL_He_burn_relative_limit = lgL_He_burn_relative_limit

 s% ctrl% delta_lgL_z_limit = delta_lgL_z_limit
 s% ctrl% delta_lgL_z_hard_limit = delta_lgL_z_hard_limit
 s% ctrl% lgL_z_burn_min = lgL_z_burn_min
 s% ctrl% lgL_z_drop_factor = lgL_z_drop_factor
 s% ctrl% lgL_z_burn_relative_limit = lgL_z_burn_relative_limit

 s% ctrl% delta_lgL_power_photo_limit = delta_lgL_power_photo_limit
 s% ctrl% delta_lgL_power_photo_hard_limit = delta_lgL_power_photo_hard_limit
 s% ctrl% lgL_power_photo_burn_min = lgL_power_photo_burn_min
 s% ctrl% lgL_power_photo_drop_factor = lgL_power_photo_drop_factor
 s% ctrl% min_lgT_for_lgL_power_photo_limit = min_lgT_for_lgL_power_photo_limit

 s% ctrl% delta_lgL_nuc_limit = delta_lgL_nuc_limit
 s% ctrl% delta_lgL_nuc_hard_limit = delta_lgL_nuc_hard_limit
 s% ctrl% delta_lgL_nuc_at_high_T_limit = delta_lgL_nuc_at_high_T_limit
 s% ctrl% delta_lgL_nuc_at_high_T_hard_limit = delta_lgL_nuc_at_high_T_hard_limit
 s% ctrl% delta_lgL_nuc_at_high_T_limit_lgT_min = delta_lgL_nuc_at_high_T_limit_lgT_min
 
 s% ctrl% max_lgT_for_lgL_nuc_limit = max_lgT_for_lgL_nuc_limit
 s% ctrl% lgL_nuc_burn_min = lgL_nuc_burn_min
 s% ctrl% lgL_nuc_drop_factor = lgL_nuc_drop_factor

 s% ctrl% delta_lgRho_cntr_limit = delta_lgRho_cntr_limit
 s% ctrl% delta_lgRho_cntr_hard_limit = delta_lgRho_cntr_hard_limit

 s% ctrl% delta_lgP_cntr_limit = delta_lgP_cntr_limit
 s% ctrl% delta_lgP_cntr_hard_limit = delta_lgP_cntr_hard_limit

 s% ctrl% delta_lgT_cntr_limit = delta_lgT_cntr_limit
 s% ctrl% delta_lgT_cntr_hard_limit = delta_lgT_cntr_hard_limit
 s% ctrl% delta_lgT_cntr_limit_only_after_near_zams = delta_lgT_cntr_limit_only_after_near_zams

 s% ctrl% delta_lgT_max_limit = delta_lgT_max_limit
 s% ctrl% delta_lgT_max_hard_limit = delta_lgT_max_hard_limit
 s% ctrl% delta_lgT_max_limit_lgT_min = delta_lgT_max_limit_lgT_min
 s% ctrl% delta_lgT_max_limit_only_after_near_zams = delta_lgT_max_limit_only_after_near_zams

 s% ctrl% delta_lgT_max_at_high_T_limit = delta_lgT_max_at_high_T_limit
 s% ctrl% delta_lgT_max_at_high_T_hard_limit = delta_lgT_max_at_high_T_hard_limit
 s% ctrl% delta_lgT_max_at_high_T_limit_lgT_min = delta_lgT_max_at_high_T_limit_lgT_min

 s% ctrl% delta_log_eps_nuc_limit = delta_log_eps_nuc_limit
 s% ctrl% delta_log_eps_nuc_hard_limit = delta_log_eps_nuc_hard_limit

 s% ctrl% delta_dX_div_X_cntr_min = delta_dX_div_X_cntr_min
 s% ctrl% delta_dX_div_X_cntr_max = delta_dX_div_X_cntr_max
 s% ctrl% delta_dX_div_X_cntr_limit = delta_dX_div_X_cntr_limit
 s% ctrl% delta_dX_div_X_cntr_hard_limit = delta_dX_div_X_cntr_hard_limit

 s% ctrl% delta_dX_div_X_drop_only = delta_dX_div_X_drop_only
 s% ctrl% delta_lg_XH_drop_only = delta_lg_XH_drop_only
 s% ctrl% delta_lg_XHe_drop_only = delta_lg_XHe_drop_only
 s% ctrl% delta_lg_XC_drop_only = delta_lg_XC_drop_only
 s% ctrl% delta_lg_XNe_drop_only = delta_lg_XNe_drop_only
 s% ctrl% delta_lg_XO_drop_only = delta_lg_XO_drop_only
 s% ctrl% delta_lg_XSi_drop_only = delta_lg_XSi_drop_only
 s% ctrl% delta_XH_drop_only = delta_XH_drop_only
 s% ctrl% delta_XHe_drop_only = delta_XHe_drop_only
 s% ctrl% delta_XC_drop_only = delta_XC_drop_only
 s% ctrl% delta_XNe_drop_only = delta_XNe_drop_only
 s% ctrl% delta_XO_drop_only = delta_XO_drop_only
 s% ctrl% delta_XSi_drop_only = delta_XSi_drop_only

 s% ctrl% delta_lg_XH_cntr_min = delta_lg_XH_cntr_min
 s% ctrl% delta_lg_XH_cntr_max = delta_lg_XH_cntr_max
 s% ctrl% delta_lg_XH_cntr_limit = delta_lg_XH_cntr_limit
 s% ctrl% delta_lg_XH_cntr_hard_limit = delta_lg_XH_cntr_hard_limit

 s% ctrl% delta_lg_XHe_cntr_min = delta_lg_XHe_cntr_min
 s% ctrl% delta_lg_XHe_cntr_max = delta_lg_XHe_cntr_max
 s% ctrl% delta_lg_XHe_cntr_limit = delta_lg_XHe_cntr_limit
 s% ctrl% delta_lg_XHe_cntr_hard_limit = delta_lg_XHe_cntr_hard_limit

 s% ctrl% delta_lg_XC_cntr_min = delta_lg_XC_cntr_min
 s% ctrl% delta_lg_XC_cntr_max = delta_lg_XC_cntr_max
 s% ctrl% delta_lg_XC_cntr_limit = delta_lg_XC_cntr_limit
 s% ctrl% delta_lg_XC_cntr_hard_limit = delta_lg_XC_cntr_hard_limit

 s% ctrl% delta_lg_XNe_cntr_limit = delta_lg_XNe_cntr_limit
 s% ctrl% delta_lg_XNe_cntr_hard_limit = delta_lg_XNe_cntr_hard_limit
 s% ctrl% delta_lg_XNe_cntr_min = delta_lg_XNe_cntr_min
 s% ctrl% delta_lg_XNe_cntr_max = delta_lg_XNe_cntr_max

 s% ctrl% delta_lg_XO_cntr_limit = delta_lg_XO_cntr_limit
 s% ctrl% delta_lg_XO_cntr_hard_limit = delta_lg_XO_cntr_hard_limit
 s% ctrl% delta_lg_XO_cntr_min = delta_lg_XO_cntr_min
 s% ctrl% delta_lg_XO_cntr_max = delta_lg_XO_cntr_max

 s% ctrl% delta_lg_XSi_cntr_limit = delta_lg_XSi_cntr_limit
 s% ctrl% delta_lg_XSi_cntr_hard_limit = delta_lg_XSi_cntr_hard_limit
 s% ctrl% delta_lg_XSi_cntr_min = delta_lg_XSi_cntr_min
 s% ctrl% delta_lg_XSi_cntr_max = delta_lg_XSi_cntr_max

 s% ctrl% delta_XH_cntr_limit = delta_XH_cntr_limit
 s% ctrl% delta_XH_cntr_hard_limit = delta_XH_cntr_hard_limit
 s% ctrl% delta_XHe_cntr_limit = delta_XHe_cntr_limit
 s% ctrl% delta_XHe_cntr_hard_limit = delta_XHe_cntr_hard_limit
 s% ctrl% delta_XC_cntr_limit = delta_XC_cntr_limit
 s% ctrl% delta_XC_cntr_hard_limit = delta_XC_cntr_hard_limit
 s% ctrl% delta_XNe_cntr_limit = delta_XNe_cntr_limit
 s% ctrl% delta_XNe_cntr_hard_limit = delta_XNe_cntr_hard_limit
 s% ctrl% delta_XO_cntr_limit = delta_XO_cntr_limit
 s% ctrl% delta_XO_cntr_hard_limit = delta_XO_cntr_hard_limit
 s% ctrl% delta_XSi_cntr_limit = delta_XSi_cntr_limit
 s% ctrl% delta_XSi_cntr_hard_limit = delta_XSi_cntr_hard_limit

 s% ctrl% delta_lgTeff_limit = delta_lgTeff_limit
 s% ctrl% delta_lgTeff_hard_limit = delta_lgTeff_hard_limit

 s% ctrl% delta_lgL_limit = delta_lgL_limit
 s% ctrl% delta_lgL_limit_L_min = delta_lgL_limit_L_min
 s% ctrl% delta_lgL_hard_limit = delta_lgL_hard_limit

 s% ctrl% delta_HR_ds_L = delta_HR_ds_L
 s% ctrl% delta_HR_ds_Teff = delta_HR_ds_Teff
 s% ctrl% delta_HR_limit = delta_HR_limit
 s% ctrl% delta_HR_hard_limit = delta_HR_hard_limit

 s% ctrl% delta_lg_star_mass_limit = delta_lg_star_mass_limit
 s% ctrl% delta_lg_star_mass_hard_limit = delta_lg_star_mass_hard_limit

 s% ctrl% delta_mdot_atol = delta_mdot_atol
 s% ctrl% delta_mdot_rtol = delta_mdot_rtol
 s% ctrl% delta_mdot_limit = delta_mdot_limit
 s% ctrl% delta_mdot_hard_limit = delta_mdot_hard_limit

 s% ctrl% adjust_J_q_limit = adjust_J_q_limit
 s% ctrl% adjust_J_q_hard_limit = adjust_J_q_hard_limit
 s% ctrl% never_skip_hard_limits = never_skip_hard_limits
 s% ctrl% relax_hard_limits_after_retry = relax_hard_limits_after_retry
 s% ctrl% report_dt_hard_limit_retries = report_dt_hard_limit_retries
 s% ctrl% report_min_dr_div_cs = report_min_dr_div_cs
 s% ctrl% report_solver_dt_info = report_solver_dt_info

 s% ctrl% limit_for_rel_error_in_energy_conservation = limit_for_rel_error_in_energy_conservation
 s% ctrl% hard_limit_for_rel_error_in_energy_conservation = hard_limit_for_rel_error_in_energy_conservation

 s% ctrl% min_chem_eqn_scale = min_chem_eqn_scale

 s% ctrl% trace_evolve = trace_evolve


 ! misc
 s% ctrl% zams_filename = zams_filename
 s% ctrl% set_rho_to_dm_div_dV = set_rho_to_dm_div_dV

 s% ctrl% use_other_mlt_results = use_other_mlt_results
 s% ctrl% use_other_surface_PT = use_other_surface_PT
 s% ctrl% use_other_kap = use_other_kap
 s% ctrl% use_other_diffusion = use_other_diffusion
 s% ctrl% use_other_diffusion_factor = use_other_diffusion_factor
 s% ctrl% use_other_adjust_mdot = use_other_adjust_mdot
 s% ctrl% use_other_j_for_adjust_J_lost = use_other_j_for_adjust_J_lost
 s% ctrl% use_other_alpha_mlt = use_other_alpha_mlt
 s% ctrl% use_other_am_mixing = use_other_am_mixing
 s% ctrl% use_other_brunt = use_other_brunt
 s% ctrl% use_other_brunt_smoothing = use_other_brunt_smoothing
 s% ctrl% use_other_solver_monitor = use_other_solver_monitor
 s% ctrl% use_other_build_initial_model = use_other_build_initial_model
 s% ctrl% use_other_cgrav = use_other_cgrav
 s% ctrl% use_other_mesh_delta_coeff_factor = use_other_mesh_delta_coeff_factor
 s% ctrl% use_other_energy_implicit = use_other_energy_implicit
 s% ctrl% use_other_remove_surface = use_other_remove_surface
 s% ctrl% use_other_momentum = use_other_momentum
 s% ctrl% use_other_momentum_implicit = use_other_momentum_implicit
 s% ctrl% use_other_pressure = use_other_pressure
 s% ctrl% use_other_energy = use_other_energy
 s% ctrl% use_other_mesh_functions = use_other_mesh_functions
 s% ctrl% use_other_eps_grav = use_other_eps_grav
 s% ctrl% use_other_gradr_factor = use_other_gradr_factor
 s% ctrl% use_other_D_mix = use_other_D_mix
 s% ctrl% use_other_neu = use_other_neu
 s% ctrl% use_other_net_get = use_other_net_get
 s% ctrl% use_other_opacity_factor = use_other_opacity_factor
 s% ctrl% use_other_diffusion_coefficients = use_other_diffusion_coefficients
 s% ctrl% use_other_pgstar_plots = use_other_pgstar_plots
 s% ctrl% use_other_eval_fp_ft = use_other_eval_fp_ft
 s% ctrl% use_other_eval_i_rot = use_other_eval_i_rot
 s% ctrl% use_other_torque = use_other_torque
 s% ctrl% use_other_torque_implicit = use_other_torque_implicit
 s% ctrl% use_other_wind = use_other_wind
 s% ctrl% use_other_accreting_state = use_other_accreting_state
 s% ctrl% use_other_after_struct_burn_mix = use_other_after_struct_burn_mix
 s% ctrl% use_other_before_struct_burn_mix = use_other_before_struct_burn_mix
 s% ctrl% use_other_astero_freq_corr = use_other_astero_freq_corr
 s% ctrl% use_other_timestep_limit = use_other_timestep_limit
 s% ctrl% use_other_set_pgstar_controls = use_other_set_pgstar_controls
 s% ctrl% use_other_screening = use_other_screening

 s% ctrl% x_ctrl = x_ctrl
 s% ctrl% x_integer_ctrl = x_integer_ctrl
 s% ctrl% x_logical_ctrl = x_logical_ctrl
 s% ctrl% x_character_ctrl = x_character_ctrl

 ! info for debugging
 s% ctrl% stop_for_bad_nums = stop_for_bad_nums
 s% ctrl% report_ierr = report_ierr
 s% ctrl% report_bad_negative_xa = report_bad_negative_xa

 s% ctrl% diffusion_dump_call_number = diffusion_dump_call_number

 s% ctrl% surface_accel_div_grav_limit = surface_accel_div_grav_limit
 s% ctrl% steps_before_start_stress_test = steps_before_start_stress_test
 s% ctrl% stress_test_relax = stress_test_relax

 end subroutine store_controls


 subroutine set_controls_for_writing(s, ierr)
 use star_private_def
 use chem_def ! categories
 type (star_info), pointer :: s
 integer, intent(out) :: ierr

 ierr = 0

 ! where to start
 initial_mass = s% ctrl% initial_mass
 initial_z = s% ctrl% initial_z
 initial_y = s% ctrl% initial_y
 initial_he3 = s% ctrl% initial_he3

 ! definition of core boundaries
 he_core_boundary_h1_fraction = s% ctrl% he_core_boundary_h1_fraction
 co_core_boundary_he4_fraction = s% ctrl% co_core_boundary_he4_fraction
 one_core_boundary_he4_c12_fraction = s% ctrl% one_core_boundary_he4_c12_fraction
 fe_core_boundary_si28_fraction = s% ctrl% fe_core_boundary_si28_fraction
 neutron_rich_core_boundary_Ye_max = s% ctrl% neutron_rich_core_boundary_Ye_max
 min_boundary_fraction = s% ctrl% min_boundary_fraction

 ! when to stop
 max_model_number = s% ctrl% max_model_number
 max_abs_rel_run_E_err = s% ctrl% max_abs_rel_run_E_err
 max_number_retries = s% ctrl% max_number_retries
 relax_max_number_retries = s% ctrl% relax_max_number_retries
 max_age = s% ctrl% max_age
 max_age_in_days = s% ctrl% max_age_in_days
 max_age_in_seconds = s% ctrl% max_age_in_seconds
 num_adjusted_dt_steps_before_max_age = s% ctrl% num_adjusted_dt_steps_before_max_age
 dt_years_for_steps_before_max_age = s% ctrl% dt_years_for_steps_before_max_age
 reduction_factor_for_max_timestep = s% ctrl% reduction_factor_for_max_timestep
 when_to_stop_rtol = s% ctrl% when_to_stop_rtol
 when_to_stop_atol = s% ctrl% when_to_stop_atol
 gamma_center_limit = s% ctrl% gamma_center_limit
 eta_center_limit = s% ctrl% eta_center_limit
 log_center_temp_limit = s% ctrl% log_center_temp_limit
 log_max_temp_upper_limit = s% ctrl% log_max_temp_upper_limit
 log_max_temp_lower_limit = s% ctrl% log_max_temp_lower_limit
 log_center_temp_lower_limit = s% ctrl% log_center_temp_lower_limit
 log_center_density_limit = s% ctrl% log_center_density_limit
 log_center_density_lower_limit = s% ctrl% log_center_density_lower_limit
 min_timestep_limit = s% ctrl% min_timestep_limit

 center_entropy_limit = s% ctrl% center_entropy_limit
 center_entropy_lower_limit = s% ctrl% center_entropy_lower_limit
 max_entropy_limit = s% ctrl% max_entropy_limit
 max_entropy_lower_limit = s% ctrl% max_entropy_lower_limit

 fe_core_infall_limit = s% ctrl% fe_core_infall_limit
 center_Ye_lower_limit = s% ctrl% center_Ye_lower_limit
 center_R_lower_limit = s% ctrl% center_R_lower_limit
 non_fe_core_infall_limit = s% ctrl% non_fe_core_infall_limit
 non_fe_core_rebound_limit = s% ctrl% non_fe_core_rebound_limit
 v_div_csound_surf_limit = s% ctrl% v_div_csound_surf_limit
 v_div_csound_max_limit = s% ctrl% v_div_csound_max_limit
 Lnuc_div_L_upper_limit = s% ctrl% Lnuc_div_L_upper_limit
 Lnuc_div_L_lower_limit = s% ctrl% Lnuc_div_L_lower_limit
 v_surf_div_v_kh_upper_limit = s% ctrl% v_surf_div_v_kh_upper_limit
 v_surf_div_v_kh_lower_limit = s% ctrl% v_surf_div_v_kh_lower_limit
 v_surf_div_v_esc_limit = s% ctrl% v_surf_div_v_esc_limit
 v_surf_kms_limit = s% ctrl% v_surf_kms_limit

 stop_near_zams = s% ctrl% stop_near_zams
 stop_at_phase_PreMS = s% ctrl% stop_at_phase_PreMS
 stop_at_phase_ZAMS = s% ctrl% stop_at_phase_ZAMS
 stop_at_phase_IAMS = s% ctrl% stop_at_phase_IAMS
 stop_at_phase_TAMS = s% ctrl% stop_at_phase_TAMS
 stop_at_phase_He_Burn = s% ctrl% stop_at_phase_He_Burn
 stop_at_phase_ZACHeB = s% ctrl% stop_at_phase_ZACHeB
 stop_at_phase_TACHeB = s% ctrl% stop_at_phase_TACHeB
 stop_at_phase_TP_AGB = s% ctrl% stop_at_phase_TP_AGB
 stop_at_phase_C_Burn = s% ctrl% stop_at_phase_C_Burn
 stop_at_phase_Ne_Burn = s% ctrl% stop_at_phase_Ne_Burn
 stop_at_phase_O_Burn = s% ctrl% stop_at_phase_O_Burn
 stop_at_phase_Si_Burn = s% ctrl% stop_at_phase_Si_Burn
 stop_at_phase_WDCS = s% ctrl% stop_at_phase_WDCS
 Lnuc_div_L_zams_limit = s% ctrl% Lnuc_div_L_zams_limit
 Pgas_div_P_limit = s% ctrl% Pgas_div_P_limit
 Pgas_div_P_limit_max_q = s% ctrl% Pgas_div_P_limit_max_q
 gamma1_limit = s% ctrl% gamma1_limit
 gamma1_limit_max_q = s% ctrl% gamma1_limit_max_q
 gamma1_limit_max_v_div_vesc = s% ctrl% gamma1_limit_max_v_div_vesc
 peak_burn_vconv_div_cs_limit = s% ctrl% peak_burn_vconv_div_cs_limit
 omega_div_omega_crit_limit = s% ctrl% omega_div_omega_crit_limit
 delta_nu_lower_limit = s% ctrl% delta_nu_lower_limit
 delta_nu_upper_limit = s% ctrl% delta_nu_upper_limit
 delta_Pg_lower_limit = s% ctrl% delta_Pg_lower_limit
 delta_Pg_upper_limit = s% ctrl% delta_Pg_upper_limit
 shock_mass_upper_limit = s% ctrl% shock_mass_upper_limit
 mach1_mass_upper_limit = s% ctrl% mach1_mass_upper_limit
 stop_when_reach_this_cumulative_extra_heating = s% ctrl% stop_when_reach_this_cumulative_extra_heating

 xa_central_lower_limit_species = s% ctrl% xa_central_lower_limit_species
 xa_central_lower_limit = s% ctrl% xa_central_lower_limit

 xa_central_upper_limit_species = s% ctrl% xa_central_upper_limit_species
 xa_central_upper_limit = s% ctrl% xa_central_upper_limit

 xa_surface_lower_limit_species = s% ctrl% xa_surface_lower_limit_species
 xa_surface_lower_limit = s% ctrl% xa_surface_lower_limit

 xa_surface_upper_limit_species = s% ctrl% xa_surface_upper_limit_species
 xa_surface_upper_limit = s% ctrl% xa_surface_upper_limit

 xa_average_lower_limit_species = s% ctrl% xa_average_lower_limit_species
 xa_average_lower_limit = s% ctrl% xa_average_lower_limit

 xa_average_upper_limit_species = s% ctrl% xa_average_upper_limit_species
 xa_average_upper_limit = s% ctrl% xa_average_upper_limit

 HB_limit = s% ctrl% HB_limit

 star_mass_max_limit = s% ctrl% star_mass_max_limit
 star_mass_min_limit = s% ctrl% star_mass_min_limit
 ejecta_mass_max_limit = s% ctrl% ejecta_mass_max_limit
 remnant_mass_min_limit = s% ctrl% remnant_mass_min_limit
 
 star_species_mass_min_limit = s% ctrl% star_species_mass_min_limit
 star_species_mass_min_limit_iso = s% ctrl% star_species_mass_min_limit_iso
 star_species_mass_max_limit = s% ctrl% star_species_mass_max_limit
 star_species_mass_max_limit_iso = s% ctrl% star_species_mass_max_limit_iso

 xmstar_min_limit = s% ctrl% xmstar_min_limit
 xmstar_max_limit = s% ctrl% xmstar_max_limit
 envelope_mass_limit = s% ctrl% envelope_mass_limit
 envelope_fraction_left_limit = s% ctrl% envelope_fraction_left_limit

 he_core_mass_limit = s% ctrl% he_core_mass_limit
 co_core_mass_limit = s% ctrl% co_core_mass_limit
 one_core_mass_limit = s% ctrl% one_core_mass_limit
 fe_core_mass_limit = s% ctrl% fe_core_mass_limit
 neutron_rich_core_mass_limit = s% ctrl% neutron_rich_core_mass_limit

 he_layer_mass_lower_limit = s% ctrl% he_layer_mass_lower_limit
 abs_diff_lg_LH_lg_Ls_limit = s% ctrl% abs_diff_lg_LH_lg_Ls_limit
 Teff_upper_limit = s% ctrl% Teff_upper_limit
 Teff_lower_limit = s% ctrl% Teff_lower_limit
 photosphere_m_upper_limit = s% ctrl% photosphere_m_upper_limit
 photosphere_m_lower_limit = s% ctrl% photosphere_m_lower_limit
 photosphere_m_sub_M_center_limit = s% ctrl% photosphere_m_sub_M_center_limit
 photosphere_r_upper_limit = s% ctrl% photosphere_r_upper_limit
 photosphere_r_lower_limit = s% ctrl% photosphere_r_lower_limit
 log_Teff_upper_limit = s% ctrl% log_Teff_upper_limit
 log_Teff_lower_limit = s% ctrl% log_Teff_lower_limit
 log_Tsurf_upper_limit = s% ctrl% log_Tsurf_upper_limit
 log_Tsurf_lower_limit = s% ctrl% log_Tsurf_lower_limit
 log_Rsurf_upper_limit = s% ctrl% log_Rsurf_upper_limit
 log_Rsurf_lower_limit = s% ctrl% log_Rsurf_lower_limit
 log_Psurf_upper_limit = s% ctrl% log_Psurf_upper_limit
 log_Psurf_lower_limit = s% ctrl% log_Psurf_lower_limit
 log_Dsurf_upper_limit = s% ctrl% log_Dsurf_upper_limit
 log_Dsurf_lower_limit = s% ctrl% log_Dsurf_lower_limit
 log_L_upper_limit = s% ctrl% log_L_upper_limit
 log_L_lower_limit = s% ctrl% log_L_lower_limit
 log_g_upper_limit = s% ctrl% log_g_upper_limit
 log_g_lower_limit = s% ctrl% log_g_lower_limit

 power_nuc_burn_upper_limit = s% ctrl% power_nuc_burn_upper_limit
 power_h_burn_upper_limit = s% ctrl% power_h_burn_upper_limit
 power_he_burn_upper_limit = s% ctrl% power_he_burn_upper_limit
 power_z_burn_upper_limit = s% ctrl% power_z_burn_upper_limit
 power_nuc_burn_lower_limit = s% ctrl% power_nuc_burn_lower_limit
 power_h_burn_lower_limit = s% ctrl% power_h_burn_lower_limit
 power_he_burn_lower_limit = s% ctrl% power_he_burn_lower_limit
 power_z_burn_lower_limit = s% ctrl% power_z_burn_lower_limit


 ! output of "snapshots" for restarts
 photo_interval = s% ctrl% photo_interval
 photo_digits = s% ctrl% photo_digits
 photo_directory = s% ctrl% photo_directory
 ! output of history and profiles.
 do_history_file = s% ctrl% do_history_file
 history_interval = s% ctrl% history_interval

 write_header_frequency = s% ctrl% write_header_frequency
 terminal_interval = s% ctrl% terminal_interval
 terminal_show_age_units = s% ctrl% terminal_show_age_units
 terminal_show_timestep_units = s% ctrl% terminal_show_timestep_units
 terminal_show_log_dt = s% ctrl% terminal_show_log_dt
 terminal_show_log_age = s% ctrl% terminal_show_log_age
 extra_terminal_output_file = s% ctrl% extra_terminal_output_file
 num_trace_history_values = s% ctrl% num_trace_history_values
 trace_history_value_name = s% ctrl% trace_history_value_name

 log_directory = s% ctrl% log_directory

 star_history_name = s% ctrl% star_history_name
 star_history_header_name = s% ctrl% star_history_header_name
 star_history_dbl_format = s% ctrl% star_history_dbl_format
 star_history_int_format = s% ctrl% star_history_int_format
 star_history_txt_format = s% ctrl% star_history_txt_format

 profiles_index_name = s% ctrl% profiles_index_name
 profile_data_prefix = s% ctrl% profile_data_prefix
 profile_data_suffix = s% ctrl% profile_data_suffix
 profile_data_header_suffix = s% ctrl% profile_data_header_suffix
 profile_int_format = s% ctrl% profile_int_format
 profile_txt_format = s% ctrl% profile_txt_format
 profile_dbl_format = s% ctrl% profile_dbl_format
 profile_header_include_sys_details = s% ctrl% profile_header_include_sys_details
 write_profiles_flag = s% ctrl% write_profiles_flag
 profile_interval = s% ctrl% profile_interval
 priority_profile_interval = s% ctrl% priority_profile_interval
 profile_model = s% ctrl% profile_model
 max_num_profile_models = s% ctrl% max_num_profile_models
 max_num_profile_zones = s% ctrl% max_num_profile_zones

 write_controls_info_with_profile = s% ctrl% write_controls_info_with_profile
 controls_data_prefix = s% ctrl% controls_data_prefix
 controls_data_suffix = s% ctrl% controls_data_suffix

 write_pulse_data_with_profile = s% ctrl% write_pulse_data_with_profile
 pulse_data_format = s% ctrl% pulse_data_format
 add_atmosphere_to_pulse_data = s% ctrl% add_atmosphere_to_pulse_data
 add_center_point_to_pulse_data = s% ctrl% add_center_point_to_pulse_data
 keep_surface_point_for_pulse_data = s% ctrl% keep_surface_point_for_pulse_data
 add_double_points_to_pulse_data = s% ctrl% add_double_points_to_pulse_data
 interpolate_rho_for_pulse_data = s% ctrl% interpolate_rho_for_pulse_data
 threshold_grad_mu_for_double_point = s% ctrl% threshold_grad_mu_for_double_point
 max_number_of_double_points = s% ctrl% max_number_of_double_points

 fgong_header = s% ctrl% fgong_header
 fgong_ivers = s% ctrl% fgong_ivers
 
 max_num_gyre_points = s% ctrl% max_num_gyre_points
 format_for_OSC_data = s% ctrl% format_for_OSC_data
 fgong_zero_A_inside_r = s% ctrl% fgong_zero_A_inside_r
 use_other_export_pulse_data = s% ctrl% use_other_export_pulse_data
 use_other_get_pulse_data = s% ctrl% use_other_get_pulse_data
 use_other_edit_pulse_data = s% ctrl% use_other_edit_pulse_data

 write_model_with_profile = s% ctrl% write_model_with_profile
 model_data_prefix = s% ctrl% model_data_prefix
 model_data_suffix = s% ctrl% model_data_suffix

 mixing_D_limit_for_log = s% ctrl% mixing_D_limit_for_log
 trace_mass_location = s% ctrl% trace_mass_location
 min_tau_for_max_abs_v_location = s% ctrl% min_tau_for_max_abs_v_location
 min_q_for_inner_mach1_location = s% ctrl% min_q_for_inner_mach1_location
 max_q_for_outer_mach1_location = s% ctrl% max_q_for_outer_mach1_location
 
 mass_depth_for_L_surf = s% ctrl% mass_depth_for_L_surf
 conv_core_gap_dq_limit = s% ctrl% conv_core_gap_dq_limit

 ! burn zone eps definitions for use in logs and profiles
 burn_min1 = s% ctrl% burn_min1
 burn_min2 = s% ctrl% burn_min2

 max_conv_vel_div_csound_maxq = s% ctrl% max_conv_vel_div_csound_maxq
 width_for_limit_conv_vel = s% ctrl% width_for_limit_conv_vel
 max_q_for_limit_conv_vel = s% ctrl% max_q_for_limit_conv_vel
 max_mass_in_gm_for_limit_conv_vel = s% ctrl% max_mass_in_gm_for_limit_conv_vel
 max_r_in_cm_for_limit_conv_vel = s% ctrl% max_r_in_cm_for_limit_conv_vel

 ! for reported average values
 surface_avg_abundance_dq = s% ctrl% surface_avg_abundance_dq
 center_avg_value_dq = s% ctrl% center_avg_value_dq

 ! mixing parameters
 min_convective_gap = s% ctrl% min_convective_gap
 min_thermohaline_gap = s% ctrl% min_thermohaline_gap
 min_semiconvection_gap = s% ctrl% min_semiconvection_gap
 min_thermohaline_dropout = s% ctrl% min_thermohaline_dropout
 max_dropout_gradL_sub_grada = s% ctrl% max_dropout_gradL_sub_grada
 remove_embedded_semiconvection = s% ctrl% remove_embedded_semiconvection
 recalc_mix_info_after_evolve = s% ctrl% recalc_mix_info_after_evolve
 remove_mixing_glitches = s% ctrl% remove_mixing_glitches
 okay_to_remove_mixing_singleton = s% ctrl% okay_to_remove_mixing_singleton
 prune_bad_cz_min_Hp_height = s% ctrl% prune_bad_cz_min_Hp_height
 prune_bad_cz_min_log_eps_nuc = s% ctrl% prune_bad_cz_min_log_eps_nuc
 redo_conv_for_dr_lt_mixing_length = s% ctrl% redo_conv_for_dr_lt_mixing_length

 alpha_semiconvection = s% ctrl% alpha_semiconvection
 semiconvection_option = s% ctrl% semiconvection_option
 use_Ledoux_criterion = s% ctrl% use_Ledoux_criterion
 num_cells_for_smooth_gradL_composition_term = s% ctrl% num_cells_for_smooth_gradL_composition_term
 threshold_for_smooth_gradL_composition_term = s% ctrl% threshold_for_smooth_gradL_composition_term
 clip_D_limit = s% ctrl% clip_D_limit
 fix_eps_grav_transition_to_grid = s% ctrl% fix_eps_grav_transition_to_grid

 okay_to_reduce_gradT_excess = s% ctrl% okay_to_reduce_gradT_excess
 gradT_excess_f1 = s% ctrl% gradT_excess_f1
 gradT_excess_f2 = s% ctrl% gradT_excess_f2
 gradT_excess_max_center_h1 = s% ctrl% gradT_excess_max_center_h1
 gradT_excess_min_center_he4 = s% ctrl% gradT_excess_min_center_he4
 gradT_excess_max_logT = s% ctrl% gradT_excess_max_logT
 gradT_excess_min_log_tau_full_on = s% ctrl% gradT_excess_min_log_tau_full_on
 gradT_excess_max_log_tau_full_off = s% ctrl% gradT_excess_max_log_tau_full_off
 gradT_excess_lambda1 = s% ctrl% gradT_excess_lambda1
 gradT_excess_beta1 = s% ctrl% gradT_excess_beta1
 gradT_excess_lambda2 = s% ctrl% gradT_excess_lambda2
 gradT_excess_beta2 = s% ctrl% gradT_excess_beta2
 gradT_excess_dlambda = s% ctrl% gradT_excess_dlambda
 gradT_excess_dbeta = s% ctrl% gradT_excess_dbeta
 
 D_mix_zero_region_bottom_q = s% ctrl% D_mix_zero_region_bottom_q
 D_mix_zero_region_top_q = s% ctrl% D_mix_zero_region_top_q
 dq_D_mix_zero_at_H_He_crossover = s% ctrl% dq_D_mix_zero_at_H_He_crossover
 dq_D_mix_zero_at_H_C_crossover = s% ctrl% dq_D_mix_zero_at_H_C_crossover

 use_superad_reduction = s% ctrl% use_superad_reduction
 superad_reduction_gamma_limit = s% ctrl% superad_reduction_gamma_limit
 superad_reduction_gamma_limit_scale = s% ctrl% superad_reduction_gamma_limit_scale
 superad_reduction_gamma_inv_scale = s% ctrl% superad_reduction_gamma_inv_scale
 superad_reduction_diff_grads_limit = s% ctrl% superad_reduction_diff_grads_limit
 superad_reduction_limit = s% ctrl% superad_reduction_limit
 
 max_logT_for_mlt = s% ctrl% max_logT_for_mlt
 mlt_make_surface_no_mixing = s% ctrl% mlt_make_surface_no_mixing
 do_normalize_dqs_as_part_of_set_qs = s% ctrl% do_normalize_dqs_as_part_of_set_qs

 thermohaline_coeff = s% ctrl% thermohaline_coeff
 thermohaline_option = s% ctrl% thermohaline_option
 mixing_length_alpha = s% ctrl% mixing_length_alpha
 remove_small_D_limit = s% ctrl% remove_small_D_limit
 alt_scale_height_flag = s% ctrl% alt_scale_height_flag
 Henyey_MLT_y_param = s% ctrl% Henyey_MLT_y_param
 Henyey_MLT_nu_param = s% ctrl% Henyey_MLT_nu_param
 make_gradr_sticky_in_solver_iters = s% ctrl% make_gradr_sticky_in_solver_iters
 min_logT_for_make_gradr_sticky_in_solver_iters = s% ctrl% min_logT_for_make_gradr_sticky_in_solver_iters
 no_MLT_below_shock = s% ctrl% no_MLT_below_shock
 MLT_option = s% ctrl% MLT_option
 steps_before_use_TDC = s% ctrl% steps_before_use_TDC
 mlt_use_rotation_correction = s% ctrl% mlt_use_rotation_correction
 mlt_Pturb_factor = s% ctrl% mlt_Pturb_factor

 burn_z_mix_region_logT = s% ctrl% burn_z_mix_region_logT
 burn_he_mix_region_logT = s% ctrl% burn_he_mix_region_logT
 burn_h_mix_region_logT = s% ctrl% burn_h_mix_region_logT
 max_Y_for_burn_z_mix_region = s% ctrl% max_Y_for_burn_z_mix_region
 max_X_for_burn_he_mix_region = s% ctrl% max_X_for_burn_he_mix_region
 
 limit_overshoot_Hp_using_size_of_convection_zone = s% ctrl% limit_overshoot_Hp_using_size_of_convection_zone

 predictive_mix = s% ctrl% predictive_mix
 predictive_superad_thresh = s% ctrl% predictive_superad_thresh
 predictive_avoid_reversal = s% ctrl% predictive_avoid_reversal
 predictive_limit_ingestion = s% ctrl% predictive_limit_ingestion
 predictive_ingestion_factor = s% ctrl% predictive_ingestion_factor
 predictive_zone_type = s% ctrl% predictive_zone_type
 predictive_zone_loc = s% ctrl% predictive_zone_loc
 predictive_bdy_loc = s% ctrl% predictive_bdy_loc
 predictive_bdy_q_min = s% ctrl% predictive_bdy_q_min
 predictive_bdy_q_max = s% ctrl% predictive_bdy_q_max

 do_conv_premix = s% ctrl% do_conv_premix
 conv_premix_avoid_increase = s% ctrl% conv_premix_avoid_increase
 conv_premix_time_factor = s% ctrl% conv_premix_time_factor
 conv_premix_fix_pgas = s% ctrl% conv_premix_fix_pgas
 conv_premix_dump_snapshots = s% ctrl% conv_premix_dump_snapshots
 do_premix_heating = s% ctrl% do_premix_heating

 overshoot_f = s% ctrl% overshoot_f
 overshoot_f0 = s% ctrl% overshoot_f0
 overshoot_D0 = s% ctrl% overshoot_D0
 overshoot_Delta0 = s% ctrl% overshoot_Delta0
 overshoot_mass_full_on = s% ctrl% overshoot_mass_full_on
 overshoot_mass_full_off = s% ctrl% overshoot_mass_full_off
 overshoot_scheme = s% ctrl% overshoot_scheme
 overshoot_zone_type = s% ctrl% overshoot_zone_type
 overshoot_zone_loc = s% ctrl% overshoot_zone_loc
 overshoot_bdy_loc = s% ctrl% overshoot_bdy_loc
 overshoot_D_min = s% ctrl% overshoot_D_min
 overshoot_brunt_B_max = s% ctrl% overshoot_brunt_B_max

 max_conv_vel_div_csound = s% ctrl% max_conv_vel_div_csound
 max_v_for_convection = s% ctrl% max_v_for_convection
 max_q_for_convection_with_hydro_on = s% ctrl% max_q_for_convection_with_hydro_on
 max_v_div_cs_for_convection = s% ctrl% max_v_div_cs_for_convection
 max_abs_du_div_cs_for_convection = s% ctrl% max_abs_du_div_cs_for_convection

 calculate_Brunt_B = s% ctrl% calculate_Brunt_B
 calculate_Brunt_N2 = s% ctrl% calculate_Brunt_N2
 brunt_N2_coefficient = s% ctrl% brunt_N2_coefficient
 threshold_for_smooth_brunt_B = s% ctrl% threshold_for_smooth_brunt_B
 min_magnitude_brunt_B = s% ctrl% min_magnitude_brunt_B

 min_overshoot_q = s% ctrl% min_overshoot_q
 overshoot_alpha = s% ctrl% overshoot_alpha
 
   RSP_max_num_periods = s% ctrl% RSP_max_num_periods
   RSP_target_steps_per_cycle = s% ctrl% RSP_target_steps_per_cycle
   RSP_min_max_R_for_periods = s% ctrl% RSP_min_max_R_for_periods
   RSP_min_deltaR_for_periods = s% ctrl% RSP_min_deltaR_for_periods
   RSP_default_PERIODLIN = s% ctrl% RSP_default_PERIODLIN
   RSP_min_PERIOD_div_PERIODLIN = s% ctrl% RSP_min_PERIOD_div_PERIODLIN
   RSP_GREKM_avg_abs_frac_new = s% ctrl% RSP_GREKM_avg_abs_frac_new
   RSP_GREKM_avg_abs_limit = s% ctrl% RSP_GREKM_avg_abs_limit
   RSP_theta = s% ctrl% RSP_theta
   RSP_thetat = s% ctrl% RSP_thetat
   RSP_thetau = s% ctrl% RSP_thetau
   RSP_thetae = s% ctrl% RSP_thetae
   RSP_thetaq = s% ctrl% RSP_thetaq
   RSP_wtr = s% ctrl% RSP_wtr
   RSP_wtc = s% ctrl% RSP_wtc
   RSP_wtt = s% ctrl% RSP_wtt
   RSP_gam = s% ctrl% RSP_gam
   RSP_alfa = s% ctrl% RSP_alfa
   RSP_alfap = s% ctrl% RSP_alfap
   RSP_alfam = s% ctrl% RSP_alfam
   RSP_alfat = s% ctrl% RSP_alfat
   RSP_alfas = s% ctrl% RSP_alfas
   RSP_alfac = s% ctrl% RSP_alfac
   RSP_alfad = s% ctrl% RSP_alfad
   RSP_gammar = s% ctrl% RSP_gammar
   RSP_efl0 = s% ctrl% RSP_efl0
   RSP_min_tau_for_turbulent_flux = s% ctrl% RSP_min_tau_for_turbulent_flux
   RSP_cq = s% ctrl% RSP_cq
   RSP_zsh = s% ctrl% RSP_zsh
   RSP_Qvisc_quadratic = s% ctrl% RSP_Qvisc_quadratic
   RSP_Qvisc_linear = s% ctrl% RSP_Qvisc_linear
   RSP_Qvisc_linear_static = s% ctrl% RSP_Qvisc_linear_static
   RSP_tol_max_corr = s% ctrl% RSP_tol_max_corr
   RSP_tol_max_resid = s% ctrl% RSP_tol_max_resid
   RSP_max_iters_per_try = s% ctrl% RSP_max_iters_per_try
   RSP_max_retries_per_step = s% ctrl% RSP_max_retries_per_step
   RSP_nz_div_IBOTOM = s% ctrl% RSP_nz_div_IBOTOM
   RSP_kick_vsurf_km_per_sec = s% ctrl% RSP_kick_vsurf_km_per_sec
   RSP_fraction_1st_overtone = s% ctrl% RSP_fraction_1st_overtone
   RSP_fraction_2nd_overtone = s% ctrl% RSP_fraction_2nd_overtone
   RSP_Avel = s% ctrl% RSP_Avel
   RSP_Arnd = s% ctrl% RSP_Arnd
   RSP_mode_for_setting_PERIODLIN = s% ctrl% RSP_mode_for_setting_PERIODLIN
   RSP_initial_dt_factor = s% ctrl% RSP_initial_dt_factor
   RSP_v_div_cs_threshold_for_dt_limit = s% ctrl% RSP_v_div_cs_threshold_for_dt_limit
   RSP_max_dt_times_min_dr_div_cs = s% ctrl% RSP_max_dt_times_min_dr_div_cs
   RSP_max_dt_times_min_rad_diff_time = s% ctrl% RSP_max_dt_times_min_rad_diff_time
   RSP_max_dt = s% ctrl% RSP_max_dt
   RSP_testing = s% ctrl% RSP_testing
   RSP_report_limit_dt = s% ctrl% RSP_report_limit_dt
   RSP_use_Prad_for_Psurf = s% ctrl% RSP_use_Prad_for_Psurf
   RSP_report_undercorrections = s% ctrl% RSP_report_undercorrections
   RSP_use_atm_grey_with_kap_for_Psurf = s% ctrl% RSP_use_atm_grey_with_kap_for_Psurf
   use_other_RSP_linear_analysis = s% ctrl% use_other_RSP_linear_analysis
   use_other_RSP_build_model = s% ctrl% use_other_RSP_build_model
   RSP_kap_density_factor = s% ctrl% RSP_kap_density_factor
   RSP_fixed_Psurf = s% ctrl% RSP_fixed_Psurf
   RSP_hydro_only = s% ctrl% RSP_hydro_only
   RSP_tau_surf_for_atm_grey_with_kap = s% ctrl% RSP_tau_surf_for_atm_grey_with_kap
   RSP_Psurf = s% ctrl% RSP_Psurf
   set_RSP_Psurf_to_multiple_of_initial_P1 = s% ctrl% set_RSP_Psurf_to_multiple_of_initial_P1
   RSP_surface_tau = s% ctrl% RSP_surface_tau
   RSP_write_map = s% ctrl% RSP_write_map
   RSP_trace_RSP_build_model = s% ctrl% RSP_trace_RSP_build_model
   RSP_map_filename = s% ctrl% RSP_map_filename
   RSP_map_columns_filename = s% ctrl% RSP_map_columns_filename
   RSP_map_history_filename = s% ctrl% RSP_map_history_filename
   RSP_map_first_period = s% ctrl% RSP_map_first_period
   RSP_map_last_period = s% ctrl% RSP_map_last_period
   RSP_map_zone_interval = s% ctrl% RSP_map_zone_interval
   RSP_nmodes = s% ctrl% RSP_nmodes
   RSP_work_period = s% ctrl% RSP_work_period
   RSP_work_filename = s% ctrl% RSP_work_filename
   RSP_nz_outer = s% ctrl% RSP_nz_outer
   RSP_max_outer_dm_tries = s% ctrl% RSP_max_outer_dm_tries
   RSP_max_inner_scale_tries = s% ctrl% RSP_max_inner_scale_tries
   RSP_relax_max_tries = s% ctrl% RSP_relax_max_tries
   RSP_T_anchor_tolerance = s% ctrl% RSP_T_anchor_tolerance
   RSP_T_inner_tolerance = s% ctrl% RSP_T_inner_tolerance
   RSP_relax_dm_tolerance = s% ctrl% RSP_relax_dm_tolerance
   RSP_dq_1_factor = s% ctrl% RSP_dq_1_factor
   use_RSP_new_start_scheme = s% ctrl% use_RSP_new_start_scheme
   RSP_do_check_omega = s% ctrl% RSP_do_check_omega
   RSP_report_check_omega_changes = s% ctrl% RSP_report_check_omega_changes
   RSP_nz = s% ctrl% RSP_nz
   RSP_T_anchor = s% ctrl% RSP_T_anchor
   RSP_T_inner = s% ctrl% RSP_T_inner
   RSP_relax_initial_model = s% ctrl% RSP_relax_initial_model
   RSP_relax_alfap_before_alfat = s% ctrl% RSP_relax_alfap_before_alfat
   RSP_relax_adjust_inner_mass_distribution = s% ctrl% RSP_relax_adjust_inner_mass_distribution
   RSP_Teff = s% ctrl% RSP_Teff
   RSP_mass = s% ctrl% RSP_mass
   RSP_L = s% ctrl% RSP_L
   RSP_X = s% ctrl% RSP_X
   RSP_Z = s% ctrl% RSP_Z

 RTI_smooth_mass = s% ctrl% RTI_smooth_mass
 RTI_smooth_iterations = s% ctrl% RTI_smooth_iterations
 RTI_smooth_fraction = s% ctrl% RTI_smooth_fraction

 alpha_RTI_diffusion_factor = s% ctrl% alpha_RTI_diffusion_factor
 dudt_RTI_diffusion_factor = s% ctrl% dudt_RTI_diffusion_factor
 dedt_RTI_diffusion_factor = s% ctrl% dedt_RTI_diffusion_factor
 dlnddt_RTI_diffusion_factor = s% ctrl% dlnddt_RTI_diffusion_factor
 composition_RTI_diffusion_factor = s% ctrl% composition_RTI_diffusion_factor
 max_M_RTI_factors_full_on = s% ctrl% max_M_RTI_factors_full_on
 min_M_RTI_factors_full_off = s% ctrl% min_M_RTI_factors_full_off

 alpha_RTI_src_min_v_div_cs = s% ctrl% alpha_RTI_src_min_v_div_cs
 alpha_RTI_src_max_q = s% ctrl% alpha_RTI_src_max_q
 alpha_RTI_src_min_q = s% ctrl% alpha_RTI_src_min_q

 T_mix_limit = s% ctrl% T_mix_limit
 mlt_gradT_fraction = s% ctrl% mlt_gradT_fraction

 ! atmosphere -- surface boundary conditions
 atm_option = s% ctrl% atm_option
 atm_off_table_option = s% ctrl% atm_off_table_option
 Pextra_factor = s% ctrl% Pextra_factor
 atm_fixed_Teff = s% ctrl% atm_fixed_Teff
 atm_fixed_Psurf = s% ctrl% atm_fixed_Psurf
 atm_fixed_Tsurf = s% ctrl% atm_fixed_Tsurf

 atm_T_tau_relation = s% ctrl% atm_T_tau_relation
 atm_T_tau_opacity = s% ctrl% atm_T_tau_opacity
 atm_T_tau_errtol = s% ctrl% atm_T_tau_errtol
 atm_T_tau_max_iters = s% ctrl% atm_T_tau_max_iters
 atm_T_tau_max_steps = s% ctrl% atm_T_tau_max_steps

 atm_table = s% ctrl% atm_table

 atm_irradiated_opacity = s% ctrl% atm_irradiated_opacity
 atm_irradiated_errtol = s% ctrl% atm_irradiated_errtol
 atm_irradiated_T_eq = s% ctrl% atm_irradiated_T_eq
 atm_irradiated_kap_v = s% ctrl% atm_irradiated_kap_v
 atm_irradiated_kap_v_div_kap_th = s% ctrl% atm_irradiated_kap_v_div_kap_th
 atm_irradiated_P_surf = s% ctrl% atm_irradiated_P_surf
 atm_irradiated_max_iters = s% ctrl% atm_irradiated_max_iters

 use_compression_outer_BC = s% ctrl% use_compression_outer_BC
 use_momentum_outer_BC = s% ctrl% use_momentum_outer_BC
 Tsurf_factor = s% ctrl% Tsurf_factor
 use_zero_Pgas_outer_BC = s% ctrl% use_zero_Pgas_outer_BC
 fixed_vsurf = s% ctrl% fixed_vsurf
 use_fixed_vsurf_outer_BC = s% ctrl% use_fixed_vsurf_outer_BC
 fixed_Psurf = s% ctrl% fixed_Psurf
 use_fixed_Psurf_outer_BC = s% ctrl% use_fixed_Psurf_outer_BC

 atm_build_tau_outer = s% ctrl% atm_build_tau_outer
 atm_build_dlogtau = s% ctrl% atm_build_dlogtau
 atm_build_errtol = s% ctrl% atm_build_errtol
 
 use_T_tau_gradr_factor = s% ctrl% use_T_tau_gradr_factor

 ! extra heat near surface to model irradiation
 irradiation_flux = s% ctrl% irradiation_flux
 column_depth_for_irradiation = s% ctrl% column_depth_for_irradiation

 ! extra heat
 inject_uniform_extra_heat = s% ctrl% inject_uniform_extra_heat
 min_q_for_uniform_extra_heat = s% ctrl% min_q_for_uniform_extra_heat
 max_q_for_uniform_extra_heat = s% ctrl% max_q_for_uniform_extra_heat
 inject_extra_ergs_sec = s% ctrl% inject_extra_ergs_sec
 base_of_inject_extra_ergs_sec = s% ctrl% base_of_inject_extra_ergs_sec
 total_mass_for_inject_extra_ergs_sec = s% ctrl% total_mass_for_inject_extra_ergs_sec
 start_time_for_inject_extra_ergs_sec = s% ctrl% start_time_for_inject_extra_ergs_sec
 duration_for_inject_extra_ergs_sec = s% ctrl% duration_for_inject_extra_ergs_sec
 inject_until_reach_model_with_total_energy = s% ctrl% inject_until_reach_model_with_total_energy

 ! mass gain or loss
 mass_change = s% ctrl% mass_change
 mass_change_full_off_dt = s% ctrl% mass_change_full_off_dt
 mass_change_full_on_dt = s% ctrl% mass_change_full_on_dt
 trace_dt_control_mass_change = s% ctrl% trace_dt_control_mass_change
 no_wind_if_no_rotation = s% ctrl% no_wind_if_no_rotation

 min_wind = s% ctrl% min_wind
 max_wind = s% ctrl% max_wind
 use_accreted_material_j = s% ctrl% use_accreted_material_j
 accreted_material_j = s% ctrl% accreted_material_j
 D_omega_mixing_rate = s% ctrl% D_omega_mixing_rate
 D_omega_mixing_across_convection_boundary = s% ctrl% D_omega_mixing_across_convection_boundary
 max_q_for_D_omega_zero_in_convection_region = s% ctrl% max_q_for_D_omega_zero_in_convection_region
 nu_omega_mixing_rate = s% ctrl% nu_omega_mixing_rate
 nu_omega_mixing_across_convection_boundary = s% ctrl% nu_omega_mixing_across_convection_boundary
 max_q_for_nu_omega_zero_in_convection_region = s% ctrl% max_q_for_nu_omega_zero_in_convection_region
 
 mdot_omega_power = s% ctrl% mdot_omega_power
 max_rotational_mdot_boost = s% ctrl% max_rotational_mdot_boost
 max_mdot_jump_for_rotation = s% ctrl% max_mdot_jump_for_rotation
 lim_trace_rotational_mdot_boost = s% ctrl% lim_trace_rotational_mdot_boost
 rotational_mdot_boost_fac = s% ctrl% rotational_mdot_boost_fac
 rotational_mdot_kh_fac = s% ctrl% rotational_mdot_kh_fac
 surf_avg_tau = s% ctrl% surf_avg_tau
 surf_avg_tau_min = s% ctrl% surf_avg_tau_min

 super_eddington_scaling_factor = s% ctrl% super_eddington_scaling_factor
 super_eddington_wind_Ledd_factor = s% ctrl% super_eddington_wind_Ledd_factor
 wind_boost_full_off_L_div_Ledd = s% ctrl% wind_boost_full_off_L_div_Ledd
 wind_boost_full_on_L_div_Ledd = s% ctrl% wind_boost_full_on_L_div_Ledd
 super_eddington_wind_max_boost = s% ctrl% super_eddington_wind_max_boost
 trace_super_eddington_wind_boost = s% ctrl% trace_super_eddington_wind_boost
 
 max_tries_for_implicit_wind = s% ctrl% max_tries_for_implicit_wind
 iwind_tolerance = s% ctrl% iwind_tolerance
 iwind_lambda = s% ctrl% iwind_lambda

 rlo_scaling_factor = s% ctrl% rlo_scaling_factor
 rlo_wind_min_L = s% ctrl% rlo_wind_min_L
 rlo_wind_max_Teff = s% ctrl% rlo_wind_max_Teff
 rlo_wind_roche_lobe_radius = s% ctrl% rlo_wind_roche_lobe_radius
 roche_lobe_xfer_full_on = s% ctrl% roche_lobe_xfer_full_on
 roche_lobe_xfer_full_off = s% ctrl% roche_lobe_xfer_full_off
 rlo_wind_base_mdot = s% ctrl% rlo_wind_base_mdot
 rlo_wind_scale_height = s% ctrl% rlo_wind_scale_height

 cool_wind_RGB_scheme = s% ctrl% cool_wind_RGB_scheme
 cool_wind_AGB_scheme = s% ctrl% cool_wind_AGB_scheme
 RGB_to_AGB_wind_switch = s% ctrl% RGB_to_AGB_wind_switch
 Reimers_scaling_factor = s% ctrl% Reimers_scaling_factor
 Blocker_scaling_factor = s% ctrl% Blocker_scaling_factor
 de_Jager_scaling_factor = s% ctrl% de_Jager_scaling_factor
 van_Loon_scaling_factor = s% ctrl% van_Loon_scaling_factor
 Nieuwenhuijzen_scaling_factor = s% ctrl% Nieuwenhuijzen_scaling_factor
 Vink_scaling_factor = s% ctrl% Vink_scaling_factor
 Dutch_scaling_factor = s% ctrl% Dutch_scaling_factor
 Dutch_wind_lowT_scheme = s% ctrl% Dutch_wind_lowT_scheme

 wind_H_envelope_limit = s% ctrl% wind_H_envelope_limit
 wind_H_He_envelope_limit = s% ctrl% wind_H_He_envelope_limit
 wind_He_layer_limit = s% ctrl% wind_He_layer_limit

 max_logT_for_k_below_const_q = s% ctrl% max_logT_for_k_below_const_q
 max_q_for_k_below_const_q = s% ctrl% max_q_for_k_below_const_q
 min_q_for_k_below_const_q = s% ctrl% min_q_for_k_below_const_q
 max_logT_for_k_const_mass = s% ctrl% max_logT_for_k_const_mass
 min_q_for_k_const_mass = s% ctrl% min_q_for_k_const_mass
 max_q_for_k_const_mass = s% ctrl% max_q_for_k_const_mass

 ! composition of added mass
 accrete_same_as_surface = s% ctrl% accrete_same_as_surface

 accrete_given_mass_fractions = s% ctrl% accrete_given_mass_fractions
 num_accretion_species = s% ctrl% num_accretion_species
 accretion_species_id = s% ctrl% accretion_species_id
 accretion_species_xa = s% ctrl% accretion_species_xa

 accretion_h1 = s% ctrl% accretion_h1
 accretion_h2 = s% ctrl% accretion_h2
 accretion_he3 = s% ctrl% accretion_he3
 accretion_he4 = s% ctrl% accretion_he4
 accretion_zfracs = s% ctrl% accretion_zfracs
 accretion_dump_missing_metals_into_heaviest = s% ctrl% accretion_dump_missing_metals_into_heaviest

 ! special list of z fractions
 z_fraction_li = s% ctrl% z_fraction_li
 z_fraction_be = s% ctrl% z_fraction_be
 z_fraction_b = s% ctrl% z_fraction_b
 z_fraction_c = s% ctrl% z_fraction_c
 z_fraction_n = s% ctrl% z_fraction_n
 z_fraction_o = s% ctrl% z_fraction_o
 z_fraction_f = s% ctrl% z_fraction_f
 z_fraction_ne = s% ctrl% z_fraction_ne
 z_fraction_na = s% ctrl% z_fraction_na
 z_fraction_mg = s% ctrl% z_fraction_mg
 z_fraction_al = s% ctrl% z_fraction_al
 z_fraction_si = s% ctrl% z_fraction_si
 z_fraction_p = s% ctrl% z_fraction_p
 z_fraction_s = s% ctrl% z_fraction_s
 z_fraction_cl = s% ctrl% z_fraction_cl
 z_fraction_ar = s% ctrl% z_fraction_ar
 z_fraction_k = s% ctrl% z_fraction_k
 z_fraction_ca = s% ctrl% z_fraction_ca
 z_fraction_sc = s% ctrl% z_fraction_sc
 z_fraction_ti = s% ctrl% z_fraction_ti
 z_fraction_v = s% ctrl% z_fraction_v
 z_fraction_cr = s% ctrl% z_fraction_cr
 z_fraction_mn = s% ctrl% z_fraction_mn
 z_fraction_fe = s% ctrl% z_fraction_fe
 z_fraction_co = s% ctrl% z_fraction_co
 z_fraction_ni = s% ctrl% z_fraction_ni
 z_fraction_cu = s% ctrl% z_fraction_cu
 z_fraction_zn = s% ctrl% z_fraction_zn

 lgT_lo_for_set_new_abundances = s% ctrl% lgT_lo_for_set_new_abundances
 lgT_hi_for_set_new_abundances = s% ctrl% lgT_hi_for_set_new_abundances

 ! automatic stops for mass loss/gain
 max_star_mass_for_gain = s% ctrl% max_star_mass_for_gain
 min_star_mass_for_loss = s% ctrl% min_star_mass_for_loss
 max_T_center_for_any_mass_loss = s% ctrl% max_T_center_for_any_mass_loss
 max_T_center_for_full_mass_loss = s% ctrl% max_T_center_for_full_mass_loss

 ! relaxation parameters
 extra_power_source = s% ctrl% extra_power_source
 relax_dlnZ = s% ctrl% relax_dlnZ
 relax_dY = s% ctrl% relax_dY

 ! mesh adjustment
 show_mesh_changes = s% ctrl% show_mesh_changes
 okay_to_remesh = s% ctrl% okay_to_remesh
 restore_mesh_on_retry = s% ctrl% restore_mesh_on_retry
 num_steps_to_hold_mesh_after_retry = s% ctrl% num_steps_to_hold_mesh_after_retry
 trace_mesh_adjust_error_in_conservation = s% ctrl% trace_mesh_adjust_error_in_conservation
 max_rel_delta_IE_for_mesh_total_energy_balance = s% ctrl% max_rel_delta_IE_for_mesh_total_energy_balance
 max_allowed_nz = s% ctrl% max_allowed_nz
 mesh_max_allowed_ratio = s% ctrl% mesh_max_allowed_ratio
 remesh_max_allowed_logT = s% ctrl% remesh_max_allowed_logT
 max_delta_x_for_merge = s% ctrl% max_delta_x_for_merge

 mesh_ok_to_merge = s% ctrl% mesh_ok_to_merge
 mesh_max_k_old_for_split = s% ctrl% mesh_max_k_old_for_split
 mesh_min_k_old_for_split = s% ctrl% mesh_min_k_old_for_split
 mesh_adjust_get_T_from_E = s% ctrl% mesh_adjust_get_T_from_E

 max_dq = s% ctrl% max_dq
 min_dq = s% ctrl% min_dq
 min_dq_for_split = s% ctrl% min_dq_for_split
 min_dq_for_xa = s% ctrl% min_dq_for_xa
 min_dq_for_xa_convective = s% ctrl% min_dq_for_xa_convective
 min_dq_for_logT = s% ctrl% min_dq_for_logT

 mesh_min_dlnR = s% ctrl% mesh_min_dlnR
 merge_if_dlnR_too_small = s% ctrl% merge_if_dlnR_too_small

 mesh_min_dr_div_dRstar = s% ctrl% mesh_min_dr_div_dRstar
 merge_if_dr_div_dRstar_too_small = s% ctrl% merge_if_dr_div_dRstar_too_small

 mesh_min_dr_div_cs = s% ctrl% mesh_min_dr_div_cs
 merge_if_dr_div_cs_too_small = s% ctrl% merge_if_dr_div_cs_too_small

 max_center_cell_dq = s% ctrl% max_center_cell_dq
 max_surface_cell_dq = s% ctrl% max_surface_cell_dq
 max_num_subcells = s% ctrl% max_num_subcells
 max_num_merge_cells = s% ctrl% max_num_merge_cells

 mesh_delta_coeff = s% ctrl% mesh_delta_coeff
 mesh_delta_coeff_for_highT = s% ctrl% mesh_delta_coeff_for_highT
 logT_max_for_standard_mesh_delta_coeff = s% ctrl% logT_max_for_standard_mesh_delta_coeff
 logT_min_for_highT_mesh_delta_coeff = s% ctrl% logT_min_for_highT_mesh_delta_coeff
 mesh_Pgas_div_P_exponent = s% ctrl% mesh_Pgas_div_P_exponent

 remesh_dt_limit = s% ctrl% remesh_dt_limit

 E_function_weight = s% ctrl% E_function_weight
 E_function_param = s% ctrl% E_function_param
 P_function_weight = s% ctrl% P_function_weight

 mesh_logX_species = s% ctrl% mesh_logX_species
 mesh_logX_min_for_extra = s% ctrl% mesh_logX_min_for_extra
 mesh_dlogX_dlogP_extra = s% ctrl% mesh_dlogX_dlogP_extra
 mesh_dlogX_dlogP_full_on = s% ctrl% mesh_dlogX_dlogP_full_on
 mesh_dlogX_dlogP_full_off = s% ctrl% mesh_dlogX_dlogP_full_off

 mesh_dlog_eps_min_for_extra = s% ctrl% mesh_dlog_eps_min_for_extra
 mesh_dlog_eps_dlogP_full_on = s% ctrl% mesh_dlog_eps_dlogP_full_on
 mesh_dlog_eps_dlogP_full_off = s% ctrl% mesh_dlog_eps_dlogP_full_off

 mesh_dlog_pp_dlogP_extra = s% ctrl% mesh_dlog_pp_dlogP_extra
 mesh_dlog_cno_dlogP_extra = s% ctrl% mesh_dlog_cno_dlogP_extra
 mesh_dlog_3alf_dlogP_extra = s% ctrl% mesh_dlog_3alf_dlogP_extra

 mesh_dlog_burn_c_dlogP_extra = s% ctrl% mesh_dlog_burn_c_dlogP_extra
 mesh_dlog_burn_n_dlogP_extra = s% ctrl% mesh_dlog_burn_n_dlogP_extra
 mesh_dlog_burn_o_dlogP_extra = s% ctrl% mesh_dlog_burn_o_dlogP_extra
 mesh_dlog_burn_ne_dlogP_extra = s% ctrl% mesh_dlog_burn_ne_dlogP_extra
 mesh_dlog_burn_na_dlogP_extra = s% ctrl% mesh_dlog_burn_na_dlogP_extra
 mesh_dlog_burn_mg_dlogP_extra = s% ctrl% mesh_dlog_burn_mg_dlogP_extra
 mesh_dlog_burn_si_dlogP_extra = s% ctrl% mesh_dlog_burn_si_dlogP_extra
 mesh_dlog_burn_s_dlogP_extra = s% ctrl% mesh_dlog_burn_s_dlogP_extra
 mesh_dlog_burn_ar_dlogP_extra = s% ctrl% mesh_dlog_burn_ar_dlogP_extra
 mesh_dlog_burn_ca_dlogP_extra = s% ctrl% mesh_dlog_burn_ca_dlogP_extra
 mesh_dlog_burn_ti_dlogP_extra = s% ctrl% mesh_dlog_burn_ti_dlogP_extra
 mesh_dlog_burn_cr_dlogP_extra = s% ctrl% mesh_dlog_burn_cr_dlogP_extra
 mesh_dlog_burn_fe_dlogP_extra = s% ctrl% mesh_dlog_burn_fe_dlogP_extra

 mesh_dlog_cc_dlogP_extra = s% ctrl% mesh_dlog_cc_dlogP_extra
 mesh_dlog_co_dlogP_extra = s% ctrl% mesh_dlog_co_dlogP_extra
 mesh_dlog_oo_dlogP_extra = s% ctrl% mesh_dlog_oo_dlogP_extra

 mesh_dlog_pnhe4_dlogP_extra = s% ctrl% mesh_dlog_pnhe4_dlogP_extra
 mesh_dlog_photo_dlogP_extra = s% ctrl% mesh_dlog_photo_dlogP_extra
 mesh_dlog_other_dlogP_extra = s% ctrl% mesh_dlog_other_dlogP_extra
 
 mesh_delta_coeff_factor_smooth_iters = s% ctrl% mesh_delta_coeff_factor_smooth_iters

 T_function1_weight = s% ctrl% T_function1_weight
 T_function2_weight = s% ctrl% T_function2_weight
 T_function2_param = s% ctrl% T_function2_param

 R_function_weight = s% ctrl% R_function_weight
 R_function_param = s% ctrl% R_function_param

 R_function2_weight = s% ctrl% R_function2_weight
 R_function2_param1 = s% ctrl% R_function2_param1
 R_function2_param2 = s% ctrl% R_function2_param2

 R_function3_weight = s% ctrl% R_function3_weight

 M_function_weight = s% ctrl% M_function_weight
 M_function_param = s% ctrl% M_function_param

 gradT_function_weight = s% ctrl% gradT_function_weight
 log_tau_function_weight = s% ctrl% log_tau_function_weight
 log_kap_function_weight = s% ctrl% log_kap_function_weight
 omega_function_weight = s% ctrl% omega_function_weight

 gam_function_weight = s% ctrl% gam_function_weight
 gam_function_param1 = s% ctrl% gam_function_param1
 gam_function_param2 = s% ctrl% gam_function_param2

 xa_function_species = s% ctrl% xa_function_species
 xa_function_weight = s% ctrl% xa_function_weight
 xa_function_param = s% ctrl% xa_function_param
 xa_mesh_delta_coeff = s% ctrl% xa_mesh_delta_coeff
 
 use_split_merge_amr = s% ctrl% use_split_merge_amr
 split_merge_amr_nz_baseline = s% ctrl% split_merge_amr_nz_baseline
 split_merge_amr_nz_r_core = s% ctrl% split_merge_amr_nz_r_core
 split_merge_amr_nz_r_core_fraction = s% ctrl% split_merge_amr_nz_r_core_fraction
 split_merge_amr_mesh_delta_coeff = s% ctrl% split_merge_amr_mesh_delta_coeff
 split_merge_amr_log_zoning = s% ctrl% split_merge_amr_log_zoning
 split_merge_amr_hybrid_zoning = s% ctrl% split_merge_amr_hybrid_zoning
 split_merge_amr_flipped_hybrid_zoning = s% ctrl% split_merge_amr_flipped_hybrid_zoning
 split_merge_amr_logtau_zoning = s% ctrl% split_merge_amr_logtau_zoning
 split_merge_amr_okay_to_split_nz = s% ctrl% split_merge_amr_okay_to_split_nz
 split_merge_amr_okay_to_split_1 = s% ctrl% split_merge_amr_okay_to_split_1
 merge_amr_inhibit_at_jumps = s% ctrl% merge_amr_inhibit_at_jumps
 split_merge_amr_MaxLong = s% ctrl% split_merge_amr_MaxLong
 split_merge_amr_MaxShort = s% ctrl% split_merge_amr_MaxShort
 merge_amr_max_abs_du_div_cs = s% ctrl% merge_amr_max_abs_du_div_cs
 merge_amr_ignore_surface_cells = s% ctrl% merge_amr_ignore_surface_cells
 merge_amr_du_div_cs_limit_only_for_compression = s% ctrl% merge_amr_du_div_cs_limit_only_for_compression
 split_merge_amr_avoid_repeated_remesh = s% ctrl% split_merge_amr_avoid_repeated_remesh
 merge_amr_k_for_ignore_surface_cells = s% ctrl% merge_amr_k_for_ignore_surface_cells
 split_merge_amr_dq_min = s% ctrl% split_merge_amr_dq_min
 split_merge_amr_dq_max = s% ctrl% split_merge_amr_dq_max
 split_merge_amr_r_core_cm = s% ctrl% split_merge_amr_r_core_cm
 split_merge_amr_max_iters = s% ctrl% split_merge_amr_max_iters
 trace_split_merge_amr = s% ctrl% trace_split_merge_amr
 equal_split_density_amr = s% ctrl% equal_split_density_amr

 ! nuclear reaction parameters
 screening_mode = s% ctrl% screening_mode
 default_net_name = s% ctrl% default_net_name

 net_logTcut_lo = s% ctrl% net_logTcut_lo
 net_logTcut_lim = s% ctrl% net_logTcut_lim

 eps_nuc_factor = s% ctrl% eps_nuc_factor
 op_split_burn_eps_nuc_infall_limit = s% ctrl% op_split_burn_eps_nuc_infall_limit
 eps_WD_sedimentation_factor = s% ctrl% eps_WD_sedimentation_factor
 max_abs_eps_nuc = s% ctrl% max_abs_eps_nuc
 dxdt_nuc_factor = s% ctrl% dxdt_nuc_factor
 max_abar_for_burning = s% ctrl% max_abar_for_burning
 fe56ec_fake_factor = s% ctrl% fe56ec_fake_factor
 min_T_for_fe56ec_fake_factor = s% ctrl% min_T_for_fe56ec_fake_factor
 weak_rate_factor = s% ctrl% weak_rate_factor

 mix_factor = s% ctrl% mix_factor

 sig_term_limit = s% ctrl% sig_term_limit

 sig_min_factor_for_high_Tcenter = s% ctrl% sig_min_factor_for_high_Tcenter
 Tcenter_min_for_sig_min_factor_full_on = s% ctrl% Tcenter_min_for_sig_min_factor_full_on
 Tcenter_max_for_sig_min_factor_full_off = s% ctrl% Tcenter_max_for_sig_min_factor_full_off
 max_delta_m_to_bdy_for_sig_min_factor = s% ctrl% max_delta_m_to_bdy_for_sig_min_factor
 delta_m_lower_for_sig_min_factor = s% ctrl% delta_m_lower_for_sig_min_factor
 delta_m_upper_for_sig_min_factor = s% ctrl% delta_m_upper_for_sig_min_factor

 am_sig_term_limit = s% ctrl% am_sig_term_limit
 am_D_mix_factor = s% ctrl% am_D_mix_factor
 am_gradmu_factor = s% ctrl% am_gradmu_factor
 am_nu_factor = s% ctrl% am_nu_factor

 D_visc_factor = s% ctrl% D_visc_factor
 D_DSI_factor = s% ctrl% D_DSI_factor
 D_SH_factor = s% ctrl% D_SH_factor
 D_SSI_factor = s% ctrl% D_SSI_factor
 D_ES_factor = s% ctrl% D_ES_factor
 D_GSF_factor = s% ctrl% D_GSF_factor
 D_ST_factor = s% ctrl% D_ST_factor

 am_nu_non_rotation_factor = s% ctrl% am_nu_non_rotation_factor
 skip_rotation_in_convection_zones = s% ctrl% skip_rotation_in_convection_zones
 am_nu_DSI_factor = s% ctrl% am_nu_DSI_factor
 am_nu_SH_factor = s% ctrl% am_nu_SH_factor
 am_nu_SSI_factor = s% ctrl% am_nu_SSI_factor
 am_nu_ES_factor = s% ctrl% am_nu_ES_factor
 am_nu_GSF_factor = s% ctrl% am_nu_GSF_factor
 am_nu_ST_factor = s% ctrl% am_nu_ST_factor
 am_nu_visc_factor = s% ctrl% am_nu_visc_factor

 am_nu_omega_rot_factor = s% ctrl% am_nu_omega_rot_factor
 am_nu_omega_non_rot_factor = s% ctrl% am_nu_omega_non_rot_factor
 am_nu_j_rot_factor = s% ctrl% am_nu_j_rot_factor
 am_nu_j_non_rot_factor = s% ctrl% am_nu_j_non_rot_factor

 smooth_nu_ST = s% ctrl% smooth_nu_ST
 smooth_D_ST = s% ctrl% smooth_D_ST
 smooth_D_DSI = s% ctrl% smooth_D_DSI
 smooth_D_SSI = s% ctrl% smooth_D_SSI
 smooth_D_SH = s% ctrl% smooth_D_SH
 smooth_D_GSF = s% ctrl% smooth_D_GSF
 smooth_D_ES = s% ctrl% smooth_D_ES
 smooth_D_omega = s% ctrl% smooth_D_omega
 smooth_am_nu_rot = s% ctrl% smooth_am_nu_rot
 ST_angsmt = s% ctrl% ST_angsmt
 ST_angsml = s% ctrl% ST_angsml

 simple_i_rot_flag = s% ctrl% simple_i_rot_flag
 do_adjust_J_lost = s% ctrl% do_adjust_J_lost
 premix_omega = s% ctrl% premix_omega
 angular_momentum_error_warn = s% ctrl% angular_momentum_error_warn
 angular_momentum_error_retry = s% ctrl% angular_momentum_error_retry
 recalc_mixing_info_each_substep = s% ctrl% recalc_mixing_info_each_substep
 adjust_J_fraction = s% ctrl% adjust_J_fraction
 min_q_for_adjust_J_lost = s% ctrl% min_q_for_adjust_J_lost
 min_J_div_delta_J = s% ctrl% min_J_div_delta_J
 max_mdot_redo_cnt = s% ctrl% max_mdot_redo_cnt
 mdot_revise_factor = s% ctrl% mdot_revise_factor
 implicit_mdot_boost = s% ctrl% implicit_mdot_boost
 min_years_dt_for_redo_mdot = s% ctrl% min_years_dt_for_redo_mdot
 surf_omega_div_omega_crit_limit = s% ctrl% surf_omega_div_omega_crit_limit
 surf_omega_div_omega_crit_tol = s% ctrl% surf_omega_div_omega_crit_tol
 w_div_wcrit_max = s% ctrl% w_div_wcrit_max
 w_div_wcrit_max2 = s% ctrl% w_div_wcrit_max2
 fp_min = s% ctrl% fp_min
 ft_min = s% ctrl% ft_min
 fp_error_limit = s% ctrl% fp_error_limit
 ft_error_limit = s% ctrl% ft_error_limit

 D_mix_rotation_max_logT_full_on = s% ctrl% D_mix_rotation_max_logT_full_on
 D_mix_rotation_min_logT_full_off = s% ctrl% D_mix_rotation_min_logT_full_off

 set_uniform_am_nu_non_rot = s% ctrl% set_uniform_am_nu_non_rot
 uniform_am_nu_non_rot = s% ctrl% uniform_am_nu_non_rot

 set_min_am_nu_non_rot = s% ctrl% set_min_am_nu_non_rot
 min_am_nu_non_rot = s% ctrl% min_am_nu_non_rot
 min_center_Ye_for_min_am_nu_non_rot = s% ctrl% min_center_Ye_for_min_am_nu_non_rot

 set_min_D_mix = s% ctrl% set_min_D_mix
 mass_lower_limit_for_min_D_mix = s% ctrl% mass_lower_limit_for_min_D_mix
 mass_upper_limit_for_min_D_mix = s% ctrl% mass_upper_limit_for_min_D_mix
 min_D_mix = s% ctrl% min_D_mix
 set_min_D_mix_below_Tmax = s% ctrl% set_min_D_mix_below_Tmax
 min_D_mix_below_Tmax = s% ctrl% min_D_mix_below_Tmax
 set_min_D_mix_in_H_He = s% ctrl% set_min_D_mix_in_H_He
 min_D_mix_in_H_He = s% ctrl% min_D_mix_in_H_He
 min_center_Ye_for_min_D_mix = s% ctrl% min_center_Ye_for_min_D_mix
 reaction_neuQs_factor = s% ctrl% reaction_neuQs_factor
 nonlocal_NiCo_kap_gamma = s% ctrl% nonlocal_NiCo_kap_gamma
 nonlocal_NiCo_decay_heat = s% ctrl% nonlocal_NiCo_decay_heat
 dtau_gamma_NiCo_decay_heat = s% ctrl% dtau_gamma_NiCo_decay_heat
 max_logT_for_net = s% ctrl% max_logT_for_net
 smooth_outer_xa_big = s% ctrl% smooth_outer_xa_big
 smooth_outer_xa_small = s% ctrl% smooth_outer_xa_small

 ! element diffusion parameters
 diffusion_use_iben_macdonald = s% ctrl% diffusion_use_iben_macdonald
 diffusion_use_paquette = s% ctrl% diffusion_use_paquette
 diffusion_use_cgs_solver = s% ctrl% diffusion_use_cgs_solver
 diffusion_use_full_net = s% ctrl% diffusion_use_full_net
 do_WD_sedimentation_heating = s% ctrl% do_WD_sedimentation_heating
 min_xa_for_WD_sedimentation_heating = s% ctrl% min_xa_for_WD_sedimentation_heating
 do_diffusion_heating = s% ctrl% do_diffusion_heating
 do_element_diffusion = s% ctrl% do_element_diffusion
 cgs_thermal_diffusion_eta_full_on = s% ctrl% cgs_thermal_diffusion_eta_full_on
 cgs_thermal_diffusion_eta_full_off = s% ctrl% cgs_thermal_diffusion_eta_full_off
 diffusion_min_dq_at_surface = s% ctrl% diffusion_min_dq_at_surface
 diffusion_min_T_at_surface = s% ctrl% diffusion_min_T_at_surface
 diffusion_min_dq_ratio_at_surface = s% ctrl% diffusion_min_dq_ratio_at_surface
 diffusion_dt_limit = s% ctrl% diffusion_dt_limit

 do_phase_separation = s% ctrl% do_phase_separation
 do_phase_separation_heating = s% ctrl% do_phase_separation_heating
 phase_separation_mixing_use_brunt = s% ctrl% phase_separation_mixing_use_brunt
 phase_separation_no_diffusion = s% ctrl% phase_separation_no_diffusion

 diffusion_min_X_hard_limit = s% ctrl% diffusion_min_X_hard_limit
 diffusion_X_total_atol = s% ctrl% diffusion_X_total_atol
 diffusion_X_total_rtol = s% ctrl% diffusion_X_total_rtol
 diffusion_upwind_abs_v_limit = s% ctrl% diffusion_upwind_abs_v_limit
 diffusion_dt_div_timescale = s% ctrl% diffusion_dt_div_timescale
 diffusion_min_num_substeps = s% ctrl% diffusion_min_num_substeps
 diffusion_max_iters_per_substep = s% ctrl% diffusion_max_iters_per_substep
 diffusion_max_retries_per_substep = s% ctrl% diffusion_max_retries_per_substep
 diffusion_v_max = s% ctrl% diffusion_v_max
 diffusion_gamma_full_off = s% ctrl% diffusion_gamma_full_off
 diffusion_gamma_full_on = s% ctrl% diffusion_gamma_full_on
 diffusion_T_full_off = s% ctrl% diffusion_T_full_off
 D_mix_ignore_diffusion = s% ctrl% D_mix_ignore_diffusion
 diffusion_T_full_on = s% ctrl% diffusion_T_full_on
 diffusion_calculates_ionization = s% ctrl% diffusion_calculates_ionization
 diffusion_nsmooth_typical_charge = s% ctrl% diffusion_nsmooth_typical_charge
 diffusion_tol_correction_max = s% ctrl% diffusion_tol_correction_max
 diffusion_tol_correction_norm = s% ctrl% diffusion_tol_correction_norm

 diffusion_AD_dm_full_on = s% ctrl% diffusion_AD_dm_full_on
 diffusion_AD_dm_full_off = s% ctrl% diffusion_AD_dm_full_off
 diffusion_AD_boost_factor = s% ctrl% diffusion_AD_boost_factor

 diffusion_SIG_factor = s% ctrl% diffusion_SIG_factor
 diffusion_GT_factor = s% ctrl% diffusion_GT_factor

 diffusion_Vlimit_dm_full_on = s% ctrl% diffusion_Vlimit_dm_full_on
 diffusion_Vlimit_dm_full_off = s% ctrl% diffusion_Vlimit_dm_full_off
 diffusion_Vlimit = s% ctrl% diffusion_Vlimit

 diffusion_max_T_for_radaccel = s% ctrl% diffusion_max_T_for_radaccel
 diffusion_min_T_for_radaccel = s% ctrl% diffusion_min_T_for_radaccel
 diffusion_max_Z_for_radaccel = s% ctrl% diffusion_max_Z_for_radaccel
 diffusion_min_Z_for_radaccel = s% ctrl% diffusion_min_Z_for_radaccel
 diffusion_screening_for_radaccel = s% ctrl% diffusion_screening_for_radaccel
 op_mono_data_path = s% ctrl% op_mono_data_path
 op_mono_data_cache_filename = s% ctrl% op_mono_data_cache_filename

 show_diffusion_info = s% ctrl% show_diffusion_info
 show_diffusion_substep_info = s% ctrl% show_diffusion_substep_info
 show_diffusion_timing = s% ctrl% show_diffusion_timing

 diffusion_num_classes = s% ctrl% diffusion_num_classes
 diffusion_class_representative = s% ctrl% diffusion_class_representative
 diffusion_class_A_max = s% ctrl% diffusion_class_A_max
 diffusion_class_typical_charge = s% ctrl% diffusion_class_typical_charge
 diffusion_class_factor = s% ctrl% diffusion_class_factor

 diffusion_use_isolve = s% ctrl% diffusion_use_isolve
 diffusion_rtol_for_isolve = s% ctrl% diffusion_rtol_for_isolve
 diffusion_atol_for_isolve = s% ctrl% diffusion_atol_for_isolve
 diffusion_maxsteps_for_isolve = s% ctrl% diffusion_maxsteps_for_isolve
 diffusion_isolve_solver = s% ctrl% diffusion_isolve_solver

 ! eos controls
 fix_d_eos_dxa_partials = s% ctrl% fix_d_eos_dxa_partials
 
 ! opacity controls
 use_simple_es_for_kap = s% ctrl% use_simple_es_for_kap
 use_starting_composition_for_kap = s% ctrl% use_starting_composition_for_kap
 min_kap_for_dPrad_dm_eqn = s% ctrl% min_kap_for_dPrad_dm_eqn

 low_logT_op_mono_full_off = s% ctrl% low_logT_op_mono_full_off
 low_logT_op_mono_full_on = s% ctrl% low_logT_op_mono_full_on
 high_logT_op_mono_full_off = s% ctrl% high_logT_op_mono_full_off
 high_logT_op_mono_full_on = s% ctrl% high_logT_op_mono_full_on
 op_mono_min_X_to_include = s% ctrl% op_mono_min_X_to_include
 use_op_mono_alt_get_kap = s% ctrl% use_op_mono_alt_get_kap

 include_L_in_correction_limits = s% ctrl% include_L_in_correction_limits
 include_v_in_correction_limits = s% ctrl% include_v_in_correction_limits
 include_u_in_correction_limits = s% ctrl% include_u_in_correction_limits
 include_w_in_correction_limits = s% ctrl% include_w_in_correction_limits

 ! asteroseismology controls

 get_delta_nu_from_scaled_solar = s% ctrl% get_delta_nu_from_scaled_solar
 nu_max_sun = s% ctrl% nu_max_sun
 delta_nu_sun = s% ctrl% delta_nu_sun
 Teff_sun = s% ctrl% Teff_sun
 delta_Pg_mode_freq = s% ctrl% delta_Pg_mode_freq

 ! hydro parameters
 energy_eqn_option = s% ctrl% energy_eqn_option
 opacity_max = s% ctrl% opacity_max
 opacity_factor = s% ctrl% opacity_factor
 min_logT_for_opacity_factor_off = s% ctrl% min_logT_for_opacity_factor_off
 min_logT_for_opacity_factor_on = s% ctrl% min_logT_for_opacity_factor_on
 max_logT_for_opacity_factor_on = s% ctrl% max_logT_for_opacity_factor_on
 max_logT_for_opacity_factor_off = s% ctrl% max_logT_for_opacity_factor_off

 non_nuc_neu_factor = s% ctrl% non_nuc_neu_factor
 use_time_centered_eps_grav = s% ctrl% use_time_centered_eps_grav
 no_dedt_form_during_relax = s% ctrl% no_dedt_form_during_relax
 dedt_eqn_r_scale = s% ctrl% dedt_eqn_r_scale
 use_mass_corrections = s% ctrl% use_mass_corrections
 use_gravity_rotation_correction = s% ctrl% use_gravity_rotation_correction
 eps_grav_factor = s% ctrl% eps_grav_factor
 eps_mdot_factor = s% ctrl% eps_mdot_factor
 include_composition_in_eps_grav = s% ctrl% include_composition_in_eps_grav
 max_abs_rel_change_surf_lnS = s% ctrl% max_abs_rel_change_surf_lnS
 max_num_surf_revisions = s% ctrl% max_num_surf_revisions
 Gamma_lnS_eps_grav_full_off = s% ctrl% Gamma_lnS_eps_grav_full_off
 Gamma_lnS_eps_grav_full_on = s% ctrl% Gamma_lnS_eps_grav_full_on

 use_dPrad_dm_form_of_T_gradient_eqn = s% ctrl% use_dPrad_dm_form_of_T_gradient_eqn
 use_gradT_actual_vs_gradT_MLT_for_T_gradient_eqn = s% ctrl% use_gradT_actual_vs_gradT_MLT_for_T_gradient_eqn
 steps_before_use_velocity_time_centering = s% ctrl% steps_before_use_velocity_time_centering
 include_P_in_velocity_time_centering = s% ctrl% include_P_in_velocity_time_centering
 include_L_in_velocity_time_centering = s% ctrl% include_L_in_velocity_time_centering
 P_theta_for_velocity_time_centering = s% ctrl% P_theta_for_velocity_time_centering
 L_theta_for_velocity_time_centering = s% ctrl% L_theta_for_velocity_time_centering
 use_P_d_1_div_rho_form_of_work_when_time_centering_velocity = s% ctrl% use_P_d_1_div_rho_form_of_work_when_time_centering_velocity

 RTI_A = s% ctrl% RTI_A
 RTI_B = s% ctrl% RTI_B
 RTI_C = s% ctrl% RTI_C
 RTI_D = s% ctrl% RTI_D
 RTI_max_alpha = s% ctrl% RTI_max_alpha
 RTI_C_X_factor = s% ctrl% RTI_C_X_factor
 RTI_C_X0_frac = s% ctrl% RTI_C_X0_frac
 RTI_dm_for_center_eta_nondecreasing = s% ctrl% RTI_dm_for_center_eta_nondecreasing
 RTI_min_dm_behind_shock_for_full_on = s% ctrl% RTI_min_dm_behind_shock_for_full_on
 RTI_energy_floor = s% ctrl% RTI_energy_floor
 RTI_D_mix_floor = s% ctrl% RTI_D_mix_floor
 RTI_min_m_for_D_mix_floor = s% ctrl% RTI_min_m_for_D_mix_floor
 RTI_log_max_boost = s% ctrl% RTI_log_max_boost 
 RTI_m_full_boost = s% ctrl% RTI_m_full_boost
 RTI_m_no_boost = s% ctrl% RTI_m_no_boost

 velocity_logT_lower_bound = s% ctrl% velocity_logT_lower_bound
 max_dt_yrs_for_velocity_logT_lower_bound = s% ctrl% max_dt_yrs_for_velocity_logT_lower_bound
 velocity_q_upper_bound = s% ctrl% velocity_q_upper_bound

 retry_for_v_above_clight = s% ctrl% retry_for_v_above_clight

 ! solvers

 tol_correction_norm = s% ctrl% tol_correction_norm
 tol_max_correction = s% ctrl% tol_max_correction
 correction_xa_limit = s% ctrl% correction_xa_limit

 tol_correction_high_T_limit = s% ctrl% tol_correction_high_T_limit
 tol_correction_norm_high_T = s% ctrl% tol_correction_norm_high_T
 tol_max_correction_high_T = s% ctrl% tol_max_correction_high_T

 tol_correction_extreme_T_limit = s% ctrl% tol_correction_extreme_T_limit
 tol_correction_norm_extreme_T = s% ctrl% tol_correction_norm_extreme_T
 tol_max_correction_extreme_T = s% ctrl% tol_max_correction_extreme_T
 
 tol_bad_max_correction = s% ctrl% tol_bad_max_correction
 bad_max_correction_series_limit = s% ctrl% bad_max_correction_series_limit

 tol_residual_norm1 = s% ctrl% tol_residual_norm1
 tol_max_residual1 = s% ctrl% tol_max_residual1
 tol_residual_norm2 = s% ctrl% tol_residual_norm2
 tol_max_residual2 = s% ctrl% tol_max_residual2
 tol_residual_norm3 = s% ctrl% tol_residual_norm3
 tol_max_residual3 = s% ctrl% tol_max_residual3
 warning_limit_for_max_residual = s% ctrl% warning_limit_for_max_residual
 trace_solver_damping = s% ctrl% trace_solver_damping
 
 relax_use_gold_tolerances = s% ctrl% relax_use_gold_tolerances
 relax_tol_correction_norm = s% ctrl% relax_tol_correction_norm
 relax_tol_max_correction = s% ctrl% relax_tol_max_correction
 relax_solver_iters_timestep_limit = s% ctrl% relax_solver_iters_timestep_limit
 relax_iter_for_resid_tol2 = s% ctrl% relax_iter_for_resid_tol2
 relax_tol_residual_norm1 = s% ctrl% relax_tol_residual_norm1
 relax_tol_max_residual1 = s% ctrl% relax_tol_max_residual1
 relax_iter_for_resid_tol3 = s% ctrl% relax_iter_for_resid_tol3
 relax_tol_residual_norm2 = s% ctrl% relax_tol_residual_norm2
 relax_tol_max_residual2 = s% ctrl% relax_tol_max_residual2
 relax_tol_residual_norm3 = s% ctrl% relax_tol_residual_norm3
 relax_tol_max_residual3 = s% ctrl% relax_tol_max_residual3
 relax_maxT_for_gold_tolerances = s% ctrl% relax_maxT_for_gold_tolerances
 
 use_gold_tolerances = s% ctrl% use_gold_tolerances
 gold_solver_iters_timestep_limit = s% ctrl% gold_solver_iters_timestep_limit 
 maxT_for_gold_tolerances = s% ctrl% maxT_for_gold_tolerances
 gold_tol_residual_norm1 = s% ctrl% gold_tol_residual_norm1
 gold_tol_max_residual1 = s% ctrl% gold_tol_max_residual1
 gold_iter_for_resid_tol2 = s% ctrl% gold_iter_for_resid_tol2
 gold_tol_residual_norm2 = s% ctrl% gold_tol_residual_norm2
 gold_tol_max_residual2 = s% ctrl% gold_tol_max_residual2
 gold_iter_for_resid_tol3 = s% ctrl% gold_iter_for_resid_tol3
 gold_tol_residual_norm3 = s% ctrl% gold_tol_residual_norm3
 gold_tol_max_residual3 = s% ctrl% gold_tol_max_residual3
 steps_before_use_gold_tolerances = s% ctrl% steps_before_use_gold_tolerances
 
 use_gold2_tolerances = s% ctrl% use_gold2_tolerances
 gold2_solver_iters_timestep_limit = s% ctrl% gold2_solver_iters_timestep_limit 
 gold2_tol_residual_norm1 = s% ctrl% gold2_tol_residual_norm1
 gold2_tol_max_residual1 = s% ctrl% gold2_tol_max_residual1
 gold2_iter_for_resid_tol2 = s% ctrl% gold2_iter_for_resid_tol2
 gold2_tol_residual_norm2 = s% ctrl% gold2_tol_residual_norm2
 gold2_tol_max_residual2 = s% ctrl% gold2_tol_max_residual2
 gold2_iter_for_resid_tol3 = s% ctrl% gold2_iter_for_resid_tol3
 gold2_tol_residual_norm3 = s% ctrl% gold2_tol_residual_norm3
 gold2_tol_max_residual3 = s% ctrl% gold2_tol_max_residual3
 steps_before_use_gold2_tolerances = s% ctrl% steps_before_use_gold2_tolerances
 
 include_rotation_in_total_energy = s% ctrl% include_rotation_in_total_energy

 convergence_ignore_equL_residuals = s% ctrl% convergence_ignore_equL_residuals
 convergence_ignore_alpha_RTI_residuals = s% ctrl% convergence_ignore_alpha_RTI_residuals

 iter_for_resid_tol2 = s% ctrl% iter_for_resid_tol2
 iter_for_resid_tol3 = s% ctrl% iter_for_resid_tol3

 solver_itermin = s% ctrl% solver_itermin
 solver_itermin_until_reduce_min_corr_coeff = s% ctrl% solver_itermin_until_reduce_min_corr_coeff
 solver_reduced_min_corr_coeff = s% ctrl% solver_reduced_min_corr_coeff
 do_solver_damping_for_neg_xa = s% ctrl% do_solver_damping_for_neg_xa
 scale_max_correction_for_negative_surf_lum = s% ctrl% scale_max_correction_for_negative_surf_lum
 max_frac_for_negative_surf_lum = s% ctrl% max_frac_for_negative_surf_lum
 hydro_mtx_max_allowed_abs_dlogT = s% ctrl% hydro_mtx_max_allowed_abs_dlogT
 hydro_mtx_max_allowed_abs_dlogRho = s% ctrl% hydro_mtx_max_allowed_abs_dlogRho
 min_logT_for_hydro_mtx_max_allowed = s% ctrl% min_logT_for_hydro_mtx_max_allowed
 hydro_mtx_max_allowed_logT = s% ctrl% hydro_mtx_max_allowed_logT
 hydro_mtx_max_allowed_logRho = s% ctrl% hydro_mtx_max_allowed_logRho
 hydro_mtx_min_allowed_logT = s% ctrl% hydro_mtx_min_allowed_logT
 hydro_mtx_min_allowed_logRho = s% ctrl% hydro_mtx_min_allowed_logRho
 
 use_DGESVX_in_bcyclic = s% ctrl% use_DGESVX_in_bcyclic
 use_equilibration_in_DGESVX = s% ctrl% use_equilibration_in_DGESVX
 report_min_rcond_from_DGESXV = s% ctrl% report_min_rcond_from_DGESXV
 
 op_split_burn = s% ctrl% op_split_burn
 op_split_burn_min_T = s% ctrl% op_split_burn_min_T
 op_split_burn_eps = s% ctrl% op_split_burn_eps
 op_split_burn_odescal = s% ctrl% op_split_burn_odescal
 op_split_burn_min_T_for_variable_T_solver = s% ctrl% op_split_burn_min_T_for_variable_T_solver

 tiny_corr_coeff_limit = s% ctrl% tiny_corr_coeff_limit
 scale_correction_norm = s% ctrl% scale_correction_norm
 num_times_solver_reuse_mtx = s% ctrl% num_times_solver_reuse_mtx
 corr_param_factor = s% ctrl% corr_param_factor
 scale_max_correction = s% ctrl% scale_max_correction
 ignore_min_corr_coeff_for_scale_max_correction = s% ctrl% ignore_min_corr_coeff_for_scale_max_correction
 ignore_too_large_correction = s% ctrl% ignore_too_large_correction
 ignore_species_in_max_correction = s% ctrl% ignore_species_in_max_correction

 corr_norm_jump_limit = s% ctrl% corr_norm_jump_limit
 max_corr_jump_limit = s% ctrl% max_corr_jump_limit
 resid_norm_jump_limit = s% ctrl% resid_norm_jump_limit
 max_resid_jump_limit = s% ctrl% max_resid_jump_limit

 corr_coeff_limit = s% ctrl% corr_coeff_limit
 tiny_corr_factor = s% ctrl% tiny_corr_factor

 solver_max_tries_before_reject = s% ctrl% solver_max_tries_before_reject
 max_tries1 = s% ctrl% max_tries1
 max_tries_for_retry = s% ctrl% max_tries_for_retry
 max_tries_after_5_retries = s% ctrl% max_tries_after_5_retries
 max_tries_after_10_retries = s% ctrl% max_tries_after_10_retries
 max_tries_after_20_retries = s% ctrl% max_tries_after_20_retries
 retry_limit = s% ctrl% retry_limit
 redo_limit = s% ctrl% redo_limit

 use_Pvsc_art_visc = s% ctrl% use_Pvsc_art_visc
 Pvsc_cq = s% ctrl% Pvsc_cq
 Pvsc_zsh = s% ctrl% Pvsc_zsh

 min_xa_hard_limit = s% ctrl% min_xa_hard_limit
 min_xa_hard_limit_for_highT = s% ctrl% min_xa_hard_limit_for_highT
 logT_max_for_min_xa_hard_limit = s% ctrl% logT_max_for_min_xa_hard_limit
 logT_min_for_min_xa_hard_limit_for_highT = s% ctrl% logT_min_for_min_xa_hard_limit_for_highT

 sum_xa_hard_limit = s% ctrl% sum_xa_hard_limit
 sum_xa_hard_limit_for_highT = s% ctrl% sum_xa_hard_limit_for_highT
 logT_max_for_sum_xa_hard_limit = s% ctrl% logT_max_for_sum_xa_hard_limit
 logT_min_for_sum_xa_hard_limit_for_highT = s% ctrl% logT_min_for_sum_xa_hard_limit_for_highT

 xa_clip_limit = s% ctrl% xa_clip_limit
 report_solver_progress = s% ctrl% report_solver_progress
 solver_test_partials_call_number = s% ctrl% solver_test_partials_call_number
 solver_test_partials_iter_number = s% ctrl% solver_test_partials_iter_number
 solver_epsder_chem = s% ctrl% solver_epsder_chem
 solver_epsder_struct = s% ctrl% solver_epsder_struct
 solver_numerical_jacobian = s% ctrl% solver_numerical_jacobian
 solver_jacobian_nzlo = s% ctrl% solver_jacobian_nzlo
 solver_jacobian_nzhi = s% ctrl% solver_jacobian_nzhi
 solver_check_everything = s% ctrl% solver_check_everything
 energy_conservation_dump_model_number = s% ctrl% energy_conservation_dump_model_number
 solver_inspect_soln_flag = s% ctrl% solver_inspect_soln_flag
 solver_test_partials_dx_0 = s% ctrl% solver_test_partials_dx_0
 solver_test_partials_k = s% ctrl% solver_test_partials_k
 solver_test_partials_k_low = s% ctrl% solver_test_partials_k_low
 solver_test_partials_k_high = s% ctrl% solver_test_partials_k_high
 solver_show_correction_info = s% ctrl% solver_show_correction_info
 solver_test_partials_write_eos_call_info = s% ctrl% solver_test_partials_write_eos_call_info
 solver_test_eos_partials = s% ctrl% solver_test_eos_partials
 solver_test_kap_partials = s% ctrl% solver_test_kap_partials
 solver_test_net_partials = s% ctrl% solver_test_net_partials
 solver_test_atm_partials = s% ctrl% solver_test_atm_partials
solver_test_partials_var_name = s% ctrl% solver_test_partials_var_name
solver_test_partials_sink_name = s% ctrl% solver_test_partials_sink_name
 solver_test_partials_equ_name = s% ctrl% solver_test_partials_equ_name
 solver_test_partials_show_dx_var_name = s% ctrl% solver_test_partials_show_dx_var_name
 solver_save_photo_call_number = s% ctrl% solver_save_photo_call_number
 fill_arrays_with_NaNs = s% ctrl% fill_arrays_with_NaNs
 zero_when_allocate = s% ctrl% zero_when_allocate
 warn_when_large_rel_run_E_err = s% ctrl% warn_when_large_rel_run_E_err
 warn_when_large_virial_thm_rel_err = s% ctrl% warn_when_large_virial_thm_rel_err
 warn_when_get_a_bad_eos_result = s% ctrl% warn_when_get_a_bad_eos_result
 warn_rates_for_high_temp = s% ctrl% warn_rates_for_high_temp
 max_safe_logT_for_rates = s% ctrl% max_safe_logT_for_rates
 eps_mdot_leak_frac_factor = s% ctrl% eps_mdot_leak_frac_factor

 alpha_TDC_DAMP = s% ctrl% alpha_TDC_DAMP
 alpha_TDC_DAMPR = s% ctrl% alpha_TDC_DAMPR
 alpha_TDC_PtdVdt = s% ctrl% alpha_TDC_PtdVdt
 compare_TDC_to_MLT = s% ctrl% compare_TDC_to_MLT

 RSP2_alfap= s% ctrl% RSP2_alfap
 RSP2_alfad = s% ctrl% RSP2_alfad
 RSP2_alfat= s% ctrl% RSP2_alfat 
 RSP2_alfam= s% ctrl% RSP2_alfam
 RSP2_alfar= s% ctrl% RSP2_alfar
 RSP2_min_Lt_div_L_for_overshooting_mixing_type = s% ctrl% RSP2_min_Lt_div_L_for_overshooting_mixing_type
 RSP2_min_Lc_div_L_for_convective_mixing_type = s% ctrl% RSP2_min_Lc_div_L_for_convective_mixing_type
 RSP2_Lsurf_factor= s% ctrl% RSP2_Lsurf_factor
 RSP2_use_Stellingwerf_Lr = s% ctrl% RSP2_use_Stellingwerf_Lr
 RSP2_remesh_when_load = s% ctrl% RSP2_remesh_when_load
 RSP2_use_L_eqn_at_surface = s% ctrl% RSP2_use_L_eqn_at_surface
 RSP2_report_adjust_w = s% ctrl% RSP2_report_adjust_w
 RSP2_assume_HSE = s% ctrl% RSP2_assume_HSE
 RSP2_use_RSP_eqn_for_Y_face = s% ctrl% RSP2_use_RSP_eqn_for_Y_face
 RSP2_use_mass_interp_face_values = s% ctrl% RSP2_use_mass_interp_face_values
 RSP2_num_outermost_cells_forced_nonturbulent = s% ctrl% RSP2_num_outermost_cells_forced_nonturbulent
 RSP2_num_innermost_cells_forced_nonturbulent = s% ctrl% RSP2_num_innermost_cells_forced_nonturbulent
 RSP2_T_anchor = s% ctrl% RSP2_T_anchor
 RSP2_dq_1_factor = s% ctrl% RSP2_dq_1_factor
 RSP2_nz = s% ctrl% RSP2_nz
 RSP2_nz_outer = s% ctrl% RSP2_nz_outer
 RSP2_nz_div_IBOTOM = s% ctrl% RSP2_nz_div_IBOTOM
 RSP2_target_steps_per_cycle = s% ctrl% RSP2_target_steps_per_cycle
 RSP2_max_num_periods = s% ctrl% RSP2_max_num_periods
 RSP2_work_period = s% ctrl% RSP2_work_period
 RSP2_map_first_period = s% ctrl% RSP2_map_first_period
 RSP2_map_last_period = s% ctrl% RSP2_map_last_period
 RSP2_min_max_R_for_periods = s% ctrl% RSP2_min_max_R_for_periods
 RSP2_GREKM_avg_abs_frac_new = s% ctrl% RSP2_GREKM_avg_abs_frac_new
 RSP2_GREKM_avg_abs_limit = s% ctrl% RSP2_GREKM_avg_abs_limit
 RSP2_map_zone_interval = s% ctrl% RSP2_map_zone_interval
 RSP2_work_filename = s% ctrl% RSP2_work_filename
 RSP2_map_columns_filename = s% ctrl% RSP2_map_columns_filename
 RSP2_map_filename = s% ctrl% RSP2_map_filename
 RSP2_map_history_filename = s% ctrl% RSP2_map_history_filename
 RSP2_write_map = s% ctrl% RSP2_write_map
 RSP2_w_min_for_damping = s% ctrl% RSP2_w_min_for_damping
 RSP2_source_seed = s% ctrl% RSP2_source_seed
 RSP2_w_fix_if_neg = s% ctrl% RSP2_w_fix_if_neg

 max_X_for_conv_timescale = s% ctrl% max_X_for_conv_timescale
 min_X_for_conv_timescale = s% ctrl% min_X_for_conv_timescale
 max_q_for_conv_timescale = s% ctrl% max_q_for_conv_timescale
 min_q_for_conv_timescale = s% ctrl% min_q_for_conv_timescale
 max_q_for_QHSE_timescale = s% ctrl% max_q_for_QHSE_timescale
 min_q_for_QHSE_timescale = s% ctrl% min_q_for_QHSE_timescale

 ! timestep
 max_timestep = s% ctrl% max_timestep
 max_years_for_timestep = s% ctrl% max_years_for_timestep

 hi_T_max_years_for_timestep = s% ctrl% hi_T_max_years_for_timestep
 max_timestep_hi_T_limit = s% ctrl% max_timestep_hi_T_limit

 min_timestep_factor = s% ctrl% min_timestep_factor
 max_timestep_factor = s% ctrl% max_timestep_factor
 max_timestep_factor_at_high_T = s% ctrl% max_timestep_factor_at_high_T
 min_logT_for_max_timestep_factor_at_high_T = s% ctrl% min_logT_for_max_timestep_factor_at_high_T
 time_delta_coeff = s% ctrl% time_delta_coeff
 timestep_factor_for_retries = s% ctrl% timestep_factor_for_retries
 retry_hold = s% ctrl% retry_hold
 neg_mass_fraction_hold = s% ctrl% neg_mass_fraction_hold
 timestep_dt_factor = s% ctrl% timestep_dt_factor
 use_dt_low_pass_controller = s% ctrl% use_dt_low_pass_controller
 
 force_timestep_min = s% ctrl% force_timestep_min
 force_timestep_min_years = s% ctrl% force_timestep_min_years
 force_timestep_min_factor = s% ctrl% force_timestep_min_factor
 force_timestep = s% ctrl% force_timestep
 force_timestep_years = s% ctrl% force_timestep_years

 varcontrol_target = s% ctrl% varcontrol_target
 min_allowed_varcontrol_target = s% ctrl% min_allowed_varcontrol_target
 varcontrol_dt_limit_ratio_hard_max = s% ctrl% varcontrol_dt_limit_ratio_hard_max
 xa_scale = s% ctrl% xa_scale

 solver_iters_timestep_limit = s% ctrl% solver_iters_timestep_limit

 burn_steps_limit = s% ctrl% burn_steps_limit
 burn_steps_hard_limit = s% ctrl% burn_steps_hard_limit

 diffusion_steps_limit = s% ctrl% diffusion_steps_limit
 diffusion_steps_hard_limit = s% ctrl% diffusion_steps_hard_limit
 diffusion_iters_limit = s% ctrl% diffusion_iters_limit
 diffusion_iters_hard_limit = s% ctrl% diffusion_iters_hard_limit

 dt_div_dt_cell_collapse_limit = s% ctrl% dt_div_dt_cell_collapse_limit
 dt_div_dt_cell_collapse_hard_limit = s% ctrl% dt_div_dt_cell_collapse_hard_limit
 dt_div_min_dr_div_cs_limit = s% ctrl% dt_div_min_dr_div_cs_limit
 dt_div_min_dr_div_cs_hard_limit = s% ctrl% dt_div_min_dr_div_cs_hard_limit
 
 min_abs_du_div_cs_for_dt_div_min_dr_div_cs_limit = s% ctrl% min_abs_du_div_cs_for_dt_div_min_dr_div_cs_limit
 min_abs_u_div_cs_for_dt_div_min_dr_div_cs_limit = s% ctrl% min_abs_u_div_cs_for_dt_div_min_dr_div_cs_limit
 min_k_for_dt_div_min_dr_div_cs_limit = s% ctrl% min_k_for_dt_div_min_dr_div_cs_limit
 min_q_for_dt_div_min_dr_div_cs_limit = s% ctrl% min_q_for_dt_div_min_dr_div_cs_limit
 max_q_for_dt_div_min_dr_div_cs_limit = s% ctrl% max_q_for_dt_div_min_dr_div_cs_limit
 check_remnant_only_for_dt_div_min_dr_div_cs_limit = s% ctrl% check_remnant_only_for_dt_div_min_dr_div_cs_limit

 dX_mix_dist_limit = s% ctrl% dX_mix_dist_limit

 dH_limit_min_H = s% ctrl% dH_limit_min_H
 dH_limit = s% ctrl% dH_limit
 dH_hard_limit = s% ctrl% dH_hard_limit
 dH_div_H_limit_min_H = s% ctrl% dH_div_H_limit_min_H
 dH_div_H_limit = s% ctrl% dH_div_H_limit
 dH_div_H_hard_limit = s% ctrl% dH_div_H_hard_limit
 dH_decreases_only = s% ctrl% dH_decreases_only

 dHe_limit_min_He = s% ctrl% dHe_limit_min_He
 dHe_limit = s% ctrl% dHe_limit
 dHe_hard_limit = s% ctrl% dHe_hard_limit
 dHe_div_He_limit_min_He = s% ctrl% dHe_div_He_limit_min_He
 dHe_div_He_limit = s% ctrl% dHe_div_He_limit
 dHe_div_He_hard_limit = s% ctrl% dHe_div_He_hard_limit
 dHe_decreases_only = s% ctrl% dHe_decreases_only

 dHe3_limit_min_He3 = s% ctrl% dHe3_limit_min_He3
 dHe3_limit = s% ctrl% dHe3_limit
 dHe3_hard_limit = s% ctrl% dHe3_hard_limit
 dHe3_div_He3_limit_min_He3 = s% ctrl% dHe3_div_He3_limit_min_He3
 dHe3_div_He3_limit = s% ctrl% dHe3_div_He3_limit
 dHe3_div_He3_hard_limit = s% ctrl% dHe3_div_He3_hard_limit
 dHe3_decreases_only = s% ctrl% dHe3_decreases_only

 dX_limit_min_X = s% ctrl% dX_limit_min_X
 dX_limit = s% ctrl% dX_limit
 dX_hard_limit = s% ctrl% dX_hard_limit
 dX_div_X_limit_min_X = s% ctrl% dX_div_X_limit_min_X
 dX_div_X_limit = s% ctrl% dX_div_X_limit
 dX_div_X_hard_limit = s% ctrl% dX_div_X_hard_limit
 dX_div_X_at_high_T_limit = s% ctrl% dX_div_X_at_high_T_limit
 dX_div_X_at_high_T_hard_limit = s% ctrl% dX_div_X_at_high_T_hard_limit
 dX_div_X_at_high_T_limit_lgT_min = s% ctrl% dX_div_X_at_high_T_limit_lgT_min
 dX_decreases_only = s% ctrl% dX_decreases_only

 dX_nuc_drop_min_X_limit = s% ctrl% dX_nuc_drop_min_X_limit
 dX_nuc_drop_max_A_limit = s% ctrl% dX_nuc_drop_max_A_limit
 dX_nuc_drop_limit = s% ctrl% dX_nuc_drop_limit
 dX_nuc_drop_limit_at_high_T = s% ctrl% dX_nuc_drop_limit_at_high_T
 dX_nuc_drop_hard_limit = s% ctrl% dX_nuc_drop_hard_limit
 dX_nuc_drop_min_yrs_for_dt = s% ctrl% dX_nuc_drop_min_yrs_for_dt

 dL_div_L_limit_min_L = s% ctrl% dL_div_L_limit_min_L
 dL_div_L_limit = s% ctrl% dL_div_L_limit
 dL_div_L_hard_limit = s% ctrl% dL_div_L_hard_limit

 delta_lgP_limit = s% ctrl% delta_lgP_limit
 delta_lgP_hard_limit = s% ctrl% delta_lgP_hard_limit
 delta_lgP_limit_min_lgP = s% ctrl% delta_lgP_limit_min_lgP

 delta_lgRho_limit = s% ctrl% delta_lgRho_limit
 delta_lgRho_hard_limit = s% ctrl% delta_lgRho_hard_limit
 delta_lgRho_limit_min_lgRho = s% ctrl% delta_lgRho_limit_min_lgRho

 delta_lgT_limit = s% ctrl% delta_lgT_limit
 delta_lgT_hard_limit = s% ctrl% delta_lgT_hard_limit
 delta_lgT_limit_min_lgT = s% ctrl% delta_lgT_limit_min_lgT

 delta_lgE_limit = s% ctrl% delta_lgE_limit
 delta_lgE_hard_limit = s% ctrl% delta_lgE_hard_limit
 delta_lgE_limit_min_lgE = s% ctrl% delta_lgE_limit_min_lgE

 delta_lgR_limit = s% ctrl% delta_lgR_limit
 delta_lgR_hard_limit = s% ctrl% delta_lgR_hard_limit
 delta_lgR_limit_min_lgR = s% ctrl% delta_lgR_limit_min_lgR

 delta_Ye_highT_limit = s% ctrl% delta_Ye_highT_limit
 delta_Ye_highT_hard_limit = s% ctrl% delta_Ye_highT_hard_limit
 minT_for_highT_Ye_limit = s% ctrl% minT_for_highT_Ye_limit

 delta_lgL_nuc_cat_limit = s% ctrl% delta_lgL_nuc_cat_limit
 delta_lgL_nuc_cat_hard_limit = s% ctrl% delta_lgL_nuc_cat_hard_limit
 lgL_nuc_cat_burn_min = s% ctrl% lgL_nuc_cat_burn_min
 lgL_nuc_mix_dist_limit = s% ctrl% lgL_nuc_mix_dist_limit

 delta_lgL_H_limit = s% ctrl% delta_lgL_H_limit
 delta_lgL_H_hard_limit = s% ctrl% delta_lgL_H_hard_limit
 lgL_H_burn_min = s% ctrl% lgL_H_burn_min
 lgL_H_drop_factor = s% ctrl% lgL_H_drop_factor
 lgL_H_burn_relative_limit = s% ctrl% lgL_H_burn_relative_limit

 delta_lgL_He_limit = s% ctrl% delta_lgL_He_limit
 delta_lgL_He_hard_limit = s% ctrl% delta_lgL_He_hard_limit
 lgL_He_burn_min = s% ctrl% lgL_He_burn_min
 lgL_He_drop_factor = s% ctrl% lgL_He_drop_factor
 lgL_He_burn_relative_limit = s% ctrl% lgL_He_burn_relative_limit

 delta_lgL_z_limit = s% ctrl% delta_lgL_z_limit
 delta_lgL_z_hard_limit = s% ctrl% delta_lgL_z_hard_limit
 lgL_z_burn_min = s% ctrl% lgL_z_burn_min
 lgL_z_drop_factor = s% ctrl% lgL_z_drop_factor
 lgL_z_burn_relative_limit = s% ctrl% lgL_z_burn_relative_limit

 delta_lgL_power_photo_limit = s% ctrl% delta_lgL_power_photo_limit
 delta_lgL_power_photo_hard_limit = s% ctrl% delta_lgL_power_photo_hard_limit
 lgL_power_photo_burn_min = s% ctrl% lgL_power_photo_burn_min
 lgL_power_photo_drop_factor = s% ctrl% lgL_power_photo_drop_factor
 min_lgT_for_lgL_power_photo_limit = s% ctrl% min_lgT_for_lgL_power_photo_limit

 delta_lgL_nuc_limit = s% ctrl% delta_lgL_nuc_limit
 delta_lgL_nuc_hard_limit = s% ctrl% delta_lgL_nuc_hard_limit
 delta_lgL_nuc_at_high_T_limit = s% ctrl% delta_lgL_nuc_at_high_T_limit
 delta_lgL_nuc_at_high_T_hard_limit = s% ctrl% delta_lgL_nuc_at_high_T_hard_limit
 delta_lgL_nuc_at_high_T_limit_lgT_min = s% ctrl% delta_lgL_nuc_at_high_T_limit_lgT_min
 
 max_lgT_for_lgL_nuc_limit = s% ctrl% max_lgT_for_lgL_nuc_limit
 lgL_nuc_burn_min = s% ctrl% lgL_nuc_burn_min
 lgL_nuc_drop_factor = s% ctrl% lgL_nuc_drop_factor

 delta_lgRho_cntr_limit = s% ctrl% delta_lgRho_cntr_limit
 delta_lgRho_cntr_hard_limit = s% ctrl% delta_lgRho_cntr_hard_limit

 delta_lgT_cntr_limit = s% ctrl% delta_lgT_cntr_limit
 delta_lgT_cntr_hard_limit = s% ctrl% delta_lgT_cntr_hard_limit
 delta_lgT_cntr_limit_only_after_near_zams = s% ctrl% delta_lgT_cntr_limit_only_after_near_zams

 delta_lgP_cntr_limit = s% ctrl% delta_lgP_cntr_limit
 delta_lgP_cntr_hard_limit = s% ctrl% delta_lgP_cntr_hard_limit

 delta_lgT_max_limit = s% ctrl% delta_lgT_max_limit
 delta_lgT_max_hard_limit = s% ctrl% delta_lgT_max_hard_limit
 delta_lgT_max_limit_lgT_min = s% ctrl% delta_lgT_max_limit_lgT_min
 delta_lgT_max_limit_only_after_near_zams = s% ctrl% delta_lgT_max_limit_only_after_near_zams

 delta_lgT_max_at_high_T_limit = s% ctrl% delta_lgT_max_at_high_T_limit
 delta_lgT_max_at_high_T_hard_limit = s% ctrl% delta_lgT_max_at_high_T_hard_limit
 delta_lgT_max_at_high_T_limit_lgT_min = s% ctrl% delta_lgT_max_at_high_T_limit_lgT_min

 delta_log_eps_nuc_limit = s% ctrl% delta_log_eps_nuc_limit
 delta_log_eps_nuc_hard_limit = s% ctrl% delta_log_eps_nuc_hard_limit

 delta_dX_div_X_cntr_min = s% ctrl% delta_dX_div_X_cntr_min
 delta_dX_div_X_cntr_max = s% ctrl% delta_dX_div_X_cntr_max
 delta_dX_div_X_cntr_limit = s% ctrl% delta_dX_div_X_cntr_limit
 delta_dX_div_X_cntr_hard_limit = s% ctrl% delta_dX_div_X_cntr_hard_limit

 delta_dX_div_X_drop_only = s% ctrl% delta_dX_div_X_drop_only
 delta_lg_XH_drop_only = s% ctrl% delta_lg_XH_drop_only
 delta_lg_XHe_drop_only = s% ctrl% delta_lg_XHe_drop_only
 delta_lg_XC_drop_only = s% ctrl% delta_lg_XC_drop_only
 delta_lg_XNe_drop_only = s% ctrl% delta_lg_XNe_drop_only
 delta_lg_XO_drop_only = s% ctrl% delta_lg_XO_drop_only
 delta_lg_XSi_drop_only = s% ctrl% delta_lg_XSi_drop_only
 delta_XH_drop_only = s% ctrl% delta_XH_drop_only
 delta_XHe_drop_only = s% ctrl% delta_XHe_drop_only
 delta_XC_drop_only = s% ctrl% delta_XC_drop_only
 delta_XNe_drop_only = s% ctrl% delta_XNe_drop_only
 delta_XO_drop_only = s% ctrl% delta_XO_drop_only
 delta_XSi_drop_only = s% ctrl% delta_XSi_drop_only

 delta_lg_XH_cntr_min = s% ctrl% delta_lg_XH_cntr_min
 delta_lg_XH_cntr_max = s% ctrl% delta_lg_XH_cntr_max
 delta_lg_XH_cntr_limit = s% ctrl% delta_lg_XH_cntr_limit
 delta_lg_XH_cntr_hard_limit = s% ctrl% delta_lg_XH_cntr_hard_limit

 delta_lg_XHe_cntr_min = s% ctrl% delta_lg_XHe_cntr_min
 delta_lg_XHe_cntr_max = s% ctrl% delta_lg_XHe_cntr_max
 delta_lg_XHe_cntr_limit = s% ctrl% delta_lg_XHe_cntr_limit
 delta_lg_XHe_cntr_hard_limit = s% ctrl% delta_lg_XHe_cntr_hard_limit

 delta_lg_XC_cntr_min = s% ctrl% delta_lg_XC_cntr_min
 delta_lg_XC_cntr_max = s% ctrl% delta_lg_XC_cntr_max
 delta_lg_XC_cntr_limit = s% ctrl% delta_lg_XC_cntr_limit
 delta_lg_XC_cntr_hard_limit = s% ctrl% delta_lg_XC_cntr_hard_limit

 delta_lg_XNe_cntr_limit = s% ctrl% delta_lg_XNe_cntr_limit
 delta_lg_XNe_cntr_hard_limit = s% ctrl% delta_lg_XNe_cntr_hard_limit
 delta_lg_XNe_cntr_min = s% ctrl% delta_lg_XNe_cntr_min
 delta_lg_XNe_cntr_max = s% ctrl% delta_lg_XNe_cntr_max

 delta_lg_XO_cntr_limit = s% ctrl% delta_lg_XO_cntr_limit
 delta_lg_XO_cntr_hard_limit = s% ctrl% delta_lg_XO_cntr_hard_limit
 delta_lg_XO_cntr_min = s% ctrl% delta_lg_XO_cntr_min
 delta_lg_XO_cntr_max = s% ctrl% delta_lg_XO_cntr_max

 delta_lg_XSi_cntr_limit = s% ctrl% delta_lg_XSi_cntr_limit
 delta_lg_XSi_cntr_hard_limit = s% ctrl% delta_lg_XSi_cntr_hard_limit
 delta_lg_XSi_cntr_min = s% ctrl% delta_lg_XSi_cntr_min
 delta_lg_XSi_cntr_max = s% ctrl% delta_lg_XSi_cntr_max

 delta_XH_cntr_limit = s% ctrl% delta_XH_cntr_limit
 delta_XH_cntr_hard_limit = s% ctrl% delta_XH_cntr_hard_limit
 delta_XHe_cntr_limit = s% ctrl% delta_XHe_cntr_limit
 delta_XHe_cntr_hard_limit = s% ctrl% delta_XHe_cntr_hard_limit
 delta_XC_cntr_limit = s% ctrl% delta_XC_cntr_limit
 delta_XC_cntr_hard_limit = s% ctrl% delta_XC_cntr_hard_limit
 delta_XNe_cntr_limit = s% ctrl% delta_XNe_cntr_limit
 delta_XNe_cntr_hard_limit = s% ctrl% delta_XNe_cntr_hard_limit
 delta_XO_cntr_limit = s% ctrl% delta_XO_cntr_limit
 delta_XO_cntr_hard_limit = s% ctrl% delta_XO_cntr_hard_limit
 delta_XSi_cntr_limit = s% ctrl% delta_XSi_cntr_limit
 delta_XSi_cntr_hard_limit = s% ctrl% delta_XSi_cntr_hard_limit

 delta_lgTeff_limit = s% ctrl% delta_lgTeff_limit
 delta_lgTeff_hard_limit = s% ctrl% delta_lgTeff_hard_limit

 delta_lgL_limit = s% ctrl% delta_lgL_limit
 delta_lgL_limit_L_min = s% ctrl% delta_lgL_limit_L_min
 delta_lgL_hard_limit = s% ctrl% delta_lgL_hard_limit

 delta_HR_ds_L = s% ctrl% delta_HR_ds_L
 delta_HR_ds_Teff = s% ctrl% delta_HR_ds_Teff
 delta_HR_limit = s% ctrl% delta_HR_limit
 delta_HR_hard_limit = s% ctrl% delta_HR_hard_limit

 delta_lg_star_mass_limit = s% ctrl% delta_lg_star_mass_limit
 delta_lg_star_mass_hard_limit = s% ctrl% delta_lg_star_mass_hard_limit

 delta_mdot_atol = s% ctrl% delta_mdot_atol
 delta_mdot_rtol = s% ctrl% delta_mdot_rtol
 delta_mdot_limit = s% ctrl% delta_mdot_limit
 delta_mdot_hard_limit = s% ctrl% delta_mdot_hard_limit

 adjust_J_q_limit = s% ctrl% adjust_J_q_limit
 adjust_J_q_hard_limit = s% ctrl% adjust_J_q_hard_limit
 never_skip_hard_limits = s% ctrl% never_skip_hard_limits
 relax_hard_limits_after_retry = s% ctrl% relax_hard_limits_after_retry
 report_dt_hard_limit_retries = s% ctrl% report_dt_hard_limit_retries
 report_min_dr_div_cs = s% ctrl% report_min_dr_div_cs
 report_solver_dt_info = s% ctrl% report_solver_dt_info

 limit_for_rel_error_in_energy_conservation = s% ctrl% limit_for_rel_error_in_energy_conservation
 hard_limit_for_rel_error_in_energy_conservation = s% ctrl% hard_limit_for_rel_error_in_energy_conservation

 min_chem_eqn_scale = s% ctrl% min_chem_eqn_scale

 ! controls for the evolve routine
 trace_evolve = s% ctrl% trace_evolve


 ! misc
 zams_filename = s% ctrl% zams_filename
 set_rho_to_dm_div_dV = s% ctrl% set_rho_to_dm_div_dV

 use_other_mlt_results = s% ctrl% use_other_mlt_results
 use_other_surface_PT = s% ctrl% use_other_surface_PT
 use_other_kap = s% ctrl% use_other_kap
 use_other_diffusion = s% ctrl% use_other_diffusion
 use_other_diffusion_factor = s% ctrl% use_other_diffusion_factor
 use_other_adjust_mdot = s% ctrl% use_other_adjust_mdot
 use_other_j_for_adjust_J_lost = s% ctrl% use_other_j_for_adjust_J_lost
 use_other_alpha_mlt = s% ctrl% use_other_alpha_mlt
 use_other_am_mixing = s% ctrl% use_other_am_mixing
 use_other_brunt = s% ctrl% use_other_brunt
 use_other_brunt_smoothing = s% ctrl% use_other_brunt_smoothing
 use_other_solver_monitor = s% ctrl% use_other_solver_monitor
 use_other_build_initial_model = s% ctrl% use_other_build_initial_model
 use_other_cgrav = s% ctrl% use_other_cgrav
 use_other_mesh_delta_coeff_factor = s% ctrl% use_other_mesh_delta_coeff_factor
 use_other_energy_implicit = s% ctrl% use_other_energy_implicit
 use_other_momentum_implicit = s% ctrl% use_other_momentum_implicit
 use_other_momentum = s% ctrl% use_other_momentum
 use_other_remove_surface = s% ctrl% use_other_remove_surface
 use_other_energy = s% ctrl% use_other_energy
 use_other_pressure = s% ctrl% use_other_pressure
 use_other_mesh_functions = s% ctrl% use_other_mesh_functions
 use_other_eps_grav = s% ctrl% use_other_eps_grav
 use_other_gradr_factor = s% ctrl% use_other_gradr_factor
 use_other_D_mix = s% ctrl% use_other_D_mix
 use_other_neu = s% ctrl% use_other_neu
 use_other_net_get = s% ctrl% use_other_net_get
 use_other_opacity_factor = s% ctrl% use_other_opacity_factor
 use_other_diffusion_coefficients = s% ctrl% use_other_diffusion_coefficients
 use_other_pgstar_plots = s% ctrl% use_other_pgstar_plots
 use_other_eval_fp_ft = s% ctrl% use_other_eval_fp_ft
 use_other_eval_i_rot = s% ctrl% use_other_eval_i_rot
 use_other_torque = s% ctrl% use_other_torque
 use_other_torque_implicit = s% ctrl% use_other_torque_implicit
 use_other_wind = s% ctrl% use_other_wind
 use_other_accreting_state = s% ctrl% use_other_accreting_state
 use_other_after_struct_burn_mix = s% ctrl% use_other_after_struct_burn_mix
 use_other_before_struct_burn_mix = s% ctrl% use_other_before_struct_burn_mix
 use_other_astero_freq_corr = s% ctrl% use_other_astero_freq_corr
 use_other_timestep_limit = s% ctrl% use_other_timestep_limit
 use_other_set_pgstar_controls = s% ctrl% use_other_set_pgstar_controls
 use_other_screening = s% ctrl% use_other_screening

 x_ctrl = s% ctrl% x_ctrl
 x_integer_ctrl = s% ctrl% x_integer_ctrl
 x_logical_ctrl = s% ctrl% x_logical_ctrl
 x_character_ctrl = s% ctrl% x_character_ctrl

 ! info for debugging
 stop_for_bad_nums = s% ctrl% stop_for_bad_nums
 report_ierr = s% ctrl% report_ierr
 report_bad_negative_xa = s% ctrl% report_bad_negative_xa

 diffusion_dump_call_number = s% ctrl% diffusion_dump_call_number

 surface_accel_div_grav_limit = s% ctrl% surface_accel_div_grav_limit
 gradT_excess_age_fraction = s% ctrl% gradT_excess_age_fraction
 gradT_excess_max_change = s% ctrl% gradT_excess_max_change
 hot_wind_scheme = s% ctrl% hot_wind_scheme
 cool_wind_full_on_T = s% ctrl% cool_wind_full_on_T
 hot_wind_full_on_T = s% ctrl% hot_wind_full_on_T
 num_cells_for_smooth_brunt_B = s% ctrl% num_cells_for_smooth_brunt_B
 steps_before_start_stress_test = s% ctrl% steps_before_start_stress_test
 stress_test_relax = s% ctrl% stress_test_relax
 


 end subroutine set_controls_for_writing

   subroutine get_control(s, name, val, ierr)
      use utils_lib, only: StrUpCase
      type (star_info), pointer :: s
      character(len=*),intent(in) :: name
      character(len=*), intent(out) :: val
      integer, intent(out) :: ierr

      character(len(name)+1) :: upper_name
      character(len=512) :: str
      integer :: iounit,iostat,ind,i


      ! First save current controls
      call set_controls_for_writing(s, ierr)
      if(ierr/=0) return

      ! Write namelist to temporay file
      open(newunit=iounit,status='scratch')
      write(iounit,nml=controls)
      rewind(iounit)

      ! Namelists get written in captials
      upper_name = trim(StrUpCase(name))//'='
      val = ''
      ! Search for name inside namelist
      do 
         read(iounit,'(A)',iostat=iostat) str
         ind = index(trim(str),trim(upper_name))
         if( ind /= 0 ) then
            val = str(ind+len_trim(upper_name):len_trim(str)-1) ! Remove final comma and starting =
            do i=1,len(val)
               if(val(i:i)=='"') val(i:i) = ' '
            end do
            exit
         end if
         if(is_iostat_end(iostat)) exit
      end do   

      if(len_trim(val) == 0 .and. ind==0 ) ierr = -1

      close(iounit)

   end subroutine get_control

   subroutine set_control(s, name, val, ierr)
      type (star_info), pointer :: s
      character(len=*), intent(in) :: name, val
      character(len=len(name)+len(val)+13) :: tmp
      integer, intent(out) :: ierr

      ! First save current controls
      call set_controls_for_writing(s, ierr)
      if(ierr/=0) return

      tmp=''
      tmp = '&controls '//trim(name)//'='//trim(val)//' /'

      ! Load into namelist
      read(tmp, nml=controls)

      ! Add to star
      call store_controls(s, ierr)
      if(ierr/=0) return

   end subroutine set_control


 end module ctrls_io

